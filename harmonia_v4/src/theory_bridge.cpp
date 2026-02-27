#include "theory_bridge.h"
#include <cstring>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <errno.h>
#include <sys/select.h>

TheoryBridge::TheoryBridge()  {}
TheoryBridge::~TheoryBridge() { stop(); }

bool TheoryBridge::start(const std::string& theory_dir) {
    if (pipe(to_py_) < 0 || pipe(from_py_) < 0) { perror("pipe"); return false; }
    py_pid_ = fork();
    if (py_pid_ < 0) { perror("fork"); return false; }
    if (py_pid_ == 0) {
        dup2(to_py_[0], STDIN_FILENO); dup2(from_py_[1], STDOUT_FILENO);
        close(to_py_[0]); close(to_py_[1]);
        close(from_py_[0]); close(from_py_[1]);
        std::string srv = theory_dir + "/server.py";
        execl("/usr/bin/python3","python3",srv.c_str(),nullptr);
        perror("execl"); _exit(1);
    }
    close(to_py_[0]); close(from_py_[1]);
    to_py_[0] = from_py_[1] = -1;
    running_.store(true);
    reader_thread_ = std::thread(&TheoryBridge::readerLoop, this);
    fprintf(stderr,"TheoryBridge: Python pid=%d\n",py_pid_);
    return true;
}

void TheoryBridge::stop() {
    running_.store(false);
    if (py_pid_ > 0) { kill(py_pid_, SIGTERM); waitpid(py_pid_,nullptr,0); py_pid_=-1; }
    if (to_py_[1]   >= 0) { close(to_py_[1]);   to_py_[1]=-1; }
    if (from_py_[0] >= 0) { close(from_py_[0]); from_py_[0]=-1; }
    if (reader_thread_.joinable()) reader_thread_.join();
}

void TheoryBridge::readerLoop() {
    char buf[65536]; std::string acc;
    while (running_.load()) {
        fd_set fds; FD_ZERO(&fds); FD_SET(from_py_[0],&fds);
        struct timeval tv{0,50000};
        if (select(from_py_[0]+1,&fds,nullptr,nullptr,&tv) <= 0) continue;
        ssize_t n = read(from_py_[0],buf,sizeof(buf)-1);
        if (n<=0) { running_.store(false); break; }
        buf[n] = '\0'; acc += buf;
        size_t pos;
        while ((pos = acc.find('\n')) != std::string::npos) {
            std::string line = acc.substr(0,pos); acc.erase(0,pos+1);
            if (!line.empty()) {
                std::lock_guard<std::mutex> lk(response_mutex_);
                responses_.push_back(std::move(line));
            }
        }
    }
}

void TheoryBridge::sendRaw(const std::string& json) {
    if (to_py_[1]<0 || !running_.load()) return;
    std::string msg = json + "\n";
    ssize_t n = write(to_py_[1],msg.c_str(),msg.size()); (void)n;
}

void TheoryBridge::poll() {
    std::deque<std::string> local;
    { std::lock_guard<std::mutex> lk(response_mutex_); local.swap(responses_); }
    for (auto& resp : local) {
        std::string tag = jStr(resp,"tag");
        if (tag.empty() || tag == "?") tag = jStr(resp,"result");

        RawJsonCb cb = nullptr;
        {
            std::lock_guard<std::mutex> lk(pending_mutex_);
            for (auto it = pending_cbs_.begin(); it != pending_cbs_.end(); ++it) {
                if (it->result_tag == tag) {
                    cb = it->raw_cb;
                    pending_cbs_.erase(it);
                    break;
                }
            }
        }
        if (cb) { cb(resp); }
    }
}

std::string TheoryBridge::sendSync(const std::string& json, int timeout_ms) {
    sendRaw(json); int waited = 0;
    while (waited < timeout_ms) {
        usleep(5000); waited += 5;
        std::lock_guard<std::mutex> lk(response_mutex_);
        if (!responses_.empty()) { auto r=responses_.front(); responses_.pop_front(); return r; }
    }
    return "{}";
}

// ────────────────────────────────────────────────────────────────────────────
//  Public structural queries
// ────────────────────────────────────────────────────────────────────────────
void TheoryBridge::queryAnalyzeChord(const std::vector<int>& pcs, int key, int edo,
                                      FunctionCb cb) {
    std::string arr = "[";
    for (size_t i=0;i<pcs.size();i++) { arr+=std::to_string(pcs[i]); if(i+1<pcs.size()) arr+=","; }
    arr += "]";
    std::string json = "{\"cmd\":\"analyze_chord\",\"tag\":\"analyze_chord\",\"pcs\":" + arr +
                       ",\"key\":" + std::to_string(key) + ",\"edo\":" + std::to_string(edo) + "}";
    auto raw_cb = [cb](const std::string& resp) {
        FunctionalAnalysis fa;
        TonnetzTension tt;
        fa.function        = jStr(resp,"function");
        fa.function_reason = jStr(resp,"function_reason");
        fa.tension_level   = jInt(resp,"tension_level");
        fa.has_tritone     = jBool(resp,"has_tritone");
        fa.in_diatonic     = jBool(resp,"in_diatonic");
        fa.tritone         = jIntArr(resp,"tritone");
        fa.tritone_names   = jStrArr(resp,"tritone_names");
        fa.chromatic_tones = jIntArr(resp,"chromatic_tones");
        // Parse tendency tones
        size_t ta = resp.find("\"tendency_tones\"");
        if (ta != std::string::npos) {
            size_t as = resp.find('[',ta), ae = resp.find(']',as);
            if (as!=std::string::npos&&ae!=std::string::npos) {
                std::string arr_str = resp.substr(as+1,ae-as-1);
                for (auto& obj : splitJsonArray(arr_str)) {
                    TendencyTone t;
                    t.pc         = jInt(obj,"pc");
                    t.name       = jStr(obj,"name");
                    t.role       = jStr(obj,"role");
                    t.tendency   = jStr(obj,"tendency");
                    t.force      = jStr(obj,"force");
                    t.derivation = jStr(obj,"derivation");
                    fa.tendency_tones.push_back(t);
                }
            }
        }
        // Parse tension block
        tt.tension       = jFloat(resp,"tension");
        tt.tension_label = jStr(resp,"tension_label");
        tt.distance      = jFloat(resp,"distance");
        if (cb) cb(fa, tt);
    };
    { std::lock_guard<std::mutex> lk(pending_mutex_);
      pending_cbs_.push_back({"analyze_chord", raw_cb}); }
    sendRaw(json);
}

void TheoryBridge::queryPLRNeighbors(int root, const std::string& quality,
                                      PLRNeighborsCb cb) {
    std::string json = "{\"cmd\":\"plr_neighbors\",\"tag\":\"plr_neighbors\",\"root\":" + std::to_string(root) +
                       ",\"quality\":\"" + quality + "\"}";
    auto raw_cb = [cb](const std::string& resp) {
        std::vector<PLRTransform> nbrs;
        size_t na = resp.find("\"neighbors\"");
        if (na == std::string::npos) { if(cb) cb(nbrs); return; }
        size_t as = resp.find('[',na), ae = resp.rfind(']');
        if (as==std::string::npos||ae==std::string::npos) { if(cb) cb(nbrs); return; }
        std::string arr_str = resp.substr(as+1,ae-as-1);
        for (auto& obj : splitJsonArray(arr_str))
            nbrs.push_back(parsePLRTransform(obj));
        if (cb) cb(nbrs);
    };
    { std::lock_guard<std::mutex> lk(pending_mutex_);
      pending_cbs_.push_back({"plr_neighbors", raw_cb}); }
    sendRaw(json);
}

void TheoryBridge::queryPLRPath(int ra, const std::string& qa,
                                 int rb, const std::string& qb,
                                 PLRPathCb cb) {
    std::string json = "{\"cmd\":\"plr_path\",\"tag\":\"plr_path\",\"root_a\":" + std::to_string(ra) +
                       ",\"quality_a\":\"" + qa + "\",\"root_b\":" + std::to_string(rb) +
                       ",\"quality_b\":\"" + qb + "\"}";
    auto raw_cb = [cb](const std::string& resp) {
        bool reachable = jBool(resp,"reachable",true);
        auto path = jStrArr(resp,"path");
        if (cb) cb(path, reachable);
    };
    { std::lock_guard<std::mutex> lk(pending_mutex_);
      pending_cbs_.push_back({"plr_path", raw_cb}); }
    sendRaw(json);
}

void TheoryBridge::queryResolutionPaths(int root, const std::string& quality, int key, int edo,
                                         ResolutionsCb cb) {
    std::string json = "{\"cmd\":\"resolution_paths\",\"tag\":\"resolution_paths\",\"root\":" + std::to_string(root) +
                       ",\"quality\":\"" + quality + "\",\"key\":" + std::to_string(key) + ",\"edo\":" + std::to_string(edo) + "}";
    auto raw_cb = [cb](const std::string& resp) {
        std::vector<ResolutionPath> paths;
        size_t pa = resp.find("\"paths\"");
        if (pa == std::string::npos) { if(cb) cb(paths); return; }
        size_t as = resp.find('[',pa), ae = resp.rfind(']');
        if (as==std::string::npos||ae==std::string::npos) { if(cb) cb(paths); return; }
        std::string arr_str = resp.substr(as+1,ae-as-1);
        for (auto& obj : splitJsonArray(arr_str))
            paths.push_back(parseResolutionPath(obj));
        if (cb) cb(paths);
    };
    { std::lock_guard<std::mutex> lk(pending_mutex_);
      pending_cbs_.push_back({"resolution_paths", raw_cb}); }
    sendRaw(json);
}

void TheoryBridge::queryPsychoacoustic(const std::vector<int>& pcs, int key,
                                       float c4_hz, int octave, float level_db,
                                       PsychoacousticCb cb) {
    std::string arr = "[";
    for (size_t i=0;i<pcs.size();i++) { arr+=std::to_string(pcs[i]); if(i+1<pcs.size()) arr+=","; }
    arr += "]";
    std::string json = "{\"cmd\":\"psychoacoustic_analysis\",\"tag\":\"psychoacoustic\",\"pcs\":" + arr +
                       ",\"key\":" + std::to_string(key) + ",\"c4_hz\":" + std::to_string(c4_hz) +
                       ",\"octave\":" + std::to_string(octave) + ",\"level_db\":" + std::to_string(level_db) + "}";
    auto raw_cb = [cb](const std::string& resp) {
        if (cb) cb(parsePsychoacoustic(resp));
    };
    { std::lock_guard<std::mutex> lk(pending_mutex_);
      pending_cbs_.push_back({"psychoacoustic", raw_cb}); }
    sendRaw(json);
}

void TheoryBridge::querySuggestCompletion(const std::vector<int>& pcs, int key, int edo,
                                           CompletionCb cb) {
    std::string arr = "[";
    for (size_t i=0;i<pcs.size();i++) { arr+=std::to_string(pcs[i]); if(i+1<pcs.size()) arr+=","; }
    arr += "]";
    std::string json = "{\"cmd\":\"suggest_completion\",\"tag\":\"suggest_completion\",\"pitch_classes\":" + arr +
                       ",\"key\":" + std::to_string(key) + ",\"edo\":" + std::to_string(edo) + "}";
    auto raw_cb = [cb](const std::string& resp) {
        std::vector<CompletionSuggestion> sugs;
        size_t sa = resp.find("\"suggestions\"");
        if (sa == std::string::npos) { if(cb) cb(sugs); return; }
        size_t as = resp.find('[',sa), ae = resp.rfind(']');
        if (as==std::string::npos||ae==std::string::npos) { if(cb) cb(sugs); return; }
        std::string arr_str = resp.substr(as+1,ae-as-1);
        for (auto& obj : splitJsonArray(arr_str))
            sugs.push_back(parseCompletion(obj));
        if (cb) cb(sugs);
    };
    { std::lock_guard<std::mutex> lk(pending_mutex_);
      pending_cbs_.push_back({"suggest_completion", raw_cb}); }
    sendRaw(json);
}

void TheoryBridge::queryOrbifoldDistance(const std::vector<int>& a,
                                          const std::vector<int>& b,
                                          OrbifoldCb cb) {
    auto mkArr = [](const std::vector<int>& v) {
        std::string s="["; for(size_t i=0;i<v.size();i++){s+=std::to_string(v[i]);if(i+1<v.size())s+=",";}; return s+"]";
    };
    std::string json = "{\"cmd\":\"orbifold_distance\",\"tag\":\"orbifold_distance\",\"chord_a\":" + mkArr(a) +
                       ",\"chord_b\":" + mkArr(b) + "}";
    auto raw_cb = [cb](const std::string& resp) {
        if (cb) cb(parseVoiceLeading(resp));
    };
    { std::lock_guard<std::mutex> lk(pending_mutex_);
      pending_cbs_.push_back({"orbifold_distance", raw_cb}); }
    sendRaw(json);
}

void TheoryBridge::queryPivotSearch(int key_from, int key_to, PivotSearchCb cb) {
    std::string json = "{\"cmd\":\"pivot_search\",\"tag\":\"pivot_search\",\"key_from\":" + std::to_string(key_from) +
                       ",\"key_to\":" + std::to_string(key_to) + "}";
    auto raw_cb = [this, cb](const std::string& resp) {
        PivotSearchResult res;
        res.key_from         = jInt(resp, "key_from");
        res.key_to           = jInt(resp, "key_to");
        res.key_from_name    = jStr(resp, "key_from_name");
        res.key_to_name      = jStr(resp, "key_to_name");
        res.cof_distance     = jInt(resp, "cof_distance");
        res.cof_relationship = jStr(resp, "cof_relationship");
        res.common_scale_tones = jIntArr(resp, "common_scale_tones");

        // Parse all_pivots
        size_t ap = resp.find("\"all_pivots\"");
        if (ap != std::string::npos) {
            size_t as = resp.find('[', ap), ae = resp.find(']', as);
            // This simple parser might fail on nested arrays, but splitJsonArray helps
            if (as != std::string::npos) {
                // Find matching ] for the array
                int depth = 0; size_t end = as;
                for (size_t i=as; i<resp.size(); i++) {
                    if (resp[i] == '[') depth++;
                    else if (resp[i] == ']') {
                        depth--;
                        if (depth == 0) { end = i; break; }
                    }
                }
                std::string arr_str = resp.substr(as + 1, end - as - 1);
                for (auto& obj : splitJsonArray(arr_str)) {
                    PivotChord pc;
                    pc.type           = jStr(obj, "type");
                    pc.pcs            = jIntArr(obj, "pcs");
                    pc.label          = jStr(obj, "label");
                    pc.root           = jInt(obj, "root");
                    pc.quality        = jStr(obj, "quality");
                    pc.roman_from     = jStr(obj, "roman_from");
                    pc.roman_to       = jStr(obj, "roman_to");
                    pc.function_from  = jStr(obj, "function_from");
                    pc.function_to    = jStr(obj, "function_to");
                    pc.pivot_score    = jFloat(obj, "pivot_score");
                    pc.interpretation = jStr(obj, "interpretation");
                    res.all_pivots.push_back(pc);
                }
            }
        }

        // Parse modulation_path
        size_t mp = resp.find("\"modulation_path\"");
        if (mp != std::string::npos) {
            size_t as = resp.find('[', mp);
            if (as != std::string::npos) {
                int depth = 0; size_t end = as;
                for (size_t i=as; i<resp.size(); i++) {
                    if (resp[i] == '[') depth++;
                    else if (resp[i] == ']') {
                        depth--;
                        if (depth == 0) { end = i; break; }
                    }
                }
                std::string arr_str = resp.substr(as + 1, end - as - 1);
                for (auto& obj : splitJsonArray(arr_str)) {
                    ModulationStep s;
                    s.root           = jInt(obj, "root");
                    s.quality        = jStr(obj, "quality");
                    s.label          = jStr(obj, "label");
                    s.roman          = jStr(obj, "roman");
                    s.key_context    = jStr(obj, "key_context");
                    s.function       = jStr(obj, "function");
                    s.vl_from_prev   = jInt(obj, "vl_from_prev");
                    s.cumulative_cost= jInt(obj, "cumulative_cost");
                    s.pivot_note     = jStr(obj, "pivot_note");
                    s.roman_as_pivot_from = jStr(obj, "roman_as_pivot_from");
                    s.roman_as_pivot_to   = jStr(obj, "roman_as_pivot_to");
                    res.modulation_path.push_back(s);
                }
            }
        }

        if (cb) cb(res);
    };
    { std::lock_guard<std::mutex> lk(pending_mutex_);
      pending_cbs_.push_back({"pivot_search", raw_cb}); }
    sendRaw(json);
}

void TheoryBridge::queryDetectSequence(
    const std::vector<std::pair<int,std::string>>& chords, SequenceCb cb) {
    std::string arr = "[";
    for (size_t i=0;i<chords.size();i++) {
        arr += "{\"root\":" + std::to_string(chords[i].first) +
               ",\"quality\":\"" + chords[i].second + "\"}";
        if (i+1 < chords.size()) arr += ",";
    }
    arr += "]";
    std::string json = "{\"cmd\":\"detect_sequence\",\"tag\":\"detect_sequence\",\"chords\":" + arr + "}";
    auto raw_cb = [cb](const std::string& resp) {
        SequenceInfo si;
        si.type                    = jStr(resp,"sequence_type");
        si.interval                = jInt(resp,"interval");
        si.interval_name           = jStr(resp,"interval_name");
        si.description             = jStr(resp,"description");
        si.algebraic_explanation   = jStr(resp,"algebraic_explanation");
        si.period                  = jInt(resp,"period",0);
        if (cb) cb(si);
    };
    { std::lock_guard<std::mutex> lk(pending_mutex_);
      pending_cbs_.push_back({"detect_sequence", raw_cb}); }
    sendRaw(json);
}

void TheoryBridge::queryEDOAnalysis(int edo, RawJsonCb cb) {
    std::string json = "{\"cmd\":\"edo_analysis\",\"tag\":\"edo_analysis\",\"edo\":" + std::to_string(edo) + "}";
    auto raw_cb = [cb](const std::string& resp) { if(cb) cb(resp); };
    { std::lock_guard<std::mutex> lk(pending_mutex_);
      pending_cbs_.push_back({"edo_analysis", raw_cb}); }
    sendRaw(json);
}

void TheoryBridge::queryRaw(const std::string& json, const std::string& tag, RawJsonCb cb) {
    auto raw_cb = [cb](const std::string& r) { if(cb) cb(r); };
    { std::lock_guard<std::mutex> lk(pending_mutex_);
      pending_cbs_.push_back({tag, raw_cb}); }
    sendRaw(json);
}

// ────────────────────────────────────────────────────────────────────────────
//  JSON helpers (minimal, no external deps)
// ────────────────────────────────────────────────────────────────────────────
static size_t findValPos(const std::string& j, const std::string& k) {
    std::string needle = "\"" + k + "\"";
    size_t p = j.find(needle);
    if (p == std::string::npos) return std::string::npos;
    p = j.find(':', p + needle.size());
    if (p == std::string::npos) return std::string::npos;
    p++; // skip colon
    while (p < j.size() && (j[p] == ' ' || j[p] == '\t' || j[p] == '\n' || j[p] == '\r')) p++;
    return p;
}

int TheoryBridge::jInt(const std::string& j, const std::string& k, int d) {
    size_t p = findValPos(j, k);
    if (p == std::string::npos) return d;
    try { return std::stoi(j.substr(p)); } catch(...) { return d; }
}

float TheoryBridge::jFloat(const std::string& j, const std::string& k, float d) {
    size_t p = findValPos(j, k);
    if (p == std::string::npos) return d;
    try { return std::stof(j.substr(p)); } catch(...) { return d; }
}

std::string TheoryBridge::jStr(const std::string& j, const std::string& k,
                                const std::string& d) {
    size_t p = findValPos(j, k);
    if (p == std::string::npos || j[p] != '"') return d;
    p++; // skip quote
    size_t e = p;
    while (e < j.size() && !(j[e] == '"' && (e == 0 || j[e-1] != '\\'))) e++;
    return (e < j.size()) ? j.substr(p, e - p) : d;
}

bool TheoryBridge::jBool(const std::string& j, const std::string& k, bool d) {
    size_t p = findValPos(j, k);
    if (p == std::string::npos) return d;
    if (j[p] == 't') return true;
    if (j[p] == 'f') return false;
    return d;
}

std::vector<int> TheoryBridge::jIntArr(const std::string& j, const std::string& k) {
    std::vector<int> result;
    size_t p = findValPos(j, k);
    if (p == std::string::npos || j[p] != '[') return result;
    size_t as = p, ae = j.find(']', as);
    if (ae == std::string::npos) return result;
    std::string arr = j.substr(as + 1, ae - as - 1);
    std::istringstream iss(arr);
    std::string tok;
    while (std::getline(iss, tok, ',')) {
        tok.erase(std::remove_if(tok.begin(), tok.end(), [](char c){ return isspace(c); }), tok.end());
        if (!tok.empty()) try { result.push_back(std::stoi(tok)); } catch(...) {}
    }
    return result;
}

std::vector<float> TheoryBridge::jFloatArr(const std::string& j, const std::string& k) {
    std::vector<float> result;
    size_t p = findValPos(j, k);
    if (p == std::string::npos || j[p] != '[') return result;
    size_t as = p, ae = j.find(']', as);
    if (ae == std::string::npos) return result;
    std::string arr = j.substr(as + 1, ae - as - 1);
    std::istringstream iss(arr);
    std::string tok;
    while (std::getline(iss, tok, ',')) {
        tok.erase(std::remove_if(tok.begin(), tok.end(), [](char c){ return isspace(c); }), tok.end());
        if (!tok.empty()) try { result.push_back(std::stof(tok)); } catch(...) {}
    }
    return result;
}

std::vector<std::string> TheoryBridge::jStrArr(const std::string& j, const std::string& k) {
    std::vector<std::string> result;
    size_t p = findValPos(j, k);
    if (p == std::string::npos || j[p] != '[') return result;
    size_t as = p, ae = j.find(']', as);
    if (ae == std::string::npos) return result;
    std::string arr = j.substr(as + 1, ae - as - 1);
    size_t pos = 0;
    while (pos < arr.size()) {
        size_t qs = arr.find('"',pos); if(qs==std::string::npos) break;
        size_t qe = arr.find('"',qs+1); if(qe==std::string::npos) break;
        result.push_back(arr.substr(qs+1,qe-qs-1));
        pos = qe+1;
    }
    return result;
}

std::vector<std::string> TheoryBridge::splitJsonArray(const std::string& arr) {
    std::vector<std::string> result;
    int depth = 0; size_t start = std::string::npos;
    for (size_t i=0;i<arr.size();i++) {
        if (arr[i]=='{') {
            if(depth==0) start=i;
            depth++;
        } else if (arr[i]=='}') {
            depth--;
            if(depth==0 && start!=std::string::npos)
                result.push_back(arr.substr(start,i-start+1));
        }
    }
    return result;
}

VoiceLeadingInfo TheoryBridge::parseVoiceLeading(const std::string& j) {
    VoiceLeadingInfo vl;
    vl.distance          = jInt(j,"distance");
    vl.smoothness        = jStr(j,"smoothness");
    vl.description       = jStr(j,"description");
    vl.contrary_motions  = jInt(j,"contrary_motions");
    vl.parallel_motions  = jInt(j,"parallel_motions");
    vl.oblique_motions   = jInt(j,"oblique_motions");
    // parse motions array
    size_t ma = j.find("\"motions\"");
    if (ma != std::string::npos) {
        size_t as = j.find('[',ma), ae = j.find(']',as);
        if (as!=std::string::npos && ae!=std::string::npos) {
            std::string arr = j.substr(as+1,ae-as-1);
            for (auto& obj : splitJsonArray(arr)) {
                VoiceMotion vm;
                vm.from_pc   = jInt(obj,"from_pc");
                vm.to_pc     = jInt(obj,"to_pc");
                vm.semitones = jInt(obj,"semitones");
                vm.from_name = jStr(obj,"from_name");
                vm.to_name   = jStr(obj,"to_name");
                vl.motions.push_back(vm);
            }
        }
    }
    return vl;
}

PLRTransform TheoryBridge::parsePLRTransform(const std::string& j) {
    PLRTransform t;
    t.op             = jStr(j,"op");
    t.op_name        = jStr(j,"op_name");
    t.op_description = jStr(j,"op_description");
    t.vl_cost        = jInt(j,"vl_cost");
    t.common_tone_names = jStrArr(j,"common_tone_names");
    // Parse "to" sub-object
    size_t tp = j.find("\"to\"");
    if (tp != std::string::npos) {
        size_t ts = j.find('{',tp), te = j.find('}',ts);
        if (ts!=std::string::npos && te!=std::string::npos) {
            std::string sub = j.substr(ts,te-ts+1);
            t.to_root    = jInt(sub,"root");
            t.to_quality = jStr(sub,"quality");
            t.to_label   = jStr(sub,"label");
            t.to_pcs     = jIntArr(sub,"pcs");
        }
    }
    // Parse voice_motions
    size_t vm = j.find("\"voice_motions\"");
    if (vm != std::string::npos) {
        size_t as = j.find('[',vm), ae = j.find(']',as);
        if (as!=std::string::npos && ae!=std::string::npos) {
            std::string arr = j.substr(as+1,ae-as-1);
            for (auto& obj : splitJsonArray(arr)) {
                VoiceMotion m;
                m.from_pc   = jInt(obj,"from_pc");
                m.to_pc     = jInt(obj,"to_pc");
                m.semitones = jInt(obj,"semitones");
                m.from_name = jStr(obj,"from_name");
                m.to_name   = jStr(obj,"to_name");
                t.voice_motions.push_back(m);
            }
        }
    }
    return t;
}

ResolutionPath TheoryBridge::parseResolutionPath(const std::string& j) {
    ResolutionPath rp;
    rp.target_pcs     = jIntArr(j,"target_pcs");
    rp.target_label   = jStr(j,"target_label");
    rp.target_root    = jInt(j,"target_root");
    rp.target_quality = jStr(j,"target_quality");
    rp.rule           = jStr(j,"rule");
    rp.rule_class     = jStr(j,"rule_class");
    rp.explanation    = jStr(j,"explanation");
    rp.priority       = jInt(j,"priority");
    // voice_leading sub-object
    size_t vp = j.find("\"voice_leading\"");
    if (vp != std::string::npos) {
        size_t ts = j.find('{',vp), te = j.rfind('}');
        if (ts!=std::string::npos && te!=std::string::npos)
            rp.voice_leading = parseVoiceLeading(j.substr(ts,te-ts+1));
    }
    return rp;
}

CompletionSuggestion TheoryBridge::parseCompletion(const std::string& j) {
    CompletionSuggestion cs;
    cs.pc                = jInt(j,"pc");
    cs.name              = jStr(j,"name");
    cs.score             = jFloat(j,"score");
    cs.roughness_delta   = jFloat(j,"roughness_delta");
    cs.cents_ji          = jFloat(j,"cents_ji");
    cs.in_key            = jBool(j,"in_key");
    cs.function_if_added = jStr(j,"function_if_added");
    cs.tension_if_added  = jFloat(j,"tension_if_added");
    cs.completes_triad   = jBool(j,"completes_triad");
    cs.structural_reasons = jStrArr(j,"structural_reasons");
    return cs;
}

PsychoacousticAnalysis TheoryBridge::parsePsychoacoustic(const std::string& j) {
    PsychoacousticAnalysis pa;

    // Helper for finding nested objects
    auto getSub = [](const std::string& json, const std::string& key) {
        size_t p = json.find("\"" + key + "\"");
        if (p == std::string::npos) return std::string("");
        size_t s = json.find('{', p), e = s;
        if (s == std::string::npos) return std::string("");
        int depth = 1;
        while (depth > 0 && ++e < json.size()) {
            if (json[e] == '{') depth++;
            else if (json[e] == '}') depth--;
        }
        return json.substr(s, e - s + 1);
    };

    std::string l1 = getSub(j, "level_1_peripheral");
    if (!l1.empty()) {
        pa.level1.roughness_normalized = jFloat(l1, "roughness_normalized");
        pa.level1.consonance_score = jFloat(l1, "consonance_score");
        pa.level1.masked_tones = jIntArr(l1, "masked_tones");
        pa.level1.masked_tone_names = jStrArr(l1, "masked_tone_names");
        pa.level1.note_frequencies_hz = jFloatArr(l1, "note_frequencies_hz");
        // Parse roughness_per_pair
        size_t rpa = l1.find("\"roughness_per_pair\"");
        if (rpa != std::string::npos) {
            size_t as = l1.find('[', rpa), ae = l1.find(']', as);
            if (as != std::string::npos && ae != std::string::npos) {
                for (auto& obj : splitJsonArray(l1.substr(as + 1, ae - as - 1))) {
                    RoughnessPair rp;
                    rp.pc_a = jInt(obj, "pc_a"); rp.pc_b = jInt(obj, "pc_b");
                    rp.name_a = jStr(obj, "name_a"); rp.name_b = jStr(obj, "name_b");
                    rp.roughness = jFloat(obj, "roughness");
                    rp.interval_semitones = jInt(obj, "interval_semitones");
                    pa.level1.roughness_per_pair.push_back(rp);
                }
            }
        }
        // Parse cb_overlaps
        size_t cba = l1.find("\"cb_overlaps\"");
        if (cba != std::string::npos) {
            size_t as = l1.find('[', cba), ae = l1.find(']', as);
            if (as != std::string::npos && ae != std::string::npos) {
                for (auto& obj : splitJsonArray(l1.substr(as + 1, ae - as - 1))) {
                    CBOverlap cbo;
                    cbo.pc_a = jInt(obj, "pc_a"); cbo.pc_b = jInt(obj, "pc_b");
                    cbo.interval_semitones = jInt(obj, "interval_semitones");
                    cbo.delta_hz = jFloat(obj, "delta_hz");
                    cbo.cb_overlap_fraction = jFloat(obj, "cb_overlap_fraction");
                    cbo.bm_distance_mm = jFloat(obj, "bm_distance_mm");
                    pa.level1.cb_overlaps.push_back(cbo);
                }
            }
        }
    }

    std::string l2 = getSub(j, "level_2_brainstem");
    if (!l2.empty()) {
        pa.level2.virtual_pitch_hz = jFloat(l2, "virtual_pitch_hz");
        pa.level2.virtual_pitch_name = jStr(l2, "virtual_pitch_name");
        pa.level2.harmonicity = jFloat(l2, "harmonicity");
        pa.level2.harmonicity_label = jStr(l2, "harmonicity_label");
        pa.level2.an_peak_cf_hz = jFloat(l2, "an_peak_cf_hz");
        pa.level2.an_peak_firing_rate = jFloat(l2, "an_peak_firing_rate");
        // vp_top_candidates
        size_t vpa = l2.find("\"vp_top_candidates\"");
        if (vpa != std::string::npos) {
            size_t as = l2.find('[', vpa), ae = l2.find(']', as);
            if (as != std::string::npos && ae != std::string::npos) {
                for (auto& obj : splitJsonArray(l2.substr(as + 1, ae - as - 1))) {
                    VPCandidate vpc;
                    vpc.f0_hz = jFloat(obj, "f0_hz"); vpc.pc = jInt(obj, "pc");
                    vpc.name = jStr(obj, "name"); vpc.score = jFloat(obj, "score");
                    pa.level2.vp_top_candidates.push_back(vpc);
                }
            }
        }
    }

    std::string l3 = getSub(j, "level_3_cortical");
    if (!l3.empty()) {
        pa.level3.kk_tonal_stability = jFloat(l3, "kk_tonal_stability");
        pa.level3.kk_stability_label = jStr(l3, "kk_stability_label");
        pa.level3.kk_raw_score = jFloat(l3, "kk_raw_score");
        pa.level3.spectral_centroid_hz = jFloat(l3, "spectral_centroid_hz");
        pa.level3.spectral_centroid_pc = jInt(l3, "spectral_centroid_pc");
        pa.level3.spectral_centroid_name = jStr(l3, "spectral_centroid_name");
    }

    pa.perceptual_tension = jFloat(j, "perceptual_tension");
    pa.perceptual_tension_label = jStr(j, "perceptual_tension_label");

    std::string tb = getSub(j, "tension_breakdown");
    if (!tb.empty()) {
        pa.tension_breakdown.roughness_component = jFloat(tb, "roughness_component");
        pa.tension_breakdown.harmonicity_component = jFloat(tb, "harmonicity_component");
        pa.tension_breakdown.tonal_component = jFloat(tb, "tonal_component");
    }

    return pa;
}
