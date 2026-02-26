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
        std::string tag = jStr(resp,"result");
        std::lock_guard<std::mutex> lk(pending_mutex_);
        for (auto it = pending_cbs_.begin(); it != pending_cbs_.end(); ++it) {
            if (it->result_tag == tag) {
                auto cb = it->raw_cb; pending_cbs_.erase(it);
                if (cb) { cb(resp); }
                break;
            }
        }
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
void TheoryBridge::queryAnalyzeChord(const std::vector<int>& pcs, int key,
                                      FunctionCb cb) {
    std::string arr = "[";
    for (size_t i=0;i<pcs.size();i++) { arr+=std::to_string(pcs[i]); if(i+1<pcs.size()) arr+=","; }
    arr += "]";
    std::string json = "{\"cmd\":\"analyze_chord\",\"pcs\":" + arr +
                       ",\"key\":" + std::to_string(key) + "}";
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
    std::string json = "{\"cmd\":\"plr_neighbors\",\"root\":" + std::to_string(root) +
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
    std::string json = "{\"cmd\":\"plr_path\",\"root_a\":" + std::to_string(ra) +
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

void TheoryBridge::queryResolutionPaths(int root, const std::string& quality, int key,
                                         ResolutionsCb cb) {
    std::string json = "{\"cmd\":\"resolution_paths\",\"root\":" + std::to_string(root) +
                       ",\"quality\":\"" + quality + "\",\"key\":" + std::to_string(key) + "}";
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

void TheoryBridge::querySuggestCompletion(const std::vector<int>& pcs, int key,
                                           CompletionCb cb) {
    std::string arr = "[";
    for (size_t i=0;i<pcs.size();i++) { arr+=std::to_string(pcs[i]); if(i+1<pcs.size()) arr+=","; }
    arr += "]";
    std::string json = "{\"cmd\":\"suggest_completion\",\"pitch_classes\":" + arr +
                       ",\"key\":" + std::to_string(key) + "}";
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
    std::string json = "{\"cmd\":\"orbifold_distance\",\"chord_a\":" + mkArr(a) +
                       ",\"chord_b\":" + mkArr(b) + "}";
    auto raw_cb = [cb](const std::string& resp) {
        if (cb) cb(parseVoiceLeading(resp));
    };
    { std::lock_guard<std::mutex> lk(pending_mutex_);
      pending_cbs_.push_back({"orbifold_distance", raw_cb}); }
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
    std::string json = "{\"cmd\":\"detect_sequence\",\"chords\":" + arr + "}";
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
    std::string json = "{\"cmd\":\"edo_analysis\",\"edo\":" + std::to_string(edo) + "}";
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
int TheoryBridge::jInt(const std::string& j, const std::string& k, int d) {
    std::string needle = "\"" + k + "\"";
    size_t p = j.find(needle); if(p==std::string::npos) return d;
    p = j.find_first_of("0123456789-", p+needle.size()+1);
    if(p==std::string::npos) return d;
    return std::stoi(j.substr(p));
}

float TheoryBridge::jFloat(const std::string& j, const std::string& k, float d) {
    std::string needle = "\"" + k + "\"";
    size_t p = j.find(needle); if(p==std::string::npos) return d;
    p = j.find_first_of("0123456789-.", p+needle.size()+1);
    if(p==std::string::npos) return d;
    try { return std::stof(j.substr(p)); } catch(...) { return d; }
}

std::string TheoryBridge::jStr(const std::string& j, const std::string& k,
                                const std::string& d) {
    std::string needle = "\"" + k + "\":\"";
    size_t p = j.find(needle); if(p==std::string::npos) return d;
    p += needle.size();
    size_t e = p;
    while (e < j.size() && !(j[e]=='"' && (e==0||j[e-1]!='\\'))) e++;
    return (e<j.size()) ? j.substr(p,e-p) : d;
}

bool TheoryBridge::jBool(const std::string& j, const std::string& k, bool d) {
    std::string needle = "\"" + k + "\":";
    size_t p = j.find(needle); if(p==std::string::npos) return d;
    p += needle.size();
    while (p<j.size() && j[p]==' ') p++;
    if (p<j.size()) {
        if (j[p]=='t') return true;
        if (j[p]=='f') return false;
    }
    return d;
}

std::vector<int> TheoryBridge::jIntArr(const std::string& j, const std::string& k) {
    std::vector<int> result;
    std::string needle = "\"" + k + "\"";
    size_t p = j.find(needle); if(p==std::string::npos) return result;
    size_t as = j.find('[',p), ae = j.find(']',as);
    if(as==std::string::npos||ae==std::string::npos) return result;
    std::string arr = j.substr(as+1,ae-as-1);
    std::istringstream iss(arr);
    std::string tok;
    while (std::getline(iss,tok,',')) {
        tok.erase(std::remove_if(tok.begin(),tok.end(),[](char c){return c==' ';}),tok.end());
        if (!tok.empty()) try { result.push_back(std::stoi(tok)); } catch(...) {}
    }
    return result;
}

std::vector<std::string> TheoryBridge::jStrArr(const std::string& j, const std::string& k) {
    std::vector<std::string> result;
    std::string needle = "\"" + k + "\"";
    size_t p = j.find(needle); if(p==std::string::npos) return result;
    size_t as = j.find('[',p), ae = j.find(']',as);
    if(as==std::string::npos||ae==std::string::npos) return result;
    std::string arr = j.substr(as+1,ae-as-1);
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
