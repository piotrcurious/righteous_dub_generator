#pragma once
#include <string>
#include <vector>
#include <functional>
#include <thread>
#include <mutex>
#include <atomic>
#include <deque>
#include <unistd.h>
#include <sys/wait.h>
#include <signal.h>

// ────────────────────────────────────────────────────────────────────────────
//  THEORY BRIDGE — structural edition
//
//  All queries now return algebraic structure, not probabilities.
//  The API mirrors the Python server's command set.
// ────────────────────────────────────────────────────────────────────────────

// ── Interval-class vector
struct IntervalVector {
    int ic[6];   // ic[0]=ic1 count ... ic[5]=ic6 count
    std::string description;  // e.g. "1xM2/m7, 1xm3/M6, 1xM3/m6"
};

// ── Set-class info
struct SetClassInfo {
    std::vector<int> prime_form;
    IntervalVector   icv;
    std::string      forte;         // e.g. "3-11"
    std::string      common_name;   // e.g. "major / minor triad"
    std::string      z_related;     // e.g. "4-Z29" or ""
    bool             inv_symmetric;
    std::vector<int> t_symmetry;    // transpositional symmetry levels
};

// ── Functional analysis result
struct TendencyTone {
    int         pc;
    std::string name;
    std::string role;         // "leading tone (7)", "subdominant (4)"
    std::string tendency;     // "resolves UP 1 semitone to C"
    std::string force;        // "strong", "moderate"
    std::string derivation;   // WHY (algebraic justification)
};

struct FunctionalAnalysis {
    std::string function;         // "TONIC", "DOMINANT", "SUBDOMINANT", etc.
    std::string function_reason;  // Full algebraic explanation
    int         tension_level;    // 0-4
    bool        has_tritone;
    std::vector<int> tritone;     // [lt_pc, sd_pc]
    std::vector<std::string> tritone_names;
    bool        in_diatonic;
    std::vector<int> chromatic_tones;
    std::vector<TendencyTone> tendency_tones;
};

// ── Tonnetz tension
struct TonnetzTension {
    float       tension;          // 0=rest, 1=maximal
    std::string tension_label;    // "at rest", "mild", "tension", etc.
    float       distance;
    float       chord_centroid[2];
    float       tonic_centroid[2];
};

// ── PLR transformation
struct VoiceMotion {
    int from_pc, to_pc, semitones;
    std::string from_name, to_name;
};

struct PLRTransform {
    std::string op, op_name, op_description;
    int         to_root;
    std::string to_quality, to_label;
    std::vector<int>         to_pcs;
    std::vector<int>         common_tones;
    std::vector<std::string> common_tone_names;
    std::vector<VoiceMotion> voice_motions;
    int                      vl_cost;
};

// ── Resolution path
struct VoiceLeadingInfo {
    int         distance;
    std::string smoothness;
    std::string description;
    std::vector<VoiceMotion> motions;
    int         contrary_motions, parallel_motions, oblique_motions;
};

struct ResolutionPath {
    std::vector<int> target_pcs;
    std::string      target_label;
    int              target_root;
    std::string      target_quality;
    std::string      rule;         // "TRITONE_RESOLUTION", "PLR_P", etc.
    std::string      rule_class;   // "necessity", "minimal_motion", "scalar"
    std::string      explanation;  // Full algebraic justification
    VoiceLeadingInfo voice_leading;
    int              priority;     // 1=necessity, 2=cadence, 3=PLR, 4=scalar
};

// ── Completion suggestion
struct CompletionSuggestion {
    int         pc;
    std::string name;
    float       score;
    float       roughness_delta;
    float       cents_ji;
    bool        in_key;
    std::string function_if_added;
    float       tension_if_added;
    std::vector<std::string> structural_reasons;  // WHY to add this note
    bool        completes_triad;
};

// ── Sequence detection
struct SequenceInfo {
    std::string type;           // "transposition", "PLR_P_sequence", "circle_of_fifths"
    int         interval;       // for transposition sequences
    std::string interval_name;
    std::string description;
    std::string algebraic_explanation;
    int         period;
};

// ─── Callbacks
using SetClassCb    = std::function<void(const SetClassInfo&)>;
using FunctionCb    = std::function<void(const FunctionalAnalysis&, const TonnetzTension&)>;
using PLRNeighborsCb= std::function<void(const std::vector<PLRTransform>&)>;
using PLRPathCb     = std::function<void(const std::vector<std::string>&, bool reachable)>;
using ResolutionsCb = std::function<void(const std::vector<ResolutionPath>&)>;
using CompletionCb  = std::function<void(const std::vector<CompletionSuggestion>&)>;
using OrbifoldCb    = std::function<void(const VoiceLeadingInfo&)>;
using SequenceCb    = std::function<void(const SequenceInfo&)>;
using RawJsonCb     = std::function<void(const std::string&)>;

// ────────────────────────────────────────────────────────────────────────────
class TheoryBridge {
public:
    TheoryBridge();
    ~TheoryBridge();

    bool start(const std::string& theory_dir);
    void stop();
    bool isRunning() const { return running_.load(); }

    // ── Structural queries (all async, result via callback)

    // Full structural analysis of current chord
    void queryAnalyzeChord(const std::vector<int>& pcs, int key,
                           FunctionCb cb);

    // PLR neighborhood (3 direct neighbors with justification)
    void queryPLRNeighbors(int root, const std::string& quality,
                            PLRNeighborsCb cb);

    // Shortest PLR path between two triads
    void queryPLRPath(int root_a, const std::string& qual_a,
                      int root_b, const std::string& qual_b,
                      PLRPathCb cb);

    // Resolution paths (algebraically justified, not probabilistic)
    void queryResolutionPaths(int root, const std::string& quality, int key,
                               ResolutionsCb cb);

    // Completion suggestions with structural reasons
    void querySuggestCompletion(const std::vector<int>& pcs, int key,
                                 CompletionCb cb);

    // Orbifold voice-leading
    void queryOrbifoldDistance(const std::vector<int>& a,
                                const std::vector<int>& b,
                                OrbifoldCb cb);

    // Sequence detection
    void queryDetectSequence(const std::vector<std::pair<int,std::string>>& chords,
                              SequenceCb cb);

    // EDO analysis
    void queryEDOAnalysis(int edo, RawJsonCb cb);

    // Raw JSON passthrough
    void queryRaw(const std::string& json, const std::string& result_tag, RawJsonCb cb);

    // ── Call from UI idle handler
    void poll();

    std::string sendSync(const std::string& json, int timeout_ms = 150);

    // ── JSON helpers (public for UI use in main.cpp)
    static int         jInt  (const std::string& j, const std::string& k, int d=0);
    static float       jFloat(const std::string& j, const std::string& k, float d=0.f);
    static std::string jStr  (const std::string& j, const std::string& k,
                               const std::string& d="");
    static bool        jBool (const std::string& j, const std::string& k, bool d=false);

private:
    pid_t py_pid_{-1};
    int   to_py_[2]{-1,-1};
    int   from_py_[2]{-1,-1};

    std::thread      reader_thread_;
    std::atomic<bool> running_{false};

    struct PendingCb {
        std::string result_tag;
        std::function<void(const std::string&)> raw_cb;
    };
    mutable std::mutex pending_mutex_;
    std::deque<PendingCb> pending_cbs_;

    mutable std::mutex response_mutex_;
    std::deque<std::string> responses_;

    void readerLoop();
    void sendRaw(const std::string& json);

    // Parse arrays
    static std::vector<int> jIntArr(const std::string& j, const std::string& k);
    static std::vector<std::string> jStrArr(const std::string& j, const std::string& k);

    // Parse structured objects from JSON
    static VoiceLeadingInfo parseVoiceLeading(const std::string& j);
    static PLRTransform     parsePLRTransform(const std::string& obj_json);
    static ResolutionPath   parseResolutionPath(const std::string& obj_json);
    static CompletionSuggestion parseCompletion(const std::string& obj_json);

    // Split a JSON array string into object strings
    static std::vector<std::string> splitJsonArray(const std::string& arr);
};
