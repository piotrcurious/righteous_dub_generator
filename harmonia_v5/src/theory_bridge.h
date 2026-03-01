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
//  THEORY BRIDGE — Psychoacoustic Edition
//
//  Integrates algebraic geometry with psychoacoustic neural models.
// ────────────────────────────────────────────────────────────────────────────

// ── Psychoacoustic data structures
struct RoughnessPair {
    int pc_a, pc_b;
    std::string name_a, name_b;
    float roughness;
    int interval_semitones;
};

struct CBOverlap {
    int pc_a, pc_b;
    int interval_semitones;
    float delta_hz;
    float cb_overlap_fraction;
    float bm_distance_mm;
};

struct ToneAudibility {
    int pc;
    std::string name;
    float freq_hz;
    float masking_threshold_db;
    float sensation_level_db;
    bool audible;
};

struct VPCandidate {
    float f0_hz;
    int pc;
    std::string name;
    float score;
};

struct PsychoacousticAnalysis {
    struct {
        float roughness_normalized;
        float consonance_score;
        std::vector<RoughnessPair> roughness_per_pair;
        std::vector<int> masked_tones;
        std::vector<std::string> masked_tone_names;
        std::vector<CBOverlap> cb_overlaps;
        std::vector<float> note_frequencies_hz;
    } level1;

    struct {
        float virtual_pitch_hz;
        int   virtual_pitch_pc;
        std::string virtual_pitch_name;
        float harmonicity;
        std::string harmonicity_label;
        std::vector<VPCandidate> vp_top_candidates;
        float an_peak_cf_hz;
        float an_peak_firing_rate;
    } level2;

    struct {
        float kk_tonal_stability;
        std::string kk_stability_label;
        float kk_raw_score;
        float spectral_centroid_hz;
        int spectral_centroid_pc;
        std::string spectral_centroid_name;
    } level3;

    float perceptual_tension;
    std::string perceptual_tension_label;
    struct {
        float roughness_component;
        float harmonicity_component;
        float tonal_component;
    } tension_breakdown;
};

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

// ── Modulation & Pivot Search
struct PivotChord {
    std::string type;
    std::vector<int> pcs;
    std::string label;
    int root;
    std::string quality;
    std::string roman_from, roman_to;
    std::string function_from, function_to;
    float pivot_score;
    std::string interpretation;
};

struct ModulationStep {
    int root;
    std::string quality, label, roman, key_context, function;
    int vl_from_prev, cumulative_cost;
    // For pivot steps
    std::string pivot_note;
    std::string roman_as_pivot_from, roman_as_pivot_to;
};

struct PivotSearchResult {
    int key_from, key_to;
    std::string key_from_name, key_to_name;
    int cof_distance;
    std::string cof_relationship;
    std::vector<int> common_scale_tones;
    std::vector<PivotChord> all_pivots;
    std::vector<ModulationStep> modulation_path;
};

// ─── Callbacks
using SetClassCb    = std::function<void(const SetClassInfo&)>;
using FunctionCb    = std::function<void(const FunctionalAnalysis&, const TonnetzTension&)>;
using PsychoacousticCb = std::function<void(const PsychoacousticAnalysis&)>;
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

    // ── Queries (all async, result via callback)

    // Full structural analysis of current chord
    void queryAnalyzeChord(const std::vector<int>& pcs, int key, int edo,
                           FunctionCb cb);

    // Full psychoacoustic neural model analysis
    void queryPsychoacoustic(const std::vector<int>& pcs, int key, int edo,
                             float c4_hz, int octave, float level_db,
                             PsychoacousticCb cb);

    // PLR neighborhood (3 direct neighbors with justification)
    void queryPLRNeighbors(int root, const std::string& quality,
                            PLRNeighborsCb cb);

    // Shortest PLR path between two triads
    void queryPLRPath(int root_a, const std::string& qual_a,
                      int root_b, const std::string& qual_b,
                      PLRPathCb cb);

    // Resolution paths (algebraically justified, not probabilistic)
    void queryResolutionPaths(int root, const std::string& quality, int key, int edo,
                               ResolutionsCb cb);

    // Completion suggestions with structural reasons
    void querySuggestCompletion(const std::vector<int>& pcs, int key, int edo,
                                 CompletionCb cb);

    // Orbifold voice-leading
    void queryOrbifoldDistance(const std::vector<int>& a,
                                const std::vector<int>& b,
                                int edo,
                                OrbifoldCb cb);

    // Sequence detection
    void queryDetectSequence(const std::vector<std::pair<int,std::string>>& chords,
                              int edo,
                              SequenceCb cb);

    // Pivot search for modulation
    using PivotSearchCb = std::function<void(const PivotSearchResult&)>;
    void queryPivotSearch(int key_from, int key_to, int edo, PivotSearchCb cb);

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

    std::atomic<int> request_counter_{0};

    void readerLoop();
    void sendRaw(const std::string& json);

    // Parse arrays
    static std::vector<int> jIntArr(const std::string& j, const std::string& k);
    static std::vector<float> jFloatArr(const std::string& j, const std::string& k);
    static std::vector<std::string> jStrArr(const std::string& j, const std::string& k);

    // Parse structured objects from JSON
    static VoiceLeadingInfo parseVoiceLeading(const std::string& j);
    static PLRTransform     parsePLRTransform(const std::string& obj_json);
    static ResolutionPath   parseResolutionPath(const std::string& obj_json);
    static CompletionSuggestion parseCompletion(const std::string& obj_json);
    static PsychoacousticAnalysis parsePsychoacoustic(const std::string& j);

    // Split a JSON array string into object strings
    static std::vector<std::string> splitJsonArray(const std::string& arr);
};
