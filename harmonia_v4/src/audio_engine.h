#pragma once
#include "voice.h"
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <functional>
#include <deque>
#include <alsa/asoundlib.h>

// ────────────────────────────────────────────────────────────────────────────
//  AUDIO ENGINE — ALSA-backed additive synthesis across N voices.
//
//  Architecture:
//    UI thread   → modifies voices (mutex-protected)
//    audio thread→ reads voices, synthesizes samples, writes to ALSA PCM
//
//  Signal path per frame:
//    sample[t] = Σ_voices  Σₖ voice.env · voice.amp[k] · sin(phase[k])
//    phase[k] += 2π·k·f / Fs
//
//  Roughness is computed every ROUGHNESS_UPDATE_FRAMES frames and
//  stored for the UI to read without blocking the audio thread.
// ────────────────────────────────────────────────────────────────────────────

constexpr int    AUDIO_CHANNELS      = 2;
constexpr int    AUDIO_PERIOD_FRAMES = 512;
constexpr int    MAX_VOICES          = 16;
constexpr int    ROUGHNESS_UPDATE_FRAMES = 4096;
constexpr int    SPECTRUM_SIZE       = 4096;   // for display FFT

// ── Pairwise roughness record
struct RoughnessRecord {
    int voice_a, voice_b;
    float roughness;          // 0=consonant, 1=maximally rough
    float virtual_pitch_hz;   // shared virtual pitch
};

// ── Abstract object: a perceptual gestalt across voices
struct AbstractObject {
    std::string chord_name;       // e.g. "Cmaj7", "G7", "Am"
    std::string quality;          // e.g. "maj", "min"
    int         root_pc;          // 0-11
    float       confidence;       // 0-1
    float       roughness_total;  // aggregate
    float       virtual_pitch_hz; // lowest implied fundamental
    std::vector<int> voice_ids;   // contributing voices
    int         tonnetz_x, tonnetz_y; // centroid in Tonnetz
};

// ── Spectrum snapshot (for display)
struct SpectrumSnapshot {
    std::array<float, SPECTRUM_SIZE/2> magnitude;
    double time;
};

class AudioEngine {
public:
    AudioEngine();
    ~AudioEngine();

    bool init(const char* device = "default");
    void shutdown();

    // ── voice management (call from UI thread, mutex-protected)
    int  addVoice(int midi_note, TimbrePreset timbre = TimbrePreset::SINE);
    void removeVoice(int voice_id);
    void noteOn(int voice_id);
    void noteOff(int voice_id);
    void setVoiceFrequency(int voice_id, double hz);
    void setVoiceAmplitude(int voice_id, float amp);
    void setVoiceTimbre(int voice_id, TimbrePreset t);
    void setVoiceDetune(int voice_id, float cents);

    // ── read-only access (UI thread)
    std::vector<Voice> getVoiceSnapshot() const;
    std::vector<RoughnessRecord> getRoughnessSnapshot() const;
    AbstractObject getAbstractObject() const;
    SpectrumSnapshot getSpectrumSnapshot() const;

    // ── settings
    void setEDO(int edo);
    int  getEDO() const { return edo_.load(); }

    // ── master
    void setMasterVolume(float v) { master_volume_.store(v); }
    float getMasterVolume() const  { return master_volume_.load(); }

    // ── callbacks
    using RoughnessCallback = std::function<void(const std::vector<RoughnessRecord>&, const AbstractObject&)>;
    void setRoughnessCallback(RoughnessCallback cb) { roughness_cb_ = std::move(cb); }

    bool isRunning() const { return running_.load(); }

private:
    // ── ALSA state
    snd_pcm_t*   pcm_{nullptr};
    int16_t*     alsa_buf_{nullptr};

    // ── voices (audio thread owns write access; UI reads snapshot)
    mutable std::mutex voices_mutex_;
    std::vector<Voice> voices_;
    int next_voice_id_{0};

    // ── synthesis state
    float  sample_buf_[AUDIO_PERIOD_FRAMES * AUDIO_CHANNELS];

    // ── audio thread
    std::thread audio_thread_;
    std::atomic<bool> running_{false};
    std::atomic<float> master_volume_{0.5f};
    std::atomic<int>   edo_{12};

    // ── shared results (written by audio thread, read by UI)
    mutable std::mutex results_mutex_;
    std::vector<RoughnessRecord> roughness_records_;
    AbstractObject abstract_object_;
    SpectrumSnapshot spectrum_snap_;

    // ── internal
    int roughness_counter_{0};
    RoughnessCallback roughness_cb_;

    void audioLoop();
    void synthesizeFrame(float* out, int frames);
    void updateEnvelope(Voice& v, int frames);
    void computeRoughness();
    void computeVirtualPitch(AbstractObject& obj);
    void identifyChord(AbstractObject& obj);
    float pairRoughness(const Voice& a, const Voice& b);
    float criticalBandwidth(float hz);
    void updateSpectrum(const float* buf, int frames);
};
