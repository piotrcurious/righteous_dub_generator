#include "audio_engine.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <map>
#include <sstream>
#include <stdio.h>

// ────────────────────────────────────────────────────────────────────────────
//  Chord template matching.
//  Each entry: { chord_name, { intervals in semitones from root } }
// ────────────────────────────────────────────────────────────────────────────
struct ChordTemplate {
    const char* name;
    std::vector<int> intervals;  // semitones from root (not including 0)
};

static const ChordTemplate CHORD_TEMPLATES[] = {
    {"maj",   {4, 7}},
    {"min",   {3, 7}},
    {"dim",   {3, 6}},
    {"aug",   {4, 8}},
    {"sus2",  {2, 7}},
    {"sus4",  {5, 7}},
    {"7",     {4, 7, 10}},
    {"maj7",  {4, 7, 11}},
    {"min7",  {3, 7, 10}},
    {"dim7",  {3, 6, 9}},
    {"m7b5",  {3, 6, 10}},
    {"9",     {4, 7, 10, 14}},
    {"add9",  {4, 7, 14}},
    {"6",     {4, 7, 9}},
    {"m6",    {3, 7, 9}},
};

// ────────────────────────────────────────────────────────────────────────────
AudioEngine::AudioEngine() {}

AudioEngine::~AudioEngine() { shutdown(); }

bool AudioEngine::init(const char* device) {
    int err;
    snd_pcm_hw_params_t* hw_params = nullptr;

    if ((err = snd_pcm_open(&pcm_, device, SND_PCM_STREAM_PLAYBACK, 0)) < 0) {
        fprintf(stderr, "ALSA open error: %s\n", snd_strerror(err));
        return false;
    }

    snd_pcm_hw_params_alloca(&hw_params);
    snd_pcm_hw_params_any(pcm_, hw_params);
    snd_pcm_hw_params_set_access(pcm_, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED);
    snd_pcm_hw_params_set_format(pcm_, hw_params, SND_PCM_FORMAT_S16_LE);

    unsigned int rate = static_cast<unsigned int>(SAMPLE_RATE);
    snd_pcm_hw_params_set_rate_near(pcm_, hw_params, &rate, nullptr);
    snd_pcm_hw_params_set_channels(pcm_, hw_params, AUDIO_CHANNELS);

    snd_pcm_uframes_t period = AUDIO_PERIOD_FRAMES;
    snd_pcm_hw_params_set_period_size_near(pcm_, hw_params, &period, nullptr);

    snd_pcm_uframes_t buffer_size = period * 4;
    snd_pcm_hw_params_set_buffer_size_near(pcm_, hw_params, &buffer_size);

    if ((err = snd_pcm_hw_params(pcm_, hw_params)) < 0) {
        fprintf(stderr, "ALSA hw_params error: %s\n", snd_strerror(err));
        snd_pcm_close(pcm_); pcm_ = nullptr;
        return false;
    }

    snd_pcm_prepare(pcm_);

    alsa_buf_ = new int16_t[AUDIO_PERIOD_FRAMES * AUDIO_CHANNELS];

    running_.store(true);
    audio_thread_ = std::thread(&AudioEngine::audioLoop, this);
    fprintf(stderr, "AudioEngine: started at %.0f Hz, period %d\n",
            SAMPLE_RATE, AUDIO_PERIOD_FRAMES);
    return true;
}

void AudioEngine::shutdown() {
    running_.store(false);
    if (audio_thread_.joinable()) audio_thread_.join();
    if (pcm_) { snd_pcm_drain(pcm_); snd_pcm_close(pcm_); pcm_ = nullptr; }
    delete[] alsa_buf_; alsa_buf_ = nullptr;
}

// ────────────────────────────────────────────────────────────────────────────
//  Voice management
// ────────────────────────────────────────────────────────────────────────────
int AudioEngine::addVoice(int midi_note, TimbrePreset timbre) {
    std::lock_guard<std::mutex> lk(voices_mutex_);
    if ((int)voices_.size() >= MAX_VOICES) return -1;
    Voice v;
    v.id = next_voice_id_++;
    v.edo = edo_.load();
    v.name = std::string(NOTE_NAMES[midi_note % 12]) + std::to_string(midi_note/12-1);
    v.setTimbre(timbre);
    v.setMidiNote(midi_note);
    v.active.store(false);
    v.env_stage = Voice::EnvStage::IDLE;
    v.env_value = 0.f;
    // assign a colour by pitch class
    static const float PC_COLORS[12][3] = {
        {1,.3,.3},{1,.5,.2},{1,.8,.2},{.7,1,.2},{.3,1,.3},{.2,.9,.5},
        {.2,.8,.9},{.3,.5,1},{.5,.3,1},{.8,.3,1},{1,.2,.8},{1,.2,.5}
    };
    int pc = midi_note % 12;
    v.color[0] = PC_COLORS[pc][0];
    v.color[1] = PC_COLORS[pc][1];
    v.color[2] = PC_COLORS[pc][2];
    voices_.push_back(std::move(v));
    return voices_.back().id;
}

void AudioEngine::removeVoice(int voice_id) {
    std::lock_guard<std::mutex> lk(voices_mutex_);
    voices_.erase(std::remove_if(voices_.begin(), voices_.end(),
        [voice_id](const Voice& v){ return v.id == voice_id; }), voices_.end());
}

void AudioEngine::noteOn(int voice_id) {
    std::lock_guard<std::mutex> lk(voices_mutex_);
    for (auto& v : voices_) if (v.id == voice_id) {
        v.note_on.store(true);
        v.active.store(true);
        v.env_stage = Voice::EnvStage::ATTACK;
        break;
    }
}

void AudioEngine::noteOff(int voice_id) {
    std::lock_guard<std::mutex> lk(voices_mutex_);
    for (auto& v : voices_) if (v.id == voice_id) {
        v.note_on.store(false);
        if (v.env_stage != Voice::EnvStage::IDLE)
            v.env_stage = Voice::EnvStage::RELEASE;
        break;
    }
}

void AudioEngine::setVoiceFrequency(int voice_id, double hz) {
    std::lock_guard<std::mutex> lk(voices_mutex_);
    for (auto& v : voices_) if (v.id == voice_id) { v.setFrequency(hz); break; }
}

void AudioEngine::setVoiceAmplitude(int voice_id, float amp) {
    std::lock_guard<std::mutex> lk(voices_mutex_);
    for (auto& v : voices_) if (v.id == voice_id) { v.amplitude = amp; break; }
}

void AudioEngine::setVoiceTimbre(int voice_id, TimbrePreset t) {
    std::lock_guard<std::mutex> lk(voices_mutex_);
    for (auto& v : voices_) if (v.id == voice_id) { v.setTimbre(t); break; }
}

void AudioEngine::setVoiceDetune(int voice_id, float cents) {
    std::lock_guard<std::mutex> lk(voices_mutex_);
    for (auto& v : voices_) if (v.id == voice_id) {
        v.detune_cents = cents;
        v.edo = edo_.load();
        v.setFrequency(v.frequency);
        break;
    }
}

std::vector<Voice> AudioEngine::getVoiceSnapshot() const {
    std::lock_guard<std::mutex> lk(voices_mutex_);
    return voices_;
}

std::vector<RoughnessRecord> AudioEngine::getRoughnessSnapshot() const {
    std::lock_guard<std::mutex> lk(results_mutex_);
    return roughness_records_;
}

AbstractObject AudioEngine::getAbstractObject() const {
    std::lock_guard<std::mutex> lk(results_mutex_);
    return abstract_object_;
}

SpectrumSnapshot AudioEngine::getSpectrumSnapshot() const {
    std::lock_guard<std::mutex> lk(results_mutex_);
    return spectrum_snap_;
}

// ────────────────────────────────────────────────────────────────────────────
//  Audio thread
// ────────────────────────────────────────────────────────────────────────────
void AudioEngine::audioLoop() {
    while (running_.load()) {
        synthesizeFrame(sample_buf_, AUDIO_PERIOD_FRAMES);

        // Convert float to S16LE
        float mv = master_volume_.load();
        for (int i = 0; i < AUDIO_PERIOD_FRAMES; i++) {
            float s = std::max(-1.f, std::min(1.f, sample_buf_[i] * mv));
            int16_t s16 = static_cast<int16_t>(s * 32767.f);
            alsa_buf_[i * AUDIO_CHANNELS + 0] = s16;
            alsa_buf_[i * AUDIO_CHANNELS + 1] = s16;
        }

        int err = snd_pcm_writei(pcm_, alsa_buf_, AUDIO_PERIOD_FRAMES);
        if (err == -EPIPE) {
            snd_pcm_prepare(pcm_);
        } else if (err < 0) {
            fprintf(stderr, "ALSA write error: %s\n", snd_strerror(err));
        }

        roughness_counter_ += AUDIO_PERIOD_FRAMES;
        if (roughness_counter_ >= ROUGHNESS_UPDATE_FRAMES) {
            roughness_counter_ = 0;
            computeRoughness();
        }
        updateSpectrum(sample_buf_, AUDIO_PERIOD_FRAMES);
    }
}

// ────────────────────────────────────────────────────────────────────────────
//  Additive synthesis + ADSR
// ────────────────────────────────────────────────────────────────────────────
void AudioEngine::updateEnvelope(Voice& v, int frames) {
    float dt = (float)frames / SAMPLE_RATE;
    switch (v.env_stage) {
    case Voice::EnvStage::ATTACK:
        v.env_value += dt / (v.attack_ms * 0.001f);
        if (v.env_value >= 1.f) { v.env_value = 1.f; v.env_stage = Voice::EnvStage::DECAY; }
        break;
    case Voice::EnvStage::DECAY:
        v.env_value -= dt / (v.decay_ms * 0.001f) * (1.f - v.sustain_level);
        if (v.env_value <= v.sustain_level) { v.env_value = v.sustain_level; v.env_stage = Voice::EnvStage::SUSTAIN; }
        break;
    case Voice::EnvStage::SUSTAIN:
        v.env_value = v.sustain_level;
        break;
    case Voice::EnvStage::RELEASE:
        v.env_value -= dt / (v.release_ms * 0.001f) * v.sustain_level;
        if (v.env_value <= 0.f) { v.env_value = 0.f; v.env_stage = Voice::EnvStage::IDLE; v.active.store(false); }
        break;
    default: break;
    }
}

void AudioEngine::synthesizeFrame(float* out, int frames) {
    memset(out, 0, sizeof(float) * frames);
    std::lock_guard<std::mutex> lk(voices_mutex_);

    for (auto& v : voices_) {
        if (!v.active.load()) continue;
        updateEnvelope(v, frames);
        float gain = v.amplitude * v.env_value;
        if (gain < 1e-6f) continue;

        for (int i = 0; i < frames; i++) {
            float sample = 0.f;
            for (int k = 0; k < MAX_HARMONICS; k++) {
                if (v.harmonic_amp[k] < 1e-6f) continue;
                sample += v.harmonic_amp[k] * std::sin(v.harmonic_phase[k]);
                v.harmonic_phase[k] += v.harmonic_phase_inc[k];
                // Keep phase in [-π, π]
                if (v.harmonic_phase[k] > M_PI) v.harmonic_phase[k] -= 2.f * M_PI;
            }
            out[i] += gain * sample;
        }
    }

    // Soft clipping
    for (int i = 0; i < frames; i++) {
        float x = out[i] * 0.3f;  // scale down for multiple voices
        out[i] = x / (1.f + std::abs(x));
    }
}

// ────────────────────────────────────────────────────────────────────────────
//  Psychoacoustics: Sethares roughness model (1993)
//
//  For each pair of partials (fᵢ, aᵢ), (fⱼ, aⱼ):
//    R_pair = aᵢ·aⱼ · s(|fᵢ-fⱼ|, mean(fᵢ,fⱼ))
//
//  Where s() is the Plomp-Levelt curve scaled by critical bandwidth.
// ────────────────────────────────────────────────────────────────────────────
float AudioEngine::criticalBandwidth(float hz) {
    // Moore & Glasberg (1990) ERB approximation
    return 24.7f * (4.37f * hz / 1000.f + 1.f);
}

float AudioEngine::pairRoughness(const Voice& a, const Voice& b) {
    if (!a.active.load() || !b.active.load()) return 0.f;
    float total = 0.f;
    // iterate over all partial pairs
    for (int ka = 1; ka <= MAX_HARMONICS; ka++) {
        if (a.harmonic_amp[ka-1] < 1e-6f) continue;
        float fa = static_cast<float>(a.frequency) * ka;
        float aa = a.harmonic_amp[ka-1] * a.amplitude * a.env_value;

        for (int kb = 1; kb <= MAX_HARMONICS; kb++) {
            if (b.harmonic_amp[kb-1] < 1e-6f) continue;
            float fb = static_cast<float>(b.frequency) * kb;
            float ab = b.harmonic_amp[kb-1] * b.amplitude * b.env_value;

            float fmean  = 0.5f * (fa + fb);
            float cbw    = criticalBandwidth(fmean);
            float df     = std::abs(fa - fb);
            float x      = df / cbw;
            // Plomp-Levelt: peak at x≈0.25, zero at x=0 and x≥1
            float roughness = (x < 1.f)
                ? aa * ab * x * std::exp(1.f - x)
                : 0.f;
            total += roughness;
        }
    }
    return std::min(1.f, total);
}

// ────────────────────────────────────────────────────────────────────────────
//  Virtual pitch via autocorrelation of summed spectrum
//
//  Build a "neural activity pattern" (NAP): sum all active harmonic amplitudes
//  at their frequencies. Compute autocorrelation → peak = virtual fundamental.
// ────────────────────────────────────────────────────────────────────────────
void AudioEngine::computeVirtualPitch(AbstractObject& obj) {
    // Build amplitude spectrum on a log-frequency grid
    constexpr int N = 1024;
    constexpr float F_MIN = 30.f, F_MAX = 4000.f;
    static float spectrum[N] = {};
    memset(spectrum, 0, sizeof(spectrum));

    float log_min = std::log2(F_MIN), log_max = std::log2(F_MAX);

    for (auto& v : voices_) {
        if (!v.active.load()) continue;
        for (int k = 1; k <= MAX_HARMONICS; k++) {
            if (v.harmonic_amp[k-1] < 1e-6f) continue;
            float f = static_cast<float>(v.frequency) * k;
            if (f < F_MIN || f > F_MAX) continue;
            float logf = std::log2(f);
            int idx = static_cast<int>((logf - log_min) / (log_max - log_min) * N);
            idx = std::max(0, std::min(N-1, idx));
            spectrum[idx] += v.harmonic_amp[k-1] * v.amplitude * v.env_value;
        }
    }

    // Autocorrelation in log-frequency space
    // A lag of 'tau' corresponds to multiplying frequency by 2^(tau*(log_max-log_min)/N)
    float best_corr = -1.f;
    int best_lag = 0;
    for (int lag = 5; lag < N/4; lag++) {
        float corr = 0.f;
        for (int i = 0; i + lag < N; i++) corr += spectrum[i] * spectrum[i + lag];
        if (corr > best_corr) { best_corr = corr; best_lag = lag; }
    }

    // Convert lag to Hz
    float log_period = best_lag * (log_max - log_min) / N;
    obj.virtual_pitch_hz = F_MIN * std::pow(2.f, log_period);
    // clamp to reasonable range and quantize to nearest semitone
    obj.virtual_pitch_hz = std::max(30.f, std::min(1000.f, obj.virtual_pitch_hz));
}

// ────────────────────────────────────────────────────────────────────────────
//  Chord identification by pitch-class set matching
// ────────────────────────────────────────────────────────────────────────────
void AudioEngine::identifyChord(AbstractObject& obj) {
    // Collect active pitch classes
    std::vector<int> pcs;
    for (auto& v : voices_) {
        if (!v.active.load()) continue;
        if (std::find(pcs.begin(), pcs.end(), v.pitch_class) == pcs.end())
            pcs.push_back(v.pitch_class);
    }
    if (pcs.empty()) { obj.chord_name = "—"; obj.quality = ""; obj.confidence = 0.f; return; }
    if (pcs.size() == 1) {
        obj.chord_name = std::string(NOTE_NAMES[pcs[0]]);
        obj.quality = "maj";
        obj.root_pc = pcs[0]; obj.confidence = 1.f; return;
    }

    // Try every rotation (potential root) against every template
    float best_score = -1.f;
    int best_root = pcs[0];
    const ChordTemplate* best_tpl = nullptr;

    for (int root : pcs) {
        // intervals relative to this root
        std::vector<int> intervals;
        for (int pc : pcs) if (pc != root) intervals.push_back(((pc - root) + 12) % 12);
        std::sort(intervals.begin(), intervals.end());

        for (auto& tpl : CHORD_TEMPLATES) {
            // count matching intervals
            int matches = 0;
            for (int iv : tpl.intervals)
                if (std::find(intervals.begin(), intervals.end(), iv % 12) != intervals.end())
                    matches++;
            float score = (float)matches / std::max((int)tpl.intervals.size(), (int)intervals.size());
            if (score > best_score) {
                best_score = score; best_root = root; best_tpl = &tpl;
            }
        }
    }

    obj.root_pc = best_root;
    obj.confidence = best_score;
    if (best_tpl) {
        obj.chord_name = std::string(NOTE_NAMES[best_root]) + best_tpl->name;
        obj.quality = best_tpl->name;
    } else {
        obj.chord_name = std::string(NOTE_NAMES[best_root]) + "?";
        obj.quality = "maj";
    }

    // Tonnetz centroid
    int sx = 0, sy = 0;
    for (auto& v : voices_) { if (!v.active.load()) continue; sx += v.tonnetz_x; sy += v.tonnetz_y; }
    int n = (int)voices_.size();
    obj.tonnetz_x = n ? sx/n : 0;
    obj.tonnetz_y = n ? sy/n : 0;
}

// ────────────────────────────────────────────────────────────────────────────
//  Master roughness + abstract object update
// ────────────────────────────────────────────────────────────────────────────
void AudioEngine::computeRoughness() {
    std::vector<RoughnessRecord> records;
    AbstractObject obj;
    obj.voice_ids.clear();
    float total_roughness = 0.f;

    {
        std::lock_guard<std::mutex> lk(voices_mutex_);
        // pairwise roughness
        for (size_t i = 0; i < voices_.size(); i++) {
            if (!voices_[i].active.load()) continue;
            obj.voice_ids.push_back(voices_[i].id);
            for (size_t j = i+1; j < voices_.size(); j++) {
                if (!voices_[j].active.load()) continue;
                RoughnessRecord r;
                r.voice_a = voices_[i].id;
                r.voice_b = voices_[j].id;
                r.roughness = pairRoughness(voices_[i], voices_[j]);
                r.virtual_pitch_hz = 0.f;
                records.push_back(r);
                total_roughness += r.roughness;
            }
        }

        obj.roughness_total = total_roughness;
        computeVirtualPitch(obj);
        identifyChord(obj);
    }

    {
        std::lock_guard<std::mutex> lk(results_mutex_);
        roughness_records_ = records;
        abstract_object_  = obj;
    }

    if (roughness_cb_) roughness_cb_(records, obj);
}

// ────────────────────────────────────────────────────────────────────────────
//  Running spectrum snapshot (simple summed magnitude for display)
// ────────────────────────────────────────────────────────────────────────────
void AudioEngine::updateSpectrum(const float* /*buf*/, int /*frames*/) {
    // Very lightweight: accumulate squared amplitude in freq bins
    (void)0;
    static int   count = 0;

    // DFT on individual partials (cheaper than FFT on audio buf)
    if (++count < 8) return;
    count = 0;

    SpectrumSnapshot snap;
    snap.magnitude.fill(0.f);
    snap.time = 0.0;

    std::lock_guard<std::mutex> lk(voices_mutex_);
    for (auto& v : voices_) {
        if (!v.active.load()) continue;
        for (int k = 1; k <= MAX_HARMONICS; k++) {
            if (v.harmonic_amp[k-1] < 1e-6f) continue;
            float f = static_cast<float>(v.frequency) * k;
            float amp = v.harmonic_amp[k-1] * v.amplitude * v.env_value;
            // map freq to bin: linear 0..Nyquist → 0..SPECTRUM_SIZE/2
            int bin = static_cast<int>(f / (SAMPLE_RATE/2) * (SPECTRUM_SIZE/2));
            if (bin >= 0 && bin < SPECTRUM_SIZE/2) snap.magnitude[bin] += amp;
        }
    }

    std::lock_guard<std::mutex> lk2(results_mutex_);
    spectrum_snap_ = snap;
}
