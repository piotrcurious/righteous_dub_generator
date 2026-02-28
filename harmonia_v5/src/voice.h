#pragma once
#include <string>
#include <array>
#include <vector>
#include <cmath>
#include <atomic>

// ────────────────────────────────────────────────────────────────────────────
//  VOICE — a single additive-synthesis voice with a harmonic series.
//
//  Signal model:
//    y(t) = Σₖ  amplitude[k] · sin(2π · k · frequency · t + phase[k])
//
//  k=1 is the fundamental; k=2..MAX_HARMONICS are overtones.
//  The harmonic profile (amplitude[k]) encodes timbre.
// ────────────────────────────────────────────────────────────────────────────

constexpr int MAX_HARMONICS = 16;
constexpr double SAMPLE_RATE = 44100.0;

enum class TimbrePreset {
    SINE,       // only fundamental
    SAWTOOTH,   // 1/k rolloff
    SQUARE,     // odd harmonics 1/k
    STRINGS,    // 1/k with slight inharmonicity
    BRASS,      // rising then falling
    FLUTE,      // fundamental + weak second harmonic
    CUSTOM
};

struct Voice {
    // ── identity
    std::string name;
    int         id;           // unique across session

    // ── pitch
    double frequency;         // fundamental Hz (concert pitch A4=440)
    int    pitch_class;       // 0..edo-1
    int    octave;            // MIDI octave (4 = middle C octave)
    int    edo{12};
    float  detune_cents;      // fine tuning in cents (−50 to +50)

    // ── amplitude / envelope
    float  amplitude;         // master gain 0..1
    float  attack_ms;         // ADSR
    float  decay_ms;
    float  sustain_level;
    float  release_ms;

    // ── harmonic series
    TimbrePreset timbre;
    std::array<float, MAX_HARMONICS> harmonic_amp;   // amplitude of kth harmonic
    std::array<float, MAX_HARMONICS> harmonic_phase;  // current phase (radians)
    std::array<float, MAX_HARMONICS> harmonic_phase_inc; // phase increment per sample

    // ── state
    std::atomic<bool> active{false};
    std::atomic<bool> note_on{false};
    float env_value{0.f};     // current envelope level
    enum class EnvStage { IDLE, ATTACK, DECAY, SUSTAIN, RELEASE } env_stage{EnvStage::IDLE};

    // ── Tonnetz coordinates (5-limit JI lattice)
    int tonnetz_x{0};   // fifths axis (×3/2)
    int tonnetz_y{0};   // thirds axis (×5/4)

    static void pcColorHSV(int pc, int edo_val, float& r, float& g, float& b) {
        r = g = b = 1.0f;
        if (edo_val <= 0) return;
        float h = (float)(pc % edo_val) / edo_val;
        float s = 0.75f, v = 0.95f;
        int i = (int)(h * 6);
        float f = h * 6 - i;
        float p = v * (1 - s);
        float q = v * (1 - f * s);
        float t = v * (1 - (1 - f) * s);
        switch (i % 6) {
            case 0: r = v; g = t; b = p; break;
            case 1: r = q; g = v; b = p; break;
            case 2: r = p; g = v; b = t; break;
            case 3: r = p; g = q; b = v; break;
            case 4: r = t; g = p; b = v; break;
            case 5: r = v; g = p; b = q; break;
        }
    }

    // ── display
    float color[3]{0.5f, 0.7f, 1.0f};   // RGB for UI

    // ────── construction
    Voice() {
        harmonic_amp.fill(0.f);
        harmonic_phase.fill(0.f);
        harmonic_phase_inc.fill(0.f);
        amplitude = 0.6f;
        attack_ms = 20.f; decay_ms = 80.f; sustain_level = 0.7f; release_ms = 300.f;
        timbre = TimbrePreset::SINE;
        setTimbre(timbre);
    }

    // std::atomic is not copyable — provide explicit copy/move
    Voice(const Voice& o)
        : name(o.name), id(o.id), frequency(o.frequency), pitch_class(o.pitch_class),
          octave(o.octave), detune_cents(o.detune_cents), amplitude(o.amplitude),
          attack_ms(o.attack_ms), decay_ms(o.decay_ms), sustain_level(o.sustain_level),
          release_ms(o.release_ms), timbre(o.timbre),
          harmonic_amp(o.harmonic_amp), harmonic_phase(o.harmonic_phase),
          harmonic_phase_inc(o.harmonic_phase_inc),
          active(o.active.load()), note_on(o.note_on.load()),
          env_value(o.env_value), env_stage(o.env_stage),
          tonnetz_x(o.tonnetz_x), tonnetz_y(o.tonnetz_y) {
        color[0]=o.color[0]; color[1]=o.color[1]; color[2]=o.color[2];
    }

    Voice& operator=(const Voice& o) {
        if (this == &o) return *this;
        name=o.name; id=o.id; frequency=o.frequency; pitch_class=o.pitch_class;
        octave=o.octave; detune_cents=o.detune_cents; amplitude=o.amplitude;
        attack_ms=o.attack_ms; decay_ms=o.decay_ms; sustain_level=o.sustain_level;
        release_ms=o.release_ms; timbre=o.timbre;
        harmonic_amp=o.harmonic_amp; harmonic_phase=o.harmonic_phase;
        harmonic_phase_inc=o.harmonic_phase_inc;
        active.store(o.active.load()); note_on.store(o.note_on.load());
        env_value=o.env_value; env_stage=o.env_stage;
        tonnetz_x=o.tonnetz_x; tonnetz_y=o.tonnetz_y;
        color[0]=o.color[0]; color[1]=o.color[1]; color[2]=o.color[2];
        return *this;
    }

    Voice(Voice&& o) noexcept
        : name(std::move(o.name)), id(o.id), frequency(o.frequency),
          pitch_class(o.pitch_class), octave(o.octave), detune_cents(o.detune_cents),
          amplitude(o.amplitude), attack_ms(o.attack_ms), decay_ms(o.decay_ms),
          sustain_level(o.sustain_level), release_ms(o.release_ms), timbre(o.timbre),
          harmonic_amp(o.harmonic_amp), harmonic_phase(o.harmonic_phase),
          harmonic_phase_inc(o.harmonic_phase_inc),
          active(o.active.load()), note_on(o.note_on.load()),
          env_value(o.env_value), env_stage(o.env_stage),
          tonnetz_x(o.tonnetz_x), tonnetz_y(o.tonnetz_y) {
        color[0]=o.color[0]; color[1]=o.color[1]; color[2]=o.color[2];
    }

    Voice& operator=(Voice&& o) noexcept { return operator=(o); }

    // ────── set frequency and recompute phase increments
    void setFrequency(double freq, bool update_coords = true) {
        frequency = freq;
        double f_actual = freq * std::pow(2.0, detune_cents / 1200.0);
        for (int k = 1; k <= MAX_HARMONICS; k++) {
            harmonic_phase_inc[k-1] = static_cast<float>(
                2.0 * M_PI * k * f_actual / SAMPLE_RATE
            );
        }

        // For pitch class, we use the pure logf without detune for theoretical labeling
        double logf_label = std::log2(freq / 261.63);
        logf_label = logf_label - std::floor(logf_label); // mod octave
        pitch_class = (int)std::round(logf_label * edo) % edo;

        if (update_coords) computeTonnetzCoords();
    }

    // ────── set pitch by MIDI note number
    void setMidiNote(int midi_note) {
        octave     = midi_note / 12 - 1;
        setFrequency(440.0 * std::pow(2.0, (midi_note - 69) / 12.0));
    }

    // ────── apply a timbre preset to harmonic amplitudes
    void setTimbre(TimbrePreset t) {
        timbre = t;
        harmonic_amp.fill(0.f);
        switch (t) {
        case TimbrePreset::SINE:
            harmonic_amp[0] = 1.0f;
            break;
        case TimbrePreset::SAWTOOTH:
            for (int k = 1; k <= MAX_HARMONICS; k++)
                harmonic_amp[k-1] = 1.0f / k;
            break;
        case TimbrePreset::SQUARE:
            for (int k = 1; k <= MAX_HARMONICS; k += 2)
                harmonic_amp[k-1] = 1.0f / k;
            break;
        case TimbrePreset::STRINGS:
            for (int k = 1; k <= MAX_HARMONICS; k++)
                harmonic_amp[k-1] = std::exp(-0.15f * (k-1)) / k;
            harmonic_amp[1] *= 1.4f;
            harmonic_amp[2] *= 1.1f;
            break;
        case TimbrePreset::BRASS:
            for (int k = 1; k <= MAX_HARMONICS; k++) {
                float env = (k <= 6) ? k/6.f : std::exp(-0.3f*(k-6));
                harmonic_amp[k-1] = env / std::sqrt(static_cast<float>(k));
            }
            break;
        case TimbrePreset::FLUTE:
            harmonic_amp[0] = 1.0f;
            harmonic_amp[1] = 0.3f;
            harmonic_amp[2] = 0.05f;
            break;
        case TimbrePreset::CUSTOM:
            break;
        }
        // normalize
        float peak = 0.f;
        for (float a : harmonic_amp) peak = std::max(peak, a);
        if (peak > 0.f) for (float& a : harmonic_amp) a /= peak;
    }

    // ────── map frequency to nearest Tonnetz coords (5-limit JI lattice)
    void computeTonnetzCoords() {
        // Project log2(f/C4) onto JI fifth (log2(3/2)) and third (log2(5/4)) axes
        double logf = std::log2(frequency / 261.63);
        logf = logf - std::floor(logf); // mod octave

        double fifth_ji = std::log2(1.5);
        double third_ji = std::log2(1.25);

        double best = 1e9;
        for (int a = -10; a <= 10; a++) {
            for (int b = -6; b <= 6; b++) {
                double val = a * fifth_ji + b * third_ji;
                val = val - std::floor(val);
                double diff = std::min(std::abs(val - logf),
                             std::min(std::abs(val - logf + 1),
                                      std::abs(val - logf - 1)));
                if (diff < best) { best = diff; tonnetz_x = a; tonnetz_y = b; }
            }
        }
    }
};

// ────── note name utilities
static const char* NOTE_NAMES[12] = {
    "C","D♭","D","E♭","E","F","G♭","G","A♭","A","B♭","B"
};

inline std::string noteName(int pitch_class, int octave) {
    return std::string(NOTE_NAMES[pitch_class]) + std::to_string(octave);
}

inline double midiToHz(int midi_note) {
    return 440.0 * std::pow(2.0, (midi_note - 69) / 12.0);
}

inline int hzToMidi(double hz) {
    return static_cast<int>(std::round(69.0 + 12.0 * std::log2(hz / 440.0)));
}
