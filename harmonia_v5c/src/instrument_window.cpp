#include "instrument_window.h"
#include <FL/fl_draw.H>
#include <FL/fl_ask.H>
#include <sstream>
#include <algorithm>

static const Fl_Color COL_BG     = fl_rgb_color(18,20,28);
static const Fl_Color COL_ACCENT = fl_rgb_color(200,160,40);

InstrumentWindow::InstrumentWindow(int W, int H, const char* L)
    : Fl_Double_Window(W, H, L) {
    color(COL_BG);
    buildUI();
    end(); // Ensure window is closed

    // Default chords and patterns
    for(int i=0; i<10; i++) {
        chords_[i].name = "Empty " + std::to_string((i+1)%10);
        arp_patterns_[i].name = "Pattern " + std::to_string((i+1)%10);
        arp_patterns_[i].steps = {0, 1, 2, 1}; // Default pattern
    }
    updateChordList();
}

void InstrumentWindow::buildUI() {
    Fl_Tabs* tabs = new Fl_Tabs(10, 10, w()-20, h()-20);

    // Chord Tab
    {
        Fl_Group* g = new Fl_Group(10, 35, w()-20, h()-45, "Chords");
        browser_chords_ = new Fl_Browser(20, 45, w()-40, h()-120);
        browser_chords_->type(FL_HOLD_BROWSER);
        browser_chords_->callback(cbChordBrowser, this);

        btn_clear_slot_ = new Fl_Button(20, h()-70, 100, 25, "Clear Slot");
        btn_sustain_ = new Fl_Check_Button(130, h()-70, 100, 25, "Sustain");
        g->end();
    }

    // Arp Tab
    {
        Fl_Group* g = new Fl_Group(10, 35, w()-20, h()-45, "Arpeggio");
        browser_arps_ = new Fl_Browser(20, 45, 150, h()-120);
        browser_arps_->type(FL_HOLD_BROWSER);
        browser_arps_->callback(cbArpBrowser, this);

        new Fl_Box(180, 45, 100, 20, "Steps (comma separated):");
        input_arp_steps_ = new Fl_Input(180, 70, w()-200, 25);
        input_arp_steps_->callback(cbArpSteps, this);
        input_arp_steps_->when(FL_WHEN_CHANGED);

        new Fl_Box(180, 105, 100, 20, "Division (s):");
        sp_arp_division_ = new Fl_Spinner(180, 130, 80, 25);
        sp_arp_division_->type(FL_FLOAT_INPUT);
        sp_arp_division_->minimum(0.01);
        sp_arp_division_->maximum(2.0);
        sp_arp_division_->step(0.01);
        sp_arp_division_->value(0.125);

        g->end();
    }

    tabs->end();

    for(int i=0; i<10; i++) {
        char buf[32];
        snprintf(buf, 32, "Slot %d: (empty)", (i+1)%10);
        browser_chords_->add(buf);

        snprintf(buf, 32, "Pattern %d", (i+1)%10);
        browser_arps_->add(buf);
    }
}

void InstrumentWindow::updateChordList() {
    browser_chords_->clear();
    for(int i=0; i<10; i++) {
        char buf[256];
        snprintf(buf, 256, "Slot %d: %s (%zu notes)", (i+1)%10, chords_[i].name.c_str(), chords_[i].notes.size());
        browser_chords_->add(buf);
    }
}

void InstrumentWindow::setChord(int slot, const InstrumentChord& chord) {
    if(slot < 0 || slot >= 10) return;
    chords_[slot] = chord;
    updateChordList();
}

int InstrumentWindow::handle(int event) {
    if (handlePerformanceKeys(event)) return 1;
    return Fl_Double_Window::handle(event);
}

int InstrumentWindow::handlePerformanceKeys(int event) {
    if (event == FL_KEYDOWN) {
        if (Fl::focus() && (dynamic_cast<Fl_Input*>(Fl::focus()) || dynamic_cast<Fl_Spinner*>(Fl::focus())))
            return 0;
        int k = Fl::event_key();
        if (k >= '1' && k <= '9') {
            int idx = k - '1';
            if (play_chord_cb_) play_chord_cb_(chords_[idx], btn_sustain_->value());
            return 1;
        }
        if (k == '0') {
            if (play_chord_cb_) play_chord_cb_(chords_[9], btn_sustain_->value());
            return 1;
        }

        // Arp keys: zxcvbnm,./
        const char* arp_keys = "zxcvbnm,./";
        const char* p = strchr(arp_keys, k);
        if (p) {
            int idx = p - arp_keys;
            playArp(idx);
            return 1;
        }
    }
    if (event == FL_KEYUP) {
        if (Fl::focus() && (dynamic_cast<Fl_Input*>(Fl::focus()) || dynamic_cast<Fl_Spinner*>(Fl::focus())))
            return 0;
        int k = Fl::event_key();
        if (k >= '1' && k <= '9') {
            int idx = k - '1';
            if (!btn_sustain_->value() && release_chord_cb_) release_chord_cb_(chords_[idx]);
            return 1;
        }
        if (k == '0') {
            if (!btn_sustain_->value() && release_chord_cb_) release_chord_cb_(chords_[9]);
            return 1;
        }
    }
    return 0;
}

void InstrumentWindow::cbChordBrowser(Fl_Widget* w, void* d) {
    (void)w; (void)d;
    // Maybe highlight or something
}

void InstrumentWindow::cbArpBrowser(Fl_Widget* w, void* d) {
    (void)w;
    InstrumentWindow* self = (InstrumentWindow*)d;
    int idx = self->browser_arps_->value() - 1;
    if (idx >= 0 && idx < 10) {
        self->current_arp_pattern_idx_ = idx;
        std::stringstream ss;
        for(size_t i=0; i<self->arp_patterns_[idx].steps.size(); i++) {
            ss << self->arp_patterns_[idx].steps[i] << (i == self->arp_patterns_[idx].steps.size()-1 ? "" : ",");
        }
        self->input_arp_steps_->value(ss.str().c_str());
        self->sp_arp_division_->value(self->arp_patterns_[idx].division);
    }
}

void InstrumentWindow::cbArpSteps(Fl_Widget* w, void* d) {
    (void)w;
    InstrumentWindow* self = (InstrumentWindow*)d;
    int idx = self->current_arp_pattern_idx_;
    std::string s = self->input_arp_steps_->value();
    std::vector<int> steps;
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, ',')) {
        try {
            steps.push_back(std::stoi(item));
        } catch(...) {}
    }
    self->arp_patterns_[idx].steps = steps;
}

struct ArpState {
    InstrumentWindow* win;
    InstrumentChord chord;
    ArpPattern pattern;
    size_t current_step;
    std::function<void(double, int, int)> play_note_cb;
};

void arp_timer_cb(void* d) {
    ArpState* state = (ArpState*)d;
    if (state->current_step >= state->pattern.steps.size() || state->chord.notes.empty()) {
        delete state;
        return;
    }

    int note_idx = state->pattern.steps[state->current_step];
    if (note_idx >= 0 && note_idx < (int)state->chord.notes.size()) {
        const auto& n = state->chord.notes[note_idx];
        if (state->play_note_cb) state->play_note_cb(n.freq, n.pc, n.oct);
    }

    state->current_step++;
    if (state->current_step < state->pattern.steps.size()) {
        Fl::repeat_timeout(state->pattern.division, arp_timer_cb, state);
    } else {
        delete state;
    }
}

void InstrumentWindow::playArp(int pattern_idx) {
    int chord_idx = browser_chords_->value() - 1;
    if (chord_idx < 0) chord_idx = 0;

    const auto& chord = chords_[chord_idx];
    if (chord.notes.empty()) return;

    const auto& pattern = arp_patterns_[pattern_idx];
    if (pattern.steps.empty()) return;

    ArpState* state = new ArpState{this, chord, pattern, 0, play_note_cb_};
    state->pattern.division = (float)sp_arp_division_->value();
    Fl::add_timeout(0.0, arp_timer_cb, state);
}
