#pragma once
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Browser.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Pack.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Spinner.H>
#include <FL/Fl_Choice.H>
#include <vector>
#include <string>
#include <functional>
#include "voice.h"

struct InstrumentNote {
    int pc;
    double freq;
    int oct;
};

struct InstrumentChord {
    std::string name;
    std::vector<InstrumentNote> notes;
};

struct ArpPattern {
    std::string name;
    std::vector<int> steps; // Indices into the current chord's notes
    float division{0.125f}; // 1/8 note
};

class InstrumentWindow : public Fl_Double_Window {
public:
    InstrumentWindow(int W, int H, const char* L = "Harmonia Instrument Builder");

    void setChord(int slot, const InstrumentChord& chord);
    const InstrumentChord& getChord(int slot) const { return chords_[slot % 10]; }

    void setPlayChordCallback(std::function<void(const InstrumentChord&, bool)> cb) { play_chord_cb_ = cb; }
    void setReleaseChordCallback(std::function<void(const InstrumentChord&)> cb) { release_chord_cb_ = cb; }
    void setPlayNoteCallback(std::function<void(double freq, int pc, int oct)> cb) { play_note_cb_ = cb; }

    int handle(int event) override;

    void updateChordList();

private:
    InstrumentChord chords_[10];
    ArpPattern arp_patterns_[10];
    int current_arp_pattern_idx_{0};

    Fl_Browser* browser_chords_{nullptr};
    Fl_Button* btn_clear_slot_{nullptr};
    Fl_Check_Button* btn_sustain_{nullptr};

    // Arp editor
    Fl_Browser* browser_arps_{nullptr};
    Fl_Input* input_arp_steps_{nullptr};
    Fl_Spinner* sp_arp_division_{nullptr};

    std::function<void(const InstrumentChord&, bool)> play_chord_cb_;
    std::function<void(const InstrumentChord&)> release_chord_cb_;
    std::function<void(double freq, int pc, int oct)> play_note_cb_;

    void buildUI();
    static void cbChordBrowser(Fl_Widget* w, void* d);
    static void cbArpBrowser(Fl_Widget* w, void* d);
    static void cbArpSteps(Fl_Widget* w, void* d);

    void playArp(int pattern_idx);
};
