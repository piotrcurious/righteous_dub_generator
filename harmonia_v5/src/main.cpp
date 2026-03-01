/*
 * harmonia/src/main.cpp
 * ─────────────────────────────────────────────────────────────────────────────
 *  HARMONIA V5 — Psychoacoustic Multi-Voice Music Coder
 *  Restored full V4 feature set with Circular Tonal Space visualization.
 * ─────────────────────────────────────────────────────────────────────────────
 */

#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Browser.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Spinner.H>
#include <FL/fl_draw.H>
#include <FL/fl_ask.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>
#include <GL/gl.h>
#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include <functional>
#include <map>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <unistd.h>
#include <algorithm>

#include "voice.h"
#include "audio_engine.h"
#include "tonal_space_widget.h"
#include "theory_bridge.h"

// ─────────────────────────────────────────────────────────────────────────────
//  FLTK colour helpers
// ─────────────────────────────────────────────────────────────────────────────
static const Fl_Color COL_BG     = fl_rgb_color(18,20,28);
static const Fl_Color COL_PANEL  = fl_rgb_color(26,30,40);
static const Fl_Color COL_ACCENT = fl_rgb_color(200,160,40);
static const Fl_Color COL_TEXT   = fl_rgb_color(220,215,200);
static const Fl_Color COL_DIM    = fl_rgb_color(100,100,110);

// ─────────────────────────────────────────────────────────────────────────────
//  SPECTRUM DISPLAY
// ─────────────────────────────────────────────────────────────────────────────
class SpectrumWidget : public Fl_Gl_Window {
public:
    SpectrumWidget(int X,int Y,int W,int H)
        : Fl_Gl_Window(X,Y,W,H)
    { mode(FL_RGB|FL_DOUBLE); }

    void update(const SpectrumSnapshot& snap, const std::vector<Voice>& voices,
                float roughness) {
        snap_ = snap; voices_ = voices; roughness_ = roughness; redraw();
    }

    void draw() override {
        if (!valid()) {
            glViewport(0,0,w(),h());
            glMatrixMode(GL_PROJECTION); glLoadIdentity();
            glOrtho(0,w(),h(),0,-1,1);
            glMatrixMode(GL_MODELVIEW); glLoadIdentity();
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        }
        glClearColor(0.07f,0.08f,0.11f,1.f);
        glClear(GL_COLOR_BUFFER_BIT);

        // Draw spectrum bars
        int sw = w(), sh = h()-20;
        float bin_w = (float)sw / (SPECTRUM_SIZE/2);
        for (int i = 0; i < SPECTRUM_SIZE/2; i++) {
            float mag = snap_.magnitude[i];
            if (mag < 1e-5f) continue;
            float bh = std::min((float)sh, mag * sh * 3.f);
            float t = (float)i / (SPECTRUM_SIZE/2);
            float r = 0.2f + 0.8f*t, g = 0.7f - 0.5f*t, b = 1.f - 0.7f*t;
            glColor4f(r,g,b, 0.85f);
            glBegin(GL_QUADS);
            glVertex2f(i*bin_w, sh);
            glVertex2f((i+1)*bin_w, sh);
            glVertex2f((i+1)*bin_w, sh-bh);
            glVertex2f(i*bin_w, sh-bh);
            glEnd();
        }

        // Roughness bar at bottom
        float rw = roughness_ * sw;
        float rg = 1.f - roughness_;
        glColor4f(roughness_, rg, 0.f, 0.8f);
        glBegin(GL_QUADS);
        glVertex2f(0,sh); glVertex2f(rw,sh);
        glVertex2f(rw,h()); glVertex2f(0,h());
        glEnd();

        gl_font(FL_HELVETICA, 9);
        glColor3f(.6f,.6f,.6f);
        gl_draw("Roughness", 2, h()-3);
        char rbuf[32]; snprintf(rbuf,sizeof(rbuf),"%.2f", roughness_);
        gl_draw(rbuf, (int)(rw + 4), h()-3);
    }

private:
    SpectrumSnapshot snap_;
    std::vector<Voice> voices_;
    float roughness_{0.f};
};

// ─────────────────────────────────────────────────────────────────────────────
//  VOICE STRIP — one row per voice in the left panel
// ─────────────────────────────────────────────────────────────────────────────
struct VoiceStrip {
    int voice_id{-1};
    bool manual{false};
    Fl_Group*       group{nullptr};
    Fl_Light_Button* btn_on{nullptr};
    Fl_Box*          lbl_note{nullptr};
    Fl_Value_Slider* sl_amp{nullptr};
    Fl_Value_Slider* sl_detune{nullptr};
    Fl_Choice*       ch_timbre{nullptr};
    Fl_Button*       btn_remove{nullptr};
    Fl_Box*          roughness_box{nullptr};
    float            roughness_val{0.f};
    float            color[3]{0.5f,0.7f,1.f};
};

// ─────────────────────────────────────────────────────────────────────────────
//  MAIN WINDOW
// ─────────────────────────────────────────────────────────────────────────────
class HarmoniaApp {
public:
    HarmoniaApp();
    ~HarmoniaApp();
    void run();

private:
    static constexpr int WIN_W = 1280, WIN_H = 820;
    static constexpr int TOOLBAR_H = 38;
    static constexpr int VOICE_PANEL_W = 260;
    static constexpr int THEORY_PANEL_W = 320;
    static constexpr int STATUS_H = 120;
    static constexpr int TONAL_W = WIN_W - VOICE_PANEL_W - THEORY_PANEL_W;
    static constexpr int TONAL_H = WIN_H - TOOLBAR_H - STATUS_H;

    // ── FLTK widgets
    Fl_Double_Window* win_{nullptr};

    // Toolbar
    Fl_Button*       btn_add_{nullptr};
    Fl_Button*       btn_play_all_{nullptr};
    Fl_Button*       btn_stop_all_{nullptr};
    Fl_Choice*       ch_key_{nullptr};
    Fl_Spinner*      sp_edo_{nullptr};
    Fl_Spinner*      sp_base_octave_{nullptr};
    Fl_Value_Slider* sl_master_{nullptr};
    Fl_Output*       out_object_{nullptr};

    // Left panel – voice list
    Fl_Scroll*        scroll_voices_{nullptr};
    std::vector<VoiceStrip*> strips_;
    int               voice_panel_y_{0};

    // Centre – Tonal Space
    TonalSpaceWidget* tonal_space_{nullptr};

    // Right – Side Panel
    Fl_Tabs*     tabs_right_{nullptr};
    Fl_Group*    grp_theory_{nullptr};
    Fl_Group*    grp_modulation_{nullptr};
    Fl_Group*    grp_instrument_{nullptr};

    // Right – structural theory panel
    Fl_Browser*  browser_function_{nullptr};
    Fl_Browser*  browser_resolutions_{nullptr};
    Fl_Browser*  browser_completion_{nullptr};
    Fl_Browser*  browser_psycho_{nullptr};
    Fl_Browser*  browser_edo_{nullptr};
    Fl_Button*   btn_analyze_{nullptr};
    Fl_Button*   btn_resolve_{nullptr};
    Fl_Button*   btn_complete_{nullptr};
    Fl_Button*   btn_psycho_{nullptr};
    Fl_Button*   btn_edo_{nullptr};
    Fl_Output*   out_orbifold_{nullptr};

    // Right – modulation panel
    Fl_Choice*   ch_target_key_{nullptr};
    Fl_Button*   btn_find_pivots_{nullptr};
    Fl_Browser*  browser_pivots_{nullptr};
    PivotSearchResult last_pivot_res_;

    // Instrument Builder
    struct ChordKey {
        std::string name;
        struct Note { int pc; double freq; int oct; };
        std::vector<Note> notes;
    };
    std::vector<ChordKey> instrument_keyboard_;
    Fl_Browser* browser_keys_{nullptr};
    Fl_Button* btn_add_chord_{nullptr};
    Fl_Button* btn_clear_keys_{nullptr};

    // Bottom – spectrum + roughness
    SpectrumWidget*   spectrum_{nullptr};
    Fl_Box*           box_object_label_{nullptr};

    // ── engines
    std::unique_ptr<AudioEngine>  audio_;
    std::unique_ptr<TheoryBridge> theory_;

    // ── state
    int  current_key_{0};
    int  current_edo_{12};
    int  last_clicked_root_{-1};
    FunctionalAnalysis  last_func_;
    TonnetzTension      last_tension_;
    std::vector<ResolutionPath> last_resolutions_;

    // ── methods
    void buildUI();
    void setupCallbacks();

    void addVoice(int midi_note = 60, TimbrePreset t = TimbrePreset::SINE);
    void removeVoice(int voice_id);
    VoiceStrip* findStrip(int voice_id);
    void relayoutStrips();
    std::vector<VoiceStrip*> pool_strips_;

    void updateFunctionalAnalysis();
    void updateResolutionPaths();
    void updateCompletionSuggestions();
    void updatePsychoAnalysis();
    void updateEDOAnalysis();
    void updatePivotSearch();
    void updateTheory();

    static void onIdle(void* data);
    static void cbAddVoice(Fl_Widget*, void* d) { ((HarmoniaApp*)d)->addVoice(); }
    static void cbPlayAll(Fl_Widget*, void* d);
    static void cbStopAll(Fl_Widget*, void* d);
    static void cbMasterVol(Fl_Widget*, void* d);
    static void cbKey(Fl_Widget*, void* d);
    static void cbEDO(Fl_Widget*, void* d);
    static void cbAnalyze(Fl_Widget*, void* d);
    static void cbResolve(Fl_Widget*, void* d);
    static void cbSuggestCompletion(Fl_Widget*, void* d);
    static void cbPsycho(Fl_Widget*, void* d);
    static void cbAnalyseEDO(Fl_Widget*, void* d);
    static void cbFindPivots(Fl_Widget*, void* d);
    static void cbBrowserPivots(Fl_Widget*, void* d);
    static void cbBrowserFunction(Fl_Widget*, void* d);
    static void cbBrowserPsycho(Fl_Widget*, void* d);
    static void cbBrowserResolutions(Fl_Widget* w, void* d);
    static void cbBrowserCompletion(Fl_Widget* w, void* d);
    static void cbBrowserKeys(Fl_Widget* w, void* d);
    static void cbAddChordToKeys(Fl_Widget*, void* d);
    static void cbClearKeys(Fl_Widget*, void* d);

    void playChord(const ChordKey& chord, bool sustain);
    void releaseChord(const ChordKey& chord);
    int handle(int event);

    void onTonalSpaceClick(int pc, int oct);
    std::pair<int,int> highlightNodeForPC(int pc, int oct = -100);

    std::string theoryDir_;
};

HarmoniaApp::HarmoniaApp() {
    audio_  = std::make_unique<AudioEngine>();
    theory_ = std::make_unique<TheoryBridge>();
    buildUI();
    setupCallbacks();
}

HarmoniaApp::~HarmoniaApp() {
    audio_->shutdown();
    theory_->stop();
}

void HarmoniaApp::run() {
    if (!audio_->init()) fl_message("ALSA init failed — audio disabled");
    char exepath[512] = {};
    ssize_t n = readlink("/proc/self/exe", exepath, sizeof(exepath)-1);
    theoryDir_ = (n > 0) ? std::string(exepath, n).substr(0, std::string(exepath, n).rfind('/')) + "/../theory" : "./theory";
    if (!theory_->start(theoryDir_)) fl_message("Python theory server not started");

    addVoice(60); addVoice(64); addVoice(67);
    for (auto& s : strips_) audio_->noteOn(s->voice_id);

    updateTheory();
    Fl::add_idle(onIdle, this);
    Fl::run();
}

void HarmoniaApp::buildUI() {
    Fl::scheme("gtk+");
    Fl::background(18,20,28);
    Fl::foreground(220,215,200);

    class HarmoniaWindow : public Fl_Double_Window {
        HarmoniaApp* app_;
    public:
        HarmoniaWindow(int W, int H, const char* L, HarmoniaApp* a) : Fl_Double_Window(W, H, L), app_(a) {}
        int handle(int e) override { if (app_->handle(e)) return 1; return Fl_Double_Window::handle(e); }
    };

    win_ = new HarmoniaWindow(WIN_W, WIN_H, "Harmonia V5 — Circular Tonal Space", this);
    win_->color(COL_BG);

    // ── TOOLBAR ──────────────────────────────────────────────────────────────
    {
        Fl_Group* tb = new Fl_Group(0, 0, WIN_W, TOOLBAR_H);
        tb->box(FL_FLAT_BOX); tb->color(fl_rgb_color(14,16,22));
        int x = 8, y = 6, bh = 26;
        Fl_Box* title = new Fl_Box(x,y,100,bh,"HARMONIA V5");
        title->labelcolor(COL_ACCENT); title->labelfont(FL_HELVETICA_BOLD); x += 110;
        new Fl_Box(x,y,30,bh,"Key:"); x+=30;
        ch_key_ = new Fl_Choice(x,y,70,bh); for(int i=0;i<12;i++) ch_key_->add(NOTE_NAMES[i]); ch_key_->value(0); x += 76;
        new Fl_Box(x,y,30,bh,"EDO:"); x+=32;
        sp_edo_ = new Fl_Spinner(x,y,50,bh); sp_edo_->minimum(5); sp_edo_->maximum(72); sp_edo_->value(12); x += 56;
        new Fl_Box(x,y,30,bh,"Oct:"); x+=32;
        sp_base_octave_ = new Fl_Spinner(x,y,45,bh); sp_base_octave_->minimum(0); sp_base_octave_->maximum(8); sp_base_octave_->value(4); x += 50;
        new Fl_Box(x,y,50,bh,"Master:"); x+=52;
        sl_master_ = new Fl_Value_Slider(x,y,100,bh); sl_master_->type(FL_HORIZONTAL); sl_master_->value(0.5); x += 110;
        btn_add_ = new Fl_Button(x,y,90,bh,"+ Add Voice"); x += 98;
        btn_play_all_ = new Fl_Button(x,y,70,bh,"▶ All"); x += 76;
        btn_stop_all_ = new Fl_Button(x,y,70,bh,"■ Stop"); x += 76;
        out_object_ = new Fl_Output(x,y, WIN_W-x-8, bh);
        out_object_->box(FL_FLAT_BOX); out_object_->color(fl_rgb_color(20,24,32)); out_object_->textcolor(COL_ACCENT);
        tb->end();
    }

    int content_y = TOOLBAR_H;
    int content_h = WIN_H - TOOLBAR_H - STATUS_H;

    // ── VOICE PANEL (left) ────────────────────────────────────────────────────
    {
        Fl_Group* panel = new Fl_Group(0, content_y, VOICE_PANEL_W, content_h);
        panel->box(FL_FLAT_BOX); panel->color(fl_rgb_color(22,26,35));
        new Fl_Box(0,content_y, VOICE_PANEL_W,22,"VOICES");
        scroll_voices_ = new Fl_Scroll(0, content_y+22, VOICE_PANEL_W, content_h-22);
        scroll_voices_->type(Fl_Scroll::VERTICAL_ALWAYS); scroll_voices_->end();
        panel->end();
    }

    // ── TONAL SPACE (centre) ──────────────────────────────────────────────────
    tonal_space_ = new TonalSpaceWidget(VOICE_PANEL_W, content_y, TONAL_W, content_h);

    // ── SIDE PANEL (right) ────────────────────────────────────────────────────
    {
        tabs_right_ = new Fl_Tabs(VOICE_PANEL_W + TONAL_W, content_y, THEORY_PANEL_W, content_h);

        // Tab 1: Theory
        grp_theory_ = new Fl_Group(tabs_right_->x(), content_y + 25, THEORY_PANEL_W, content_h - 25, "Theory");
        int py = content_y + 30, px = tabs_right_->x() + 4, pw = THEORY_PANEL_W - 8;
        new Fl_Box(px, py, pw, 16, "HARMONIC FUNCTION"); py += 18;
        browser_function_ = new Fl_Browser(px, py, pw, 120); browser_function_->textsize(9); py += 122;
        btn_analyze_ = new Fl_Button(px, py, pw, 20, "Analyse Structure"); py += 24;
        new Fl_Box(px, py, pw, 16, "RESOLUTION PATHS"); py += 18;
        browser_resolutions_ = new Fl_Browser(px, py, pw, 120); browser_resolutions_->textsize(9); py += 122;
        btn_resolve_ = new Fl_Button(px, py, pw, 20, "Show Resolutions"); py += 24;
        new Fl_Box(px, py, pw, 16, "COMPLETION"); py += 18;
        browser_completion_ = new Fl_Browser(px, py, pw, 80); browser_completion_->textsize(9); py += 82;
        btn_complete_ = new Fl_Button(px, py, pw, 20, "Suggest Voice"); py += 24;
        new Fl_Box(px, py, pw, 16, "PSYCHOACOUSTICS"); py += 18;
        browser_psycho_ = new Fl_Browser(px, py, pw, 80); browser_psycho_->textsize(9); py += 82;
        btn_psycho_ = new Fl_Button(px, py, pw, 20, "Neural Analysis");
        grp_theory_->end();

        // Tab 2: Modulate
        grp_modulation_ = new Fl_Group(tabs_right_->x(), content_y + 25, THEORY_PANEL_W, content_h - 25, "Modulate");
        py = content_y + 30; new Fl_Box(px, py, pw, 16, "PIVOT SEARCH"); py += 20;
        ch_target_key_ = new Fl_Choice(px+80, py, 70, 20, "Target:"); for(int i=0;i<12;i++) ch_target_key_->add(NOTE_NAMES[i]); ch_target_key_->value(7); py += 25;
        btn_find_pivots_ = new Fl_Button(px, py, pw, 22, "Find Paths"); py += 26;
        browser_pivots_ = new Fl_Browser(px, py, pw, content_h - (py-content_y) - 10); browser_pivots_->textsize(9);
        grp_modulation_->end();

        // Tab 3: Instrument
        grp_instrument_ = new Fl_Group(tabs_right_->x(), content_y + 25, THEORY_PANEL_W, content_h - 25, "Instrument");
        py = content_y + 30; new Fl_Box(px, py, pw, 16, "KEYBOARD BUILDER"); py += 20;
        browser_keys_ = new Fl_Browser(px, py, pw, content_h - 120); py += content_h - 118;
        btn_add_chord_ = new Fl_Button(px, py, pw, 26, "+ Add Current"); py += 30;
        btn_clear_keys_ = new Fl_Button(px, py, pw, 26, "Clear Keys");
        grp_instrument_->end();

        tabs_right_->end();
    }

    // Bottom
    {
        spectrum_ = new SpectrumWidget(0, WIN_H - STATUS_H, WIN_W - 250, STATUS_H);
        box_object_label_ = new Fl_Box(WIN_W - 248, WIN_H - STATUS_H, 248, STATUS_H, "");
        box_object_label_->labelcolor(COL_ACCENT); box_object_label_->labelfont(FL_HELVETICA_BOLD); box_object_label_->labelsize(16);
    }

    win_->end(); win_->show();
}

void HarmoniaApp::setupCallbacks() {
    btn_add_->callback(cbAddVoice, this);
    btn_analyze_->callback(cbAnalyze, this);
    btn_resolve_->callback(cbResolve, this);
    btn_complete_->callback(cbSuggestCompletion, this);
    btn_psycho_->callback(cbPsycho, this);
    btn_find_pivots_->callback(cbFindPivots, this);
    btn_add_chord_->callback(cbAddChordToKeys, this);
    btn_clear_keys_->callback(cbClearKeys, this);
    btn_play_all_->callback(cbPlayAll, this);
    btn_stop_all_->callback(cbStopAll, this);
    sl_master_->callback(cbMasterVol, this);
    ch_key_->callback(cbKey, this);
    sp_edo_->callback(cbEDO, this);

    browser_function_->callback(cbBrowserFunction, this);
    browser_resolutions_->callback(cbBrowserResolutions, this);
    browser_completion_->callback(cbBrowserCompletion, this);
    browser_psycho_->callback(cbBrowserPsycho, this);
    browser_pivots_->callback(cbBrowserPivots, this);
    browser_keys_->callback(cbBrowserKeys, this);

    tonal_space_->setNodeClickCallback([this](int pc, int tx, int ty){ onTonalSpaceClick(pc, ty); (void)tx; });
}

void HarmoniaApp::addVoice(int midi_note, TimbrePreset t) {
    int id = audio_->addVoice(midi_note, t); if (id < 0) return;
    VoiceStrip* s = nullptr;
    if (!pool_strips_.empty()) { s = pool_strips_.back(); pool_strips_.pop_back(); s->group->show(); scroll_voices_->add(s->group); }
    else {
        scroll_voices_->begin(); s = new VoiceStrip();
        s->group = new Fl_Group(4, 0, VOICE_PANEL_W-20, 76); s->group->box(FL_FLAT_BOX); s->group->color(fl_rgb_color(28,32,44));
        s->lbl_note = new Fl_Box(s->group->x()+6, s->group->y()+4, 60, 26, "—"); s->lbl_note->labelfont(FL_HELVETICA_BOLD); s->lbl_note->labelsize(15);
        s->roughness_box = new Fl_Box(s->group->x()+70, s->group->y()+4, 46, 13, "R:0.00"); s->roughness_box->labelsize(8);
        s->btn_on = new Fl_Light_Button(s->group->x()+70, s->group->y()+18, 46, 14, "PLAY"); s->btn_on->labelsize(8);
        s->btn_remove = new Fl_Button(s->group->x()+s->group->w()-22, s->group->y()+4, 18, 18, "×");
        s->ch_timbre = new Fl_Choice(s->group->x()+118, s->group->y()+4, s->group->w()-142, 18);
        s->ch_timbre->add("Sine|Saw|Square|Strings|Brass|Flute");
        s->sl_amp = new Fl_Value_Slider(s->group->x()+30, s->group->y()+36, s->group->w()-34, 13); s->sl_amp->type(FL_HORIZONTAL);
        s->sl_detune = new Fl_Value_Slider(s->group->x()+30, s->group->y()+56, s->group->w()-34, 13); s->sl_detune->type(FL_HORIZONTAL);
        s->group->end(); scroll_voices_->end();
    }
    s->voice_id = id; s->manual = false;
    s->btn_on->user_data((void*)(intptr_t)id); s->btn_on->callback([](Fl_Widget* w, void* d){ int id=(int)(intptr_t)w->user_data(); if(((Fl_Light_Button*)w)->value()) ((HarmoniaApp*)d)->audio_->noteOn(id); else ((HarmoniaApp*)d)->audio_->noteOff(id); }, this);
    s->btn_remove->user_data((void*)(intptr_t)id); s->btn_remove->callback([](Fl_Widget* w, void* d){ ((HarmoniaApp*)d)->removeVoice((int)(intptr_t)w->user_data()); }, this);
    s->ch_timbre->user_data((void*)(intptr_t)id); s->ch_timbre->callback([](Fl_Widget* w, void* d){ ((HarmoniaApp*)d)->audio_->setVoiceTimbre((int)(intptr_t)w->user_data(), (TimbrePreset)((Fl_Choice*)w)->value()); }, this);
    s->sl_amp->user_data((void*)(intptr_t)id); s->sl_amp->callback([](Fl_Widget* w, void* d){ ((HarmoniaApp*)d)->audio_->setVoiceAmplitude((int)(intptr_t)w->user_data(), (float)((Fl_Value_Slider*)w)->value()); }, this);
    s->sl_detune->user_data((void*)(intptr_t)id); s->sl_detune->callback([](Fl_Widget* w, void* d){ ((HarmoniaApp*)d)->audio_->setVoiceDetune((int)(intptr_t)w->user_data(), (float)((Fl_Value_Slider*)w)->value()); }, this);
    s->sl_amp->value(0.6); s->sl_detune->value(0); s->ch_timbre->value((int)t);
    strips_.push_back(s); relayoutStrips(); win_->redraw();
}

void HarmoniaApp::removeVoice(int id) {
    audio_->noteOff(id); audio_->removeVoice(id);
    auto it = std::find_if(strips_.begin(), strips_.end(), [id](VoiceStrip* s){ return s->voice_id == id; });
    if (it != strips_.end()) { VoiceStrip* s = *it; strips_.erase(it); s->group->hide(); scroll_voices_->remove(s->group); pool_strips_.push_back(s); relayoutStrips(); }
}

void HarmoniaApp::relayoutStrips() { int y = 5; for(auto* s : strips_) { s->group->position(s->group->x(), y); y += 79; } scroll_voices_->init_sizes(); }

VoiceStrip* HarmoniaApp::findStrip(int id) { for(auto* s:strips_) if(s->voice_id==id) return s; return nullptr; }

void HarmoniaApp::onIdle(void* data) {
    HarmoniaApp* app = (HarmoniaApp*)data; app->theory_->poll();
    auto snap = app->audio_->getSpectrumSnapshot(); auto rr = app->audio_->getRoughnessSnapshot(); auto obj = app->audio_->getAbstractObject();
    float tr = 0; for(auto& r : rr) tr += r.roughness; app->spectrum_->update(snap, app->audio_->getVoiceSnapshot(), std::min(1.f, tr));
    auto voices = app->audio_->getVoiceSnapshot();
    app->tonal_space_->setVoices(voices); app->tonal_space_->setAbstractObject(obj); app->tonal_space_->setRoughnessRecords(rr);
    if (obj.confidence > 0.2f) { app->out_object_->value(obj.chord_name.c_str()); app->box_object_label_->copy_label(obj.chord_name.c_str()); }
    for(auto& v : voices) {
        if (auto* s = app->findStrip(v.id)) {
            s->btn_on->value(v.note_on.load() ? 1 : 0);
            s->lbl_note->copy_label(noteName(v.pitch_class, v.octave, app->current_edo_).c_str());
            if (v.env_stage == Voice::EnvStage::IDLE && !v.note_on.load() && !s->manual) { app->removeVoice(v.id); break; }
        }
    }
}

void HarmoniaApp::onTonalSpaceClick(int pc, int oct) {
    // Standardize frequency calculation:
    // f = C4 * 2^((pc + (oct-4)*edo) / edo)
    double total_steps = (double)pc + (double)(oct - 4) * current_edo_;
    double freq = C4_HZ * std::pow(2.0, total_steps / current_edo_);

    auto voices = audio_->getVoiceSnapshot(); int eid = -1;
    for (auto& v : voices) {
        // Use a small epsilon in log space for robust comparison
        double diff_steps = std::abs(std::log2(v.frequency / freq) * current_edo_);
        if (diff_steps < 0.1) { eid = v.id; break; }
    }

    if (eid != -1) removeVoice(eid);
    else {
        int id = audio_->addVoice(60);
        if (id >= 0) {
            audio_->setVoiceFrequency(id, freq);
            audio_->noteOn(id);
            if (auto* s = findStrip(id)) s->manual = true;
        }
    }
    updateTheory();
}

std::pair<int,int> HarmoniaApp::highlightNodeForPC(int pc, int oct) {
    if (pc < 0) { tonal_space_->setHighlightedPC(-1); tonal_space_->setHighlightedNode(-1, -1); }
    else { if (oct != -100) tonal_space_->setHighlightedNode(pc, oct); else tonal_space_->setHighlightedPC(pc); }
    return {pc, oct};
}

void HarmoniaApp::updateFunctionalAnalysis() {
    auto voices = audio_->getVoiceSnapshot(); std::vector<int> pcs;
    for (auto& v : voices) if(v.active) pcs.push_back(v.pitch_class);
    if (pcs.empty()) { browser_function_->clear(); return; }
    theory_->queryAnalyzeChord(pcs, current_key_, current_edo_, [this](const FunctionalAnalysis& fa, const TonnetzTension& tt) {
        browser_function_->clear(); char buf[256]; snprintf(buf,256,"@bFUNCTION: %s  [T=%.2f]", fa.function.c_str(), tt.tension);
        browser_function_->add(buf); browser_function_->add(fa.function_reason.c_str());
        for(auto& t:fa.tendency_tones){ snprintf(buf,256,"  %s: %s", t.name.c_str(), t.tendency.c_str()); browser_function_->add(buf); browser_function_->data(browser_function_->size(), (void*)(intptr_t)(t.pc+1000)); }
    });
}
void HarmoniaApp::updateResolutionPaths() {
    auto obj = audio_->getAbstractObject(); int root = last_clicked_root_>=0 ? last_clicked_root_ : obj.root_pc;
    theory_->queryResolutionPaths(root, "maj", current_key_, current_edo_, [this](const std::vector<ResolutionPath>& paths){
        last_resolutions_ = paths; browser_resolutions_->clear();
        for(size_t i=0; i<paths.size(); i++){
            char buf[256]; snprintf(buf,256,"%s -> %s", paths[i].rule.c_str(), paths[i].target_label.c_str());
            browser_resolutions_->add(buf); browser_resolutions_->data(browser_resolutions_->size(), (void*)(intptr_t)(i+1));
        }
    });
}
void HarmoniaApp::updateCompletionSuggestions() {
    auto voices = audio_->getVoiceSnapshot(); std::vector<int> pcs; for(auto& v:voices) if(v.active) pcs.push_back(v.pitch_class);
    theory_->querySuggestCompletion(pcs, current_key_, current_edo_, [this](const std::vector<CompletionSuggestion>& ns){
        browser_completion_->clear(); for(auto& s:ns){ char buf[128]; snprintf(buf,128,"@b%s (score %.2f)", s.name.c_str(), s.score); browser_completion_->add(buf); browser_completion_->data(browser_completion_->size(), (void*)(intptr_t)(s.pc+1000)); }
    });
}
void HarmoniaApp::updatePsychoAnalysis() {
    auto voices = audio_->getVoiceSnapshot(); std::vector<int> pcs; for(auto& v:voices) if(v.active) pcs.push_back(v.pitch_class); if(pcs.empty()) return;
    theory_->queryPsychoacoustic(pcs, current_key_, current_edo_, (float)C4_HZ, 4, 70.0f, [this](const PsychoacousticAnalysis& pa){
        browser_psycho_->clear(); char buf[128]; snprintf(buf,128,"@bTENSION: %s (%.2f)", pa.perceptual_tension_label.c_str(), pa.perceptual_tension);
        browser_psycho_->add(buf); snprintf(buf,128,"  Consonance: %.2f", pa.level1.consonance_score); browser_psycho_->add(buf);
        snprintf(buf,128,"  Harmonicity: %.2f", pa.level2.harmonicity); browser_psycho_->add(buf);
        if(!pa.level2.virtual_pitch_name.empty()){
            snprintf(buf,128,"  Virtual Root: %s", pa.level2.virtual_pitch_name.c_str()); browser_psycho_->add(buf);
            browser_psycho_->data(browser_psycho_->size(), (void*)(intptr_t)(pa.level2.virtual_pitch_pc+1000));
        }
    });
}
void HarmoniaApp::updatePivotSearch() {
    theory_->queryPivotSearch(current_key_, ch_target_key_->value(), current_edo_, [this](const PivotSearchResult& res){
        last_pivot_res_ = res; browser_pivots_->clear();
        for(size_t i=0; i<res.all_pivots.size(); i++){
            char buf[256]; snprintf(buf,256,"[%zu] %s: %s -> %s", i, res.all_pivots[i].label.c_str(), res.all_pivots[i].roman_from.c_str(), res.all_pivots[i].roman_to.c_str());
            browser_pivots_->add(buf); browser_pivots_->data(browser_pivots_->size(), (void*)(intptr_t)(i+1));
        }
    });
}
void HarmoniaApp::updateTheory() { updateFunctionalAnalysis(); updateResolutionPaths(); updateCompletionSuggestions(); updatePsychoAnalysis(); }
void HarmoniaApp::playChord(const ChordKey& chord, bool sustain) {
    (void)sustain;
    for(const auto& n:chord.notes){
        int id = audio_->addVoice(60);
        if(id>=0){
            audio_->setVoiceFrequency(id, n.freq);
            audio_->noteOn(id);
            if(auto* s=findStrip(id)) s->manual=true;
        }
    }
    updateTheory();
}
void HarmoniaApp::releaseChord(const ChordKey& chord) { (void)chord; }

int HarmoniaApp::handle(int event) {
    if(event == FL_KEYDOWN){
        int k = Fl::event_key(); if(k >= '1' && k <= '9'){ int i = k-'1'; if(i<(int)instrument_keyboard_.size()) playChord(instrument_keyboard_[i], false); return 1; }
    }
    return 0;
}

void HarmoniaApp::cbPlayAll(Fl_Widget*, void* d) { for(auto* s:((HarmoniaApp*)d)->strips_) ((HarmoniaApp*)d)->audio_->noteOn(s->voice_id); }
void HarmoniaApp::cbStopAll(Fl_Widget*, void* d) { for(auto* s:((HarmoniaApp*)d)->strips_) ((HarmoniaApp*)d)->audio_->noteOff(s->voice_id); }
void HarmoniaApp::cbMasterVol(Fl_Widget* w, void* d) { ((HarmoniaApp*)d)->audio_->setMasterVolume((float)((Fl_Value_Slider*)w)->value()); }
void HarmoniaApp::cbKey(Fl_Widget* w, void* d) { ((HarmoniaApp*)d)->current_key_ = ((Fl_Choice*)w)->value(); ((HarmoniaApp*)d)->updateTheory(); }
void HarmoniaApp::cbEDO(Fl_Widget* w, void* d) {
    auto* a = (HarmoniaApp*)d; a->current_edo_ = (int)a->sp_edo_->value();
    a->tonal_space_->setEDO(a->current_edo_); a->audio_->setEDO(a->current_edo_); a->updateTheory();
    (void)w;
}
void HarmoniaApp::cbAnalyze(Fl_Widget*, void* d) { ((HarmoniaApp*)d)->updateFunctionalAnalysis(); }
void HarmoniaApp::cbResolve(Fl_Widget*, void* d) { ((HarmoniaApp*)d)->updateResolutionPaths(); }
void HarmoniaApp::cbSuggestCompletion(Fl_Widget*, void* d) { ((HarmoniaApp*)d)->updateCompletionSuggestions(); }
void HarmoniaApp::cbPsycho(Fl_Widget*, void* d) { ((HarmoniaApp*)d)->updatePsychoAnalysis(); }
void HarmoniaApp::cbFindPivots(Fl_Widget*, void* d) { ((HarmoniaApp*)d)->updatePivotSearch(); }
void HarmoniaApp::cbBrowserFunction(Fl_Widget* w, void* d) {
    auto* b=(Fl_Browser*)w; int l=b->value(); if(l<=0) return;
    int v=(int)(intptr_t)b->data(l); if(v>=1000) ((HarmoniaApp*)d)->highlightNodeForPC(v-1000);
}
void HarmoniaApp::cbBrowserPsycho(Fl_Widget* w, void* d) {
    auto* b=(Fl_Browser*)w; int l=b->value(); if(l<=0) return;
    int v=(int)(intptr_t)b->data(l); if(v>=1000) ((HarmoniaApp*)d)->highlightNodeForPC(v-1000);
}
void HarmoniaApp::cbBrowserResolutions(Fl_Widget* w, void* d) {
    auto* a=(HarmoniaApp*)d; auto* b=(Fl_Browser*)w; int l=b->value(); if(l<=0) return;
    int idx=(int)(intptr_t)b->data(l)-1; if(idx>=0 && idx<(int)a->last_resolutions_.size()){
        const auto& rp = a->last_resolutions_[idx]; a->highlightNodeForPC(rp.target_root);
        if(Fl::event_clicks()){ ChordKey ck; ck.name=rp.target_label; for(int pc:rp.target_pcs) ck.notes.push_back({pc, C4_HZ*std::pow(2.0, (double)pc/a->current_edo_), 4}); a->playChord(ck, true); }
    }
}
void HarmoniaApp::cbBrowserCompletion(Fl_Widget* w, void* d) {
    auto* b=(Fl_Browser*)w; int l=b->value(); if(l<=0) return;
    int v=(int)(intptr_t)b->data(l); if(v>=1000) ((HarmoniaApp*)d)->highlightNodeForPC(v-1000);
}
void HarmoniaApp::cbBrowserPivots(Fl_Widget* w, void* d) {
    auto* a=(HarmoniaApp*)d; auto* b=(Fl_Browser*)w; int l=b->value(); if(l<=0) return;
    int idx=(int)(intptr_t)b->data(l)-1; if(idx>=0 && idx<(int)a->last_pivot_res_.all_pivots.size()){
        const auto& p = a->last_pivot_res_.all_pivots[idx]; a->highlightNodeForPC(p.root);
        if(Fl::event_clicks()){ ChordKey ck; ck.name=p.label; for(int pc:p.pcs) ck.notes.push_back({pc, C4_HZ*std::pow(2.0, (double)pc/a->current_edo_), 4}); a->playChord(ck, true); }
    }
}
void HarmoniaApp::cbBrowserKeys(Fl_Widget* w, void* d) { if(Fl::event_clicks()){ auto* b=(Fl_Browser*)w; if(b->value()>0) ((HarmoniaApp*)d)->playChord(((HarmoniaApp*)d)->instrument_keyboard_[b->value()-1], true); } }
void HarmoniaApp::cbAddChordToKeys(Fl_Widget*, void* d) {
    auto* a = (HarmoniaApp*)d; auto voices = a->audio_->getVoiceSnapshot(); if(voices.empty()) return;
    ChordKey ck; ck.name = "Chord " + std::to_string(a->instrument_keyboard_.size()+1);
    for(auto& v:voices) if(v.note_on.load()){ ck.notes.push_back({v.pitch_class, v.frequency, v.octave}); }
    a->instrument_keyboard_.push_back(ck); a->browser_keys_->add(ck.name.c_str());
}
void HarmoniaApp::cbClearKeys(Fl_Widget*, void* d) { ((HarmoniaApp*)d)->instrument_keyboard_.clear(); ((HarmoniaApp*)d)->browser_keys_->clear(); }

int main() { Fl::lock(); HarmoniaApp app; app.run(); return 0; }
