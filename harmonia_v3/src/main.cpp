/*
 * harmonia/src/main.cpp
 * ─────────────────────────────────────────────────────────────────────────────
 *  HARMONIA — Multi-voice psychoacoustic music coder
 *
 *  Window layout (1280 × 820):
 *  ┌──────────────────────────────────────────────────────────────────────────┐
 *  │  [Harmonia]  Key: C▼  EDO: 12▼  Master:━━━━  [Add Voice] [Play All]    │ ← toolbar
 *  ├─────────────────┬────────────────────────────┬───────────────────────────┤
 *  │                 │                            │                           │
 *  │   VOICE LIST    │      TONNETZ (OpenGL)      │   THEORY PANEL            │
 *  │   + controls    │                            │   • Next chords           │
 *  │                 │  Interactive harmonic       │   • Completion            │
 *  │  [voice strip]  │  lattice. Click=add note   │   • Orbifold dist         │
 *  │  [voice strip]  │  Middle-drag=pan           │   • EDO analysis          │
 *  │  ...            │  Scroll=zoom               │                           │
 *  │                 │                            │                           │
 *  ├─────────────────┴────────────────────────────┴───────────────────────────┤
 *  │  SPECTRUM + ROUGHNESS bar            OBJECT DISPLAY                      │ ← status
 *  └──────────────────────────────────────────────────────────────────────────┘
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

#include "voice.h"
#include "audio_engine.h"
#include "tonnetz_widget.h"
#include "theory_bridge.h"

// ─────────────────────────────────────────────────────────────────────────────
//  FLTK colour helpers
// ─────────────────────────────────────────────────────────────────────────────
static Fl_Color flColor(float r, float g, float b) {
    return fl_rgb_color((uchar)(r*255),(uchar)(g*255),(uchar)(b*255));
}
static const Fl_Color COL_BG     = fl_rgb_color(18,20,28);
static const Fl_Color COL_PANEL  = fl_rgb_color(26,30,40);
static const Fl_Color COL_ACCENT = fl_rgb_color(200,160,40);
static const Fl_Color COL_TEXT   = fl_rgb_color(220,215,200);
static const Fl_Color COL_DIM    = fl_rgb_color(100,100,110);

// ─────────────────────────────────────────────────────────────────────────────
//  SPECTRUM DISPLAY — OpenGL widget showing partials per voice
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
            // color by frequency
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

        // Freq axis labels
        glColor3f(.4f,.4f,.5f);
        static const int FREQ_MARKS[] = {100,200,500,1000,2000,4000,8000,0};
        for (int i = 0; FREQ_MARKS[i]; i++) {
            float f = FREQ_MARKS[i];
            float x = (f / (SAMPLE_RATE/2)) * sw;
            glBegin(GL_LINES); glVertex2f(x,sh); glVertex2f(x,sh-4); glEnd();
            char lbl[16]; snprintf(lbl,sizeof(lbl),"%d",FREQ_MARKS[i]);
            gl_draw(lbl, (int)x-8, sh-6);
        }
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
    static constexpr int THEORY_PANEL_W = 300;
    static constexpr int STATUS_H = 120;
    static constexpr int TONNETZ_W = WIN_W - VOICE_PANEL_W - THEORY_PANEL_W;
    static constexpr int TONNETZ_H = WIN_H - TOOLBAR_H - STATUS_H;

    // ── FLTK widgets
    Fl_Double_Window* win_{nullptr};

    // Toolbar
    Fl_Button*       btn_add_{nullptr};
    Fl_Button*       btn_play_all_{nullptr};
    Fl_Button*       btn_stop_all_{nullptr};
    Fl_Choice*       ch_key_{nullptr};
    Fl_Spinner*      sp_edo_{nullptr};
    Fl_Value_Slider* sl_master_{nullptr};
    Fl_Output*       out_object_{nullptr};

    // Left panel – voice list
    Fl_Scroll*        scroll_voices_{nullptr};
    std::vector<VoiceStrip*> strips_;
    int               voice_panel_y_{0};  // next y for new strip

    // Centre – Tonnetz
    TonnetzWidget*    tonnetz_{nullptr};

    // Right – structural theory panel
    Fl_Browser*  browser_function_{nullptr};    // function + algebraic reasons
    Fl_Browser*  browser_resolutions_{nullptr}; // resolution paths with justification
    Fl_Browser*  browser_completion_{nullptr};  // completion suggestions + reasons
    Fl_Browser*  browser_psycho_{nullptr};      // psychoacoustic neural model
    Fl_Browser*  browser_edo_{nullptr};
    Fl_Button*   btn_analyze_{nullptr};
    Fl_Button*   btn_resolve_{nullptr};
    Fl_Button*   btn_complete_{nullptr};
    Fl_Button*   btn_psycho_{nullptr};
    Fl_Button*   btn_edo_{nullptr};
    Fl_Output*   out_orbifold_{nullptr};

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
    std::string last_chord_type_{"maj"};
    FunctionalAnalysis  last_func_;
    TonnetzTension      last_tension_;
    std::vector<ResolutionPath> last_resolutions_;

    // ── methods
    void buildUI();
    void setupCallbacks();

    // Voice management
    void addVoice(int midi_note = 60, TimbrePreset t = TimbrePreset::STRINGS);
    void removeVoice(int voice_id);
    VoiceStrip* findStrip(int voice_id);
    void relayoutStrips();

    // Structural theory updates
    void updateFunctionalAnalysis();
    void updateResolutionPaths();
    void updateCompletionSuggestions();
    void updatePsychoAnalysis();
    void updateEDOAnalysis();

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

    // Browser callbacks
    static void cbBrowserResolutions(Fl_Widget* w, void* d);
    static void cbBrowserCompletion(Fl_Widget* w, void* d);

    // Voice strip callbacks
    static void cbVoiceOn(Fl_Widget* w, void* d);
    static void cbVoiceAmp(Fl_Widget* w, void* d);
    static void cbVoiceDetune(Fl_Widget* w, void* d);
    static void cbVoiceTimbre(Fl_Widget* w, void* d);
    static void cbVoiceRemove(Fl_Widget* w, void* d);

    // Tonnetz click
    void onTonnetzClick(int pitch_class, int tx, int ty);

    std::string theoryDir_;
};

// ─────────────────────────────────────────────────────────────────────────────
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
    // ── start audio (ALSA)
    if (!audio_->init()) {
        fl_message("ALSA init failed — audio disabled");
    }

    // ── find theory server
    char exepath[512] = {};
    ssize_t n = readlink("/proc/self/exe", exepath, sizeof(exepath)-1);
    if (n > 0) {
        std::string dir(exepath, n);
        size_t pos = dir.rfind('/');
        theoryDir_ = dir.substr(0, pos) + "/../theory";
    } else {
        theoryDir_ = "./theory";
    }

    if (!theory_->start(theoryDir_)) {
        fl_message("Python theory server not started — theory features disabled");
    }

    // ── seed with a C major chord
    addVoice(60, TimbrePreset::STRINGS);  // C4
    addVoice(64, TimbrePreset::STRINGS);  // E4
    addVoice(67, TimbrePreset::STRINGS);  // G4

    // ── play them
    for (auto& s : strips_) audio_->noteOn(s->voice_id);

    // ── initial structural analysis
    updateFunctionalAnalysis();
    updatePsychoAnalysis();

    Fl::add_idle(onIdle, this);
    Fl::run();
}

// ─────────────────────────────────────────────────────────────────────────────
//  UI CONSTRUCTION
// ─────────────────────────────────────────────────────────────────────────────
void HarmoniaApp::buildUI() {
    Fl::scheme("gtk+");
    Fl::background(18,20,28);
    Fl::foreground(220,215,200);

    win_ = new Fl_Double_Window(WIN_W, WIN_H, "Harmonia — Psychoacoustic Music Coder");
    win_->color(COL_BG);

    // ── TOOLBAR ──────────────────────────────────────────────────────────────
    {
        Fl_Group* tb = new Fl_Group(0, 0, WIN_W, TOOLBAR_H);
        tb->box(FL_FLAT_BOX); tb->color(fl_rgb_color(14,16,22));

        int x = 8, y = 6, bh = 26;

        Fl_Box* title = new Fl_Box(x,y,100,bh,"HARMONIA");
        title->labelcolor(COL_ACCENT); title->labelfont(FL_HELVETICA_BOLD);
        title->labelsize(14); title->box(FL_NO_BOX);
        x += 108;

        new Fl_Box(x,y,30,bh,"Key:"); x+=30;
        ch_key_ = new Fl_Choice(x,y,70,bh);
        for (int i=0;i<12;i++) ch_key_->add(NOTE_NAMES[i]);
        ch_key_->value(0); x += 76;

        new Fl_Box(x,y,30,bh,"EDO:"); x+=32;
        sp_edo_ = new Fl_Spinner(x,y,60,bh);
        sp_edo_->minimum(5); sp_edo_->maximum(72); sp_edo_->value(12);
        x += 68;

        new Fl_Box(x,y,50,bh,"Master:"); x+=52;
        sl_master_ = new Fl_Value_Slider(x,y,100,bh);
        sl_master_->type(FL_HORIZONTAL); sl_master_->minimum(0); sl_master_->maximum(1);
        sl_master_->value(0.5); sl_master_->color(fl_rgb_color(40,44,58));
        x += 110;

        btn_add_ = new Fl_Button(x,y,90,bh,"+ Add Voice");
        btn_add_->color(fl_rgb_color(30,80,50)); btn_add_->labelcolor(fl_rgb_color(120,220,100));
        x += 98;

        btn_play_all_ = new Fl_Button(x,y,70,bh,"▶ All");
        btn_play_all_->color(fl_rgb_color(30,60,100)); btn_play_all_->labelcolor(fl_rgb_color(100,180,255));
        x += 76;

        btn_stop_all_ = new Fl_Button(x,y,70,bh,"■ Stop");
        btn_stop_all_->color(fl_rgb_color(80,30,30)); btn_stop_all_->labelcolor(fl_rgb_color(255,120,100));
        x += 76;

        out_object_ = new Fl_Output(x,y, WIN_W-x-8, bh);
        out_object_->box(FL_FLAT_BOX); out_object_->color(fl_rgb_color(20,24,32));
        out_object_->textcolor(COL_ACCENT); out_object_->textsize(13);
        out_object_->value("Object: —");

        tb->end();
    }

    int content_y = TOOLBAR_H;
    int content_h = WIN_H - TOOLBAR_H - STATUS_H;

    // ── VOICE PANEL (left) ────────────────────────────────────────────────────
    {
        Fl_Group* panel = new Fl_Group(0, content_y, VOICE_PANEL_W, content_h);
        panel->box(FL_FLAT_BOX); panel->color(fl_rgb_color(22,26,35));

        Fl_Box* hdr = new Fl_Box(0,content_y, VOICE_PANEL_W,22,"VOICES");
        hdr->box(FL_FLAT_BOX); hdr->color(fl_rgb_color(30,34,46));
        hdr->labelcolor(COL_DIM); hdr->labelfont(FL_HELVETICA); hdr->labelsize(10);

        scroll_voices_ = new Fl_Scroll(0, content_y+22, VOICE_PANEL_W, content_h-22);
        scroll_voices_->box(FL_FLAT_BOX); scroll_voices_->color(fl_rgb_color(22,26,35));
        scroll_voices_->type(Fl_Scroll::VERTICAL_ALWAYS);
        scroll_voices_->end();

        voice_panel_y_ = content_y + 22;
        panel->end();
    }

    // ── TONNETZ (centre) ──────────────────────────────────────────────────────
    {
        tonnetz_ = new TonnetzWidget(VOICE_PANEL_W, content_y, TONNETZ_W, content_h);
    }

    // ── THEORY PANEL (right) — STRUCTURAL EDITION ────────────────────────────
    {
        int tx = VOICE_PANEL_W + TONNETZ_W;
        Fl_Group* panel = new Fl_Group(tx, content_y, THEORY_PANEL_W, content_h);
        panel->box(FL_FLAT_BOX); panel->color(fl_rgb_color(22,26,35));

        int py = content_y+4, pw = THEORY_PANEL_W-8, px = tx+4;

        // ─ Functional Analysis
        auto* h1 = new Fl_Box(px,py,pw,16,"HARMONIC FUNCTION");
        h1->labelfont(FL_HELVETICA_BOLD); h1->labelsize(9);
        h1->labelcolor(COL_ACCENT); h1->box(FL_NO_BOX); py+=18;

        browser_function_ = new Fl_Browser(px,py,pw,130);
        browser_function_->textsize(9); browser_function_->color(fl_rgb_color(18,22,30));
        browser_function_->textcolor(COL_TEXT); py+=132;

        // Set callback inline so it cannot be missed by setupCallbacks()
        btn_analyze_ = new Fl_Button(px,py,pw,20,"Analyse chord structure");
        btn_analyze_->labelsize(9); btn_analyze_->color(fl_rgb_color(30,50,80));
        btn_analyze_->callback(cbAnalyze, this); py+=24;

        // ─ Resolution Paths
        auto* h2 = new Fl_Box(px,py,pw,16,"RESOLUTION PATHS");
        h2->labelfont(FL_HELVETICA_BOLD); h2->labelsize(9);
        h2->labelcolor(fl_rgb_color(200,120,60)); h2->box(FL_NO_BOX); py+=18;

        browser_resolutions_ = new Fl_Browser(px,py,pw,140);
        browser_resolutions_->textsize(9); browser_resolutions_->color(fl_rgb_color(18,22,30));
        browser_resolutions_->textcolor(fl_rgb_color(220,200,170)); py+=142;

        btn_resolve_ = new Fl_Button(px,py,pw,20,"Show resolution paths");
        btn_resolve_->labelsize(9); btn_resolve_->color(fl_rgb_color(50,30,20));
        btn_resolve_->callback(cbResolve, this); py+=24;

        // ─ Voice Completion
        auto* h3 = new Fl_Box(px,py,pw,16,"COMPLETION SUGGESTIONS");
        h3->labelfont(FL_HELVETICA_BOLD); h3->labelsize(9);
        h3->labelcolor(fl_rgb_color(100,200,140)); h3->box(FL_NO_BOX); py+=18;

        browser_completion_ = new Fl_Browser(px,py,pw,100);
        browser_completion_->textsize(9); browser_completion_->color(fl_rgb_color(18,22,30));
        browser_completion_->textcolor(fl_rgb_color(160,220,180)); py+=102;

        btn_complete_ = new Fl_Button(px,py,pw,20,"Suggest next voice");
        btn_complete_->labelsize(9); btn_complete_->color(fl_rgb_color(20,50,30));
        btn_complete_->callback(cbSuggestCompletion, this); py+=24;

        // ─ Psychoacoustic Neural Model
        auto* h5 = new Fl_Box(px,py,pw,16,"AUDITORY CORTEX MODEL");
        h5->labelfont(FL_HELVETICA_BOLD); h5->labelsize(9);
        h5->labelcolor(fl_rgb_color(200,100,200)); h5->box(FL_NO_BOX); py+=18;

        browser_psycho_ = new Fl_Browser(px,py,pw,120);
        browser_psycho_->textsize(9); browser_psycho_->color(fl_rgb_color(18,22,30));
        browser_psycho_->textcolor(fl_rgb_color(220,180,220)); py+=122;

        btn_psycho_ = new Fl_Button(px,py,pw,20,"Neural psychoanalysis");
        btn_psycho_->labelsize(9); btn_psycho_->color(fl_rgb_color(60,20,60));
        btn_psycho_->callback(cbPsycho, this); py+=24;

        // ─ Orbifold distance
        out_orbifold_ = new Fl_Output(px,py,pw,30);
        out_orbifold_->textsize(9); out_orbifold_->color(fl_rgb_color(18,22,30));
        out_orbifold_->textcolor(fl_rgb_color(150,200,150)); py+=34;

        // ─ EDO analysis
        auto* h4 = new Fl_Box(px,py,pw,16,"EDO ANALYSIS");
        h4->labelfont(FL_HELVETICA_BOLD); h4->labelsize(9);
        h4->labelcolor(fl_rgb_color(140,100,220)); h4->box(FL_NO_BOX); py+=18;

        browser_edo_ = new Fl_Browser(px,py,pw, content_h-(py-content_y)-28);
        browser_edo_->textsize(9); browser_edo_->color(fl_rgb_color(18,22,30));
        browser_edo_->textcolor(fl_rgb_color(180,160,220));
        py = content_y+content_h-26;

        btn_edo_ = new Fl_Button(px,py,pw,22,"Analyse EDO");
        btn_edo_->labelsize(9); btn_edo_->color(fl_rgb_color(40,25,70));
        btn_edo_->callback(cbAnalyseEDO, this);

        panel->end();
    }

    // ── STATUS / SPECTRUM BAR (bottom) ────────────────────────────────────────
    {
        int sy = WIN_H - STATUS_H;
        Fl_Group* sb = new Fl_Group(0,sy, WIN_W, STATUS_H);
        sb->box(FL_FLAT_BOX); sb->color(fl_rgb_color(14,16,22));

        spectrum_ = new SpectrumWidget(0, sy, WIN_W - 250, STATUS_H);

        box_object_label_ = new Fl_Box(WIN_W-248, sy, 248, STATUS_H, "");
        box_object_label_->box(FL_FLAT_BOX); box_object_label_->color(fl_rgb_color(18,22,30));
        box_object_label_->labelcolor(COL_ACCENT); box_object_label_->labelfont(FL_HELVETICA_BOLD);
        box_object_label_->labelsize(16); box_object_label_->align(FL_ALIGN_CENTER|FL_ALIGN_INSIDE);

        sb->end();
    }

    win_->end();
    win_->show();
}

// ─────────────────────────────────────────────────────────────────────────────
void HarmoniaApp::setupCallbacks() {
    btn_add_->callback(cbAddVoice, this);
    btn_play_all_->callback(cbPlayAll, this);
    btn_stop_all_->callback(cbStopAll, this);
    sl_master_->callback(cbMasterVol, this);
    ch_key_->callback(cbKey, this);
    sp_edo_->callback(cbEDO, this);
    browser_resolutions_->callback(cbBrowserResolutions, this);
    browser_completion_->callback(cbBrowserCompletion, this);
    tonnetz_->setNodeClickCallback([this](int pc, int tx, int ty){
        onTonnetzClick(pc, tx, ty);
    });
}

// ─────────────────────────────────────────────────────────────────────────────
//  VOICE MANAGEMENT
// ─────────────────────────────────────────────────────────────────────────────
void HarmoniaApp::addVoice(int midi_note, TimbrePreset t) {
    int id = audio_->addVoice(midi_note, t);
    if (id < 0) { fl_message("Max voices reached"); return; }

    static const int STRIP_H = 76;
    int sw = VOICE_PANEL_W - 12;

    scroll_voices_->begin();
    VoiceStrip* strip = new VoiceStrip();
    strip->voice_id = id;

    // Initial position doesn't matter much as relayoutStrips() will fix it
    int gx = scroll_voices_->x() + 4;
    int gy = scroll_voices_->y() + 4;
    strip->group = new Fl_Group(gx, gy, sw, STRIP_H);
    strip->group->box(FL_FLAT_BOX);
    strip->group->color(fl_rgb_color(28,32,44));

    // ── Row 1: colour bar + note name + mute + remove ──────────────────
    // Thin colour accent bar on left edge
    Fl_Box* colorbar = new Fl_Box(gx, gy, 4, STRIP_H, "");
    colorbar->box(FL_FLAT_BOX);
    colorbar->color(fl_rgb_color(80,120,200));

    // Note name — large and readable
    strip->lbl_note = new Fl_Box(gx+6, gy+4, 60, 26, "—");
    strip->lbl_note->box(FL_FLAT_BOX);
    strip->lbl_note->color(fl_rgb_color(18,22,32));
    strip->lbl_note->labelcolor(fl_rgb_color(200,220,255));
    strip->lbl_note->labelfont(FL_HELVETICA_BOLD);
    strip->lbl_note->labelsize(15);
    strip->lbl_note->align(FL_ALIGN_CENTER | FL_ALIGN_INSIDE);

    // Roughness indicator — compact pill next to note name
    strip->roughness_box = new Fl_Box(gx+70, gy+4, 46, 13, "R:0.00");
    strip->roughness_box->box(FL_FLAT_BOX);
    strip->roughness_box->color(fl_rgb_color(20,26,36));
    strip->roughness_box->labelcolor(fl_rgb_color(160,200,140));
    strip->roughness_box->labelsize(8);

    // Mute/unmute (light button, green when active)
    strip->btn_on = new Fl_Light_Button(gx+70, gy+18, 46, 14, "PLAY");
    strip->btn_on->labelsize(8);
    strip->btn_on->value(0);
    strip->btn_on->selection_color(fl_rgb_color(40,180,80));
    strip->btn_on->user_data((void*)(intptr_t)id);
    strip->btn_on->callback(cbVoiceOn, this);

    // Remove button — right-aligned, small
    strip->btn_remove = new Fl_Button(gx+sw-22, gy+4, 18, 18, "×");
    strip->btn_remove->labelsize(12);
    strip->btn_remove->box(FL_FLAT_BOX);
    strip->btn_remove->color(fl_rgb_color(60,18,18));
    strip->btn_remove->labelcolor(fl_rgb_color(220,80,80));
    strip->btn_remove->user_data((void*)(intptr_t)id);
    strip->btn_remove->callback(cbVoiceRemove, this);

    // Timbre dropdown — spans top right area
    strip->ch_timbre = new Fl_Choice(gx+118, gy+4, sw-142, 18);
    strip->ch_timbre->textsize(9);
    strip->ch_timbre->add("Sine|Saw|Square|Strings|Brass|Flute");
    strip->ch_timbre->value(3);
    strip->ch_timbre->user_data((void*)(intptr_t)id);
    strip->ch_timbre->callback(cbVoiceTimbre, this);

    // ── Row 2: Amp slider ───────────────────────────────────────────────
    auto* lbl_a = new Fl_Box(gx+6, gy+36, 22, 12, "Vol");
    lbl_a->labelsize(8); lbl_a->labelcolor(fl_rgb_color(140,150,170)); lbl_a->box(FL_NO_BOX);
    strip->sl_amp = new Fl_Value_Slider(gx+30, gy+36, sw-34, 13);
    strip->sl_amp->type(FL_HORIZONTAL);
    strip->sl_amp->minimum(0); strip->sl_amp->maximum(1); strip->sl_amp->value(0.6);
    strip->sl_amp->textsize(8);
    strip->sl_amp->color(fl_rgb_color(22,26,38));
    strip->sl_amp->selection_color(fl_rgb_color(40,100,180));
    strip->sl_amp->user_data((void*)(intptr_t)id);
    strip->sl_amp->callback(cbVoiceAmp, this);

    // ── Row 3: Detune slider ────────────────────────────────────────────
    auto* lbl_d = new Fl_Box(gx+6, gy+56, 22, 12, "Det");
    lbl_d->labelsize(8); lbl_d->labelcolor(fl_rgb_color(140,150,170)); lbl_d->box(FL_NO_BOX);
    strip->sl_detune = new Fl_Value_Slider(gx+30, gy+56, sw-34, 13);
    strip->sl_detune->type(FL_HORIZONTAL);
    strip->sl_detune->minimum(-50); strip->sl_detune->maximum(50); strip->sl_detune->value(0);
    strip->sl_detune->textsize(8);
    strip->sl_detune->color(fl_rgb_color(22,26,38));
    strip->sl_detune->selection_color(fl_rgb_color(100,60,160));
    strip->sl_detune->user_data((void*)(intptr_t)id);
    strip->sl_detune->callback(cbVoiceDetune, this);

    strip->group->end();
    scroll_voices_->end();

    // Update note label and voice colour
    auto voices = audio_->getVoiceSnapshot();
    for (auto& v : voices) {
        if (v.id == id) {
            strip->color[0] = v.color[0]; strip->color[1] = v.color[1]; strip->color[2] = v.color[2];
            char buf[16];
            snprintf(buf,sizeof(buf),"%s%d", NOTE_NAMES[v.pitch_class], v.octave);
            strip->lbl_note->copy_label(buf);
            // Tint: group bg slightly, colour bar fully
            Fl_Color accent = fl_rgb_color(
                (uchar)(v.color[0]*200), (uchar)(v.color[1]*200), (uchar)(v.color[2]*200));
            strip->group->color(flColor(v.color[0]*0.12f, v.color[1]*0.12f, v.color[2]*0.12f));
            // The first child of group is the colour bar
            if (strip->group->children() > 0)
                strip->group->child(0)->color(accent);
        }
    }

    strips_.push_back(strip);
    relayoutStrips();
    win_->redraw();
}

void HarmoniaApp::removeVoice(int voice_id) {
    audio_->noteOff(voice_id);
    audio_->removeVoice(voice_id);

    auto it = std::find_if(strips_.begin(), strips_.end(),
        [voice_id](VoiceStrip* s){ return s->voice_id == voice_id; });
    if (it == strips_.end()) return;

    VoiceStrip* strip = *it;
    strips_.erase(it);   // remove from our list first

    // Hide the group immediately so it disappears on the next redraw.
    // Then detach from scroll's child list and queue the actual deletion.
    // This sequence is safe even when called from the × button's own callback
    // because Fl::delete_widget() defers destruction to after event dispatch.
    strip->group->hide();
    scroll_voices_->remove(strip->group);
    Fl::delete_widget(strip->group);  // FLTK now owns and will free
    strip->group = nullptr;
    delete strip;                     // free our wrapper (no FLTK widgets inside)

    relayoutStrips();
    win_->redraw();
}

void HarmoniaApp::relayoutStrips() {
    static const int STRIP_H = 76;
    // Save scroll position
    int sx = scroll_voices_->xposition();
    int sy = scroll_voices_->yposition();
    // Reset scroll to 0 to position widgets correctly in absolute coordinates
    scroll_voices_->scroll_to(0, 0);

    int y = scroll_voices_->y() + 4;
    for (auto* s : strips_) {
        // Use position() not resize() — resize cascades to children and can
        // reset slider values; we only need to move the group vertically.
        s->group->position(s->group->x(), y);
        y += STRIP_H + 3;
    }
    // Recompute the scroll's virtual canvas so it knows the new total height.
    scroll_voices_->init_sizes();
    // Restore scroll position
    scroll_voices_->scroll_to(sx, sy);
    scroll_voices_->redraw();
}

VoiceStrip* HarmoniaApp::findStrip(int voice_id) {
    for (auto* s : strips_) if (s->voice_id == voice_id) return s;
    return nullptr;
}

// ─────────────────────────────────────────────────────────────────────────────
//  IDLE — update UI from audio thread results
// ─────────────────────────────────────────────────────────────────────────────
void HarmoniaApp::onIdle(void* data) {
    HarmoniaApp* app = (HarmoniaApp*)data;

    // Throttle to ~30 fps. NEVER call Fl::check() inside an idle handler —
    // it re-enters the FLTK event loop recursively and breaks all click delivery.
    static double last_t = 0.0;
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    double now = ts.tv_sec + ts.tv_nsec * 1e-9;
    if (now - last_t < 0.033) return;
    last_t = now;

    // Poll theory bridge (fires queued callbacks on the main thread, safe)
    app->theory_->poll();

    // Grab audio snapshots (mutex-protected, no FLTK calls inside lock)
    auto snap   = app->audio_->getSpectrumSnapshot();
    auto voices = app->audio_->getVoiceSnapshot();
    auto rr     = app->audio_->getRoughnessSnapshot();
    auto obj    = app->audio_->getAbstractObject();

    // Spectrum widget
    float total_rough = 0.f;
    for (auto& r : rr) total_rough += r.roughness;
    app->spectrum_->update(snap, voices, std::min(1.f, total_rough));

    // Tonnetz
    app->tonnetz_->setVoices(voices);
    app->tonnetz_->setAbstractObject(obj);
    app->tonnetz_->setRoughnessRecords(rr);

    // Abstract object label
    if (obj.confidence > 0.2f) {
        char buf[256];
        snprintf(buf, sizeof(buf), "%s  [vp=%.0fHz]  R=%.2f",
                 obj.chord_name.c_str(), obj.virtual_pitch_hz, obj.roughness_total);
        app->out_object_->value(buf);
        app->box_object_label_->copy_label(obj.chord_name.c_str());
    }

    // Roughness indicators on voice strips (only redraw if changed)
    std::map<int,float> voice_rough;
    for (auto& r : rr) {
        voice_rough[r.voice_a] += r.roughness;
        voice_rough[r.voice_b] += r.roughness;
    }
    bool changed = false;
    for (auto* s : app->strips_) {
        float rv = voice_rough.count(s->voice_id) ? voice_rough[s->voice_id] : 0.f;
        if (std::fabs(rv - s->roughness_val) < 0.005f) continue;
        s->roughness_val = rv;
        float rn = std::min(1.f, rv);
        uchar rc = (uchar)(std::min(1.f, rn * 2.0f) * 200);
        uchar gc = (uchar)(std::max(0.f, 1.f - rn * 1.5f) * 160);
        char rbuf[16]; snprintf(rbuf, sizeof(rbuf), "R %.2f", rv);
        s->roughness_box->copy_label(rbuf);
        s->roughness_box->color(fl_rgb_color(rc > 20 ? rc : 20, gc > 20 ? gc : 20, 36));
        changed = true;
    }
    if (changed) app->scroll_voices_->redraw();

    // Sync PLAY button visual state with actual note_on state
    for (auto& v : voices) {
        if (auto* s = app->findStrip(v.id)) {
            int want = v.note_on.load() ? 1 : 0;
            if (s->btn_on->value() != want) s->btn_on->value(want);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  TONNETZ CLICK — add voice at clicked pitch class
// ─────────────────────────────────────────────────────────────────────────────
void HarmoniaApp::onTonnetzClick(int pitch_class, int /*tx*/, int /*ty*/) {
    // pitch_class is 0..edo-1. Map to nearest MIDI note + detune.
    double cents = (double)pitch_class * 1200.0 / current_edo_;
    int closest_midi = 60 + (int)std::round(cents / 100.0);
    float detune = (float)(cents - (closest_midi - 60) * 100.0);

    last_clicked_root_ = pitch_class;
    auto voices = audio_->getVoiceSnapshot();

    int existing_id = -1;
    for (auto& v : voices) {
        if (v.pitch_class == pitch_class) {
            existing_id = v.id; break;
        }
    }

    if (existing_id != -1) {
        removeVoice(existing_id);
    } else {
        addVoice(closest_midi);
        if (!strips_.empty()) {
            int newid = strips_.back()->voice_id;
            audio_->setVoiceDetune(newid, detune);
            if (auto* s = findStrip(newid)) {
                s->sl_detune->value(detune);
                s->btn_on->value(1);
            }
            audio_->noteOn(newid);
        }
    }

    updateFunctionalAnalysis();
    updateCompletionSuggestions();
    updatePsychoAnalysis();
}

// ─────────────────────────────────────────────────────────────────────────────
//  STRUCTURAL THEORY UPDATES
// ─────────────────────────────────────────────────────────────────────────────

void HarmoniaApp::updateFunctionalAnalysis() {
    auto voices = audio_->getVoiceSnapshot();
    std::vector<int> pcs;
    for (auto& v : voices) pcs.push_back(v.pitch_class);
    if (pcs.empty()) return;

    theory_->queryAnalyzeChord(pcs, current_key_, current_edo_,
        [this](const FunctionalAnalysis& fa, const TonnetzTension& tt) {
            last_func_    = fa;
            last_tension_ = tt;

            browser_function_->clear();

            // ── Function + reason
            char buf[256];
            // function badge
            snprintf(buf,sizeof(buf),"@b%s  [tension=%d]  dist=%.2f",
                     fa.function.c_str(), fa.tension_level, tt.distance);
            browser_function_->add(buf);

            // truncated reason (split into lines)
            std::string reason = fa.function_reason;
            while (reason.size() > 50) {
                size_t split = reason.rfind(' ', 50);
                if (split == std::string::npos) split = 50;
                browser_function_->add(("  " + reason.substr(0,split)).c_str());
                reason = reason.substr(split+1);
            }
            if (!reason.empty()) browser_function_->add(("  " + reason).c_str());

            // ── Set class
            // We get this from the same analyze_chord call via raw query
            browser_function_->add("  —");

            // ── Tritone
            if (fa.has_tritone && fa.tritone.size() >= 2) {
                snprintf(buf,sizeof(buf),"  Tritone: %s-%s  (resolves inward)",
                         fa.tritone_names.size()>=2
                             ? fa.tritone_names[0].c_str() : "?",
                         fa.tritone_names.size()>=2
                             ? fa.tritone_names[1].c_str() : "?");
                browser_function_->add(buf);
            }

            // ── Tendency tones
            if (!fa.tendency_tones.empty()) {
                browser_function_->add("  Tendency tones:");
                for (auto& t : fa.tendency_tones) {
                    snprintf(buf,sizeof(buf),"   %s (%s) -> %s [%s]",
                             t.name.c_str(), t.role.c_str(),
                             t.tendency.c_str(), t.force.c_str());
                    browser_function_->add(buf);
                }
            }

            // ── Tension bar (text)
            snprintf(buf,sizeof(buf),"  Tension: %s (%.0f%%)",
                     tt.tension_label.c_str(), tt.tension*100.f);
            browser_function_->add(buf);

            browser_function_->redraw();
        });

    // Also request set-class info for ICV
    std::string pcs_json = "[";
    for (size_t i=0;i<pcs.size();i++) {
        pcs_json += std::to_string(pcs[i]);
        if (i+1<pcs.size()) pcs_json += ",";
    }
    pcs_json += "]";
    theory_->queryRaw(
        "{\"cmd\":\"analyze_chord\",\"tag\":\"analyze_chord_icv\",\"pcs\":" + pcs_json +
        ",\"key\":" + std::to_string(current_key_) + ",\"edo\":" + std::to_string(current_edo_) + "}",
        "analyze_chord_icv",  // separate tag so it doesn't conflict
        [this](const std::string& resp) {
            // Extract ICV and forte number
            // find ic_description
            size_t p = resp.find("\"ic_description\"");
            if (p == std::string::npos) return;
            size_t colon = resp.find(':', p + 16);
            if (colon == std::string::npos) return;
            size_t qs = resp.find('"', colon + 1);
            if (qs == std::string::npos) return;
            size_t qe = resp.find('"', qs+1);
            if (qe == std::string::npos) return;
            std::string ic_desc = resp.substr(qs+1, qe-qs-1);
            std::string forte   = TheoryBridge::jStr(resp,"forte","?");
            std::string cname   = TheoryBridge::jStr(resp,"common_name","");
            char buf[128];
            snprintf(buf,sizeof(buf),"  Forte: %s  (%s)", forte.c_str(), cname.c_str());
            if (browser_function_) browser_function_->insert(3, buf);
            snprintf(buf,sizeof(buf),"  ICV:   [%s]", ic_desc.c_str());
            if (browser_function_) browser_function_->insert(4, buf);
            browser_function_->redraw();
        });
}

void HarmoniaApp::updateResolutionPaths() {
    auto obj = audio_->getAbstractObject();
    if (obj.root_pc < 0 && last_clicked_root_ < 0) return;
    int root = last_clicked_root_ >= 0 ? last_clicked_root_ : obj.root_pc;
    std::string qual = obj.quality;
    if (qual == "" || qual == "?") qual = "maj";

    theory_->queryResolutionPaths(root, qual, current_key_, current_edo_,
        [this](const std::vector<ResolutionPath>& paths) {
            last_resolutions_ = paths;
            browser_resolutions_->clear();

            // Build Tonnetz progression path for visualization
            int path_idx = 0;
            for (auto& rp : paths) {
                char buf[256];
                // Rule class badge
                const char* col = (rp.rule_class == "necessity") ? "@C1" :
                                  (rp.rule_class == "minimal_motion") ? "@C4" : "@C6";
                snprintf(buf,sizeof(buf),"%s[P%d] %s -> %s (%s, VL=%d)",
                         col, rp.priority,
                         rp.rule.c_str(),
                         rp.target_label.c_str(),
                         rp.voice_leading.smoothness.c_str(),
                         rp.voice_leading.distance);
                browser_resolutions_->add(buf);
                browser_resolutions_->data(browser_resolutions_->size(), (void*)(intptr_t)path_idx);

                // Explanation (word-wrapped)
                std::string ex = rp.explanation;
                if (ex.size() > 48) ex = ex.substr(0,48) + "...";
                browser_resolutions_->add(("  " + ex).c_str());

                // Voice motions
                for (auto& m : rp.voice_leading.motions) {
                    if (m.semitones == 0) continue;
                    snprintf(buf,sizeof(buf),"   %s->%s (%+d)",
                             m.from_name.c_str(), m.to_name.c_str(), m.semitones);
                    browser_resolutions_->add(buf);
                    browser_resolutions_->data(browser_resolutions_->size(), (void*)(intptr_t)path_idx);
                }
                browser_resolutions_->add(" ");
                path_idx++;
            }
            browser_resolutions_->redraw();
        });
}

void HarmoniaApp::updateCompletionSuggestions() {
    auto voices = audio_->getVoiceSnapshot();
    std::vector<int> pcs;
    for (auto& v : voices) pcs.push_back(v.pitch_class);

    theory_->querySuggestCompletion(pcs, current_key_, current_edo_,
        [this](const std::vector<CompletionSuggestion>& ns) {
            browser_completion_->clear();
            for (auto& s : ns) {
                char buf[128];
                // Score bar
                int bars = (int)(s.score * 8);
                std::string bar(bars, '|');
                snprintf(buf,sizeof(buf),"@b%-3s  %s  R+%.2f  %s",
                         s.name.c_str(), bar.c_str(),
                         s.roughness_delta,
                         s.in_key ? "(diatonic)" : "(chromatic)");
                browser_completion_->add(buf);
                browser_completion_->data(browser_completion_->size(), (void*)(intptr_t)s.pc);

                // Structural reasons
                for (auto& r : s.structural_reasons) {
                    browser_completion_->add(("   " + r).c_str());
                    browser_completion_->data(browser_completion_->size(), (void*)(intptr_t)s.pc);
                }
                if (!s.structural_reasons.empty())
                    browser_completion_->add(" ");
            }
            browser_completion_->redraw();
        });
}

void HarmoniaApp::updatePsychoAnalysis() {
    auto voices = audio_->getVoiceSnapshot();
    std::vector<int> pcs;
    for (auto& v : voices) pcs.push_back(v.pitch_class);
    if (pcs.empty()) return;

    theory_->queryPsychoacoustic(pcs, current_key_, 261.63f, 4, 70.0f,
        [this](const PsychoacousticAnalysis& pa) {
            browser_psycho_->clear();
            char buf[256];

            snprintf(buf, sizeof(buf), "@bTENSION: %s (%.2f)",
                     pa.perceptual_tension_label.c_str(), pa.perceptual_tension);
            browser_psycho_->add(buf);

            snprintf(buf, sizeof(buf), "  Breakdown: R:%.2f H:%.2f T:%.2f",
                     pa.tension_breakdown.roughness_component,
                     pa.tension_breakdown.harmonicity_component,
                     pa.tension_breakdown.tonal_component);
            browser_psycho_->add(buf);

            browser_psycho_->add("@bL1: PERIPHERAL");
            snprintf(buf, sizeof(buf), "  Consonance: %.2f", pa.level1.consonance_score);
            browser_psycho_->add(buf);
            if (!pa.level1.masked_tone_names.empty()) {
                std::string mt = "  Masked: ";
                for (auto& s : pa.level1.masked_tone_names) mt += s + " ";
                browser_psycho_->add(mt.c_str());
            }

            browser_psycho_->add("@bL2: BRAINSTEM");
            snprintf(buf, sizeof(buf), "  Virtual Root: %s (%.2fHz)",
                     pa.level2.virtual_pitch_name.c_str(), pa.level2.virtual_pitch_hz);
            browser_psycho_->add(buf);
            snprintf(buf, sizeof(buf), "  Harmonicity: %s (%.2f)",
                     pa.level2.harmonicity_label.c_str(), pa.level2.harmonicity);
            browser_psycho_->add(buf);

            browser_psycho_->add("@bL3: CORTICAL");
            snprintf(buf, sizeof(buf), "  Tonal Stability: %s (%.2f)",
                     pa.level3.kk_stability_label.c_str(), pa.level3.kk_tonal_stability);
            browser_psycho_->add(buf);
            snprintf(buf, sizeof(buf), "  Centroid: %s (%.0fHz)",
                     pa.level3.spectral_centroid_name.c_str(), pa.level3.spectral_centroid_hz);
            browser_psycho_->add(buf);

            browser_psycho_->redraw();
        });
}

void HarmoniaApp::updateEDOAnalysis() {
    theory_->queryEDOAnalysis(current_edo_,
        [this](const std::string& raw) {
            browser_edo_->clear();
            char buf[128];
            snprintf(buf,sizeof(buf),"EDO-%d  step=%.2f¢",
                     current_edo_, TheoryBridge::jFloat(raw,"step_cents"));
            browser_edo_->add(buf);
            snprintf(buf,sizeof(buf),"5-limit score: %.1f%%",
                     TheoryBridge::jFloat(raw,"consonance_score"));
            browser_edo_->add(buf);
            snprintf(buf,sizeof(buf),"7-limit score: %.1f%%",
                     TheoryBridge::jFloat(raw,"seven_limit_score"));
            browser_edo_->add(buf);
            browser_edo_->add("  — Prime errors from JI —");

            static const struct { const char* prime; const char* name; } PRIMES[] = {
                {"3","P5/P4"},{"5","M3/m6"},{"7","h7"},{"11","11th"},{"13","13th"},{nullptr,nullptr}
            };
            for (int i=0; PRIMES[i].prime; i++) {
                size_t p = raw.find(std::string("\"")+PRIMES[i].prime+"\"");
                if (p==std::string::npos) continue;
                size_t ep = raw.find("\"error_cents\"",p);
                if (ep==std::string::npos) continue;
                size_t colon = raw.find(':', ep + 13);
                if (colon == std::string::npos) continue;
                size_t val_start = raw.find_first_of("0123456789-.", colon + 1);
                if (val_start == std::string::npos) continue;
                size_t ee = raw.find_first_not_of("0123456789-.", val_start);
                std::string err = raw.substr(val_start, ee-val_start);
                float ev = std::stof(err);
                snprintf(buf,sizeof(buf),"  prime %s (%s):  %+.2f¢",
                         PRIMES[i].prime, PRIMES[i].name, ev);
                browser_edo_->add(buf);
            }
            browser_edo_->redraw();
        });
}

// ─────────────────────────────────────────────────────────────────────────────
//  STATIC CALLBACKS
// ─────────────────────────────────────────────────────────────────────────────
void HarmoniaApp::cbPlayAll(Fl_Widget*, void* d) {
    HarmoniaApp* app = (HarmoniaApp*)d;
    for (auto* s : app->strips_) {
        app->audio_->noteOn(s->voice_id);
        s->btn_on->value(1);
    }
}

void HarmoniaApp::cbStopAll(Fl_Widget*, void* d) {
    HarmoniaApp* app = (HarmoniaApp*)d;
    for (auto* s : app->strips_) {
        app->audio_->noteOff(s->voice_id);
        s->btn_on->value(0);
    }
}

void HarmoniaApp::cbMasterVol(Fl_Widget* w, void* d) {
    ((HarmoniaApp*)d)->audio_->setMasterVolume((float)((Fl_Value_Slider*)w)->value());
}

void HarmoniaApp::cbKey(Fl_Widget*, void* d) {
    HarmoniaApp* app = (HarmoniaApp*)d;
    app->current_key_ = app->ch_key_->value();
    app->updateFunctionalAnalysis();
    app->updatePsychoAnalysis();
}

void HarmoniaApp::cbEDO(Fl_Widget*, void* d) {
    HarmoniaApp* app = (HarmoniaApp*)d;
    app->current_edo_ = (int)app->sp_edo_->value();
    app->tonnetz_->setEDO(app->current_edo_);
    app->audio_->setEDO(app->current_edo_);

    // Refresh analysis for the new EDO
    app->updateFunctionalAnalysis();
    app->updateCompletionSuggestions();
    app->updatePsychoAnalysis();
    app->updateEDOAnalysis();
}

void HarmoniaApp::cbAnalyze(Fl_Widget*, void* d) {
    ((HarmoniaApp*)d)->updateFunctionalAnalysis();
}

void HarmoniaApp::cbResolve(Fl_Widget*, void* d) {
    ((HarmoniaApp*)d)->updateResolutionPaths();
}

void HarmoniaApp::cbSuggestCompletion(Fl_Widget*, void* d) {
    ((HarmoniaApp*)d)->updateCompletionSuggestions();
}

void HarmoniaApp::cbPsycho(Fl_Widget*, void* d) {
    ((HarmoniaApp*)d)->updatePsychoAnalysis();
}

void HarmoniaApp::cbAnalyseEDO(Fl_Widget*, void* d) {
    ((HarmoniaApp*)d)->updateEDOAnalysis();
}

void HarmoniaApp::cbBrowserResolutions(Fl_Widget* w, void* d) {
    HarmoniaApp* app = (HarmoniaApp*)d;
    Fl_Browser* b = (Fl_Browser*)w;
    int line = b->value();
    if (line <= 0) {
        app->tonnetz_->setHighlightedPC(-1);
        return;
    }

    int path_idx = (int)(intptr_t)b->data(line);
    if (path_idx < 0 || path_idx >= (int)app->last_resolutions_.size()) {
        app->tonnetz_->setHighlightedPC(-1);
        return;
    }

    const auto& rp = app->last_resolutions_[path_idx];
    app->tonnetz_->setHighlightedPC(rp.target_root);

    // Highlight nodes for the target chord
    if (Fl::event_clicks()) {
        for (int pc : rp.target_pcs) {
            // Add voice if not present
            bool present = false;
            auto voices = app->audio_->getVoiceSnapshot();
            for (auto& v : voices) if (v.pitch_class == pc) { present = true; break; }
            if (!present) {
                double cents = (double)pc * 1200.0 / app->current_edo_;
                int closest_midi = 60 + (int)std::round(cents / 100.0);
                float detune = (float)(cents - (closest_midi - 60) * 100.0);

                app->addVoice(closest_midi);
                if (!app->strips_.empty()) {
                    int newid = app->strips_.back()->voice_id;
                    app->audio_->setVoiceDetune(newid, detune);
                    if (auto* s = app->findStrip(newid)) {
                        s->sl_detune->value(detune);
                        s->btn_on->value(1);
                    }
                    app->audio_->noteOn(newid);
                }
            }
        }
    }
}

void HarmoniaApp::cbBrowserCompletion(Fl_Widget* w, void* d) {
    HarmoniaApp* app = (HarmoniaApp*)d;
    Fl_Browser* b = (Fl_Browser*)w;
    int line = b->value();
    if (line <= 0) {
        app->tonnetz_->setHighlightedPC(-1);
        return;
    }

    int pc = (int)(intptr_t)b->data(line);
    if (pc < 0) {
        app->tonnetz_->setHighlightedPC(-1);
        return;
    }

    app->tonnetz_->setHighlightedPC(pc);

    // Single click for highlighting (already done above), double click to add
    if (Fl::event_clicks()) {
        app->onTonnetzClick(pc, 0, 0);
    }
}

void HarmoniaApp::cbVoiceOn(Fl_Widget* w, void* d) {
    HarmoniaApp* app = (HarmoniaApp*)d;
    int id = (int)(intptr_t)w->user_data();
    Fl_Light_Button* btn = (Fl_Light_Button*)w;
    if (btn->value()) app->audio_->noteOn(id);
    else              app->audio_->noteOff(id);
}

void HarmoniaApp::cbVoiceAmp(Fl_Widget* w, void* d) {
    HarmoniaApp* app = (HarmoniaApp*)d;
    int id = (int)(intptr_t)w->user_data();
    app->audio_->setVoiceAmplitude(id, (float)((Fl_Value_Slider*)w)->value());
}

void HarmoniaApp::cbVoiceDetune(Fl_Widget* w, void* d) {
    HarmoniaApp* app = (HarmoniaApp*)d;
    int id = (int)(intptr_t)w->user_data();
    app->audio_->setVoiceDetune(id, (float)((Fl_Value_Slider*)w)->value());
}

void HarmoniaApp::cbVoiceTimbre(Fl_Widget* w, void* d) {
    HarmoniaApp* app = (HarmoniaApp*)d;
    int id = (int)(intptr_t)w->user_data();
    int t  = ((Fl_Choice*)w)->value();
    static const TimbrePreset PRESETS[] = {
        TimbrePreset::SINE, TimbrePreset::SAWTOOTH, TimbrePreset::SQUARE,
        TimbrePreset::STRINGS, TimbrePreset::BRASS, TimbrePreset::FLUTE
    };
    if (t >= 0 && t < 6) app->audio_->setVoiceTimbre(id, PRESETS[t]);
}

void HarmoniaApp::cbVoiceRemove(Fl_Widget* w, void* d) {
    HarmoniaApp* app = (HarmoniaApp*)d;
    int id = (int)(intptr_t)w->user_data();
    app->removeVoice(id);
}

// ─────────────────────────────────────────────────────────────────────────────
int main(int /*argc*/, char** /*argv*/) {
    // Request FLTK lock for threaded use
    Fl::lock();

    HarmoniaApp app;
    app.run();
    return 0;
}
