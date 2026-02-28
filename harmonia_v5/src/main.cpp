/*
 * harmonia/src/main.cpp
 * ─────────────────────────────────────────────────────────────────────────────
 *  HARMONIA V5 — Psychoacoustic Music Coder with Circular Tonal Space
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

static const Fl_Color COL_BG     = fl_rgb_color(18,20,28);
static const Fl_Color COL_ACCENT = fl_rgb_color(200,160,40);

// ─────────────────────────────────────────────────────────────────────────────
//  SPECTRUM DISPLAY
// ─────────────────────────────────────────────────────────────────────────────
class SpectrumWidget : public Fl_Gl_Window {
public:
    SpectrumWidget(int X,int Y,int W,int H) : Fl_Gl_Window(X,Y,W,H) { mode(FL_RGB|FL_DOUBLE); }
    void update(const SpectrumSnapshot& snap, float roughness) { snap_ = snap; roughness_ = roughness; redraw(); }
    void draw() override {
        if (!valid()) {
            glViewport(0,0,w(),h()); glMatrixMode(GL_PROJECTION); glLoadIdentity();
            glOrtho(0,w(),h(),0,-1,1); glMatrixMode(GL_MODELVIEW); glLoadIdentity();
            glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        }
        glClearColor(0.07f,0.08f,0.11f,1.f); glClear(GL_COLOR_BUFFER_BIT);
        int sw = w(), sh = h()-20; float bin_w = (float)sw / (SPECTRUM_SIZE/2);
        for (int i = 0; i < SPECTRUM_SIZE/2; i++) {
            float mag = snap_.magnitude[i]; if (mag < 1e-5f) continue;
            float bh = std::min((float)sh, mag * sh * 3.f);
            float t = (float)i / (SPECTRUM_SIZE/2);
            glColor4f(0.2f + 0.8f*t, 0.7f - 0.5f*t, 1.f - 0.7f*t, 0.85f);
            glBegin(GL_QUADS); glVertex2f(i*bin_w, sh); glVertex2f((i+1)*bin_w, sh);
            glVertex2f((i+1)*bin_w, sh-bh); glVertex2f(i*bin_w, sh-bh); glEnd();
        }
        float rw = roughness_ * sw; glColor4f(roughness_, 1.f-roughness_, 0.f, 0.8f);
        glBegin(GL_QUADS); glVertex2f(0,sh); glVertex2f(rw,sh); glVertex2f(rw,h()); glVertex2f(0,h()); glEnd();
    }
private:
    SpectrumSnapshot snap_; float roughness_{0.f};
};

// ─────────────────────────────────────────────────────────────────────────────
//  VOICE STRIP
// ─────────────────────────────────────────────────────────────────────────────
struct VoiceStrip {
    int voice_id{-1}; bool manual{false};
    Fl_Group* group{nullptr}; Fl_Light_Button* btn_on{nullptr};
    Fl_Box* lbl_note{nullptr}; Fl_Value_Slider* sl_amp{nullptr};
    Fl_Choice* ch_timbre{nullptr};
};

// ─────────────────────────────────────────────────────────────────────────────
//  MAIN APP
// ─────────────────────────────────────────────────────────────────────────────
class HarmoniaApp {
public:
    HarmoniaApp(); ~HarmoniaApp(); void run();
private:
    static constexpr int WIN_W = 1280, WIN_H = 820, TOOLBAR_H = 38, VOICE_W = 260, THEORY_W = 320, STATUS_H = 120;
    Fl_Double_Window* win_; Fl_Button *btn_add_, *btn_play_all_, *btn_stop_all_;
    Fl_Choice *ch_key_; Fl_Spinner *sp_edo_; Fl_Value_Slider *sl_master_; Fl_Output *out_object_;
    Fl_Scroll* scroll_voices_; std::vector<VoiceStrip*> strips_;
    TonalSpaceWidget* tonal_space_; Fl_Tabs* tabs_right_;
    Fl_Browser *browser_function_, *browser_psycho_;
    SpectrumWidget* spectrum_; Fl_Box* box_object_label_;
    std::unique_ptr<AudioEngine> audio_; std::unique_ptr<TheoryBridge> theory_;
    int current_key_{0}, current_edo_{12};
    std::vector<VoiceStrip*> pool_strips_;

    void buildUI(); void setupCallbacks();
    void addVoice(int midi_note = 60, TimbrePreset t = TimbrePreset::SINE);
    void removeVoice(int voice_id); void relayoutStrips();
    void updateTheory();
    static void onIdle(void* data);
    void onTonalSpaceClick(int pc, int oct);
    std::string theoryDir_;
};

HarmoniaApp::HarmoniaApp() {
    audio_ = std::make_unique<AudioEngine>(); theory_ = std::make_unique<TheoryBridge>();
    buildUI(); setupCallbacks();
}
HarmoniaApp::~HarmoniaApp() { audio_->shutdown(); theory_->stop(); }

void HarmoniaApp::run() {
    if (!audio_->init()) fl_message("ALSA init failed");
    char exepath[512] = {}; ssize_t n = readlink("/proc/self/exe", exepath, sizeof(exepath)-1);
    theoryDir_ = (n > 0) ? std::string(exepath, n).substr(0, std::string(exepath, n).rfind('/')) + "/../theory" : "./theory";
    if (!theory_->start(theoryDir_)) fl_message("Theory server failed");
    addVoice(60); addVoice(64); addVoice(67);
    for (auto& s : strips_) audio_->noteOn(s->voice_id);
    updateTheory();
    Fl::add_idle(onIdle, this); Fl::run();
}

void HarmoniaApp::buildUI() {
    Fl::scheme("gtk+"); win_ = new Fl_Double_Window(WIN_W, WIN_H, "Harmonia V5 — Circular Tonal Space"); win_->color(COL_BG);
    {
        Fl_Group* tb = new Fl_Group(0, 0, WIN_W, TOOLBAR_H); tb->box(FL_FLAT_BOX); tb->color(fl_rgb_color(14,16,22));
        int x = 8, y = 6, bh = 26;
        Fl_Box* title = new Fl_Box(x,y,100,bh,"HARMONIA V5"); title->labelcolor(COL_ACCENT); title->labelfont(FL_HELVETICA_BOLD); x += 110;
        ch_key_ = new Fl_Choice(x+35,y,70,bh,"Key:"); for(int i=0;i<12;i++) ch_key_->add(NOTE_NAMES[i]); ch_key_->value(0); x+=110;
        sp_edo_ = new Fl_Spinner(x+35,y,50,bh,"EDO:"); sp_edo_->minimum(5); sp_edo_->maximum(72); sp_edo_->value(12); x+=90;
        sl_master_ = new Fl_Value_Slider(x+50,y,100,bh,"Vol:"); sl_master_->type(FL_HORIZONTAL); sl_master_->value(0.5); x+=160;
        btn_add_ = new Fl_Button(x,y,90,bh,"+ Voice"); x+=95;
        btn_play_all_ = new Fl_Button(x,y,60,bh,"Play"); x+=65;
        btn_stop_all_ = new Fl_Button(x,y,60,bh,"Stop"); x+=65;
        out_object_ = new Fl_Output(x,y, WIN_W-x-8, bh); out_object_->color(fl_rgb_color(20,24,32)); out_object_->textcolor(COL_ACCENT);
        tb->end();
    }
    int cy = TOOLBAR_H, ch = WIN_H - TOOLBAR_H - STATUS_H;
    scroll_voices_ = new Fl_Scroll(0, cy, VOICE_W, ch); scroll_voices_->type(Fl_Scroll::VERTICAL_ALWAYS); scroll_voices_->end();
    tonal_space_ = new TonalSpaceWidget(VOICE_W, cy, WIN_W-VOICE_W-THEORY_W, ch);
    tabs_right_ = new Fl_Tabs(WIN_W-THEORY_W, cy, THEORY_W, ch);
    {
        Fl_Group* gt = new Fl_Group(WIN_W-THEORY_W, cy+25, THEORY_W, ch-25, "Theory");
        browser_function_ = new Fl_Browser(tabs_right_->x()+4, cy+45, THEORY_W-8, ch/2-50);
        browser_psycho_ = new Fl_Browser(tabs_right_->x()+4, cy+ch/2+10, THEORY_W-8, ch/2-50);
        gt->end();
    }
    tabs_right_->end();
    spectrum_ = new SpectrumWidget(0, WIN_H-STATUS_H, WIN_W-250, STATUS_H);
    box_object_label_ = new Fl_Box(WIN_W-248, WIN_H-STATUS_H, 248, STATUS_H, "");
    win_->end(); win_->show();
}

void HarmoniaApp::setupCallbacks() {
    btn_add_->callback([](Fl_Widget*, void* d){ ((HarmoniaApp*)d)->addVoice(); }, this);
    btn_play_all_->callback([](Fl_Widget*, void* d){ auto* a = (HarmoniaApp*)d; for(auto* s:a->strips_) a->audio_->noteOn(s->voice_id); }, this);
    btn_stop_all_->callback([](Fl_Widget*, void* d){ auto* a = (HarmoniaApp*)d; for(auto* s:a->strips_) a->audio_->noteOff(s->voice_id); }, this);
    sp_edo_->callback([](Fl_Widget*, void* d){
        auto* a = (HarmoniaApp*)d; a->current_edo_ = (int)a->sp_edo_->value();
        a->tonal_space_->setEDO(a->current_edo_); a->audio_->setEDO(a->current_edo_); a->updateTheory();
    }, this);
    ch_key_->callback([](Fl_Widget*, void* d){ ((HarmoniaApp*)d)->current_key_ = ((Fl_Choice*)((HarmoniaApp*)d)->ch_key_)->value(); ((HarmoniaApp*)d)->updateTheory(); }, this);
    tonal_space_->setNodeClickCallback([this](int pc, int tx, int ty){ onTonalSpaceClick(pc, ty); });
}

void HarmoniaApp::addVoice(int midi_note, TimbrePreset t) {
    int id = audio_->addVoice(midi_note, t); if (id < 0) return;
    VoiceStrip* s = new VoiceStrip(); s->voice_id = id;
    s->group = new Fl_Group(4, 0, VOICE_W-20, 60); s->group->box(FL_FLAT_BOX); s->group->color(fl_rgb_color(30,35,45));
    s->lbl_note = new Fl_Box(s->group->x()+5, s->group->y()+5, 80, 25, "Note"); s->lbl_note->labelcolor(COL_ACCENT);
    s->btn_on = new Fl_Light_Button(s->group->x()+90, s->group->y()+5, 50, 25, "PLAY");
    s->btn_on->callback([](Fl_Widget* w, void* d){
        int id = (int)(intptr_t)w->user_data(); auto* a = (HarmoniaApp*)d;
        if (((Fl_Light_Button*)w)->value()) a->audio_->noteOn(id); else a->audio_->noteOff(id);
    }, this); s->btn_on->user_data((void*)(intptr_t)id);
    auto voices = audio_->getVoiceSnapshot();
    for(auto& v : voices) if(v.id == id) {
        char buf[16];
        if (current_edo_ == 12) snprintf(buf,16,"%s%d", NOTE_NAMES[v.pitch_class%12], v.octave);
        else snprintf(buf,16,"%d/%d", v.pitch_class, v.octave);
        s->lbl_note->copy_label(buf);
    }
    s->group->end(); scroll_voices_->add(s->group); strips_.push_back(s); relayoutStrips();
}

void HarmoniaApp::relayoutStrips() {
    int y = 5; for(auto* s : strips_) { s->group->position(s->group->x(), y); y += 65; }
    scroll_voices_->init_sizes();
}

void HarmoniaApp::removeVoice(int id) {
    audio_->noteOff(id); audio_->removeVoice(id);
    auto it = std::find_if(strips_.begin(), strips_.end(), [id](VoiceStrip* s){ return s->voice_id == id; });
    if (it != strips_.end()) { scroll_voices_->remove((*it)->group); strips_.erase(it); relayoutStrips(); }
}

void HarmoniaApp::onIdle(void* data) {
    HarmoniaApp* app = (HarmoniaApp*)data; app->theory_->poll();
    auto snap = app->audio_->getSpectrumSnapshot(); auto rr = app->audio_->getRoughnessSnapshot(); auto obj = app->audio_->getAbstractObject();
    float tr = 0; for(auto& r : rr) tr += r.roughness; app->spectrum_->update(snap, std::min(1.f, tr));
    auto voices = app->audio_->getVoiceSnapshot();
    app->tonal_space_->setVoices(voices); app->tonal_space_->setAbstractObject(obj); app->tonal_space_->setRoughnessRecords(rr);
    if (obj.confidence > 0.2f) { app->out_object_->value(obj.chord_name.c_str()); app->box_object_label_->copy_label(obj.chord_name.c_str()); }
    for(auto& v : voices) {
        auto it = std::find_if(app->strips_.begin(), app->strips_.end(), [&](VoiceStrip* s){ return s->voice_id == v.id; });
        if (it != app->strips_.end()) {
            (*it)->btn_on->value(v.note_on.load() ? 1 : 0);
            char buf[16];
            if (app->current_edo_ == 12) snprintf(buf,16,"%s%d", NOTE_NAMES[v.pitch_class%12], v.octave);
            else snprintf(buf,16,"%d/%d", v.pitch_class, v.octave);
            (*it)->lbl_note->copy_label(buf);
        }
    }
}

void HarmoniaApp::onTonalSpaceClick(int pc, int oct) {
    double freq = 261.625565 * std::pow(2.0, (double)pc / current_edo_ + (oct - 4));
    auto voices = audio_->getVoiceSnapshot(); int eid = -1;
    for (auto& v : voices) if (v.pitch_class == pc && v.octave == oct) { eid = v.id; break; }
    if (eid != -1) removeVoice(eid);
    else {
        int id = audio_->addVoice(60); if (id >= 0) {
            audio_->setVoiceFrequency(id, freq); audio_->noteOn(id);
            auto it = std::find_if(strips_.begin(), strips_.end(), [&](VoiceStrip* s){ return s->voice_id == id; });
            if (it != strips_.end()) (*it)->manual = true;
        }
    }
    updateTheory();
}

void HarmoniaApp::updateTheory() {
    auto voices = audio_->getVoiceSnapshot(); std::vector<int> pcs;
    for (auto& v : voices) if(v.active) pcs.push_back(v.pitch_class); if (pcs.empty()) return;
    theory_->queryAnalyzeChord(pcs, current_key_, current_edo_, [this](const FunctionalAnalysis& fa, const TonnetzTension& tt) {
        browser_function_->clear(); char buf[256]; snprintf(buf,256,"@bFUNCTION: %s", fa.function.c_str()); browser_function_->add(buf);
        browser_function_->add(fa.function_reason.c_str());
    });
    theory_->queryPsychoacoustic(pcs, current_key_, 261.63f, 4, 70.0f, [this](const PsychoacousticAnalysis& pa) {
        browser_psycho_->clear(); char buf[256]; snprintf(buf,256,"@bTENSION: %s (%.2f)", pa.perceptual_tension_label.c_str(), pa.perceptual_tension);
        browser_psycho_->add(buf); snprintf(buf,256,"Consonance: %.2f", pa.level1.consonance_score); browser_psycho_->add(buf);
    });
}

int main() { Fl::lock(); HarmoniaApp app; app.run(); return 0; }
