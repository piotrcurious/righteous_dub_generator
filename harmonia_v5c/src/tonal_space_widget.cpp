#include "tonal_space_widget.h"
#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <FL/gl.h>
#include <FL/glu.h>
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const char* NOTE_NAMES_FLAT[12] = {
    "C","D♭","D","E♭","E","F","G♭","G","A♭","A","B♭","B"
};

static std::string pcLabel(int pc, int edo) {
    pc = ((pc % edo) + edo) % edo;
    if (edo == 12) return NOTE_NAMES_FLAT[pc];

    double exact_semitones = pc * 12.0 / edo;
    int idx12 = (int)std::round(exact_semitones) % 12;
    int dev_cents = (int)std::round((exact_semitones - std::round(exact_semitones)) * 100);
    std::string base = NOTE_NAMES_FLAT[idx12];
    if (dev_cents == 0) return base;
    std::string sign = (dev_cents > 0) ? "+" : "";
    return base + sign + std::to_string(dev_cents) + "¢";
}

// ────────────────────────────────────────────────────────────────────────────
TonalSpaceWidget::TonalSpaceWidget(int X,int Y,int W,int H,const char* l)
    : Fl_Gl_Window(X,Y,W,H,l)
{
    mode(FL_RGB | FL_DOUBLE | FL_DEPTH);
    buildSpace();
}

void TonalSpaceWidget::buildSpace() {
    nodes_.clear();
    float base_radius = 80.0f;
    float ring_step = 45.0f;
    for (int oct = 2; oct <= 6; oct++) {
        for (int pc = 0; pc < edo_; pc++) {
            TonalNode n;
            n.pitch_class = pc;
            n.octave = oct;
            // Angle: PC 0 is at PI/2 (top). progression is clockwise.
            n.angle = (float)(M_PI / 2.0 - ((double)pc / edo_) * 2.0 * M_PI);
            n.radius = base_radius + (oct - 2) * ring_step;
            n.label = (oct == 4) ? pcLabel(pc, edo_) : "";
            n.cx = n.cy = 0.f;
            // Expected frequency: C4 is the reference (pc=0, oct=4).
            // Each EDO step is one equal-tempered step in the current tuning.
            n.freq = (float)(C4_HZ * std::pow(2.0,
                        (n.pitch_class + (n.octave - 4.0) * edo_) / edo_));
            nodes_.push_back(n);
        }
    }
}

void TonalSpaceWidget::setVoices(const std::vector<Voice>& v) {
    voices_ = v; redraw();
}
void TonalSpaceWidget::setAbstractObject(const AbstractObject& o) {
    abs_obj_ = o; redraw();
}
void TonalSpaceWidget::setRoughnessRecords(const std::vector<RoughnessRecord>& rr) {
    roughness_ = rr; redraw();
}
void TonalSpaceWidget::setEDO(int edo) {
    if (edo_ == edo) return;
    edo_ = edo;
    buildSpace();
    redraw();
}

// ────────────────────────────────────────────────────────────────────────────
//  Coordinate transforms
// ────────────────────────────────────────────────────────────────────────────
void TonalSpaceWidget::worldToScreen(float wx, float wy, float& sx, float& sy) {
    float center_x = w() * 0.5f, center_y = h() * 0.5f;
    sx = center_x + (wx + pan_x_) * zoom_;
    sy = center_y - (wy + pan_y_) * zoom_;
}

void TonalSpaceWidget::screenToWorld(int sx, int sy, float& wx, float& wy) {
    float center_x = w() * 0.5f, center_y = h() * 0.5f;
    wx = ((float)sx - center_x) / zoom_ - pan_x_;
    wy = (center_y - (float)sy) / zoom_ - pan_y_;
}

TonalNode* TonalSpaceWidget::nodeAt(int sx, int sy) {
    updateNodeCoordinates();
    TonalNode* best = nullptr;
    float best_dist = 25.0f; // pixel radius

    for (auto& n : nodes_) {
        float dist = std::hypot((float)sx - n.cx, (float)sy - n.cy);
        if (dist < best_dist) {
            best_dist = dist;
            best = &n;
        }
    }
    return best;
}

// ────────────────────────────────────────────────────────────────────────────
//  OpenGL helpers
// ────────────────────────────────────────────────────────────────────────────
void TonalSpaceWidget::drawFilledCircle(float cx,float cy,float r,
                                     float red,float g,float b,float a) {
    glColor4f(red,g,b,a);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(cx, cy);
    for (int i = 0; i <= 32; i++) {
        float a2 = i * 2.f * (float)M_PI / 32.f;
        glVertex2f(cx + r*cosf(a2), cy + r*sinf(a2));
    }
    glEnd();
}

void TonalSpaceWidget::drawRing(float cx,float cy,float r,float thick,
                              float red,float g,float b,float a) {
    glLineWidth(thick);
    glColor4f(red,g,b,a);
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < 48; i++) {
        float ang = i * 2.f * (float)M_PI / 48.f;
        glVertex2f(cx + r*cosf(ang), cy + r*sinf(ang));
    }
    glEnd();
    glLineWidth(1.f);
}

void TonalSpaceWidget::drawLine(float x1,float y1,float x2,float y2,
                              float r,float g,float b,float lw) {
    glLineWidth(lw);
    glColor4f(r,g,b,1.f);
    glBegin(GL_LINES);
    glVertex2f(x1,y1); glVertex2f(x2,y2);
    glEnd();
    glLineWidth(1.f);
}

// ────────────────────────────────────────────────────────────────────────────
//  DRAW
// ────────────────────────────────────────────────────────────────────────────
void TonalSpaceWidget::updateNodeCoordinates() {
    for (auto& n : nodes_) {
        float wx = n.radius * std::cos(n.angle + rotation_);
        float wy = n.radius * std::sin(n.angle + rotation_);
        worldToScreen(wx, wy, n.cx, n.cy);
    }
}

void TonalSpaceWidget::draw() {
    if (!valid()) {
        glViewport(0, 0, w(), h());
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0, w(), h(), 0, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    }

    glClearColor(0.08f, 0.09f, 0.12f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT);

    updateNodeCoordinates();

    drawConnections();
    drawRoughnessHeat();
    drawAbstractObject();
    drawNodes();
    drawLabels();

    // Legend
    gl_font(FL_HELVETICA, 10);
    glColor3f(.6f,.6f,.6f);
    gl_draw("Circular Pitch Class Space", 8, h()-14);
    char buf[64];
    snprintf(buf, sizeof(buf), "EDO: %d", edo_);
    gl_draw(buf, 8, h()-28);
}

void TonalSpaceWidget::drawConnections() {
    int fifth_steps = (int)std::round(std::log2(1.5) * edo_);
    int third_steps = (int)std::round(std::log2(1.25) * edo_);

    for (int oct = 2; oct <= 6; oct++) {
        int offset = (oct - 2) * edo_;
        for (int i = 0; i < edo_; i++) {
            int next_fifth = (i + fifth_steps) % edo_;
            int next_third = (i + third_steps) % edo_;

            drawLine(nodes_[offset + i].cx, nodes_[offset + i].cy,
                     nodes_[offset + next_fifth].cx, nodes_[offset + next_fifth].cy,
                     0.2f, 0.25f, 0.35f, 1.0f);
            drawLine(nodes_[offset + i].cx, nodes_[offset + i].cy,
                     nodes_[offset + next_third].cx, nodes_[offset + next_third].cy,
                     0.25f, 0.2f, 0.3f, 0.8f);

            // Octave connections
            if (oct < 6) {
                int next_oct_offset = (oct - 1) * edo_;
                drawLine(nodes_[offset + i].cx, nodes_[offset + i].cy,
                         nodes_[next_oct_offset + i].cx, nodes_[next_oct_offset + i].cy,
                         0.15f, 0.15f, 0.2f, 0.5f);
            }
        }
    }
}

void TonalSpaceWidget::drawNodes() {
    for (auto& n : nodes_) {
        bool has_voice = false;
        float vr=0, vg=0, vb=0;
        for (auto& v : voices_) {
            // Match by frequency proximity (within half an EDO step) rather
            // than by pitch_class integer.  The audio engine may store
            // pitch_class in 12-EDO units regardless of the active tuning,
            // so an integer comparison would select the wrong node for any
            // EDO != 12.
            if (v.active && v.frequency > 0.f && n.freq > 0.f) {
                float ratio = std::abs((float)std::log2(v.frequency / n.freq));
                if (ratio < 0.5f / edo_) {
                    has_voice = true;
                    vr = v.color[0]; vg = v.color[1]; vb = v.color[2];
                    break;
                }
            }
        }

        float pr, pg, pb;
        Voice::pcColorHSV(n.pitch_class, edo_, pr, pg, pb);

        bool is_highlighted = (std::find(highlighted_pcs_.begin(), highlighted_pcs_.end(), n.pitch_class) != highlighted_pcs_.end()) ||
                              (n.pitch_class == highlighted_pc_ && (highlighted_oct_ == -1 || n.octave == highlighted_oct_));

        if (has_voice) {
            drawFilledCircle(n.cx, n.cy, 12.f, vr, vg, vb, 1.f);
            drawRing(n.cx, n.cy, 15.f, 2.f, 1.f, 1.f, 1.f, 0.7f);
            if (is_highlighted) drawRing(n.cx, n.cy, 18.f, 3.f, 1.f, 1.0f, 0.4f, 0.9f);
        } else if (is_highlighted) {
            drawFilledCircle(n.cx, n.cy, 10.f, pr, pg, pb, 0.8f);
            drawRing(n.cx, n.cy, 12.f, 2.f, 1.0f, 1.0f, 1.0f, 0.9f);
        } else {
            drawFilledCircle(n.cx, n.cy, 8.f, pr*.4f, pg*.4f, pb*.4f, 1.f);
            drawRing(n.cx, n.cy, 8.f, 1.f, pr, pg, pb, 0.5f);
        }
    }
}

void TonalSpaceWidget::drawLabels() {
    gl_font(FL_HELVETICA_BOLD, 10);
    for (auto& n : nodes_) {
        if (n.label.empty()) continue;
        bool has_voice = false;
        for (auto& v : voices_) {
            if (v.active && v.frequency > 0.f && n.freq > 0.f) {
                float ratio = std::abs((float)std::log2(v.frequency / n.freq));
                if (ratio < 0.5f / edo_) { has_voice = true; break; }
            }
        }

        if (has_voice) glColor3f(0.05f,0.05f,0.05f);
        else           glColor3f(0.85f,0.85f,0.85f);

        int tw = (int)gl_width(n.label.c_str());
        gl_draw(n.label.c_str(), (int)n.cx - tw/2, (int)n.cy + 4);
    }
}

void TonalSpaceWidget::drawAbstractObject() {
    if (abs_obj_.confidence < 0.2f) return;

    float center_x = w() * 0.5f, center_y = h() * 0.5f;
    static float pulse = 0.f; pulse += 0.04f;
    float r = 50.f + 10.f * sinf(pulse);

    drawRing(center_x, center_y, r, 2.f, 1.f, 0.85f, 0.2f, abs_obj_.confidence);
}

void TonalSpaceWidget::drawRoughnessHeat() {
    for (const auto& r : roughness_) {
        if (r.roughness < 0.01f) continue;
        const Voice *va = nullptr, *vb = nullptr;
        for (const auto& v : voices_) {
            if (v.id == r.voice_a) va = &v;
            if (v.id == r.voice_b) vb = &v;
        }
        if (!va || !vb) continue;

        auto findNode = [&](float vfreq) -> TonalNode* {
            // Find the node whose expected frequency is closest to vfreq,
            // within half an EDO step.  Frequency-based lookup is correct for
            // any EDO because it doesn't assume the voice's pitch_class field
            // is stored in the same EDO units as the widget's node indices.
            TonalNode* best = nullptr;
            float best_ratio = 0.5f / edo_;
            for (auto& node : nodes_) {
                if (node.freq <= 0.f) continue;
                float ratio = std::abs((float)std::log2(vfreq / node.freq));
                if (ratio < best_ratio) {
                    best_ratio = ratio;
                    best = &node;
                }
            }
            return best;
        };

        TonalNode *na = findNode(va->frequency);
        TonalNode *nb = findNode(vb->frequency);
        if (!na || !nb) continue;

        float intensity = std::min(1.0f, r.roughness * 2.0f);
        drawLine(na->cx, na->cy, nb->cx, nb->cy, 1.0f, 0.4f, 0.2f, 1.0f + intensity * 6.0f);
    }
}

int TonalSpaceWidget::handle(int event) {
    // For an Fl_Gl_Window, events are local to the widget's origin.
    int lx = Fl::event_x();
    int ly = Fl::event_y();

    switch (event) {
    case FL_PUSH:
        drag_x_ = lx; drag_y_ = ly;
        dragging_ = (Fl::event_button() == FL_MIDDLE_MOUSE);
        if (Fl::event_button() == FL_LEFT_MOUSE) {
            TonalNode* n = nodeAt(lx, ly);
            if (n && node_click_cb_) {
                // Return PC and octave as surrogates
                node_click_cb_(n->pitch_class, n->pitch_class, n->octave);
            }
        }
        return 1;
    case FL_DRAG:
        if (dragging_) {
            // Dragging left-to-right (positive lx delta) rotates the wheel clockwise (decreasing angle).
            rotation_ -= (lx - drag_x_) * 0.01f;
            drag_x_ = lx; drag_y_ = ly;
            redraw();
        }
        return 1;
    case FL_RELEASE:
        dragging_ = false;
        return 1;
    case FL_MOUSEWHEEL:
        zoom_ *= (Fl::event_dy() < 0) ? 1.1f : 0.91f;
        zoom_ = std::max(0.3f, std::min(4.f, zoom_));
        redraw();
        return 1;
    default:
        return Fl_Gl_Window::handle(event);
    }
}
