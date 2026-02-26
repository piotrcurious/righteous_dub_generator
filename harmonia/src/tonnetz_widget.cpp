#include "tonnetz_widget.h"
#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <FL/gl.h>
#include <FL/glu.h>
#include <FL/Fl_Gl_Window.H>
#include <cmath>
#include <cstring>
#include <algorithm>

// ────────────────────────────────────────────────────────────────────────────
//  Pitch class → RGB colour (hue wheel by semitone)
// ────────────────────────────────────────────────────────────────────────────
static void pcColor(int pc, float& r, float& g, float& b) {
    static const float COLS[12][3] = {
        {1.0f,.35f,.35f}, {1.0f,.55f,.20f}, {1.0f,.85f,.15f}, {.75f,1.0f,.20f},
        {.30f,1.0f,.30f}, {.20f,.95f,.55f}, {.15f,.85f,.95f}, {.20f,.50f,1.0f},
        {.50f,.25f,1.0f}, {.80f,.20f,1.0f}, {1.0f,.20f,.80f}, {1.0f,.20f,.50f}
    };
    r = COLS[pc][0]; g = COLS[pc][1]; b = COLS[pc][2];
}

static const char* NOTE_NAMES_FLAT[12] = {
    "C","D♭","D","E♭","E","F","G♭","G","A♭","A","B♭","B"
};

static std::string pcLabel(int pc, int edo) {
    if (edo == 12) return NOTE_NAMES_FLAT[pc % 12];
    return std::to_string(pc);
}

// ────────────────────────────────────────────────────────────────────────────
TonnetzWidget::TonnetzWidget(int X,int Y,int W,int H,const char* l)
    : Fl_Gl_Window(X,Y,W,H,l)
{
    mode(FL_RGB | FL_DOUBLE | FL_DEPTH);
    buildLattice();
}

void TonnetzWidget::buildLattice() {
    nodes_.clear();
    for (int gy = 0; gy < GRID_Y; gy++) {
        for (int gx = 0; gx < GRID_X; gx++) {
            TonnetzNode n;
            n.x = gx - GRID_X/2;
            n.y = gy - GRID_Y/2;
            n.pitch_class = pcFromCoord(n.x, n.y);
            n.label = pcLabel(n.pitch_class, edo_);
            n.cx = n.cy = 0.f;
            nodes_.push_back(n);
        }
    }
}

void TonnetzWidget::setVoices(const std::vector<Voice>& v) {
    voices_ = v; redraw();
}
void TonnetzWidget::setAbstractObject(const AbstractObject& o) {
    abs_obj_ = o; redraw();
}
void TonnetzWidget::setRoughnessRecords(const std::vector<RoughnessRecord>& rr) {
    roughness_ = rr; redraw();
}
void TonnetzWidget::setEDO(int edo) {
    if (edo_ == edo) return;
    edo_ = edo;
    buildLattice();
    redraw();
}

// ────────────────────────────────────────────────────────────────────────────
//  Coordinate transforms
// ────────────────────────────────────────────────────────────────────────────
void TonnetzWidget::worldToScreen(float wx, float wy, float& sx, float& sy) {
    float cx = w() * 0.5f, cy = h() * 0.5f;
    sx = cx + (wx + pan_x_) * zoom_;
    sy = cy - (wy + pan_y_) * zoom_;   // OpenGL Y-up but window Y-down
}

void TonnetzWidget::screenToWorld(int sx, int sy, float& wx, float& wy) {
    float cx = w() * 0.5f, cy = h() * 0.5f;
    wx = (sx - cx) / zoom_ - pan_x_;
    wy = -(sy - cy) / zoom_ + pan_y_;
}

TonnetzNode* TonnetzWidget::nodeAt(int sx, int sy) {
    float wx, wy;
    screenToWorld(sx, sy, wx, wy);
    TonnetzNode* best = nullptr;
    float bestd = 22.f * zoom_;
    for (auto& n : nodes_) {
        // node world position
        float nx = n.x * node_dx_ + n.y * node_dx_ * 0.5f;
        float ny = n.y * node_dy_;
        float d = std::sqrt((wx-nx)*(wx-nx)+(wy-ny)*(wy-ny));
        if (d < bestd) { bestd = d; best = &n; }
    }
    return best;
}

// ────────────────────────────────────────────────────────────────────────────
//  OpenGL helpers
// ────────────────────────────────────────────────────────────────────────────
void TonnetzWidget::drawFilledCircle(float cx,float cy,float r,
                                     float red,float g,float b,float a) {
    glColor4f(red,g,b,a);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(cx, cy);
    for (int i = 0; i <= 32; i++) {
        float a2 = i * 2.f * M_PI / 32.f;
        glVertex2f(cx + r*cosf(a2), cy + r*sinf(a2));
    }
    glEnd();
}

void TonnetzWidget::drawRing(float cx,float cy,float r,float thick,
                              float red,float g,float b,float a) {
    glLineWidth(thick);
    glColor4f(red,g,b,a);
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < 48; i++) {
        float ang = i * 2.f * M_PI / 48.f;
        glVertex2f(cx + r*cosf(ang), cy + r*sinf(ang));
    }
    glEnd();
    glLineWidth(1.f);
}

void TonnetzWidget::drawTriangle(float x1,float y1,float x2,float y2,
                                  float x3,float y3,
                                  float r,float g,float b,float a) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(r,g,b,a);
    glBegin(GL_TRIANGLES);
    glVertex2f(x1,y1); glVertex2f(x2,y2); glVertex2f(x3,y3);
    glEnd();
}

void TonnetzWidget::drawLine(float x1,float y1,float x2,float y2,
                              float r,float g,float b,float lw) {
    glLineWidth(lw);
    glColor4f(r,g,b,1.f);
    glBegin(GL_LINES);
    glVertex2f(x1,y1); glVertex2f(x2,y2);
    glEnd();
    glLineWidth(1.f);
}

void TonnetzWidget::drawText(float x, float y, const char* text,
                              float r, float g, float b) {
    // Use FLTK overlay text (GL text requires GLUT or freetype)
    // We'll draw via gl_draw at transformed position
    float sx, sy;
    worldToScreen(x, y, sx, sy);
    glColor3f(r, g, b);
    gl_color(fl_rgb_color((uchar)(r*255),(uchar)(g*255),(uchar)(b*255)));
    gl_draw(text, (int)sx - (int)(strlen(text)*3), (int)sy + 4);
}

// ────────────────────────────────────────────────────────────────────────────
//  Node position in world-space
// ────────────────────────────────────────────────────────────────────────────
static void nodeWorldPos(const TonnetzNode& n, float& wx, float& wy, float dx, float dy) {
    // Hexagonal layout: x-axis = fifths, diagonal offset for thirds
    wx = n.x * dx + n.y * dx * 0.5f;
    wy = n.y * dy;
}

// ────────────────────────────────────────────────────────────────────────────
//  DRAW
// ────────────────────────────────────────────────────────────────────────────
void TonnetzWidget::draw() {
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

    // Clear
    glClearColor(0.08f, 0.09f, 0.12f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT);

    // Precompute node screen coords (store in node.cx, cy for hit-testing)
    for (auto& n : nodes_) {
        float wx, wy;
        nodeWorldPos(n, wx, wy, node_dx_, node_dy_);
        worldToScreen(wx, wy, n.cx, n.cy);
    }

    drawTriads();
    drawEdges();
    drawRoughnessHeat();
    drawAbstractObject();
    drawProgressionPath();
    drawNodes();
    drawLabels();

    // Legend
    gl_font(FL_HELVETICA, 10);
    glColor3f(.6f,.6f,.6f);
    gl_draw("Horizontal: Perfect Fifths (×3/2)", 8, h()-28);
    gl_draw("Diagonal:   Major Thirds (×5/4)",   8, h()-14);
    gl_draw("△=major  ▽=minor", 8, h()-42);
}

void TonnetzWidget::drawEdges() {
    for (int gy = 0; gy < GRID_Y; gy++) {
        for (int gx = 0; gx < GRID_X; gx++) {
            const TonnetzNode* n = nullptr;
            for (auto& nd : nodes_)
                if (nd.x == gx - GRID_X/2 && nd.y == gy - GRID_Y/2) { n = &nd; break; }
            if (!n) continue;

            // right neighbour (fifth)
            if (gx+1 < GRID_X) {
                const TonnetzNode* nr = nullptr;
                for (auto& nd : nodes_)
                    if (nd.x == gx+1-GRID_X/2 && nd.y == gy-GRID_Y/2) { nr = &nd; break; }
                if (nr) drawLine(n->cx, n->cy, nr->cx, nr->cy, .25f,.28f,.35f, 1.2f);
            }
            // diagonal upper-right (third)
            if (gx < GRID_X && gy+1 < GRID_Y) {
                const TonnetzNode* nu = nullptr;
                for (auto& nd : nodes_)
                    if (nd.x == gx-GRID_X/2 && nd.y == gy+1-GRID_Y/2) { nu = &nd; break; }
                if (nu) drawLine(n->cx, n->cy, nu->cx, nu->cy, .25f,.28f,.35f, 1.f);
            }
        }
    }
}

void TonnetzWidget::drawTriads() {
    // Major triads: (gx,gy),(gx+1,gy),(gx,gy+1)
    // Minor triads: (gx+1,gy),(gx,gy+1),(gx+1,gy+1)
    auto findNode = [&](int gx, int gy) -> const TonnetzNode* {
        for (auto& n : nodes_)
            if (n.x == gx && n.y == gy) return &n;
        return nullptr;
    };

    for (int gy = -GRID_Y/2; gy < GRID_Y/2; gy++) {
        for (int gx = -GRID_X/2; gx < GRID_X/2; gx++) {
            auto* a = findNode(gx, gy);
            auto* b = findNode(gx+1, gy);
            auto* c = findNode(gx, gy+1);
            auto* d = findNode(gx+1, gy+1);
            if (!a||!b||!c) continue;

            // Check if any voice is active with matching pitch classes
            auto voiceActive = [&](int pc) {
                for (auto& v : voices_) if (v.active && v.pitch_class==pc) return true;
                return false;
            };

            bool maj_lit = voiceActive(a->pitch_class) && voiceActive(b->pitch_class) && voiceActive(c->pitch_class);
            float alpha = maj_lit ? 0.55f : 0.08f;
            float pr,pg,pb; pcColor(a->pitch_class, pr,pg,pb);
            drawTriangle(a->cx,a->cy, b->cx,b->cy, c->cx,c->cy, pr,pg,pb, alpha);

            if (!d) continue;
            bool min_lit = voiceActive(b->pitch_class) && voiceActive(c->pitch_class) && voiceActive(d->pitch_class);
            alpha = min_lit ? 0.45f : 0.06f;
            pcColor(c->pitch_class, pr,pg,pb);
            drawTriangle(b->cx,b->cy, c->cx,c->cy, d->cx,d->cy, pr*0.6f,pg*0.6f,pb*0.6f, alpha);
        }
    }
}

void TonnetzWidget::drawAbstractObject() {
    if (abs_obj_.confidence < 0.2f || abs_obj_.voice_ids.empty()) return;
    // Draw a pulsing halo at the centroid
    // Find all active voices, average their screen positions
    float sx = 0, sy = 0; int cnt = 0;
    for (auto& v : voices_) {
        if (!v.active) continue;
        // find node
        float wx = v.tonnetz_x * node_dx_ + v.tonnetz_y * node_dx_ * 0.5f;
        float wy = v.tonnetz_y * node_dy_;
        float nsx, nsy; worldToScreen(wx, wy, nsx, nsy);
        sx += nsx; sy += nsy; cnt++;
    }
    if (!cnt) return;
    sx /= cnt; sy /= cnt;

    // pulsing radius
    static float pulse = 0.f;
    pulse += 0.04f;
    float r = 40.f + 8.f * sinf(pulse);

    // glow ring
    float conf = abs_obj_.confidence;
    drawRing(sx,sy, r,     4.f, 1.f,.85f,.2f, conf*0.9f);
    drawRing(sx,sy, r+8.f, 2.f, 1.f,.85f,.2f, conf*0.4f);
    drawRing(sx,sy, r+16.f,1.f, 1.f,.85f,.2f, conf*0.15f);

    // label
    gl_font(FL_HELVETICA_BOLD, 14);
    glColor4f(1.f,.9f,.3f, conf);
    std::string label = abs_obj_.chord_name;
    if (abs_obj_.virtual_pitch_hz > 30.f) {
        char buf[64];
        snprintf(buf, sizeof(buf), " [%.0fHz]", abs_obj_.virtual_pitch_hz);
        label += buf;
    }
    int tw = (int)gl_width(label.c_str());
    gl_draw(label.c_str(), (int)sx - tw/2, (int)sy - (int)r - 8);
}

void TonnetzWidget::drawProgressionPath() {
    if (progression_path_.size() < 2) return;
    auto findNode = [&](int gx, int gy) -> const TonnetzNode* {
        for (auto& n : nodes_) if (n.x==gx && n.y==gy) return &n;
        return nullptr;
    };
    for (size_t i = 0; i+1 < progression_path_.size(); i++) {
        auto* a = findNode(progression_path_[i].first,   progression_path_[i].second);
        auto* b = findNode(progression_path_[i+1].first, progression_path_[i+1].second);
        if (!a||!b) continue;
        float t = (float)i / (progression_path_.size()-1);
        drawLine(a->cx,a->cy, b->cx,b->cy, .2f+.6f*t, .7f, 1.f-t*.5f, 2.5f);
        // Arrow tip
        float dx = b->cx-a->cx, dy = b->cy-a->cy;
        float len = sqrtf(dx*dx+dy*dy);
        if (len > 0.f) { dx/=len; dy/=len; }
        drawFilledCircle(b->cx-dx*22.f, b->cy-dy*22.f, 5.f, .2f,.9f,.8f, 0.8f);
    }
}

void TonnetzWidget::drawNodes() {
    gl_font(FL_HELVETICA, 11);

    for (auto& n : nodes_) {
        // Is any voice playing this pitch class?
        bool has_voice = false;
        float vr=0,vg=0,vb=0;
        for (auto& v : voices_) {
            if (v.active && v.pitch_class == n.pitch_class) {
                has_voice = true; vr=v.color[0]; vg=v.color[1]; vb=v.color[2];
            }
        }

        float pr,pg,pb;
        pcColor(n.pitch_class, pr,pg,pb);

        if (has_voice) {
            // Active voice: bright node with glow
            drawFilledCircle(n.cx, n.cy, 20.f, vr*.3f,vg*.3f,vb*.3f, 0.5f);
            drawFilledCircle(n.cx, n.cy, 18.f, vr,vg,vb, 1.f);
            drawRing(n.cx, n.cy, 22.f, 2.f, 1.f,1.f,1.f, 0.7f);
        } else {
            // Inactive: dark node
            drawFilledCircle(n.cx, n.cy, 16.f, pr*.15f,pg*.15f,pb*.15f, 1.f);
            drawFilledCircle(n.cx, n.cy, 14.f, pr*.4f, pg*.4f, pb*.4f,  1.f);
            drawRing(n.cx, n.cy, 16.f, 1.f, pr,pg,pb, 0.5f);
        }
    }
}

void TonnetzWidget::drawLabels() {
    gl_font(FL_HELVETICA_BOLD, 11);
    for (auto& n : nodes_) {
        bool has_voice = false;
        for (auto& v : voices_)
            if (v.active && v.pitch_class == n.pitch_class) { has_voice = true; break; }

        if (has_voice) glColor3f(0.05f,0.05f,0.05f);
        else           glColor3f(0.85f,0.85f,0.85f);

        int tw = (int)gl_width(n.label.c_str());
        gl_draw(n.label.c_str(), (int)n.cx - tw/2, (int)n.cy + 4);
    }
}

void TonnetzWidget::drawRoughnessHeat() {
    for (const auto& r : roughness_) {
        if (r.roughness < 0.01f) continue;

        const Voice *va = nullptr, *vb = nullptr;
        for (const auto& v : voices_) {
            if (v.id == r.voice_a) va = &v;
            if (v.id == r.voice_b) vb = &v;
        }
        if (!va || !vb) continue;

        // Get world positions
        float wax = va->tonnetz_x * node_dx_ + va->tonnetz_y * node_dx_ * 0.5f;
        float way = va->tonnetz_y * node_dy_;
        float wbx = vb->tonnetz_x * node_dx_ + vb->tonnetz_y * node_dx_ * 0.5f;
        float wby = vb->tonnetz_y * node_dy_;

        float sax, say, sbx, sby;
        worldToScreen(wax, way, sax, say);
        worldToScreen(wbx, wby, sbx, sby);

        // Draw glowing edge for roughness
        float intensity = std::min(1.0f, r.roughness * 2.0f);
        drawLine(sax, say, sbx, sby, 1.0f, 0.4f, 0.2f, 1.0f + intensity * 6.0f);
    }
}

// ────────────────────────────────────────────────────────────────────────────
//  Event handling
// ────────────────────────────────────────────────────────────────────────────
int TonnetzWidget::handle(int event) {
    switch (event) {
    case FL_PUSH:
        drag_x_ = Fl::event_x(); drag_y_ = Fl::event_y();
        dragging_ = (Fl::event_button() == FL_MIDDLE_MOUSE);
        if (Fl::event_button() == FL_LEFT_MOUSE) {
            TonnetzNode* n = nodeAt(drag_x_, drag_y_);
            if (n && node_click_cb_) {
                node_click_cb_(n->pitch_class, n->x, n->y);
            }
        }
        return 1;
    case FL_DRAG:
        if (dragging_) {
            pan_x_ += (Fl::event_x() - drag_x_) / zoom_;
            pan_y_ -= (Fl::event_y() - drag_y_) / zoom_;
            drag_x_ = Fl::event_x(); drag_y_ = Fl::event_y();
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
