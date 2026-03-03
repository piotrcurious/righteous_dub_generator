#include "psycho_graph_widget.h"
#include <FL/gl.h>
#include <FL/fl_draw.H>
#include <cmath>
#include <algorithm>

PsychoGraphWidget::PsychoGraphWidget(int X, int Y, int W, int H, const char* label)
    : Fl_Gl_Window(X, Y, W, H, label) {
    mode(FL_RGB | FL_DOUBLE);
}

void PsychoGraphWidget::update(const PsychoacousticAnalysis& pa) {
    pa_ = pa;
    data_available_ = true;
    redraw();
}

void PsychoGraphWidget::draw() {
    if (!valid()) {
        glViewport(0, 0, w(), h());
        glMatrixMode(GL_PROJECTION); glLoadIdentity();
        glOrtho(0, w(), h(), 0, -1, 1);
        glMatrixMode(GL_MODELVIEW); glLoadIdentity();
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }

    glClearColor(0.05f, 0.06f, 0.08f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    if (!data_available_) {
        gl_font(FL_HELVETICA, 12);
        glColor3f(0.4f, 0.4f, 0.4f);
        gl_draw("No analysis data", w()/2 - 50, h()/2);
        return;
    }

    drawRoughnessPairs();
    drawStabilityBar();
    drawVirtualPitchCandidates();

    // Overall score
    char buf[64];
    snprintf(buf, sizeof(buf), "Perceptual Tension: %.2f (%s)", pa_.perceptual_tension, pa_.perceptual_tension_label.c_str());
    gl_font(FL_HELVETICA_BOLD, 12);
    glColor3f(0.8f, 0.7f, 0.3f);
    gl_draw(buf, 10, 20);
}

void PsychoGraphWidget::drawRoughnessPairs() {
    int x = 10, y = 50;
    gl_font(FL_HELVETICA, 10);
    glColor3f(0.7f, 0.7f, 0.7f);
    gl_draw("Roughness Breakdown:", x, y);
    y += 15;

    float bar_max_w = w() - 100;
    int count = 0;
    for (const auto& rp : pa_.level1.roughness_per_pair) {
        if (count++ > 5) break;
        char buf[64];
        snprintf(buf, sizeof(buf), "%s - %s", rp.name_a.c_str(), rp.name_b.c_str());
        glColor3f(0.6f, 0.6f, 0.6f);
        gl_draw(buf, x, y + 10);

        float bw = rp.roughness * bar_max_w * 2.0f; // Scale for visibility
        bw = std::min(bw, bar_max_w);
        drawFilledRect(x + 80, y, bw, 12, 0.8f, 0.3f, 0.2f, 0.8f);
        y += 18;
    }
}

void PsychoGraphWidget::drawStabilityBar() {
    int x = 10, y = h() - 60;
    gl_font(FL_HELVETICA, 10);
    glColor3f(0.7f, 0.7f, 0.7f);
    gl_draw("Tonal Stability (KK):", x, y);
    y += 10;

    float bw = pa_.level3.kk_tonal_stability * (w() - 20);
    drawFilledRect(x, y, bw, 15, 0.2f, 0.6f, 0.8f, 0.8f);

    char buf[64];
    snprintf(buf, sizeof(buf), "%s (%.2f)", pa_.level3.kk_stability_label.c_str(), pa_.level3.kk_tonal_stability);
    glColor3f(1, 1, 1);
    gl_draw(buf, x + 5, y + 12);
}

void PsychoGraphWidget::drawVirtualPitchCandidates() {
    int x = 10, y = 180;
    gl_font(FL_HELVETICA, 10);
    glColor3f(0.7f, 0.7f, 0.7f);
    gl_draw("Virtual Pitch Candidates (Terhardt):", x, y);
    y += 15;

    int count = 0;
    for (const auto& c : pa_.level2.vp_top_candidates) {
        if (count++ > 3) break;
        char buf[64];
        snprintf(buf, sizeof(buf), "%s: %.1f Hz (score %.2f)", c.name.c_str(), c.f0_hz, c.score);
        float alpha = 0.3f + 0.7f * (c.score / pa_.level2.vp_top_candidates[0].score);
        glColor4f(0.3f, 0.8f, 0.4f, alpha);
        gl_draw(buf, x + 10, y);
        y += 15;
    }
}

void PsychoGraphWidget::drawFilledRect(float x, float y, float w, float h, float r, float g, float b, float a) {
    glColor4f(r, g, b, a);
    glBegin(GL_QUADS);
    glVertex2f(x, y);
    glVertex2f(x + w, y);
    glVertex2f(x + w, y + h);
    glVertex2f(x, y + h);
    glEnd();
}
