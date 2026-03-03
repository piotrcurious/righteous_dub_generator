#pragma once
#include <FL/Fl_Gl_Window.H>
#include <GL/gl.h>
#include "theory_bridge.h"
#include <vector>
#include <string>

// ────────────────────────────────────────────────────────────────────────────
//  PSYCHO GRAPH WIDGET
//
//  Visualizes the internal state of the psychoacoustic model:
//  - Consonance / Roughness landscape
//  - Tonal hierarchy stability
//  - Virtual pitch candidates
// ────────────────────────────────────────────────────────────────────────────

class PsychoGraphWidget : public Fl_Gl_Window {
public:
    PsychoGraphWidget(int X, int Y, int W, int H, const char* label = nullptr);

    void update(const PsychoacousticAnalysis& pa);

    void draw() override;

private:
    PsychoacousticAnalysis pa_;
    bool data_available_{false};

    void drawRoughnessPairs();
    void drawStabilityBar();
    void drawVirtualPitchCandidates();

    void drawText(float x, float y, const char* text, float r, float g, float b);
    void drawFilledRect(float x, float y, float w, float h, float r, float g, float b, float a=1.f);
};
