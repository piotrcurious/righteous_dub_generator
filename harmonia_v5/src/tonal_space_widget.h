#pragma once
#include <FL/Fl_Gl_Window.H>
#include <GL/gl.h>
#include <GL/glu.h>
#include "voice.h"
#include "audio_engine.h"
#include <functional>
#include <vector>
#include <string>

// ────────────────────────────────────────────────────────────────────────────
//  TONAL SPACE WIDGET — Circular Pitch Class Space
//
//  Layout:
//    Nodes are arranged in a circle representing the pitch class space.
//    Distance from center represents octave or spectral centroid.
//    Angles represent pitch classes (0 to EDO-1).
//
//  Interactive:
//    Left-click  node  → toggle voice on/off
//    Middle-drag       → pan/rotate view
//    Scroll           → zoom
// ────────────────────────────────────────────────────────────────────────────

struct TonalNode {
    int pitch_class;   // 0 to EDO-1
    int octave;        // 0-8
    float angle;       // radians
    float radius;      // distance from center
    float cx, cy;      // screen coords (updated each draw)
    std::string label;
};

class TonalSpaceWidget : public Fl_Gl_Window {
public:
    TonalSpaceWidget(int X, int Y, int W, int H, const char* label = nullptr);

    // ── inject live data before each draw
    void setVoices(const std::vector<Voice>& voices);
    void setAbstractObject(const AbstractObject& obj);
    void setRoughnessRecords(const std::vector<RoughnessRecord>& rr);
    void setEDO(int edo);

    // ── callback: user clicked a node → suggested note
    using NodeClickCb = std::function<void(int pitch_class, int tx, int ty)>;
    void setNodeClickCallback(NodeClickCb cb) { node_click_cb_ = std::move(cb); }

    void setProgressionPath(const std::vector<std::pair<int,int>>& path) {
        progression_path_ = path; redraw();
    }

    void setHighlightedPC(int pc) {
        highlighted_pcs_.clear();
        if (pc >= 0) highlighted_pcs_.push_back(pc);
        redraw();
    }
    void setHighlightedPCs(const std::vector<int>& pcs) {
        highlighted_pcs_ = pcs;
        redraw();
    }
    void setHighlightedNode(int tx, int ty) {
        // tx/ty are legacy Tonnetz coords; for circular space we map to PC
        // In this version, we'll just use PC for highlighting.
        highlighted_pcs_.clear();
        // Since we don't have a 1:1 mapping back to PC easily without EDO info here,
        // we'll rely on the app to pass PC if possible or we calculate it.
        redraw();
    }

    const std::vector<TonalNode>& nodes() const { return nodes_; }

    void draw() override;
    int  handle(int event) override;

private:
    // ── space geometry
    std::vector<TonalNode> nodes_;

    // ── display state
    int                       edo_{12};
    std::vector<Voice>        voices_;
    AbstractObject            abs_obj_;
    std::vector<RoughnessRecord> roughness_;
    std::vector<std::pair<int,int>> progression_path_;
    std::vector<int> highlighted_pcs_;

    // ── camera
    float pan_x_{0.f}, pan_y_{0.f};
    float rotation_{0.f};
    float zoom_{1.f};
    int   drag_x_{0},  drag_y_{0};
    bool  dragging_{false};

    // ── callbacks
    NodeClickCb node_click_cb_;

    // ── private drawing
    void initGL();
    void buildSpace();
    void drawBackground();
    void drawConnections();
    void drawNodes();
    void drawAbstractObject();
    void drawRoughnessHeat();
    void drawLabels();

    void worldToScreen(float wx, float wy, float& sx, float& sy);
    void screenToWorld(int sx, int sy, float& wx, float& wy);
    TonalNode* nodeAt(int sx, int sy);
    void drawText(float x, float y, const char* text, float r, float g, float b);
    void drawFilledCircle(float cx, float cy, float r, float red, float g, float b, float a=1.f);
    void drawRing(float cx, float cy, float r, float thick, float red, float g, float b, float a=1.f);
    void drawLine(float x1,float y1, float x2,float y2, float r,float g,float b,float lw=1.f);

    bool gl_inited_{false};
};
