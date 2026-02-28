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
//  TONNETZ WIDGET — OpenGL rendering of the 5-limit harmonic lattice.
//
//  Layout:
//    Horizontal axis  = perfect fifths (×3/2, +7 semitones)
//    Diagonal axis    = major thirds   (×5/4, +4 semitones)
//
//  Triangles:
//    Upward   △  = major triad   (root, M3, P5)
//    Downward ▽  = minor triad   (root, m3, P5)
//
//  Interactive:
//    Left-click  node  → toggle voice on/off
//    Right-click node  → context menu (timbre, detune)
//    Middle-drag       → pan view
//    Scroll           → zoom
// ────────────────────────────────────────────────────────────────────────────

struct TonnetzNode {
    int x, y;          // lattice coordinates
    int pitch_class;   // 0-11
    float cx, cy;      // screen coords (updated each draw)
    std::string label;
};

class TonnetzWidget : public Fl_Gl_Window {
public:
    TonnetzWidget(int X, int Y, int W, int H, const char* label = nullptr);

    // ── inject live data before each draw
    void setVoices(const std::vector<Voice>& voices);
    void setAbstractObject(const AbstractObject& obj);
    void setRoughnessRecords(const std::vector<RoughnessRecord>& rr);
    void setEDO(int edo);

    // ── callback: user clicked a node → suggested note
    using NodeClickCb = std::function<void(int pitch_class, int tonnetz_x, int tonnetz_y)>;
    void setNodeClickCallback(NodeClickCb cb) { node_click_cb_ = std::move(cb); }

    // ── highlight a path (progression suggestion)
    void setProgressionPath(const std::vector<std::pair<int,int>>& path) {
        progression_path_ = path; redraw();
    }

    void setHighlightedPC(int pc) {
        highlighted_pcs_.clear();
        highlighted_node_x_ = highlighted_node_y_ = -100;
        if (pc >= 0) highlighted_pcs_.push_back(pc);
        redraw();
    }
    void setHighlightedPCs(const std::vector<int>& pcs) {
        highlighted_pcs_ = pcs;
        highlighted_node_x_ = highlighted_node_y_ = -100;
        redraw();
    }
    void setHighlightedNode(int tx, int ty) {
        highlighted_pcs_.clear();
        highlighted_node_x_ = tx;
        highlighted_node_y_ = ty;
        redraw();
    }

    const std::vector<TonnetzNode>& nodes() const { return nodes_; }

    void draw() override;
    int  handle(int event) override;

private:
    // ── lattice geometry
    static constexpr int GRID_X = 11;  // columns (fifths)
    static constexpr int GRID_Y = 7;   // rows    (thirds)
    std::vector<TonnetzNode> nodes_;

    // ── display state
    int                       edo_{12};
    std::vector<Voice>        voices_;
    AbstractObject            abs_obj_;
    std::vector<RoughnessRecord> roughness_;
    std::vector<std::pair<int,int>> progression_path_;
    std::vector<int> highlighted_pcs_;
    int              highlighted_node_x_{-100};
    int              highlighted_node_y_{-100};

    // ── camera
    float pan_x_{0.f}, pan_y_{0.f};
    float zoom_{1.f};
    int   drag_x_{0},  drag_y_{0};
    bool  dragging_{false};

    // ── node spacing
    float node_dx_{52.f}, node_dy_{46.f};

    // ── callbacks
    NodeClickCb node_click_cb_;

    // ── private drawing
    void initGL();
    void buildLattice();
    void drawBackground();
    void drawTriads();
    void drawEdges();
    void drawNodes();
    void drawAbstractObject();
    void drawProgressionPath();
    void drawRoughnessHeat();
    void drawLabels();

    void worldToScreen(float wx, float wy, float& sx, float& sy);
    void screenToWorld(int sx, int sy, float& wx, float& wy);
    TonnetzNode* nodeAt(int sx, int sy);
    void drawText(float x, float y, const char* text, float r, float g, float b);
    void drawFilledCircle(float cx, float cy, float r, float red, float g, float b, float a=1.f);
    void drawRing(float cx, float cy, float r, float thick, float red, float g, float b, float a=1.f);
    void drawTriangle(float x1,float y1, float x2,float y2, float x3,float y3,
                      float r, float g, float b, float a=0.3f);
    void drawLine(float x1,float y1, float x2,float y2, float r,float g,float b,float lw=1.f);

    int pcFromCoord(int gx, int gy) const {
        // Fixed 12-tone topology: Fifth = 7 steps, Third = 4 steps
        int pc12 = ((gx * 7 + gy * 4) % 12 + 12) % 12;
        // Map 12-tone pitch class to the nearest step in current EDO
        return (int)std::round(pc12 * edo_ / 12.0) % edo_;
    }

    bool gl_inited_{false};
};
