// CS 247 - Scientific Visualization, KAUST
//
// Programming Assignment #5
#include <cstring>
#include <algorithm>
#include <vector>
#include <glm/glm.hpp>
#include "CS247_prog.h"

using glm::vec2;
using glm::vec3;
using glm::vec4;

namespace {

vec2 gridCellToNDC(int ix, int iy, int W, int H)
{
    float fx = (static_cast<float>(ix) + 0.5f) / static_cast<float>(W); // cell center x in [0,1]
    float fy = (static_cast<float>(iy) + 0.5f) / static_cast<float>(H); // cell center y in [0,1]
    float ndcX = 2.0f * fx - 1.0f;  // [-1,1] left..right
    // glTexImage2D puts memory row y=0 at texture t=0 (bottom); quad maps window top to t=1 (iy=H-1).
    float ndcY = -1.0f + 2.0f * fy;
    return vec2(ndcX, ndcY);
}

void pushLineVertex(std::vector<float>& buf, const vec3& p)
{
    buf.push_back(p.x);
    buf.push_back(p.y);
    buf.push_back(p.z);
    buf.push_back(0.0f);
    buf.push_back(0.0f);
    buf.push_back(0.0f);
}

void pushLine(std::vector<float>& buf, const vec3& a, const vec3& b)
{
    pushLineVertex(buf, a);
    pushLineVertex(buf, b);
}

} // namespace

static std::vector<float> g_streamlineVertexData;
static std::vector<GLint> g_streamlineFirst;
static std::vector<GLsizei> g_streamlineCount;

static std::vector<float> g_pathlineVertexData;
static std::vector<GLint> g_pathlineFirst;
static std::vector<GLsizei> g_pathlineCount;

static void sampleVectorBilinear(float gx, float gy, int W, int H, int tstep, float* outVx, float* outVy)
{
    gx = std::max(0.f, std::min(static_cast<float>(W - 1), gx));
    gy = std::max(0.f, std::min(static_cast<float>(H - 1), gy));
    const int ix0 = std::min(W - 2, std::max(0, static_cast<int>(std::floor(gx))));
    const int iy0 = std::min(H - 2, std::max(0, static_cast<int>(std::floor(gy))));
    const float fx = gx - static_cast<float>(ix0);
    const float fy = gy - static_cast<float>(iy0);

    auto cellV = [&](int ix, int iy) {
        ix = std::max(0, std::min(W - 1, ix));
        iy = std::max(0, std::min(H - 1, iy));
        const int j = iy * W + ix;
        const int b = 3 * j + 3 * tstep * data_size;
        return glm::vec2(vector_array[b], vector_array[b + 1]);
    };

    const glm::vec2 v00 = cellV(ix0, iy0);
    const glm::vec2 v10 = cellV(ix0 + 1, iy0);
    const glm::vec2 v01 = cellV(ix0, iy0 + 1);
    const glm::vec2 v11 = cellV(ix0 + 1, iy0 + 1);
    const glm::vec2 vx0 = glm::mix(v00, v10, fx);
    const glm::vec2 vx1 = glm::mix(v01, v11, fx);
    const glm::vec2 v = glm::mix(vx0, vx1, fy);
    *outVx = v.x;
    *outVy = v.y;
}

// bilinear in (gx, gy) and linear between consecutive time slices (integer k and k+1).
static void sampleVectorTrilinear(float gx, float gy, float tFloat, int W, int H, int numTs, float* outVx, float* outVy)
{
    if (numTs <= 1) {
        sampleVectorBilinear(gx, gy, W, H, 0, outVx, outVy);
        return;
    }
    tFloat = std::max(0.f, std::min(static_cast<float>(numTs - 1), tFloat));
    const int k0 = std::min(numTs - 2, std::max(0, static_cast<int>(std::floor(tFloat))));
    const float alphaT = tFloat - static_cast<float>(k0);
    float vx0 = 0.f;
    float vy0 = 0.f;
    float vx1 = 0.f;
    float vy1 = 0.f;
    sampleVectorBilinear(gx, gy, W, H, k0, &vx0, &vy0);
    sampleVectorBilinear(gx, gy, W, H, k0 + 1, &vx1, &vy1);
    *outVx = glm::mix(vx0, vx1, alphaT);
    *outVy = glm::mix(vy0, vy1, alphaT);
}

static glm::vec3 gridFloatToNDC(float gx, float gy, int W, int H)
{
    const float ndcX = 2.f * (gx + 0.5f) / static_cast<float>(W) - 1.f;
    const float ndcY = -1.f + 2.f * (gy + 0.5f) / static_cast<float>(H); // same y convention as gridCellToNDC / texture
    return glm::vec3(ndcX, ndcY, 0.f);
}

static const float kStreamlineZeroEps = 1.0e-5f;
static const int kStreamlineMaxSteps = 50000;

// euler in grid space: p += (vx, vy) * dt * direction sign (grid iy matches texture after gridCellToNDC fix).
static std::vector<glm::vec2> 
integrateEulerDirection(const glm::vec2& seed, int W, int H, int tstep, float stepDt, int directionSign)
{
    // output polyline, state at seed, arc-length sum, and max length before we stop
    std::vector<glm::vec2> pts;
    glm::vec2 p = seed;
    pts.push_back(p);
    float accumLen = 0.f;
    const float maxPath = 4.f * static_cast<float>(W + H);

    // euler steps until domain exit, stagnation, length cap, or step limit
    for (int step = 0; step < kStreamlineMaxSteps; ++step) {
        if (p.x < 0.f || p.x > static_cast<float>(W - 1) || p.y < 0.f || p.y > static_cast<float>(H - 1))
            break;

        // interpolate (vx, vy) at p; bail if speed is ~zero so we do not stall in a dead region
        float vx = 0.f;
        float vy = 0.f;
        sampleVectorBilinear(p.x, p.y, W, H, tstep, &vx, &vy); // bilinear (vx, vy)
        const glm::vec2 vel(vx, vy);
        const float speed = glm::length(vel);
        if (speed < kStreamlineZeroEps)
            break;

        // euler update: move along velocity scaled by dt and integration direction (+/-)
        const glm::vec2 delta = vel * (stepDt * static_cast<float>(directionSign));
        const glm::vec2 pNext = p + delta;
        accumLen += glm::length(delta);

        // stopping: max accumulated length, or next foot outside the valid cell range
        if (accumLen > maxPath)
            break;
        if (pNext.x < 0.f || pNext.x > static_cast<float>(W - 1) || pNext.y < 0.f || pNext.y > static_cast<float>(H - 1))
            break;

        p = pNext;
        pts.push_back(p);
    }
    return pts;
}

// RK2 midpoint in grid space: full step uses velocity at p + (dt*sign/2)*V(p); same bounds/stagnation rules as euler.
static std::vector<glm::vec2>
integrateRK2Direction(const glm::vec2& seed, int W, int H, int tstep, float stepDt, int directionSign)
{
    // output polyline, state at seed, arc-length sum, and max length before we stop
    std::vector<glm::vec2> pts;
    glm::vec2 p = seed;
    pts.push_back(p);
    float accumLen = 0.f;
    const float maxPath = 4.f * static_cast<float>(W + H);
    const float halfH = 0.5f * stepDt * static_cast<float>(directionSign);
    const float fullH = stepDt * static_cast<float>(directionSign);

    // midpoint (rk2) steps until domain exit, stagnation, length cap, or step limit
    for (int step = 0; step < kStreamlineMaxSteps; ++step) {
        if (p.x < 0.f || p.x > static_cast<float>(W - 1) || p.y < 0.f || p.y > static_cast<float>(H - 1))
            break;

        // interpolate (vx, vy) at p; bail if speed is ~zero so we do not stall in a dead region
        float vx1 = 0.f;
        float vy1 = 0.f;
        sampleVectorBilinear(p.x, p.y, W, H, tstep, &vx1, &vy1); // bilinear (vx, vy)
        const glm::vec2 v1(vx1, vy1);
        const float s1 = glm::length(v1);
        if (s1 < kStreamlineZeroEps)
            break;

        // midpoint probe: half-step along v(p); bail if that point leaves the valid cell range
        const glm::vec2 pMid = p + v1 * halfH;
        if (pMid.x < 0.f || pMid.x > static_cast<float>(W - 1) || pMid.y < 0.f || pMid.y > static_cast<float>(H - 1))
            break;

        // second sample at midpoint for rk2 slope; bail if speed is ~zero at the probe
        float vx2 = 0.f;
        float vy2 = 0.f;
        sampleVectorBilinear(pMid.x, pMid.y, W, H, tstep, &vx2, &vy2); // bilinear (vx, vy)
        const glm::vec2 v2(vx2, vy2);
        const float s2 = glm::length(v2);
        if (s2 < kStreamlineZeroEps)
            break;

        // rk2 update: move along v(p_mid) scaled by full dt and integration direction (+/-)
        const glm::vec2 delta = v2 * fullH;
        const glm::vec2 pNext = p + delta;
        accumLen += glm::length(delta);

        // stopping: max accumulated length, or next foot outside the valid cell range
        if (accumLen > maxPath)
            break;
        if (pNext.x < 0.f || pNext.x > static_cast<float>(W - 1) || pNext.y < 0.f || pNext.y > static_cast<float>(H - 1))
            break;

        p = pNext;
        pts.push_back(p);
    }
    return pts;
}

static void appendVec3Interleaved(std::vector<float>& buf, const glm::vec3& p)
{
    buf.push_back(p.x);
    buf.push_back(p.y);
    buf.push_back(p.z);
    buf.push_back(0.f);
    buf.push_back(0.f);
    buf.push_back(0.f);
}

// pathline euler in grid space: p += (vx, vy)*dt*sign; velocity from trilinear(gx, gy, t) when num_ts>1 else bilinear on slice 0; t advances by dt*sign in parallel.
static std::vector<glm::vec2>
integratePathEulerDirection(const glm::vec2& seed_xy, float seedT, int W, int H, int numTs, float stepDt, int directionSign)
{
    // output polyline (spatial trace only), state (p, tF), arc-length sum, and max length before we stop
    std::vector<glm::vec2> pts;
    glm::vec2 p = seed_xy;
    const float tMax = static_cast<float>(std::max(1, numTs) - 1);
    float tF = std::max(0.f, std::min(tMax, seedT));
    pts.push_back(p);
    float accumLen = 0.f;
    const float maxPath = 4.f * static_cast<float>(W + H);
    const float dir = static_cast<float>(directionSign);
    const bool lockT = (numTs <= 1);

    // euler steps until domain exit (space or time), stagnation, length cap, or step limit
    for (int step = 0; step < kStreamlineMaxSteps; ++step) {
        if (p.x < 0.f || p.x > static_cast<float>(W - 1) || p.y < 0.f || p.y > static_cast<float>(H - 1))
            break;
        if (!lockT && (tF < 0.f || tF > tMax))
            break;

        // trilinear across slices, or bilinear on slice 0 if only one timestep
        float vx = 0.f;
        float vy = 0.f;
        if (lockT)
            sampleVectorBilinear(p.x, p.y, W, H, 0, &vx, &vy); // bilinear (vx, vy)
        else
            sampleVectorTrilinear(p.x, p.y, tF, W, H, numTs, &vx, &vy); // bilinear in space + linear in time
        const glm::vec2 vel(vx, vy);
        if (glm::length(vel) < kStreamlineZeroEps)
            break;

        // spatial step along v; advance fractional time index when multiple timesteps exist
        const glm::vec2 delta = vel * (stepDt * dir);
        const glm::vec2 pNext = p + delta;
        const float tNext = lockT ? tF : (tF + stepDt * dir);
        accumLen += glm::length(delta);

        // max accumulated length, next foot outside spatial box, or next time outside [0, tMax]
        if (accumLen > maxPath)
            break;
        if (pNext.x < 0.f || pNext.x > static_cast<float>(W - 1) || pNext.y < 0.f || pNext.y > static_cast<float>(H - 1))
            break;
        if (!lockT && (tNext < 0.f || tNext > tMax))
            break;

        p = pNext;
        tF = tNext;
        pts.push_back(p);
    }
    return pts;
}

// pathline rk2 midpoint: same (p, t) advance pattern as spatial rk2; v from trilinear at (p, tF) then at midpoint (p_mid, t_mid).
static std::vector<glm::vec2>
integratePathRK2Direction(const glm::vec2& seed_xy, float seedT, int W, int H, int numTs, float stepDt, int directionSign)
{
    // output polyline (spatial trace only), state (p, tF), arc-length sum, and max length before we stop
    std::vector<glm::vec2> pts;
    glm::vec2 p = seed_xy;
    const float tMax = static_cast<float>(std::max(1, numTs) - 1);
    float tF = std::max(0.f, std::min(tMax, seedT));
    pts.push_back(p);
    float accumLen = 0.f;
    const float maxPath = 4.f * static_cast<float>(W + H);
    const float halfH = 0.5f * stepDt * static_cast<float>(directionSign);
    const float fullH = stepDt * static_cast<float>(directionSign);
    const bool lockT = (numTs <= 1);

    // midpoint (rk2) steps until domain exit (space or time), stagnation, length cap, or step limit
    for (int step = 0; step < kStreamlineMaxSteps; ++step) {
        if (p.x < 0.f || p.x > static_cast<float>(W - 1) || p.y < 0.f || p.y > static_cast<float>(H - 1))
            break;
        if (!lockT && (tF < 0.f || tF > tMax))
            break;

        // interpolate (vx, vy) at (p, tF); bail if speed is ~zero so we do not stall in a dead region
        float vx1 = 0.f;
        float vy1 = 0.f;
        if (lockT)
            sampleVectorBilinear(p.x, p.y, W, H, 0, &vx1, &vy1); // bilinear (vx, vy)
        else
            sampleVectorTrilinear(p.x, p.y, tF, W, H, numTs, &vx1, &vy1); // bilinear in space + linear in time
        const glm::vec2 v1(vx1, vy1);
        if (glm::length(v1) < kStreamlineZeroEps)
            break;

        // midpoint probe in (p, t): half-step in space; half-step in time when num_ts > 1; bail if probe leaves valid range
        const glm::vec2 pMid = p + v1 * halfH;
        const float tMid = lockT ? tF : (tF + halfH);
        if (pMid.x < 0.f || pMid.x > static_cast<float>(W - 1) || pMid.y < 0.f || pMid.y > static_cast<float>(H - 1))
            break;
        if (!lockT && (tMid < 0.f || tMid > tMax))
            break;

        // second sample at (p_mid, t_mid) for rk2 slope; bail if speed is ~zero at the probe
        float vx2 = 0.f;
        float vy2 = 0.f;
        if (lockT)
            sampleVectorBilinear(pMid.x, pMid.y, W, H, 0, &vx2, &vy2); // bilinear (vx, vy)
        else
            sampleVectorTrilinear(pMid.x, pMid.y, tMid, W, H, numTs, &vx2, &vy2); // bilinear in space + linear in time
        const glm::vec2 v2(vx2, vy2);
        if (glm::length(v2) < kStreamlineZeroEps)
            break;

        // full spatial step along v(p_mid, t_mid); advance time by full dt*sign when multiple timesteps exist
        const glm::vec2 delta = v2 * fullH;
        const glm::vec2 pNext = p + delta;
        const float tNext = lockT ? tF : (tF + fullH);
        accumLen += glm::length(delta);

        // max accumulated length, next foot outside spatial box, or next time outside [0, tMax]
        if (accumLen > maxPath)
            break;
        if (pNext.x < 0.f || pNext.x > static_cast<float>(W - 1) || pNext.y < 0.f || pNext.y > static_cast<float>(H - 1))
            break;
        if (!lockT && (tNext < 0.f || tNext > tMax))
            break;

        p = pNext;
        tF = tNext;
        pts.push_back(p);
    }
    return pts;
}

// rebuild all streamlines for current timestep: integrate each seed forward+backward, merge to ndc strips, pack vbo for multidraw.
void rebuildAllStreamlines()
{
    //  clear batch buffers; bail if no field, no grid, or no seeds
    g_streamlineVertexData.clear();
    g_streamlineFirst.clear();
    g_streamlineCount.clear();
    if (!vector_array || !grid_data_loaded || streamline_seeds.empty())
        return;

    // grid size and time slice shared by every integration in this rebuild
    const int W = static_cast<int>(vol_dim[0]);
    const int H = static_cast<int>(vol_dim[1]);
    const int tstep = loaded_timestep;
    GLint vertexBase = 0;

    // one open line strip per seed (backward branch + forward branch, single polyline)
    for (const glm::vec2& seed : streamline_seeds) {
        std::vector<glm::vec2> backward;
        std::vector<glm::vec2> forward;
        switch (static_cast<int>(streamline_use_rk2)) {
        case 0:
            backward = integrateEulerDirection(seed, W, H, tstep, dt, -1);
            forward = integrateEulerDirection(seed, W, H, tstep, dt, +1);
            break;
        case 1:
            backward = integrateRK2Direction(seed, W, H, tstep, dt, -1);
            forward = integrateRK2Direction(seed, W, H, tstep, dt, +1);
            break;
        }
        std::vector<glm::vec3> stripNdc;

        // reverse backward (far end -> seed), then forward from first step after seed (avoids duplicate seed)
        for (int i = static_cast<int>(backward.size()) - 1; i >= 0; --i)
            stripNdc.push_back(gridFloatToNDC(backward[i].x, backward[i].y, W, H));
        for (size_t i = 1; i < forward.size(); ++i)
            stripNdc.push_back(gridFloatToNDC(forward[i].x, forward[i].y, W, H));
        if (stripNdc.size() < 2)
            continue; // skip degenerate strip (single point)

        // record first vertex index + vertex count, append interleaved verts (same layout as glyphs)
        g_streamlineFirst.push_back(vertexBase);
        g_streamlineCount.push_back(static_cast<GLsizei>(stripNdc.size()));
        for (const glm::vec3& p : stripNdc)
            appendVec3Interleaved(g_streamlineVertexData, p);
        vertexBase += static_cast<GLint>(stripNdc.size()); // next strip starts after this strip's vertex count
    }
}

// upload rebuilt streamline batches and draw each seed strip as a line strip (one multidraw call).
void drawStreamlines()
{
    // feature off or nothing to draw after rebuild
    if (!en_streamline || g_streamlineFirst.empty() || g_streamlineVertexData.empty())
        return;

    // bind streamline vao/vbo and push interleaved cpu buffer to gpu (stream = replace each frame)
    glBindVertexArray(streamlineVAO);
    glBindBuffer(GL_ARRAY_BUFFER, streamlineVBO);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(g_streamlineVertexData.size() * sizeof(float)),
                 g_streamlineVertexData.data(), GL_STREAM_DRAW);

    // one line strip per seed using first[]/count[] from rebuild (same primitive as integration polyline)
    glLineWidth(2.0f);
    glMultiDrawArrays(GL_LINE_STRIP, g_streamlineFirst.data(), g_streamlineCount.data(),
                      static_cast<GLsizei>(g_streamlineFirst.size()));

    // restore default vao/vbo bindings
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// similar to rebuildstreamlines
void rebuildAllPathlines()
{
    g_pathlineVertexData.clear();
    g_pathlineFirst.clear();
    g_pathlineCount.clear();
    if (!vector_array || !grid_data_loaded || pathline_seeds.empty())
        return;

    const int W = static_cast<int>(vol_dim[0]);
    const int H = static_cast<int>(vol_dim[1]);
    const int numTs = num_timesteps;
    GLint vertexBase = 0;

    for (const glm::vec3& seed : pathline_seeds) {
        const glm::vec2 xy(seed.x, seed.y);
        const float t0 = seed.z;
        std::vector<glm::vec2> backward;
        std::vector<glm::vec2> forward;
        switch (static_cast<int>(streamline_use_rk2)) {
        case 0:
            backward = integratePathEulerDirection(xy, t0, W, H, numTs, dt, -1);
            forward = integratePathEulerDirection(xy, t0, W, H, numTs, dt, +1);
            break;
        case 1:
            backward = integratePathRK2Direction(xy, t0, W, H, numTs, dt, -1);
            forward = integratePathRK2Direction(xy, t0, W, H, numTs, dt, +1);
            break;
        }
        std::vector<glm::vec3> stripNdc;
        for (int i = static_cast<int>(backward.size()) - 1; i >= 0; --i)
            stripNdc.push_back(gridFloatToNDC(backward[i].x, backward[i].y, W, H));
        for (size_t i = 1; i < forward.size(); ++i)
            stripNdc.push_back(gridFloatToNDC(forward[i].x, forward[i].y, W, H));
        if (stripNdc.size() < 2)
            continue;

        g_pathlineFirst.push_back(vertexBase);
        g_pathlineCount.push_back(static_cast<GLsizei>(stripNdc.size()));
        for (const glm::vec3& p : stripNdc)
            appendVec3Interleaved(g_pathlineVertexData, p);
        vertexBase += static_cast<GLint>(stripNdc.size());
    }
}

// similar to drawStreamlines with different toggles and buffers
void drawPathlines()
{
    if (!en_pathline || g_pathlineFirst.empty() || g_pathlineVertexData.empty())
        return;

    glBindVertexArray(pathlineVAO);
    glBindBuffer(GL_ARRAY_BUFFER, pathlineVBO);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(g_pathlineVertexData.size() * sizeof(float)),
                 g_pathlineVertexData.data(), GL_STREAM_DRAW);
    glLineWidth(2.0f);
    glMultiDrawArrays(GL_LINE_STRIP, g_pathlineFirst.data(), g_pathlineCount.data(),
                      static_cast<GLsizei>(g_pathlineFirst.size()));
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}


// cycle clear colors
static void nextClearColor()
{
    clearColor = (++clearColor) % 3;

    switch(clearColor)
    {
        case 0:
            glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
            break;
        case 1:
            glClearColor(0.2f, 0.2f, 0.3f, 1.0f);
            break;
        default:
            glClearColor(0.7f, 0.7f, 0.7f, 1.0f);
            break;
    }
}


// callbacks
// framebuffer to fix viewport
void frameBufferCallback(GLFWwindow* window, int width, int height)
{
    view_width = width;
    view_height = height;
    glViewport(0, 0, width, height);
}

// key callback to take user inputs for both windows
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action != GLFW_RELEASE) {
        char* status[ 2 ];
        status[ 0 ] = "disabled";
        status[ 1 ] = "enabled";

        switch (key) {
            case '1':
                toggle_xy = 0;
                LoadData( filenames[ 0 ] );
                loaded_file = 0;
                fprintf( stderr, "Loading " );
                fprintf( stderr, filenames[ 0 ] );
                fprintf( stderr, " dataset.\n");
                break;
            case '2':
                toggle_xy = 0;
                LoadData(filenames[ 1 ] );
                loaded_file = 1;
                fprintf( stderr, "Loading " );
                fprintf( stderr, filenames[ 1 ] );
                fprintf( stderr, " dataset.\n");
                break;
            case '3':
                toggle_xy = 1;
                LoadData( filenames[ 2 ] );
                loaded_file = 2;
                fprintf( stderr, "Loading " );
                fprintf( stderr, filenames[ 2 ] );
                fprintf( stderr, " dataset.\n");
                break;
            case '0':
                if( num_timesteps > 1 ){
                    loadNextTimestep();
                    fprintf( stderr, "Timestep %d.\n", loaded_timestep );
                }
                break;
            case GLFW_KEY_A:
                en_arrow = !en_arrow;
                fprintf(stderr, "%s drawing arrows.\n", en_arrow? "enabling" : "disabling");
                break;
            case GLFW_KEY_L:
                glyph_length_by_magnitude = !glyph_length_by_magnitude;
                fprintf(stderr, "Glyph length: %s\n", glyph_length_by_magnitude ? "proportional to speed" : "constant");
                break;
            case GLFW_KEY_S:
                current_scalar_field = (current_scalar_field + 1)%num_scalar_fields;
                DownloadScalarFieldAsTexture();
                fprintf( stderr, "Scalar field changed.\n");
                break;
            case GLFW_KEY_B:
                nextClearColor();
                fprintf( stderr, "Next clear color.\n");
                break;
            case GLFW_KEY_EQUAL:
                sampling_rate = std::min(sampling_rate + 5, 100);
                fprintf(stderr, "Increasing sampling rate to %d.\n", sampling_rate);
                break;
            case GLFW_KEY_MINUS:
                sampling_rate = std::max(sampling_rate - 5, 5);
                fprintf(stderr, "Decreasing sampling rate to: %d.\n", sampling_rate);
                break;
            case GLFW_KEY_I:
                dt = std::min(dt + 0.005, 1.);
                fprintf(stderr, "Increase dt: %.2f\n", dt);
                if (!streamline_seeds.empty())
                    rebuildAllStreamlines();
                if (!pathline_seeds.empty())
                    rebuildAllPathlines();
                break;
            case GLFW_KEY_K:
                dt = std::max(dt - 0.005, 0.0001);
                fprintf(stderr, "Decrease dt: %.2f\n", dt);
                if (!streamline_seeds.empty())
                    rebuildAllStreamlines();
                if (!pathline_seeds.empty())
                    rebuildAllPathlines();
                break;
            case GLFW_KEY_T:
                en_streamline = !en_streamline;
                fprintf(stderr, "%s drawing streamlines.\n", en_streamline? "enabling" : "disabling");
                break;
            case GLFW_KEY_P:
                en_pathline = !en_pathline;
                fprintf(stderr, "%s drawing pathlines.\n", en_pathline? "enabling" : "disabling");
                break;
            case GLFW_KEY_R:
                streamline_use_rk2 = !streamline_use_rk2;
                fprintf(stderr, "Streamline/pathline integration: %s\n", streamline_use_rk2 ? "RK2 (midpoint)" : "Euler");
                if (!streamline_seeds.empty())
                    rebuildAllStreamlines();
                if (!pathline_seeds.empty())
                    rebuildAllPathlines();
                break;
            // TODO: add keyboard controls for:
            //   - toggle colormap mode (cycle off/rainbow/cool-warm)
            //   - adjust blend factor (increase/decrease between 0.0 and 1.0)
            //   - clear all seeds
            case GLFW_KEY_C:
                colormap_mode = (colormap_mode + 1) % 3;
                fprintf(stderr, "Scalar colormap: %s\n",
                        colormap_mode == 0 ? "off (grayscale texture)" :
                        (colormap_mode == 1 ? "rainbow" : "cool-warm"));
                break;
            case GLFW_KEY_LEFT_BRACKET:
                scalar_colormap_blend = std::max(0.f, scalar_colormap_blend - 0.05f);
                fprintf(stderr, "Colormap blend (0=gray, 1=color): %.2f\n", scalar_colormap_blend);
                break;
            case GLFW_KEY_RIGHT_BRACKET:
                scalar_colormap_blend = std::min(1.f, scalar_colormap_blend + 0.05f);
                fprintf(stderr, "Colormap blend (0=gray, 1=color): %.2f\n", scalar_colormap_blend);
                break;
            case GLFW_KEY_X:
                streamline_seeds.clear();
                pathline_seeds.clear();
                rebuildAllStreamlines();
                rebuildAllPathlines();
                fprintf(stderr, "Cleared all streamline and pathline seeds.\n");
                break;
            case GLFW_KEY_Q:
            case GLFW_KEY_ESCAPE:
                exit( 0 );
                break;
            default:
                fprintf( stderr, "\nKeyboard commands:\n\n"
                                 "1, load %s dataset\n"
                                 "2, load %s dataset\n"
                                 "3, load %s dataset\n"
                                 "0, cycle through timesteps\n"
                                 "b, switch backgropund color\n"
                                 "a, en-/disable arrows.\n"
                                 "l, toggle arrow length: constant vs proportional to speed.\n"
                                 "t, en-/disable streamlines.\n"
                                 "p, en-/disable pathlines.\n"
                                 "r, toggle streamline/pathline integration: Euler vs RK2 (midpoint).\n"
                                 "+, increase sampling rate.\n"
                                 "-, decrease sampling rate.\n"
                                 "i, increase dt.\n"
                                 "k, decrease dt.\n"
                                 "c, cycle scalar colormap: off / rainbow / cool-warm.\n"
                                 "[ / ], decrease / increase grayscale-vs-colormap blend (0..1).\n"
                                 "x, clear all streamline and pathline seeds.\n"
                                 "q, <esc> - Quit\n",
                         filenames[0], filenames[1], filenames[2] );
                break;
        }
    }
}

// mouse callback to seed streamlines/pathlines
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
    if(button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        double xpos, ypos;
        //getting cursor position
        glfwGetCursorPos(window, &xpos, &ypos);
        // TODO: seed streamlines & pathlines using mouse clicks
        // Hint: convert screen coords to grid coords (y-flip needed),
        //       then call computeStreamline/computePathline when enabled
        if (grid_data_loaded && vector_array) {
            int fbw = 0;
            int fbh = 0;
            glfwGetFramebufferSize(window, &fbw, &fbh);
            if (fbw < 1)
                fbw = 1;
            if (fbh < 1)
                fbh = 1;
            const int W = static_cast<int>(vol_dim[0]);
            const int H = static_cast<int>(vol_dim[1]);
            float gx = static_cast<float>(xpos / static_cast<double>(fbw)) * static_cast<float>(W - 1);
            float gy = (1.f - static_cast<float>(ypos / static_cast<double>(fbh))) * static_cast<float>(H - 1);
            gx = std::max(0.f, std::min(static_cast<float>(W - 1), gx));
            gy = std::max(0.f, std::min(static_cast<float>(H - 1), gy));
            if (en_streamline)
                computeStreamline(gx, gy);
            if (en_pathline)
                computePathline(gx, gy);
        }
        // TODO: seed pathlines when en_pathline
    }
}

// glfw error callback
static void errorCallback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}

// data

void loadNextTimestep( void )
{
    loaded_timestep = ( loaded_timestep + 1 ) % num_timesteps;
    DownloadScalarFieldAsTexture();
    rebuildAllStreamlines();
}


/*
 * load .gri dataset
 * This only reads the header information and calls the dat loader
 * For now we ignore the grid data and assume a rectangular grid
 */
void LoadData( char* base_filename )
{
    //reset
    reset_rendering_props();

    char filename[ 80 ];
    strcpy( filename, base_filename );
    strcat( filename, ".gri");

    fprintf( stderr, "loading grid file %s\n", filename );

    // open grid file, read only, binary mode
    FILE* fp = fopen( filename, "rb" );
    if ( fp == NULL ) {
        fprintf( stderr, "Cannot open file %s for reading.\n", filename );
        return;
    }

    // read header
    char header[ 40 ];
    fread( header, sizeof( char ), 40, fp );
    sscanf( header, "SN4DB %d %d %d %d %d %f",
            &vol_dim[ 0 ], &vol_dim[ 1 ], &vol_dim[ 2 ],
            &num_scalar_fields, &num_timesteps ,&timestep );

    fprintf( stderr, "dimensions: x: %d, y: %d, z: %d.\n", vol_dim[ 0 ], vol_dim[ 1 ], vol_dim[ 2 ] );
    fprintf( stderr, "additional info: # scalar fields: %d, # timesteps: %d, timestep: %f.\n", num_scalar_fields, num_timesteps, timestep );

    // read data
    char dat_filename[ 80 ];
    strcpy( dat_filename, base_filename );

    if( num_timesteps <= 1 ){

        strcat( dat_filename, ".dat");

    } else {

        strcat( dat_filename, ".00000.dat");

    }

    loaded_timestep = 0;
    LoadVectorData( base_filename );

    glfwSetWindowSize(window, vol_dim[ 0 ], vol_dim[ 1 ] );
    grid_data_loaded = true;
}

/*
 * load .dat dataset
 * loads vector and scalar fields
 */
void LoadVectorData( const char* filename )
{
    fprintf( stderr, "loading scalar file %s\n", filename );

    // open data file, read only, binary mode
    char ts_name[ 80 ];
    if( num_timesteps > 1 )
    {
        sprintf( ts_name, "%s.%.5d.dat", filename, 0 );
    }
    else
        sprintf( ts_name, "%s.dat",filename);

    FILE* fp = fopen( ts_name, "rb" );
    if ( fp == NULL ) {
        fprintf( stderr, "Cannot open file %s for reading.\n", filename );
        return;
    }
    else
    {
        fclose( fp );
    }

    data_size = vol_dim[ 0 ] * vol_dim[ 1 ] * vol_dim[ 2 ];

    if (!vector_array) {
        delete[] vector_array;
        vector_array = NULL;
    }
    // dim.xyz * vector.xyz * timesteps
    vector_array = new float[ data_size * 3 * num_timesteps];

    // read data
    if (!scalar_fields) {
        delete[] scalar_fields;
        scalar_fields = NULL;
        delete[] scalar_bounds;
        scalar_bounds = NULL;
    }
    // dim.xyz * scalarfields(2) * timesteps
    scalar_fields = new float[ data_size * num_scalar_fields * num_timesteps ];
    scalar_bounds = new float[ 2 * num_scalar_fields * num_timesteps ];

    int num_total_fields = num_scalar_fields + 3; // scalar fields + vec.xyz
    float *tmp = new float[ data_size * num_total_fields * num_timesteps ];

    for( int k = 0 ; k < num_timesteps; k++ )
    {
        char times_name[ 80 ];
        if( num_timesteps > 1 )
        {
            sprintf( times_name, "%s.%.5d.dat", filename, k );
            fp = fopen( times_name, "rb" );
        }
        else
            fp = fopen( ts_name, "rb" );
        // read scalar data
        fread( &tmp[k*data_size*num_total_fields], sizeof( float ), ( data_size * num_total_fields ), fp );

        // close file
        fclose( fp );

        // copy and scan for min and max values
        for( int  i = 0; i < num_scalar_fields; i++ ){

            float min_val = 99999.9f;
            float max_val = 0.0f;

            float avg = 0.0;

            int offset = i * data_size * num_timesteps;

            for( int j = 0; j < data_size; j++ ){

                float val = tmp[ j * num_total_fields + 3 + i + k*data_size*num_total_fields ];

                scalar_fields[ j + k*data_size + offset ] = val;

                if( toggle_xy ) {
                    // overwrite
                    if( i == 0 ){
                        vector_array[ 3*j + 0 + 3*k*data_size ] = tmp[ j * num_total_fields + 1 + k*data_size*num_total_fields ];//toggle x and y components in vector field
                        vector_array[ 3*j + 1 + 3*k*data_size ] = tmp[ j * num_total_fields + 0 + k*data_size*num_total_fields ];
                        vector_array[ 3*j + 2 + 3*k*data_size ] = tmp[ j * num_total_fields + 2 + k*data_size*num_total_fields ];
                    }
                } else {
                    // overwrite
                    if( i == 0 ){
                        vector_array[ 3*j + 0 + 3*k*data_size ] = tmp[ j * num_total_fields + 0 + k*data_size*num_total_fields ];
                        vector_array[ 3*j + 1 + 3*k*data_size ] = tmp[ j * num_total_fields + 1 + k*data_size*num_total_fields ];
                        vector_array[ 3*j + 2 + 3*k*data_size ] = tmp[ j * num_total_fields + 2 + k*data_size*num_total_fields ];
                    }
                }

                min_val = std::min( val, min_val );
                max_val = std::max( val, max_val );

                avg += scalar_fields[ offset + j + k*data_size ] / data_size;
            }
            scalar_bounds[ 2 * i     + k*num_scalar_fields*2 ] = min_val;
            scalar_bounds[ 2 * i + 1 + k*num_scalar_fields*2 ] = max_val;
        }

        // normalize
        for( int  i = 0; i < num_scalar_fields; i++ ){

            int offset = i * data_size * num_timesteps;

            float lower_bound = scalar_bounds[ 2 * i     + k*num_scalar_fields*2 ];
            float upper_bound = scalar_bounds[ 2 * i + 1 + k*num_scalar_fields*2 ];

            // scale between [0..1] where 1 is original zero
            // the boundary of the bigger abs border will be used to scale
            // meaning one boundary will likely not be hit i.e real scale
            // for -50..100 will be [0.25..1.0]
            if( lower_bound < 0.0 && upper_bound > 0.0 ){

                float scale = 0.5f / std::max( -lower_bound, upper_bound );

                for( int j = 0; j < data_size; j++ ){

                    scalar_fields[ offset + j + k*data_size ] = 0.5f + scalar_fields[ offset + j + k*data_size ] * scale;
                }
                scalar_bounds[ 2 * i     + k*num_scalar_fields*2 ] = 0.5f + scalar_bounds[ 2 * i     + k*num_scalar_fields*2 ] * scale;
                scalar_bounds[ 2 * i + 1 + k*num_scalar_fields*2 ] = 0.5f + scalar_bounds[ 2 * i + 1 + k*num_scalar_fields*2 ] * scale;


                // scale between [0..1]
            } else {

                float sign = upper_bound <= 0.0 ? -1.0f : 1.0f;

                float scale = 1.0f / ( upper_bound - lower_bound ) * sign;

                for( int j = 0; j < data_size; j++ ){

                    scalar_fields[ offset + j + k*data_size ] = ( scalar_fields[ offset + j + k*data_size ] - lower_bound ) * scale;
                }
                scalar_bounds[ 2 * i     + k*num_scalar_fields*2 ] = ( scalar_bounds[ 2 * i     + k*num_scalar_fields*2 ] + lower_bound ) * scale; //should be 0.0
                scalar_bounds[ 2 * i + 1 + k*num_scalar_fields*2 ] = ( scalar_bounds[ 2 * i + 1 + k*num_scalar_fields*2 ] + lower_bound ) * scale; //should be 1.0
            }
        }
    }
    delete[] tmp;
    DownloadScalarFieldAsTexture();

    scalar_data_loaded = true;
}


void DownloadScalarFieldAsTexture( void )
{
    fprintf( stderr, "downloading scalar field to 2D texture\n" );

    glEnable( GL_TEXTURE_2D );

    // One GL texture name for the lifetime of the app; refresh pixels with glTexImage2D.
    if (scalar_field_texture == 0) {
        glGenTextures(1, &scalar_field_texture);
    }
    glBindTexture(GL_TEXTURE_2D, scalar_field_texture);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    int datasize = vol_dim[0] * vol_dim[1];
    // GL_RED / GL_R32F: reliable single-channel float for modern GL (legacy GL_LUMINANCE can mis-sample).
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, vol_dim[0], vol_dim[1], 0, GL_RED, GL_FLOAT,
                 &scalar_fields[(loaded_timestep + current_scalar_field * num_timesteps) * datasize]);


    glDisable( GL_TEXTURE_2D );
}

char *textFileRead( char *fn ){

    FILE *fp;
    char *content = NULL;

    int count=0;

    if (!fn) {
        fp = fopen(fn,"rt");

        if (!fp) {

            fseek(fp, 0, SEEK_END);
            count = ftell(fp);
            rewind(fp);

            if (count > 0) {
                content = (char *)malloc(sizeof(char) * (count+1));
                count = fread(content,sizeof(char),count,fp);
                content[count] = '\0';
            }
            fclose(fp);
        }
    }
    return content;
}


// initializations
// init application
bool initApplication(int argc, char **argv)
{

    std::string version((const char *)glGetString(GL_VERSION));
    std::stringstream stream(version);
    unsigned major, minor;
    char dot;

    stream >> major >> dot >> minor;

    assert(dot == '.');
    if (major > 3 || (major == 2 && minor >= 0)) {
        std::cout << "OpenGL Version " << major << "." << minor << std::endl;
    } else {
        std::cout << "The minimum required OpenGL version is not supported on this machine. Supported is only " << major << "." << minor << std::endl;
        return false;
    }

    return true;
}

void reset_rendering_props( void )
{
    num_scalar_fields = 0;
    streamline_seeds.clear();
    g_streamlineVertexData.clear();
    g_streamlineFirst.clear();
    g_streamlineCount.clear();
    pathline_seeds.clear();
    g_pathlineVertexData.clear();
    g_pathlineFirst.clear();
    g_pathlineCount.clear();
}

// set up the scene
void setup() {
    LoadData( filenames[ 0 ] );
    loaded_file = 0;

    DownloadScalarFieldAsTexture();


    // compile & link shader 
    vectorProgram.compileShader("../../shaders/vertex.vs");
    vectorProgram.compileShader("../../shaders/fragment.fs");
    vectorProgram.link();

    // make quad to render texture
    // see: vboquad.h and vboquad.cpp
    quad.init();

    // TODO: glyph/streamlines/pathlines VAO and VBO
    glGenVertexArrays(1, &glyphVAO);
    glGenBuffers(1, &glyphVBO);
    glBindVertexArray(glyphVAO);
    glBindBuffer(GL_ARRAY_BUFFER, glyphVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 4, nullptr, GL_STREAM_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glGenVertexArrays(1, &streamlineVAO);
    glGenBuffers(1, &streamlineVBO);
    glBindVertexArray(streamlineVAO);
    glBindBuffer(GL_ARRAY_BUFFER, streamlineVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 4, nullptr, GL_STREAM_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glGenVertexArrays(1, &pathlineVAO);
    glGenBuffers(1, &pathlineVBO);
    glBindVertexArray(pathlineVAO);
    glBindBuffer(GL_ARRAY_BUFFER, pathlineVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 4, nullptr, GL_STREAM_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// rendering
void render() {
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glEnable( GL_TEXTURE_2D );

    // draw the texture
    glBindTexture(GL_TEXTURE_2D, scalar_field_texture);
    vectorProgram.use();

    model = mat4(1);

    vectorProgram.setUniform("drawMode", 0); // add draw mode
    vectorProgram.setUniform("vertexColor", glm::vec4(0));
    vectorProgram.setUniform("model", model);

    // TODO: pass colormap uniforms to shader before drawing the quad
    // Hint: set colormapMode and blendFactor uniforms here
    vectorProgram.setUniform("colormapMode", colormap_mode);
    vectorProgram.setUniform("blendFactor", scalar_colormap_blend);

    quad.render();
    glDisable( GL_TEXTURE_2D );

    // TODO: reset colormap mode to 0 before drawing overlays
    // so that glyphs/streamlines/pathlines use solid colors
    // (Overlays use drawMode==1; fragment.fs ignores colormap for drawMode==1.)

    // TODO: draw glyphs, streamlines, pathlines
    if (en_arrow && grid_data_loaded && vector_array) {
        vectorProgram.setUniform("drawMode", 1);
        vectorProgram.setUniform("vertexColor", vec4(1.0f, 0.95f, 0.15f, 1.0f));
        vectorProgram.setUniform("model", model);
        drawGlyphs();
    }

    if (en_streamline && grid_data_loaded && vector_array) {
        vectorProgram.setUniform("drawMode", 1);
        vectorProgram.setUniform("vertexColor", vec4(0.25f, 0.95f, 1.0f, 1.0f));
        vectorProgram.setUniform("model", model);
        drawStreamlines();
    }

    if (en_pathline && grid_data_loaded && vector_array) {
        vectorProgram.setUniform("drawMode", 1);
        vectorProgram.setUniform("vertexColor", vec4(1.0f, 0.45f, 0.85f, 1.0f));
        vectorProgram.setUniform("model", model);
        drawPathlines();
    }

}

// entry point
int main(int argc, char** argv)
{
    // init variables
    view_width = 0;
    view_height = 0;

    toggle_xy = 0;

    en_arrow = false;
    en_streamline = false;
    en_pathline = false;
    sampling_rate = 15;
    dt = 0.1;
    glyph_length_by_magnitude = false;
    streamline_use_rk2 = false;
    colormap_mode = 0;
    scalar_colormap_blend = 1.0f;

    reset_rendering_props();

    vector_array = NULL;
    scalar_fields = NULL;
    scalar_bounds = NULL;
    grid_data_loaded = false;
    scalar_data_loaded = false;
    current_scalar_field = 0;
    clearColor = 0;

    filenames[0] = (char*)"../../data/block/c_block";
    filenames[1] = (char*)"../../data/tube/tube";
    filenames[2] = (char*)"../../data/hurricane/hurricane_p_tc";



    // set glfw error callback
    glfwSetErrorCallback(errorCallback);

    // init glfw
    if (!glfwInit()) {
        exit(EXIT_FAILURE);
    }

    // init glfw window
    window = glfwCreateWindow(gWindowWidth, gWindowHeight, "AMCS/CS247 Scientific Visualization", nullptr, nullptr);
    if (!window)
    {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    // set GLFW callback functions
    glfwSetKeyCallback(window, keyCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetFramebufferSizeCallback(window, frameBufferCallback);

    // make context current (once is sufficient)
    glfwMakeContextCurrent(window);

    // get the frame buffer size
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    // init the OpenGL API (we need to do this once before any calls to the OpenGL API)
    gladLoadGL();

    // init our application
    if (!initApplication(argc, argv)) {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }


    // set up the scene
    setup();

    // print menu
    keyCallback(window, GLFW_KEY_BACKSLASH, 0, GLFW_PRESS, 0);

    // start traversing the main loop
    // loop until the user closes the window
    while (!glfwWindowShouldClose(window))
    {
        // clear frame buffer
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // render one frame
        render();

        // swap front and back buffers
        glfwSwapBuffers(window);

        // poll and process input events (keyboard, mouse, window, ...)
        glfwPollEvents();
    }

    glfwTerminate();
    return EXIT_SUCCESS;
}

// TODO: define any useful functions you might need.
//  e.g., indexing, linear interpolation ..etc


void computeStreamline(float gx, float gy)
{
    // TODO: compute streamlines starting from x,y position. enable switching between euler and runge kutta
    // Hint: implement bilinear interpolation of vectors, forward+backward integration,
    //       and stopping conditions (boundary, zero vector, max accumulated length)

    // TODO: set any useful uniforms & update VBO & draw

    if (!grid_data_loaded || !vector_array)
        return;
    if (!en_streamline)
        return;

    fprintf(stderr, "loading streamline...\n");
    fflush(stderr);

    streamline_seeds.push_back(glm::vec2(gx, gy));
    rebuildAllStreamlines();

    fprintf(stderr, "done streamline.\n");
    fflush(stderr);
}

void computePathline(float gx, float gy)
{
    if (!grid_data_loaded || !vector_array)
        return;
    if (!en_pathline)
        return;

    fprintf(stderr, "loading pathline...\n");
    fflush(stderr);

    pathline_seeds.push_back(glm::vec3(gx, gy, static_cast<float>(loaded_timestep)));
    rebuildAllPathlines();

    fprintf(stderr, "done pathline.\n");
    fflush(stderr);
}

void drawGlyphs()
{
    // TODO: draw arrows/glyphs
    // Hint: iterate over grid with sampling_rate stride, compute arrow geometry
    //       (shaft + arrowhead) in NDC, upload to VBO, draw with GL_LINES

    if (!vector_array || !grid_data_loaded)
        return;

    const int W = static_cast<int>(vol_dim[0]);
    const int H = static_cast<int>(vol_dim[1]);
    if (W <= 0 || H <= 0)
        return;

    const int stride = std::max(1, sampling_rate); // subsample: one glyph every stride cells
    const int t = loaded_timestep;
    const float eps = 1.0e-8f; // treat smaller |v| as zero (skip / avoid div-by-zero)

    const float cell = std::min(2.0f / static_cast<float>(W), 2.0f / static_cast<float>(H)); // NDC size of ~one cell
    const float baseLen = std::max(6.0f * cell, 0.012f); // constant-mode length; max length in magnitude mode

    float maxSpeed = 0.0f;
    if (glyph_length_by_magnitude) {
        for (int iy = 0; iy < H; iy += stride) {
            for (int ix = 0; ix < W; ix += stride) {
                const int j = iy * W + ix;
                const int base = 3 * j + 3 * t * data_size; // index into vector_array (xyz per cell, per timestep)
                const float vx = vector_array[base + 0];
                const float vy = vector_array[base + 1];
                const float s = glm::length(vec2(vx, vy));
                if (s > maxSpeed)
                    maxSpeed = s; // normalize lengths by fastest vector on this subsampled set
            }
        }
        if (maxSpeed < eps)
            return;
    }

    std::vector<float> verts;
    verts.reserve(static_cast<size_t>((W / stride + 2) * (H / stride + 2) * 6 * 6));

    for (int iy = 0; iy < H; iy += stride) {
        for (int ix = 0; ix < W; ix += stride) {
            const int j = iy * W + ix;
            const int base = 3 * j + 3 * t * data_size;
            const float vx = vector_array[base + 0];
            const float vy = vector_array[base + 1];

            vec2 dir(vx, vy);
            const float len = glm::length(dir); // scalar speed at this cell
            if (len < eps)
                continue;
            dir /= len; // unit direction in NDC-consistent 2D

            // set arrowlen based on current mode
            float arrowLen = baseLen;
            if (glyph_length_by_magnitude) {
                const float speedFraction = len / maxSpeed; // 1 = fastest on subsampled grid, <1 = slower
                arrowLen = baseLen * speedFraction;
            }

            const float headLen = 0.35f * arrowLen;   // arrowhead depth along -dir from tip
            const float headHalfW = 0.22f * arrowLen; // half-width of head perpendicular to dir

            const vec2 c = gridCellToNDC(ix, iy, W, H);
            const vec3 tip(c.x + dir.x * arrowLen, c.y + dir.y * arrowLen, 0.0f);
            const vec3 tail(c.x - dir.x * arrowLen, c.y - dir.y * arrowLen, 0.0f);
            const vec3 headBase(tip.x - dir.x * headLen, tip.y - dir.y * headLen, 0.0f); // inward from tip along shaft
            const vec2 perp(-dir.y, dir.x); // unit normal in 2D (dir already unit)
            const vec3 wingA(headBase.x + perp.x * headHalfW, headBase.y + perp.y * headHalfW, 0.0f);
            const vec3 wingB(headBase.x - perp.x * headHalfW, headBase.y - perp.y * headHalfW, 0.0f);

            pushLine(verts, tail, tip);
            pushLine(verts, tip, wingA);
            pushLine(verts, tip, wingB);
        }
    }

    if (verts.empty())
        return;

    glBindVertexArray(glyphVAO);
    glBindBuffer(GL_ARRAY_BUFFER, glyphVBO);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(verts.size() * sizeof(float)), verts.data(), GL_STREAM_DRAW);

    glLineWidth(1.5f);
    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(verts.size() / 6));
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}
