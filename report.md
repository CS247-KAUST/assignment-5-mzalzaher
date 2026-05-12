# CS247 Assignment 5 — Report

## 1. Glyphs and related visualization

### Arrows on the slice (`drawGlyphs`, `CS247_prog.cpp`)

Glyphs are built each frame from **`vector_array`** at **`loaded_timestep`**: index **`3 * (iy * W + ix) + 3 * t * data_size`** for **`(vx, vy)`** (third component unused here). **`drawGlyphs`** loops **`ix`, `iy`** with **`stride = max(1, sampling_rate)`**; cells with **`length((vx,vy)) < eps`** are skipped. **`gridCellToNDC(ix, iy, W, H)`** places the tail anchor at the cell center in NDC, with **`fy = (iy + 0.5)/H`** and **`ndcY = -1 + 2*fy`** so row **0** aligns with the texture top. Direction uses **`normalize(vec2(vx, vy))`**; shaft runs **`tail → tip`** and two head edges **`tip → wingA/B`** with head size proportional to **`arrowLen`**.

**Assignment requirements**

- **Draw glyphs:** **`en_arrow`** toggled with **`A`**; **`render()`** sets **`drawMode == 1`**, yellow **`vertexColor`**, then **`drawGlyphs()`**. Geometry is uploaded to **`glyphVBO`** / **`glyphVAO`** (interleaved position + texCoord slots for the shared shader) and drawn with **`GL_LINES`**, **`glLineWidth(1.5f)`**.

- **Length mode (constant vs speed):** **`glyph_length_by_magnitude`** toggled with **`L`**. Constant mode: **`arrowLen = baseLen`** with **`baseLen`** from cell size (**`~6 * min(2/W, 2/H)`**, floor **0.012**). Magnitude mode: first pass on the same stride grid finds **`maxSpeed`**, then **`arrowLen = baseLen * (|v| / maxSpeed)`** per cell.

- **Subsampling:** **`sampling_rate`** controls stride; **`+` / `-`** ( **`GLFW_KEY_EQUAL` / `GLFW_KEY_MINUS`** ) change it in steps of **5**, clamped **5…100**.

---

## 2. Streamlines

### Euler vs RK2 (`CS247_prog.cpp`)

Tracing uses continuous grid coords in `[0, W-1] × [0, H-1]`, **`sampleVectorBilinear`** on **`vector_array`** at **`loaded_timestep`**, step **`dt`**, sign ±1.

- **Euler** (`integrateEulerDirection`): one bilinear sample at `p`; `pNext = p + (vx,vy) * (dt * sign)`.
- **RK2** (`integrateRK2Direction`): sample at `p` → `v1`, midpoint `pMid = p + v1 * (0.5f * dt * sign)` (must stay in-bounds), second sample at `pMid` → `v2`; `pNext = p + v2 * (dt * sign)`. Two bilinear calls per accepted step.

**Assignment requirements**

- **Keep streamlines when adding seeds:** Seeds are cleared only in **`reset_rendering_props`** when loading a new dataset (`LoadData`).

- **Recalculate when the time slice changes:** **`loadNextTimestep`** advances **`loaded_timestep`**, refreshes the scalar texture, then calls **`rebuildAllStreamlines()`**, which re-integrates every entry in **`streamline_seeds`** at the new timestep.

- **Bilinear vector sampling (no nearest-cell snap):** **`sampleVectorBilinear`** uses the cell under **`floor(gx), floor(gy)`**, fractional weights **`fx`, `fy`**, and **`glm::mix`** across the four corner vectors for each step of Euler and RK2.

- **Stopping conditions:** both integrators stop if **`length(v) < kStreamlineZeroEps`** (near-zero field), if **`accumLen > 4(W+H)`** (arc-length cap), if **`p`** or **`pNext`** leaves **`[0, W-1] × [0, H-1]`** (boundary), or after **`kStreamlineMaxSteps`** iterations. RK2 also stops if **`pMid`** leaves that box.

- **Backward and forward:** **`rebuildAllStreamlines`** runs the chosen integrator with **`directionSign`** **`-1`** and **`+1`** per seed, then concatenates **reversed backward** with **forward from index 1** so the seed is not duplicated—one polyline per seed.

---

## 3. Pathlines

### Euler vs RK2 (`CS247_prog.cpp`)

Tracing uses the same spatial grid **`[0, W-1] × [0, H-1]`** and step **`dt`** with sign ±1, but velocity is sampled in **space and time**: **`sampleVectorTrilinear`** when **`num_timesteps > 1`** (linear mix between two slices after bilinear in each slice); when **`num_timesteps <= 1`**, **`sampleVectorBilinear`** on slice **0** only and the fractional time index does not move (**`lockT`**).

- **Euler** (`integratePathEulerDirection`): one sample at **`(p, tF)`**; **`pNext = p + (vx, vy) * (dt * sign)`**; **`tNext = tF + dt * sign`** unless **`lockT`**.
- **RK2** (`integratePathRK2Direction`): sample at **`(p, tF)`** → **`v1`**, spatial midpoint **`pMid = p + v1 * (0.5f * dt * sign)`**, time midpoint **`tMid = tF + 0.5f * dt * sign`** when time is unlocked; second sample at **`(pMid, tMid)`** → **`v2`**; **`pNext = p + v2 * (dt * sign)`**, **`tNext`** advanced by **`dt * sign`** when not **`lockT`**. Two trilinear (or bilinear) calls per accepted step.

**Assignment requirements**

- **Keep pathlines when adding seeds:** **`pathline_seeds`** (`std::vector<glm::vec3>`: **`gx`, `gy`, `t0`**) is only appended in **`computePathline`** (**`t0`** = **`loaded_timestep`** at click); it is not cleared on new clicks. It is cleared in **`reset_rendering_props`** when a new dataset is loaded (**`LoadData`**).

- **Trilinear vector sampling (no snap to nearest cell or time slice):** **`sampleVectorTrilinear`** clamps **`tFloat`** to **`[0, num_timesteps - 1]`**, picks **`k0 = floor(tFloat)`** capped to **`[0, num_timesteps - 2]`**, **`alphaT = tFloat - k0`**, then **`glm::mix`** between **`sampleVectorBilinear`** at slices **`k0`** and **`k0 + 1`**.

- **Stopping conditions:** same stagnation threshold **`kStreamlineZeroEps`**, arc cap **`accumLen > 4(W+H)`**, and **`kStreamlineMaxSteps`** as streamlines; spatial exit if **`p`**, **`pNext`**, or (RK2) **`pMid`** leaves **`[0, W-1] × [0, H-1]`**; when time is active, also stop if **`tF`**, **`tNext`**, or (RK2) **`tMid`** leaves **`[0, tMax]`** with **`tMax = num_timesteps - 1`**.

- **Backward and forward:** **`rebuildAllPathlines`** uses the same **`streamline_use_rk2`** **`switch`** as streamlines, runs **`integratePath*Direction`** with **`-1`** and **`+1`**, then merges **reversed backward** and **forward from index 1** into one strip per seed and packs **`g_pathline*`** for **`drawPathlines`** / **`glMultiDrawArrays(GL_LINE_STRIP, …)`**.
