#version 450

in vec2 texCoord;

uniform sampler2D txtr;
uniform vec4 vertexColor;
uniform int drawMode; // 0 = scalar quad (texture), 1 = solid overlay (glyphs, lines)

// TODO: add uniform for colormap mode (int) and blend factor (float)
uniform int colormapMode;   // 0 = original grayscale texture, 1 = rainbow, 2 = cool-warm
uniform float blendFactor;  // 0 = grayscale only, 1 = full colormap (when colormapMode != 0)

out vec4 fragColor;

// TODO: implement colormap functions
// Hint: define at least two colormaps (e.g. rainbow and cool-warm)
//   - Rainbow: interpolate between blue -> cyan -> green -> yellow -> red
//     using 4 segments at thresholds 0.25, 0.50, 0.75 with mix()
//   - Cool-warm: interpolate blue -> white -> red
//     using 2 segments at threshold 0.5 with mix()
// Full-spectrum rainbow: hue maps linearly with t (piecewise blue..red keeps G=1 for half of [0,1],
//   so mid-range scalars look "all green"; HSV avoids that).

vec3 rainbow(float t) {
    t = clamp(t, 0.0, 1.0);
    vec3 c = vec3(t, 1.0, 1.0); // hue, saturation, value
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec3 coolWarm(float t) {
    t = clamp(t, 0.0, 1.0);
    if (t < 0.5)
        return mix(vec3(0.0, 0.0, 1.0), vec3(1.0, 1.0, 1.0), t / 0.5);
    return mix(vec3(1.0, 1.0, 1.0), vec3(1.0, 0.0, 0.0), (t - 0.5) / 0.5);
}

void main() {
    // TODO: apply colormap to the scalar field texture value
    // Hint: sample the texture to get a scalar value in [0,1],
    //   use it to look up the colormap color, then blend between
    //   grayscale and the mapped color using the blend factor.
    //   When colormap mode is off (0), use the original behavior below.
    //   Make sure overlays (glyphs, streamlines) still work by
    //   setting colormapMode to 0 from C++ before drawing them.

    if (drawMode == 1) {
        fragColor = vertexColor;
        return;
    }

    vec4 texSample = texture(txtr, texCoord);
    if (colormapMode == 0) {
        // Scalar lives in .r (GL_RED upload); match old luminance look as neutral gray.
        float s0 = texSample.r;
        fragColor = vertexColor + vec4(vec3(s0), 1.0);
        return;
    }

    // Single-channel scalar: .r from GL_RED; luminance-style uploads often replicate to RGB.
    float s = max(texSample.r, max(texSample.g, texSample.b));
    vec3 gray = vec3(s);
    vec3 mapped = (colormapMode == 1) ? rainbow(s) : coolWarm(s);
    float b = clamp(blendFactor, 0.0, 1.0);
    vec3 rgb = mix(gray, mapped, b);
    fragColor = vec4(rgb, 1.0) + vertexColor;
}
