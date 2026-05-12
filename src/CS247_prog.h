#ifndef CS247_PROG_H
#define CS247_PROG_H

#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>

// framework includes
#include "glslprogram.h"
#include "vboquad.h"


////////////////
// Structures //
////////////////

// window size
const unsigned int gWindowWidth = 512;
const unsigned int gWindowHeight = 512;

int current_scalar_field;
int data_size;
bool en_arrow;
bool en_streamline;
bool en_pathline;

int sampling_rate;
float dt;

// false = constant NDC arrow length; true = length scales with |v| (normalized by max on sampled grid)
bool glyph_length_by_magnitude;

// false = explicit Euler streamlines; true = RK2 midpoint (Heun-style 2nd order) in grid space
bool streamline_use_rk2;

// TODO: define colormap variables
// Hint: you need a colormap mode (off/rainbow/cool-warm) and a blend factor
// scalar quad: 0 = original texture (grayscale), 1 = rainbow, 2 = cool-warm (see fragment.fs)
int colormap_mode;
// 0 = all grayscale, 1 = full colormap (mix in shader when colormap_mode != 0)
float scalar_colormap_blend;




//////////////////////
//  Global defines  //
//////////////////////
#define TIMER_FREQUENCY_MILLIS  50

//////////////////////
// Global variables //
//////////////////////

// Handle of the window we're rendering to
static GLFWwindow* window;

char bmModifiers;	// keyboard modifiers (e.g. ctrl,...)

int clearColor;

// data handling
char* filenames[ 3 ];
bool grid_data_loaded;
bool scalar_data_loaded;
unsigned short vol_dim[ 3 ];
float* vector_array;
float* scalar_fields;
float* scalar_bounds;

GLuint scalar_field_texture;

GLuint glyphVAO;
GLuint glyphVBO;

GLuint streamlineVAO;
GLuint streamlineVBO;

GLuint pathlineVAO;
GLuint pathlineVBO;

// streamline seeds in continuous grid coordinates
std::vector<glm::vec2> streamline_seeds;

// pathline seeds
std::vector<glm::vec3> pathline_seeds;

int num_scalar_fields;
int num_timesteps; //stores number of time steps

int loaded_file;
int loaded_timestep;
float timestep;

int view_width, view_height; // height and width of entire view

GLuint displayList_idx;

int toggle_xy;

////////////////
// Prototypes //
////////////////

void drawGlyphs();

void rebuildAllStreamlines();
void drawStreamlines();

void rebuildAllPathlines();
void drawPathlines();

void computeStreamline(float gx, float gy);

void computePathline(float gx, float gy);

void loadNextTimestep( void );

void LoadData( char* base_filename );
void LoadVectorData( const char* filename );

void DownloadScalarFieldAsTexture( void );
void initGL( void );

void reset_rendering_props( void );

// TODO: define data arrays, VAO and VBO
// Hint: you need one for the glyphs, streamlines, pathlines

// make quad to load texture to
VBOQuad quad;

// GLSL
GLSLProgram vectorProgram;
glm::mat4 model;


#endif //CS247_PROG_H
