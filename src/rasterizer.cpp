#include "rasterizer.h"
#include <time.h>

using namespace std;

namespace CGL {

RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
                                       size_t width, size_t height,
                                       unsigned int sample_rate) {
  this->psm = psm;
  this->lsm = lsm;
  this->width = width;
  this->height = height;
  this->sample_rate = sample_rate;

  supersample_buffer.resize(width * height * sample_rate, Color::White);
}

// Used by rasterize_point and rasterize_line
void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
    for (size_t i = 0; i < sample_rate; ++i) {
        fill_supersample(x, y, i, c);
    }   
}

// Optional helper function to add a sample to the supersample_buffer
void RasterizerImp::fill_supersample(size_t x, size_t y, size_t s, Color c) {
  // TODO: Task 2: You may want to implement this function. Hint: our solution uses one line
    this->supersample_buffer[y * width * sample_rate + x * sample_rate + s] = c;
};

// Rasterize a point: simple example to help you start familiarizing
// yourself with the starter code.
//
void RasterizerImp::rasterize_point(float x, float y, Color color) {
  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= width) return;
  if (sy < 0 || sy >= height) return;

  fill_pixel(sx, sy, color);
  return;
}

// Rasterize a line.
void RasterizerImp::rasterize_line(float x0, float y0,
                                    float x1, float y1,  Color color) {
    if (x0 > x1) {
    swap(x0, x1);
    swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;

    if (steep) {
    dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
    dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
    rasterize_point(pt[0], pt[1], color);
    pt[0] += dpt[0]; pt[1] += dpt[1];
    }
}

// Rasterize a triangle.
void RasterizerImp::rasterize_triangle(float x0, float y0,
                                       float x1, float y1,
                                       float x2, float y2,
                                       Color color) {
  // TODO: Task 1: Implement basic triangle rasterization here, no supersampling

    // find x,y boundries to prevent reduntant computations
    int smallestX = (int)max(0.0f, floor(min(min(x0, x1), x2))),
        smallestY = (int)max(0.0f, floor(min(min(y0, y1), y2))),
        largestX = (int)min(width * 1.0f, ceil(max(max(x0, x1), x2))),
        largestY = (int)min(height * 1.0f, ceil(max(max(y0, y1), y2)));
 
    size_t nSamples = sqrt(sample_rate);
    size_t denominator = nSamples * 2;

    for (size_t row = smallestX; row < largestX; ++row) {
        for (size_t col = smallestY; col < largestY; ++col) {            
            for (size_t ny = 0; ny < nSamples; ++ny) {
                for (size_t nx = 0; nx < nSamples; ++nx) {
                    if (inside(x0, y0, x1, y1, x2, y2, row + (1.0 + 2.0 * ny) / denominator, col + (1.0 + 2.0 * nx) / denominator)) {
                        fill_supersample(row, col, ny * nSamples + nx, color);
                    }
                }
            }
        }
    }

  
  // TODO: Task 2: Update to implement super-sampled rasterization
  
}


void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
                                                          float x1, float y1, Color c1,
                                                          float x2, float y2, Color c2)
{
  // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    int smallestX = (int)max(0.0f, floor(min(min(x0, x1), x2))),
        smallestY = (int)max(0.0f, floor(min(min(y0, y1), y2))),
        largestX = (int)min(width * 1.0f, ceil(max(max(x0, x1), x2))),
        largestY = (int)min(height * 1.0f, ceil(max(max(y0, y1), y2)));

    size_t nSamples = sqrt(sample_rate);
    size_t denominator = nSamples * 2;
    float alpha = 0.0, beta = 0.0, epsilon = 0.0;

    for (int row = smallestX; row < largestX; ++row) {
        for (int col = smallestY; col < largestY; ++col) {
            for (size_t ny = 0; ny < nSamples; ++ny) {
                for (size_t nx = 0; nx < nSamples; ++nx) {
                    if (baycentricCoords(x0, y0, x1, y1, x2, y2, row + (1.0 + 2.0 * ny) / denominator, col + (1.0 + 2.0 * nx) / denominator, alpha, beta, epsilon)) {
                        fill_supersample(row, col, ny * nSamples + nx, alpha * c0 + beta * c1 + epsilon * c2);
                    }
                }
            }
        }
    }

  // Hint: You can reuse code from rasterize_triangle

}

void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
{
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    
    SampleParams params;
    Vector2D uv, dx, dy;
    int smallestX = (int)max(0.0f, floor(min(min(x0, x1), x2))),
        smallestY = (int)max(0.0f, floor(min(min(y0, y1), y2))),
        largestX = (int)min(width * 1.0f, ceil(max(max(x0, x1), x2))),
        largestY = (int)min(height * 1.0f, ceil(max(max(y0, y1), y2)));

    size_t nSamples = sqrt(sample_rate);
    size_t denominator = nSamples * 2;
    float alpha = 0.0, beta = 0.0, epsilon = 0.0;

    for (int row = smallestX; row < largestX; ++row) {
        for (int col = smallestY; col < largestY; ++col) {
            for (size_t ny = 0; ny < nSamples; ++ny) {
                for (size_t nx = 0; nx < nSamples; ++nx) {
                    if (baycentricCoords(x0, y0, x1, y1, x2, y2, row + (1.0 + 2.0 * ny) / denominator, col + (1.0 + 2.0 * nx) / denominator, alpha, beta, epsilon)) {
                        // find (u, v)
                        uv.x = alpha * u0 + beta * u1 + epsilon * u2;
                        uv.y = alpha * v0 + beta * v1 + epsilon * v2;
                        
                        // find (du/dx, dv/dx)
                        baycentricCoords(x0, y0, x1, y1, x2, y2, row + 1.0 + (1.0 + 2.0 * ny) / denominator, col + (1.0 + 2.0 * nx) / denominator, alpha, beta, epsilon);
                        dx.x = alpha * u0 + beta * u1 + epsilon * u2 - uv.x;
                        dx.y = alpha * v0 + beta * v1 + epsilon * v2 - uv.y;

                        // find (du/dy, dv/dy)
                        baycentricCoords(x0, y0, x1, y1, x2, y2, row + (1.0 + 2.0 * ny) / denominator, col + 1.0 + (1.0 + 2.0 * nx) / denominator, alpha, beta, epsilon);
                        dy.x = alpha * u0 + beta * u1 + epsilon * u2 - uv.x;
                        dy.y = alpha * v0 + beta * v1 + epsilon * v2 - uv.y;

                        dy *= tex.height;
                        dx *= tex.width;

                        params.p_uv = uv;
                        params.p_dx_uv = dx;
                        params.p_dy_uv = dy;
                        params.psm = psm;
                        params.lsm = lsm;
                        
                        fill_supersample(row, col, ny * nSamples + nx, tex.sample(params));
                    }
                }
            }
        }
    }

    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle

}

void RasterizerImp::set_sample_rate(unsigned int rate) {
  // TODO: Task 2: You may want to update this function for supersampling support

  this->sample_rate = rate;
  this->supersample_buffer.resize(width * height * sample_rate, Color::White);

}


void RasterizerImp::set_framebuffer_target( unsigned char* rgb_framebuffer,
                                                size_t width, size_t height )
{
  // TODO: Task 2: You may want to update this function for supersampling support

  this->width = width;
  this->height = height;
  this->rgb_framebuffer_target = rgb_framebuffer;
  
}


void RasterizerImp::clear_buffers() {
  // TODO: Task 2: You may want to update this function for supersampling support
  // Hint: With supersampling, you have an additional buffer to take care of

  std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
  this->supersample_buffer.assign(width * height * sample_rate, Color::White);
}


// This function is called at the end of rasterizing all elements of the
// SVG file.  If you use a supersample buffer to rasterize SVG elements
// for antialising, you could use this call to fill the target framebuffer
// pixels from the supersample buffer data.
//
void RasterizerImp::resolve_to_framebuffer() {
  // TODO: Task 2: You will likely want to update this function for supersampling support
    Color color(0.0, 0.0, 0.0);

    for (size_t y = 0; y < height; ++y) {
        for (size_t x = 0; x < width; ++x) {
            color = supersample_buffer[y * width * sample_rate + x * sample_rate];
            for (size_t s = 1; s < sample_rate; ++s) {
                color += supersample_buffer[y * width * sample_rate + x * sample_rate + s];
            }
            color *= (1.0 / sample_rate);
            
            rgb_framebuffer_target[3 * (y * width + x)] = (unsigned char)(color.r * 255);
            rgb_framebuffer_target[3 * (y * width + x) + 1] = (unsigned char)(color.g * 255);
            rgb_framebuffer_target[3 * (y * width + x) + 2] = (unsigned char)(color.b * 255);
        }
    }  
}

Rasterizer::~Rasterizer() { }


}// CGL
