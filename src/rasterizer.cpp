#include "rasterizer.h"
#include <list>
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

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
    // sample_buffer[y * width + x] = c;
      for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&c.r)[k] * 255;
      }
  }

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

    // fill_pixel(sx, sy, color);
      for (int i = 0; i < sample_rate; ++i) {
          sample_buffer[(sy * width + sx) * sample_rate + i] = color;
      }
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
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

bool is_inside(float x, float y, float x0, float y0, float x1, float y1, float x2, float y2) {
    // vertices of the triangle  point
    Vector2D p = Vector2D(x, y);
    Vector2D v1 = Vector2D(x0, y0);
    Vector2D v2 = Vector2D(x1, y1);
    Vector2D v3 = Vector2D(x2, y2);
    // edges of the triangle
    Vector2D e1 = v2 - v1;
    Vector2D e2 = v3 - v2;
    Vector2D e3 = v1 - v3;
    // normal vectors
    Vector2D n1 = Vector2D(-e1.y, e1.x);
    Vector2D n2 = Vector2D(-e2.y, e2.x);
    Vector2D n3 = Vector2D(-e3.y, e3.x);
    // point vectors
    Vector2D p1 = p - v1;
    Vector2D p2 = p - v2;
    Vector2D p3 = p - v3;
    // compute dot product
    bool inside = dot(p1, n1) >= 0 && dot(p2, n2) >= 0 && dot(p3, n3) >= 0;
    bool outside = dot(p1, n1) <= 0 && dot(p2, n2) <= 0 && dot(p3, n3) <= 0;
    
    return (inside || outside);
}

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
      // Task 1: Implement basic triangle rasterization here, no supersampling
      // define the x-y values of bounding box for the triangles
      float minX = std::min({x0, x1, x2});
      float minY = std::min({y0, y1, y2});
      float maxX = std::max({x0, x1, x2});
      float maxY = std::max({y0, y1, y2});
      // loop through the bouding box and fill pixel if is_inside test passes
//      for (int x = minX; x < maxX; ++x) {
//          for (int y = minY; y < maxY; ++y) {
//              if (is_inside(x-0.5, y-0.5, x0, y0, x1, y1, x2, y2)) {
//                  fill_pixel(x, y, color);
//              }
//          }
//      }
    // TODO: Task 2: Update to implement super-sampled rasterization
      int sqrt_r = (int)floor(sqrt(sample_rate));
      for (int y = minY; y < maxY; ++y) {
          for (int x = minX; x < maxX; ++x) {
              for (int j = 0; j < sqrt_r; ++j) {
                  for (int i = 0; i < sqrt_r; ++i) {
                      if (is_inside(x + ((float(i) + 0.5) / sqrt_r), y + ((float(j) + 0.5) / sqrt_r), x0, y0, x1, y1, x2, y2)) {
                          sample_buffer[(y * width + x) * sample_rate + i + j * sqrt_r] = color;
                      }
                  }
              }
          }
      }
  }

Vector3D barycentric(float x, float y, float x0, float y0, float x1, float y1, float x2, float y2) {
    float alpha = (-(x-x1)*(y2-y1)+(y-y1)*(x2-x1))/(-(x0-x1)*(y2-y1)+(y0-y1)*(x2-x1));
    float beta = (-(x-x2)*(y0-y2)+(y-y2)*(x0-x2))/(-(x1-x2)*(y0-y2)+(y1-y2)*(x0-x2));
    float gamma = 1 - alpha - beta;
    return Vector3D(alpha, beta, gamma);
}


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
      float minX = std::min({x0, x1, x2});
      float minY = std::min({y0, y1, y2});
      float maxX = std::max({x0, x1, x2});
      float maxY = std::max({y0, y1, y2});
      int sqrt_r = (int)floor(sqrt(sample_rate));
      for (int y = minY; y < maxY; ++y) {
          for (int x = minX; x < maxX; ++x) {
              for (int j = 0; j < sqrt_r; ++j) {
                  for (int i = 0; i < sqrt_r; ++i) {
                      if (is_inside(x + ((float(i) + 0.5) / sqrt_r), y + ((float(j) + 0.5) / sqrt_r), x0, y0, x1, y1, x2, y2)) {
                          Vector3D bcoord = barycentric(x, y, x0, y0, x1, y1, x2, y2);
                          Color interp_color = bcoord.x * c0 + bcoord.y * c1 + bcoord.z * c2;
                          sample_buffer[(y * width + x) * sample_rate + i + j * sqrt_r] = interp_color;
                      }
                  }
              }
          }
      }
  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5:
      float minX = std::min({x0, x1, x2});
      float minY = std::min({y0, y1, y2});
      float maxX = std::max({x0, x1, x2});
      float maxY = std::max({y0, y1, y2});
      SampleParams sample;
      sample.lsm = lsm;
      sample.psm = psm;
      double data[] = {u0, u1, u2, v0, v1, v2, 1, 1, 1};
      Matrix3x3 m_texture = Matrix3x3(data);
      int sqrt_r = (int)floor(sqrt(sample_rate));
      for (int y = minY; y < maxY; ++y) {
          for (int x = minX; x < maxX; ++x) {
              for (int j = 0; j < sqrt_r; ++j) {
                  for (int i = 0; i < sqrt_r; ++i) {
                      if (is_inside(x + ((float(i) + 0.5) / sqrt_r), y + ((float(j) + 0.5) / sqrt_r), x0, y0, x1, y1, x2, y2)) {
                          Vector3D b_coord = barycentric(x, y, x0, y0, x1, y1, x2, y2);
                          Vector3D b_coord_dx = barycentric(x+1, y, x0, y0, x1, y1, x2, y2);
                          Vector3D b_coord_dy = barycentric(x, y+1, x0, y0, x1, y1, x2, y2);
                          Vector3D ts = m_texture * b_coord;
                          Vector3D ts_dx = m_texture * b_coord_dx;
                          Vector3D ts_dy = m_texture * b_coord_dy;
                          Vector2D uv = Vector2D(ts.x, ts.y);
                          Vector2D uv_dx = Vector2D(ts_dx.x, ts_dy.y);
                          Vector2D uv_dy = Vector2D(ts_dy.x, ts_dy.y);
                          // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
                          sample.p_uv = uv;
                          sample.p_dx_uv = uv_dx - uv;
                          sample.p_dy_uv = uv_dy - uv;
                          //                          if (psm == P_NEAREST) {
                          //                              c = tex.sample_nearest(uv, 0);
                          //                          } else if (psm == P_LINEAR) {
                          //                              c = tex.sample_bilinear(uv, 0);
                          //                          }
                          Color c = tex.sample(sample);
                          sample_buffer[(y * width + x) * sample_rate + i + j * sqrt_r] = c;
                      }
                  }
              }
          }
      }
  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;
    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support


    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
          Color avg = Color();
          for (int i = 0; i < sample_rate; ++i) {
              Color curSample = sample_buffer[(y * width + x) * sample_rate + i];
              avg += curSample;
          }
          avg.r /= sample_rate;
          avg.g /= sample_rate;
          avg.b /= sample_rate;
          for (int k = 0; k < 3; ++k) {
              this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&avg.r)[k] * 255;
          }
      }
    }
  }

  Rasterizer::~Rasterizer() { }


}// CGL
