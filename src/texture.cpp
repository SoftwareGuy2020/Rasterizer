#include "texture.h"
#include "CGL/color.h"

#include <algorithm>
#include <cmath>

namespace CGL {

Color Texture::sample(const SampleParams &sp) {
    
  // TODO: Task 6: Fill this in.
  // return magenta for invalid level
    float lev = get_level(sp);
    if (lev < 0.0) {
        lev = 0.0;
    }
    if (sp.psm) {
        if (sp.lsm == LevelSampleMethod::L_LINEAR) {
            float alpha = min((size_t)lev, mipmap.size() - 1), beta = min((size_t)lev + 1, mipmap.size() - 1);
            
            Color d1 = sample_bilinear(sp.p_uv, alpha),
                d2 = sample_bilinear(sp.p_uv, beta);

            Color c = (lev - alpha) * d1 + (beta - lev) * d2;
            c.r = min(1.0f, max(0.0f, c.r));
            c.g = min(1.0f, max(0.0f, c.g));
            c.b = min(1.0f, max(0.0f, c.b));
            return c;

        }
        return sample_bilinear(sp.p_uv, (int)min((size_t)lev, mipmap.size() - 1));
    }
    else {
        if (sp.lsm == LevelSampleMethod::L_LINEAR)
            return sample_bilinear(sp.p_uv, (int)min((size_t)lev, mipmap.size() - 1));
        else
            return sample_nearest(sp.p_uv, (int)min((size_t)lev, mipmap.size() - 1));
    }
}

float Texture::get_level(const SampleParams &sp) {
  // TODO: Task 6: Fill this in.
    if (sp.lsm == LevelSampleMethod::L_ZERO)
        return 0.0;

    else if (sp.lsm == LevelSampleMethod::L_NEAREST) {
        float r = floor(max(sqrt(pow(sp.p_dx_uv.x, 2.0) + pow(sp.p_dx_uv.y, 2.0)), sqrt(pow(sp.p_dy_uv.x, 2.0) + pow(sp.p_dy_uv.y, 2.0)))); 
        return (r == 0.0) ? r : log2(r);
    }
    else {
        float r = max(sqrt(pow(sp.p_dx_uv.x, 2.0) + pow(sp.p_dx_uv.y, 2.0)), sqrt(pow(sp.p_dy_uv.x, 2.0) + pow(sp.p_dy_uv.y, 2.0)));
        return log2(r);
    }
    //return 0.0;
}

Color MipLevel::get_texel(int tx, int ty) {
  return Color(&texels[tx * 3 + ty * width * 3]);
}

Color Texture::sample_nearest(Vector2D uv, int level) {
  // TODO: Task 5: Fill this in.
  // return magenta for invalid level
    if (level < 0 || level >= kMaxMipLevels)
        return Color(1, 0, 1);

    MipLevel* m = &mipmap[level];
    size_t x = max(min(m->width - 1, (size_t)floor(m->width * uv.x)), size_t(0)),
           y = max(min(m->height - 1, (size_t)floor(m->height * uv.y)), size_t(0));
    return m->get_texel(x, y);

}
// helper function to calculate a lerp
Color lerp(float x, Color v1, Color v2) {
    Color c = x * v2 + v1 * (1 - x);
    c.r = min(1.0f, max(0.0f, c.r));
    c.g = min(1.0f, max(0.0f, c.g));
    c.b = min(1.0f, max(0.0f, c.b));
    return c;
}

Color Texture::sample_bilinear(Vector2D uv, int level) {
  // TODO: Task 5: Fill this in.
  // return magenta for invalid level
    if (level < 0 || level >= kMaxMipLevels)
        return Color(1, 0, 1);

    MipLevel* m = &mipmap[level];
    float x = x = max(min((float) m->width - 1, (float) floor(m->width * uv.x)), 0.0f),
        y = max(min((float )m->height - 1, (float)floor(m->height * uv.y)), 0.0f);
    float t = m->height * uv.y - y, s = m->width * uv.x - x;

    Color u1 = lerp(s, m->get_texel((int)x, (int)y), m->get_texel((int)min(x + 1.0, m->width - 1.0), (int)y)),
        u2 = lerp(s, m->get_texel((int)min(x + 1.0, m->width - 1.0), (int)min(y + 1.0, m->height - 1.0)),
            m->get_texel((int)min(x + 1.0, m->width - 1.0), (int)min(y + 1.0, m->height - 1.0)));

    return lerp(t, u1, u2);
}

/****************************************************************************/

// Helpers



inline void uint8_to_float(float dst[3], unsigned char *src) {
  uint8_t *src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
}

inline void float_to_uint8(unsigned char *dst, float src[3]) {
  uint8_t *dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
  dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
  dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
}

void Texture::generate_mips(int startLevel) {

  // make sure there's a valid texture
  if (startLevel >= mipmap.size()) {
    std::cerr << "Invalid start level";
  }

  // allocate sublevels
  int baseWidth = mipmap[startLevel].width;
  int baseHeight = mipmap[startLevel].height;
  int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  mipmap.resize(startLevel + numSubLevels + 1);

  int width = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel &level = mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width = max(1, width / 2);
    // assert (width > 0);
    height = max(1, height / 2);
    // assert (height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(3 * width * height);
  }

  // create mips
  int subLevels = numSubLevels - (startLevel + 1);
  for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
       mipLevel++) {

    MipLevel &prevLevel = mipmap[mipLevel - 1];
    MipLevel &currLevel = mipmap[mipLevel];

    int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
    int currLevelPitch = currLevel.width * 3; // 32 bit RGB

    unsigned char *prevLevelMem;
    unsigned char *currLevelMem;

    currLevelMem = (unsigned char *)&currLevel.texels[0];
    prevLevelMem = (unsigned char *)&prevLevel.texels[0];

    float wDecimal, wNorm, wWeight[3];
    int wSupport;
    float hDecimal, hNorm, hWeight[3];
    int hSupport;

    float result[3];
    float input[3];

    // conditional differentiates no rounding case from round down case
    if (prevLevel.width & 1) {
      wSupport = 3;
      wDecimal = 1.0f / (float)currLevel.width;
    } else {
      wSupport = 2;
      wDecimal = 0.0f;
    }

    // conditional differentiates no rounding case from round down case
    if (prevLevel.height & 1) {
      hSupport = 3;
      hDecimal = 1.0f / (float)currLevel.height;
    } else {
      hSupport = 2;
      hDecimal = 0.0f;
    }

    wNorm = 1.0f / (2.0f + wDecimal);
    hNorm = 1.0f / (2.0f + hDecimal);

    // case 1: reduction only in horizontal size (vertical size is 1)
    if (currLevel.height == prevLevel.height) {
      // assert (currLevel.height == 1);

      for (int i = 0; i < currLevel.width; i++) {
        wWeight[0] = wNorm * (1.0f - wDecimal * i);
        wWeight[1] = wNorm * 1.0f;
        wWeight[2] = wNorm * wDecimal * (i + 1);

        result[0] = result[1] = result[2] = 0.0f;

        for (int ii = 0; ii < wSupport; ii++) {
          uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
          result[0] += wWeight[ii] * input[0];
          result[1] += wWeight[ii] * input[1];
          result[2] += wWeight[ii] * input[2];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (3 * i), result);
      }

      // case 2: reduction only in vertical size (horizontal size is 1)
    } else if (currLevel.width == prevLevel.width) {
      // assert (currLevel.width == 1);

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        result[0] = result[1] = result[2] = 0.0f;
        for (int jj = 0; jj < hSupport; jj++) {
          uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
          result[0] += hWeight[jj] * input[0];
          result[1] += hWeight[jj] * input[1];
          result[2] += hWeight[jj] * input[2];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (currLevelPitch * j), result);
      }

      // case 3: reduction in both horizontal and vertical size
    } else {

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          // convolve source image with a trapezoidal filter.
          // in the case of no rounding this is just a box filter of width 2.
          // in the general case, the support region is 3x3.
          for (int jj = 0; jj < hSupport; jj++)
            for (int ii = 0; ii < wSupport; ii++) {
              float weight = hWeight[jj] * wWeight[ii];
              uint8_to_float(input, prevLevelMem +
                                        prevLevelPitch * (2 * j + jj) +
                                        3 * (2 * i + ii));
              result[0] += weight * input[0];
              result[1] += weight * input[1];
              result[2] += weight * input[2];
            }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
        }
      }
    }
  }
}

} // namespace CGL
