#include "vec2.h"
#include <math.h>
#include <xmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>

Vec2 vec2_zero(void) {
    return (Vec2){0.0f, 0.0f};
}

Vec2 vec2_new(float x, float y) {
    return (Vec2){x, y};
}

Vec2 vec2_add(Vec2 a, Vec2 b) {
    __m128 va = _mm_set_ps(0, 0, a.y, a.x);  // Pack a into [0, 0, y, x]
    __m128 vb = _mm_set_ps(0, 0, b.y, b.x);  // Pack b into [0, 0, y, x]
    __m128 result = _mm_add_ps(va, vb);     // Add
    float res[4];
    _mm_storeu_ps(res, result);
    return (Vec2){res[0], res[1]};
}

Vec2 vec2_sub(Vec2 a, Vec2 b) {
    __m128 va = _mm_set_ps(0, 0, a.y, a.x);
    __m128 vb = _mm_set_ps(0, 0, b.y, b.x);
    __m128 result = _mm_sub_ps(va, vb);
    float res[4];
    _mm_storeu_ps(res, result);
    return (Vec2){res[0], res[1]};
}

Vec2 vec2_mul(Vec2 v, float s) {
    __m128 vv = _mm_set_ps(0, 0, v.y, v.x);
    __m128 ss = _mm_set1_ps(s);  // Broadcast scalar
    __m128 result = _mm_mul_ps(vv, ss);
    float res[4];
    _mm_storeu_ps(res, result);
    return (Vec2){res[0], res[1]};
}

float vec2_mag(Vec2 v) {
    return sqrtf(v.x * v.x + v.y * v.y);
}

Vec2 vec2_normalize(Vec2 v) {
    float mag = vec2_mag(v);
    if (mag == 0.0f) return vec2_zero();
    return vec2_mul(v, 1.0f / mag);
} 