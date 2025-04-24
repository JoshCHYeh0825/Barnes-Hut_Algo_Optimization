#include "vec2.h"
#include <immintrin.h>  // AVX/SSE intrinsics
#include <math.h>

// Scalar fallback for single ops
Vec2 vec2_zero(void) {
    return (Vec2){0.0f, 0.0f};
}

Vec2 vec2_new(float x, float y) {
    return (Vec2){x, y};
}

Vec2 vec2_add(Vec2 a, Vec2 b) {
    return (Vec2){a.x + b.x, a.y + b.y};
}

Vec2 vec2_sub(Vec2 a, Vec2 b) {
    return (Vec2){a.x - b.x, a.y - b.y};
}

Vec2 vec2_mul(Vec2 v, float s) {
    return (Vec2){v.x * s, v.y * s};
}

float vec2_mag(Vec2 v) {
    return sqrtf(v.x * v.x + v.y * v.y);
}

Vec2 vec2_normalize(Vec2 v) {
    float mag = vec2_mag(v);
    if (mag == 0.0f) return vec2_zero();
    return vec2_mul(v, 1.0f / mag);
}

// SIMD batch operations: expects arrays of Vec2s
void vec2_array_add(Vec2* out, const Vec2* a, const Vec2* b, int count) {
    for (int i = 0; i < count; i += 2) {
        __m128 va = _mm_loadu_ps((float*)&a[i]);  // Load 2 Vec2s (4 floats)
        __m128 vb = _mm_loadu_ps((float*)&b[i]);
        __m128 res = _mm_add_ps(va, vb);
        _mm_storeu_ps((float*)&out[i], res);
    }
}

void vec2_array_mul_scalar(Vec2* out, const Vec2* a, float scalar, int count) {
    __m128 ss = _mm_set1_ps(scalar);
    for (int i = 0; i < count; i += 2) {
        __m128 va = _mm_loadu_ps((float*)&a[i]);
        __m128 res = _mm_mul_ps(va, ss);
        _mm_storeu_ps((float*)&out[i], res);
    }
}
