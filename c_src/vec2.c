#include "vec2.h"
#include <math.h>

// pretty good, not much need for optimization
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