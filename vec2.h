#ifndef VEC2_H
#define VEC2_H

#ifdef __CUDACC__
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif

typedef struct {
    float x, y;
} Vec2;

// Declare functions
HOST_DEVICE Vec2 vec2_zero(void);
HOST_DEVICE Vec2 vec2_new(float x, float y);
HOST_DEVICE Vec2 vec2_add(Vec2 a, Vec2 b);
HOST_DEVICE Vec2 vec2_sub(Vec2 a, Vec2 b);
HOST_DEVICE Vec2 vec2_mul(Vec2 v, float s);
HOST_DEVICE float vec2_mag(Vec2 v);
HOST_DEVICE Vec2 vec2_normalize(Vec2 v);

#endif // VEC2_H
