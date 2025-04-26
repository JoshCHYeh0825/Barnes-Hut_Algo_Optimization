#include <xmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>
#include <cuda_runtime.h>

#ifndef VEC2_H
#define VEC2_H

typedef struct {
    float x;
    float y;
} Vec2;

__host__ __device__ Vec2 vec2_zero(void);
__host__ __device__ Vec2 vec2_new(float x, float y);
__host__ __device__ Vec2 vec2_add(Vec2 a, Vec2 b);
__host__ __device__ Vec2 vec2_sub(Vec2 a, Vec2 b);
__host__ __device__ Vec2 vec2_mul(Vec2 v, float s);
__host__ __device__ float vec2_mag(Vec2 v);
__host__ __device__ Vec2 vec2_normalize(Vec2 v);

#endif // VEC2_H 