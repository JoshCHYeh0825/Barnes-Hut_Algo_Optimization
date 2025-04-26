#include "vec2.h"
#include "body.h"
#include "quadtree.h"

#include <cuda_runtime.h>
#include <math.h>

__device__ Vec2 quadtree_acc(Quadtree* qt, Vec2 pos) {
    float acc_x = 0.0f, acc_y = 0.0f;

    unsigned int stack[16384]; // fixed stack depth
    int stack_size = 0;
    stack[stack_size++] = 0; // Root

    while (stack_size > 0) {
        unsigned int node_index = stack[--stack_size];
        Node* n = &qt->nodes[node_index];

        float dx = n->pos.x - pos.x;
        float dy = n->pos.y - pos.y;
        float d_sq = dx * dx + dy * dy;

        if (n->mass <= 0.0f || d_sq < 0.0001f) {
            continue;
        }

        if ((n->children[0] == UINT_MAX && n->children[1] == UINT_MAX &&
            n->children[2] == UINT_MAX && n->children[3] == UINT_MAX) || (n->size * n->size < d_sq * qt->t_sq)) {

            float denom = powf(d_sq + qt->e_sq, 1.5f);
            if (denom > 0.0f) {
                float f = n->mass / denom;
                acc_x += dx * f;
                acc_y += dy * f;
            }
        } else {
            // Add all non-null children
            for (int i = 0; i < 4; i++) {
                if (n->children[i] != UINT_MAX) {
                    stack[stack_size++] = n->children[i];
                }
            }
        }
    }

    return vec2_new(acc_x, acc_y);
}

__global__ void update_bodies_kernel(Body* bodies, Quadtree* quadtree, int num_bodies, float dt, float G) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_bodies) return;

    Vec2 acc = quadtree_acc(quadtree, bodies[idx].pos);
    bodies[idx].acc = vec2_mul(acc, G);

    bodies[idx].vel = vec2_add(bodies[idx].vel, vec2_mul(bodies[idx].acc, dt * 0.5f));
    bodies[idx].pos = vec2_add(bodies[idx].pos, vec2_mul(bodies[idx].vel, dt));
}

void update_simulation_gpu(Body* d_bodies, Quadtree* d_quadtree, int num_bodies, float dt, float G) {
    int threadsPerBlock = 256;
    int blocksPerGrid = (num_bodies + threadsPerBlock - 1) / threadsPerBlock;
    update_bodies_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_bodies, d_quadtree, num_bodies, dt, G);
    cudaDeviceSynchronize();
}