#include "vec2.h"
#include "body.h"
#include "quadtree.h"

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

__device__ Vec2 quadtree_acc(Quadtree* qt, Vec2 pos) {
    Vec2 acc = vec2_zero();
    if (qt->node_count == 0 || qt->nodes[ROOT].mass == 0.0f) return acc;

    unsigned int node = ROOT;
    while (node < qt->node_count) {
        Node* n = &qt->nodes[node];

        if (n->mass <= 0.0f) {
            node = (n->next == 0) ? qt->node_count : n->next;
            continue;
        }

        Vec2 d = vec2_sub(n->pos, pos);
        float d_sq = d.x * d.x + d.y * d.y;
        if (d_sq < 0.0001f) {
            node = (n->next == 0) ? qt->node_count : n->next;
            continue;
        }

        float s2 = n->quad.size * n->quad.size;
        if (node_is_leaf(n) || s2 < d_sq * qt->t_sq) {
            float denom = powf(d_sq + qt->e_sq, 1.5f);
            if (denom > 0.0f)
                acc = vec2_add(acc, vec2_mul(d, n->mass / denom));
            node = (n->next == 0) ? qt->node_count : n->next;
        } else {
            node = n->children;
        }
    }

    return acc;
}

__global__ void update_bodies_kernel(Body* bodies, Quadtree* quadtree, int num_bodies, float dt, float G) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_bodies) return;

    Vec2 acc = quadtree_acc(quadtree, bodies[idx].pos);
    bodies[idx].acc = vec2_mul(acc, G);

    bodies[idx].vel = vec2_add(bodies[idx].vel, vec2_mul(bodies[idx].acc, dt * 0.5f));
    bodies[idx].pos = vec2_add(bodies[idx].pos, vec2_mul(bodies[idx].vel, dt));
}

void initialize_simulation(int num_bodies);
void update_simulation(float dt, int num_bodies);

void initialize_simulation(int num_bodies) {
    srand((unsigned int)time(NULL));
    bodies = (Body*)malloc(num_bodies * sizeof(Body));

    float center_x = WINDOW_WIDTH / 2.0f;
    float center_y = WINDOW_HEIGHT / 2.0f;
    float max_radius = fminf(WINDOW_WIDTH, WINDOW_HEIGHT) * 0.4f;

    for (int i = 0; i < num_bodies; i++) {
        float angle = ((float)rand() / RAND_MAX) * 2.0f * M_PI;
        float distance = ((float)rand() / RAND_MAX) * max_radius;
        float x = center_x + cosf(angle) * distance;
        float y = center_y + sinf(angle) * distance;

        float vel_angle = angle + M_PI / 2.0f + (((float)rand() / RAND_MAX) - 0.5f) * 0.5f;
        float vel_magnitude = MAX_VELOCITY * sqrtf(distance / max_radius);
        float vx = cosf(vel_angle) * vel_magnitude;
        float vy = sinf(vel_angle) * vel_magnitude;

        float radius = 2.0f + ((float)rand() / RAND_MAX) * 6.0f;
        float mass = radius * radius;

        bodies[i] = body_new(vec2_new(x, y), vec2_new(vx, vy), mass, radius);
    }

    quadtree = quadtree_new(THETA, EPSILON);
}

void update_simulation(float dt, int num_bodies) {
    dt *= TIME_SCALE;

    Quad quad = quad_new_containing(bodies, num_bodies);
    quadtree_clear(quadtree, quad);

    for (int i = 0; i < num_bodies; i++)
        quadtree_insert(quadtree, bodies[i].pos, bodies[i].mass);
    quadtree_propagate(quadtree);

    for (int i = 0; i < num_bodies; i++) {
        bodies[i].acc = vec2_mul(quadtree_acc(quadtree, bodies[i].pos), G);
    }

    for (int i = 0; i < num_bodies; i++) {
        bodies[i].vel = vec2_add(bodies[i].vel, vec2_mul(bodies[i].acc, dt * 0.5f));
        bodies[i].pos = vec2_add(bodies[i].pos, vec2_mul(bodies[i].vel, dt));
        handle_wall_collisions(&bodies[i]);
    }

    quad = quad_new_containing(bodies, num_bodies);
    quadtree_clear(quadtree, quad);
    for (int i = 0; i < num_bodies; i++)
        quadtree_insert(quadtree, bodies[i].pos, bodies[i].mass);
    quadtree_propagate(quadtree);

    for (int i = 0; i < num_bodies; i++) {
        Vec2 new_acc = vec2_mul(quadtree_acc(quadtree, bodies[i].pos), G);
        bodies[i].vel = vec2_add(bodies[i].vel, vec2_mul(new_acc, dt * 0.5f));
        bodies[i].acc = new_acc;
    }
}

void handle_wall_collisions(Body* body) {
    float damping = 0.8f;

    if (body->pos.x - body->radius < 0) {
        body->pos.x = body->radius;
        body->vel.x = fabsf(body->vel.x) * damping;
    } else if (body->pos.x + body->radius > WINDOW_WIDTH) {
        body->pos.x = WINDOW_WIDTH - body->radius;
        body->vel.x = -fabsf(body->vel.x) * damping;
    }

    if (body->pos.y - body->radius < 0) {
        body->pos.y = body->radius;
        body->vel.y = fabsf(body->vel.y) * damping;
    } else if (body->pos.y + body->radius > WINDOW_HEIGHT) {
        body->pos.y = WINDOW_HEIGHT - body->radius;
        body->vel.y = -fabsf(body->vel.y) * damping;
    }
}

// CPU-callable function
void update_simulation_gpu(Body* d_bodies, Quadtree* d_quadtree, int num_bodies, float dt, float G) {
    int threadsPerBlock = 256;
    int blocksPerGrid = (num_bodies + threadsPerBlock - 1) / threadsPerBlock;
    update_bodies_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_bodies, d_quadtree, num_bodies, dt, G);
    cudaDeviceSynchronize();
}