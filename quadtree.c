#include "quadtree.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>

#define ROOT 0
#define QUADTREE_INITIAL_CAPACITY 256

// Create a new quadrant containing all bodies
Quad quad_new_containing(Body *bodies, int count) {
    float min_x = FLT_MAX;
    float min_y = FLT_MAX;
    float max_x = -FLT_MAX;
    float max_y = -FLT_MAX;

    for (int i = 0; i < count; i++) {
        min_x = fminf(min_x, bodies[i].pos.x);
        min_y = fminf(min_y, bodies[i].pos.y);
        max_x = fmaxf(max_x, bodies[i].pos.x);
        max_y = fmaxf(max_y, bodies[i].pos.y);
    }

    Vec2 center = vec2_mul(vec2_new(min_x + max_x, min_y + max_y), 0.5f);
    float size = fmaxf(max_x - min_x, max_y - min_y);

    Quad quad = {center, size};
    return quad;
}

// Find which quadrant a position falls into (0-3)
unsigned int quad_find_quadrant(/* Quad* quad */ Vec2 center, Vec2 pos) {
    return ((pos.y > center.y) << 1) | (pos.x > center.x);
}

// Convert a quadrant into a subquadrant
Quad quad_into_quadrant(Quad quad, unsigned int quadrant) {
    quad.size *= 0.5f;
    quad.center.x += ((quadrant & 1) - 0.5f) * quad.size;
    quad.center.y += ((quadrant >> 1) - 0.5f) * quad.size;
    return quad;
}

// Subdivide a quadrant into four subquadrants
void quad_subdivide(Quad *quad, Quad *subquads) {
    for (int i = 0; i < 4; i++) {
        subquads[i] = quad_into_quadrant(*quad, i);
    }
}

// Create a new node
/* Node node_new(unsigned int next, Quad quad) {
    Node node;
    node.children = 0;
    node.next = next;
    node.pos = vec2_zero();
    node.mass = 0.0f;
    node.quad = quad;
    return node;
}
    */
Node node_new(void) {
    Node node;
    for (int i = 0; i < 4; i++)
        node.children[i] = UINT_MAX;
    node.pos = vec2_zero();
    node.mass = 0.0f;
    node.size = 0.0f;
    return node;
}

// Check if a node is a leaf (has no children)
int node_is_leaf(Node *node) {
    return (node->children[0] == UINT_MAX && node->children[1] == UINT_MAX && node->children[2] == UINT_MAX &&
            node->children[3] == UINT_MAX);
}

// Check if a node is a branch (has children)
int node_is_branch(Node *node) {
    return node->children != 0;
}

// Check if a node is empty (has no mass)
int node_is_empty(Node *node) {
    return node->mass == 0.0f;
}

// Create a new quadtree
Quadtree *quadtree_new(float theta, float epsilon) {
    Quadtree *qt = (Quadtree *)malloc(sizeof(Quadtree));
    qt->t_sq = theta * theta;
    qt->e_sq = epsilon * epsilon;
    qt->capacity = QUADTREE_INITIAL_CAPACITY;
    qt->nodes = (Node *)malloc(qt->capacity * sizeof(Node));
    qt->parents = (unsigned int *)malloc(qt->capacity * sizeof(unsigned int));
    qt->node_count = 0;
    qt->parent_count = 0;
    return qt;
}

// Clear the quadtree and initialize with a new quadrant
void quadtree_clear(Quadtree *qt, Quad quad) {
    qt->node_count = 1;
    qt->parent_count = 0;
    // qt->nodes[ROOT] = node_new(0, quad);
    qt->nodes[ROOT] = node_new();
    qt->nodes[ROOT].pos = quad.center;
    qt->nodes[ROOT].size = quad.size;
}

// Ensure the quadtree has enough capacity
void quadtree_ensure_capacity(Quadtree *qt, unsigned int needed) {
    if (qt->capacity >= needed)
        return;

    while (qt->capacity < needed) {
        qt->capacity *= 2;
    }

    qt->nodes = (Node *)realloc(qt->nodes, qt->capacity * sizeof(Node));
    qt->parents = (unsigned int *)realloc(qt->parents, qt->capacity * sizeof(unsigned int));
}

// Subdivide a node in the quadtree
unsigned int quadtree_subdivide(Quadtree *qt, unsigned int node_index) {
    // Ensure capacity for 4 new nodes
    quadtree_ensure_capacity(qt, qt->node_count + 4);

    // Record this as a parent
    qt->parents[qt->parent_count++] = node_index;

    // unsigned int children = qt->node_count;
    // qt->nodes[node_index].children = children;
    // qt->parents[qt->parent_count++] = node_index;
    Node *parent = &qt->nodes[node_index];
    float half_size = parent->size * 0.5f;

    /*
    unsigned int nexts[4] = {
        children + 1,
        children + 2,
        children + 3,
        qt->nodes[node_index].next
    };

    Quad quads[4];
    quad_subdivide(&qt->nodes[node_index].quad, quads);

    for (int i = 0; i < 4; i++) {
        qt->nodes[qt->node_count++] = node_new(nexts[i], quads[i]);
    }
    */

    for (int i = 0; i < 4; i++) {
        Node child = node_new();
        child.size = half_size;
        // Position center adjustment based on quadrant
        child.pos.x = parent->pos.x + ((i & 1) ? 0.5f : -0.5f) * half_size;
        child.pos.y = parent->pos.y + ((i >> 1) ? 0.5f : -0.5f) * half_size;

        parent->children[i] = qt->node_count;
        qt->nodes[qt->node_count++] = child;
    }

    return node_index;
}

// Insert a position and mass into the quadtree
void quadtree_insert(Quadtree *qt, Vec2 pos, float mass) {
    unsigned int node = ROOT;

    while (/* node_is_branch */ !node_is_leaf(&qt->nodes[node])) {
        unsigned int quadrant = quad_find_quadrant(qt->nodes[node].pos, pos);
        node = qt->nodes[node].children[quadrant];
    }

    if (node_is_empty(&qt->nodes[node])) {
        qt->nodes[node].pos = pos;
        qt->nodes[node].mass = mass;
        return;
    }

    Vec2 p = qt->nodes[node].pos;
    float m = qt->nodes[node].mass;

    // If positions are identical, just add mass
    if (p.x == pos.x && p.y == pos.y) {
        qt->nodes[node].mass += mass;
        return;
    }

    // Otherwise, subdivide until we can separate them
    while (1) {
        if (qt->nodes[node].size < 1e-6f) { // Safeguard against infinite subdivisions
            qt->nodes[node].mass += mass;
            return;
        }
        quadtree_subdivide(qt, node);

        unsigned int q1 = quad_find_quadrant(qt->nodes[node].pos, p);
        unsigned int q2 = quad_find_quadrant(qt->nodes[node].pos, pos);

        if (q1 == q2) {
            node = qt->nodes[node].children[q1];
        } else {
            unsigned int n1 = qt->nodes[node].children[q1];
            unsigned int n2 = qt->nodes[node].children[q2];

            qt->nodes[n1].pos = p;
            qt->nodes[n1].mass = m;
            qt->nodes[n2].pos = pos;
            qt->nodes[n2].mass = mass;

            return;
        }
    }
}

// Propagate center of mass calculations up the tree
// Changed to loop unrolling
// Trying OpenMP
void quadtree_propagate(Quadtree *qt) {

#pragma omp parallel for
    for (int i = qt->parent_count - 1; i >= 0; i--) {
        unsigned int node = qt->parents[i];
        Node *parent = &qt->nodes[node];

        // Declared these new floats
        float total_mass = 0.0f;
        float com_x = 0.0f;
        float com_y = 0.0f;

        /*
        Vec2 com = vec2_zero();
        float total_mass = 0.0f;

        for (int j = 0; j < 4; j++) {
            unsigned int child = children + j;
            float child_mass = qt->nodes[child].mass;

            com = vec2_add(com, vec2_mul(qt->nodes[child].pos, child_mass));
            total_mass += child_mass;
        }
        */

        // Manual unrolling for 4 children
        for (int j = 0; j < 4; j++) {
            unsigned int child_idx = parent->children[j];
            if (child_idx == UINT_MAX)
                continue;
            Node *child_node = &qt->nodes[child_idx];

            float mass = child_node->mass;
            total_mass += mass;
            com_x += child_node->pos.x * mass;
            com_y += child_node->pos.y * mass;
        }

        qt->nodes[node].mass = total_mass;
        if (total_mass > 0.0f) {
            qt->nodes[node].pos.x = com_x / total_mass;
            qt->nodes[node].pos.y = com_y / total_mass;
        }
    }
}

// Calculate acceleration due to gravity at a position
// Flattened flat tree
Vec2 quadtree_acc(Quadtree *qt, Vec2 pos) {
    float acc_x = 0.0f, acc_y = 0.0f;

    unsigned int stack[16384]; // fixed stack depth
    int stack_size = 0;
    stack[stack_size++] = ROOT;

    while (stack_size > 0) {
        unsigned int node_index = stack[--stack_size];
        Node *n = &qt->nodes[node_index];

// Prefetch all children
#pragma omp simd
        for (int j = 0; j < 4; j++) {
            if (n->children[j] != UINT_MAX) {
                __builtin_prefetch(&qt->nodes[n->children[j]], 0, 1);
            }
        }

        float dx = n->pos.x - pos.x;
        float dy = n->pos.y - pos.y;
        float d_sq = dx * dx + dy * dy;

        if (n->mass <= 0.0f || d_sq < 0.0001f) {
            continue;
        }

        if (node_is_leaf(n) || (n->size * n->size < d_sq * qt->t_sq)) {
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

// Free the quadtree
void quadtree_free(Quadtree *qt) {
    free(qt->nodes);
    free(qt->parents);
    free(qt);
}