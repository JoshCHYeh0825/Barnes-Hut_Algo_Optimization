#include "quadtree_tree.h"
#include <float.h>
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
unsigned int quad_find_quadrant(Quad *quad, Vec2 pos) {
    return ((pos.y > quad->center.y) << 1) | (pos.x > quad->center.x);
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
Node node_new(unsigned int next, Quad quad) {
    Node node;
    node.children = 0;
    node.next = next;
    node.pos = vec2_zero();
    node.mass = 0.0f;
    node.quad = quad;
    return node;
}

// Check if a node is a leaf (has no children)
int node_is_leaf(Node *node) {
    return node->children == 0;
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
    qt->nodes[ROOT] = node_new(0, quad);
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

    unsigned int children = qt->node_count;
    qt->nodes[node_index].children = children;

    unsigned int nexts[4] = {children + 1, children + 2, children + 3, qt->nodes[node_index].next};

    Quad quads[4];
    quad_subdivide(&qt->nodes[node_index].quad, quads);

    for (int i = 0; i < 4; i++) {
        qt->nodes[qt->node_count++] = node_new(nexts[i], quads[i]);
    }

    return children;
}

// Insert a position and mass into the quadtree
void quadtree_insert(Quadtree *qt, Vec2 pos, float mass) {
    unsigned int node = ROOT;

    while (node_is_branch(&qt->nodes[node])) {
        unsigned int quadrant = quad_find_quadrant(&qt->nodes[node].quad, pos);
        node = qt->nodes[node].children + quadrant;
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
        unsigned int children = quadtree_subdivide(qt, node);

        unsigned int q1 = quad_find_quadrant(&qt->nodes[node].quad, p);
        unsigned int q2 = quad_find_quadrant(&qt->nodes[node].quad, pos);

        if (q1 == q2) {
            node = children + q1;
        } else {
            unsigned int n1 = children + q1;
            unsigned int n2 = children + q2;

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
        unsigned int children = qt->nodes[node].children;

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
            Node *child_node = &qt->nodes[children + j];
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
// Optimized memory and accumulators
Vec2 quadtree_acc(Quadtree *qt, Vec2 pos) {
    // Vec2 acc = vec2_zero();
    float acc_x = 0.0f, acc_y = 0.0f;
    unsigned int node = ROOT;

    // Prefetch quadtree
    __builtin_prefetch(&qt->nodes[node], 0, 1);
    while (node < qt->node_count) {
        Node *n = &qt->nodes[node];

        // Skip nodes with no mass
        // Include node = n->next into the condition
        if (n->mass <= 0.0f) {
            if (n->next == 0)
                break;
            node = n->next;
            continue;
        }

        float dx = n->pos.x - pos.x;
        float dy = n->pos.y - pos.y;
        float d_sq = dx * dx + dy * dy;

        // Vec2 d = vec2_sub(n->pos, pos);
        // More efficient way to calculate magnitude squared
        // float d_sq = d.x * d.x + d.y * d.y;

        // Skip self-node (prevent self-gravity)
        // Include node = n->next into the condition
        if (d_sq < 0.0001f) {
            if (n->next == 0)
                break;
            node = n->next;
            continue;
        }

        if (node_is_leaf(n) || (n->quad.size * n->quad.size < d_sq * qt->t_sq)) {
            // Use approximation if far enough or a leaf
            float denom = powf(d_sq + qt->e_sq, 1.5f);
            if (denom > 0.0f) {
                /*
                float force = n->mass / denom;
                acc = vec2_add(acc, vec2_mul(d, force));
                */

                // Scalar float accumulators, replacing temporary vec2 structs
                // OpenMP and SIMD looped the calculations

                float f = n->mass / denom;
                acc_x += dx * f;
                acc_y += dy * f;
                // replaces acc to stop vec2 function calls
            }

            // Move to next node at same level
            if (n->next == 0)
                break;
            node = n->next;
        } else {
            // Descend into children if too close
            node = n->children;
        }
    }

    return /* acc */ vec2_new(acc_x, acc_y);
}

// Free the quadtree
void quadtree_free(Quadtree *qt) {
    free(qt->nodes);
    free(qt->parents);
    free(qt);
}