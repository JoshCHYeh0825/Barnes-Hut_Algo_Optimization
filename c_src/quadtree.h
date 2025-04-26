#ifndef QUADTREE_H
#define QUADTREE_H

#include "vec2.h"
#include "body.h"

// Structure to represent a quadrant
typedef struct {
    Vec2 center;
    float size;
} Quad;

// Structure to represent a node in the quadtree
// Flattening the tree, processor version first

/*
typedef struct {
    unsigned int children;  // Index of first child (0 if leaf)
    unsigned int next;      // Index of next node at same level
    Vec2 pos;               // Center of mass
    float mass;             // Total mass
    Quad quad;              // The quadrant this node represents
} Node;*/
typedef struct {
    unsigned int children[4]; // Instead of 1 child index, store 4
    Vec2 pos;                 // Center of mass
    float mass;               // Total mass
    float size;               // Size of quadrant
} Node;


// Structure to represent the quadtree
typedef struct {
    float t_sq;             // Theta squared (for approximation)
    float e_sq;             // Epsilon squared (for softening)
    Node* nodes;            // Array of nodes
    unsigned int* parents;  // Array of parent indices
    unsigned int node_count;    // Number of nodes
    unsigned int parent_count;  // Number of parent nodes
    unsigned int capacity;      // Capacity of nodes array
} Quadtree;

// Create a new quadrant containing all bodies
Quad quad_new_containing(Body* bodies, int count);

// Find which quadrant a position falls into (0-3)
unsigned int quad_find_quadrant(Vec2 center, Vec2 pos);

// Convert a quadrant into a subquadrant
Quad quad_into_quadrant(Quad quad, unsigned int quadrant);

// Subdivide a quadrant into four subquadrants
void quad_subdivide(Quad* quad, Quad* subquads);

// Create a new node
Node node_new(void);

// Check if a node is a leaf (has no children)
int node_is_leaf(Node* node);

// Check if a node is a branch (has children)
int node_is_branch(Node* node);

// Check if a node is empty (has no mass)
int node_is_empty(Node* node);

// Create a new quadtree
Quadtree* quadtree_new(float theta, float epsilon);

// Clear the quadtree and initialize with a new quadrant
void quadtree_clear(Quadtree* qt, Quad quad);

// Subdivide a node in the quadtree
unsigned int quadtree_subdivide(Quadtree* qt, unsigned int node_index);

// Insert a position and mass into the quadtree
void quadtree_insert(Quadtree* qt, Vec2 pos, float mass);

// Propagate center of mass calculations up the tree
void quadtree_propagate(Quadtree* qt);

// Calculate acceleration due to gravity at a position
Vec2 quadtree_acc(Quadtree* qt, Vec2 pos);

// Free the quadtree
void quadtree_free(Quadtree* qt);

#endif // QUADTREE_H 