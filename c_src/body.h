#ifndef BODY_H
#define BODY_H

#include "vec2.h"

typedef struct {
    Vec2 pos;
    Vec2 vel;
    Vec2 acc;
    float mass;
    float radius;
} Body;

Body body_new(Vec2 pos, Vec2 vel, float mass, float radius);
void body_update(Body* body, float dt);

#endif // BODY_H 