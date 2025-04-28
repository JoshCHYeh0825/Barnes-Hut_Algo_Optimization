#include "body.h"

Body body_new(Vec2 pos, Vec2 vel, float mass, float radius) {
    Body body = {.pos = pos, .vel = vel, .acc = vec2_zero(), .mass = mass, .radius = radius};
    return body;
}

void body_update(Body *body, float dt) {
    body->vel = vec2_add(body->vel, vec2_mul(body->acc, dt));
    body->pos = vec2_add(body->pos, vec2_mul(body->vel, dt));
}