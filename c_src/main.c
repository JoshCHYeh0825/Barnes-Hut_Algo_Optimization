#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "body.h"
#include "quadtree.h"

#define WINDOW_WIDTH 900
#define WINDOW_HEIGHT 900
#define G 6.67430e-2f
#define THETA 0.5f
#define EPSILON 10.0f
#define NUM_BODIES 1000
#define NUM_TRIALS 5

#define MAX_VELOCITY 1.0f
#define TIME_SCALE 1.0f

SDL_Window* window = NULL;
SDL_Renderer* renderer = NULL;
bool running = true;
Body* bodies = NULL;
Quadtree* quadtree = NULL;

void initialize_simulation(void);
void update_simulation(float dt);
void cleanup_simulation(void);

void initialize_simulation(void) {
    srand((unsigned int)time(NULL));
    bodies = (Body*)malloc(NUM_BODIES * sizeof(Body));

    float center_x = WINDOW_WIDTH / 2.0f;
    float center_y = WINDOW_HEIGHT / 2.0f;
    float max_radius = fminf(WINDOW_WIDTH, WINDOW_HEIGHT) * 0.4f;

    for (int i = 0; i < NUM_BODIES; i++) {
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

// pretty straightforward already
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

// can optimize anything in this function
void update_simulation(float dt) {
    dt *= TIME_SCALE;

    Quad quad = quad_new_containing(bodies, NUM_BODIES);
    quadtree_clear(quadtree, quad);

    for (int i = 0; i < NUM_BODIES; i++)
        quadtree_insert(quadtree, bodies[i].pos, bodies[i].mass);
    quadtree_propagate(quadtree);

    // can optimize this - multithreading / gpu
    for (int i = 0; i < NUM_BODIES; i++) {
        bodies[i].acc = vec2_mul(quadtree_acc(quadtree, bodies[i].pos), G);
    }

    for (int i = 0; i < NUM_BODIES; i++) {
        bodies[i].vel = vec2_add(bodies[i].vel, vec2_mul(bodies[i].acc, dt * 0.5f));
        bodies[i].pos = vec2_add(bodies[i].pos, vec2_mul(bodies[i].vel, dt));
        handle_wall_collisions(&bodies[i]);
    }

    quad = quad_new_containing(bodies, NUM_BODIES);
    quadtree_clear(quadtree, quad);
    for (int i = 0; i < NUM_BODIES; i++)
        quadtree_insert(quadtree, bodies[i].pos, bodies[i].mass);
    quadtree_propagate(quadtree);

    // can optimize this - multithreading / gpu
    for (int i = 0; i < NUM_BODIES; i++) {
        Vec2 new_acc = vec2_mul(quadtree_acc(quadtree, bodies[i].pos), G);
        bodies[i].vel = vec2_add(bodies[i].vel, vec2_mul(new_acc, dt * 0.5f));
        bodies[i].acc = new_acc;
    }
}

void cleanup_simulation(void) {
    free(bodies);
    quadtree_free(quadtree);
}

void render(void) {
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    for (int i = 0; i < NUM_BODIES; i++) {
        SDL_FRect rect = {
            bodies[i].pos.x - bodies[i].radius,
            bodies[i].pos.y - bodies[i].radius,
            bodies[i].radius * 2.0f,
            bodies[i].radius * 2.0f
        };
        SDL_RenderFillRectF(renderer, &rect);
    }
    SDL_RenderPresent(renderer);
}

int main(int argc, char* argv[]) {
        
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL_Init Error: %s\n", SDL_GetError());
        return 1;
    }

    window = SDL_CreateWindow("Barnes-Hut Simulation",
                            SDL_WINDOWPOS_CENTERED,
                            SDL_WINDOWPOS_CENTERED,
                            WINDOW_WIDTH, WINDOW_HEIGHT, 0);
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    initialize_simulation();

    Uint32 last_time = SDL_GetTicks();

    // optimize here!!!!
    Uint32 start = SDL_GetTicks();
    int iterations = 0;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT)
                running = false;
        }

        Uint32 current_time = SDL_GetTicks();
        float dt = (current_time - last_time) / 1000.0f;
        last_time = current_time;

        if (dt > 0.05f) dt = 0.05f;

        update_simulation(dt);
        render();
        SDL_Delay(16);

        iterations++;

        // if its been 10 seconds, exit
        if (current_time - start > 10000) {
            running = false;
        }
    }

    cleanup_simulation();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    // print the average iterations per second
    printf("Average iterations per second: %f\n", (float)iterations / 10.0f);
    
    return 0;
}

