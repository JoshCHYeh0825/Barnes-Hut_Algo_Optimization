#include "body.h"
#include "quadtree.h"
#include <SDL2/SDL.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define WINDOW_WIDTH 450
#define WINDOW_HEIGHT 450
#define G 6.67430e-2f
#define THETA 0.5f
#define EPSILON 10.0f
#define NUM_TRIALS 10

#define MAX_VELOCITY 1.0f
#define TIME_SCALE 1.0f

SDL_Window *window = NULL;
SDL_Renderer *renderer = NULL;
Body *bodies = NULL;
Quadtree *quadtree = NULL;

void initialize_simulation(int num_bodies);
void update_simulation(float dt, int num_bodies);
void cleanup_simulation(void);

void initialize_simulation(int num_bodies) {
    srand((unsigned int)time(NULL));
    bodies = (Body *)malloc(num_bodies * sizeof(Body));

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

// pretty straightforward already
void handle_wall_collisions(Body *body) {
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
void update_simulation(float dt, int num_bodies) {
    dt *= TIME_SCALE;

    Quad quad = quad_new_containing(bodies, num_bodies);
    quadtree_clear(quadtree, quad);

    for (int i = 0; i < num_bodies; i++)
        quadtree_insert(quadtree, bodies[i].pos, bodies[i].mass);
    quadtree_propagate(quadtree);

    // can optimize this - multithreading / gpu
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

    // can optimize this - multithreading / gpu
    for (int i = 0; i < num_bodies; i++) {
        Vec2 new_acc = vec2_mul(quadtree_acc(quadtree, bodies[i].pos), G);
        bodies[i].vel = vec2_add(bodies[i].vel, vec2_mul(new_acc, dt * 0.5f));
        bodies[i].acc = new_acc;
    }
}

void cleanup_simulation(void) {
    free(bodies);
    quadtree_free(quadtree);
}

void render(int num_bodies) {
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    for (int i = 0; i < num_bodies; i++) {
        SDL_FRect rect = {bodies[i].pos.x - bodies[i].radius, bodies[i].pos.y - bodies[i].radius,
                          bodies[i].radius * 2.0f, bodies[i].radius * 2.0f};
        SDL_RenderFillRectF(renderer, &rect);
    }
    SDL_RenderPresent(renderer);
}

int main(int argc, char *argv[]) {

    int a = 1000;
    int b = 500;
    int c = 100;

    int num_bodies = 0;
    int i = 0;

    // logging
    int num_bodies_log[NUM_TRIALS];
    int num_iterations_log[NUM_TRIALS];

    for (i = 0; i < NUM_TRIALS; i++) {
        bodies = NULL;
        quadtree = NULL;

        num_bodies = a * i * i + b * i + c;
        printf("num_bodies: %d\n", num_bodies);

        if (SDL_Init(SDL_INIT_VIDEO) != 0) {
            fprintf(stderr, "SDL_Init Error: %s\n", SDL_GetError());
            return 1;
        }

        window = SDL_CreateWindow("Barnes-Hut Simulation", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WINDOW_WIDTH,
                                  WINDOW_HEIGHT, 0);
        renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

        initialize_simulation(num_bodies);

        Uint32 last_time = SDL_GetTicks();

        // optimize here!!!!
        Uint32 start = SDL_GetTicks();
        int iterations = 0;
        bool running = true;
        while (running) {
            SDL_Event e;
            while (SDL_PollEvent(&e)) {
                if (e.type == SDL_QUIT)
                    running = false;
            }

            Uint32 current_time = SDL_GetTicks();
            float dt = (current_time - last_time) / 1000.0f;
            last_time = current_time;

            if (dt > 0.05f)
                dt = 0.05f;

            update_simulation(dt, num_bodies);
            render(num_bodies);
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

        // update logging
        num_bodies_log[i] = num_bodies;
        num_iterations_log[i] = iterations;
    }

    printf("num_bodies,num_iterations\n");
    for (i = 0; i < NUM_TRIALS; i++) {
        printf("%d %.3f\n", num_bodies_log[i], (float)num_iterations_log[i] / 10.0);
    }

    return 0;
}
