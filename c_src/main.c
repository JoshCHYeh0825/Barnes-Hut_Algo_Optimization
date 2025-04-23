#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "body.h"
#include "quadtree.h"

#define WINDOW_WIDTH 900
#define WINDOW_HEIGHT 900
#define G 6.67430e-2f       // Gravitational constant (scaled for visualization)
#define THETA 0.5f          // Barnes-Hut theta parameter (accuracy vs. speed)
#define EPSILON 10.0f       // Softening parameter to prevent singularities
#define NUM_BODIES 50       // Number of bodies in simulation
#define MAX_VELOCITY 1.0f   // Maximum initial velocity
#define TIME_SCALE 1.0f     // Time scale factor (1.0 = normal speed)

HWND hwnd;
HDC hdc;
bool running = true;
Body* bodies = NULL;
Quadtree* quadtree = NULL;

// Forward declarations
LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);
void initialize_simulation(void);
void update_simulation(float dt);
void cleanup_simulation(void);

// Initialize the simulation
void initialize_simulation(void) {
    srand((unsigned int)time(NULL));
    
    // Allocate bodies
    bodies = (Body*)malloc(NUM_BODIES * sizeof(Body));
    
    // Create bodies in a disc with random velocities
    float center_x = WINDOW_WIDTH / 2.0f;
    float center_y = WINDOW_HEIGHT / 2.0f;
    float max_radius = fminf(WINDOW_WIDTH, WINDOW_HEIGHT) * 0.4f;
    
    for (int i = 0; i < NUM_BODIES; i++) {
        // Random position in a disc
        float angle = ((float)rand() / RAND_MAX) * 2.0f * 3.14159f;
        float distance = ((float)rand() / RAND_MAX) * max_radius;
        float x = center_x + cosf(angle) * distance;
        float y = center_y + sinf(angle) * distance;
        
        // Random velocity, perpendicular to radius with some randomness
        float vel_angle = angle + 3.14159f / 2.0f + (((float)rand() / RAND_MAX) - 0.5f) * 0.5f;
        float vel_magnitude = MAX_VELOCITY * sqrtf(distance / max_radius);
        float vx = cosf(vel_angle) * vel_magnitude;
        float vy = sinf(vel_angle) * vel_magnitude;
        
        // Mass proportional to radius
        float radius = 2.0f + ((float)rand() / RAND_MAX) * 6.0f;
        float mass = radius * radius;
        
        bodies[i] = body_new(vec2_new(x, y), vec2_new(vx, vy), mass, radius);
    }
    
    // Create quadtree
    quadtree = quadtree_new(THETA, EPSILON);
}

// Handle wall collisions with improved stability
void handle_wall_collisions(Body* body) {
    float damping = 0.8f;  // Energy loss factor (0.8 = 20% energy loss on collision)
    
    // Check if the body is crossing any wall boundary
    bool collision_x = false;
    bool collision_y = false;
    
    // X-boundary collisions
    if (body->pos.x - body->radius < 0) {
        body->pos.x = body->radius;
        body->vel.x = fabsf(body->vel.x) * damping;
        collision_x = true;
    } else if (body->pos.x + body->radius > WINDOW_WIDTH) {
        body->pos.x = WINDOW_WIDTH - body->radius;
        body->vel.x = -fabsf(body->vel.x) * damping;
        collision_x = true;
    }
    
    // Y-boundary collisions
    if (body->pos.y - body->radius < 0) {
        body->pos.y = body->radius;
        body->vel.y = fabsf(body->vel.y) * damping;
        collision_y = true;
    } else if (body->pos.y + body->radius > WINDOW_HEIGHT) {
        body->pos.y = WINDOW_HEIGHT - body->radius;
        body->vel.y = -fabsf(body->vel.y) * damping;
        collision_y = true;
    }
    
    // If there was a collision, slightly reduce the perpendicular velocity component
    // to simulate friction with the wall
    if (collision_x) {
        body->vel.y *= 0.95f;
    }
    if (collision_y) {
        body->vel.x *= 0.95f;
    }
}

// Update the simulation
void update_simulation(float dt) {
    // Scale dt by the time scale factor
    dt *= TIME_SCALE;
    
    // Create a quad containing all bodies
    Quad quad = quad_new_containing(bodies, NUM_BODIES);
    quadtree_clear(quadtree, quad);
    
    // Insert all bodies into the quadtree
    for (int i = 0; i < NUM_BODIES; i++) {
        quadtree_insert(quadtree, bodies[i].pos, bodies[i].mass);
    }
    
    // Propagate center of mass calculations
    quadtree_propagate(quadtree);
    
    // Calculate accelerations, PARALLELIZED
    #pragma omp parallel for
    for (int i = 0; i < NUM_BODIES; i++) {
        bodies[i].acc = quadtree_acc(quadtree, bodies[i].pos);
        // Scale by G
        bodies[i].acc = vec2_mul(bodies[i].acc, G);
    }
    
    // Update positions and velocities using a more stable integrator (velocity Verlet)
    #pragma omp parallel for
    for (int i = 0; i < NUM_BODIES; i++) {
        // Save old acceleration for integrator
        Vec2 old_acc = bodies[i].acc;
        
        // First half of velocity update
        bodies[i].vel = vec2_add(bodies[i].vel, vec2_mul(old_acc, dt * 0.5f));
        
        // Position update
        bodies[i].pos = vec2_add(bodies[i].pos, vec2_mul(bodies[i].vel, dt));
        
        // Handle wall collisions
        handle_wall_collisions(&bodies[i]);
    }
    
    // Rebuild quadtree with new positions
    quad = quad_new_containing(bodies, NUM_BODIES);
    quadtree_clear(quadtree, quad);
    
    for (int i = 0; i < NUM_BODIES; i++) {
        quadtree_insert(quadtree, bodies[i].pos, bodies[i].mass);
    }
    
    quadtree_propagate(quadtree);
    
    // Second half of velocity update
    for (int i = 0; i < NUM_BODIES; i++) {
        Vec2 new_acc = quadtree_acc(quadtree, bodies[i].pos);
        new_acc = vec2_mul(new_acc, G);
        
        bodies[i].vel = vec2_add(bodies[i].vel, vec2_mul(new_acc, dt * 0.5f));
        bodies[i].acc = new_acc;
    }
}

// Cleanup the simulation
void cleanup_simulation(void) {
    free(bodies);
    quadtree_free(quadtree);
}

// Initialize the Win32 window
bool init_window(HINSTANCE hInstance) {
    // Register the window class
    WNDCLASS wc = {0};
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = "BarnesHutClass";
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    
    RegisterClass(&wc);
    
    // Create the window
    hwnd = CreateWindowEx(
        0,                              // Optional window styles
        "BarnesHutClass",               // Window class
        "Barnes-Hut Simulation",        // Window text
        WS_OVERLAPPEDWINDOW,            // Window style
        CW_USEDEFAULT, CW_USEDEFAULT,   // Position
        WINDOW_WIDTH, WINDOW_HEIGHT,    // Size
        NULL,                           // Parent window    
        NULL,                           // Menu
        hInstance,                      // Instance handle
        NULL                            // Additional application data
    );
    
    if (hwnd == NULL) {
        MessageBox(NULL, "Window Creation Failed!", "Error", MB_ICONEXCLAMATION | MB_OK);
        return false;
    }
    
    // Get device context for drawing
    hdc = GetDC(hwnd);
    if (hdc == NULL) {
        MessageBox(NULL, "Failed to get device context!", "Error", MB_ICONEXCLAMATION | MB_OK);
        return false;
    }
    
    // Show window
    ShowWindow(hwnd, SW_SHOWDEFAULT);
    UpdateWindow(hwnd);
    
    return true;
}

// Window procedure
LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
    switch (uMsg) {
        case WM_CLOSE:
            running = false;
            PostQuitMessage(0);
            return 0;
        case WM_DESTROY:
            running = false;
            PostQuitMessage(0);
            return 0;
    }
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

// Process Windows messages
void process_messages() {
    MSG msg;
    while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
        if (msg.message == WM_QUIT) {
            running = false;
        }
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
}

// Render the bodies
void render(Body* bodies, int count) {
    // Clear the window
    RECT clientRect;
    GetClientRect(hwnd, &clientRect);
    HBRUSH blackBrush = CreateSolidBrush(RGB(0, 0, 0));
    FillRect(hdc, &clientRect, blackBrush);
    DeleteObject(blackBrush);
    
    // Draw each body as a white circle
    HBRUSH whiteBrush = CreateSolidBrush(RGB(255, 255, 255));
    HPEN nullPen = (HPEN)GetStockObject(NULL_PEN);
    
    SelectObject(hdc, nullPen);
    SelectObject(hdc, whiteBrush);
    
    for (int i = 0; i < count; i++) {
        int x = (int)(bodies[i].pos.x - bodies[i].radius);
        int y = (int)(bodies[i].pos.y - bodies[i].radius);
        int width = (int)(bodies[i].radius * 2);
        int height = (int)(bodies[i].radius * 2);
        
        Ellipse(hdc, x, y, x + width, y + height);
    }
    
    DeleteObject(whiteBrush);
}

// Main entry point
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow) {
    (void)hPrevInstance;
    (void)lpCmdLine;
    (void)nCmdShow;
    
    if (!init_window(hInstance)) {
        return 1;
    }
    
    initialize_simulation();
    
    // Main loop
    DWORD last_time = GetTickCount();
    while (running) {
        process_messages();
        
        DWORD current_time = GetTickCount();
        float dt = (current_time - last_time) / 1000.0f;
        last_time = current_time;
        
        // Limit dt to avoid instability
        if (dt > 0.05f) dt = 0.05f;
        
        update_simulation(dt);
        render(bodies, NUM_BODIES);
        
        Sleep(16); // ~60 FPS
    }
    
    cleanup_simulation();
    ReleaseDC(hwnd, hdc);
    return 0;
} 