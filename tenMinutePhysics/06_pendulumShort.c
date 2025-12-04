#include <stdlib.h>
#include <math.h>

#include <raylib.h>
#define RAYMATH_IMPLEMENTATION
#include <raymath.h>

#define STB_DS_IMPLEMENTATION
#include "../stb_ds.h"


float lengths[] = {0.2, 0.2, 0.2};
float masses[] = {1.0, 0.5, 0.3};
float angles[] = {0.5 * M_PI, M_PI, M_PI};

int window_innerWidth;
int window_innerHeight;
int canvas_width;
int canvas_height;
float simMinWidth;
float cScale;

float cX(float x) { return canvas_width / 2 + x * cScale; }
float cY(float y) { return 0.4 * canvas_height - y * cScale; }

typedef struct
{
    float *masses;
    float *lengths;
    double *posX;
    double *posY;
    double *prevPosX;
    double *prevPosY;
    double *velX;
    double *velY;
} Pendulum;
typedef struct
{
    float gravity;
    float dt;
    int numSubSteps;
    Pendulum pendulum;
} Scene;
Scene scene = {0};
Pendulum Pendulum_constructor()
{
    Pendulum this = {0};

    int numMasses = sizeof(masses) / sizeof(masses[0]);

    arrput(this.masses, 0.0);
    arrput(this.lengths, 0.0);
    arrput(this.posX, 0.0);
    arrput(this.posY, 0.0);
    arrput(this.prevPosX, 0.0);
    arrput(this.prevPosY, 0.0);
    arrput(this.velX, 0.0);
    arrput(this.velY, 0.0);
    double x = 0.0, y = 0.0;
    for (int i = 0; i < numMasses; i++) {
        arrput(this.masses, masses[i]);
        arrput(this.lengths, lengths[i]);
        x += lengths[i] * sinf(angles[i]);
        y += lengths[i] * -cosf(angles[i]);
        arrput(this.posX, x);
        arrput(this.posY, y);
        arrput(this.prevPosX, x);
        arrput(this.prevPosY, y);
        arrput(this.velX, 0.0);
        arrput(this.velY, 0.0);
    }

    return this;
}
void Pendulum_simulate(Pendulum *this, float dt, float gravity)
{
    Pendulum *p = this;
    for (int i = 1; i < arrlen(p->masses); i++) {
        p->velY[i] += dt * scene.gravity;
        p->prevPosX[i] = p->posX[i];
        p->prevPosY[i] = p->posY[i];
        p->posX[i] += p->velX[i] * dt;
        p->posY[i] += p->velY[i] * dt;
    }
    for (int i = 1; i < arrlen(p->masses); i++) {
        float dx = p->posX[i] - p->posX[i-1];
        float dy = p->posY[i] - p->posY[i-1];
        float d = sqrt(dx * dx + dy * dy);
        float w0 = p->masses[i - 1] > 0.0 ? 1.0 / p->masses[i - 1] : 0.0;
        float w1 = p->masses[i] > 0.0 ? 1.0 / p->masses[i] : 0.0;
        float corr = (p->lengths[i] - d) / d / (w0 + w1);
        p->posX[i - 1] -= w0 * corr * dx;
        p->posY[i - 1] -= w0 * corr * dy;
        p->posX[i] += w1 * corr * dx;
        p->posY[i] += w1 * corr * dy;
    }
    for (int i = 1; i < arrlen(p->masses); i++) {
        p->velX[i] = (p->posX[i] - p->prevPosX[i]) / dt;
        p->velY[i] = (p->posY[i] - p->prevPosY[i]) / dt;
    }
}
void Pendulum_draw(Pendulum *this)
{
    Pendulum *p = this;
    Color c = {0x30, 0x30, 0x30, 0xFF};
    int lineWidth = 10;
    Vector2 prev = {cX(p->posX[0]), cY(p->posY[0])};
    Vector2 curr = prev;
    for (int i = 0; i < arrlen(p->masses); i++) {
        curr = (Vector2){cX(p->posX[i]), cY(p->posY[i])};
        DrawLineEx(prev, curr, lineWidth, c);
        prev = curr;
    }
    lineWidth = 1;

    c = (Color){0xFF, 0x00, 0x00, 0xFF};
    for (int i = 1; i < arrlen(p->masses); i++) {
        float r = 0.05 * sqrtf(p->masses[i]);
        DrawCircle(cX(p->posX[i]), cY(p->posY[i]), cScale * r, c);
    }
}

void draw()
{
    ClearBackground(BLACK);
    Pendulum_draw(&scene.pendulum);
}

void simulate()
{
    float sdt = scene.dt / scene.numSubSteps;
    for (int step = 0; step < scene.numSubSteps; step++)
        Pendulum_simulate(&scene.pendulum, sdt, scene.gravity);
}

void update()
{
    simulate();
    draw();
}

int main()
{
    window_innerWidth = 1200;
    window_innerHeight = 700;
    canvas_width  = (window_innerWidth - 20);
    canvas_height = (window_innerHeight - 20);
    simMinWidth = 1.0;
    cScale = fminf(canvas_width, canvas_height) / simMinWidth;

    scene = (Scene){
        gravity : -10.0,
        dt : 0.01,
        numSubSteps : 10000,
        pendulum : Pendulum_constructor()
    };

    InitWindow(canvas_width, canvas_height, "pendulum short");
    SetTargetFPS(30);

    while (!WindowShouldClose()) {
        BeginDrawing(); {
            update();
        } EndDrawing();
    }

    arrfree(scene.pendulum.masses);
    arrfree(scene.pendulum.lengths);
    arrfree(scene.pendulum.posX);
    arrfree(scene.pendulum.posY);
    arrfree(scene.pendulum.prevPosX);
    arrfree(scene.pendulum.prevPosY);
    arrfree(scene.pendulum.velX);
    arrfree(scene.pendulum.velY);
    
    CloseWindow();
}
