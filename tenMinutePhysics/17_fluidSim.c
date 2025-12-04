#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <raylib.h>
#define RAYGUI_IMPLEMENTATION
#include "../raygui.h"

#include "../stb_ds.h"

/*
Translated from: Euler Fluid (JavaScript/HTML) by Matthias Müller (Ten Minute Physics)
Target: single-file C program using raylib + raygui for visualization and UI,
and stb_ds.h for any dynamic arrays (not strictly needed here but available).
Notes:
 - Keeps structure/variable names close to original JS where reasonable.
 - Uses simple float arrays allocated with malloc; stb_ds.h is included per request.
 - Compile with: gcc -O2 -o euler_fluid euler_fluid.c -lraylib -lm -lpthread -ldl -lrt
   (Linking may vary by platform. On Windows/MSYS use appropriate raylib linking.)
 - Run and use GUI buttons/checkboxes on screen. Mouse drags place obstacle.
*/

#define U_FIELD 0
#define V_FIELD 1
#define S_FIELD 2

// screen / simulation globals (set after window init)
static int screenWidth = 1280;
static int screenHeight = 720;
static float simHeight = 1.1f;
static float cScale = 1.0f;    // pixels per simulation unit (set in init)
static float simWidth = 1.0f;

static int cnt = 0;

/* Utility coordinate transforms like cX / cY in JS */
static inline float cX(float x) { return x * cScale; }
static inline float cY(float y) { return (float)screenHeight - y * cScale; }

/* Fluid struct mirroring the JS class */
typedef struct
{
    float density;
    int numX;      // includes ghost/solid boundary cells (numX + 2 in JS constructor)
    int numY;
    int numCells;
    float h;
    float *u;
    float *v;
    float *newU;
    float *newV;
    float *p;
    float *s;      // solid flag (1 fluid, 0 solid)
    float *m;      // smoke / marker
    float *newM;
} Fluid;

static Fluid* Fluid_Create(float density, int numX_in, int numY_in, float h) {
    Fluid *f = (Fluid*)malloc(sizeof(Fluid));
    // JS added +2 in constructor; the passed numX/numY were already computed as floor(domainWidth/h)
    // Here we mimic JS: store numX+2, numY+2
    f->density = density;
    f->numX = numX_in + 2;
    f->numY = numY_in + 2;
    f->numCells = f->numX * f->numY;
    f->h = h;
    f->u = (float*)calloc(f->numCells, sizeof(float));
    f->v = (float*)calloc(f->numCells, sizeof(float));
    f->newU = (float*)calloc(f->numCells, sizeof(float));
    f->newV = (float*)calloc(f->numCells, sizeof(float));
    f->p = (float*)calloc(f->numCells, sizeof(float));
    f->s = (float*)calloc(f->numCells, sizeof(float));
    f->m = (float*)calloc(f->numCells, sizeof(float));
    f->newM = (float*)calloc(f->numCells, sizeof(float));
    // initialize marker m to 1.0 (JS: this.m.fill(1.0))
    for (int i = 0; i < f->numCells; ++i) f->m[i] = 1.0f;
    return f;
}

static void Fluid_Destroy(Fluid *f)
{
    if (!f) return;
    free(f->u); free(f->v); free(f->newU); free(f->newV);
    free(f->p); free(f->s); free(f->m); free(f->newM);
    free(f);
}

/* Forward: scene struct (global) */
typedef struct
{
    float gravity;
    float dt;
    int numIters;
    int frameNr;
    float overRelaxation;
    float obstacleX;
    float obstacleY;
    float obstacleRadius;
    bool paused;
    int sceneNr;
    bool showObstacle;
    bool showStreamlines;
    bool showVelocities;
    bool showPressure;
    bool showSmoke;
    float obstacleVX;
    float obstacleVY;
    Fluid *fluid;
} Scene;

static Scene scene;

/* Fluid methods ported from JS */

static void Fluid_integrate(Fluid *f, float dt, float gravity)
{
    int n = f->numY;
    for (int i = 1; i < f->numX; ++i) {
        for (int j = 1; j < f->numY - 1; ++j) {
            if (f->s[i*n + j] != 0.0f && f->s[i*n + j-1] != 0.0f) {
                f->v[i*n + j] += gravity * dt;
            }
        }
    }
}

static void Fluid_solveIncompressibility(Fluid *f, int numIters, float dt)
{
    int n = f->numY;
    float cp = f->density * f->h / dt;

    for (int iter = 0; iter < numIters; ++iter) {
        for (int i = 1; i < f->numX - 1; ++i) {
            for (int j = 1; j < f->numY - 1; ++j) {

                if (f->s[i*n + j] == 0.0f) continue;

                float sx0 = f->s[(i-1)*n + j];
                float sx1 = f->s[(i+1)*n + j];
                float sy0 = f->s[i*n + j-1];
                float sy1 = f->s[i*n + j+1];
                float ssum = sx0 + sx1 + sy0 + sy1;
                if (ssum == 0.0f) continue;

                float div = f->u[(i+1)*n + j] - f->u[i*n + j] +
                            f->v[i*n + j+1] - f->v[i*n + j];

                float p = -div / ssum;
                p *= scene.overRelaxation;
                f->p[i*n + j] += cp * p;

                f->u[i*n + j] -= sx0 * p;
                f->u[(i+1)*n + j] += sx1 * p;
                f->v[i*n + j] -= sy0 * p;
                f->v[i*n + j+1] += sy1 * p;
            }
        }
    }
}

static void Fluid_extrapolate(Fluid *f)
{
    int n = f->numY;
    for (int i = 0; i < f->numX; ++i) {
        f->u[i*n + 0] = f->u[i*n + 1];
        f->u[i*n + f->numY-1] = f->u[i*n + f->numY-2];
    }
    for (int j = 0; j < f->numY; ++j) {
        f->v[0*f->numY + j] = f->v[1*f->numY + j];
        f->v[(f->numX-1)*f->numY + j] = f->v[(f->numX-2)*f->numY + j];
    }
}

static float Fluid_sampleField(Fluid *f, float x, float y, int field)
{
    int n = f->numY;
    float h = f->h;
    float h1 = 1.0f / h;
    float h2 = 0.5f * h;

    if (x < h) x = h;
    if (x > f->numX * h) x = f->numX * h;
    if (y < h) y = h;
    if (y > f->numY * h) y = f->numY * h;

    float dx = 0.0f, dy = 0.0f;
    float *arr = NULL;

    switch (field) {
        case U_FIELD: arr = f->u; dy = h2; break;
        case V_FIELD: arr = f->v; dx = h2; break;
        case S_FIELD: arr = f->m; dx = h2; dy = h2; break;
        default: arr = f->u; break;
    }

    int x0 = (int)floorf((x-dx)*h1);
    if (x0 > f->numX-1) x0 = f->numX-1;
    if (x0 < 0) x0 = 0;
    float tx = ((x-dx) - x0*h) * h1;
    int x1 = x0 + 1;
    if (x1 > f->numX-1) x1 = f->numX-1;

    int y0 = (int)floorf((y-dy)*h1);
    if (y0 > f->numY-1) y0 = f->numY-1;
    if (y0 < 0) y0 = 0;
    float ty = ((y-dy) - y0*h) * h1;
    int y1 = y0 + 1;
    if (y1 > f->numY-1) y1 = f->numY-1;

    float sx = 1.0f - tx;
    float sy = 1.0f - ty;

    float val = sx*sy * arr[x0*n + y0] +
                tx*sy * arr[x1*n + y0] +
                tx*ty * arr[x1*n + y1] +
                sx*ty * arr[x0*n + y1];
    return val;
}

static float Fluid_avgU(Fluid *f, int i, int j)
{
    int n = f->numY;
    float u = (f->u[i*n + j-1] + f->u[i*n + j] +
               f->u[(i+1)*n + j-1] + f->u[(i+1)*n + j]) * 0.25f;
    return u;
}

static float Fluid_avgV(Fluid *f, int i, int j)
{
    int n = f->numY;
    float v = (f->v[(i-1)*n + j] + f->v[i*n + j] +
               f->v[(i-1)*n + j+1] + f->v[i*n + j+1]) * 0.25f;
    return v;
}

static void Fluid_advectVel(Fluid *f, float dt)
{
    // newU/newV <- u/v then advect into them
    memcpy(f->newU, f->u, sizeof(float) * f->numCells);
    memcpy(f->newV, f->v, sizeof(float) * f->numCells);

    int n = f->numY;
    float h = f->h;
    float h2 = 0.5f * h;

    for (int i = 1; i < f->numX; ++i) {
        for (int j = 1; j < f->numY; ++j) {

            cnt++;

            if (f->s[i*n + j] != 0.0f && f->s[(i-1)*n + j] != 0.0f && j < f->numY - 1) {
                float x = i*h;
                float y = j*h + h2;
                float u = f->u[i*n + j];
                float v = Fluid_avgV(f, i, j);
                x = x - dt*u;
                y = y - dt*v;
                u = Fluid_sampleField(f, x, y, U_FIELD);
                f->newU[i*n + j] = u;
            }
            if (f->s[i*n + j] != 0.0f && f->s[i*n + j-1] != 0.0f && i < f->numX - 1) {
                float x = i*h + h2;
                float y = j*h;
                float u = Fluid_avgU(f, i, j);
                float v = f->v[i*n + j];
                x = x - dt*u;
                y = y - dt*v;
                v = Fluid_sampleField(f, x, y, V_FIELD);
                f->newV[i*n + j] = v;
            }
        }
    }

    memcpy(f->u, f->newU, sizeof(float) * f->numCells);
    memcpy(f->v, f->newV, sizeof(float) * f->numCells);
}

static void Fluid_advectSmoke(Fluid *f, float dt)
{
    memcpy(f->newM, f->m, sizeof(float) * f->numCells);

    int n = f->numY;
    float h = f->h;
    float h2 = 0.5f * h;

    for (int i = 1; i < f->numX - 1; ++i) {
        for (int j = 1; j < f->numY - 1; ++j) {
            if (f->s[i*n + j] != 0.0f) {
                float u = (f->u[i*n + j] + f->u[(i+1)*n + j]) * 0.5f;
                float v = (f->v[i*n + j] + f->v[i*n + j+1]) * 0.5f;
                float x = i*h + h2 - dt*u;
                float y = j*h + h2 - dt*v;
                f->newM[i*n + j] = Fluid_sampleField(f, x, y, S_FIELD);
            }
        }
    }
    memcpy(f->m, f->newM, sizeof(float) * f->numCells);
}

/* simulate wrapper */
static void Fluid_simulate(Fluid *f, float dt, float gravity, int numIters)
{
    Fluid_integrate(f, dt, gravity);
    // zero pressure
    for (int i = 0; i < f->numCells; ++i) f->p[i] = 0.0f;
    Fluid_solveIncompressibility(f, numIters, dt);
    Fluid_extrapolate(f);
    Fluid_advectVel(f, dt);
    Fluid_advectSmoke(f, dt);
}

/* Scene setup and utilities ported from JS */

static void setObstacle(float x, float y, bool reset);

/* Setup scene similar to setupScene(sceneNr) in JS */
static void setupScene(int sceneNr)
{
    scene.sceneNr = sceneNr;
    scene.obstacleRadius = 0.15f;
    scene.overRelaxation = 1.9f;
    scene.dt = 1.0f / 60.0f;
    scene.numIters = 40;
    int res = 100;
    if (sceneNr == 0) res = 50;
    else if (sceneNr == 3) res = 200;

    float domainHeight = 1.0f;
    float domainWidth = domainHeight / simHeight * simWidth;
    float h = domainHeight / res;

    int numX = (int)floorf(domainWidth / h);
    int numY = (int)floorf(domainHeight / h);

    float density = 1000.0f;

    if (scene.fluid)
        Fluid_Destroy(scene.fluid);
    scene.fluid = Fluid_Create(density, numX, numY, h);
    Fluid *f = scene.fluid;
    int n = f->numY;

    if (sceneNr == 0) {
        // tank
        for (int i = 0; i < f->numX; ++i) {
            for (int j = 0; j < f->numY; ++j) {
                float s = 1.0f;
                if (i == 0 || i == f->numX-1 || j == 0) s = 0.0f;
                f->s[i*n + j] = s;
            }
        }
        scene.gravity = -9.81f;
        scene.showPressure = true;
        scene.showSmoke = false;
        scene.showStreamlines = false;
        scene.showVelocities = false;
    } else if (sceneNr == 1 || sceneNr == 3) {
        float inVel = 2.0f;
        for (int i = 0; i < f->numX; ++i) {
            for (int j = 0; j < f->numY; ++j) {
                float s = 1.0f;
                if (i == 0 || j == 0 || j == f->numY-1) s = 0.0f;
                f->s[i*n + j] = s;
                if (i == 1) f->u[i*n + j] = inVel;
            }
        }

        float pipeH = 0.1f * f->numY;
        int minJ = (int)floorf(0.5f * f->numY - 0.5f * pipeH);
        int maxJ = (int)floorf(0.5f * f->numY + 0.5f * pipeH);
        for (int j = minJ; j < maxJ; ++j)
            f->m[j] = 0.0f;

        setObstacle(0.4f, 0.5f, true);

        scene.gravity = 0.0f;
        scene.showPressure = false;
        scene.showSmoke = true;
        scene.showStreamlines = false;
        scene.showVelocities = false;

        if (sceneNr == 3) {
            scene.dt = 1.0f / 120.0f;
            scene.numIters = 100;
            scene.showPressure = true;
        }
    } else if (sceneNr == 2) {
        scene.gravity = 0.0f;
        scene.overRelaxation = 1.0f;
        scene.showPressure = false;
        scene.showSmoke = true;
        scene.showStreamlines = false;
        scene.showVelocities = false;
        scene.obstacleRadius = 0.1f;
    }
}

/* Color utilities (getSciColor) */
static void getSciColor(float val, float minVal, float maxVal, unsigned char out[4])
{
    if (val < minVal)
        val = minVal;
    if (val >= maxVal)
        val = maxVal - 0.0001f;
    float d = maxVal - minVal;
    float v = (d == 0.0f) ? 0.5f : (val - minVal) / d;
    float m = 0.25f;
    int num = (int)floorf(v / m);
    float s = (v - num * m) / m;
    float r = 0;
    float g = 0;
    float b = 0;
    switch (num) {
        case 0:
            r = 0.0f;
            g = s;
            b = 1.0f;
            break;
        case 1:
            r = 0.0f;
            g = 1.0f;
            b = 1.0f - s;
            break;
        case 2:
            r = s;
            g = 1.0f;
            b = 0.0f;
            break;
        case 3:
            r = 1.0f;
            g = 1.0f - s;
            b = 0.0f;
            break;
        default:
            r = g = b = 0.0f;
            break;
    }
    out[0] = (unsigned char)(255.0f * r);
    out[1] = (unsigned char)(255.0f * g);
    out[2] = (unsigned char)(255.0f * b);
    out[3] = 255;
}

static void drawScene(Image img, Texture2D tex, Color *px)
{
    Fluid *f = scene.fluid;
    if (!f) return;
    int n = f->numY;
    float h = f->h;
    float cellScale = 1.1f;

    // find min/max pressure
    float minP = f->p[0], maxP = f->p[0];
    for (int i = 0; i < f->numCells; ++i) {
        if (f->p[i] < minP)
            minP = f->p[i];
        if (f->p[i] > maxP)
            maxP = f->p[i];
    }

    // initialize with black
    for (int i = 0; i < screenWidth * screenHeight; ++i)
        px[i] = BLACK;

    for (int i = 0; i < f->numX; ++i) {
        for (int j = 0; j < f->numY; ++j) {
            unsigned char col[4] = { 0, 0, 0, 255 };

            if (scene.showPressure) {
                float p = f->p[i*n + j];
                float s = f->m[i*n + j];
                getSciColor(p, minP, maxP, col);
                if (scene.showSmoke) {
                    int r = col[0] > (unsigned char)(255.0f * s) ? col[0] - (unsigned char)(255.0f * s) : 0;
                    int g = col[1] > (unsigned char)(255.0f * s) ? col[1] - (unsigned char)(255.0f * s) : 0;
                    int b = col[2] > (unsigned char)(255.0f * s) ? col[2] - (unsigned char)(255.0f * s) : 0;
                    col[0] = (unsigned char)r;
                    col[1] = (unsigned char)g;
                    col[2] = (unsigned char)b;
                }
            } else if (scene.showSmoke) {
                float s = f->m[i*n + j];
                if (scene.sceneNr == 2) {
                    getSciColor(s, 0.0f, 1.0f, col);
                } else {
                    unsigned char v = (unsigned char)(255.0f * s);
                    col[0] = col[1] = col[2] = v;
                }
            } else if (f->s[i*n + j] == 0.0f) {
                col[0] = col[1] = col[2] = 0;
            }

            int x = (int)floorf(cX(i * h));
            int y = (int)floorf(cY((j+1) * h));
            int cx = (int)floorf(cScale * cellScale * h) + 1;
            int cy = (int)floorf(cScale * cellScale * h) + 1;

            for (int yi = y; yi < y + cy; ++yi) {
                if (yi < 0 || yi >= screenHeight)
                    continue;
                for (int xi = x; xi < x + cx; ++xi) {
                    if (xi < 0 || xi >= screenWidth) continue;
                    int idx = yi * screenWidth + xi;
                    px[idx].r = col[0];
                    px[idx].g = col[1];
                    px[idx].b = col[2];
                    px[idx].a = col[3];
                }
            }
        }
    }

    // draw cells to an Image buffer for performance
    // push pixels into image and draw
    memcpy(img.data, px, screenWidth*screenHeight*sizeof(Color));

    UpdateTexture(tex, img.data);
    DrawTexture(tex, 0, 0, WHITE);

    // draw velocities
    if (scene.showVelocities) {
        DrawRectangleLines(0,0, screenWidth, screenHeight, BLACK);
        float scale = 0.02f;
        for (int i = 0; i < f->numX; ++i) {
            for (int j = 0; j < f->numY; ++j) {
                float u = f->u[i*n + j];
                float v = f->v[i*n + j];
                float x0 = cX(i * h);
                float x1 = cX(i * h + u * scale);
                float y = cY((j + 0.5f) * h);
                DrawLineV((Vector2){x0,y}, (Vector2){x1,y}, BLACK);

                float x = cX((i + 0.5f) * h);
                float y0 = cY(j * h);
                float y1 = cY(j * h + v * scale);
                DrawLineV((Vector2){x,y0}, (Vector2){x,y1}, BLACK);
            }
        }
    }

    // streamlines
    if (scene.showStreamlines) {
        float segLen = f->h * 0.2f;
        int numSegs = 15;
        for (int i = 1; i < f->numX - 1; i += 5) {
            for (int j = 1; j < f->numY - 1; j += 5) {
                float x = (i + 0.5f) * f->h;
                float y = (j + 0.5f) * f->h;
                Vector2 prev = { cX(x), cY(y) };
                for (int ns = 0; ns < numSegs; ++ns) {
                    float u = Fluid_sampleField(f, x, y, U_FIELD);
                    float v = Fluid_sampleField(f, x, y, V_FIELD);
                    float l = sqrtf(u*u + v*v);
                    x += u * 0.01f;
                    y += v * 0.01f;
                    if (x > f->numX * f->h) break;
                    Vector2 cur = { cX(x), cY(y) };
                    DrawLineV(prev, cur, BLACK);
                    prev = cur;
                }
            }
        }
    }

    // obstacle
    if (scene.showObstacle) {
        float r = scene.obstacleRadius + scene.fluid->h;
        Color fill = scene.showPressure ? BLACK : (Color){0xDD,0xDD,0xDD,255};
        DrawCircle((int)cX(scene.obstacleX), (int)cY(scene.obstacleY), cScale * r+3, BLACK);
        DrawCircleV((Vector2){ cX(scene.obstacleX), cY(scene.obstacleY) }, cScale * r, fill);
//        DrawCircleLines((int)cX(scene.obstacleX), (int)cY(scene.obstacleY), cScale * r, BLACK);
    }

    if (scene.showPressure) {
        char buf[128];
        sprintf(buf, "pressure: %.0f - %.0f N/m", minP, maxP);
        DrawText(buf, 10, 45, 16, BLACK);
    }
}

static void setObstacle(float x, float y, bool reset)
{
    float vx = 0.0f, vy = 0.0f;
    if (!reset) {
        vx = (x - scene.obstacleX) / scene.dt;
        vy = (y - scene.obstacleY) / scene.dt;
    }
    scene.obstacleX = x; scene.obstacleY = y;
    float r = scene.obstacleRadius;
    Fluid *f = scene.fluid;
    int n = f->numY;

    for (int i = 1; i < f->numX - 2; ++i) {
        for (int j = 1; j < f->numY - 2; ++j) {
            f->s[i*n + j] = 1.0f;
            float dx = (i + 0.5f) * f->h - x;
            float dy = (j + 0.5f) * f->h - y;
            if (dx*dx + dy*dy < r*r) {
                f->s[i*n + j] = 0.0f;
                if (scene.sceneNr == 2)
                    f->m[i*n + j] = 0.5f + 0.5f * sinf(0.1f * scene.frameNr);
                else
                    f->m[i*n + j] = 1.0f;
                f->u[i*n + j] = vx;
                f->u[(i+1)*n + j] = vx;
                f->v[i*n + j] = vy;
                f->v[i*n + j+1] = vy;
            }
        }
    }
    scene.showObstacle = true;
}

static void simulateStep()
{
    if (!scene.paused && scene.fluid) {
        Fluid_simulate(scene.fluid, scene.dt, scene.gravity, scene.numIters);
    }
    scene.frameNr++;
}

/* Input handling: mouse drag to place obstacle */
static bool mouseDown = false;
static void StartDrag(float mx, float my)
{
    mouseDown = true;
    float x = mx / cScale;
    float y = ((float)screenHeight - my) / cScale;
    setObstacle(x, y, true);
}
static void Drag(float mx, float my)
{
    if (!mouseDown) return;
    float x = mx / cScale;
    float y = ((float)screenHeight - my) / cScale;
    setObstacle(x, y, false);
}
static void EndDrag() { mouseDown = false; }

int main()
{
    // Init window
    screenWidth = GetScreenWidth();
    screenHeight = GetScreenHeight();
    // If default values not set, set to 1280x720
    if (screenWidth == 0) screenWidth = 1280;
    if (screenHeight == 0) screenHeight = 720;

    InitWindow(screenWidth, screenHeight, "Euler Fluid - translated");
    SetTargetFPS(60);

    // compute scales
    simHeight = 1.1f;
    cScale = (float)screenHeight / simHeight;
    simWidth = (float)screenWidth / cScale;

    // initialize scene defaults
    scene.gravity = -9.81f;
    scene.dt = 1.0f / 120.0f;
    scene.numIters = 100;
    scene.frameNr = 0;
    scene.overRelaxation = 1.9f;
    scene.obstacleX = 0.0f; scene.obstacleY = 0.0f;
    scene.obstacleRadius = 0.15f;
    scene.paused = false;
    scene.sceneNr = 0;
    scene.showObstacle = false;
    scene.showStreamlines = false;
    scene.showVelocities = false;
    scene.showPressure = false;
    scene.showSmoke = true;
    scene.fluid = NULL;

    setupScene(1);

    // UI variables
    bool showStreamlines = scene.showStreamlines;
    bool showVelocities = scene.showVelocities;
    bool showPressure = scene.showPressure;
    bool showSmoke = scene.showSmoke;
    bool overrelax = scene.overRelaxation > 1.0f;

    Image img = {0};
    img.format = PIXELFORMAT_UNCOMPRESSED_R8G8B8A8;
    img.width = screenWidth;
    img.height = screenHeight;
    img.mipmaps = 1;
    img.data = MemAlloc(screenWidth*screenHeight*sizeof(Color));
//    Image img = GenImageColor(screenWidth, screenHeight, BLACK);
//    ImageFormat(&img, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);

    Texture tex = LoadTextureFromImage(img);

    Color *px = (Color*)MemAlloc(screenWidth * screenHeight * sizeof(Color));

    while (!WindowShouldClose()) {
        // mouse/touch handling
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            Vector2 m = GetMousePosition();
            StartDrag(m.x, m.y);
        }
        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) {
            EndDrag();
        }
        if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
            Vector2 m = GetMousePosition();
            Drag(m.x, m.y);
        }

        // keyboard
        if (IsKeyPressed(KEY_P)) scene.paused = !scene.paused;
        if (IsKeyPressed(KEY_M)) { scene.paused = false; simulateStep(); scene.paused = true; }

        simulateStep();
        BeginDrawing(); {
            ClearBackground(BLACK);

            drawScene(img, tex, px);

            // UI (simple top-left controls)
            if (GuiButton((Rectangle){10,10,100,24}, "Wind Tunnel")) { setupScene(1); }
            if (GuiButton((Rectangle){120,10,100,24}, "Hires Tunnel")) { setupScene(3); }
            if (GuiButton((Rectangle){230,10,100,24}, "Tank")) { setupScene(0); }
            if (GuiButton((Rectangle){340,10,100,24}, "Paint")) { setupScene(2); }

            if (GuiToggle((Rectangle){450,10,120,24}, "Streamlines", &scene.showStreamlines)) {
                scene.showStreamlines = !scene.showStreamlines;
            }
            if (GuiToggle((Rectangle){580,10,100,24}, "Velocities", &scene.showVelocities)) {
                scene.showVelocities = !scene.showVelocities;
            }
            if (GuiToggle((Rectangle){690,10,80,24}, "Pressure", &scene.showPressure)) {
                scene.showPressure = !scene.showPressure;
            }
            if (GuiToggle((Rectangle){780,10,80,24}, "Smoke", &scene.showSmoke)) {
                scene.showSmoke = !scene.showSmoke;
            }
            if (GuiToggle((Rectangle){870,10,100,24}, "Overrelax", &overrelax)) {
                overrelax = !overrelax;
                scene.overRelaxation = overrelax ? 1.9f : 1.0f;
            }
        } EndDrawing();
    }

    MemFree(px);
    UnloadTexture(tex);
    UnloadImage(img);

    // cleanup
    if (scene.fluid)
        Fluid_Destroy(scene.fluid);

    CloseWindow();
    return 0;
}
