#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <raylib.h>
#define RAYGUI_IMPLEMENTATION
#include "../raygui.h"

static const int TRAIL_CAP = 1000; // Int32Array(1000) storing x,y pairs
static const double SIM_MIN_WIDTH = 1.0f;

static int winWidth = 800;
static int winHeight = 600;
static double cScale = 1.0f;

static double cX(double x) { return winWidth / 2.0f + x * cScale; }
static double cY(double y) { return 0.4f * winHeight - y * cScale; }

typedef struct
{
    bool usePBD;
    Color color;
    double *masses;   // dynamic arrays, 1-based in JS but we store 0..n-1 with a leading zero mimic
    double *lengths;
    double *posX;
    double *posY;
    double *prevPosX;
    double *prevPosY;
    double *velX;
    double *velY;
    double *theta;
    double *omega;
    int count; // number of masses (including a dummy 0 index to mirror JS array where masses[0]=0)
    int *trail; // int array of length TRAIL_CAP
    int trailFirst;
    int trailLast;
} Pendulum;

typedef struct
{
    double gravity;
    double dt;
    int numSubSteps;
    bool paused;
    Pendulum *pendulumPBD;
    Pendulum *pendulumAnalytic;
} Scene;

static Scene scene = {
        .gravity = -10.0f,
        .dt = 0.01f,
        .numSubSteps = 10000,
        .paused = true,
        .pendulumPBD = NULL,
        .pendulumAnalytic = NULL
};

static int sceneNr = 0;

static Color ColorFromHex(const char *hex)
{
    int r=255,g=0,b=0;
    if (hex && hex[0]=='#' && strlen(hex)>=7) {
        char rs[3]={hex[1],hex[2],0}, gs[3]={hex[3],hex[4],0}, bs[3]={hex[5],hex[6],0};
        r = (int)strtol(rs, NULL, 16);
        g = (int)strtol(gs, NULL, 16);
        b = (int)strtol(bs, NULL, 16);
    }
    return (Color){ r, g, b, 255 };
}

static Pendulum *Pendulum_Create(bool usePBD, const char *hexColor, double *masses_in, int masses_n, double *lengths_in, int lengths_n, double *angles_in, int angles_n)
{
    Pendulum *p = calloc(1, sizeof(Pendulum));
    p->usePBD = usePBD;
    p->color = ColorFromHex(hexColor);

    // js arrays has a leading zero element at index 0, then actual masses at 1..n
    int n = masses_n; // number of actual pendulum links
    p->count = n + 1; // include dummy 0 index
    p->masses = NULL;
    p->lengths = NULL;
    p->posX = NULL;
    p->posY = NULL;
    p->prevPosX = NULL;
    p->prevPosY = NULL;
    p->velX = NULL;
    p->velY = NULL;
    p->theta = NULL;
    p->omega = NULL;

    p->masses = calloc(p->count, sizeof(double));
    p->lengths = calloc(p->count, sizeof(double));
    p->posX = calloc(p->count, sizeof(double));
    p->posY = calloc(p->count, sizeof(double));
    p->prevPosX = calloc(p->count, sizeof(double));
    p->prevPosY = calloc(p->count, sizeof(double));
    p->velX = calloc(p->count, sizeof(double));
    p->velY = calloc(p->count, sizeof(double));
    p->theta = calloc(p->count, sizeof(double));
    p->omega = calloc(p->count, sizeof(double));

    // trail
    p->trail = calloc(TRAIL_CAP, sizeof(int));
    p->trailFirst = 0;
    p->trailLast = 0;

    // initialize leading zero as in js
    p->masses[0] = 0.0f;
    p->lengths[0] = 0.0f;
    p->posX[0] = 0.0f;
    p->posY[0] = 0.0f;
    p->prevPosX[0] = 0.0f;
    p->prevPosY[0] = 0.0f;
    p->velX[0] = 0.0f;
    p->velY[0] = 0.0f;
    p->theta[0] = 0.0f;
    p->omega[0] = 0.0f;

    double x = 0.0f, y = 0.0f;
    for (int i = 0; i < n; ++i) {
        p->masses[i+1] = masses_in[i];
        p->lengths[i+1] = lengths_in[i];
        p->theta[i+1] = angles_in[i];
        p->omega[i+1] = 0.0f;
        x += lengths_in[i] * sinf(angles_in[i]);
        y += lengths_in[i] * -cosf(angles_in[i]);
        p->posX[i+1] = x;
        p->posY[i+1] = y;
        p->prevPosX[i+1] = x;
        p->prevPosY[i+1] = y;
        p->velX[i+1] = 0.0f;
        p->velY[i+1] = 0.0f;
    }

    return p;
}

static void Pendulum_Free(Pendulum *p)
{
    if (!p) return;
    free(p->masses);
    free(p->lengths);
    free(p->posX);
    free(p->posY);
    free(p->prevPosX);
    free(p->prevPosY);
    free(p->velX);
    free(p->velY);
    free(p->theta);
    free(p->omega);
    free(p->trail);
    free(p);
}

// simulate analytic (only works for up to 3 masses as in JS)
static void Pendulum_simulateAnalytic(Pendulum *p, double dt, double gravity)
{
    // only valid when count-1 <= 3 (masses length <= 3)
    if (p->count-1 > 3)
        return;
    double g = -gravity;
    double m1 = p->masses[1];
    double m2 = p->masses[2];
    double m3 = p->masses[3];
    double l1 = p->lengths[1];
    double l2 = p->lengths[2];
    double l3 = p->lengths[3];
    double t1 = p->theta[1];
    double t2 = p->theta[2];
    double t3 = p->theta[3];
    double w1 = p->omega[1];
    double w2 = p->omega[2];
    double w3 = p->omega[3];

    double b1 =
        g*l1*m1*sinf(t1) + g*l1*m2*sinf(t1)+g*l1*m3*sinf(t1) + m2*l1*l2*sinf(t1-t2)*w1*w2 +
        m3*l1*l3*sinf(t1-t3)*w1*w3       +   m3*l1*l2*sinf(t1-t2)*w1*w2  +
        m2*l1*l2*sinf(t2-t1)*(w1-w2)*w2  +
        m3*l1*l2*sinf(t2-t1)*(w1-w2)*w2  +
        m3*l1*l3*sinf(t3-t1)*(w1-w3)*w3;

    double a11 = l1*l1*(m1+m2+m3);
    double a12 = m2*l1*l2*cosf(t1-t2) + m3*l1*l2*cosf(t1-t2);
    double a13 = m3*l1*l3*cosf(t1-t3);

    double b2 =
        g*l2*m2*sinf(t2) + g*l2*m3*sinf(t2) + w1*w2*l1*l2*sinf(t2-t1)*(m2 + m3) +
        m3*l2*l3*sinf(t2-t3)*w2*w3               +
        (m2 + m3)*l1*l2*sinf(t2-t1)*(w1-w2)*w1   +
        m3*l2*l3*sinf(t3-t2)*(w2-w3)*w3;

    double a21 = (m2 + m3)*l1*l2*cosf(t2-t1);
    double a22 = l2*l2*(m2+m3);
    double a23 = m3*l2*l3*cosf(t2-t3);

    double b3 =
        m3*g*l3*sinf(t3) - m3*l2*l3*sinf(t2-t3)*w2*w3 - m3*l1*l3*sinf(t1-t3)*w1*w3 +
        m3*l1*l3*sinf(t3-t1)*(w1-w3)*w1    +
        m3*l2*l3*sinf(t3-t2)*(w2-w3)*w2;

    double a31 = m3*l1*l3*cosf(t1-t3);
    double a32 = m3*l2*l3*cosf(t2-t3);
    double a33 = m3*l3*l3;

    b1 = -b1;
    b2 = -b2;
    b3 = -b3;

    double det = a11 * (a22 * a33 - a23 * a32) + a21 * (a32 * a13 - a33 * a12) + a31 * (a12 * a23 - a13 * a22);
    if (det == 0.0f)
        return;

    double A1 = b1 * (a22 * a33 - a23 * a32) + b2 * (a32 * a13 - a33 * a12) + b3 * (a12 * a23 - a13 * a22);
    double A2 = b1 * (a23 * a31 - a21 * a33) + b2 * (a33 * a11 - a31 * a13) + b3 * (a13 * a21 - a11 * a23);
    double A3 = b1 * (a21 * a32 - a22 * a31) + b2 * (a31 * a12 - a32 * a11) + b3 * (a11 * a22 - a12 * a21);

    A1 /= det;
    A2 /= det;
    A3 /= det;

    p->omega[1] += A1 * dt;
    p->omega[2] += A2 * dt;
    p->omega[3] += A3 * dt;
    p->theta[1] += p->omega[1] * dt;
    p->theta[2] += p->omega[2] * dt;
    p->theta[3] += p->omega[3] * dt;

    double x=0.0f, y=0.0f;
    for (int i = 1; i < p->count; i++) {
        x += p->lengths[i] * sinf(p->theta[i]);
        y += p->lengths[i] * -cosf(p->theta[i]);
        p->posX[i] = x;
        p->posY[i] = y;
    }
}

static void Pendulum_updateTrail(Pendulum *p)
{
    int idx = p->trailLast;
    // store as int of screen coords
    int sx = (int)roundf(cX(p->posX[p->count-1]));
    int sy = (int)roundf(cY(p->posY[p->count-1]));
    p->trail[idx] = sx;
    p->trail[(idx+1)%TRAIL_CAP] = sy;
    p->trailLast = (p->trailLast + 2) % TRAIL_CAP;
    if (p->trailLast == p->trailFirst)
        p->trailFirst = (p->trailFirst + 2) % TRAIL_CAP;
}

static void Pendulum_simulatePBD(Pendulum *p, double dt, double gravity)
{
    // integrate velocities and positions
    for (int i = 1; i < p->count; ++i) {
        p->velY[i] += dt * gravity;
        p->prevPosX[i] = p->posX[i];
        p->prevPosY[i] = p->posY[i];
        p->posX[i] += p->velX[i] * dt;
        p->posY[i] += p->velY[i] * dt;
    }
    // constraints
    for (int i = 1; i < p->count; ++i) {
        double dx = p->posX[i] - p->posX[i-1];
        double dy = p->posY[i] - p->posY[i-1];
        double d = sqrtf(dx*dx + dy*dy);
        double w0 = (p->masses[i-1] > 0.0f) ? 1.0f / p->masses[i-1] : 0.0f;
        double w1 = (p->masses[i] > 0.0f) ? 1.0f / p->masses[i] : 0.0f;
        if (d < 1e-8f)
            d = 1e-8f;
        double corr = (p->lengths[i] - d) / d / (w0 + w1);
        p->posX[i-1] -= w0 * corr * dx;
        p->posY[i-1] -= w0 * corr * dy;
        p->posX[i] += w1 * corr * dx;
        p->posY[i] += w1 * corr * dy;
    }
    // update velocities
    for (int i = 1; i < p->count; ++i) {
        p->velX[i] = (p->posX[i] - p->prevPosX[i]) / dt;
        p->velY[i] = (p->posY[i] - p->prevPosY[i]) / dt;
    }
}

static void Pendulum_simulate(Pendulum *p, double dt, double gravity)
{
    if (p->usePBD)
        Pendulum_simulatePBD(p, dt, gravity);
    else
        Pendulum_simulateAnalytic(p, dt, gravity);
}

static void Pendulum_draw(Pendulum *p)
{
    // draw trail
    if (p->trailLast != p->trailFirst) {
        int i = p->trailFirst;
        Vector2 prev = { (double)p->trail[i], (double)p->trail[i+1] };
        i = (i + 2) % TRAIL_CAP;
        while (i != p->trailLast) {
            Vector2 cur = { (double)p->trail[i], (double)p->trail[(i+1)%TRAIL_CAP] };
            DrawLineV(prev, cur, p->color);
            prev = cur;
            i = (i + 2) % TRAIL_CAP;
        }
    }
    // draw rods
    DrawLineEx((Vector2){ cX(p->posX[0]), cY(p->posY[0]) }, (Vector2){ cX(p->posX[1]), cY(p->posY[1]) }, 10.0f, DARKGRAY);
    for (int i = 1; i < p->count; ++i) {
        Vector2 a = { cX(p->posX[i-1]), cY(p->posY[i-1]) };
        Vector2 b = { cX(p->posX[i]), cY(p->posY[i]) };
        DrawLineEx(a, b, 10.0f, DARKGRAY);
    }
    // draw masses
    for (int i = 1; i < p->count; ++i) {
        double r = 0.03f * sqrtf(p->masses[i]);
        DrawCircleV((Vector2){ cX(p->posX[i]), cY(p->posY[i]) }, cScale * r, p->color);
    }
}

// --- Scene setup ---

static void setupScene()
{
    double anglesArr[] = { 0.5f * (double)M_PI, (double)M_PI, (double)M_PI, (double)M_PI, (double)M_PI };
    double lengths[5];
    double masses[5];
    int useN = 0;

    switch(sceneNr % 6) {
        case 0: {
            double l[] = {0.15f,0.15f,0.15f};
            double m[] = {1.0f,1.0f,1.0f};
            memcpy(lengths, l, sizeof(l));
            memcpy(masses, m, sizeof(m));
            useN=3;
            break;
        }
        case 1: {
            double l[] = {0.06f,0.15f,0.2f};
            double m[] = {1.0f,0.5f,0.1f};
            memcpy(lengths, l, sizeof(l));
            memcpy(masses, m, sizeof(m));
            useN=3;
            break;
        }
        case 2: {
            double l[] = {0.15f,0.15f,0.15f};
            double m[] = {1.0f,0.01f,1.0f};
            memcpy(lengths, l, sizeof(l));
            memcpy(masses, m, sizeof(m));
            useN=3;
            break;
        }
        case 3: {
            double l[] = {0.15f,0.15f,0.15f};
            double m[] = {0.01f,1.0f,0.01f};
            memcpy(lengths, l, sizeof(l));
            memcpy(masses, m, sizeof(m));
            useN=3;
            break;
        }
        case 4: {
            double l[] = {0.2f,0.133f,0.04f};
            double m[] = {0.3f,0.3f,0.3f};
            memcpy(lengths, l, sizeof(l));
            memcpy(masses, m, sizeof(m));
            useN=3;
            break;
        }
        case 5: {
            double l[] = {0.1f,0.12f,0.1f,0.15f,0.05f};
            double m[] = {0.2f,0.6f,0.4f,0.3f,0.2f};
            memcpy(lengths, l, sizeof(l));
            memcpy(masses, m, sizeof(m));
            useN=5;
            break;
        }
    }

    if (scene.pendulumPBD) {
        Pendulum_Free(scene.pendulumPBD);
        scene.pendulumPBD = NULL;
    }
    if (scene.pendulumAnalytic) {
        Pendulum_Free(scene.pendulumAnalytic);
        scene.pendulumAnalytic = NULL;
    }

    // create pendulums
    scene.pendulumPBD = Pendulum_Create(true, "#FF3030", masses, useN, lengths, useN, anglesArr, useN);
    if (useN <= 3) {
        scene.pendulumAnalytic = Pendulum_Create(false, "#00FF00", masses, useN, lengths, useN, anglesArr, useN);
    } else {
        scene.pendulumAnalytic = NULL;
    }

    scene.paused = true;
    sceneNr++;
}

static void drawScene()
{
    ClearBackground(BLACK);

    if (scene.pendulumPBD)
        Pendulum_draw(scene.pendulumPBD);
    if (scene.numSubSteps >= 100) {
        if (scene.pendulumAnalytic)
            Pendulum_draw(scene.pendulumAnalytic);
    }
}

static void simulateStep()
{
    if (scene.paused) return;
    float sdt = scene.dt / (float)scene.numSubSteps;
    if (sdt < 1e-6f) sdt = 1e-6f;
    int trailGap = scene.numSubSteps / 10;
    if (trailGap < 1)
        trailGap = 1;

    for (int step = 0; step < scene.numSubSteps; ++step) {
        if (scene.pendulumPBD) {
            Pendulum_simulate(scene.pendulumPBD, sdt, scene.gravity);
            if ((step % trailGap) == 0) Pendulum_updateTrail(scene.pendulumPBD);
        }
        if (scene.pendulumAnalytic) {
            Pendulum_simulate(scene.pendulumAnalytic, sdt, scene.gravity);
            if ((step % trailGap) == 0) Pendulum_updateTrail(scene.pendulumAnalytic);
        }
    }
}

// --- UI state ---
static int stepsSliderValue = 5; // corresponds to JS default value = 5
static const int stepsOptions[6] = {1,5,10,100,1000,10000};
static Rectangle btnRestart = { 10, 10, 120, 30 };
static Rectangle btnRun = { 140, 10, 80, 30 };
static Rectangle btnStep = { 230, 10, 80, 30 };
static Rectangle sliderRect = { 60, 50, 250, 20 };

int main()
{
    winWidth = 1000;
    winHeight = 700;
    InitWindow(winWidth, winHeight, "Pendulum");
    SetTargetFPS(30);

    cScale = (double)fmin(winWidth, winHeight) / SIM_MIN_WIDTH;

    setupScene();

    scene.numSubSteps = stepsOptions[stepsSliderValue];

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_S)) {
            scene.paused = false;
            simulateStep();
            scene.paused = true;
        }

        if (GuiButton(btnRestart, "Restart")) {
            setupScene();
        }
        if (GuiButton(btnRun, "Run")) {
            scene.paused = false;
        }
        if (GuiButton(btnStep, "Step")) {
            // step once as above
            bool prevPaused = scene.paused;
            scene.paused = false;
            simulateStep();
            scene.paused = true;
            scene.paused = prevPaused;
            scene.paused = true;
        }

        float stepsSliderValue_f = (float)stepsSliderValue;
        int newSlider = GuiSliderBar(sliderRect, "Sub-steps", TextFormat("%d", stepsOptions[stepsSliderValue]), &stepsSliderValue_f, 0, 5);
        if (newSlider) {
            stepsSliderValue = (int)stepsSliderValue_f;
            scene.numSubSteps = stepsOptions[stepsSliderValue];
        }

        simulateStep();

        BeginDrawing(); {
            drawScene();

            DrawText(TextFormat("Number of sub-steps: %d", scene.numSubSteps), 10, 80, 14, RAYWHITE);
            DrawText("Press 's' to Step", 10, 100, 12, GRAY);
        } EndDrawing();
    }

    if (scene.pendulumPBD)
        Pendulum_Free(scene.pendulumPBD);
    if (scene.pendulumAnalytic)
        Pendulum_Free(scene.pendulumAnalytic);

    CloseWindow();

    return 0;
}
