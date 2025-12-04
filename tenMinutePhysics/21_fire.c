#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

#include <raylib.h>
#define RAYGUI_IMPLEMENTATION
#include "../raygui.h"

typedef struct
{
    int numX;
    int numY;
    int numCells;
    float h;
    float *u;
    float *v;
    float *newU;
    float *newV;
    float *s;
    float *t;
    float *newT;
    int num;

    int numSwirls;

    int maxSwirls;
    float swirlGlobalTime;
    float *swirlX;
    float *swirlY;
    float *swirlOmega;
    float *swirlRadius;
    float *swirlTime;
} Fluid;
typedef struct
{
    float gravity;
    float dt;
    int numIters;
    int frameNr;
    float obstacleX;
    float obstacleY;
    float obstacleRadius;
    bool burningObstacle;
    bool burningFloor;
    bool paused;
    bool showPressure;
    bool showObstacle;
    bool showSwirls;
    float swirlProbability;
    float swirlMaxRadius;
    Fluid fluid;
} Scene;
Scene scene = {
        gravity : 0.00,
        dt : 1.0 / 60.0,
        numIters : 10,
        frameNr : 0,
        obstacleX : 0.0,
        obstacleY : 0.0,
        obstacleRadius: 0.2,
        burningObstacle: true,
        burningFloor: false,
        paused: false,
        showObstacle: false,
        showSwirls: false,
        swirlProbability: 50.0,
        swirlMaxRadius: 0.05,
        fluid: {0}
};


static int windowWidth;
static int windowHeight;

float cScale = 1.0;

#define U_FIELD 0
#define V_FIELD 1
#define T_FIELD 2

int cnt = 0;

static inline float cX(float x)
{
    return x * cScale;
}

static inline float cY(float y)
{
    return windowHeight - y * cScale;
}

float kernel(float r, float rmax, float leftBorder, float rightBorder)
{
    if (r >= rmax || r <= 0.0)
        return 0.0;
    else if (r < leftBorder)
        return r / leftBorder;
    else if (r <= rmax - rightBorder)
        return 1.0;
    else
        return (rmax - r) / rightBorder;
}


// ----------------- start of simulator ------------------------------

Fluid Fluid_constructor(int numX, int numY, float h)
{
    Fluid this = {0};

    this.numX = numX + 2;
    this.numY = numY + 2;
    this.numCells = this.numX * this.numY;
    this.h = h;
    this.u = calloc(this.numCells, sizeof(float));
    this.v = calloc(this.numCells, sizeof(float));
    this.newU = calloc(this.numCells, sizeof(float));
    this.newV = calloc(this.numCells, sizeof(float));
    this.s = calloc(this.numCells, sizeof(float));
    this.t = calloc(this.numCells, sizeof(float));
    this.newT = calloc(this.numCells, sizeof(float));
    this.num = numX * numY;
    for (int i = 0; i < this.numCells; ++i) {
        this.t[i] = 0.0;
        this.s[i] = 1.0;
    }

    this.numSwirls = 0;

    this.maxSwirls = 100;
    this.swirlGlobalTime = 0.0;
    this.swirlX = calloc(this.maxSwirls, sizeof(float));
    this.swirlY = calloc(this.maxSwirls, sizeof(float));
    this.swirlOmega = calloc(this.maxSwirls, sizeof(float));
    this.swirlRadius = calloc(this.maxSwirls, sizeof(float));
    this.swirlTime = calloc(this.maxSwirls, sizeof(float));
    for (int i = 0; i < this.maxSwirls; ++i) {
        this.swirlTime[i] = 0.0;
    }

    return this;
}

void Fluid_dispose(Fluid *this)
{
    free(this->u);
    free(this->v);
    free(this->newU);
    free(this->newV);
    free(this->s);
    free(this->t);
    free(this->newT);

    free(this->swirlX);
    free(this->swirlY);
    free(this->swirlOmega);
    free(this->swirlRadius);
    free(this->swirlTime);
}

void Fluid_integrate(Fluid *this, float dt, float gravity)
{
    int n = this->numY;
    for (int i = 1; i < this->numX; i++) {
        for (int j = 1; j < this->numY-1; j++) {
            if (this->s[i*n + j] != 0.0 && this->s[i*n + j-1] != 0.0)
                this->v[i*n + j] += gravity * dt;
        }
    }
}

void Fluid_solveIncompressibility(Fluid *this, int numIters, float dt)
{
    int n = this->numY;
    float overRelaxation = 1.9;

    for (int iter = 0; iter < numIters; iter++) {

        for (int i = 1; i < this->numX-1; i++) {
            for (int j = 1; j < this->numY-1; j++) {

                if (this->s[i*n + j] == 0.0)
                    continue;

                float sx0 = this->s[(i-1)*n + j];
                float sx1 = this->s[(i+1)*n + j];
                float sy0 = this->s[i*n + j-1];
                float sy1 = this->s[i*n + j+1];
                float s = sx0 + sx1 + sy0 + sy1;
                if (s == 0.0)
                    continue;

                float div = this->u[(i+1)*n + j] - this->u[i*n + j] +
                        this->v[i*n + j+1] - this->v[i*n + j];

                float p = -div / s;
                p *= overRelaxation;
                this->u[i*n + j] -= sx0 * p;
                this->u[(i+1)*n + j] += sx1 * p;
                this->v[i*n + j] -= sy0 * p;
                this->v[i*n + j+1] += sy1 * p;
            }
        }
    }
}

void Fluid_extrapolate(Fluid *this)
{
    int n = this->numY;
    for (int i = 0; i < this->numX; i++) {
        this->u[i*n + 0] = this->u[i*n + 1];
        this->u[i*n + this->numY-1] = this->u[i*n + this->numY-2];
    }
    for (int j = 0; j < this->numY; j++) {
        this->v[0*n + j] = this->v[1*n + j];
        this->v[(this->numX-1)*n + j] = this->v[(this->numX-2)*n + j];
    }
}

float Fluid_sampleField(Fluid *this, float x, float y, int field)
{
    int n = this->numY;
    float h = this->h;
    float h1 = 1.0 / h;
    float h2 = 0.5 * h;

    x = fmaxf(fminf(x, this->numX * h), h);
    y = fmaxf(fminf(y, this->numY * h), h);

    float dx = 0.0;
    float dy = 0.0;

    float *f = NULL;

    switch (field) {
        case U_FIELD: f = this->u; dy = h2; break;
        case V_FIELD: f = this->v; dx = h2; break;
        case T_FIELD: f = this->t; dx = h2; dy = h2; break;
    }

    int x0 = fminf(floorf((x-dx)*h1), this->numX-1);
    float tx = ((x-dx) - x0*h) * h1;
    int x1 = fminf(x0 + 1, this->numX-1);

    int y0 = fminf(floorf((y-dy)*h1), this->numY-1);
    float ty = ((y-dy) - y0*h) * h1;
    int y1 = fminf(y0 + 1, this->numY-1);

    float sx = 1.0 - tx;
    float sy = 1.0 - ty;

    float val = sx*sy * f[x0*n + y0] +
            tx*sy * f[x1*n + y0] +
            tx*ty * f[x1*n + y1] +
            sx*ty * f[x0*n + y1];

    return val;
}

float Fluid_avgU(Fluid *this, int i, int j)
{
    int n = this->numY;
    float u = (this->u[i*n + j-1] + this->u[i*n + j] +
            this->u[(i+1)*n + j-1] + this->u[(i+1)*n + j]) * 0.25;
    return u;
}

float Fluid_avgV(Fluid *this, int i, int j)
{
    int n = this->numY;
    float v = (this->v[(i-1)*n + j] + this->v[i*n + j] +
            this->v[(i-1)*n + j+1] + this->v[i*n + j+1]) * 0.25;
    return v;
}

void Fluid_advectVel(Fluid *this, float dt)
{
    memcpy(this->newU, this->u, this->numCells * sizeof(float));
    memcpy(this->newV, this->v, this->numCells * sizeof(float));

    int n = this->numY;
    float h = this->h;
    float h2 = 0.5 * h;

    for (int i = 1; i < this->numX; i++) {
        for (int j = 1; j < this->numY; j++) {
            cnt++;

            // u component
            if (this->s[i*n + j] != 0.0 && this->s[(i-1)*n + j] != 0.0 && j < this->numY - 1) {
                float x = i*h;
                float y = j*h + h2;
                float u = this->u[i*n + j];
                float v = Fluid_avgV(this, i, j);
                x = x - dt*u;
                y = y - dt*v;
                u = Fluid_sampleField(this, x,y, U_FIELD);
                this->newU[i*n + j] = u;
            }
            // v component
            if (this->s[i*n + j] != 0.0 && this->s[i*n + j-1] != 0.0 && i < this->numX - 1) {
                float x = i*h + h2;
                float y = j*h;
                float u = Fluid_avgU(this, i, j);
                float v = this->v[i*n + j];
                x = x - dt*u;
                y = y - dt*v;
                v = Fluid_sampleField(this, x,y, V_FIELD);
                this->newV[i*n + j] = v;
            }
        }
    }

    memcpy(this->u, this->newU, this->numCells * sizeof(float));
    memcpy(this->v, this->newV, this->numCells * sizeof(float));
}

void Fluid_advectTemperature(Fluid *this, float dt)
{
    memcpy(this->newT, this->t, this->numCells * sizeof(float));

    int n = this->numY;
    float h = this->h;
    float h2 = 0.5 * h;

    for (int i = 1; i < this->numX-1; i++) {
        for (int j = 1; j < this->numY-1; j++) {
            if (this->s[i*n + j] != 0.0) {
                float u = (this->u[i*n + j] + this->u[(i+1)*n + j]) * 0.5;
                float v = (this->v[i*n + j] + this->v[i*n + j+1]) * 0.5;
                float x = i*h + h2 - dt*u;
                float y = j*h + h2 - dt*v;

                this->newT[i*n + j] = Fluid_sampleField(this, x,y, T_FIELD);
            }
        }
    }
    memcpy(this->t, this->newT, this->numCells * sizeof(float));
}

void Fluid_updateFire(Fluid *this, float dt)
{
    float h = this->h;

    float swirlTimeSpan = 1.0;
    float swirlOmega = 20.0;
    float swirlDamping = 10.0 * dt;
    float swirlProbability = scene.swirlProbability * h * h;

    float fireCooling = 1.2 * dt;
    float smokeCooling = 0.3 * dt;
    float lift = 3.0;
    float acceleration = 6.0 * dt;
    float kernelRadius = scene.swirlMaxRadius;

    // update swirls

    int n = this->numY;
    float maxX = (this->numX - 1) * this->h;
    float maxY = (this->numY - 1) * this->h;

    // kill swirls

    int num = 0;
    for (int nr = 0; nr < this->numSwirls; nr++) {
        this->swirlTime[nr] -= dt;
        if (this->swirlTime[nr] > 0.0) {
            this->swirlTime[num] = this->swirlTime[nr];
            this->swirlX[num] = this->swirlX[nr];
            this->swirlY[num] = this->swirlY[nr];
            this->swirlOmega[num] = this->swirlOmega[nr];
            num++;
        }
    }
    this->numSwirls = num;

    // advect and modify velocity field

    for (int nr = 0; nr < this->numSwirls; nr++) {
//        float ageScale = this->swirlTime[nr] / swirlTimeSpan;
        float x = this->swirlX[nr];
        float y = this->swirlY[nr];
        float swirlU = (1.0 - swirlDamping) * Fluid_sampleField(this, x, y, U_FIELD);
        float swirlV = (1.0 - swirlDamping) * Fluid_sampleField(this, x, y, V_FIELD);
        x += swirlU * dt;
        y += swirlV * dt;
        x = fminf(fmaxf(x, h), maxX);
        y = fminf(fmaxf(y, h), maxY);

        this->swirlX[nr] = x;
        this->swirlY[nr] = y;
        float omega = this->swirlOmega[nr];

        // update surrounding velocity field

        int x0 = fmaxf(floorf((x - kernelRadius) / h), 0);
        int y0 = fmaxf(floorf((y - kernelRadius) / h), 0);
        int x1 = fminf(floorf((x + kernelRadius) / h) + 1, this->numX - 1);
        int y1 = fminf(floorf((y + kernelRadius) / h) + 1, this->numY - 1);

        for (int i = x0; i <= x1; i++) {
            for (int j = y0; j <= y1; j++) {
                for (int dim = 0; dim < 2; dim++) {
                    float vx = dim == 0 ? i * h : (i + 0.5) * h;
                    float vy = dim == 0 ? (j + 0.5) * h : j * h;

                    float rx = vx - x;
                    float ry = vy - y;
                    float r = sqrtf(rx * rx + ry * ry);

                    if (r < kernelRadius) {
                        float s = 1.0;
                        if (r > 0.8 * kernelRadius) {
                            s = 5.0 - 5.0 / kernelRadius * r;
                            // s = (kernelRadius - r) / kernelRadius * ageScale;
                        }

                        if (dim == 0) {
                            float target = ry * omega + swirlU;
                            float u = this->u[n*i + j];
                            this->u[n*i + j] = (target - u) * s;
                        } else {
                            float target = -rx * omega + swirlV;
                            float v = this->v[n*i + j];
                            this->v[n*i + j] += (target - v) * s;
                        }
                    }
                }
            }
        }
    }

    // update temperatures

    float minR = 0.85 * scene.obstacleRadius;
    float maxR = scene.obstacleRadius + h;

    for (int i = 0; i < this->numX; i++) {
        for (int j = 0; j < this->numY; j++) {
            float t = this->t[i*n + j];

            float cooling = t < 0.3 ? smokeCooling : fireCooling;
            this->t[i*n + j] = fmaxf(t - cooling, 0.0);
//            float u = this->u[i*n + j];
            float v = this->v[i*n + j];
            float targetV = t * lift;
            this->v[i*n + j] += (targetV - v) * acceleration;

            int numNewSwirls = 0;

            // obstacle burning
            if (scene.burningObstacle) {
                float dx = (i + 0.5) * this->h - scene.obstacleX;
                float dy = (j + 0.5) * this->h - scene.obstacleY - 3.0 * this->h;
                float d = dx * dx + dy * dy;
                if (minR * minR <= d && d < maxR * maxR) {
                    this->t[i*n + j] = 1.0;
                    if (drand48() < 0.5 * swirlProbability)
                        numNewSwirls++;
                }
            }

            // floor burning
            if (j < 4 && scene.burningFloor) {
                this->t[i*n + j] = 1.0;
                this->u[i*n + j] = 0.0;
                this->v[i*n + j] = 0.0;
                if (drand48() < swirlProbability)
                    numNewSwirls++;
            }

            for (int k = 0; k < numNewSwirls; k++) {
                if (this->numSwirls >= this->maxSwirls)
                    break;
                int nr = this->numSwirls;
                this->swirlX[nr] = i * h;
                this->swirlY[nr] = j * h;
                this->swirlOmega[nr] = (-1.0 + 2.0 * drand48()) * swirlOmega;
                this->swirlTime[nr] = swirlTimeSpan;
                this->numSwirls++;
            }
        }
    }

    // foo: maybe for obstacle

    // smooth temperatures

    for (int i = 1; i < this->numX - 1; i++) {
        for (int j = 1; j < this->numY - 1; j++) {
            float t = this->t[i * n + j];
            if (t == 1.0) {
                float avg = (
                        this->t[(i - 1) * n + (j - 1)] +
                        this->t[(i + 1) * n + (j - 1)] +
                        this->t[(i + 1) * n + (j + 1)] +
                        this->t[(i - 1) * n + (j + 1)]) * 0.25;
                this->t[i * n + j] = avg;
            }
        }
    }
}

// ----------------- end of simulator ------------------------------

void Fluid_simulate(Fluid *this, float dt, float gravity, int numIters)
{
    Fluid_integrate(this, dt, gravity);

    Fluid_solveIncompressibility(this, numIters, dt);
    Fluid_extrapolate(this);
    Fluid_advectVel(this, dt);
    Fluid_advectTemperature(this, dt);
    Fluid_updateFire(this, dt);
}

void setupScene()
{
//    var canvas = document.getElementById("myCanvas");
    windowWidth = 1200;
    windowHeight = 700;

    float simHeight = 1.0;
    cScale = windowHeight / simHeight;
    float simWidth = windowWidth / cScale;

    int numCells = 100000;
    float h = sqrtf(simWidth * simHeight / numCells);

    int numX = floorf(simWidth / h);
    int numY = floorf(simHeight / h);

    if (numX < numY) {
        scene.swirlProbability = 80.0;
        scene.swirlMaxRadius = 0.04;
    }

    scene.obstacleX = 0.5 * numX * h;
    scene.obstacleY = 0.3 * numY * h;
    scene.showObstacle = scene.burningObstacle;

    scene.fluid = Fluid_constructor(numX, numY, h);
}


// draw -------------------------------------------------------

void getFireColor(float val, unsigned char out[4])
{
    val = fminf(fmaxf(val, 0.0), 1.0);
    float r, g, b;

    if (val < 0.3) {
        float s = val / 0.3;
        r = 0.2 * s; g = 0.2 * s; b = 0.2 * s;
    }
    else if (val < 0.5) {
        float s = (val-0.3) / 0.2;
        r = 0.2 + 0.8 * s; g = 0.1; b = 0.1;
    }
    else  {
        float s = (val - 0.5) / 0.48;
        r = 1.0; g = s; b = 0.0;
    }
    out[0] = 255*r;
    out[1] = 255*g;
    out[2] = 255*b;
    out[3] = 255;
}

void draw(Image img, Texture2D tex, Color *px)
{
    ClearBackground(BLACK);

    Color c = {0xFF, 0x00, 0x00, 255};
    Fluid *f = &scene.fluid;
    int n = f->numY;

    float cellScale = 1.1;

    float h = f->h;

    unsigned char color[4] = { 255, 255, 255, 255 };

    for (int i = 0; i < f->numX; i++) {
        for (int j = 0; j < f->numY; j++) {
            float t = f->t[i*n + j];
            getFireColor(t, color);

            int x = floorf(cX(i * h));
            int y = floorf(cY((j+1) * h));
            int cx = (int)ceilf(cScale * cellScale * h) + 1;
            int cy = (int)ceilf(cScale * cellScale * h) + 1;

            unsigned char r = color[0];
            unsigned char g = color[1];
            unsigned char b = color[2];

            for (int yi = y; yi < y + cy; yi++) {
                if (yi < 0 || yi >= windowHeight) continue;
                for (int xi = 0; xi < cx; xi++) {
                    int xx = x + xi;
                    if (xx < 0 || xx >= windowWidth) continue;
                    int p = yi * windowWidth + xx;
                    px[p] = (Color){r, g, b, 255};
                }
            }
        }
    }

    // draw cells to an Image buffer for performance
    // push pixels into image and draw
    memcpy(img.data, px, windowWidth*windowHeight*sizeof(Color));

    UpdateTexture(tex, img.data);
    DrawTexture(tex, 0, 0, WHITE);

    if (scene.showObstacle) {
        float r = scene.obstacleRadius + f->h;
        if (scene.showPressure)
            c = BLACK;
        else
            c = (Color){0xDD, 0xDD, 0xDD, 255};
        // c.beginPath();
        // c.arc(
        //  cX(scene.obstacleX), cY(scene.obstacleY), cScale * r, 0.0, 2.0 * Math.PI);
        // c.closePath();
        // c.fill();

        float lineWidth = 20.0;
        c = (Color){0x40, 0x40, 0x40, 255};
        DrawRing((Vector2){cX(scene.obstacleX), cY(scene.obstacleY)},
                cScale*r - lineWidth/2.0,
                cScale*r + lineWidth/2.0,
                0.0,
                360.0,
                50,
                c);
        lineWidth = 1.0;
    }

    if (scene.showSwirls) {
        for (int i = 0; i < f->numSwirls; i++) {
            float x = f->swirlX[i];
            float y = f->swirlY[i];
            float r = scene.swirlMaxRadius;

            c = (Color){0x30, 0x30, 0x30, 255};
            DrawCircleLines(cX(x), cY(y), cScale * r, c);
        }
    }
}

void setObstacle(float x, float y, bool reset)
{
    float vx = 0.0;
    float vy = 0.0;

    if (!reset) {
        vx = (x - scene.obstacleX) / scene.dt;
        vy = (y - scene.obstacleY) / scene.dt;
    }

    scene.obstacleX = x;
    scene.obstacleY = y;
    scene.showObstacle = false;

    float r = scene.obstacleRadius;
    Fluid *f = &scene.fluid;
    int n = f->numY;
//    float cd = sqrtf(2) * f->h;

    for (int i = 1; i < f->numX-2; i++) {
        for (int j = 1; j < f->numY-2; j++) {
            f->s[i*n + j] = 1.0;

            float dx = (i + 0.5) * f->h - x;
            float dy = (j + 0.5) * f->h - y;

            float d2 = dx * dx + dy * dy;
            if (d2 < r * r) {
                // f->s[i*n + j] = 0.0;
                f->u[i*n + j] += 0.2 * vx;
                f->u[(i+1)*n + j] += 0.2 * vx;
                f->v[i*n + j] += 0.2 * vy;
                f->v[i*n + j+1] += 0.2 * vy;
            }
        }
    }
}

// interaction -------------------------------------------------------

static bool mouseDown = false;

void startDrag(float x, float y)
{
    mouseDown = true;

    x = x / cScale;
    y = (windowHeight - y) / cScale;

    setObstacle(x, y, true);
}

void drag(float x, float y)
{
    if (mouseDown) {
        x = x / cScale;
        y = (windowHeight - y) / cScale;
        setObstacle(x, y, false);
    }
}

void endDrag()
{
    mouseDown = false;
}

// main -------------------------------------------------------

void simulate()
{
    if (!scene.paused) {
        Fluid_simulate(&scene.fluid, scene.dt, scene.gravity, scene.numIters);
        scene.frameNr++;
    }
}

void update(Image img, Texture2D tex, Color *px)
{
    simulate();
    draw(img, tex, px);
}

int main()
{
    setupScene();

    InitWindow(windowWidth, windowHeight, "Fire Simulation");
    SetTargetFPS(30);

    Image img = {0};
    img.format = PIXELFORMAT_UNCOMPRESSED_R8G8B8A8;
    img.width = windowWidth;
    img.height = windowHeight;
    img.mipmaps = 1;
    img.data = MemAlloc(windowWidth*windowHeight*sizeof(Color));

    Texture tex = LoadTextureFromImage(img);

    Color *px = (Color*)MemAlloc(windowWidth*windowHeight*sizeof(Color));

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_F1)) 
            scene.showObstacle = !scene.showObstacle;
        if (IsKeyPressed(KEY_F2)) 
            scene.showPressure = !scene.showPressure;

        if (IsKeyPressed(KEY_P))
            scene.paused = !scene.paused;
        if (IsKeyPressed(KEY_M)) {
            scene.paused = false;
            simulate();
            scene.paused = true;
        }

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            Vector2 mp = GetMousePosition();
            startDrag(mp.x, mp.y);
        }
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) {
            endDrag();
        }
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            Vector2 mp = GetMousePosition();
            drag(mp.x, mp.y);
        }

        BeginDrawing(); {
            update(img, tex, px);

            GuiSlider((Rectangle){25, 10, 200, 20}, NULL, TextFormat("swirl probability: %.3f", scene.swirlProbability), &scene.swirlProbability, 0, 100); // Slider control
            GuiCheckBox((Rectangle){25, 45, 20, 20}, "burning floor", &scene.burningFloor);                          // Check Box control, returns true when active
            GuiCheckBox((Rectangle){25, 75, 20, 20}, "burning obstacle", &scene.burningObstacle);                          // Check Box control, returns true when active
            GuiCheckBox((Rectangle){25, 105, 20, 20}, "show swirls", &scene.showSwirls);                          // Check Box control, returns true when active
        } EndDrawing();
    }

    MemFree(px);
    UnloadTexture(tex);
    MemFree(img.data);

    Fluid_dispose(&scene.fluid);

    CloseWindow();

    return 0;
}
