#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#define GL_SILENCE_DEPRECATION
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/gl.h>

//<title>FLIP Fluid</title>
//<input type = "checkbox" checked  onclick = "scene.showParticles = !scene.showParticles">Particles &nbsp;
//<input type = "checkbox" onclick = "scene.showGrid = !scene.showGrid">Grid &nbsp;
//<input type = "checkbox" checked onclick = "scene.compensateDrift = !scene.compensateDrift">Compensate Drift &nbsp;
//<input type = "checkbox" checked onclick = "scene.separateParticles = !scene.separateParticles">Separate Particles &nbsp;
//PIC
//<input type = "range" min = "0" max = "10" value = "9" class = "slider" onchange="scene.flipRatio = 0.1 * this.value"> FLIP
//<br>
//<canvas id="myCanvas" style="border:2px solid"></canvas>


static int canvas_width;
static int canvas_height;

static float simHeight;
static float cScale;
static float simWidth;

#define U_FIELD 0
#define V_FIELD 1

#define FLUID_CELL 0
#define AIR_CELL   1
#define SOLID_CELL 2

static inline float clamp(float x, float min, float max)
{
    if (x < min)
        return min;
    else if (x > max)
        return max;
    else
        return x;
}

// ----------------- start of simulator ------------------------------

typedef struct
{
    float density;
    int fNumX;
    int fNumY;
    float h;
    float fInvSpacing;
    int fNumCells;

    float *u;
    float *v;
    float *du;
    float *dv;
    float *prevU;
    float *prevV;
    float *p;
    float *s;
    int32_t *cellType;
    float *cellColor;

    // particles

    int maxParticles;

    float *particlePos;
    float *particleColor;

    float *particleVel;
    float *particleDensity;
    float particleRestDensity;

    float particleRadius;
    float pInvSpacing;
    int pNumX;
    int pNumY;
    int pNumCells;

    int32_t *numCellParticles;
    int32_t *firstCellParticle;
    int32_t *cellParticleIds;

    int numParticles;
} FlipFluid;
typedef struct
{
    float gravity;
    float dt;
    float flipRatio;
    int numPressureIters;
    int numParticleIters;
    int frameNr;
    float overRelaxation;
    bool compensateDrift;
    bool separateParticles;
    float obstacleX;
    float obstacleY;
    float obstacleRadius;
    bool paused;
    bool showObstacle;
    float obstacleVelX;
    float obstacleVelY;
    bool showParticles;
    bool showGrid;
    FlipFluid fluid;
} Scene;
static Scene scene =
{
    gravity : -9.81,
//    gravity : 0.0,
    dt : 1.0 / 120.0,
    flipRatio : 0.9,
    numPressureIters : 100,
    numParticleIters : 2,
    frameNr : 0,
    overRelaxation : 1.9,
    compensateDrift : true,
    separateParticles : true,
    obstacleX : 0.0,
    obstacleY : 0.0,
    obstacleRadius: 0.15,
    paused: true,
    showObstacle: true,
    obstacleVelX: 0.0,
    obstacleVelY: 0.0,
    showParticles: true,
    showGrid: false,
    fluid: {0}
};

FlipFluid FlipFluid_ctor(float density, float width, float height, float spacing, float particleRadius, int maxParticles)
{
    FlipFluid this = {0};

    // fluid

    this.density = density;
    this.fNumX = floorf(width / spacing) + 1;
    this.fNumY = floorf(height / spacing) + 1;
    this.h = fmaxf(width / this.fNumX, height / this.fNumY);
    this.fInvSpacing = 1.0 / this.h;
    this.fNumCells = this.fNumX * this.fNumY;

    this.u = calloc(this.fNumCells, sizeof(float));
    this.v = calloc(this.fNumCells, sizeof(float));
    this.du = calloc(this.fNumCells, sizeof(float));
    this.dv = calloc(this.fNumCells, sizeof(float));
    this.prevU = calloc(this.fNumCells, sizeof(float));
    this.prevV = calloc(this.fNumCells, sizeof(float));
    this.p = calloc(this.fNumCells, sizeof(float));
    this.s = calloc(this.fNumCells, sizeof(float));
    this.cellType = calloc(this.fNumCells, sizeof(int32_t));
    this.cellColor = calloc(3 * this.fNumCells, sizeof(float));

    // particles

    this.maxParticles = maxParticles;

    this.particlePos = calloc(2 * this.maxParticles, sizeof(float));
    this.particleColor = calloc(3 * this.maxParticles, sizeof(float));
    for (int i = 0; i < this.maxParticles; i++)
        this.particleColor[3 * i + 2] = 1.0; // blue particles

    this.particleVel = calloc(2 * this.maxParticles, sizeof(float));
    this.particleDensity = calloc(this.fNumCells, sizeof(float));
    this.particleRestDensity = 0.0;

    this.particleRadius = particleRadius;
    this.pInvSpacing = 1.0 / (2.2 * particleRadius);
    this.pNumX = floorf(width * this.pInvSpacing) + 1;
    this.pNumY = floorf(height * this.pInvSpacing) + 1;
    this.pNumCells = this.pNumX * this.pNumY;

    this.numCellParticles = calloc(this.pNumCells, sizeof(int32_t));
    this.firstCellParticle = calloc(this.pNumCells + 1, sizeof(int32_t));
    this.cellParticleIds = calloc(this.maxParticles, sizeof(int32_t));

    this.numParticles = 0;

    return this;
}

void FlipFluid_dtor(FlipFluid *this)
{
    // particles

    free(this->particlePos);
    free(this->particleColor);

    free(this->particleVel);
    free(this->particleDensity);

    free(this->numCellParticles);
    free(this->firstCellParticle);
    free(this->cellParticleIds);

    // fluid

    free(this->u);
    free(this->v);
    free(this->du);
    free(this->dv);
    free(this->prevU);
    free(this->prevV);
    free(this->p);
    free(this->s);
    free(this->cellType);
    free(this->cellColor);
}

void FlipFluid_integrateParticles(FlipFluid *this, float dt, float gravity)
{
    for (int i = 0; i < this->numParticles; i++) {
        this->particleVel[2 * i + 1] += dt * gravity;
        this->particlePos[2 * i] += this->particleVel[2 * i] * dt;
        this->particlePos[2 * i + 1] += this->particleVel[2 * i + 1] * dt;
    }
}

void FlipFluid_pushParticlesApart(FlipFluid *this, int numIters)
{
    float colorDiffusionCoeff = 0.001;

    // count particles per cell

    memset(this->numCellParticles, 0, sizeof(int32_t) * this->pNumCells);;

    for (int i = 0; i < this->numParticles; i++) {
        float x = this->particlePos[2 * i];
        float y = this->particlePos[2 * i + 1];

        int xi = clamp(floorf(x * this->pInvSpacing), 0, this->pNumX - 1);
        int yi = clamp(floorf(y * this->pInvSpacing), 0, this->pNumY - 1);
        int cellNr = xi * this->pNumY + yi;
        this->numCellParticles[cellNr]++;
    }

    // partial sums

    int first = 0;

    for (int i = 0; i < this->pNumCells; i++) {
        first += this->numCellParticles[i];
        this->firstCellParticle[i] = first;
    }
    this->firstCellParticle[this->pNumCells] = first;     // guard

    // fill particles into cells

    for (int i = 0; i < this->numParticles; i++) {
        float x = this->particlePos[2 * i];
        float y = this->particlePos[2 * i + 1];

        int xi = clamp(floorf(x * this->pInvSpacing), 0, this->pNumX - 1);
        int yi = clamp(floorf(y * this->pInvSpacing), 0, this->pNumY - 1);
        int cellNr = xi * this->pNumY + yi;
        this->firstCellParticle[cellNr]--;
        this->cellParticleIds[this->firstCellParticle[cellNr]] = i;
    }

    // push particles apart

    float minDist = 2.0 * this->particleRadius;
    float minDist2 = minDist * minDist;

    for (int iter = 0; iter < numIters; iter++) {
        for (int i = 0; i < this->numParticles; i++) {
            float px = this->particlePos[2 * i];
            float py = this->particlePos[2 * i + 1];

            int pxi = floorf(px * this->pInvSpacing);
            int pyi = floorf(py * this->pInvSpacing);
            float x0 = fmaxf(pxi - 1, 0);
            float y0 = fmaxf(pyi - 1, 0);
            float x1 = fminf(pxi + 1, this->pNumX - 1);
            float y1 = fminf(pyi + 1, this->pNumY - 1);

            for (int xi = x0; xi <= x1; xi++) {
                for (int yi = y0; yi <= y1; yi++) {
                    int cellNr = xi * this->pNumY + yi;
                    int first = this->firstCellParticle[cellNr];
                    int last = this->firstCellParticle[cellNr + 1];
                    for (int j = first; j < last; j++) {
                        int id = this->cellParticleIds[j];
                        if (id == i)
                            continue;
                        float qx = this->particlePos[2 * id];
                        float qy = this->particlePos[2 * id + 1];

                        float dx = qx - px;
                        float dy = qy - py;
                        float d2 = dx * dx + dy * dy;
                        if (d2 > minDist2 || d2 == 0.0)
                            continue;
                        float d = sqrtf(d2);
                        float s = 0.5 * (minDist - d) / d;
                        dx *= s;
                        dy *= s;
                        this->particlePos[2 * i] -= dx;
                        this->particlePos[2 * i + 1] -= dy;
                        this->particlePos[2 * id] += dx;
                        this->particlePos[2 * id + 1] += dy;

                        // diffuse colors

                        for (int k = 0; k < 3; k++) {
                            float color0 = this->particleColor[3 * i + k];
                            float color1 = this->particleColor[3 * id + k];
                            float color = (color0 + color1) * 0.5;
                            this->particleColor[3 * i + k] = color0 + (color - color0) * colorDiffusionCoeff;
                            this->particleColor[3 * id + k] = color1 + (color - color1) * colorDiffusionCoeff;
                        }
                    }
                }
            }
        }
    }
}

void FlipFluid_handleParticleCollisions(FlipFluid *this, float obstacleX, float obstacleY, float obstacleRadius)
{
    float h = 1.0 / this->fInvSpacing;
    float r = this->particleRadius;
//    float or = obstacleRadius;
//    float or2 = or * or;
    float minDist = obstacleRadius + r;
    float minDist2 = minDist * minDist;

    float minX = h + r;
    float maxX = (this->fNumX - 1) * h - r;
    float minY = h + r;
    float maxY = (this->fNumY - 1) * h - r;

    for (int i = 0; i < this->numParticles; i++) {
        float x = this->particlePos[2 * i];
        float y = this->particlePos[2 * i + 1];

        float dx = x - obstacleX;
        float dy = y - obstacleY;
        float d2 = dx * dx + dy * dy;

        // obstacle collision

        if (d2 < minDist2) {
            //float d = sqrtf(d2);
            //float s = (minDist - d) / d;
            //x += dx * s;
            //y += dy * s;

            this->particleVel[2*i + 0] += scene.obstacleVelX;
            this->particleVel[2*i + 1] += scene.obstacleVelY;
        }

        // wall collisions

        if (x < minX) {
            x = minX;
            this->particleVel[2 * i] = 0.0;
        }
        if (x > maxX) {
            x = maxX;
            this->particleVel[2 * i] = 0.0;
        }
        if (y < minY) {
            y = minY;
            this->particleVel[2 * i + 1] = 0.0;
        }
        if (y > maxY) {
            y = maxY;
            this->particleVel[2 * i + 1] = 0.0;
        }
        this->particlePos[2 * i] = x;
        this->particlePos[2 * i + 1] = y;
    }
}

void FlipFluid_updateParticleDensity(FlipFluid *this)
{
    int n = this->fNumY;
    float h = this->h;
    float h1 = this->fInvSpacing;
    float h2 = 0.5 * h;

    FlipFluid *f = this;
    float *d = f->particleDensity;

    memset(d, 0, sizeof(float) * this->fNumCells);

    for (int i = 0; i < this->numParticles; i++) {
        float x = this->particlePos[2 * i];
        float y = this->particlePos[2 * i + 1];

        x = clamp(x, h, (this->fNumX - 1) * h);
        y = clamp(y, h, (this->fNumY - 1) * h);

        int x0 = floorf((x - h2) * h1);
        float tx = ((x - h2) - x0 * h) * h1;
        int x1 = fminf(x0 + 1, this->fNumX-2);

        int y0 = floorf((y-h2)*h1);
        float ty = ((y - h2) - y0*h) * h1;
        int y1 = fminf(y0 + 1, this->fNumY-2);

        float sx = 1.0 - tx;
        float sy = 1.0 - ty;

        if (x0 < this->fNumX && y0 < this->fNumY) d[x0 * n + y0] += sx * sy;
        if (x1 < this->fNumX && y0 < this->fNumY) d[x1 * n + y0] += tx * sy;
        if (x1 < this->fNumX && y1 < this->fNumY) d[x1 * n + y1] += tx * ty;
        if (x0 < this->fNumX && y1 < this->fNumY) d[x0 * n + y1] += sx * ty;
    }

    if (this->particleRestDensity == 0.0) {
        float sum = 0.0;
        int numFluidCells = 0;

        for (int i = 0; i < this->fNumCells; i++) {
            if (this->cellType[i] == FLUID_CELL) {
                sum += d[i];
                numFluidCells++;
            }
        }

        if (numFluidCells > 0)
            this->particleRestDensity = sum / numFluidCells;
    }

//    for (int xi = 1; xi < this->fNumX; xi++) {
//        for (int yi = 1; yi < this->fNumY; yi++) {
//            int cellNr = xi * n + yi;
//            if (this->cellType[cellNr] != FLUID_CELL)
//                continue;
//            float hx = this->h;
//            float hy = this->h;
//
//            if (this->cellType[(xi - 1) * n + yi] == SOLID_CELL || this->cellType[(xi + 1) * n + yi] == SOLID_CELL)
//                hx -= this->particleRadius;
//            if (this->cellType[xi * n + yi - 1] == SOLID_CELL || this->cellType[xi * n + yi + 1] == SOLID_CELL)
//                hy -= this->particleRadius;
//
//            float scale = this->h * this->h / (hx * hy);
//            d[cellNr] *= scale;
//        }
//    }
}

void FlipFluid_transferVelocities(FlipFluid *this, bool toGrid, float flipRatio)
{
    int n = this->fNumY;
    float h = this->h;
    float h1 = this->fInvSpacing;
    float h2 = 0.5 * h;

    if (toGrid) {
        memcpy(this->prevU, this->u, sizeof(float) * this->fNumCells);
        memcpy(this->prevV, this->v, sizeof(float) * this->fNumCells);

        memset(this->du, 0, sizeof(float) * this->fNumCells);
        memset(this->dv, 0, sizeof(float) * this->fNumCells);
        memset(this->u, 0, sizeof(float) * this->fNumCells);
        memset(this->v, 0, sizeof(float) * this->fNumCells);

        for (int i = 0; i < this->fNumCells; i++)
            this->cellType[i] = this->s[i] == 0.0 ? SOLID_CELL : AIR_CELL;

        for (int i = 0; i < this->numParticles; i++) {
            float x = this->particlePos[2 * i];
            float y = this->particlePos[2 * i + 1];
            int xi = clamp(floorf(x * h1), 0, this->fNumX - 1);
            int yi = clamp(floorf(y * h1), 0, this->fNumY - 1);
            int cellNr = xi * n + yi;
            if (this->cellType[cellNr] == AIR_CELL)
                this->cellType[cellNr] = FLUID_CELL;
        }
    }

    for (int component = 0; component < 2; component++) {
        float dx = component == 0 ? 0.0 : h2;
        float dy = component == 0 ? h2 : 0.0;

        float *f = component == 0 ? this->u : this->v;
        float *prevF = component == 0 ? this->prevU : this->prevV;
        float *d = component == 0 ? this->du : this->dv;

        for (int i = 0; i < this->numParticles; i++) {
            float x = this->particlePos[2 * i];
            float y = this->particlePos[2 * i + 1];

            x = clamp(x, h, (this->fNumX - 1) * h);
            y = clamp(y, h, (this->fNumY - 1) * h);

            int x0 = fminf(floorf((x - dx) * h1), this->fNumX - 2);
            float tx = ((x - dx) - x0 * h) * h1;
            int x1 = fminf(x0 + 1, this->fNumX-2);

            int y0 = fminf(floorf((y-dy)*h1), this->fNumY-2);
            float ty = ((y - dy) - y0*h) * h1;
            int y1 = fminf(y0 + 1, this->fNumY-2);

            float sx = 1.0 - tx;
            float sy = 1.0 - ty;

            float d0 = sx*sy;
            float d1 = tx*sy;
            float d2 = tx*ty;
            float d3 = sx*ty;

            int nr0 = x0*n + y0;
            int nr1 = x1*n + y0;
            int nr2 = x1*n + y1;
            int nr3 = x0*n + y1;

            if (toGrid) {
                float pv = this->particleVel[2 * i + component];
                f[nr0] += pv * d0;  d[nr0] += d0;
                f[nr1] += pv * d1;  d[nr1] += d1;
                f[nr2] += pv * d2;  d[nr2] += d2;
                f[nr3] += pv * d3;  d[nr3] += d3;
            }
            else {
                int offset = component == 0 ? n : 1;
                float valid0 = this->cellType[nr0] == FLUID_CELL || this->cellType[nr0 - offset] == FLUID_CELL ? 1.0 : 0.0;
                float valid1 = this->cellType[nr1] == FLUID_CELL || this->cellType[nr1 - offset] == FLUID_CELL ? 1.0 : 0.0;
                float valid2 = this->cellType[nr2] == FLUID_CELL || this->cellType[nr2 - offset] == FLUID_CELL ? 1.0 : 0.0;
                float valid3 = this->cellType[nr3] == FLUID_CELL || this->cellType[nr3 - offset] == FLUID_CELL ? 1.0 : 0.0;

                float v = this->particleVel[2 * i + component];
                float d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                if (d > 0.0) {
                    float picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d;
                    float corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) + valid1 * d1 * (f[nr1] - prevF[nr1])
                        + valid2 * d2 * (f[nr2] - prevF[nr2]) + valid3 * d3 * (f[nr3] - prevF[nr3])) / d;
                    float flipV = v + corr;

                    this->particleVel[2 * i + component] = (1.0 - flipRatio) * picV + flipRatio * flipV;
                }
            }
        }

        if (toGrid) {
            for (int i = 0; i < this->fNumCells; i++) {
                if (d[i] > 0.0)
                    f[i] /= d[i];
            }

            // restore solid cells

            for (int i = 0; i < this->fNumX; i++) {
                for (int j = 0; j < this->fNumY; j++) {
                    int solid = this->cellType[i * n + j] == SOLID_CELL;
                    if (solid || (i > 0 && this->cellType[(i - 1) * n + j] == SOLID_CELL))
                        this->u[i * n + j] = this->prevU[i * n + j];
                    if (solid || (j > 0 && this->cellType[i * n + j - 1] == SOLID_CELL))
                        this->v[i * n + j] = this->prevV[i * n + j];
                }
            }
        }
    }
}

void FlipFluid_solveIncompressibility(FlipFluid *this, int numIters, float dt, float overRelaxation, bool compensateDrift /*= true*/)
{
    memset(this->p, 0, sizeof(float) * this->fNumCells);
    memcpy(this->prevU, this->u, sizeof(float) * this->fNumCells);
    memcpy(this->prevV, this->v, sizeof(float) * this->fNumCells);

    int n = this->fNumY;
    float cp = this->density * this->h / dt;

//    for (int i = 0; i < this->fNumCells; i++) {
//        float u = this->u[i];
//        float v = this->v[i];
//    }

    for (int iter = 0; iter < numIters; iter++) {
        for (int i = 1; i < this->fNumX-1; i++) {
            for (int j = 1; j < this->fNumY-1; j++) {
                if (this->cellType[i*n + j] != FLUID_CELL)
                    continue;

                int center = i * n + j;
                int left = (i - 1) * n + j;
                int right = (i + 1) * n + j;
                int bottom = i * n + j - 1;
                int top = i * n + j + 1;

                float sx0 = this->s[left];
                float sx1 = this->s[right];
                float sy0 = this->s[bottom];
                float sy1 = this->s[top];
                float s = sx0 + sx1 + sy0 + sy1;
                if (s == 0.0)
                    continue;

                float div = this->u[right] - this->u[center] +
                        this->v[top] - this->v[center];

                if (this->particleRestDensity > 0.0 && compensateDrift) {
                    float k = 1.0;
                    float compression = this->particleDensity[i*n + j] - this->particleRestDensity;
                    if (compression > 0.0)
                        div = div - k * compression;
                }

                float p = -div / s;
                p *= overRelaxation;
                this->p[center] += cp * p;

                this->u[center] -= sx0 * p;
                this->u[right] += sx1 * p;
                this->v[center] -= sy0 * p;
                this->v[top] += sy1 * p;
            }
        }
    }
}

void FlipFluid_updateParticleColors(FlipFluid *this)
{
    // for (int i = 0; i < this->numParticles; i++) {
    //  this->particleColor[3 * i] *= 0.99;
    //  this->particleColor[3 * i + 1] *= 0.99;
    //  this->particleColor[3 * i + 2] =
    //      clamp(this->particleColor[3 * i + 2] + 0.001, 0.0, 1.0);
    // }

    // return;

    float h1 = this->fInvSpacing;

    for (int i = 0; i < this->numParticles; i++) {
        float s = 0.01;

        this->particleColor[3 * i] = clamp(this->particleColor[3 * i] - s, 0.0, 1.0);
        this->particleColor[3 * i + 1] = clamp(this->particleColor[3 * i + 1] - s, 0.0, 1.0);
        this->particleColor[3 * i + 2] = clamp(this->particleColor[3 * i + 2] + s, 0.0, 1.0);

        float x = this->particlePos[2 * i];
        float y = this->particlePos[2 * i + 1];
        int xi = clamp(floorf(x * h1), 1, this->fNumX - 1);
        int yi = clamp(floorf(y * h1), 1, this->fNumY - 1);
        int cellNr = xi * this->fNumY + yi;

        float d0 = this->particleRestDensity;

        if (d0 > 0.0) {
            float relDensity = this->particleDensity[cellNr] / d0;
            if (relDensity < 0.7) {
                float s = 0.8;
                this->particleColor[3 * i] = s;
                this->particleColor[3 * i + 1] = s;
                this->particleColor[3 * i + 2] = 1.0;
            }
        }
    }
}

void FlipFluid_setSciColor(FlipFluid *this, int cellNr, float val, float minVal, float maxVal)
{
    val = fminf(fmaxf(val, minVal), maxVal- 0.0001);
    float d = maxVal - minVal;
    val = d == 0.0 ? 0.5 : (val - minVal) / d;
    float m = 0.25;
    int num = floorf(val / m);
    float s = (val - num * m) / m;
    float r = 0.0;
    float g = 0.0;
    float b = 0.0;

    switch (num) {
        case 0 : r = 0.0; g = s; b = 1.0; break;
        case 1 : r = 0.0; g = 1.0; b = 1.0-s; break;
        case 2 : r = s; g = 1.0; b = 0.0; break;
        case 3 : r = 1.0; g = 1.0 - s; b = 0.0; break;
    }

    this->cellColor[3 * cellNr] = r;
    this->cellColor[3 * cellNr + 1] = g;
    this->cellColor[3 * cellNr + 2] = b;
}

void FlipFluid_updateCellColors(FlipFluid *this)
{
    for (int i = 0; i < this->fNumCells; ++i) {
        this->cellColor[3 * i] = 0.0;
        this->cellColor[3 * i + 1] = 0.0;
        this->cellColor[3 * i + 2] = 0.0;
    }

    for (int i = 0; i < this->fNumCells; i++) {
        if (this->cellType[i] == SOLID_CELL) {
            this->cellColor[3*i] = 0.5;
            this->cellColor[3*i + 1] = 0.5;
            this->cellColor[3*i + 2] = 0.5;
        }
        else if (this->cellType[i] == FLUID_CELL) {
            float d = this->particleDensity[i];
            if (this->particleRestDensity > 0.0)
                d /= this->particleRestDensity;
            FlipFluid_setSciColor(this, i, d, 0.0, 2.0);
        }
    }
}

void FlipFluid_simulate(FlipFluid *this, Scene *s)
{
    int numSubSteps = 1;
    float sdt = s->dt / numSubSteps;

    for (int step = 0; step < numSubSteps; step++) {
        FlipFluid_integrateParticles(this, sdt, s->gravity);
        if (s->separateParticles)
            FlipFluid_pushParticlesApart(this, s->numParticleIters);
        FlipFluid_handleParticleCollisions(this, s->obstacleX, s->obstacleY, s->obstacleRadius);
        FlipFluid_transferVelocities(this, true, 0.0f);
        FlipFluid_updateParticleDensity(this);
        FlipFluid_solveIncompressibility(this, s->numPressureIters, sdt, s->overRelaxation, s->compensateDrift);
        FlipFluid_transferVelocities(this, false, s->flipRatio);
    }

    FlipFluid_updateParticleColors(this);
    FlipFluid_updateCellColors(this);
}

// ----------------- end of simulator ------------------------------

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
    float r = scene.obstacleRadius;
    FlipFluid *f = &scene.fluid;
    int n = f->fNumY;
//    float cd = sqrtf(2) * f->h;

    for (int i = 1; i < f->fNumX-2; i++) {
        for (int j = 1; j < f->fNumY-2; j++) {
            f->s[i*n + j] = 1.0;

            float dx = (i + 0.5) * f->h - x;
            float dy = (j + 0.5) * f->h - y;

            if (dx * dx + dy * dy < r * r) {
                f->s[i*n + j] = 0.0;
                f->u[i*n + j] = vx;
                f->u[(i+1)*n + j] = vx;
                f->v[i*n + j] = vy;
                f->v[i*n + j+1] = vy;
            }
        }
    }

    scene.showObstacle = true;
    scene.obstacleVelX = vx;
    scene.obstacleVelY = vy;
}

void setupScene()
{
    scene.obstacleRadius = 0.15;
    scene.overRelaxation = 1.9;

    scene.dt = 1.0 / 60.0;
    scene.numPressureIters = 50;
    scene.numParticleIters = 2;

    int res = 100;

    float tankHeight = 1.0 * simHeight;
    float tankWidth = 1.0 * simWidth;
    float h = tankHeight / res;
    float density = 1000.0;

    float relWaterHeight = 0.8;
    float relWaterWidth = 0.6;

    // dam break

    // compute number of particles

    float r = 0.3 * h;    // particle radius w.r.t. cell size
    float dx = 2.0 * r;
    float dy = sqrtf(3.0) / 2.0 * dx;

    int numX = floorf((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
    int numY = floorf((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy);
    int maxParticles = numX * numY;

    // create fluid

    scene.fluid = FlipFluid_ctor(density, tankWidth, tankHeight, h, r, maxParticles);
    FlipFluid *f = &scene.fluid;

    // create particles

    f->numParticles = numX * numY;
    int p = 0;
    for (int i = 0; i < numX; i++) {
        for (int j = 0; j < numY; j++) {
            f->particlePos[p++] = h + r + dx * i + (j % 2 == 0 ? 0.0 : r);
            f->particlePos[p++] = h + r + dy * j;
        }
    }

    // setup grid cells for tank

    int n = f->fNumY;

    for (int i = 0; i < f->fNumX; i++) {
        for (int j = 0; j < f->fNumY; j++) {
            float s = 1.0;    // fluid
            if (i == 0 || i == f->fNumX-1 || j == 0 || j == f->fNumY-1)
                s = 0.0;    // solid
            f->s[i*n + j] = s;
        }
    }

    setObstacle(3.0, 2.0, true);
}


// draw -------------------------------------------------------

static const char *pointVertexShaderSrc =
        "#version 330 core\n"
        "layout(location = 0) in vec2 attrPosition;\n"
        "layout(location = 1) in vec3 attrColor;\n"
        "uniform vec2 domainSize;\n"
        "uniform float pointSize;\n"
        "uniform float drawDisk;\n"
        "\n"
        "out vec3 fragColor;\n"
        "out float fragDrawDisk;\n"
        "\n"
        "void main() {\n"
        "    vec4 screenTransform = vec4(2.0 / domainSize.x, 2.0 / domainSize.y, -1.0, -1.0);\n"
        "    gl_Position = vec4(attrPosition * screenTransform.xy + screenTransform.zw, 0.0, 1.0);\n"
        "\n"
        "    gl_PointSize = pointSize;\n"
        "    fragColor = attrColor;\n"
        "    fragDrawDisk = drawDisk;\n"
        "}\n";

static const char *pointFragmentShaderSrc =
        "#version 330 core\n"
        "in vec3 fragColor;\n"
        "in float fragDrawDisk;\n"
        "out vec4 outColor;\n"
        "\n"
        "void main() {\n"
        "    if (fragDrawDisk == 1.0) {\n"
        "        float rx = 0.5 - gl_PointCoord.x;\n"
        "        float ry = 0.5 - gl_PointCoord.y;\n"
        "        float r2 = rx * rx + ry * ry;\n"
        "        if (r2 > 0.25)\n"
        "            discard;\n"
        "    }\n"
        "    outColor = vec4(fragColor, 1.0);\n"
        "}\n";

static const char *meshVertexShaderSrc =
        "#version 330 core\n"
        "layout(location = 0) in vec2 attrPosition;\n"
        "uniform vec2 domainSize;\n"
        "uniform vec3 color;\n"
        "uniform vec2 translation;\n"
        "uniform float scale;\n"
        "\n"
        "out vec3 fragColor;\n"
        "\n"
        "void main() {\n"
        "    vec2 v = translation + attrPosition * scale;\n"
        "    vec4 screenTransform = vec4(2.0 / domainSize.x, 2.0 / domainSize.y, -1.0, -1.0);\n"
        "    gl_Position = vec4(v * screenTransform.xy + screenTransform.zw, 0.0, 1.0);\n"
        "    fragColor = color;\n"
        "}\n";

static const char *meshFragmentShaderSrc =
        "#version 330 core\n"
        "in vec3 fragColor;\n"
        "out vec4 outColor;\n"
        "\n"
        "void main() {\n"
        "    outColor = vec4(fragColor, 1.0);\n"
        "}\n";

// --- Helper compile/link ---
static GLuint compile_shader(GLenum type, const char *src)
{
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, NULL);
    glCompileShader(s);
    GLint ok = 0;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        GLint len = 0;
        glGetShaderiv(s, GL_INFO_LOG_LENGTH, &len);
        char *log = malloc(len + 1);
        glGetShaderInfoLog(s, len, NULL, log);
        fprintf(stderr, "Shader compile error: %s\n", log);
        free(log);
    }
    return s;
}

static GLuint create_program(const char *vs_src, const char *fs_src)
{
    GLuint vs = compile_shader(GL_VERTEX_SHADER, vs_src);
    GLuint fs = compile_shader(GL_FRAGMENT_SHADER, fs_src);
    GLuint prog = glCreateProgram();
    glAttachShader(prog, vs);
    glAttachShader(prog, fs);
    glLinkProgram(prog);
    GLint ok = 0;
    glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if (!ok) {
        GLint len = 0;
        glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &len);
        char *log = malloc(len + 1);
        glGetProgramInfoLog(prog, len, NULL, log);
        fprintf(stderr, "Program link error: %s\n", log);
        free(log);
    }
    glDeleteShader(vs);
    glDeleteShader(fs);
    return prog;
}

// --- GL resources (global) ---
static GLuint pointProgram = 0;
static GLuint meshProgram = 0;

static GLuint pointVertexBuffer = 0;
static GLuint pointColorBuffer = 0;

static GLuint gridVertBuffer = 0;
static GLuint gridColorBuffer = 0;

static GLuint diskVertBuffer = 0;
static GLuint diskIdBuffer = 0;

// --- draw() function converted ---
static float *cellCenters = NULL;
static float *diskVerts = NULL;
static unsigned short *diskIds = NULL;
void draw_scene()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glViewport(0, 0, canvas_width, canvas_height);

    if (pointProgram == 0)
        pointProgram = create_program(pointVertexShaderSrc, pointFragmentShaderSrc);
    if (meshProgram == 0)
        meshProgram = create_program(meshVertexShaderSrc, meshFragmentShaderSrc);

    // GRID
    if (gridVertBuffer == 0) {
        int fNumX = scene.fluid.fNumX;
        int fNumY = scene.fluid.fNumY;
        int fNumCells = scene.fluid.fNumCells;
        float h = scene.fluid.h;

        glGenBuffers(1, &gridVertBuffer);
        if (cellCenters == NULL)
            cellCenters = malloc(sizeof(float) * 2 * fNumCells);
        int p = 0;
        for (int i = 0; i < fNumX; ++i) {
            for (int j = 0; j < fNumY; ++j) {
                cellCenters[p++] = (i + 0.5f) * h;
                cellCenters[p++] = (j + 0.5f) * h;
            }
        }
        glBindBuffer(GL_ARRAY_BUFFER, gridVertBuffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 2 * fNumCells, cellCenters, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    if (cellCenters != NULL) {
        free(cellCenters);
        cellCenters = NULL;
    }

    if (gridColorBuffer == 0) {
        glGenBuffers(1, &gridColorBuffer);
    }

    if (scene.showGrid) {
        float cellPixel = scene.fluid.h / simWidth * (float)canvas_width;
        float pointSize = ceilf(0.98f * cellPixel); // leggermente più grande e arrotondato

        glUseProgram(pointProgram);
        GLint loc_domain = glGetUniformLocation(pointProgram, "domainSize");
        glUniform2f(loc_domain, simWidth, simHeight);
        GLint loc_pointSize = glGetUniformLocation(pointProgram, "pointSize");
        glUniform1f(loc_pointSize, pointSize);
        GLint loc_drawDisk = glGetUniformLocation(pointProgram, "drawDisk");
        glUniform1f(loc_drawDisk, 0.0f); // 1.0f per dischi pieni

        glBindBuffer(GL_ARRAY_BUFFER, gridVertBuffer);
        GLint posLoc = glGetAttribLocation(pointProgram, "attrPosition");
        glEnableVertexAttribArray((GLuint)posLoc);
        glVertexAttribPointer((GLuint)posLoc, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glBindBuffer(GL_ARRAY_BUFFER, gridColorBuffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * scene.fluid.fNumCells, scene.fluid.cellColor, GL_DYNAMIC_DRAW);
        GLint colorLoc = glGetAttribLocation(pointProgram, "attrColor");
        glEnableVertexAttribArray((GLuint)colorLoc);
        glVertexAttribPointer((GLuint)colorLoc, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glEnable(GL_PROGRAM_POINT_SIZE);
        glDrawArrays(GL_POINTS, 0, scene.fluid.fNumCells);

        glDisableVertexAttribArray((GLuint)posLoc);
        glDisableVertexAttribArray((GLuint)colorLoc);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    // PARTICLES
    if (scene.showParticles) {
        float pointSize = 2.0f * scene.fluid.particleRadius / simWidth * (float)canvas_width;

        glUseProgram(pointProgram);
        GLint loc_domain2 = glGetUniformLocation(pointProgram, "domainSize");
        glUniform2f(loc_domain2, simWidth, simHeight);
        GLint loc_pointSize2 = glGetUniformLocation(pointProgram, "pointSize");
        glUniform1f(loc_pointSize2, pointSize);
        GLint loc_drawDisk2 = glGetUniformLocation(pointProgram, "drawDisk");
        glUniform1f(loc_drawDisk2, 1.0f);

        if (pointVertexBuffer == 0) glGenBuffers(1, &pointVertexBuffer);
        if (pointColorBuffer == 0) glGenBuffers(1, &pointColorBuffer);

        glBindBuffer(GL_ARRAY_BUFFER, pointVertexBuffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 2 * scene.fluid.numParticles, scene.fluid.particlePos, GL_DYNAMIC_DRAW);
        GLint posLoc = glGetAttribLocation(pointProgram, "attrPosition");
        glEnableVertexAttribArray((GLuint)posLoc);
        glVertexAttribPointer((GLuint)posLoc, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glBindBuffer(GL_ARRAY_BUFFER, pointColorBuffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * scene.fluid.numParticles, scene.fluid.particleColor, GL_DYNAMIC_DRAW);
        GLint colorLoc = glGetAttribLocation(pointProgram, "attrColor");
        glEnableVertexAttribArray((GLuint)colorLoc);
        glVertexAttribPointer((GLuint)colorLoc, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glEnable(GL_PROGRAM_POINT_SIZE);
        glDrawArrays(GL_POINTS, 0, scene.fluid.numParticles);

        glDisableVertexAttribArray((GLuint)posLoc);
        glDisableVertexAttribArray((GLuint)colorLoc);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    // DISK
    const int numSegs = 50;
    if (diskVertBuffer == 0) {
        glGenBuffers(1, &diskVertBuffer);
        float dphi = 2.0f * (float)M_PI / (float)numSegs;
        if (diskVerts == NULL)
            diskVerts = malloc(sizeof(float) * (2 * numSegs + 2));
        int p = 0;
        diskVerts[p++] = 0.0f;
        diskVerts[p++] = 0.0f;
        for (int i = 0; i < numSegs; ++i) {
            diskVerts[p++] = cosf(i * dphi);
            diskVerts[p++] = sinf(i * dphi);
        }
        glBindBuffer(GL_ARRAY_BUFFER, diskVertBuffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * (2 * numSegs + 2), diskVerts, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        glGenBuffers(1, &diskIdBuffer);
        if (diskIds == NULL)
            diskIds = malloc(sizeof(unsigned short) * 3 * numSegs);
        p = 0;
        for (int i = 0; i < numSegs; ++i) {
            diskIds[p++] = 0;
            diskIds[p++] = 1 + i;
            diskIds[p++] = 1 + (i + 1) % numSegs;
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, diskIdBuffer);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short) * 3 * numSegs, diskIds, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }
    if (diskIds != NULL) {
        free(diskIds);
        diskIds = NULL;
    }
    if (diskVerts != NULL) {
        free(diskVerts);
        diskVerts = NULL;
    }

    glUseProgram(meshProgram);
    GLint loc_domain_mesh = glGetUniformLocation(meshProgram, "domainSize");
    glUniform2f(loc_domain_mesh, simWidth, simHeight);
    GLint loc_color_mesh = glGetUniformLocation(meshProgram, "color");
    glUniform3f(loc_color_mesh, 1.0f, 0.0f, 0.0f);
    GLint loc_trans = glGetUniformLocation(meshProgram, "translation");
    glUniform2f(loc_trans, scene.obstacleX, scene.obstacleY);
    GLint loc_scale = glGetUniformLocation(meshProgram, "scale");
    glUniform1f(loc_scale, scene.obstacleRadius + scene.fluid.particleRadius);

    GLint posLoc = glGetAttribLocation(meshProgram, "attrPosition");
    glEnableVertexAttribArray((GLuint)posLoc);
    glBindBuffer(GL_ARRAY_BUFFER, diskVertBuffer);
    glVertexAttribPointer((GLuint)posLoc, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, diskIdBuffer);
    glDrawElements(GL_TRIANGLES, 3 * numSegs, GL_UNSIGNED_SHORT, (void*)0);

    glDisableVertexAttribArray((GLuint)posLoc);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

// interaction state
static int mouseDown = 0;
static double lastMouseX = 0.0;
static double lastMouseY = 0.0;

// x_window,y_window in pixel (framebuffer coords), started = 1 se startDrag
static void setObstacleFromWindow(double wx, double wy, int started, GLFWwindow *win)
{
    int fbw, fbh;
    glfwGetFramebufferSize(win, &fbw, &fbh); // framebuffer size in pixel
    // canvas clientHeight equivalente: fbh
    // converti: mx = wx; my = wy (eventi GLFW forniscono già client coords)
    // nota: GLFW cursor pos y origin is top; converti a domain y like JS: y = (canvas.height - my)/cScale
    // qui domain va da 0..simHeight, e abbiamo usato attrPosition in dominio 0..1
    double mx = wx;
    double my = wy;
    // cScale in JS mappava coordinate canvas pixel -> domain units. qui assumiamo domain mapped to framebuffer width/height:
    // x_domain = mx / fbw * simWidth; y_domain = (fbh - my) / fbh * simHeight
    double x = (mx / (double)fbw) * simWidth;
    double y = ((fbh - my) / (double)fbh) * simHeight;

    setObstacle(x, y, started);

    scene.obstacleX = (float)x;
    scene.obstacleY = (float)y;

    if (started) {
        scene.paused = 0; // uncomment if you have scene.paused
    }
}

void startDrag(float x, float y)
{
//    let bounds = canvas.getBoundingClientRect();
//
//    let mx = x - bounds.left - canvas.clientLeft;
//    let my = y - bounds.top - canvas.clientTop;
//    mouseDown = true;
//
//    x = mx / cScale;
//    y = (canvas_height - my) / cScale;
//
//    setObstacle(x,y, true);
//    scene.paused = false;
}

void drag(float x, float y)
{
//    if (mouseDown) {
//        let bounds = canvas.getBoundingClientRect();
//        let mx = x - bounds.left - canvas.clientLeft;
//        let my = y - bounds.top - canvas.clientTop;
//        x = mx / cScale;
//        y = (canvas_height - my) / cScale;
//        setObstacle(x,y, false);
//    }
}

void endDrag()
{
    mouseDown = false;
    scene.obstacleVelX = 0.0;
    scene.obstacleVelY = 0.0;
}

//canvas.addEventListener('mousedown', event => {
//    startDrag(event.x, event.y);
//});
//
//canvas.addEventListener('mouseup', event => {
//    endDrag();
//});
//
//canvas.addEventListener('mousemove', event => {
//    drag(event.x, event.y);
//});
//
//canvas.addEventListener('touchstart', event => {
//    startDrag(event.touches[0].clientX, event.touches[0].clientY)
//});
//
//canvas.addEventListener('touchend', event => {
//    endDrag()
//});
//
//canvas.addEventListener('touchmove', event => {
//    event.preventDefault();
//    event.stopImmediatePropagation();
//    drag(event.touches[0].clientX, event.touches[0].clientY)
//}, { passive: false});
//
//
//document.addEventListener('keydown', event => {
//    switch(event.key) {
//        case 'p': scene.paused = !scene.paused; break;
//        case 'm': scene.paused = false; simulate(); scene.paused = true; break;
//    }
//});

void toggleStart()
{
//    var button = document.getElementById('startButton');
//    if (scene.paused)
//        button.innerHTML = "Stop";
//    else
//        button.innerHTML = "Start";
//    scene.paused = !scene.paused;
}

// main -------------------------------------------------------

void simulate()
{
    if (!scene.paused) {
        FlipFluid_simulate(&scene.fluid, &scene);
        scene.frameNr++;
    }
}

void update()
{
    simulate();
    draw_scene();
}

// --- GLFW error callback ---
static void error_callback(int error, const char *description)
{
    fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

// ---------- Input callbacks ----------
static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            mouseDown = 1;
            double mx, my;
            glfwGetCursorPos(window, &mx, &my);
            lastMouseX = mx;
            lastMouseY = my;
            setObstacleFromWindow(mx, my, 1, window);
        } else if (action == GLFW_RELEASE) {
            mouseDown = 0;
            // azzera velocità ostacolo come in endDrag()
            scene.obstacleVelX = 0.0f;
            scene.obstacleVelY = 0.0f;
        }
    }
}

static void cursor_pos_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (mouseDown) {
        lastMouseX = xpos;
        lastMouseY = ypos;
        setObstacleFromWindow(xpos, ypos, 0, window);
    }
}

void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS) {
        if (key == GLFW_KEY_ESCAPE) glfwSetWindowShouldClose(window, GLFW_TRUE);
        if (key == GLFW_KEY_G) scene.showGrid = !scene.showGrid;
        if (key == GLFW_KEY_P) scene.showParticles = !scene.showParticles;
        if (key == GLFW_KEY_O) scene.showObstacle = !scene.showObstacle;
        if (key == GLFW_KEY_Z) scene.compensateDrift = !scene.compensateDrift;
        if (key == GLFW_KEY_X) scene.separateParticles = !scene.separateParticles;
        if (key == GLFW_KEY_SPACE) scene.paused = !scene.paused;
        if (key == GLFW_KEY_UP) {
            if (scene.flipRatio < 1.0f)
                scene.flipRatio += 0.1f;
            printf("%f\n", scene.flipRatio);
        }
        if (key == GLFW_KEY_DOWN) {
            if (scene.flipRatio > 0.0f)
                scene.flipRatio -= 0.1f;
            printf("%f\n", scene.flipRatio);
        }
    }
}

void sceneDispose()
{
    FlipFluid_dtor(&scene.fluid);
}

int main()
{
    canvas_width = 1280;
    canvas_height = 720;

    simHeight = 3.0;
    cScale = canvas_height / simHeight;
    simWidth = canvas_width / cScale;

    glfwSetErrorCallback(error_callback);
    if (!glfwInit()) {
        fprintf(stderr, "Failed to init GLFW\n");
        return -1;
    }

    // Request OpenGL 3.3 core
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(canvas_width, canvas_height, "FLIP Fluid", NULL, NULL);
    if (!window) {
        fprintf(stderr, "Failed to create GLFW window\n");
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_pos_callback);

    // Init GLEW
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        fprintf(stderr, "GLEW init error: %s\n", glewGetErrorString(err));
        glfwTerminate();
        return -1;
    }
    GLuint vao;
    glGenVertexArrays(1,&vao);
    glBindVertexArray(vao);

    // OpenGL state
    glDisable(GL_DEPTH_TEST); //glEnable(GL_DEPTH_TEST);

    setupScene();

    // main loop
    while (!glfwWindowShouldClose(window)) {
        // update viewport size in case of resize
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        canvas_width = w;
        canvas_height = h;

        // draw
        update();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // cleanup GL resources
    if (pointProgram) glDeleteProgram(pointProgram);
    if (meshProgram) glDeleteProgram(meshProgram);

    if (pointVertexBuffer) glDeleteBuffers(1, &pointVertexBuffer);
    if (pointColorBuffer) glDeleteBuffers(1, &pointColorBuffer);
    if (gridVertBuffer) glDeleteBuffers(1, &gridVertBuffer);
    if (gridColorBuffer) glDeleteBuffers(1, &gridColorBuffer);
    if (diskVertBuffer) glDeleteBuffers(1, &diskVertBuffer);
    if (diskIdBuffer) glDeleteBuffers(1, &diskIdBuffer);

    sceneDispose();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
