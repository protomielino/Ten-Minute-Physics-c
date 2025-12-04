#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <raylib.h>

#define STB_DS_IMPLEMENTATION
#include "../stb_ds.h"

//    <button class="button" onclick="setupScene()">Restart</button>
//    <br>
//    <canvas id="myCanvas"></canvas>



// drawing -------------------------------------------------------
#define window_innerWidth 1200
#define window_innerHeight 700

double canvas_width = window_innerWidth - 20;
double canvas_height = window_innerHeight - 100;

double simMinWidth;
double cScale;
double simWidth;
double simHeight;

// vector math -------------------------------------------------------

typedef struct
{
    double x;
    double y;
} v2;
v2 v2_constructor(double x /*= 0.0*/, double y /*= 0.0*/)
{
    return (v2){x,y};
}

void v2_set(v2 *this, v2 v)
{
    this->x = v.x; this->y = v.y;
}

v2 v2_clone(v2 *this)
{
    return (v2){this->x, this->y};
}

v2 v2_add(v2 *this, v2 v, double s /*= 1.0*/)
{
    this->x += v.x * s;
    this->y += v.y * s;
    return *this;
}

v2 v2_addVectors(v2 *this, v2 a, v2 b)
{
    this->x = a.x + b.x;
    this->y = a.y + b.y;
    return *this;
}

v2 v2_subtract(v2 *this, v2 v, double s /*= 1.0*/)
{
    this->x -= v.x * s;
    this->y -= v.y * s;
    return *this;
}

v2 v2_subtractVectors(v2 *this, v2 a, v2 b)
{
    this->x = a.x - b.x;
    this->y = a.y - b.y;
    return *this;
}

double v2_length(v2 *this)
{
    return sqrt(this->x * this->x + this->y * this->y);
}

v2 v2_scale(v2 *this, double s)
{
    this->x *= s;
    this->y *= s;
    return *this;
}

double v2_dot(v2 *this, v2 v)
{
    return this->x * v.x + this->y * v.y;
}

v2 v2_perp(v2 *this)
{
    return (v2){-this->y, this->x};
}

double cX(v2 pos)
{
    return pos.x * cScale;
}

double cY(v2 pos)
{
    return canvas_height - pos.y * cScale;
}

// scene -------------------------------------------------------

typedef struct
{
    double radius;
    double mass;
    v2 pos;
    v2 prevPos;
    v2 vel;
} Bead;

typedef struct
{
    v2     gravity;
    double  dt;
    int    numSteps;
    v2     wireCenter;
    double  wireRadius;
    Bead  *beads;
} physicsScene;
physicsScene scene =
{
    gravity    : (v2){0.0, -10.0},
    dt         : 1.0 / 60.0,
    numSteps   : 100,
    wireCenter : (v2){0},
    wireRadius : 0.0,
    beads      : NULL,
};

Bead Bead_constructor(double radius, double mass, v2 pos)
{
    Bead this = {0};
    this.radius = radius;
    this.mass = mass;
    this.pos = v2_clone(&pos);
    this.prevPos = v2_clone(&pos);
    this.vel = (v2){0};
    return this;
}
void Bead_startStep(Bead *this, double dt, v2 gravity)
{
    v2_add(&this->vel, gravity, dt);
    v2_set(&this->prevPos, this->pos);
    v2_add(&this->pos, this->vel, dt);
}
double Bead_keepOnWire(Bead *this, v2 center, double radius)
{
    v2 dir = {0};
    v2_subtractVectors(&dir, this->pos, center);
    double len = v2_length(&dir);
    if (len == 0.0)
        return 0.0;
    v2_scale(&dir, 1.0 / len);
    double lambda = scene.wireRadius - len;
    v2_add(&this->pos, dir, lambda);
    return lambda;
}
void Bead_endStep(Bead *this, double dt)
{
    v2_subtractVectors(&this->vel, this->pos, this->prevPos);
    v2_scale(&this->vel, 1.0 / dt);
}

// -----------------------------------------------------

void setupScene()
{
    scene.beads = NULL;

    scene.wireCenter.x = simWidth / 2.0;
    scene.wireCenter.y = simHeight / 2.0;
    scene.wireRadius   = simMinWidth * 0.4;

    int numBeads = 5;
    double mass = 1.0;

    double r = 0.1;
    double angle = 0.0;
    for (int i = 0; i < numBeads; i++) {
        mass = M_PI * r * r;
        v2 pos = {
                scene.wireCenter.x + scene.wireRadius * cos(angle),
                scene.wireCenter.y + scene.wireRadius * sin(angle)
        };

        arrput(scene.beads, Bead_constructor(r, mass, pos));
        angle += M_PI / numBeads;
        r = 0.05 + drand48() * 0.1;
    }
}

// draw -------------------------------------------------------

void drawCircle(v2 pos, double radius, bool filled, Color col)
{
    if (filled)
        DrawCircle(cX(pos), cY(pos), cScale * radius, col);
    else
        DrawCircleLines(cX(pos), cY(pos), cScale * radius, col);
}

void draw()
{
//    ClearBackground(BLACK);

    drawCircle(scene.wireCenter, scene.wireRadius, false, WHITE);

    for (int i = 0; i < arrlen(scene.beads); i++) {
        Bead bead = scene.beads[i];
        drawCircle(bead.pos, bead.radius, true, RED);
    }
}

// --- collision handling -------------------------------------------------------

void handleBeadBeadCollision(Bead *bead1, Bead *bead2)
{
    double restitution = 1.0;
    v2 dir = {0};
    v2_subtractVectors(&dir, bead2->pos, bead1->pos);
    double d = v2_length(&dir);
    if (d == 0.0 || d > bead1->radius + bead2->radius)
        return;

    v2_scale(&dir, 1.0 / d);

    double corr = (bead1->radius + bead2->radius - d) / 2.0;
    v2_add(&bead1->pos, dir, -corr);
    v2_add(&bead2->pos, dir, corr);

    double v1 = v2_dot(&bead1->vel, dir);
    double v2 = v2_dot(&bead2->vel, dir);

    double m1 = bead1->mass;
    double m2 = bead2->mass;

    double newV1 = (m1 * v1 + m2 * v2 - m2 * (v1 - v2) * restitution) / (m1 + m2);
    double newV2 = (m1 * v1 + m2 * v2 - m1 * (v2 - v1) * restitution) / (m1 + m2);

    v2_add(&bead1->vel, dir, newV1 - v1);
    v2_add(&bead2->vel, dir, newV2 - v2);
}

// ------------------------------------------------

void simulate()
{
    double sdt = scene.dt / scene.numSteps;

    for (int step = 0; step < scene.numSteps; step++) {
        for (int i = 0; i < arrlen(scene.beads); i++)
            Bead_startStep(&scene.beads[i], sdt, scene.gravity);

        for (int i = 0; i < arrlen(scene.beads); i++) {
            Bead_keepOnWire(&scene.beads[i], scene.wireCenter, scene.wireRadius);
        }

        for (int i = 0; i < arrlen(scene.beads); i++)
            Bead_endStep(&scene.beads[i], sdt);

        for (int i = 0; i < arrlen(scene.beads); i++) {
            for (int j = 0; j < i; j++) {
                handleBeadBeadCollision(&scene.beads[i], &scene.beads[j]);
            }
        }
    }
}

// --------------------------------------------------------

void update()
{
    simulate();
    draw();
}

int main()
{
    srand48(clock());

    simMinWidth = 2.0;
    cScale = fmin(canvas_width, canvas_height) / simMinWidth;
    simWidth = canvas_width / cScale;
    simHeight = canvas_height / cScale;

    setupScene();

    InitWindow(window_innerWidth, window_innerHeight, "Constrained Dynamics - 05_manyBeads");
    SetTargetFPS(30);

    while (!WindowShouldClose()) {
        BeginDrawing(); {
            ClearBackground(BLACK);
            update();
        } EndDrawing();
    }

    arrfree(scene.beads);
    CloseWindow();
    return 0;
}
