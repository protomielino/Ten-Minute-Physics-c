#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include <raylib.h>
#define RAYMATH_IMPLEMENTATION
#include <raymath.h>
#define RAYGUI_IMPLEMENTATION
#include "../raygui.h"


// drawing -------------------------------------------------------

#define window_innerWidth 1200
#define window_innerHeight 700
#define canvas_height (window_innerHeight - 100)

#define canvas_width (window_innerWidth - 20)
#define canvas_height (window_innerHeight - 100)

static double simMinWidth;
static double cScale;
static double simWidth;
static double simHeight;

// vector math -------------------------------------------------------

typedef struct
{
    double x;
    double y;
} v2;
v2 v2_constructor(double x /*= 0.0*/, double y /*= 0.0*/)
{
    v2 this = {0};
    this.x = x;
    this.y = y;
    return this;
}

void v2_set(v2 *this, v2 v)
{
    this->x = v.x;
    this->y = v.y;
}

v2 v2_clone(v2 *this)
{
    return v2_constructor(this->x, this->y);
}

v2* v2_add(v2 *this, v2 v, double s /*= 1.0*/)
{
    this->x += v.x * s;
    this->y += v.y * s;
    return this;
}

v2* v2_addVectors(v2 *this, v2 a, v2 b)
{
    this->x = a.x + b.x;
    this->y = a.y + b.y;
    return this;
}

v2* v2_subtract(v2 *this, v2 v, double s /*= 1.0*/)
{
    this->x -= v.x * s;
    this->y -= v.y * s;
    return this;
}

v2* v2_subtractVectors(v2 *this, v2 a, v2 b)
{
    this->x = a.x - b.x;
    this->y = a.y - b.y;
    return this;
}

double v2_length(v2 *this)
{
    return sqrt(this->x * this->x + this->y * this->y);
}

v2* v2_scale(v2 *this, double s)
{
    this->x *= s;
    this->y *= s;
    return this;
}

double v2_dot(v2 *this, v2 v)
{
    return this->x * v.x + this->y * v.y;
}

v2 v2_perp(v2 *this)
{
    return v2_constructor(-this->y, this->x);
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

typedef struct Bead_s
{
    double radius;
    double mass;
    v2 pos;
    v2 prevPos;
    v2 vel;
} Bead;

typedef struct
{
    double radius;
    double beadRadius;
    double mass;
    double angle;
    double omega;
} AnalyticBead;

typedef struct physicsScene_s
{
    v2 gravity;
    double dt;
    int numSteps;
    bool paused;
    v2 wireCenter;
    double wireRadius;
    Bead bead;
    AnalyticBead analyticBead;
} physicsScene;
physicsScene scene =
{
    gravity : (v2){0.0, -10.0},
    dt : 1.0 / 60.0,
    numSteps : 1000.0,
    paused : false,
    wireCenter : (v2){0, 0},
    wireRadius : 0.0,
    bead : {0},
    analyticBead : {0}
};

// -------------------------------------------------------

Bead Bead_constructor(double radius, double mass, v2 pos)
{
    Bead this = {0};
    this.radius = radius;
    this.mass = mass;
    this.pos = v2_clone(&pos);
    this.prevPos = v2_clone(&pos);
    this.vel = (v2){0,0};
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

// -------------------------------------------------------

AnalyticBead AnalyticBead_constructor(double radius, double beadRadius, double mass, double angle)
{
    AnalyticBead this = {0};
    this.radius = radius;
    this.beadRadius = beadRadius;
    this.mass = mass;
    this.angle = angle;
    this.omega = 0.0;
    return this;
}
double AnalyticBead_simulate(AnalyticBead *this, double dt, double gravity)
{
    double acc = -gravity / this->radius * sin(this->angle);
    this->omega += acc * dt;
    this->angle += this->omega * dt;

    double centrifugalForce = this->omega * this->omega * this->radius;
    double force = centrifugalForce + cos(this->angle) * fabsf(gravity);
    return force;
}
v2 AnalyticBead_getPos(AnalyticBead *this)
{
    return (v2){
         sin(this->angle) * this->radius,
        -cos(this->angle) * this->radius
    };
}

// -----------------------------------------------------

void setupScene()
{
    scene.paused = true;

    scene.wireCenter.x = simWidth / 2.0;
	scene.wireCenter.y = simHeight / 2.0;
	scene.wireRadius = simMinWidth * 0.4;

	v2 pos = {
	        scene.wireCenter.x + scene.wireRadius,
	        scene.wireCenter.y
	};

	scene.bead = Bead_constructor(0.1, 1.0, pos);

	scene.analyticBead = AnalyticBead_constructor(
	        scene.wireRadius, 0.1, 1.0, 0.5 * M_PI);

//	document.getElementById("force").innerHTML = 0.0.toFixed(3);
//	document.getElementById("aforce").innerHTML = 0.0.toFixed(3);
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
//	ClearBackground(BLACK);

//	c.lineWidth = 2.0;
	drawCircle(scene.wireCenter, scene.wireRadius, false, WHITE);

	Bead bead = scene.bead;
	drawCircle(bead.pos, bead.radius, true, RED);

	AnalyticBead analyticBead = scene.analyticBead;
	v2 pos = AnalyticBead_getPos(&analyticBead);
	v2_add(&pos, scene.wireCenter, 1.0);
	drawCircle(pos, analyticBead.beadRadius, true, GREEN);
}

// ------------------------------------------------

void simulate()
{
	if (scene.paused)
		return;

	double sdt = scene.dt / scene.numSteps;
	double force, analyticForce;

	for (int step = 0; step < scene.numSteps; step++) {
	    Bead_startStep(&scene.bead, sdt, scene.gravity);

	    double lambda = Bead_keepOnWire(&scene.bead,
		        scene.wireCenter, scene.wireRadius);

		force = fabsf(lambda / sdt / sdt);

		Bead_endStep(&scene.bead, sdt);

		analyticForce = AnalyticBead_simulate(&scene.analyticBead, sdt, -scene.gravity.y);
	}

//	document.getElementById("force").innerHTML = force.toFixed(3);
//	document.getElementById("aforce").innerHTML = analyticForce.toFixed(3);
}

// --------------------------------------------------------

void run()
{
    scene.paused = false;
}

void step()
{
    scene.paused = false;
	simulate();
	scene.paused = true;
}

void update()
{
	simulate();
	draw();
}
	
int main()
{
	simMinWidth = 2.0;
	cScale = fmin(canvas_width, canvas_height) / simMinWidth;
	simWidth = canvas_width / cScale;
	simHeight = canvas_height / cScale;

	setupScene();

	InitWindow(window_innerWidth, window_innerHeight, "Bead");
	SetTargetFPS(30);

	while (!WindowShouldClose()) {
	    if (IsKeyPressed(KEY_SPACE)) scene.paused = !scene.paused;
	    BeginDrawing(); {
	        ClearBackground(BLACK);

	        update();
	    } EndDrawing();
	}

	CloseWindow();

	//  <button class="button" onclick="setupScene()">Restart</button>
	//  <button class="button" onclick="run()">Run</button>
	//  <button class="button" onclick="step()">Step</button>
	//  PBD <span id = "force">0</span>  &emsp; Analytic <span id ="aforce">0</span>

	return 0;
}
