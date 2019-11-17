

#pragma once

#include <iostream>
#include <vector>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include "cloth.h"
#include "sphere.h"
 

using namespace std;


class Simulator
{
public:
	Simulator(void);
	~Simulator(void);

public:
	void					Initialize(void);
 	void					Update();
	void					Render();
	void					Lighting(void);
	void					DrawGround(void);

public:
	mass_cloth			*cloth;
	vec3				ground;
	mass_sphere *sphere;
	float timsStep;
	int					  numerical_method;
	vec3          external_force;
};

