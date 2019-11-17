
#pragma once

#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <GL/glut.h>

class mass_sphere
{
public:
	//std::vector<Node*> nodes;
	//std::vector<Node*> faces;
	GLUquadric* gl_sphere;

	vec3 origin;
	double			size;
	//int node_per_2pi;
	int			drawMode;


	mass_sphere()
	{
	}
	~mass_sphere()
	{
	}
	enum DrawModeEnum {
		DRAW_MASS_NODES,
		DRAW_SPRINGS,
		DRAW_FACES
	};

public:
	void init()
	{
		gl_sphere = gluNewQuadric();
	}
	void draw();
};