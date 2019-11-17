#pragma once

#define MIN_INIT_VELOCITY 10
#define MAX_INIT_VELOCITY 100
#define LENGTH 30

#include "spring.h"
#include "Node.h"
#include "sphere.h"
#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <GL/glut.h>


class particle
{
public:

	double  mass;
	vec3    velocity, position;
	int			iteration_n;
	int			drawMode;


	particle()
	{
		mass = 1;
		velocity = vec3(0, MIN_INIT_VELOCITY - rand() % (MAX_INIT_VELOCITY - MIN_INIT_VELOCITY), 0);
		position = vec3((1 - 2 * rand_float()) * (LENGTH - 10), (1 - rand_float()) * LENGTH - 10, 0);
	}
	~particle()
	{
	}
	enum DrawModeEnum {
		DRAW_MASS_NODES,
		DRAW_SPRINGS,
		DRAW_FACES
	};

	float rand_float()
	{
		float value = rand() / float(RAND_MAX);
		return value;
	}
public:
	void init()
	{
		//Basic Implements 1. Init Nodes and Shear and Structural Springs
		//Additional Implements 1. Init Bending Spring
		/*
			Node *xp = new Node(vec3(x, y, z));

			mass_spring *sp = new mass_spring(p[Node_Index_A], p[Node_Index_B]);
			sp->spring_coef = spring_Type_coef;
			spring.push_back(sp);
		*/
		//Basic Implements 3-1. Generate Faces
		/*
			faces.push_back(p[Node_Index_A]);
			faces.push_back(p[Node_Index_C]);
			faces.push_back(p[Node_Index_B]);
		*/
		//Additional Implements 4-2. Initialize Texture Coordinates	
	}

	void computeNormal()
	{
		std::vector<vec3> face_normals;

		//Basic Implements 3-2. Compute Vertex Normal
		/*
			for(each face)
			{
				compute face normal
			}
			for(each node)
			{
				������ face�� ��� normal
			}
		*/
	}

	void add_force(vec3 additional_force)
	{
	}

	void compute_force(double dt, vec3 gravity, vec3 external_force)
	{
	}


	void integrate(double dt, vec3 gravity, vec3 external_force, int method)
	{
		switch (method)
		{
		case 0:
			//printf("euler\n");
			break;
		case 1:
			//RK2 Method
			compute_force(dt, gravity, external_force);
			//printf("rk2\n");
			break;
		default:
			break;
		}
	}

	void collision_response(vec3 ground, mass_sphere* sphere)
	{
		//Basic Implements 4. Collision Check with ground
		//Additional Implements 2. Collision Check with Sphere
		//Additional Implements 3. Collision Check with Mesh Object
		/*
			if(Collision Detection)
			{
				Collision Response
			}
		*/

	}

	void draw();
};