#pragma once

#include "particle.h"
#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <GL/glut.h>


class particle_system
{
public:
	std::vector<particle *> particles;

	int			size_x, size_y, size_z;
	double		dx, dy, dz;
	double		structural_coef;
	double		shear_coef;
	double		bending_coef;
	int			iteration_n;
	int			drawMode;


	particle_system()
	{
	}
	~particle_system()
	{
	}
	enum DrawModeEnum {
		DRAW_MASS_NODES,
		DRAW_SPRINGS,
		DRAW_FACES
	};

public:
	void init(int n)
	{
		particles.clear();
		for (int i = 0; i < n; i++)
		{
			particle* temp;
			particles.push_back(temp);
		}
		for (int i = 0; i < n; i++)
		{
			particles[i]->init();
		}
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