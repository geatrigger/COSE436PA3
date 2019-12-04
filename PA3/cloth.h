
#pragma once
#define checkImageWidth 64
#define checkImageHeight 64

#include "spring.h"
#include "Node.h"
#include "sphere.h"
#include "Particle.h"
#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <GL/glut.h>
#include <GL/stb_image.h>

class mass_cloth
{
public:

	std::vector<Node *> nodes;
	std::vector<mass_spring *> spring;
	std::vector<Node*> faces;
	std::vector<Particle*> particles;

	int			size_x, size_y, size_z;
	double		dx, dy, dz;
	double		structural_coef;
	double		shear_coef;
	double		bending_coef;
	int			iteration_n;
	int			drawMode;

	GLubyte checkImage[checkImageHeight][checkImageWidth][4];
	GLuint texName;

	mass_cloth()
	{ 	 
	}
	~mass_cloth()
	{ 
		for (int i = 0; i < nodes.size(); i++){ delete nodes[i]; }
		for (int i = 0; i < spring.size(); i++){ delete spring[i]; }
		nodes.clear();
		spring.clear();
		faces.clear();
	}
	enum DrawModeEnum{
		DRAW_MASS_NODES,
		DRAW_SPRINGS,
		DRAW_FACES
	};
 
public:
	void makeCheckImage(void)
	{
		int i, j, c;

		for (i = 0; i < checkImageHeight; i++) {
			for (j = 0; j < checkImageWidth; j++) {
				c = ((((i & 0x8) == 0) ^ ((j & 0x8)) == 0)) * 255;
				//c = 0;
				checkImage[i][j][0] = (GLubyte)c;
				checkImage[i][j][1] = (GLubyte)c;
				checkImage[i][j][2] = (GLubyte)c;
				checkImage[i][j][3] = (GLubyte)255;
			}
		}
	}
	void init()
	{
		//makeCheckImage();

		glGenTextures(1, &texName);
		glBindTexture(GL_TEXTURE_2D, texName);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
			GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
			GL_NEAREST);
		int width, height, nrChannels;
		unsigned char* data = stbi_load("1.jpg", &width, &height, &nrChannels, 0);
		if (data)
		{
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
		}
		else
		{
			std::cout << "Failed to load texture" << std::endl;
		}
		stbi_image_free(data);
		//Node 배치
		for (int i = 0; i < size_x; i++)
		{
			for (int j = 0; j < size_y; j++)
			{
				for (int k = 0; k < size_z; k++)
				{
					Node* xp = new Node(vec3((i - size_x / 2.0) * dx, (k - size_z / 2.0) * dz + (size_y / 2.0) * dy - 10, (j - size_y / 2.0) * dy));
					Particle* pp = new Particle(xp->getPosX(), xp->getPosY(), xp->getPosZ(), -1);
					//pp->mass = 0.5;
					xp->pp = pp;
					particles.push_back(pp);
					if ((j == 0 || j == size_y - 1) && (i == 0 || i == size_x - 1))
						xp->isFixed = true;
					else
						xp->isFixed = false;
					xp->tex_num_x = i;
					xp->tex_num_y = j;
					xp->tex_num_z = k;
					nodes.push_back(xp);
				}
			}
		}
		mass_spring* sp;
		int spring_index = 0;
		//structural x axis
		for (int i = 0; i < size_x - 1; i++)
		{
			for (int j = 0; j < size_y; j++)
			{
				for (int k = 0; k < size_z; k++)
				{
					int index = i * size_y * size_z + j * size_z + k;
					sp = new mass_spring(nodes[index], nodes[index + size_z * size_y * 1]);
					sp->spring_coef = structural_coef;
					spring.push_back(sp);
					nodes[index]->add_spring(spring_index);
					nodes[index + size_z * size_y * 1]->add_spring(spring_index);
					spring_index++;
				}
			}
		}
		//structural y axis
		for (int i = 0; i < size_x; i++)
		{
			for (int j = 0; j < size_y - 1; j++)
			{
				for (int k = 0; k < size_z; k++)
				{
					int index = i * size_y * size_z + j * size_z + k;
					sp = new mass_spring(nodes[index], nodes[index + size_z * 1]);
					sp->spring_coef = structural_coef;
					spring.push_back(sp);
					nodes[index]->add_spring(spring_index);
					nodes[index + size_z * 1]->add_spring(spring_index);
					spring_index++;
				}
			}
		}
		//structural z axis
		for (int i = 0; i < size_x; i++)
		{
			for (int j = 0; j < size_y; j++)
			{
				for (int k = 0; k < size_z - 1; k++)
				{
					int index = i * size_y * size_z + j * size_z + k;
					sp = new mass_spring(nodes[index], nodes[index + 1]);
					sp->spring_coef = structural_coef;
					spring.push_back(sp);
					nodes[index]->add_spring(spring_index);
					nodes[index + 1]->add_spring(spring_index);
					spring_index++;
				}
			}
		}
		//shear +x+y
		for (int i = 0; i < size_x - 1; i++)
		{
			for (int j = 0; j < size_y - 1; j++)
			{
				for (int k = 0; k < size_z; k++)
				{
					int index = i * size_y * size_z + j * size_z + k;
					sp = new mass_spring(nodes[index], nodes[index + size_y * size_z + size_z]);
					sp->spring_coef = shear_coef;
					spring.push_back(sp);
					nodes[index]->add_spring(spring_index);
					nodes[index + size_y * size_z + size_z]->add_spring(spring_index);
					spring_index++;
				}
			}
		}
		//shear +x-y
		for (int i = 0; i < size_x - 1; i++)
		{
			for (int j = 1; j < size_y; j++)
			{
				for (int k = 0; k < size_z; k++)
				{
					int index = i * size_y * size_z + j * size_z + k;
					sp = new mass_spring(nodes[index], nodes[index + size_y * size_z - size_z]);
					sp->spring_coef = shear_coef;
					spring.push_back(sp);
					nodes[index]->add_spring(spring_index);
					nodes[index + size_y * size_z - size_z]->add_spring(spring_index);
					spring_index++;
				}
			}
		}
		//bend x axis
		for (int i = 0; i < size_x - 2; i++)
		{
			for (int j = 0; j < size_y; j++)
			{
				for (int k = 0; k < size_z; k++)
				{
					int index = i * size_y * size_z + j * size_z + k;
					sp = new mass_spring(nodes[index], nodes[index + size_z * size_y * 2]);
					sp->spring_coef = structural_coef;
					spring.push_back(sp);
					nodes[index]->add_spring(spring_index);
					nodes[index + size_z * size_y * 2]->add_spring(spring_index);
					spring_index++;
				}
			}
		}
		//bend y axis
		for (int i = 0; i < size_x; i++)
		{
			for (int j = 0; j < size_y - 2; j++)
			{
				for (int k = 0; k < size_z; k++)
				{
					int index = i * size_y * size_z + j * size_z + k;
					sp = new mass_spring(nodes[index], nodes[index + size_z * 2]);
					sp->spring_coef = structural_coef;
					spring.push_back(sp);
					nodes[index]->add_spring(spring_index);
					nodes[index + size_z * 2]->add_spring(spring_index);
					spring_index++;
				}
			}
		}
		//face
		for (int i = 0; i < size_x - 1; i++)
		{
			for (int j = 0; j < size_y - 1; j++)
			{
				for (int k = 0; k < size_z; k++)
				{
					int index = i * size_y * size_z + j * size_z + k;
					//BAD
					faces.push_back(nodes[index + size_y * size_z]);
					faces.push_back(nodes[index]);
					faces.push_back(nodes[index + size_y * size_z + size_z]);
					//CDA
					faces.push_back(nodes[index + size_z]);
					faces.push_back(nodes[index + size_y * size_z + size_z]);
					faces.push_back(nodes[index]);
				}
			}
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

		for (int i = 0; i < faces.size(); i += 3)
		{
			vec3 a = faces[i + 1]->position - faces[i]->position;
			vec3 b = faces[i + 2]->position - faces[i]->position;
			face_normals.push_back(a.Cross(b).normalize());
		}
		for (int i = 0; i < nodes.size(); i++)
		{
			//FIX 191015
			vector<vec3> adjacent_faces;
			int x = i / (size_y * size_z);
			int y = (i - x * size_y * size_z) / size_z;
			int z = (i - x * size_y * size_z - y * size_z);


			if (x != size_x - 1 && y != size_y - 1)
			{
				//right_bottom, right_bottom + 1
				int right_bottom = 2 * x * (size_y - 1) * (size_z) + 2 * y * (size_z) + 2 * z;
				adjacent_faces.push_back(face_normals[right_bottom]);
				adjacent_faces.push_back(face_normals[right_bottom + 1]);
			}
			if (x != 0 && y != 0)
			{
				//left_top, left_top + 1
				int left_top = 2 * (x - 1) * (size_y - 1) * (size_z) + 2 * (y - 1) * (size_z) + 2 * z;
				adjacent_faces.push_back(face_normals[left_top]);
				adjacent_faces.push_back(face_normals[left_top + 1]);
			}
			if(x != size_x - 1 && y != 0)
			{
				//right_top + 1
				int right_top = 2 * x * (size_y - 1) * (size_z) + 2 * (y - 1) * (size_z) + 2 * z;
				adjacent_faces.push_back(face_normals[right_top + 1]);
			}
			if (x != 0 && y != size_y - 1)
			{
				//left_bottom
				int left_bottom = 2 * (x - 1) * (size_y - 1) * (size_z) + 2 * y * (size_z) + 2 * z;
				adjacent_faces.push_back(face_normals[left_bottom]);
			}
			nodes[i]->normal = vec3(0, 0, 0);
			for (int j = 0; j < adjacent_faces.size(); j++)
			{
				nodes[i]->normal += adjacent_faces[j];
			}
			nodes[i]->normal = nodes[i]->normal.normalize();
		}

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
		for (int i = 0; i < nodes.size(); i++)
		{
			nodes[i]->add_force(additional_force);
		}
	}

	void compute_force(double dt, vec3 gravity, vec3 external_force)
	{
		for (int i = 0; i < nodes.size(); i++)
		{
			nodes[i]->add_force(gravity * nodes[i]->mass);
		}
		/* Compute Force for all springs */
		for (int i = 0; i < spring.size(); i++)
		{
			spring[i]->internal_force(dt);
		}
		//external force
		for (int i = 0; i < nodes.size(); i++)
		{
			nodes[i]->add_force(external_force);
			//sph force
			vec3 sph_force = nodes[i]->pp->fpressure + nodes[i]->pp->fviscosity;
			sph_force = sph_force/30;
			//if(i == nodes.size() / 2)
			//  printf("sph force : %lf, %lf, %lf\n", sph_force.getX(), sph_force.getY(), sph_force.getZ());
			nodes[i]->add_force(sph_force);
		}
	}
	

	void integrate(double dt, vec3 gravity, vec3 external_force, int method)
	{
		/* integrate Nodes*/
		for (int i = 0; i < nodes.size(); i++)
		{
			nodes[i]->integrate_euler(dt);
		}
		switch (method)
		{
		case 0:
			//printf("euler\n");
			break;
		case 1:
			//RK2 Method
			compute_force(dt, gravity, external_force);
			for (int i = 0; i < nodes.size(); i++)
			{
				nodes[i]->integrate(dt);
			}
			//printf("rk2\n");
			break;
		default:
			break;
		}
	}
	
	void collision_response(vec3 ground)
	{
		vec3 p, n;
		double r;
		for (int i = 0; i < nodes.size(); i++)
		{
			vec3 x = nodes[i]->position;
			vec3 v = nodes[i]->velocity;
			//ground
			p = ground;
			n = vec3(0, 1, 0);
			r = 0.1;
			if ((x - p).dot(n) < r && n.dot(v) < 0)
			{
				v = v - (2 * n.dot(v) * n);
				nodes[i]->velocity = v;
				nodes[i]->position += n * r;
				continue;
			}
		}

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

	void set_particle_position()
	{
		for (int i = 0; i < nodes.size(); i++)
		{
			Particle* pp = nodes[i]->pp;
			pp->position = vec3(nodes[i]->getPosX(), nodes[i]->getPosY(), nodes[i]->getPosZ());
		}
	}

	void draw();
};