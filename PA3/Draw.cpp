#include "Node.h"
#include <GL/glut.h>
#include "cloth.h"
#include "sphere.h"
#include "Particle.h"

void Node::draw()
{
	glDisable(GL_LIGHTING);
	glColor3f(0.97, 0.95, 0.15);
	glPointSize(2.0);

	glBegin(GL_POINTS);	
	glVertex3f(getPosX(), getPosY(), getPosZ());
	glEnd();
	glEnable(GL_LIGHTING);
}

void mass_spring::draw()
{
	glDisable(GL_LIGHTING);
	glColor3f(1.0, 1.0, 1.0);
	glLineWidth(2.0);

 	glBegin(GL_LINES);
	glVertex3f(p1->position.x, p1->position.y, p1->position.z);
	glVertex3f(p2->position.x, p2->position.y, p2->position.z);
  glEnd();	 
	glEnable(GL_LIGHTING);

}


void mass_cloth::draw()
{	
	switch (drawMode)
	{
	case DRAW_MASS_NODES:
		glDisable(GL_LIGHTING);
		for (int i = 0; i < nodes.size(); i++)
			nodes[i]->draw();
		glEnable(GL_LIGHTING);
		break;
	case DRAW_SPRINGS:
		glDisable(GL_LIGHTING);
		for (int i = 0; i < spring.size(); i++)
			spring[i]->draw();
		glEnable(GL_LIGHTING);
		break;
	case DRAW_FACES:
		glEnable(GL_TEXTURE_2D);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		glBindTexture(GL_TEXTURE_2D, texName);;
		
		for (int i = 0; i < faces.size(); i += 3)
		{
			float clean_x[3];
			for (int j = 0; j < 3; j++)
			{
				clean_x[j] = 0;
				if (faces[i+j]->isCleaned)
				{
					clean_x[j] = 0.5f;
				}
			}
			glBegin(GL_TRIANGLES);
			//glColor4f(0.8f, 0.6f, 1.0f, 1.0f);
			glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
			glTexCoord2f(faces[i]->tex_num_x / (float)size_x / 2.0f + clean_x[0], faces[i]->tex_num_y / (float)size_y);
			glNormal3f(faces[i]->normal.x, faces[i]->normal.y, faces[i]->normal.z);
			glVertex3f(faces[i]->position.x, faces[i]->position.y, faces[i]->position.z);
			
			glTexCoord2f(faces[i + 1]->tex_num_x / (float)size_x / 2.0f + clean_x[1], faces[i + 1]->tex_num_y / (float)size_y);
			glNormal3f(faces[i + 1]->normal.x, faces[i + 1]->normal.y, faces[i + 1]->normal.z);
			glVertex3f(faces[i + 1]->position.x, faces[i + 1]->position.y, faces[i + 1]->position.z);
			
			glTexCoord2f(faces[i + 2]->tex_num_x / (float)size_x / 2.0f + clean_x[2], faces[i + 2]->tex_num_y / (float)size_y);
			glNormal3f(faces[i + 2]->normal.x, faces[i + 2]->normal.y, faces[i + 2]->normal.z);
			glVertex3f(faces[i + 2]->position.x, faces[i + 2]->position.y, faces[i + 2]->position.z);
			
			glEnd();
		}
		for (int i = 0; i < nodes.size(); i++)
		{
			if(nodes[i]->hasParticle)
			  nodes[i]->pp->draw();
		}
		glDisable(GL_TEXTURE_2D);
		//Basic Implements 3-3. Draw Call for Cloth
		//Additional Implements 4-3. Texture Coordinate Mapping
		break;
	default:
		break;
	}
	glPopMatrix();
}


void mass_sphere::draw()
{
	switch (drawMode)
	{
	case DRAW_MASS_NODES:
		break;
	case DRAW_SPRINGS:
		break;
	case DRAW_FACES:
		glColor3f(0.5, 1.0, 0.5);
		glTranslatef(origin.x, origin.y, origin.z);
		gluSphere(gl_sphere, size, 50, 10);
		break;
	default:
		break;
	}
}


void Particle::draw()
{
	if (idx < 0)
	{
		glDisable(GL_LIGHTING);
		glColor3f(0.0f, 0.0f, 0.0f);
		glPointSize(5.0f);
		glBegin(GL_POINTS);
		glVertex3f(getPosX(), getPosY(), getPosZ());
		glEnd();
		glEnable(GL_LIGHTING);
		return;
	}
	//glDisable(GL_LIGHTING);
	glPointSize(25.0f);
	glEnable(GL_POINT_SMOOTH);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA_SATURATE);
	glBegin(GL_POINTS);
	glColor4f(0.3f, 0.7f, 0.9f, 0.2f);
	glNormal3f(normal.x, normal.y, normal.z);
	glVertex3f(getPosX(), getPosY(), getPosZ());
	glEnd();
	int fake_particle_num = 1;
	for (int i = 0; i < neighbor_positions.size(); i++)
	{
		for (int j = 1; j < fake_particle_num; j++)
		{
			glBegin(GL_POINTS);
			glColor4f(0.5f, 0.5f, 1.0f, 0.3f);
			glNormal3f((normal.x * j + neighbor_normals[i].x * (fake_particle_num - j)), (normal.y * j + neighbor_normals[i].y * (fake_particle_num - j)), (normal.z * j + neighbor_normals[i].z * (fake_particle_num - j)));
			glVertex3f(getPosX() * j + neighbor_positions[i].x * (fake_particle_num - j), getPosY() * j + neighbor_positions[i].y * (fake_particle_num - j), getPosZ() * j + neighbor_positions[i].z * (fake_particle_num - j));
			glEnd();
		}

	}
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glEnable(GL_LIGHTING);
}