

#pragma once
#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "vector.h"
#include "Particle.h"

class Node
{
private:
	vec3	_position;
	vec3	_velocity;
	vec3	_acceleration;
	vec3	position_euler;
	vec3	velocity_euler;
	vec3	acceleration_euler;
public:
	std::vector<int> connect_springs;
	Particle* pp;

	double	mass;
	vec3	force;
	vec3	position;
	vec3	velocity;
	vec3	acceleration;
	vec3	normal;
	bool	isFixed;
	bool  hasParticle;
	bool  isCleaned;

	int tex_num_x;
	int tex_num_y;
	int tex_num_z;

public:
	Node(void)
	{
		isFixed = false;
		//position = vec3(0, 0, 0);
		//velocity = vec3(0, 0, 0);
		//acceleration = vec3(0, 0, 0);
		mass = 1.0;
	}
	Node(vec3 init_pos)
	{
		isFixed = false;
		position = init_pos;
		//velocity = vec3(0, 0, 0);
		//acceleration = vec3(0, 0, 0);
		mass = 1.0;
	}
 
	~Node(void)
	{
	}

	double	getPosX(void) { return position.getX(); }
	double	getPosY(void) { return position.getY(); }
	double	getPosZ(void){ return position.getZ(); }

	void add_spring(int spring_index)
	{
		connect_springs.push_back(spring_index);
	}

	void add_force(vec3 additional_force)
	{
		force += additional_force;
	}

	void integrate(double dt)
	{
		if (!isFixed)
		{
			acceleration = force / mass;
			//RK2 Method
			velocity = _velocity + (acceleration_euler + acceleration) / 2.0 * dt;
			position = _position + (velocity_euler + velocity) / 2.0 * dt;
			/*
			velocity += force / mass * dt;
			position += velocity * dt;
			*/
			//Basic Implements 2-2. Integration
		}
		/*initialize Force*/
		force.x = force.y = force.z= 0.0;
	}
	//find position_euler, velocity_euler
	void integrate_euler(double dt)
	{
		if (!isFixed)
		{
			_position = position;
			_velocity = velocity;
			_acceleration = force / mass;

			acceleration = force / mass;
			velocity += acceleration * dt;
			position += velocity * dt;

			position_euler = position;
			velocity_euler = velocity;
			acceleration_euler = acceleration;
			/*
			velocity += force / mass * dt;
			position += velocity * dt;
			*/
			//Basic Implements 2-2. Integration
		}
		/*initialize Force*/
		force.x = force.y = force.z = 0.0;
	}

	void draw();
};
