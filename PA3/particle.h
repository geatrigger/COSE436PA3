#pragma once
#include <math.h>
#include <vector>
#include <algorithm>
#include "vector.h"
class Particle
{
public:
	double	mass;
	vec3	position;
	vec3	velocity;
	vec3	acceleration;
	double	density;
	int		idx;
	vec3	fpressure;
	vec3	fviscosity;
	double  restitution;
public:
	Particle(void)
	{
	}
	Particle(double _x, double _y, double _z)
	{
		position = vec3(_x, _y, _z);
		velocity = vec3(0.0, 0.0, 0.0);
		mass = 1.0;
		restitution = 0.5;
	}

	Particle(double _x, double _y, double _z, int _idx) : position(_x, _y, _z), velocity(0.0, 0.0, 0.0), acceleration(0.0, 0.0, 0.0), mass(1.0)
	{
		fpressure = vec3(0.0, 0.0, 0.0);
		fviscosity = vec3(0.0, 0.0, 0.0);
		density = 0.0;
		idx = _idx;
		restitution = 0.5;
	}
	~Particle(void)
	{
	}

	double	getPosX(void) { return position.getX(); }
	double	getPosY(void) { return position.getY(); }
	double	getPosZ(void) { return position.getZ(); }

	void integrate(double dt, vec3 gravity)
	{
		vec3 fgrav = gravity * mass;

		// Update velocity and position
		acceleration = (fpressure + fviscosity) / density + fgrav;
		velocity = velocity + acceleration * dt;
		position = position + velocity * dt;

		// Boundary condition
		if (position.x < -10.0 && velocity.x < 0.0)
		{
			velocity.x *= -restitution;
			position.x = -10.0+0.1;
		}
		if (position.x > 10.0 && velocity.x > 0.0)
		{
			velocity.x *= -restitution;
			position.x = 10.0-0.1;
		}
		if (position.y < -10.0 && velocity.y < 0.0)
		{
			velocity.y *= -restitution;
			position.y = -10.0+0.1;
		}
		if (position.y > 10.0 && velocity.y > 0.0)
		{
			velocity.y *= -restitution;
			position.y = 10.0-0.1;
		}
		if (position.z < -10.0 && velocity.z < 0.0)
		{
			velocity.z *= -restitution;
			position.z = -10.0 + 0.1;
		}
		if (position.z > 10.0 && velocity.z > 0.0)
		{
			velocity.z *= -restitution;
			position.z = 10.0 - 0.1;
		}
	}

	void draw();
};
