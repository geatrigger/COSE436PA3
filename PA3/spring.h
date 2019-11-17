
#pragma once
 
#include <iostream>
#include "Node.h"
#include <vector>
using namespace std;

class mass_spring{
public:
	double		spring_coef;
	double		damping_coef;
	Node	*p1;
	Node	*p2;
	double		initial_length;
 
public:
 
	mass_spring(Node *p1, Node *p2)
	{
		damping_coef = 5.0;
		this->p1 = p1;
		this->p2 = p2;
		init();
	}

	void init()
	{
		vec3 S_length = (p2->position - p1->position);
		initial_length = S_length.length();
	}

	void internal_force(double dt)
	{
		vec3 xj_xi = p1->position - p2->position;
		vec3 vj_vi = p1->velocity - p2->velocity;
		vec3 spring_force = xj_xi.normalize() * (spring_coef*(xj_xi.length() - initial_length) +
			damping_coef * vj_vi.dot(xj_xi.normalize()));
		p2->add_force(spring_force);
		p1->add_force(vec3(0,0,0)-spring_force);
	}
	void draw();

};