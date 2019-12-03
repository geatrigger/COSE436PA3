#include "SPH.h"
#include <time.h>
#include <omp.h>

SPH::SPH()
{

}

SPH::SPH(int numparticle)
{
	MaxParticle = numparticle;
	index = 0;
	//rest_density = 4.2;
	rest_density = 2.0;
	k = 50.0; //gas
	mu = 0.8; //viscosity
	h = 1.0; //kernel radius
	/*
	cpoly6 = 4.0 / PI;
	cpoly6G = -24.0 / PI;
	cpoly6L = -48.0 / PI;
	cspikyG = -30.0 / PI;
	cviscosityL = 40 / PI;
	*/
}

SPH::~SPH()
{
	while (!particles.empty())
	{
		particles.pop_back();
	}
	particles.clear();
}

void SPH::resetParticle()
{
	index = 0;
	while (!particles.empty())
	{
		particles.pop_back();
	}
}

void SPH::init()
{
	resetParticle();
	damBreaking();
}

void SPH::damBreaking()
{
	for (double z = 0.0; z < 2.0; z += 0.5) {
		for (double y = -9.0; y < -2.0; y += 0.5) {
			for (double x = 0.0; x < 2.0; x += 0.5) {
				if (particles.size() < MaxParticle)
				{
					Particle* p = new Particle(x, y, z, index++);
					particles.push_back(p);
				}
			}
		}
	}
	cout << "SPH" << particles.size() << " Paricles" << endl;
}

void SPH::pouring(float dt)
{
	static float pouring_time = 0.0f;
	pouring_time += dt;
	if (pouring_time < 0.1f)
	{
		return;
	}
	if (particles.size() >= MaxParticle)
		return;

	for (double z = -1.0; z < 0.0; z += 0.4) {
		for (double y = -5.0; y < -4.0; y += 0.4) {
			for (double x = -10.0; x < -10.0 + 2.0; x += 0.4) {
				if (particles.size() < MaxParticle)
				{
					Particle* p = new Particle(x, y, z, index++);
					p->velocity.x = 5.0f;
					p->velocity.y = -2.0f;
					particles.push_back(p);
					//printf("x : %lf, y : %lf, z : %lf, thread : %d\n", x, y, z, omp_get_thread_num());
				}
			}
		}
	}
	pouring_time = 0.0f;
	cout << "SPH" << particles.size() << " Paricles" << endl;
}

void SPH::update(float dt, vec3 gravity)
{
	pouring(dt);
	makeHashTable();
	computeDensity();
	computeForce();
	integrate(dt, gravity);
}

void SPH::draw()
{
	for (int i = 0; i < particles.size(); i++)
	{
		Particle *p = particles[i];
		p->draw();
	}
}

// need 3d?
double SPH::poly6Kernel(vec3 rij, double h)
{
	double temp = 0.0;
	double norm = rij.length();
	if (norm == 0)
		return 0;
	temp = h*h - norm * norm;

	return 4.0 / (SPH_PI*h*h*h*h*h*h*h*h) * temp*temp*temp;
}

// need 3d?
vec3 SPH::spikygradientKernel(vec3 rij, double q)
{
	double temp = 0.0;
	temp = 1.0 - q*q;
	temp = -30.0 / (SPH_PI*h*h*h*h) * (temp*temp)/q;

	return  vec3(temp*rij.x, temp*rij.y, temp * rij.z);
}

double SPH::viscositylaplacianKernel(vec3 rij, double q)
{
	return 40.0 / (SPH_PI*h*h*h*h) * (1.0 - q);
}

void SPH::computeDensity()
{
  #pragma omp parallel for collapse(5)
	for (int x = 0; x < GRIDSIZE; x++)
	{
		for (int y = 0; y < GRIDSIZE; y++)
		{
			for (int z = 0; z < GRIDSIZE; z++)
			{
				vector<Particle*> ris;
				vector<Particle*> rjs = getNeighbor(x, y, z, h, ris);
				//printf("x : %d, y : %d, z : %d, thread : %d\n", x, y, z, omp_get_thread_num());

				for (int i = 0; i < ris.size(); i++)
				{
					Particle* pi = ris[i];
					pi->density = 0.0;	//compute with poly6Kernel

					/*Implements - Compute Density 작성, 아래 density 초기화는 지울 것*/
					for (int j = 0; j < rjs.size(); j++)
					{
						Particle* pj = rjs[j];
						vec3 rij = pi->position - pj->position;
						double q = rij.length() / h;
						if (0.0 <= q && q < 1.0)
						{
							pi->density = pi->density + pj->mass * poly6Kernel(rij, h);
							//printf("pi->density : %lf\n", pi->density);
						}
					}

				}
			}
		}
	}
}

void SPH::computeForce() // Compute Pressure and Viscosity
{
#pragma omp parallel for collapse(5)
	for (int x = 0; x < GRIDSIZE; x++)
	{
		for (int y = 0; y < GRIDSIZE; y++)
		{
			for (int z = 0; z < GRIDSIZE; z++)
			{
				vector<Particle*> ris;
				vector<Particle*> rjs = getNeighbor(x, y, z, h, ris);

				for (int i = 0; i < ris.size(); i++)
				{
					Particle* pi = ris[i];
					pi->fpressure = vec3(0, 0, 0);
					pi->fviscosity = vec3(0, 0, 0);
					//printf("pi->fpressur : %lf, pi->fviscosity : %lf\n", pi->fpressure, pi->fviscosity);
					for (int j = 0; j < rjs.size(); j++)
					{
						Particle* pj = rjs[j];
						vec3 rij = pi->position - pj->position;
						double q = rij.length() / h;
						if (0.0 < q && q < 1.0)
						{
							pi->fpressure = pi->fpressure + pj->mass * (k * ((pi->density - rest_density) + (pj->density - rest_density))
								/ (2.0 * pj->density)) * spikygradientKernel(rij, q);//compute with spikygradientKernel
							pi->fviscosity = pi->fviscosity + pj->mass * ((pj->velocity - pi->velocity) / pj->density)
								* viscositylaplacianKernel(rij, q);//compute with viscositylaplacianKernel
							//printf("spiky : (%lf, %lf)\n", spikygradientKernel(rij, q).x, spikygradientKernel(rij, q).y);
							//printf("pj->mass : %lf, k : %lf, pi->density : %lf, pj->density : %lf\n", pj->mass, k, pi->density, pj->density);
						}
					}

					/*Implements - Compute Pressure and Viscosity Forces 작성*/
					pi->fpressure = -1.0 * pi->fpressure;
					pi->fviscosity = mu * pi->fviscosity;
					//printf("pi->fpressur : %lf, pi->fviscosity : %lf\n", pi->fpressure, pi->fviscosity);
				}
			}
		}
	}
}

void SPH::integrate(double dt, vec3 gravity)
{
	for (int i = 0; i < particles.size(); i++)
	{
		Particle *p = particles[i];
		p->integrate(dt, gravity);
		//printf("X : %lf, Y : %lf, Z : %lf\n", p->getPosX(), p->getPosY(), p->getPosZ());
	}
}

void SPH::makeHashTable()
{
#pragma omp parallel for collapse(3)
	for (int p = 0; p < GRIDSIZE; p++)
	{
		for (int q = 0; q < GRIDSIZE; q++)
		{
			for (int r = 0; r < GRIDSIZE; r++)
			{
				hashGrid[p][q][r].clear();
			}
		}
	}

	for (int i = 0; i < particles.size(); i++)
	{
		Particle *p = particles[i];
		double x = (p->getPosX() + GRIDSIZE / 2);
		double y = (p->getPosY() + GRIDSIZE / 2);
		double z = (p->getPosZ() + GRIDSIZE / 2);
		int gridx = (int)(x);					
		int gridy = (int)(y);
		int gridz = (int)(z);

		if (gridx < 0) gridx = 0;
		if (gridx > GRIDSIZE - 1) gridx = GRIDSIZE - 1;
		if (gridy < 0) gridy = 0;
		if (gridy > GRIDSIZE - 1) gridy = GRIDSIZE - 1;
		if (gridz < 0) gridz = 0;
		if (gridz > GRIDSIZE - 1) gridz = GRIDSIZE - 1;

		hashGrid[gridx][gridy][gridz].push_back(p);
		//printf("gridx : %d, gridy : %d\n", gridx, gridy);
	}
}

vector<Particle *> SPH::getNeighbor(int gridx, int gridy, int gridz, double radius, vector<Particle*> &mine)
{
	vector<Particle *>res;
	mine.clear();
#pragma omp parallel for collapse(4)
	for (int i = gridx - (int)radius; i <= gridx + (int)radius; i++)
	{
		for (int j = gridy - (int)radius; j <= gridy + (int)radius; j++)
		{
			for (int k = gridz - (int)radius; k <= gridz + (int)radius; k++)
			{
				if (i < 0 || i > GRIDSIZE - 1 || j < 0 || j > GRIDSIZE - 1 || k < 0 || k > GRIDSIZE - 1)
					continue;

				for (int l = 0; l < hashGrid[i][j][k].size(); l++)
				{
					res.push_back(hashGrid[i][j][k][l]);

					if (i == gridx && j == gridy && k == gridz)
						mine.push_back(hashGrid[i][j][k][l]);
				}
			}
		}
	}
	return res;
}
