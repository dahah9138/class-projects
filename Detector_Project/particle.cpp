#include "particle.h"

//========================================
//Constructors/destructor                |
//========================================
//default constructor
particle::particle()
{
	mass = 0;
	radius = 1; //by default radius of the particle is 1
	collisionCheck = false;
	decay_parent = false;
	decay_product = false;
	//initialize velocity and position arrays to zero
	for (int i = 0; i < 3; i++)
	{
		velocity[i] = 0.0;
		position[i] = 0.0;
	}
}

//specialized constructor
particle::particle(string name, double mass, double pos[], double vel[])
{
	this->name = name;
	this->mass = mass;
	radius = 1;
	collisionCheck = false;
	decay_product = false;
	decay_parent = false;
	for (int i = 0; i < 3; i++)
	{
		position[i] = pos[i];
		velocity[i] = vel[i];
	}
}

//destructor
//particle::~particle() {} //no dynamically allocated memory
//========================================

//========================================
//Getters                                |
//========================================
string particle::getParticleName() const
{
	return name;
}

double particle::getParticleMass() const
{
	return mass;
}

double particle::getParticleRadius() const
{
	return radius;
}

double *particle::getParticlePosition()
{
	return position;
}

double *particle::getParticleVelocity()
{
	return velocity;
}

bool particle::getCollision() const
{
	return collisionCheck;
}

bool particle::getParentDecay() const
{
	return decay_parent;
}

bool particle::getProductDecay() const
{
	return decay_product;
}
//========================================

//========================================
//Setters                                |
//========================================
void particle::setCollision(bool c)
{
	collisionCheck = c;
}

void particle::setParticleName(string name)
{
	this->name = name;
}

void particle::setParticleMass(double mass)
{
	this->mass = mass;
}

void particle::setParticleRadius(double r)
{
	this->radius = r;
}

//update the position of the particle
void particle::setParticlePosition(double pos[])
{
	for (int i = 0; i < 3; i++)
	{
		position[i] = pos[i];
	}
}

//update the velocity of the particle
void particle::setParticleVelocity(double vel[])
{
	for (int i = 0; i < 3; i++)
	{
		velocity[i] = vel[i];
	}
}

void particle::setProductDecay(bool d)
{
	decay_product = d;
}

void particle::setParentDecay(bool d)
{
	decay_parent = d;
}
//========================================

void particle::particlePrint()
{
	cout << "Particle name: " << name << endl;
	cout << "Particle mass: " << mass << endl;
	cout << "Particle position: [";
	for (int i = 0; i < 3; i++)
	{
		cout << position[i];
		if (i != 2)
		{
			cout << ", ";
		}
	}
	cout << "]" << endl;

	cout << "Particle velocity: [";
	for (int i = 0; i < 3; i++)
	{
		cout << velocity[i];
		if (i != 2)
		{
			cout << ", ";
		}
	}
	cout << "]" << endl;
}

void particle::updatePosition(double dt)
{
	for (int i = 0; i < 3; i++)
	{
		position[i] += velocity[i] * dt;
	}
}

/*
	Converts velocity to momentum.
*/
vector<double> particle::velocity_to_momentum(double m, vector<double> beta)
{
	constants cst;

	vector<double> p; //momentum of the particle
	
	beta_gamma bg = calculate_beta_gamma(beta[0], beta[1], beta[2]);

	double v = bg.beta;
	double p_norm = m * v * bg.gamma / cst.C; //divide because mass is in terms of MeV

	//Calculate components

	double theta = atan(beta[1] / beta[2]);
	double phi = acos(beta[2] / p_norm);

	double p0 = p_norm * cos(theta) * sin(phi);
	double p1 = p_norm * sin(theta) * sin(phi);
	double p2 = p_norm * cos(phi);

	p.push_back(p0);
	p.push_back(p1);
	p.push_back(p2);

	return p;
}



beta_gamma particle::calculate_beta_gamma(double bx, double by, double bz)
{
	double beta = sqrt(bx*bx + by * by + bz * bz);

	if (beta == 1) //massless particle
	{
		beta_gamma bg(beta, -1); //-1 indicates massless
		return bg;
	}
	else
	{
		double temp = 1 - beta * beta;
		if (temp < 0) { temp * -1; } //make sure temp is positive
		double gamma = 1.0 / sqrt(temp);
		beta_gamma bg(beta, gamma);
		return bg;
	}
}