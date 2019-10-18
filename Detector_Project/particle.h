#ifndef PARTICLE_H
#define PARTICLE_H
#define _USE_MATH_DEFINES
#include "constants.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
using namespace std;

/*
	This is a very basic struct that stores beta 
	and gamma to do calculations with.
*/
struct beta_gamma
{
	double beta;
	double gamma;
	beta_gamma() { beta = 0; gamma = 0; }
	beta_gamma(double b, double g) { beta = b; gamma = g; }
};

/*
	The class 'particle' is meant to keep current known velocity
	and momentum of particles. Universal features such as time
	steps and forces will be handled by higher structures. Recording
	the evolution of the particles will also be handled by a different
	struct.
*/
class particle
{
public:
	//Constructors/destructor
	//========================================
	particle(); //default constructor
	particle(string, double, double[], double[]); //specialized constructor: (name, mass, pos, vel)
	//~particle();
	//========================================
	//Getters
	//========================================
	string getParticleName() const;
	double getParticleMass() const;
	double getParticleRadius() const;
	double *getParticleVelocity(); //returns pointer to velocity array
	double *getParticlePosition(); //returns pointer to position array
	bool getCollision() const;
	bool getParentDecay() const;
	bool getProductDecay() const;
	//========================================
	//Setters
	//========================================
	void setParticleName(string);
	void setParticleMass(double);
	void setParticleRadius(double);
	void setParticleVelocity(double[]);
	void setParticlePosition(double[]);
	void setCollision(bool);
	void setParentDecay(bool);
	void setProductDecay(bool);
	//========================================

	//========================================
	//Particle functions
	//========================================
	void particlePrint(); //print the particle
	void updatePosition(double dt); //takes time step dt
	vector<double> velocity_to_momentum(double m, vector<double> beta); //mass, and vector of betas
	beta_gamma calculate_beta_gamma(double bx, double by, double bz); //returns beta_gamma object from betas
	//========================================

private:
	//Particle qualities
	//========================================
	string name;   //particle name
	double mass;   //particle mass
	bool collisionCheck; //true if added to eventQueue

	bool decay_product; //true if particle was created from decayed particle
	bool decay_parent; //true if particle decayed

	double radius; //particle radius
	double position[3]; //x,y,z
	double velocity[3]; //vx,vy,vz
	//========================================
};

#endif