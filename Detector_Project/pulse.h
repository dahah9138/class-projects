#ifndef PULSE_H
#define PULSE_H
#include <random>
#include "particle.h"

struct Mesh //used to dimensionally spread out particles in pulse
{
	double x, y, z;
	Mesh() { x = 0; y = 0; z = 0; }
	Mesh(double x, double y, double z) { this->x = x; this->y = y; this->z = z; }
};


/*
	Pulse:
	Generates a pulse of particles that are shot at the detector.

*/

class pulse
{
public:
	pulse(); //default
	pulse(int n, Mesh M); //max number of particles initialized with spacing M
	pulse(int n, Mesh Mp, Mesh Mv, string type, double v[], double p[]); //number of particles to generate of certain type with spacing M
	pulse(vector<particle> p); //special constructor to pass decay products
	~pulse();
	//========================================
	//Getters
	//========================================
	particle *getPulse(); //returns a pointer to particles
	int getSize();
	//========================================
	//Setters
	//========================================

	//========================================
	//Functions
	//========================================
	void showPulse(); //show the particle information inside pulse
	void generateParticles(string type, Mesh Mp, Mesh Mv, double vel[], double pos[]); //generates particles with some type, mesh M, velocity vel, and position pos

private:
	particle *particles; //pointer to array of particle objects
	int capacity; //maximum size of the particles initialized
	string type; //type of the pulse
	//================================
	//Helper functions
	//================================

};

#endif