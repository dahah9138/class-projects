#ifndef DETECTOR_H
#define DETECTOR_H
#include "eventQueue.h"
#include "pulse.h"
#include <vector>

/*
	This detector is an octree that stored particle objects.
	The choice to use an octree instead of a graph
	was made because when building the mesh of the graph, it
	took an incredibly large amount of computational power
	for space that may or may not have been used.

*/

struct Point3D //used to dimensionally represent non-particle objects
{
	double x, y, z;
	Point3D() { x = 0; y = 0; z = 0; }
	Point3D(double x, double y, double z) { this->x = x; this->y = y; this->z = z; }
};



class Detector
{
public:
	//========================================
	//Constructor/destructor
	//========================================
	Detector(Point3D c, Point3D halfDim); //center and half of side for the detector (halfDim.x *2, halfDim.y*2, halfDim.z*2)
	Detector(); //default is (0,0,0) and lengths (1,1,1)
	~Detector();
	//========================================
	//Tree functions
	//========================================
	int getOctant(Point3D p) const; //returns the octant of a point, not particle
	int getOctant(particle p) const; //returns the octant of a particle
	bool isLeafNode() const;
	void insertParticle(particle *p);
	void getPoints(Point3D& min, Point3D& max, vector<particle*> &discovered);
	void updateParticles(vector<particle*> arr, int size, double DT); //Update the particle stored in the tree with dt
	void updateParticles(vector<pulse*> arr, double DT); //Secondary update function for pulses
	void getParticlesInTree(vector<particle*> &stored); //Store all particles in the tree
	bool containedInTree(particle *p); //checks if particle p is stored in the tree
	void getCollisions(eventQueue &, double T); //Query that particles that have collided
	void getParticleRange(Point3D &pCenter, double r, vector<particle*> &discovered, bool firstCol); //modify vector of particles with all particles in range of center + r
	double getDist(Point3D, Point3D); //Returns distance between two points
	bool isParallel(particle *p1, particle *p2);
	void setDimensions(Point3D c, Point3D halfSize); //set the dimensions of the detector
private:
	//===============================================
	//Tree Properties
	//===============================================
	Point3D center;
	Point3D halfSize; // (l/2,w/2,h/2)
	Detector *children[8];
	particle *points; //stores particle data
	//===============================================
	//Helper Functions
	//===============================================
	void clear(); //deletes the detector
};

#endif