#ifndef EVENTQUEUE_H
#define EVENTQUEUE_H
#include "particle.h"
#include "constants.h"
#include <vector>

struct particlePair
{
	particle p1, p2; //two particles to be stored
	string eventID; //id of event ~ want to have a generator to read these into the program
	double t_c; //collision time
	particlePair() { t_c = -1; }
	particlePair(particle p1, particle p2, string id, double t) { this->p1 = p1; this->p2 = p2; 
															 this->eventID = id; this->t_c = t; }
	void setID(string id) { eventID = id; }
	void setTC(double t) { t_c = t; }
};

class eventQueue
{
public:
	//===============================================
	//Constructor/destructor						|
	//===============================================
	eventQueue(int); //initializes heap with given size
	~eventQueue();
	//===============================================

	//===============================================
	//Heap functions								|
	//===============================================
	void heapify(int);
	particlePair dequeueNextEvent(); //extract the colliding particle
	particlePair peekCollision() const; //look at the next collision
	void dequeuePair(particlePair); //remove some particle pair
	void enqueuePair(particlePair); //add some particle pair
	void printQueue(); //print the queue of particles
	double calculateCollisionTime(particle p1, particle p2); //determine collision time for p1, p2
	string assignEventID(particlePair pair); //assign an ID to the event
	string checkID(particlePair pair); //create ID to be checked. Does not update class variables
	void particle_decay(double DT); //returns vector of event decay products
	bool isEmpty(); //true if empty, false otherwise
	int getSize(); //return the current size of the queue
	void queueStats(); //prints the different types of collisions seen in queue
	void clearProducts(); //empty products vector
	vector<particle> *getProducts(); //returns a pointer to the products
	//===============================================

private:
	//===============================================
	//Heap properties  								|
	//===============================================
	particlePair *heap; //array of particle pairs
	int max_size; //max size of the heap
	int current_size; //current size of the heap


	int nnEvents; //number of neutron-neutron events
	int npEvents; //number of neutron-proton events
	int mmEvents; //number of muon-muon events
	int ppEvents; //number of proton-proton events
	int pipiEvents; //number of pion-pion events

	int nmEvents; //neutron-muon
	int npiEvents; //neutron-pion
	int ppiEvents; //proton-pion
	int mpiEvents; //muon-pion
	int mpEvents; //muon-proton

	constants cst; //struct of constants used for calculations
	vector<particle> products; //products created from collisions
	//===============================================

	//===============================================
	//Helper functions								|
	//===============================================
	void swap(particlePair &pair1, particlePair &pair2); //swaps two events in the queue
	int parent(int index) const; //returns parent index of given index
	int leftChild(int index) const; //returns left child of given index
	int rightChild(int index) const; //returns right child of given index
	void decrementTime(int index, double t_new); //decrement the collision time of the event at index
	void deletePair(int); //delete pair at given index
	int search(string id); //search for event id in the heap
	void proton_downscatter(vector<particle>& p); //Generate downscattering of protons
	vector<double> momentum_to_velocity(double mass, double mom[]); //converts momentum to velocity
	//===============================================
};

#endif