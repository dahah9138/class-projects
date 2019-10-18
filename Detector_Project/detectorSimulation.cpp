// detectorSimulation.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <iostream>
#include <stdlib.h>
#include "detector.h"
#include "pulse.h"

//Time step and number of steps used for calculations
//======================================================
#define STEPS 500 //number of steps
#define DT .000000005 //delta t
//======================================================



/*
	Prints the particles currently stored in the detector.
	Mainly used to debug if the detector is acting strange.
*/
void print(Detector &detector)
{
	vector<particle*> query; //query the detector for particles within box

	detector.getParticlesInTree(query);
	cout << endl << "================================" << endl;
	cout << "Particle Query:";
	cout << endl << "================================" << endl;

	vector<particle*>::iterator i;
	for (i = query.begin(); i != query.end(); i++)
	{
		(*i)->particlePrint();
	}
}

vector<particle> update(Detector &detector, eventQueue &queue, vector<particle*> particles, int size, double &T, vector<double>& decayTime)
{
	vector<particle> products;
	int prevSize = queue.getSize();
	detector.updateParticles(particles, size, DT);
	detector.getCollisions(queue, T); //get the collisions at time T
	if (prevSize != queue.getSize())
	{
		cout << endl << "================================" << endl;
		cout << "Particles collisions stored to queue:";
		cout << endl << "================================" << endl;

		//print the eventQueue
		queue.printQueue();

		queue.particle_decay(DT);
		products = *queue.getProducts();
		for (size_t i = 0; i < products.size(); i++)
		{
			decayTime.push_back(T);
		}
	}

	//Return decay products to be appended to total decay products
	
	return products;
}

vector<particle> testDetector(vector<double> &decayTime, string particle1, string particle2, int N)
{
	//===========================================================
	//Testing the detector implementation
	//===========================================================
	double T = 0; //initialize time at 0
	Point3D center(0, 0, 0); //coordinates of the center of the detector
	Point3D halfSize(500, 500, 500); //1000x1000x1000 detector
	Detector detector(center, halfSize);

	//const int N = 2000;

	Mesh Mp(2, 1, 1);
	Mesh Mv(.01, .0001, .0001);

	constants cst;

	double vel[] = { cst.C * .85, 0, 0};
	double pos[] = { -300,0,0 };

	double vel2[] = { -cst.C * .85,0,0 };
	double pos2[] = { 300,0,0 };

	pulse pulse1(N, Mp, Mv, particle1, vel, pos);
	pulse pulse2(N, Mp, Mv, particle2, vel2, pos2);
	particle *p = pulse1.getPulse();
	particle *p2 = pulse2.getPulse();

	vector<particle*> merged;
	merged.push_back(p);
	merged.push_back(p2);
	vector<particle*>::iterator vec;
	for (vec = merged.begin(); vec != merged.end(); vec++)
	{
		for (int k = 0; k < N; k++)
		{
			if (detector.containedInTree(&((*vec)[k])))
			{
				detector.insertParticle(&((*vec)[k]));
			}
		}
	}
	

	//create the queue
	eventQueue queue(100);

	vector<particle> totalProducts;
	
	//update 'STEPS' times
	for (int i = 0; i < STEPS; i++)
	{
		cout << "(" << i + 1 << "/" << STEPS << ") updating..." << endl;
		T += DT;
		
		vector<particle> temp = update(detector, queue, merged, N, T, decayTime);
		
		if (temp.size() != 0)
		{
			//append new products to total products
			for (vector<particle>::iterator m = temp.begin(); m != temp.end(); m++)
			{
				totalProducts.push_back(*m);
			}
		}
		//clear products stored in queue (no double counting)
		queue.clearProducts();
	}
	queue.queueStats();
	cout << "totalProducts: " << totalProducts.size() << endl;
	//Now process all the products appended
	return totalProducts;
}



/*
	Prints the decay products
*/
void printProducts(vector<double> times, vector<particle> particles)
{
	int count = 0;
	for (vector<particle>::iterator i = particles.begin(); i != particles.end(); i++)
	{

		cout << "Decay time: " << times[count++] << endl;
		i->particlePrint();
		cout << endl;
	}
}


/*
	Menu of options displayed to user.
*/

void menu()
{
	cout << "==================================" << endl;
	cout << ">> 1. Detect particles" << endl;
	cout << ">> 2. Print products" << endl;
	cout << ">> 3. Write to file" << endl;
	cout << ">> 4. Quit" << endl;
}

/*
	Checks if particles passed are valid particles
*/
bool checkValid(string p1, string p2)
{
	bool p1Check = false;
	bool p2Check = false;
	if (p1 == "Muon" || p1 == "Proton" || p1 == "Pion")
	{
		p1Check = true;
	}
	if (p1 == "Muon" || p1 == "Proton" || p1 == "Pion")
	{
		p2Check = true;
	}

	if (p1Check && p2Check)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/*
	This function writes the particles to the file passed.
*/
void writeToFile(string file, vector<particle> particles, vector<double> times)
{
	int count = 0;
	ofstream out_file(file);
	if (out_file.is_open())
	{
		for (vector<particle>::iterator i = particles.begin(); i != particles.end(); i++)
		{
			//write each individual particle
			double *pos = i->getParticlePosition();
			double *vel = i->getParticleVelocity();
			string name = i->getParticleName();
			out_file << name << endl;
			out_file << "Decay time: " << times[count++] << endl;
			out_file << "Velocity: (";
			for (int j = 0; j < 3; j++)
			{
				out_file << vel[j];
				if (j != 2)
				{
					out_file << ", ";
				}
			}
			out_file << ")" << endl;

			out_file << "Position: (";
			for (int j = 0; j < 3; j++)
			{
				out_file << pos[j];
				if (j != 2)
				{
					out_file << ", ";
				}
			}
			out_file << ")" << endl;
		}
	}
	else
	{
		cout << "Error opening file" << endl;
	}
	out_file.close();
}



/*
	If there is an overflow when initializing particles it means that the compiler was not
	able to secure enough memory in order to initialize that many particles. If there is an
	overflow, it is recommended to increase the size of the detector in testDetector() and
	increase the values in Mp (that is the sigma used to distribute particles, a higher sigma
	distributes the particles over a larger space).

	As far as I can tell the overflow occurs because the octree subdivides too many times
	causing an error.
*/


int main(int argc, char * argv[])
{
	vector <double> times;
	vector<particle> particles; 
	
	int option = 0;
	bool go = true;

	//Menu variables
	string particle1;
	string particle2;
	string file;
	int N = 0;
	
	constants cst; //struct of constants used


	//Continue displaying menu until user quits
	while (go)
	{
		menu(); //Display the menu
		cin >> option;

		switch(option)
		{
		case 1:
			while (!checkValid(particle1, particle2))
			{
				particle1 = "";
				particle2 = "";
				cout << "Enter the first particle type:" << endl;
				cin >> particle1;
				cout << "Enter the second particle type:" << endl;
				cin >> particle2;
			}
			while (N == 0)
			{
				cout << "How many particles per pulse? " << endl;
				cin >> N;
				if (N < 10)
				{
					cout << "Please enter a larger number." << endl;
					N = 0;
				}
			}

			//initialize detector with desired particles
			particles = testDetector(times, particle1, particle2, N);
			particle1 = "";
			particle2 = ""; //set back to empty to avoid false negative
			N = 0; //set N back to zero
			break;
		case 2:
			cout << endl;
			if (particles.size() != 0)
			{
				printProducts(times, particles);
			}
			else
			{
				cout << "No events occurred" << endl;
			}
			break;
		case 3:
			cout << "Enter file name (no spaces):" << endl;
			cin >> file;
			//Once called, clear particles and times
			writeToFile(file, particles, times);
			particles.clear();
			times.clear();
			cout << "Write to file complete." << endl;
			break;
		case 4:
			cout << "Goodbye!" << endl;
			go = false;
			break;
		}
	}

	return 0;
}