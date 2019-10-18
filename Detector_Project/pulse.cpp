#include "pulse.h"
//======================================================================
//Constructors and destructor
//======================================================================
/*
	Default pulse object to quickly generate a few random particles.
	Initializes an array of size n.
*/
pulse::pulse()
{
	int n = 10;
	particles = new particle[n];
	capacity = n;
	type = "Proton";
	Mesh M(1, 1, 1);
	double V[] = { 4,1,1 };
	double P[] = { 0,0,0 };

	generateParticles(type, M, M, V, P);
}

/*
	Specialized pulse object generated with mesh M and random particles.
*/
pulse::pulse(int n, Mesh M)
{
	capacity = n;
	particles = new particle[n];
	this->type = "Random";

	default_random_engine generator;
	normal_distribution<double> distribution(5.0, 1.0); //velocity of 5 with sigma = 1.0

	for (int i = 0; i < n; i++)
	{
		// Choose a random particle type
		string type;

		int j = rand() % 4;
		if (j == 0)
		{
			type = "Muon";
		}
		else if (j == 1)
		{
			type = "Proton";
		}
		else if (j == 2)
		{
			type = "Neutron";
		}
		else if (j == 3)
		{
			type = "Pion";
		}

		//Generate velocity distribution
		double v[3];
		double p[3];
		for (int k = 0; k < 3; k++)
		{
			double mesh;
			if (k == 0) { mesh = M.x; }
			else if (k == 1) { mesh = M.y; }
			else if (k == 2) { mesh = M.z; }

			v[k] = distribution(generator);
			
			p[k] = distribution(generator) + mesh * i;
		}
		//add the particle to the particle array
		particles[i] = particle(type, 1, p, v); // (name, mass, pos, vel)
	}
}

/*
	Specialized pulse object generated with mesh M and a particle
	of some specified type with velocity v and position p.
*/
pulse::pulse(int n, Mesh Mp, Mesh Mv, string type, double v[], double p[])
{
	capacity = n;
	particles = new particle[n];
	this->type = type;
	generateParticles(type, Mp, Mv, v, p);
}

/*
	Very special pulse object used only to keep track of decay products.
*/
pulse::pulse(vector<particle> p)
{
	capacity = p.size();
	type = "Products";
	particles = new particle[capacity];

	for (int i = 0; i < capacity; i++)
	{
		particles[i] = p[i];
	}
}

pulse::~pulse()
{
	delete[] particles;
}

//======================================================================
// Pulse functions
//======================================================================

/*
	Print the pulse of the generated particles.
*/
void pulse::showPulse()
{
	cout << "Pulse type: " << type << endl;
	cout << "Pulse size: " << capacity << endl;
	cout << "=================================" << endl;

	for (int i = 0; i < capacity; i++)
	{
		int state = 0;
		if (particles[i].getCollision()) { state = 1; }
		if (particles[i].getParentDecay()) { state = 2; }
		if (particles[i].getProductDecay()) { state = 3; }
		cout << "Particle[" << i << "] :: State[" << state << "]" << endl;
		particles[i].particlePrint();
		cout << endl;
	}
}

/*
	Generates particles for the pulse using the given spacing M and
	the given velocities.
*/
void pulse::generateParticles(string type, Mesh Mp, Mesh Mv, double vel[], double pos[])
{
	default_random_engine generator;
	normal_distribution<double> distributionVX(vel[0], Mv.x); //velocity of vel[0] with sigma = 1.0
	normal_distribution<double> distributionVY(vel[1], Mv.y); //velocity of vel[1] with sigma = 1.0
	normal_distribution<double> distributionVZ(vel[2], Mv.z); //velocity of vel[2] with sigma = 1.0

	normal_distribution<double> distributionPX(pos[0], Mp.x); //position of pos[0] with sigma = 1.0
	normal_distribution<double> distributionPY(pos[1], Mp.y); //position of pos[1] with sigma = 1.0
	normal_distribution<double> distributionPZ(pos[2], Mp.z); //position of pos[2] with sigma = 1.0

	for (int i = 0; i < capacity; i++)
	{
		//Generate velocity distribution
		double v[3];
		double p[3];
		for (int k = 0; k < 3; k++)
		{
			normal_distribution<double> dP, dV;
			if (k == 0) { dP = distributionPX; dV = distributionVX; }
			else if (k == 1) { dP = distributionPY; dV = distributionVY;}
			else if (k == 2) { dP = distributionPZ; dV = distributionVZ;}

			v[k] = dV(generator);
			p[k] = dP(generator);
		}
		//add the particle to the particle array
		particles[i] = particle(type, 1, p, v); // (name, mass, pos, vel)
		double r;
		if (type == "Proton") { r = .008; }
		else if (type == "Muon") { r = .0001656; }
		else if (type == "Neutron") { r = .0008; }
		else if (type == "Pion") { r = .00049; }
		particles[i].setParticleRadius(r); //initialize with size 0.1 units
	}
}

//========================================
// Getters
//========================================
/*
	Returns a pointer to the array of particles generated.
*/
particle *pulse::getPulse()
{
	return particles;
}


/*
	Get the number of particles initialized
*/
int pulse::getSize()
{
	return capacity;
}

//========================================
//Setters
//========================================


//======================================================================
// Helper Functions
//======================================================================