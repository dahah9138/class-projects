#include "eventQueue.h"

//===============================================================
//Helper functions												|
//===============================================================

/*
	This function swaps two particles in the queue
*/
void eventQueue::swap(particlePair &p1, particlePair &p2)
{
	particlePair pSwap = p1;
	p1 = p2;
	p2 = pSwap;
}

/*
	This function returns the parent of index
*/
int eventQueue::parent(int index) const
{
	return (index - 1) / 2;
}

/*
	This function returns the left child of index
*/
int eventQueue::leftChild(int index) const
{
	return (2 * index + 1);
}

/*
	This function returns the right child of index
*/
int eventQueue::rightChild(int index) const
{
	return (2 * index + 2);
}

/*
	This function decreases the collision time between two particles to t_new
*/
void eventQueue::decrementTime(int index, double t_new)
{
	heap[index].t_c = t_new;
	//repair upwards
	while (index != 0 && heap[parent(index)].t_c > heap[index].t_c)
	{
		swap(heap[index], heap[parent(index)]);
		index = parent(index);
	}
}

/*
	This function deletes the pair at the given index. This function
	assumes that the pair exists.
*/
void eventQueue::deletePair(int index)
{
	heap[index].t_c = INT_MIN;
	dequeueNextEvent();
}

/*
	Searches for the id in the heap. If it is unable to find the pair,
	the function returns an index of -1.
*/

int eventQueue::search(string id)
{
	for (int i = 0; i < current_size; i++)
	{
		if (id == heap[i].eventID)
		{
			return i;
		}
	}
	return -1; //could not find id
}
//===============================================================

//===============================================================
//Constructor/destructor										|
//===============================================================
eventQueue::eventQueue(int n)
{
	max_size = n;
	nnEvents = 0;
	mmEvents = 0;
	ppEvents = 0;
	npEvents = 0;
	pipiEvents = 0;

	nmEvents = 0;
	npiEvents = 0;
	ppiEvents = 0;
	mpiEvents = 0;
	mpEvents = 0;

	current_size = 0;
	heap = new particlePair[max_size];
}

eventQueue::~eventQueue()
{
	if (heap == 0)
	{
		return; //prevent errors
	}
	delete[] heap;
	heap = 0;
}
//===============================================================

//===============================================================
//Heap functions												|
//===============================================================
/*
	Returns the next event in the queue also removing it from the queue.
	Returns a default particlePair object if empty.
*/
particlePair eventQueue::dequeueNextEvent()
{
	if (current_size == 0)
	{
		cout << "Event queue is empty." << endl;
		particlePair defaultPair;
		return defaultPair;
	}
	else if (current_size == 1)
	{
		current_size = 0;
		return heap[0];
	}
	else
	{
		particlePair _pair = heap[0];
		heap[0] = heap[current_size - 1]; //set first to last
		current_size--;
		heapify(0); //call heapify on root
		return _pair;
	}
}

/*
	Returns the next particlePair collision.
*/
particlePair eventQueue::peekCollision() const
{
	if (current_size == 0)
	{
		particlePair p;
		return p; //return default p
	}
	//otherwise return the first value
	return heap[0];
}

/*
	Makes the heap a proper min heap
*/
void eventQueue::heapify(int i)
{
	int left = leftChild(i);
	int right = rightChild(i);
	int smallest = i;

	if (left < current_size && heap[left].t_c < heap[i].t_c)
	{
		smallest = left;
	}
	else if (right < current_size && heap[right].t_c < heap[smallest].t_c)
	{
		smallest = right;
	}

	if (smallest != i)
	{
		swap(heap[smallest], heap[i]);
		heapify(smallest);
	}
}

/*
	Dequeues the particlePair _pair from the array
*/
void eventQueue::dequeuePair(particlePair _pair)
{
	//search for the pair
	int i = search(_pair.eventID);
	if (i == -1)
	{	
		cout << "Could not find event id: " << _pair.eventID << endl;
		return; //could not find
	}
	else //found the pair
	{
		deletePair(i); //call deletePair helper function on index i
	}
}

/*
	This function enqueues a pair of pairs. It assumes that
	the time of the collisions has already been given.
*/
void eventQueue::enqueuePair(particlePair p)
{
	
	if (current_size == 0)
	{
		string ID = assignEventID(p);
		p.eventID = ID;
		heap[0] = p;
		current_size = 1;
		return;
	}
	if (search(checkID(p)) != -1)
	{
		//pair already exists
		cout << "Pair already exists." << endl;
		//want to avoid being here because it means double counting is occurring
		return;
	}
	if (max_size == current_size)
	{
		cout << "eventQueue full: Must increase size of heap." << endl;
		return; //prevent seg fault
	}
	
	string ID = assignEventID(p);
	p.eventID = ID;
	current_size++;
	int i = current_size - 1; //insert at end
	heap[i] = p;

	while (i != 0 && heap[parent(i)].t_c > heap[i].t_c)
	{
		swap(heap[i], heap[parent(i)]);
		i = parent(i);
	}
}

/*
	Prints the queue of particle collision pairs and their estimated
	collision times.
*/
void eventQueue::printQueue()
{
	if (current_size == 0)
	{
		cout << "Queue is empty." << endl;
		return;
	}

	//queue is not empty
	for (int i = 0; i < current_size; i++)
	{
		particlePair _pair = heap[i];
		cout << "[" << i << "]: Collision " << i + 1 << endl << "Pair: (";
		cout << _pair.p1.getParticleName() << ", " << _pair.p2.getParticleName() << ") " << endl;
		cout << "ID: " << _pair.eventID << endl;
		cout << "Estimated collision time: " << _pair.t_c << endl << endl;
	}
}

/*
	This function predicts the time at which the collision between
	p1 and p2 will take place. It is not necessary to use this function
	for the detector, because collisions are detected via particle radii
	in the detector.
*/
double eventQueue::calculateCollisionTime(particle p1, particle p2)
{
	//Use formula on whiteboard at home
	double t = 0;
	
	return 0;
}

/*
	This function assigns an ID to the event.

	As the simulation is able to support more
	particles, this function will increase in complexity.
*/
string eventQueue::assignEventID(particlePair pair)
{
	string IDname; //string to concatenate
	particle p1 = pair.p1;
	particle p2 = pair.p2;
	if (p1.getParticleName() == "Muon")
	{
		IDname += "m";
	}
	else if (p1.getParticleName() == "Neutron")
	{
		IDname += "n";
	}
	else if (p1.getParticleName() == "Proton")
	{
		IDname += "p";
	}
	else if (p1.getParticleName() == "Pion")
	{
		IDname += "pi";
	}

	if (p2.getParticleName() == "Muon")
	{
		IDname += "m";
	}
	else if (p2.getParticleName() == "Neutron")
	{
		IDname += "n";
	}
	else if (p2.getParticleName() == "Proton")
	{
		IDname += "p";
	}
	else if (p2.getParticleName() == "Pion")
	{
		IDname += "pi";
	}

	if (IDname == "mm")
	{
		IDname += to_string(mmEvents);
		mmEvents++; //increment number of muon-muon events
	}
	else if (IDname == "nn")
	{
		IDname += to_string(nnEvents);
		nnEvents++; //increment number of neutron-neutron events
	}
	else if (IDname == "pp")
	{
		IDname += to_string(ppEvents);
		ppEvents++; //increment number of proton-proton events
	}
	else if (IDname == "pipi")
	{
		IDname += to_string(pipiEvents);
		pipiEvents++; //increment number of pion-pion events
	}
	//account for mixtures of particles
	else if (IDname == "np" || IDname == "pn")
	{
		IDname = "np" + to_string(npEvents);
		npEvents++; //increment number of neutron-proton events
	}
	else if (IDname == "mp" || IDname == "pm")
	{
		IDname = "mp" + to_string(mpEvents);
		mpEvents++;
	}
	else if (IDname == "mpi" || IDname == "pim")
	{
		IDname = "mpi" + to_string(mpiEvents);
		mpiEvents++;
	}
	else if (IDname == "ppi" || IDname == "pip")
	{
		IDname = "pip" + to_string(ppiEvents);
		ppiEvents++;
	}
	else if (IDname == "npi" || IDname == "pin")
	{
		IDname = "npi" + to_string(npiEvents);
		npiEvents++;
	}
	else if (IDname == "nm" || IDname == "mn")
	{
		IDname = "nm" + to_string(nmEvents);
		nmEvents++;
	}
	return IDname;
}

/*
	This function behaves similarly to assignEventID
	but does not update the number of collision events
	for a certain collision pair type.
*/
string eventQueue::checkID(particlePair pair)
{
	string IDname; //string to concatenate
	particle p1 = pair.p1;
	particle p2 = pair.p2;
	if (p1.getParticleName() == "Muon")
	{
		IDname += "m";
	}
	else if (p1.getParticleName() == "Neutron")
	{
		IDname += "n";
	}
	else if (p1.getParticleName() == "Proton")
	{
		IDname += "p";
	}
	else if (p1.getParticleName() == "Pion")
	{
		IDname += "pi";
	}

	if (p2.getParticleName() == "Muon")
	{
		IDname += "m";
	}
	else if (p2.getParticleName() == "Neutron")
	{
		IDname += "n";
	}
	else if (p2.getParticleName() == "Proton")
	{
		IDname += "p";
	}
	else if (p2.getParticleName() == "Pion")
	{
		IDname += "pi";
	}

	if (IDname == "mm")
	{
		IDname += to_string(mmEvents);
	}
	else if (IDname == "nn")
	{
		IDname += to_string(nnEvents);
	}
	else if (IDname == "pp")
	{
		IDname += to_string(ppEvents);
	}
	else if (IDname == "pipi")
	{
		IDname += to_string(pipiEvents);
	}
	//account for mixtures of particles
	else if (IDname == "np" || IDname == "pn")
	{
		IDname = "np" + to_string(npEvents);
	}
	else if (IDname == "mp" || IDname == "pm")
	{
		IDname = "mp" + to_string(mpEvents);
	}
	else if (IDname == "mpi" || IDname == "pim")
	{
		IDname = "mpi" + to_string(mpiEvents);
	}
	else if (IDname == "ppi" || IDname == "pip")
	{
		IDname = "pip" + to_string(ppiEvents);
	}
	else if (IDname == "npi" || IDname == "pin")
	{
		IDname = "npi" + to_string(npiEvents);
	}
	else if (IDname == "nm" || IDname == "mn")
	{
		IDname = "nm" + to_string(nmEvents);
	}
	return IDname;
}

/*
	This function iterates through all particles stored
	and pops them from the queue after getting decay products.

*/
void eventQueue::particle_decay(double DT)
{
	default_random_engine generator;
	uniform_real_distribution<double> distribution(0, 1);

	while (current_size != 0)
	{
		// Dequeue the event
		particlePair pair = dequeueNextEvent();
		particle p[2];
		p[0] = pair.p1;
		p[1] = pair.p2;
		
		//Calculate decay for each particle
		for (int i = 0; i < 2; i++)
		{
			string name = p[i].getParticleName();
			double decay_prob = 0; //probability of decay
			if (name == "Proton")
			{
				//Proton decay func
				
				double threshold = cst.m_proton + 2 * cst.m_muon; // E = (m_proton + 2 * m_muon)c^2
				

				//Check if proton meets energy threshold
				double *vel = p[i].getParticleVelocity();
				vector<double> betas;
				betas.push_back(vel[0] / cst.C);
				betas.push_back(vel[1] / cst.C);
				betas.push_back(vel[2] / cst.C);
				
				vector<double> p_mom = p[i].velocity_to_momentum(cst.m_proton, betas);
				double E = sqrt(p_mom[0] * p_mom[0] + p_mom[1] * p_mom[1] + p_mom[2] * p_mom[2])*cst.C;
				E *= E;
				E += cst.m_proton * cst.m_proton;
				E = sqrt(E);
				//cout << "energy: " << E << endl;

				if (E > threshold)
				{
					// Interaction occurs
					beta_gamma bg = p[i].calculate_beta_gamma(betas[0], betas[1], betas[2]);
					double time_const = 1 / bg.beta * 7.006 * 10e-6;
					decay_prob = 1 - exp(-DT/time_const);
					//cout << "decay prob: " << decay_prob << endl;
				}
			}
			else if (name == "Pion")
			{
				double *vel = p[i].getParticleVelocity();
				vector<double> betas;
				betas.push_back(vel[0] / cst.C);
				betas.push_back(vel[1] / cst.C);
				betas.push_back(vel[2] / cst.C);

				beta_gamma bg = p[i].calculate_beta_gamma(betas[0], betas[1], betas[2]);
				double time_const = bg.gamma * 2.603e-8;
				//cout << "time const: " << time_const << endl;
				decay_prob = 1 - exp(-DT / time_const);
			}
			else if (name == "Muon")
			{
				double *vel = p[i].getParticleVelocity();
				vector<double> betas;
				betas.push_back(vel[0] / cst.C);
				betas.push_back(vel[1] / cst.C);
				betas.push_back(vel[2] / cst.C);

				beta_gamma bg = p[i].calculate_beta_gamma(betas[0], betas[1], betas[2]);
				double time_const = bg.gamma * 2.197e-6;
				decay_prob = 1 - exp(-DT / time_const);
			}

			//Check is particle decayed

			//Generate random double between 0 and 1

			double num = distribution(generator);

			//Compare if (num < decay_prob) particle decayed; else append particle to products
			//cout << "num: " << num << endl;
			if (num < decay_prob)
			{
				//particle decayed
				p[i].setParentDecay(true);

				if (p[i].getParticleName() == "Proton")
				{
					//Call proton down scatter
					vector<particle>scatter;
					scatter.push_back(p[i]);
					proton_downscatter(scatter);
					//products.push_back(/*particles from downscatter*/); // iterate though all particles
					for (vector<particle>::iterator l = scatter.begin(); l != scatter.end(); l++)
					{
						products.push_back(*l);
					}
				}
				else if (name == "Pion")
				{
					p[i].setParentDecay(true);

					//Generate a muon and a neutrino
					double theta = 2 * M_PI * distribution(generator);
					double phi = 2 * M_PI * distribution(generator);
					double mom = (cst.m_pion * cst.m_pion - cst.m_muon * cst.m_muon) / cst.C / (2 * cst.m_pion);

					double muon_x = mom * cos(theta) * sin(phi);
					double muon_y = mom * sin(theta) * sin(phi);
					double muon_z = mom * cos(phi);

					double muon_p[] = { muon_x, muon_y, muon_z };
					double muon_v[] = { cst.C * cst.C * muon_x / cst.m_muon, cst.C * cst.C * muon_y / cst.m_muon, cst.C * cst.C * muon_z / cst.m_muon };

					particle p1("Muon", cst.m_muon, p[i].getParticlePosition(), muon_v);
					
					double neutrino_v[] = { -muon_v[0], -muon_v[1], -muon_v[2] };
					particle p2("Neutrino", 0, p[i].getParticlePosition(), neutrino_v);

					products.push_back(p1);
					products.push_back(p2);
				}
				else if (name == "Muon")
				{
					p[i].setParentDecay(true);

					//Generates an electron and a neutrino
					double theta = 2 * M_PI * distribution(generator);
					double phi = 2 * M_PI * distribution(generator);
					double mom = (cst.m_electron * cst.m_electron - cst.m_muon * cst.m_muon) / cst.C / (2 * cst.m_electron);

					double electron_x = mom * cos(theta) * sin(phi);
					double electron_y = mom * sin(theta) * sin(phi);
					double electron_z = mom * cos(phi);

					double electron_p[] = { electron_x, electron_y, electron_z };
					double electron_v[] = { cst.C * cst.C * electron_x / cst.m_electron, cst.C * cst.C * electron_y / cst.m_electron, cst.C * cst.C * electron_z / cst.m_electron };

					particle p1("Electron", cst.m_electron, p[i].getParticlePosition(), electron_v);

					double neutrino_v[] = { -electron_v[0], -electron_v[1], -electron_v[2] };
					particle p2("Neutrino", 0, p[i].getParticlePosition(), neutrino_v);

					products.push_back(p1);
					products.push_back(p2);
				}
				
			}
			else
			{
				//No interaction
				p[i].setCollision(false);
				//products.push_back(p[i]);
			}
		}
		
	}
}

void eventQueue::clearProducts()
{
	products.clear();
}

bool eventQueue::isEmpty()
{
	if (current_size == 0) { return true; }
	else { return false; }
}

int eventQueue::getSize()
{
	return current_size;
}


void eventQueue::queueStats()
{
	cout << "(Proton, Proton): " << ppEvents << endl;
	cout << "(Proton, Neutron): " << npEvents << endl;
	cout << "(Proton, Pion): " << ppiEvents << endl;
	cout << "(Neutron, Neutron): " << nnEvents << endl;
	cout << "(Pion, Pion): " << pipiEvents << endl;
	cout << "(Neutron, Muon): " << nmEvents << endl;
	cout << "(Muon, Proton): " << mpEvents << endl;
	cout << "(Muon, Pion): " << mpiEvents << endl;
	cout << "(Muon, Muon): " << mmEvents << endl;
	cout << "(Neutron, Pion): " << npiEvents << endl;
}

vector<particle> *eventQueue::getProducts()
{
	return &products;
}

/*
	This function generates two pions and the initial
	proton that down-scattered.

	The downscattered proton is stored in the vector being passed.
*/
void eventQueue::proton_downscatter(vector<particle>&p)
{
	double b = 500; //MeV

	//Generate random numbers
	default_random_engine generator;
	uniform_real_distribution<double> distribution(0, 1);

	double r1 = distribution(generator);
	double r2 = distribution(generator);

	double randAngle1 = 2 * M_PI * distribution(generator);
	double randAngle2 = 2 * M_PI * distribution(generator);

	//Calculate energies
	double E_pi1 = cst.m_pion - (1 / b)*log(1 - r1);
	double E_pi2 = cst.m_pion - (1 / b)*log(1 - r2);

	//Calculate momentum
	double pi_1 = 1 / cst.C * sqrt(E_pi1 * E_pi1 - cst.m_pion * cst.m_pion);
	double pi_2 = 1 / cst.C * sqrt(E_pi2 * E_pi2 - cst.m_pion * cst.m_pion);
	
	//Get components
	double pi_1_x = pi_1 * cos(randAngle1) * sin(randAngle2);
	double pi_1_y = pi_1 * sin(randAngle1) * sin(randAngle2);
	double pi_1_z = pi_1 * cos(randAngle2);

	double pi_2_x = pi_2 * cos(randAngle2) * sin(randAngle1);
	double pi_2_y = pi_2 * sin(randAngle2) * sin(randAngle1);
	double pi_2_z = pi_2 * cos(randAngle1);

	double proton_x = -pi_1_x - pi_2_x;
	double proton_y = -pi_1_y - pi_2_y;
	double proton_z = -pi_1_z - pi_2_z;

	double mom1[] = { pi_1_x, pi_1_y, pi_1_z };
	double mom2[] = { pi_2_x, pi_2_y, pi_2_z };
	double mom3[] = { proton_x, proton_y, proton_z };

	vector<double> velocity1 = momentum_to_velocity(cst.m_pion, mom1);
	vector<double> velocity2 = momentum_to_velocity(cst.m_pion, mom2);
	vector<double> velocity3 = momentum_to_velocity(cst.m_proton, mom3);

	double vel1[] = { velocity1[0], velocity1[1], velocity1[2] };
	double vel2[] = { velocity2[0], velocity2[1], velocity2[2] };
	double vel3[] = { velocity3[0], velocity3[1], velocity3[2] };

	particle pi1("Pion", cst.m_pion, p[0].getParticlePosition(), vel1);
	particle pi2("Pion", cst.m_pion, p[0].getParticlePosition(), vel2);
	particle proton("Proton", cst.m_proton, p[0].getParticlePosition(), vel3);

	//insert the particles
	p.empty();
	p.push_back(proton);
	p.push_back(pi1);
	p.push_back(pi2);
}


vector<double> eventQueue::momentum_to_velocity(double mass, double mom[])
{
	vector<double> velocity;
	if (mass == cst.m_neutrino)
	{
		double theta = atan(mom[1] / mom[0]);
		double phi = atan(sqrt(mom[0] * mom[0] + mom[1] * mom[1]) / mom[2]);
		double beta_x = cos(theta) * sin(phi);
		double beta_y = sin(theta) * sin(phi);
		double beta_z = cos(phi);
		velocity.push_back(beta_x);
		velocity.push_back(beta_y);
		velocity.push_back(beta_z);
	}
	else
	{
		double theta = atan(mom[1] / mom[0]);
		double phi = atan(sqrt(mom[0] * mom[0] + mom[1] * mom[1]) / mom[2]);
		double p = sqrt(mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2]);
		double beta = (cst.C)*(p / mass) / sqrt(1 + (p / mass / cst.C));
		double beta_x = beta * cos(theta) * sin(phi);
		double beta_y = beta * sin(theta) * sin(phi);
		double beta_z = beta * cos(phi);
		velocity.push_back(beta_x);
		velocity.push_back(beta_y);
		velocity.push_back(beta_z);
	}
	return velocity;
}
//===============================================================