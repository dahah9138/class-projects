#include "detector.h"

/*
	This constructor initializes the detector with
	center c and dimensions twice that of halfDim.
*/

Detector::Detector(Point3D c, Point3D halfDim)
{
	center = c;
	halfSize = halfDim;
	points = 0;
	for (int i = 0; i < 8; i++)
	{
		children[i] = 0;
	}
}

Detector::Detector()
{
	Point3D c(0, 0, 0);
	Point3D l(1, 1, 1);
	center = c;
	halfSize = l;
	for (int i = 0; i < 8; i++)
	{
		children[i] = 0;
	}
}

/*
	This destructor destroys the children. I have
	not yet tested if there is leaked memory.
*/

Detector::~Detector()
{
	clear();
}

/*
	This function takes a point and returns the
	octant that it belongs to.
*/

int Detector::getOctant(Point3D p) const
{
	int octant;
	if (p.x < center.x && p.y < center.y && p.z < center.z) { octant = 0; }
	if (p.x < center.x && p.y < center.y && p.z >= center.z) { octant = 1; }
	if (p.x < center.x && p.y >= center.y && p.z < center.z) { octant = 2; }
	if (p.x < center.x && p.y >= center.y && p.z >= center.z) { octant = 3; }
	if (p.x >= center.x && p.y < center.y && p.z < center.z) { octant = 4; }
	if (p.x >= center.x && p.y < center.y && p.z >= center.z) { octant = 5; }
	if (p.x >= center.x && p.y >= center.y && p.z < center.z) { octant = 6; }
	if (p.x >= center.x && p.y >= center.y && p.z >= center.z) { octant = 7; }
	return octant;
}

/*
	This function takes a particle and returns the octant
	that the particle belongs to.
*/

int Detector::getOctant(particle p) const
{
	double *pos = p.getParticlePosition();
	Point3D point(pos[0], pos[1], pos[2]); //convert particle position to point
	return getOctant(point); //Calls getOctant for point
}

/*
	Checks if the Octree is a leaf or not.
	If the first child is NULL then the
	node has not divided and the current Octree
	must be a leaf.
*/

bool Detector::isLeafNode() const
{
	return children[0] == NULL;
}

void Detector::setDimensions(Point3D p, Point3D l)
{
	center = p;
	halfSize = l;
}

/*
	This functions inserts a particle into the Octree.
	If there is already a particle present, then the tree
	subdivides until there is only one particle per square.
	Divided nodes do not store children, only leafs.

	***This function does not prevent particles not stored
	in the tree from being entered into the tree. Make sure
	to check if the particle is within the boundaries of the
	tree before adding it.

	When adding particles to the tree, make sure to add them as
	new objects or insert them with a kept array (preserve the array).
*/

void Detector::insertParticle(particle *p)
{
	
	if (isLeafNode())
	{
		if (points == 0) //first data
		{
			points = p;
			return;
		}
		else //not first data
		{
			particle *prevPoint = points;
			points = 0;
			//subdivide and store points in smaller node
			for (int i = 0; i < 8; i++)
			{
				Point3D origin = center;
				origin.x += halfSize.x * (i & 4 ? .5f : -.5f);
				origin.y += halfSize.y * (i & 2 ? .5f : -.5f);
				origin.z += halfSize.z * (i & 1 ? .5f : -.5f);
				Point3D newHalf(halfSize.x *.5f, halfSize.y *.5f, halfSize.z *.5f);
				children[i] = new Detector(origin, newHalf);
			}

			//insert points
			int octPrev = getOctant(*prevPoint);
			int octP = getOctant(*p);

			children[octPrev]->insertParticle(prevPoint);
			children[octP]->insertParticle(p);
		}
	}
	else
	{
		//This node has children
		children[getOctant(*p)]->insertParticle(p);
	}
}

void Detector::getPoints(Point3D& min, Point3D& max, vector<particle*> &discovered)
{	
	if (isLeafNode())
	{
		if (points != 0)
		{
			//there is a point
			double *pos = points->getParticlePosition();
			Point3D ppoints(pos[0], pos[1], pos[2]);

			if (ppoints.x > max.x || ppoints.y > max.y || ppoints.z > max.z)
			{
				return; //out of bounds
			}
			if (ppoints.x < min.x || ppoints.y < min.y || ppoints.z < min.z)
			{
				return; //out of bounds
			}
			discovered.push_back(points); //store the actual particle
		}
	}
	else
	{
		//has children
		for (int i = 0; i < 8; i++)
		{
			//find max and min of children
			Point3D childMax, childMin;
			childMax.x = children[i]->center.x + children[i]->halfSize.x;
			childMax.y = children[i]->center.y + children[i]->halfSize.y;
			childMax.z = children[i]->center.z + children[i]->halfSize.z;

			childMin.x = children[i]->center.x - children[i]->halfSize.x;
			childMin.y = children[i]->center.y - children[i]->halfSize.y;
			childMin.z = children[i]->center.z - children[i]->halfSize.z;

			if (childMax.x < min.x || childMax.y < min.y || childMax.z < min.z)
			{
				continue;
			}
			if (childMin.x > max.x || childMin.y > max.y || childMin.z > max.z)
			{
				continue;
			}
			//The child might have a point within the box
			children[i]->getPoints(min, max, discovered);
		}
	}
}

/*
	This function traverses the tree and pushes all particles
	in the tree to the array. 
*/
void Detector::getParticlesInTree(vector<particle*> &arr)
{
	// These are the boundaries of the detector
	Point3D max(center.x + halfSize.x, center.y + halfSize.y, center.z + halfSize.z);
	Point3D min(center.x - halfSize.x, center.y - halfSize.y, center.z - halfSize.z);

	if (isLeafNode())
	{
		if (points != 0)
		{
			//there is a point
			double *pos = points->getParticlePosition();
			Point3D ppoints(pos[0], pos[1], pos[2]);

			if (ppoints.x > max.x || ppoints.y > max.y || ppoints.z > max.z)
			{
				return; //out of bounds
			}
			if (ppoints.x < min.x || ppoints.y < min.y || ppoints.z < min.z)
			{
				return; //out of bounds
			}
			
			arr.push_back(points); //particle did not interact yet
		}
	}
	else
	{
		//has children
		for (int i = 0; i < 8; i++)
		{
			//find max and min of children
			Point3D childMax, childMin;
			childMax.x = children[i]->center.x + children[i]->halfSize.x;
			childMax.y = children[i]->center.y + children[i]->halfSize.y;
			childMax.z = children[i]->center.z + children[i]->halfSize.z;

			childMin.x = children[i]->center.x - children[i]->halfSize.x;
			childMin.y = children[i]->center.y - children[i]->halfSize.y;
			childMin.z = children[i]->center.z - children[i]->halfSize.z;

			if (childMax.x < min.x || childMax.y < min.y || childMax.z < min.z)
			{
				continue;
			}
			if (childMin.x > max.x || childMin.y > max.y || childMin.z > max.z)
			{
				continue;
			}
			//The child might have a point within the box
			children[i]->getPoints(min, max, arr); //Don't recalculate min and max
		}
	}
}

/*
	This function updates the positions of the particles
	but not the velocities. Every time the particles are
	updated the tree is rebuilt.

*/
void Detector::updateParticles(vector<particle*> arr, int size, double DT)
{
	// Clear the detector
	clear();
	// Insert the particles back into the cleared tree
	vector<particle*>::iterator i;
	for (i = arr.begin(); i != arr.end(); i++)
	{
		particle *temp = *i;
		for (int j = 0; j < size; j++)
		{
			particle *p = &temp[j];
			//modify p
			p->updatePosition(DT);

			//check if particle should still be stored
			if (containedInTree(p) && !p->getCollision() && !p->getParentDecay())
			{
				insertParticle(p); //Reinsert the particle
			}
		}
	}
}


/*
	This function updates the positions of the particles
	but not the velocities. Every time the particles are
	updated the tree is rebuilt.

	Instead of taking a vector of particle pointers, this
	function takes a vector of pulse pointers.
*/
void Detector::updateParticles(vector<pulse*> arr, double DT)
{
	// Clear the detector
	//clear();
	// Insert the particles back into the cleared tree
	vector<pulse*>::iterator i;
	for (i = arr.begin(); i != arr.end(); i++)
	{
		pulse *temp = *i;
		particle *ptemp = temp->getPulse();

		for (int j = 0; j < temp->getSize(); j++)
		{
			particle *p = &ptemp[j];
			//modify p
			p->updatePosition(DT);

			//check if particle should still be stored
			if (containedInTree(p) && !p->getCollision() && !p->getParentDecay())
			{
				insertParticle(p); //Reinsert the particle
			}
		}
	}
}


void Detector::clear()
{
	for (int i = 0; i < 8; i++)
	{
		delete children[i];
		children[i] = 0;
	}
}

/*
	This function returns true or false
	depending on whether the particle is within
	the boundaries of the detector.
*/
bool Detector::containedInTree(particle *p)
{
	// These are the boundaries of the detector
	Point3D max(center.x + halfSize.x, center.y + halfSize.y, center.z + halfSize.z);
	Point3D min(center.x - halfSize.x, center.y - halfSize.y, center.z - halfSize.z);

	if (p != 0)
	{
		//there is a point
		double *pos = p->getParticlePosition();
		Point3D ppoints(pos[0], pos[1], pos[2]);

		if (ppoints.x > max.x || ppoints.y > max.y || ppoints.z > max.z)
		{
			return false; //out of bounds
		}
		if (ppoints.x < min.x || ppoints.y < min.y || ppoints.z < min.z)
		{
			return false; //out of bounds
		}
		return true;
	}
	return false; //The particle was NULL
}

/*
	This function takes an event queue, and
	pushes all particles that collided into the 
	queue

*/
void Detector::getCollisions(eventQueue &queue, double T)
{
	// Use getPoints to search for collisions

	// Iterate through particles in the tree

	vector<particle*> particles; //particles in tree

	//store particles of the tree to a vector
	getParticlesInTree(particles);

	//iterate through each particle
	vector<particle*>::iterator i;
	for (i = particles.begin(); i != particles.end(); i++)
	{
		if (!(*i)->getCollision())
		{
			//query the tree of collisions for each particle
			//create boundaries to search with
			particle *p = *i;
			double r = .016; // 2 * proton radius
			double *pos = p->getParticlePosition();
			Point3D min(pos[0] - r, pos[1] - r, pos[2] - r);
			Point3D max(pos[0] + r, pos[1] + r, pos[2] + r);

			vector<particle*> partner;

			//call getPoints
			getPoints(min, max, partner);

			vector<particle*>::iterator j;
			for (j = partner.begin(); j != partner.end(); j++)
			{
				if (*j != *i && !((*j)->getCollision()) && isParallel((*i), (*j)))
				{
					(*j)->setCollision(true);
					(*i)->setCollision(true);
					particlePair pair(**i, **j, "", T);
					//cout << "(" << (*i)->getParticleName() << ", " << (*j)->getParticleName() << ")" << endl;
					queue.enqueuePair(pair);
				}
			}
		}
	}
}

/*
	This function takes an empty vector of particle objects, and pushes all
	found particles within given spherical range to vector.
*/
void Detector::getParticleRange(Point3D &pCenter, double r, vector<particle*> &discovered, bool firstCol)
{
	// These are the boundaries of the particle passed
	Point3D max(pCenter.x + r, pCenter.y + r, pCenter.z + r);
	Point3D min(pCenter.x - r, pCenter.y - r, pCenter.z - r);

	if (isLeafNode())
	{
		if (points != 0)
		{
			//there is a point
			double *pos = points->getParticlePosition();
			Point3D ppoints(pos[0], pos[1], pos[2]);

			
			//compare squared distance between centers and radii
			double d = getDist(ppoints, pCenter); //if d = 0, it is the same particle
			double rr = points->getParticleRadius() + r;
			if (d > rr*rr || d == 0)
			{
				return; //particle is not colliding
			}
			//particle is colliding
			if (firstCol && !points->getCollision()) //particle of interest has not collided yet
			{
				points->setCollision(true);
				discovered.push_back(points); //store the collided particle
				firstCol = false; //only return the first particle found
			}
		}
	}
	else
	{
		//has children
		for (int i = 0; i < 8; i++)
		{
			//find max and min of children
			Point3D childMax, childMin;
			childMax.x = children[i]->center.x + children[i]->halfSize.x;
			childMax.y = children[i]->center.y + children[i]->halfSize.y;
			childMax.z = children[i]->center.z + children[i]->halfSize.z;

			childMin.x = children[i]->center.x - children[i]->halfSize.x;
			childMin.y = children[i]->center.y - children[i]->halfSize.y;
			childMin.z = children[i]->center.z - children[i]->halfSize.z;

			if (childMax.x < min.x || childMax.y < min.y || childMax.z < min.z)
			{
				continue;
			}
			if (childMin.x > max.x || childMin.y > max.y || childMin.z > max.z)
			{
				continue;
			}
			//The child might have a point within the sphere
			children[i]->getParticleRange(pCenter, r, discovered, firstCol);
		}
	}
}

/*
	Returns the squared distance between p1 and p2.
*/
double Detector::getDist(Point3D p1, Point3D p2)
{
	double d;
	double d_x = p1.x - p2.x;
	double d_y = p1.y - p2.y;
	double d_z = p1.z - p2.z;
	d = d_x * d_x + d_y * d_y + d_z * d_z;
	return d;
}

/*
	Checks if the two particles are heading toward each other
	to generate a collision. Returns true if the particles
	will collide and false otherwise.
*/

bool Detector::isParallel(particle *p1, particle *p2)
{
	//cout << "(" << p1->getParticleName() << ", " << p2->getParticleName() << ")" << endl;

	double *v1 = p1->getParticleVelocity();
	double *pos1 = p1->getParticlePosition();
	double *v2 = p2->getParticleVelocity();
	double *pos2 = p2->getParticlePosition();

	double dv[3]; //store change in velocity
	double dp[3]; //store change in displacement

	double magV = 0;
	double sum = 0;
	for (int i = 0; i < 3; i++)
	{
		dp[i] = pos1[i] - pos2[i];
		dv[i] = v1[i] - v2[i];

		magV += dv[i] * dv[i];
		sum += dp[i] * dv[i];
	}
	if (sum < 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}