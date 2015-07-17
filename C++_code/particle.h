#ifndef __PARTICLE_H
#define __PARTICLE_H

#include <string>
#include <cmath>
#define MAXNEIGHBOURS 30

class particle
{
public:
    particle();

    
    int type;
    // cartesian coordinates, x y z w
    double cartesian[4];
    // polar coordinates, 
    // theta [0,pi]
    // psi  [0,pi]
    // phi  [0,2pi]
    double polar[3];

    int neighbours[MAXNEIGHBOURS];
    int num_of_neighs;

    void set_cartesian(double* Values);
    void set_polar(double* Values);

    void reassign_cartesian(double radius);

};


#endif

