#ifndef __INTERACTIONS_H
#define __INTERACTIONS_H

#define EPSILON 1.
#define SIGMA 1.
#define CUT  -0.016316891136000006 //energy at cutoff 2.5


#define EPSILON_AA 1.
#define SIGMA_AA 1.
#define CUTOFF_AA -0.016316891136000006 

#define EPSILON_BB 1.
#define SIGMA_BB .833333333
#define CUTOFF_BB -0.005479441744238765 //energy at cutoff 2.5

#define EPSILON_AB 1.
#define SIGMA_AB .916666667
#define CUTOFF_AB -0.009696877283178888 //energy at cutoff 2.5

#include "particle.h"
#include <string>
#include <iostream>

double Lennard_Jones(double rij2, double epsilon, double sigma, double cutoff);

double Interact(particle &a, particle &b, double dr, int Interaction_Type);

#endif