#ifndef __UTILITIES_H
#define __UTILITIES_H


#include <random>
#include <functional>
#include <iostream>
#include <cmath>
#include "particle.h"




double rnd_real(double a, double b);
double rnd_int(double a, double b);
double mod2(double x, double y, double z, double w);
double euclidean_distance2(particle &a, particle &b);
double distance(particle &a, particle &b, double R);
double distance2(particle &a, particle &b, double R);
template<typename T> T SQ(T x);
double norm(double vec[]);
#endif
