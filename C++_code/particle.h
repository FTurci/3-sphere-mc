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


inline double sincos(double x, double &sin, double &cos){
    //always wrap input angle to -PI..PI
    if (x < -3.14159265)
        x += 6.28318531;
    else
    if (x >  3.14159265)
        x -= 6.28318531;

    //compute sine
    if (x < 0)
    {
        sin = 1.27323954 * x + .405284735 * x * x;
        
        if (sin < 0)
            sin = .225 * (sin *-sin - sin) + sin;
        else
            sin = .225 * (sin * sin - sin) + sin;
    }
    else
    {
        sin = 1.27323954 * x - 0.405284735 * x * x;
        
        if (sin < 0)
            sin = .225 * (sin *-sin - sin) + sin;
        else
            sin = .225 * (sin * sin - sin) + sin;
    }
        //compute cosine: sin(x + PI/2) = cos(x)
    x += 1.57079632;
    if (x >  3.14159265)
        x -= 6.28318531;

    if (x < 0)
    {
        cos = 1.27323954 * x + 0.405284735 * x * x;
        
        if (cos < 0)
            cos = .225 * (cos *-cos - cos) + cos;
        else
            cos = .225 * (cos * cos - cos) + cos;
    }
    else
    {
        cos = 1.27323954 * x - 0.405284735 * x * x;

        if (cos < 0)
            cos = .225 * (cos *-cos - cos) + cos;
        else
            cos = .225 * (cos * cos - cos) + cos;
    }
}
#endif

