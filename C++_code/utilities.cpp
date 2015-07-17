#include "utilities.h"

std::mt19937 engine;

double rnd_real(double a, double b){
    // random number generator, needs c++11

    std::uniform_real_distribution<double> Uniform(a,b);
    return Uniform(engine);
}

double rnd_int(double a, double b){
    // random number generator, needs c++11
    std::uniform_int_distribution<int> Uniform(a,b-1);
    return Uniform(engine);
}


template<typename T> T SQ(T x) { return x * x; }

double mod2(double x, double y, double z, double w){return (x*x)+(y*y)+(z*z)+(w*w);}

double norm(double vec[]){return sqrt(SQ(vec[0])+SQ(vec[1])+SQ(vec[2])+SQ(vec[3]) );}

double euclidean_distance2(particle &a, particle &b){
    double dx,dy,dz,dw;
    dx=a.cartesian[0]-b.cartesian[0];
    dy=a.cartesian[1]-b.cartesian[1];
    dz=a.cartesian[2]-b.cartesian[2];
    dw=a.cartesian[3]-b.cartesian[3];
    return  mod2(dx,dy,dz,dw);
}

double distance(particle &a, particle &b, double R){

    // a la Straley
    double x=sqrt(euclidean_distance2(a,b));
    // return 2*R*acos(x/R*0.5);
    return x;
    // return R * acos(1 - euclidean_distance2(a,b) / (2.*R*R) );
}

double distance2(particle &a, particle &b, double R){

    // a la Straley
    return euclidean_distance2(a,b);
    // double x= R * acos(1 - euclidean_distance2(a,b) / (2.*R*R) );
    // return x*x;
}



