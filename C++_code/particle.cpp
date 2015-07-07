#include "particle.h"


particle::particle(){

};

void particle::set_cartesian(double* Values){
    for (int i = 0; i < 3; ++i)
        this->cartesian[i]=Values[i];
}
void particle::set_polar(double* Values){
    for (int i = 0; i < 3; ++i)
        this->polar[i]=Values[i];

}
void particle::reassign_cartesian(double radius){
    double theta=this->polar[0];
    double psi=this->polar[1];
    double phi=this->polar[2];

    double cos_psi=cos(psi);
    double sin_psi=sin(psi);
    double cos_theta=cos(theta);
    double sin_theta=sin(theta);
    double cos_phi=cos(phi);
    double sin_phi=sin(phi);

    this->cartesian[0]=radius * sin_psi * cos_theta;
    this->cartesian[1]=radius * sin_psi * sin_theta * cos_phi;
    this->cartesian[2]=radius * cos_psi;
    this->cartesian[3]=radius * sin_psi * sin_theta * sin_phi;

}