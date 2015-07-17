#ifndef __model_H
#define __model_H
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <complex>
#include <sstream>
#include <iomanip>
#include "particle.h"
#include "utilities.h"
#include "interactions.h"





class model
{
public:
    model(double Radius,int N);
    
    int Npart;
    double radius;
    double Energy;

    double Acceptance;

    int Interaction;

    particle * particles;

    double * polar_coords_for_Verlet; //3N

    void add_random_particles(int N, int Type);
    void initialise_Verlet();
    void set_interaction(std::string Type);

    //output
    void write_polar_configuration(std::ofstream &Fout, int Iteration);
    void write_cartesian_configuration(std::ofstream &Fout, int Iteration);
    void load_cartesian_configuration(std::ifstream &fin, int N, int Type);
    
    void build_Verlet_lists(double Verlet_radius);
  
    void perform_a_Metropolis_move(double Max_angular_perturbation,double Verlet_radius, double Temperature);

    void perform_a_Metropolis_move(double Max_angular_perturbation, double Temperature);

    double local_energy(particle &P);
    double get_total_energy();
    double get_energy(int i);

    void control_step( double &step_value,double elapsed_time);


};

#endif