#include "model.h"


model::model(double Radius_value, int N){
        this->radius=Radius_value;
        this->Npart=N;
        this->particles=new particle [N];
        this->Acceptance=0;
};

void model::add_random_particles(int N, int Type){

    for (int i = 0; i < N; ++i)
    {
        double values[3];
        // theta
        values[0]=rnd_real(0,M_PI);
        // psi
        values[1]=rnd_real(0,M_PI);
        // phi
        values[2]=rnd_real(0,2*M_PI);

        this->particles[i].set_polar(values);
        this->particles[i].type=Type;
        this->particles[i].reassign_cartesian(this->radius);
    }
    // this->initialise_Verlet();
};

void model::load_configuration(std::ifstream &fin, int N, int Type){

    for (int i = 0; i < N; ++i)
    {
        double values[3];
        // theta
        fin>>values[0]
        // psi
        >>values[1]
        // phi
        >>values[2];

        this->particles[i].set_polar(values);
        this->particles[i].type=Type;
        this->particles[i].reassign_cartesian(this->radius);
    }

};

// void model::read_configuration(std::ifstream &Fin,int Iteration){

// }
void model::set_interaction(std::string InteractionType){
    if(InteractionType=="LJ") {
        this->Interaction=1;
         std::cout<<"* Interaction set to LJ"<<std::endl;
    }
    else if(InteractionType=="Repulsive LJ") 
        {
            this->Interaction=2;
        std::cout<<"* Interaction set to Repulsive LJ"<<std::endl;
    }
    else if(InteractionType=="Wahnstrom"){
        this->Interaction=3;
        std::cout<<"* Interaction set to Wahnstrom"<<std::endl;
    }
    else{
        std::cerr<<"* \n===> UNDEFINED INTERACTION\nInteraction \""<<InteractionType<<"\" not defined. The only available options are\n - LJ\n - Repulsive LJ\n\n Check for a typo.\n====";
        exit(-1);
    }

};

void model::write_polar_configuration(std::ofstream &Fout, int Iteration){
    // write Header
    Fout<<Npart<<"\nITERATION\t"<<Iteration<<std::endl;

    for (int i = 0; i < this->Npart; ++i)
    {   
        if(this->particles[i].type==1)
        Fout<<'A'<<'\t';
        else if (this->particles[i].type==2)
        Fout<<'B'<<'\t';
        for (int k = 0; k < 3; ++k)
            Fout<<this->particles[i].polar[k]<<'\t';
        Fout<<std::endl;
    }
}

// void model::initialise_Verlet(){
//     for (int i = 0; i < this->Npart; ++i)
//     {   for (int k = 0; k < 3; ++k)
//         {
//             this->polar_coords_for_Verlet.push_back(this->particles[i].polar[k]);
//         }
        
//     }
// }
void model::build_Verlet_lists(double Verlet_radius){
    // store position of particles
    for (int i = 0; i < this->Npart; ++i)
        {
            for (int k = 0; k < 3; ++k)
            {
                this->polar_coords_for_Verlet[i*3+k]=this->particles[i].polar[k];

            }
            this->particles[i].num_of_neighs=0;
        }

    // build local lists of neighbours
    double dist;
    for (int i = 0; i < this->Npart-1; ++i)
    {
        for (int j = i+1; j < this->Npart; ++j)
        {
            dist=distance(particles[i], particles[j],this->radius);
            if(dist<Verlet_radius)
            {
                particles[i].neighbours[particles[i].num_of_neighs]=j;
                particles[j].neighbours[particles[j].num_of_neighs]=i;
                particles[i].num_of_neighs++;
                particles[j].num_of_neighs++;
            }
        }
    }
        
    
}

double model::local_energy(particle &P){
    double dE=0,U;
    // std::cout<<"The selected particle has "<<P.num_of_neighs<<" neighbours "<<std::endl;
    for (int i = 0; i < P.num_of_neighs; ++i)
    {
        int neigh_index=P.neighbours[i];
        double dr=distance(P,this->particles[neigh_index],this->radius);
        U=Interact(P,this->particles[neigh_index],dr, this->Interaction);
        
        
        dE+=U;
    }
    return dE;
}
double model::get_total_energy(){
    
    double dr,U,dr2;
    double E=0;
    for (int i = 0; i < this->Npart-1; ++i)
    {   
      

        for (int j = i+1; j < this->Npart; ++j)
        {
    
            dr2=distance2(particles[i], particles[j],this->radius);            
            E+=Interact(particles[i], particles[j],dr2, this->Interaction);   
        }
  
    }
    // std::cout<<"\n========\n\n\n";
    this->Energy=E;
    return E;

}

double model::get_energy(int i){
    double dr,U,dr2;
    double E=0;

    for (int j = 0; j < this->Npart; ++j)
        {
            if(i!=j){
                dr2=distance2(particles[i], particles[j],this->radius);            
                E+=Interact(particles[i], particles[j],dr2, this->Interaction);
            }
        }
    return E;

}

void model::perform_a_Metropolis_move(double Max_angular_perturbation, double Temperature){
    double Eold, Enew,dE;
    // pick a particle
    int selected=rnd_int(0,this->Npart) ;
    Eold=get_energy(selected);
    // store old polar coordinates
    double old_polar[3];
    for (int k = 0; k < 3; ++k)
        old_polar[k]=this->particles[selected].polar[k];
    // perturb the coordinates
    for (int k = 0; k < 3; ++k)
        {
            particles[selected].polar[k]+=rnd_real(-Max_angular_perturbation,Max_angular_perturbation);
        }
    particles[selected].reassign_cartesian(this->radius);
    // evaluate the new energy
    Enew=get_energy(selected);
    dE=Enew-Eold;

    if(rnd_real(0,1)>exp(-dE/Temperature)){
        
        for (int k = 0; k < 3; ++k)
            particles[selected].polar[k]=old_polar[k];
        // std::cout<<"rejected"<<std::endl;
        particles[selected].reassign_cartesian(this->radius);
    }

    else {
        // if accepted
        // std::cout<<dE<<std::endl;
        this->Energy+=dE;
        this->Acceptance++;
    }
    
}


// void model::perform_a_Metropolis_move(double Max_angular_perturbation,double Verlet_radius, double Temperature){
//     double Eold, Enew,dE;
//     // pick a particle
//     int selected=rnd_int(0,this->Npart) ;
//     // check whether we need to rebuild the neighbour lists
//     // unfortunately I have to construct a particle for this scope , 
//     // in order to use the distance function
//     particle Verlet_particle(0);
//     double Verlet_polar[3];
//     for (int k = 0; k < 3; ++k) Verlet_polar[k]=this->polar_coords_for_Verlet[selected*3+k];
//     Verlet_particle.set_polar(Verlet_polar);
//     Verlet_particle.reassign_cartesian(this->radius);

//     if(distance(particles[selected],Verlet_particle, this->radius)<Verlet_radius){
//         build_Verlet_lists(Verlet_radius);
//     }

//     // evaluate old local energy
//     Eold=local_energy(particles[selected]);
//     // store old polar coordinates
//     double old_polar[3];
//     for (int k = 0; k < 3; ++k)
//         old_polar[k]=particles[selected].polar[k];
//     // perturb the coordinates
//     for (int k = 0; k < 3; ++k)
//         {
//             particles[selected].polar[k]+=rnd_real(-Max_angular_perturbation,Max_angular_perturbation);
//             particles[selected].reassign_cartesian(this->radius);

//         }
//     // evaluate the new energy
//     Enew=local_energy(particles[selected]);
//     dE=Enew-Eold;
//     // std::cout<<Eold<<" "<<Enew<<std::endl;
//     // std::cout<<"new E "<<Enew<<std::endl;
//     // reject according to Metropolis rule
//     if(rnd_real(0,1)>exp(-dE/Temperature)){
//         std::cout<"rejected\n";
//         for (int k = 0; k < 3; ++k)
//             particles[selected].polar[k]=old_polar[k];
//         dE=0;
//     }

//     this->Energy+=dE;
// }