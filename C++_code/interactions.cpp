#include "interactions.h"


double Lennard_Jones(double rij2, double epsilon, double sigma, double cutoff){
    double i_rij6=1./(rij2*rij2*rij2);
    double sigma6=sigma*sigma*sigma*sigma*sigma*sigma;

    double U=4*epsilon*(sigma6*i_rij6-1.)*sigma6*i_rij6;

    return U-cutoff;
}

double Repulsive_Lennard_Jones(double rij2){

    //a la Straley
    return 4/(rij2*rij2*rij2*rij2*rij2*rij2);
}


double Interact(particle &a, particle &b,double dr2, int Interaction_Type){
    // dr 
    if (Interaction_Type==1)
    {   
        if(dr2<2.5*2.5)
        return Lennard_Jones(dr2, EPSILON, SIGMA, CUT);
        else return 0;
  
    }
    else if (Interaction_Type==2)
    {   
    
        return Repulsive_Lennard_Jones(dr2);
  
    }
    else if(Interaction_Type==3)
    {
        if(a.type*b.type==1)
        return Lennard_Jones(dr2,EPSILON_AA, SIGMA_AA,CUTOFF_AA);
        else if(a.type*b.type==2)
        return Lennard_Jones(dr2,EPSILON_AB, SIGMA_AB,CUTOFF_AB);
        else 
        return Lennard_Jones(dr2,EPSILON_BB, SIGMA_BB,CUTOFF_BB);       
    }
    else{std::cout<<"Wrong potential type"<<std::endl;exit(-1) ;}
}

