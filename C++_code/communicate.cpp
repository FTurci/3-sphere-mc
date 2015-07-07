#include "communicate.h"

void print_info(const model &Sphere,int MCs, double stepsize){
    // print header
    int N=Sphere.Npart;
    if(MCs==0) std::cout<<std::left<<std::setw(15)<<"Iteration\tStepsize\tAcceptance\tEnergy\n";
    std::cout<<std::left<<
    std::setw(15)<<MCs<<'\t'<<
    std::setw(15)<<stepsize<<'\t'<<
    std::setw(15)<<Sphere.Acceptance/(MCs*N+N)<<'\t'<<
    std::setw(15)<<Sphere.Energy/N<<'\n' 
        ;
}

void save_info(std::ofstream &file,const model &Sphere,int MCs, double stepsize){
    // print header
    int N=Sphere.Npart;
    if(MCs==0) file<<"#Iteration\tStepsize\tAcceptance\tEnergy\n";
    file<<
        MCs<<'\t'<<
        stepsize<<'\t'<<
        Sphere.Acceptance/(MCs*N+N)<<'\t'<<
        Sphere.Energy/N<<'\n' 
        ;
    file.flush();
}