// g++ -Ofast   -std=c++11   *.cpp
#include <iostream>
#include "particle.h"
#include "model.h"
#include "communicate.h"
#include "composestream.h"
#include <ctime>

using namespace std;


int main(int argc, char const *argv[])
{   
    clock_t tStart = clock();
   
    ifstream fin(argv[1],ifstream::in);

    string word,mode,Interaction,filename="none",Folder="./", inputname;
    inputname=argv[1];
    double radius=-1, Temperature=-1,step=-1;
    int N=-1;
    long int Number_of_sweeps=-1;
    int Snapshots=500;
    while (fin.good()){
        fin>>word;
        if(word[0]=='#') getline(fin,word);
        else if(word=="Straley"){
            radius=0.5*(1.+sqrt(5));
            N=120;
            Interaction="Repulsive_LJ";
            step= M_PI/180.*2;//;//radius;//M_PI/180.*0.5;
        }
        else if (word=="Temperature") fin>>Temperature;
        else if (word=="N") fin>>N;
        else if (word=="Radius") fin>>radius;
        else if (word=="Sweeps") fin>>Number_of_sweeps;
        else if (word=="Interaction") fin>>Interaction;
        else if (word=="Step") fin>>step; //step size in length units
        else if (word=="Load") fin>>filename;
        else if (word=="Snapshots") fin>>Snapshots;
        else if (word=="Folder") fin>>Folder;

    }
    fin.close();
    system(("mkdir "+Folder).c_str());
    system(("cp "+inputname+" "+Folder).c_str());

    if(radius<0 || Temperature<0 || Number_of_sweeps<0 ||N<0 ||step<0 ) {cout<<"Radius/N/Temperature/Sweeps/Step are not set. Check your input file."<<endl; exit(0);}
    double Angular_step=step/radius;
    
   
    model Sphere(radius,N);

    ComposeStream   out;
    ofstream file(Folder+"/log.txt", ofstream::out);
    out.linkStream(std::cout);
    out.linkStream(file);
    out<<"* The radius is "<<Sphere.radius<<endl;
    out<<"* The number of particles is "<<Sphere.Npart<<endl;
    out<<"* The temperature is "<<Temperature<<endl;
    out<<"* The density is "<<N/(2*M_PI*M_PI*radius*radius*radius)<<endl;
    out<<"* The linear displacement is "<<step<<endl;
    out<<"* The angular displacement is "<<step/radius<<endl;
    out<<"* The number of sweeps is "<<Number_of_sweeps<<endl;
    out<<"* The interval between snapshots is "<<Snapshots<<endl;

  
    if(filename!="none"){
    out<<"* Loading file "<<filename<<endl;
    ifstream start(filename, ifstream::in);
    Sphere.load_cartesian_configuration( start,  N, 1);
    start.close();
    }
    else{
          // add N random particles of type 1
     Sphere.add_random_particles(N,1);
    }

    Sphere.set_interaction(Interaction);

    // modify type for half of the particles
    // for (int i = N/2; i < N; ++i)
    // {
    //     Sphere.particles[i].type=2;    
    // }
    ofstream Trajectory(Folder+"/output.tj", ofstream::out);
    
    double E=Sphere.get_total_energy();
    Sphere.Energy=E;

    out<<"\n===> The initial energy per particle is "<<E/N<<"\n\n";

    for (int i = 0; i < Number_of_sweeps; ++i)
    {
       for (int j = 0; j < N; ++j)
       {
        Sphere.perform_a_Metropolis_move(step,  Temperature);
          
       }

         
       if(i%Snapshots==0) 
       {
        
        Sphere.get_total_energy();
        print_info(Sphere,i,step);
        save_info(file,Sphere,i,step);
        Sphere.write_cartesian_configuration(Trajectory,i);
        }
       
    }

    out<<"\n===> The final energy per particle is "<<Sphere.get_total_energy()/N<<"\n\n";
    out<<"\n===> The total execution time is "<< ((double)(clock() - tStart)/CLOCKS_PER_SEC)/60<<" minutes."<<endl;
    file.close();
    return 0;
}