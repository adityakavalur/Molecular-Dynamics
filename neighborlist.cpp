#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include "loadconfig.h"
#include "initializenl.h"
#include "neighcalculation.h"
#include "constants.h"
#include "binning.h"
#include "potential_calc.h"
#include "velverlet.h"
#include "force_numerical.h"
using namespace std;


int main ()
{

bool force_check = false; //turn this boolean true when only checking for numerical differentiation of force against analytical values

int t;//number of time-steps
double dt; //value of each time-step


std::string initialfile;//simulation run parameters initial input file
/*
run_parameters.txt file description
dt - time step size in neighborlist.cpp
t - total real time of simulation in neighborlist.cpp
configuration - configuration file in loadconfig.cpp
rb - force cut-off distance in loadconfig.cpp
rc - additional distance for neighborlist calculation in loadconfig.cpp
sigma epsilon - sigma epsilon values for different atom type combinations each type combo in one line in loadconfig.cpp
*/



neighborlist_type *neighlist; //neighborlist struct defined in header initializenl.h
config_type *fromfile; //read configuration from file into config_type struct defined in loadconfig.h
potential *run_parameters;//read/calculate/take user input on potentials & related inputs, defined in loadconfig.h
bin_type *bin;//bin class defined in loadconfig.h


fromfile = new config_type;
run_parameters=new potential;
bin = new bin_type;

    initialfile="run_parameters.txt";//only place this is defined
    const char *initial_file = initialfile.c_str();


ifstream File;
    File.open(initial_file);
        if (File==NULL)
        {
        cout << "Check input file name simulation time file"; abort();
        }

    File >> dt;
    File >> t;
    File.close();

loadconfig(fromfile, run_parameters, bin, initial_file, 1);//reading configuration from user-defined file


neighlist = new neighborlist_type[fromfile->N];
binning(neighlist, fromfile, run_parameters, bin);//binning of the configuration file



//writing the output file of the neighborlist at the start of the simulation
ofstream myfile ;
myfile.open ("neighbor_list_sim_start.txt");
myfile << "No. of Atoms " << "No. of Atom Types "<< "Length of unit cells in x,y,z directions " << "Cutoff distance" << " "<< "Margin for neighbor list" << "\n" ;
myfile << fromfile->N << "\t" << fromfile->atom_types << "\t" << fromfile->L[0] << "\t" << fromfile->L[1] << "\t" << fromfile->L[2] <<"\t"<< run_parameters->rb << "\t" << run_parameters->rc << "\n" ;
myfile << "Atom," << "\t" << "Neighbor Atom," << "\t" << "Neighbor Atom Type," << "\t" << "Neighbor Atom cell direction XYZ \n";

for (long i=0; i<fromfile->N; i++)
    {
    for (int j=0; j<((neighlist+i)->countsn); j++)
        {
        int k = (neighlist+i)->numb[j];
        myfile << i << " "  << (neighlist+i)->numb[j] << " " <<fromfile->c[k] << " " << (neighlist+i)->xyz[j][0] << " " << (neighlist+i)->xyz[j][1] << " " << (neighlist+i)->xyz[j][2] <<"\n" ;
        }
    }

myfile.close();


//if only forces of numerical and analytical values are to be checked and the program exited
if (force_check==true){force_numerical(neighlist, fromfile, run_parameters);}



potential_calc(neighlist,fromfile, run_parameters);//calculating the force fields

//for printing the energy at zero time step
run_parameters->KinEnergy =0;//setting zero before calculation
for (long i=0; i<fromfile->N; i++)
    {
    run_parameters->ke[i] = 0.5*(pow(fromfile->v[i][0],2)+pow(fromfile->v[i][1],2)+pow(fromfile->v[i][2],2));
    run_parameters->KinEnergy+=run_parameters->ke[i] ;
    }

myfile.open ("energy.txt", std::ios_base::app);
myfile << 0 << " " << run_parameters->KinEnergy << " " << run_parameters->PotEnergy << " " << run_parameters->KinEnergy+run_parameters->PotEnergy  << endl;



for (int i=1; i<=t;i++)
    {
    velverlet(neighlist,fromfile, run_parameters, bin, dt);//calculating the velocities

        {

        ofstream myfile ;
        myfile.open ("energy.txt", std::ios_base::app);
        myfile << i*dt << " " << run_parameters->KinEnergy << " " << run_parameters->PotEnergy << " " << run_parameters->KinEnergy+run_parameters->PotEnergy  << endl;
        }

    }
myfile.close();


//writing the output file of the neighborlist at the end of the simulation

myfile.open ("neighbor_list_sim_end.txt");
myfile << "No. of Atoms " << "No. of Atom Types "<< "Length of unit cells in x,y,z directions " << "Cutoff distance" << " "<< "Margin for neighbor list" << "\n" ;
myfile << fromfile->N << "\t" << fromfile->atom_types << "\t" << fromfile->L[0] << "\t" << fromfile->L[1] << "\t" << fromfile->L[2] <<"\t"<< run_parameters->rb << "\t" << run_parameters->rc << "\n" ;
myfile << "Atom," << "\t" << "Neighbor Atom," << "\t" << "Neighbor Atom Type," << "\t" << "Neighbor Atom cell direction XYZ \n";

for (long i=0; i<fromfile->N; i++)
    {
    for (int j=0; j<((neighlist+i)->countsn); j++)
        {
        int k = (neighlist+i)->numb[j];
        myfile << i << " "  << (neighlist+i)->numb[j] << " " <<fromfile->c[k] << " " << (neighlist+i)->xyz[j][0] << " " << (neighlist+i)->xyz[j][1] << " " << (neighlist+i)->xyz[j][2] <<"\n" ;
        }
    }
myfile.close();


myfile.open ("end_configuration.txt");

for (long i=0; i<fromfile->N; i++)
    {
    myfile << i << " "  << fromfile->r[i][0] << " " << fromfile->r[i][1] << " " << fromfile->r[i][2] <<"\n" ;

    }

myfile.close();


for (int i=0; i<fromfile->N; i++)
    {
    initializenl(neighlist+i, 0);
    }
delete neighlist;

loadconfig(fromfile, run_parameters, bin, initial_file, 0);
delete fromfile;


}
