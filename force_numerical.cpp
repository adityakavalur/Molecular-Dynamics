#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include "initializenl.h"
#include "loadconfig.h"
#include "constants.h"
#include "potential_calc.h"
#include "force_numerical.h"
using namespace std;


void force_numerical(neighborlist_type *neighlist, config_type *fromfile, potential *run_parameters)
{
//two point method

double h=0.0001;//division for accuracy
double erg1, erg2;
double **tempr, **tempf;
tempr = new double *[fromfile->N];
tempf= new double *[fromfile->N];
for (long i=0; i<fromfile->N; i++)
    {
    tempr[i] = new double [Dim];
    tempf[i] = new double [Dim];
    }

for (long i=0; i<fromfile->N;i++)
    {
    for(int j=0; j<Dim; j++)
        {
        tempr[i][j]=fromfile->r[i][j];
        }
    }
    ofstream myfile ;
    myfile.open ("force_comparison.txt", std::ios_base::app);
for (long i=0; i<fromfile->N;i++)
    {
    for(int j=0; j<Dim; j++)
        {
        fromfile->r[i][j] = tempr[i][j] +h;
        potential_calc(neighlist,fromfile, run_parameters);//calculating the force fields for x+h
        fromfile->r[i][j] = tempr[i][j];
        erg1 = run_parameters->PotEnergy;

        fromfile->r[i][j] = tempr[i][j] -h;
        potential_calc(neighlist,fromfile, run_parameters);//calculating the force fields for x-h
        fromfile->r[i][j] = tempr[i][j];
        erg2 = run_parameters->PotEnergy;
        tempf[i][j] = -(erg1-erg2)/(2*h);

        potential_calc(neighlist,fromfile, run_parameters);//calculating the force fields for x using this for analytical values

        }
    myfile.open ("force_comparison.txt", std::ios_base::app);
    myfile <<i<<" " <<setprecision(15)<< run_parameters->force[i][0] << " " << tempf[i][0] << " " << run_parameters->force[i][1] << " " <<tempf[i][1] << " " << run_parameters->force[i][2] << " " <<tempf[i][2] << endl;

    }

myfile.close();


for (long i=0; i<fromfile->N; i++)
    {
    delete [] tempr[i] ;
    delete [] tempf[i] ;
    }
delete [] tempr;
delete [] tempf;
exit(1);


}
