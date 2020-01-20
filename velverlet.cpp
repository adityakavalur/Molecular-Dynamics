#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include "initializenl.h"
#include "loadconfig.h"
#include "constants.h"
#include "neighcalculation.h"
#include "potential_calc.h"
#include "velverlet.h"
#include "binning.h"
using namespace std;


void velverlet(neighborlist_type *neighlist, config_type *fromfile, potential *run_parameters, bin_type *bin, double dt)
{
double max_disp2[2], temp_disp[3], temp_displacement2, rc2;

rc2= pow(run_parameters->rc,2);

//2 maximum displacements
for (int i=0; i<2; i++)
    {
    max_disp2[i]=0;//i=0 stores the higher value than i=1
    }

run_parameters->KinEnergy=0.0;

//As given in the book Modeling Materials Algorithm 9.2 page no.500

//updating the position vector
for (long i=0; i<fromfile->N; i++)
    {
    for (int j=0; j<Dim; j++)
        {
        fromfile->r[i][j] += dt*fromfile->v[i][j] + 0.5*pow(dt,2)*run_parameters->force[i][j];
        }
    }


//upadting the velocity
for (long i=0; i<fromfile->N; i++)
    {
    for (int j=0; j<Dim; j++)
        {
        fromfile->v[i][j] += (0.5*dt*run_parameters->force[i][j]);
        }
    }




//calculating two maximum displacements against original binning configuration

for (long i=0; i<fromfile->N; i++)
    {
    for (int j=0; j<Dim; j++)
        {
        temp_disp[j] = fromfile->xr[i][j] - fromfile->r[i][j];
        }
    temp_displacement2 = pow(temp_disp[0],2)+pow(temp_disp[1],2)+pow(temp_disp[2],2);

    if (temp_displacement2 > max_disp2[0]) {max_disp2[1]=max_disp2[0]; max_disp2[0]=temp_displacement2;}
    else if (temp_displacement2>max_disp2[1]) {max_disp2[1]=temp_displacement2;}
    }


//check displacement requirement and the binning algorith with if condition
if (rc2<(max_disp2[0]+max_disp2[1]+2*max_disp2[0]*max_disp2[1])) {binning(neighlist, fromfile, run_parameters, bin);}//binning of the configuration file






potential_calc(neighlist,fromfile, run_parameters);//calculating the force fields


//upadting the velocity
for (long i=0; i<fromfile->N; i++)
    {
    for (int j=0; j<Dim; j++)
        {
        fromfile->v[i][j] += (0.5*dt*run_parameters->force[i][j]);
        }
    run_parameters->ke[i] = 0.5*(pow(fromfile->v[i][0],2)+pow(fromfile->v[i][1],2)+pow(fromfile->v[i][2],2));
    run_parameters->KinEnergy+=run_parameters->ke[i] ;
    }



for (long i=0; i<fromfile->N; i++)
    {
    for (int j=0; j<Dim; j++)
        {
        while (fromfile->r[i][j]<0) {fromfile->r[i][j] += fromfile->L[j];}
        while (fromfile->r[i][j]>fromfile->L[j]) {fromfile->r[i][j] -= fromfile->L[j];}
        }
    }





}
