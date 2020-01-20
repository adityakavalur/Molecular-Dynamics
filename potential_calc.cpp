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
using namespace std;


void potential_calc(neighborlist_type *neighlist, config_type *fromfile, potential *run_parameters)
{

//double PE, KE;//total system PE & KE
double *dist, disp, disp2, rb2, **potential_correction, sig_r6;
double force_temp;
dist = new double [Dim];
rb2 = pow(run_parameters->rb,2);



//calculating the potential at cut-off for each combination of atoms
potential_correction= new double *[fromfile->atom_types];
for (int i=0; i<fromfile->atom_types; i++)
    {
    potential_correction[i] = new double [fromfile->atom_types];
    }



for (int i=0; i<fromfile->atom_types; i++)
        {
        for (int j=0; j<fromfile->atom_types; j++)
            {
            potential_correction[i][j] = (4*run_parameters->epsilon[i][j]*(pow(run_parameters->sigma[i][j]/2.5,12)-pow(run_parameters->sigma[i][j]/2.5,6)));
            }
        }


//setting 0 for the new time step force to be calculated
run_parameters->PotEnergy=0.0;



for(long i=0; i<fromfile->N; i++)
    {
    run_parameters->pe[i]=0;
    run_parameters->ke[i]=0;

    for (int j=0; j<Dim; j++)
        {
        run_parameters->force[i][j]=0;
        }
    }




//calculating the force directly, without potential, since adding a constant wont affect the differential
for (long i=0; i<fromfile->N; i++)
    {
    long atom_p=i;




    for (long j=0; j<neighlist[i].countsn; j++)
        {

        long atom_s=neighlist[i].numb[j];

        for (int k=0; k<Dim; k++)
            {
            dist[k]= ((fromfile->r[atom_s][k] + neighlist[atom_p].xyz[j][k]*fromfile->L[k]) - fromfile->r[atom_p][k]);
            }
        disp2 = pow(dist[0],2)+pow(dist[1],2)+pow(dist[2],2);
        if (disp2<rb2)
            {
            disp=pow(disp2,0.5);



            int c1=fromfile->c[atom_p];
            int c2=fromfile->c[atom_s];





            sig_r6=pow(run_parameters->sigma[c1-1][c2-1]/disp,6);//(sigma/r)^6




            //calculating the force-scalar
            run_parameters->pe[atom_p] += 0.5*(4*run_parameters->epsilon[c1-1][c2-1]*((sig_r6*sig_r6)-sig_r6)) - potential_correction[c1-1][c2-1];

            force_temp=-0.5*(((48.0*run_parameters->epsilon[c1-1][c2-1])/disp2)*((sig_r6*sig_r6)-((sig_r6)*0.5)));



            for (int k=0; k<Dim;k++)
                {
                run_parameters->force[atom_p][k]+=force_temp*(dist[k]);
                run_parameters->force[atom_s][k]-=force_temp*(dist[k]);
                }

            }
        }
    }


//total potential energy
for (long i=0; i<fromfile->N; i++)
    {
    run_parameters->PotEnergy += run_parameters->pe[i];

    }

//delete potential correction at the end
for (int i=0; i<fromfile->atom_types; i++)
    {
    delete [] potential_correction[i];
    }
delete [] potential_correction;
}
