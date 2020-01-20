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
using namespace std;



bool check_bins_scanned(int **bins_scanned, int bin_counter, int s_actual, int t_actual, int u_actual)
{
for (int i=0; i<bin_counter; i++)
    {
    if (bins_scanned[i][0]==s_actual && bins_scanned[i][1]==t_actual && bins_scanned[i][2] ==u_actual)
        {return (false);}
    }

    bins_scanned[bin_counter][0]=s_actual;
    bins_scanned[bin_counter][1]=t_actual;
    bins_scanned[bin_counter][2]=u_actual;
    return (true);
}





void neighcalculation(neighborlist_type *neighlist, config_type *fromfile, potential *run_parameters, bin_type *bin)
{

int nx, ny, nz; // images to be considered - to + in each direction
int **tempxyz, tempcountsn;//temporary storage for values
int inc_neighmax=10, init_neighmax=100, neighmax=init_neighmax;
long *tempnumb;
double rcb2, disp, dist[Dim];
int **bins_scanned, total_bins_scan, bin_counter;




//number of bins/mirrors to be considered
nx = ceil ((run_parameters->rc+run_parameters->rb)/(fromfile->L[0]/bin->tot_bins[0]));
ny = ceil ((run_parameters->rc+run_parameters->rb)/(fromfile->L[1]/bin->tot_bins[1]));
nz = ceil ((run_parameters->rc+run_parameters->rb)/(fromfile->L[2]/bin->tot_bins[2]));
rcb2 = (pow(run_parameters->rc+run_parameters->rb,2));
total_bins_scan = (2*nx+1)*(2*ny+1)*(2*nz+1);



//assigning memory to temporary variables
tempxyz = new int *[neighmax];
tempnumb = new long [neighmax];
bins_scanned = new int *[total_bins_scan];
for (int i=0; i<neighmax; i++)
    {
    tempnumb[i]=-1;
    tempxyz[i] = new int [Dim];
    }
for (int i=0; i<total_bins_scan; i++)
    {
    bins_scanned[i] = new int [Dim];
    }


//calculating the neighborlist going bin-wise
for (long i=0; i<bin->tot_bins[0]; i++)
    {

    for (long j=0; j<bin->tot_bins[1]; j++)
        {
        for (long k=0; k<bin->tot_bins[2]; k++)
            {

            long atom_p=bin->bin_no[i][j][k];//assigning the first atom of the bin to the variable

            while (atom_p !=-1)
                {

                bin_counter=0;//bin counter of real bins
                tempcountsn=0;//counter for number of atoms already in the neighborlist of the current atom

                for (int s=i-nx; s<i+nx+1; s++)
                    {
                    int s_actual=(s%bin->tot_bins[0]); //storing real bin number in s_actual
                    if (s_actual<0) {s_actual += bin->tot_bins[0];}


                    for (int t=j-ny; t<j+ny+1; t++)
                        {
                        int t_actual=(t%bin->tot_bins[1]);
                        if(t_actual<0) {t_actual += bin->tot_bins[1];}


                        for (int u=k-nz; u<k+nz+1; u++)
                            {
                            int u_actual=(u%bin->tot_bins[2]);
                            if (u_actual<0) {u_actual += bin->tot_bins[2];}


                            long atom_s=bin->bin_no[s_actual][t_actual][u_actual];//assigning the first atom of the bin to the variable

                            if ( check_bins_scanned(bins_scanned, bin_counter, s_actual, t_actual, u_actual) == true )
                                {//will go into loop if true/new real bin

                                while (atom_s!= -1)
                                    {
                                    for (int e=-nx; e<nx+1; e++)
                                        {
                                        for (int f=-ny; f<ny+1; f++)
                                            {
                                            for (int g=-nz; g<nz+1; g++)
                                                {

                                                if ((atom_p != atom_s) || ((e!=0) || (f!=0) || (g!=0)))//to avoid considering the atom itself and yet considering the mirrors if necessary
                                                    {
                                                dist[0] = fromfile->r[atom_p][0]-(e*fromfile->L[0]+fromfile->r[atom_s][0]);
                                                dist[1] = fromfile->r[atom_p][1]-(f*fromfile->L[1]+fromfile->r[atom_s][1]);
                                                dist[2] = fromfile->r[atom_p][2]-(g*fromfile->L[2]+fromfile->r[atom_s][2]);
                                                disp = (pow(dist[0],2)+pow(dist[1],2)+pow(dist[2],2));

                                                if (disp < rcb2)
                                                    {

                                                    if (tempcountsn == (neighmax) -1)//if array limit i.e. neighmax of the temp variables (tempxyz,tempnumb) is reached
                                                        {
                                                        int **tempxyz2;
                                                        long *tempnumb2;


                                                        //equating temp2 to temp
                                                        tempxyz2=tempxyz;
                                                        tempnumb2 = tempnumb;

                                                        //increment neighmax value by pre-determined value
                                                        int old_neighmax=neighmax;
                                                        neighmax +=inc_neighmax;


                                                        //initializing temp
                                                        tempnumb = new long [neighmax];
                                                        tempxyz = new int *[neighmax];
                                                        for (int l=0; l<neighmax;l++)
                                                            {
                                                            tempxyz[l] = new int [Dim];
                                                            tempnumb[l]=-1;
                                                            }


                                                        //transferring values from temp2 to temp
                                                        for (int l=0; l<old_neighmax;l++)
                                                            {
                                                            tempxyz[l][0]=tempxyz2[l][0];
                                                            tempxyz[l][1]=tempxyz2[l][1];
                                                            tempxyz[l][2]=tempxyz2[l][2];
                                                            tempnumb[l]=tempnumb2[l];
                                                            }



                                                        //deleting temp2
                                                        delete []tempnumb2;
                                                        for (int l=0; l<old_neighmax; l++)
                                                            {
                                                            delete [] tempxyz2[l];
                                                            }
                                                        delete [] tempxyz2;
                                                        }

                                                    //assigning the secondary atom to the temporary neighbor array
                                                    tempnumb[tempcountsn] = atom_s;
                                                    tempxyz[tempcountsn][0] = e;
                                                    tempxyz[tempcountsn][1] = f;
                                                    tempxyz[tempcountsn][2] = g;


                                                    tempcountsn +=1;
                                                    }

                                                }
                                                }
                                            }
                                        }

                                    atom_s=bin->next[atom_s];
                                    }

                                bin_counter ++;
                                }
                            }
                        }
                    }
                initializenl(neighlist+atom_p, tempcountsn);
                (neighlist+atom_p)->countsn=tempcountsn;//transferring current atom i's neighbourlist from temp to neighlist

                for (int l=0; l<tempcountsn; l++)
                    {
                    (neighlist+atom_p)->numb[l] = tempnumb[l];
                    tempnumb[l]=-1;

                    for (int m=0; m<Dim; m++)
                        {
                        (neighlist+atom_p)->xyz[l][m]=tempxyz[l][m];
                        tempxyz[l][m]=0;
                        }
                    }
                atom_p=bin->next[atom_p];
                }

            }
        }
    }




//deleting temp
delete [] tempnumb;
for (int k=0; k<(neighmax); k++)
    {
    delete [] tempxyz[k];
    }
delete [] tempxyz;

}





