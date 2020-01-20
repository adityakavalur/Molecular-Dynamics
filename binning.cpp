#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include "loadconfig.h"
#include "constants.h"
#include "initializenl.h"
#include "neighcalculation.h"
#include "binning.h"
using namespace std;


void binning (neighborlist_type *neighlist, config_type *fromfile, potential *run_parameters, bin_type *bin)

{



int *temp_bin;
temp_bin = new int [Dim];//saves the bin number temporarily

double *bin_length;
bin_length = new double [Dim];//length of each bin

    //deleting the old bin list
    for (int i=0; i<bin->tot_bins[0]; i++)
        {
        for (int j=0; j<bin->tot_bins[1]; j++)
            {
            delete [] bin->bin_no[i][j];
            }
        delete [] bin->bin_no[i];
        }



    //identifying the number of bins, so that if there are multiple bins in the cell, they are sized such that only the nearest neighbor needs to scanned for neighborlist
    for (int i=0; i<Dim; i++)
        {
        bin->tot_bins[i]= floor(fromfile->L[i]/(run_parameters->rc+run_parameters->rb));
        if (bin->tot_bins[i]==0) {bin->tot_bins[i]=1;}
        bin_length[i] = fromfile->L[i]/bin->tot_bins[i];
        }


    //allocating memory based on the number of bins in each dimension
    bin->bin_no = new long **[bin->tot_bins[0]];
    for (int i=0; i< bin->tot_bins[0]; i++)
        {
        bin->bin_no[i] = new long *[bin->tot_bins[1]];

        for (int j=0; j<bin->tot_bins[1]; j++)
            {
            bin->bin_no[i][j] = new long [bin->tot_bins[2]];

            for (int k=0; k<bin->tot_bins[2]; k++)
                {
                bin->bin_no[i][j][k]=-1;
                }
            }
        }
        //setting the bin next to default -1 value
    for (long i=0; i<fromfile->N; i++)
        {
        bin->next[i] = -1;
        }


    //re-allocating the atoms to the primary cell
    for (long i=0; i<fromfile->N; i++)
        {
        for (int j=0; j<Dim; j++)
            {
            double quotient;


            quotient = floor(fromfile->r[i][j]/fromfile->L[j]);

            fromfile->r[i][j] -= quotient*fromfile->L[j];


            if (fromfile->r[i][j] > fromfile->L[j] || fromfile->r[i][j] < 0) {cout << fromfile->r[i][j] << " atoms not in primary cell" << endl; abort();}
            }
        }

//saving the binning configuration for future use in max displacement
for (long i=0; i<fromfile->N; i++)
    {
    for (int j=0; j<Dim; j++)
        {
        fromfile->xr[i][j]=fromfile->r[i][j];
        }
    }



    //binning the atoms
    for (long i=0; i<fromfile->N; i++)
        {
        for (int j=0; j<Dim; j++) // single atom taken
            {
            temp_bin[j]= floor(fromfile->r[i][j]/bin_length[j]);
            if (temp_bin[j]>(bin->tot_bins[j]-1) || temp_bin[j]<0) {cout << "Check binnning" << endl; abort();}//to compensate for any truncation related errors
            }


        //adding it to the identified bin
        if (bin->bin_no[temp_bin[0]][temp_bin[1]][temp_bin[2]] == -1)
            {

            bin->bin_no[temp_bin[0]][temp_bin[1]][temp_bin[2]]=i;
            }
        else
            {

            long j = bin->bin_no[temp_bin[0]][temp_bin[1]][temp_bin[2]];
            while (bin->next[j] != -1)
                {
                j = bin->next[j];
                }
            bin->next[j]=i;
            }
        }


neighcalculation(neighlist,fromfile, run_parameters, bin);//sending for neighborlist calculation
}
