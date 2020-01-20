#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include "loadconfig.h"
#include "constants.h"

using namespace std;

void loadconfig (config_type *fromfile, potential *run_parameters, bin_type *bin, const char *initial_file, int check)
{
if (check == 1)
    {
    //file run_parameters.txt is reopened a few times
    double dummy;
    char dummytext[1000];

    fromfile->L = new double [Dim];
    bin->tot_bins = new int [Dim];//total no of bins in each dimension

    long a;
    int n=0;
    //std::string initialfile;
    char configfile[1000];


    //initialfile="run_parameters.txt";
    //const char *initial_file = initialfile.c_str();


    //initial file with parameter inputs
    ifstream File;
    File.open(initial_file);
        if (File==NULL)
        {
        cout << "Check input file name run parameters file"; abort();
        }

    File >> dummy;//already read in neighborlist.cpp
    File >> dummy;//already read in neighborlist.cpp
    File >> configfile;
    File >> run_parameters->rb;
    File >> run_parameters->rc;
    /*sigma epsilon values still to be read from this file, it is read below after
    configfile and declarations since total number of atom types are known from configfile
    */
    File.close();


    //configuration file from disk
    File.open(configfile);

    if (File==NULL)
        {
        cout << "Check configuration file name"; abort();
        }
    File >> fromfile->N >> fromfile->atom_types >> fromfile->L[0] >> fromfile->L[1] >> fromfile->L[2];





    //allocating memory for reading the configuration file
    fromfile->r = new double *[fromfile->N];
    fromfile->xr = new double *[fromfile->N];
    fromfile->v = new double *[fromfile->N];
    fromfile->c = new int [fromfile->N];
    bin->next = new long [fromfile->N];
    run_parameters->force=new double *[fromfile->N];
    run_parameters->pe = new double [fromfile->N];
    run_parameters->ke = new double [fromfile->N];
    for (long i=0; i<fromfile->N; i++)
        {
        fromfile->r[i] = new double [Dim];
        fromfile->xr[i] = new double [Dim];
        fromfile->v[i] = new double [Dim];
        run_parameters->force[i] = new double [Dim];
        }


    run_parameters->sigma = new double *[fromfile->atom_types];
    run_parameters->epsilon = new double *[fromfile->atom_types];
    for (int i=0; i<fromfile->atom_types; i++)
        {
        run_parameters->sigma[i] = new double [fromfile->atom_types];
        run_parameters->epsilon[i] = new double [fromfile->atom_types];
        }

    //reading the configuration file
    while (!File.eof())
        {
        File >> a >> fromfile->c[n] >> fromfile->r[n][0] >> fromfile->r[n][1] >> fromfile->r[n][2] >> fromfile->v[n][0] >> fromfile->v[n][1] >> fromfile->v[n][2];
        if (File) n++;
        }
    File.close();

    while (n != fromfile->N)
        {
        cout << "Data file corrupted"; abort();
        }

    //initial file with run inputs
    File.open(initial_file);
        if (File==NULL)
        {
        cout << "Check input file name"; abort();
        }
    File >> dummy;
    File >> dummy;
    File >> dummytext;
    File >> dummy;
    File >> dummy;
    int values=0;
    for (int i=0; i<fromfile->atom_types; i++)
        {
        for (int j=0; j<fromfile->atom_types; j++)
            {
            File >> run_parameters->sigma[i][j] >> run_parameters->epsilon[i][j];
            if (File) values++;
            }
        }
    if (values!=fromfile->atom_types*fromfile->atom_types) {cout <<"check epsilon & sigma values \n"; abort();}
    File.close();



    for (int i=0; i<Dim; i++)
        {
        bin->tot_bins[i]=0;//setting the default as 0, so that when it goes to binning it knows that bins are not populated
        }


    }





if (check == 0)//deleting dynamic memories in fromfile
    {

    for (int i=0; i<fromfile->N; i++)
        {
        delete [] fromfile->r[i];
        delete [] fromfile->v[i];
        delete [] fromfile->xr[i];
        delete [] run_parameters->force[i];
        }
    for (int i=0; i<bin->tot_bins[0]; i++)
        {
        for (int j=0; j<bin->tot_bins[1]; j++)
            {
            delete [] bin->bin_no[i][j];
            }
        delete [] bin->bin_no[i];
        }
    for (int i=0; i<fromfile->atom_types; i++)
        {
        delete [] run_parameters->sigma[i];
        delete [] run_parameters->epsilon[i];
        }

    delete [] run_parameters->pe;
    delete [] run_parameters->ke;
    delete [] run_parameters->force;
    delete [] fromfile->r;
    delete [] fromfile->xr;
    delete [] fromfile->v;
    delete [] fromfile->L;
    delete [] bin->bin_no;
    delete [] fromfile->c;
    delete [] bin->tot_bins;
    delete [] bin->next;

    }

}
