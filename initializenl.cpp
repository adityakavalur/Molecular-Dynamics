#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include "initializenl.h"
#include "constants.h"
using namespace std;


void initializenl(neighborlist_type *neighlist, int neighmax)
{
// (re)initialising the neighborlist with the required number of neighbors
if (neighmax>0)
    {
    neighlist->xyz = new int *[neighmax];
    neighlist->numb = new long [neighmax];
    neighlist->countsn = 0;
    for (long j=0; j<neighmax; j++)
        {
        neighlist->numb[j] =-1;
        neighlist->xyz[j] = new int [Dim];
        for (int k=0; k<Dim; k++)
            {
            neighlist->xyz[j][k]=0;
            }
        }
    }



else
    {
    //deleting the dynamic memories of neighborlist

    int i;
    for (i=0; i<((neighlist)->countsn);i++)
        {
        delete [] neighlist->xyz[i];
        }

    delete [] neighlist->xyz;
    delete neighlist->numb;
    }
}
