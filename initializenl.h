#ifndef _initializenl
#define _initializenl

class neighborlist_type
{
public:
	int **xyz;
	long *numb;
	int countsn;// countsn array counter of each atom till which neighborlist is filled
};

void initializenl(neighborlist_type *neighlist, int neighmax);

#endif
