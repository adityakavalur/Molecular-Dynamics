#ifndef _loadconfig
#define _loadconfig

class config_type
{
public:
	long N;
	double **r, **v, *L, **xr;
	int *c, atom_types;
};


class bin_type
{
public:
	int *tot_bins;
	long ***bin_no, *next;
};


class potential
{
public:
	double rb, rc, PotEnergy, KinEnergy;
	double **sigma, **epsilon;
	double **force, *pe, *ke;
};


void loadconfig(config_type *fromfile, potential *run_parameters, bin_type *bin, const char *initial_file, int check);

#endif
