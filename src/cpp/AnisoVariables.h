
#ifndef ANISOVARIABLES_H_
#define ANISOVARIABLES_H_

const int N_max = 100;			// max number of iterations
const double tol_dX = 1.e-4;	// tolerance error for dX
const double tol_F = 1.e-4;		// tolerance error for F

typedef struct
{
	double lambda;
	double aT;
	double aL;
	bool did_not_find_solution;
	int number_of_iterations;

} aniso_variables;

aniso_variables find_anisotropic_variables(double e, double pl, double pt, double lambda_0, double aT_0, double aL_0);


#endif




