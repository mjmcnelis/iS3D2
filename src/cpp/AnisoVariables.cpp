#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <gsl/gsl_linalg.h>
#include "AnisoVariables.h"
#include "iS3D.h"
#include "Macros.h"
#include "Arsenal.h"

using namespace std;


void compute_F(double Ea, double PTa, double PLa, int Nparticles, double *Mass, double *Sign, double *Degeneracy, double * Baryon, double * X, double * F)
{
	double lambda = X[0];
	double aT = X[1];
	double aL = X[2];

	double aT2 = aT * aT;				// useful expressions
	double aL2 = aL * aL;
	double aT2_minus_aL2 = aT2 - aL2;
	double common_factor = aT2 * aL * lambda * lambda * lambda * lambda / four_pi2_hbarC3;

	double I_200 = 0;					// anisotropic integrals
	double I_220 = 0;
	double I_201 = 0;

	for(int n = 0; n < Nparticles; n++)	// loop over hadrons in PDG
	{
		double mass = Mass[n];			// particle info from PDG arrays
		double sign = Sign[n];
		double degeneracy = Degeneracy[n];
		double baryon = Baryon[n];

		if(mass == 0)
		{
			continue;					// skip photons
		}

		double mbar = mass / lambda;
		double mbar2 = mbar * mbar;
		double chem = 0;				// replace w/ baryon * upsilonB (anisotropic version of baryon * alphaB)

		double I_200_n = 0;				// individual hadron contribution to anisotropic integrals
		double I_220_n = 0;
		double I_201_n = 0;

		for(int i = 0; i < pbar_pts; i++)		// radial momentum integration loop
		{
			double pbar =     pbar_root_a2[i];	// pbar roots and weights for a = 2 (a = n + s)
			double weight = pbar_weight_a2[i];

			double Ebar = sqrt(pbar * pbar  +  mbar2);	// E / lambda

			double w = sqrt(aL2  +  mbar2 / (pbar * pbar));
			double z = aT2_minus_aL2 / (w * w);

			double t_200;						// hypergeometric functions
			double t_220;
			double t_201;

			if(z > delta)						// compute hypergeometric functions
			{
				double sqrtz = sqrt(z);
				double t = atan(sqrtz) / sqrtz;

				t_200 = 1.  +  (1. + z) * t;
				t_220 = (-1.  +  (1. + z) * t) / z;
				t_201 = (1.  +  (z - 1.) * t) / z;
			}
			else if(z < -delta && z > -1.)
			{
				double sqrtmz = sqrt(-z);
				double t = atanh(sqrtmz) / sqrtmz;

				t_200 = 1.  +  (1. + z) * t;
				t_220 = (-1.  +  (1. + z) * t) / z;
				t_201 = (1.  +  (z - 1.) * t) / z;
			}
			else if(fabs(z) <= delta)
			{
				double z2 =  z * z;
				double z3 = z2 * z;
				double z4 = z3 * z;
				double z5 = z4 * z;
				double z6 = z5 * z;

				t_200 = 2. + 0.6666666666666667*z - 0.1333333333333333*z2 + 0.05714285714285716*z3 - 0.031746031746031744*z4 + 0.020202020202020193*z5 - 0.013986013986013984*z6;

				t_220 = 0.6666666666666667 - 0.1333333333333333*z + 0.05714285714285716*z2 - 0.031746031746031744*z3 + 0.020202020202020193*z4 - 0.013986013986013984*z5 + 0.010256410256410262*z6;

				t_201 = 1.3333333333333333 - 0.5333333333333333*z + 0.34285714285714286*z2 - 0.25396825396825395*z3 + 0.20202020202020202*z4 - 0.16783216783216784*z5 + 0.14358974358974358*z6;
			}
			else
			{
			#ifdef FLAGS
				printf("compute_F flag: z = %lf is out of bounds\n", z);
			#endif
			}

			double common_weight = pbar * weight * exp(pbar) / (exp(Ebar + chem) + sign);

			I_200_n += common_weight * t_200 * w;
			I_220_n += common_weight * t_220 / w;
			I_201_n += common_weight * t_201 / w;
		}

		I_200_n *= degeneracy;			// multiply by degeneracy factor
		I_220_n *= degeneracy;
		I_201_n *= degeneracy;

		I_200 += I_200_n;				// add individual hadron contribution
		I_220 += I_220_n;
		I_201 += I_201_n;
	}

	I_200 *= common_factor;				// multiply by common factor (and other things)
	I_220 *= common_factor * aL2;
	I_201 *= common_factor * aT2 / 2.;

	// printf("(E, I_200) = (%lf, %lf)\n", Ea, I_200);
	// printf("(PLa, I_220) = (%lf, %lf)\n", PLa, I_220);
	// printf("(PTa, I_201) = (%lf, %lf)\n\n", PTa, I_201);
	// exit(-1);

	F[0] = I_200 - Ea;					// compute F
	F[1] = I_201 - PTa;
	F[2] = I_220 - PLa;
}


void compute_J(double Ea, double PTa, double PLa, int Nparticles, double *Mass, double *Sign, double *Degeneracy, double * Baryon, double * X, double * F, double ** J)
{
	double lambda = X[0];
	double aT = X[1];
	double aL = X[2];

	double aT2 = aT * aT;				// useful expressions
	double aL2 = aL * aL;
	double aT2_minus_aL2 = aT2 - aL2;
	double lambda2 = lambda  * lambda;
  	double lambda3 = lambda2 * lambda;
	double lambda_aT3 = lambda * aT2 * aT;
  	double lambda_aL3 = lambda * aL2 * aL;
	double common_factor = aT2 * aL * lambda2 * lambda3 / four_pi2_hbarC3;

    double J_2001 = 0;					// anisotropic integrals
    double J_2011 = 0;
    double J_2201 = 0;

    double J_402m1 = 0;
    double J_421m1 = 0;
    double J_440m1 = 0;					// this is wrong (missing quantum statistics derivative)

    for(int n = 0; n < Nparticles; n++)	// loop over hadrons in PDG
	{
		double mass = Mass[n];			// particle info from PDG arrays
		double sign = Sign[n];
		double degeneracy = Degeneracy[n];
		double baryon = Baryon[n];

		if(mass == 0)
		{
			continue;					// skip photons
		}

		double mbar = mass / lambda;
		double mbar2 = mbar * mbar;
		double chem = 0;				// replace w/ baryon * upsilonB (anisotropic version of baryon * alphaB)

		double J_2001_n = 0;			// individual hadron contribution to anisotropic integrals
		double J_2011_n = 0;
    	double J_2201_n = 0;

    	double J_402m1_n = 0;
    	double J_421m1_n = 0;
    	double J_440m1_n = 0;

    	for(int i = 0; i < pbar_pts; i++)		// radial momentum integration loop
		{
			double pbar =     pbar_root_a3[i];	// pbar roots and weights for a = 3 (a = n + s)
			double weight = pbar_weight_a3[i];

			double pbar2 = pbar * pbar;
			double Ebar = sqrt(pbar2  +  mbar2);

			double w = sqrt(aL2  +  mbar2 / pbar2);
			double z = aT2_minus_aL2 / (w * w);
			double z2 =  z * z;

			double t_200;						// hypergeometric functions
			double t_201;
			double t_220;

	  		double t_402;
	  		double t_421;
	  		double t_440;

			if(z > delta)						// compute hypergeometric functions
			{
				double sqrtz = sqrt(z);
				double t = atan(sqrtz) / sqrtz;

				t_200 = 1.  +  (1. + z) * t;
				t_220 = (-1.  +  (1. + z) * t) / z;
				t_201 = (1.  +  (z - 1.) * t) / z;

				t_402 = (3. * (z - 1.)  +  (z * (3.*z - 2.) + 3.) * t) / (4. * z2);
				t_421 = (3.  +  z  +  (1. + z) * (z - 3.) * t) / (4. * z2);
				t_440 = (-(3. + 5.*z) + 3. * (z + 1.) * (z + 1.) * t) / (4. * z2);
			}
			else if(z < -delta && z > -1.)
			{
				double sqrtmz = sqrt(-z);
				double t = atanh(sqrtmz) / sqrtmz;

				t_200 = 1.  +  (1. + z) * t;
				t_220 = (-1.  +  (1. + z) * t) / z;
				t_201 = (1.  +  (z - 1.) * t) / z;

				t_402 = (3. * (z - 1.)  +  (z * (3.*z - 2.) + 3.) * t) / (4. * z2);
				t_421 = (3.  +  z  +  (1. + z) * (z - 3.) * t) / (4. * z2);
				t_440 = (-(3. + 5.*z) + 3. * (z + 1.) * (z + 1.) * t) / (4. * z2);
			}
			else if(fabs(z) <= delta)
			{
				double z3 = z2 * z;
				double z4 = z3 * z;
				double z5 = z4 * z;
				double z6 = z5 * z;

				t_200 = 2. + 0.6666666666666667*z - 0.1333333333333333*z2 + 0.05714285714285716*z3 - 0.031746031746031744*z4 + 0.020202020202020193*z5 - 0.013986013986013984*z6;

				t_220 = 0.6666666666666667 - 0.1333333333333333*z + 0.05714285714285716*z2 - 0.031746031746031744*z3 + 0.020202020202020193*z4 - 0.013986013986013984*z5 + 0.010256410256410262*z6;

				t_201 = 1.3333333333333333 - 0.5333333333333333*z + 0.34285714285714286*z2 - 0.25396825396825395*z3 + 0.20202020202020202*z4 - 0.16783216783216784*z5 + 0.14358974358974358*z6;

				t_402 = 1.0666666666666667 - 0.4571428571428572*z + 0.3047619047619048*z2 - 0.23088023088023088*z3 + 0.1864801864801865*z4 - 0.15664335664335666*z5 + 0.13514328808446457*z6;

				t_421 = 0.2666666666666666 - 0.0761904761904762*z + 0.0380952380952381*z2 - 0.023088023088023088*z3 + 0.015540015540015537*z4 - 0.011188811188811189*z5 + 0.00844645550527904*z6;

					t_440 = 0.4 - 0.057142857142857106*z + 0.019047619047619063*z2 - 0.008658008658008663*z3 + 0.004662004662004657*z4 - 0.002797202797202792*z5 + 0.0018099547511312257*z6;
			}
			else
			{
			#ifdef FLAGS
				printf("compute J flag: z = %lf is out of bounds\n", z);
			#endif
			}

			double quantum_stat = exp(Ebar + chem) + sign;
			double common_weight = weight * exp(pbar + Ebar) / (quantum_stat * quantum_stat);

			J_2001_n += Ebar * common_weight * t_200 * w;
			J_2011_n += Ebar * common_weight * t_201 / w;
			J_2201_n += Ebar * common_weight * t_220 / w;

			J_402m1_n += pbar2 / Ebar * common_weight * t_402 / w;
	    	J_421m1_n += pbar2 / Ebar * common_weight * t_421 / w;
	    	J_440m1_n += pbar2 / Ebar * common_weight * t_440 / w;
		}

		J_2001_n *= degeneracy;			// multiply by degeneracy factor
		J_2011_n *= degeneracy;
    	J_2201_n *= degeneracy;

    	J_402m1_n *= degeneracy;
    	J_421m1_n *= degeneracy;
    	J_440m1_n *= degeneracy;


		J_2001 += J_2001_n;				// add individual hadron contribution
		J_2011 += J_2011_n;
    	J_2201 += J_2201_n;

    	J_402m1 += J_402m1_n;
    	J_421m1 += J_421m1_n;
    	J_440m1 += J_440m1_n;
    }

	J_2001 *= common_factor;			// multiply by common factor (and other things)
	J_2011 *= common_factor * aT2 / 2.;
	J_2201 *= common_factor * aL2;

	J_402m1 *= common_factor * aT2 * aT2 / 8.;
	J_421m1 *= common_factor * aT2 * aL2 / 2.;
	J_440m1 *= common_factor * aL2 * aL2;

	double Eai =  F[0] + Ea;     		// compute Eai, PTai, PLai from F
	double PTai = F[1] + PTa;
	double PLai = F[2] + PLa;

	// compute Jacobian
    J[0][0] = J_2001 / lambda2;			J[0][1] = 2. * (Eai + PTai) / aT;		J[0][2] = (Eai + PLai) / aL;
    J[1][0] = J_2011 / lambda2;			J[1][1] = 4. * J_402m1 / lambda_aT3;	J[1][2] = J_421m1 / lambda_aL3;
    J[2][0] = J_2201 / lambda2;			J[2][1] = 2. * J_421m1 / lambda_aT3;	J[2][2] = J_440m1 / lambda_aL3;
}


double line_backtrack(double Ea, double PTa, double PLa, int Nparticles, double *Mass, double *Sign, double *Degeneracy, double * Baryon, double * Xcurrent, double * dX, double dX_abs, double g0, double * F)
{
	// This line backtracking algorithm is from the book Numerical Recipes in C

	// initial data for g(l) model:
	// g0 = f(Xcurrent)                // f at Xcurrent
	// f  = f(Xcurrent + dX)           // f at full newton step Xcurrent + dX
	// gprime0 = - 2g0                 // descent derivative at Xcurrent

	double X[3];

	for(int i = 0; i < 3; i++)
	{
		X[i] = Xcurrent[i] + dX[i];                 // default newton step
	}

	compute_F(Ea, PTa, PLa, Nparticles, Mass, Sign, Degeneracy, Baryon, X, F);	// update F at least once, default = F(Xcurrent + dX)

	double f = (F[0] * F[0]  +  F[1] * F[1]  +  F[2] * F[2]) / 2.;
	double gprime0 = - 2. * g0;

	double l = 1;                               	// default value for partial step parameter
	double alpha = 0.0001;                       	// descent rate

	double lroot, lprev, fprev;

	for(int n = 0; n < partial_backtracks; n++)		// line search iterations
	{
		if((l * dX_abs) <= tol_dX)                  // check if l.|dX| within desired tolerance
		{
			return l;
		}
		else if(f <= (g0  +  l * alpha * gprime0))  // check for sufficient decrease in f
		{
			return l;
		}
		else if(n == 0)                             // compute l (start with quadratic model)
		{
			lroot = - gprime0 / (2. * (f - g0 - gprime0));
		}
		else                                        // cubic model for subsequent iterations
		{
			// fixed bug on 3/25/20
			double a = ((f  -  g0  -  l * gprime0) / (l * l)  -  (fprev  -  g0  -  lprev * gprime0) / (lprev * lprev)) / (l - lprev);
			double b = (-lprev * (f  -  g0  -  l * gprime0) / (l * l)  +  l * (fprev  -  g0  -  lprev * gprime0)  /  (lprev * lprev)) / (l - lprev);

			if(a == 0)                              // quadratic solution to dg/dl = 0
			{
				lroot = - gprime0 / (2. * b);
			}
			else
			{
				double z = b * b  -  3. * a * gprime0;

				if(z < 0)
				{
					lroot = 0.5 * l;
				}
				else if(b <= 0)
				{
					lroot = (-b + sqrt(z)) / (3. * a);
				}
				else
				{
					lroot = - gprime0 / (b + sqrt(z));   // what does this mean?...
				}
			}

			lroot = fmin(lroot, 0.5 * l);
		}

		lprev = l;                                  // store current values for the next iteration
		fprev = f;

		// l = fmax(lroot, 0.1 * l);                   // update l and f
		l = fmax(lroot, 0.5 * l);                   // update l and f

		for(int i = 0; i < 3; i++)
		{
			X[i] = Xcurrent[i]  +  l * dX[i];
		}

		compute_F(Ea, PTa, PLa, Nparticles, Mass, Sign, Degeneracy, Baryon, X, F);

		f = (F[0] * F[0]  +  F[1] * F[1]  +  F[2] * F[2]) / 2.;
	}

	return l;
}


aniso_variables find_anisotropic_variables(double E, double pl, double pt, double lambda_0, double aT_0, double aL_0, int Nparticles, double *Mass, double *Sign, double *Degeneracy, double * Baryon)
{
	gsl_vector *x = gsl_vector_alloc(3);					// holds dX
    gsl_permutation *p = gsl_permutation_alloc(3);			// permutation vector

#ifndef ABORT_GSL
    gsl_set_error_handler_off();
#endif

	double Ea = E;		// kinetic energy density
	double PTa = pt;	// kinetic transverse pressure
	double PLa = pl;	// kinetic longitudinal pressure

	if(Ea < 0 || PTa < 0 || PLa < 0)
	{
	#ifdef FLAGS
		printf("find_anisotropic_variables flag: (E, pt, pl) = (%lf, %lf, %lf) is negative\n", Ea, PTa, PLa);
	#endif

		aniso_variables variables;
		variables.lambda = lambda_0;
		variables.aT = aT_0;
		variables.aL = aL_0;
		variables.did_not_find_solution = 1;
		variables.number_of_iterations = 0;

		return variables;
	}

	double X[3] = {lambda_0, aT_0, aL_0};					// current solution
	double dX[3];											// dX iteration
  	double F[3];											// F(X)
	double **J = (double**)malloc(3 * sizeof(double*));		// J(X)

	for(int i = 0; i < 3; i++)
	{
		J[i] = (double*)malloc(3 * sizeof(double));
 	}

 	compute_F(Ea, PTa, PLa, Nparticles, Mass, Sign, Degeneracy, Baryon, X, F);	// compue F

 	// double tolmin = 1.0e-6;								// tolerance for spurious convergence to local min of f = F.F/2 (what does this mean?)

	// scaled maximum step length allowed in line searches (is this ever updated?)
	double stepmax = 100. * fmax(sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]), 3.);


	for(int n = 0; n < N_max; n++)							// newton iteration loop
	{
		// compute J and f at X
	    compute_J(Ea, PTa, PLa, Nparticles, Mass, Sign, Degeneracy, Baryon, X, F, J);

	    double f = (F[0]*F[0] + F[1]*F[1] + F[2]*F[2]) / 2.;// f = F(X).F(X) / 2

	    double J_gsl[] = {J[0][0], J[0][1], J[0][2],
                  	  	  J[1][0], J[1][1], J[1][2],
                	  	  J[2][0], J[2][1], J[2][2]};		// Jacobian matrix in gsl format

        for(int i = 0; i < 3; i++)							// change sign of F
	    {
	    	F[i] *= -1.;
	    }

	    int s;
   		gsl_matrix_view A = gsl_matrix_view_array(J_gsl, 3, 3);
    	gsl_vector_view b = gsl_vector_view_array(F, 3);
       	gsl_linalg_LU_decomp(&A.matrix, p, &s);
       	gsl_linalg_LU_solve(&A.matrix, p, &b.vector, x);	// solve matrix equations J.dX = -F

       	for(int i = 0; i < 3; i++)
       	{
       		dX[i] = gsl_vector_get(x, i);					// get dX Newton iteration
       	}

	    double dX_abs = sqrt(dX[0]*dX[0] + dX[1]*dX[1] + dX[2]*dX[2]);	// l2 norm
	    // double dX_abs = fabs(dX[0]) + fabs(dX[1]) + fabs(dX[2]);		// l1 norm variant

		if(dX_abs > stepmax)
		{
			for(int i = 0; i < 3; i++)
			{
				dX[i] *= stepmax / dX_abs;					// rescale dX if too large
			}
			dX_abs = stepmax;
		}

		// compute partial step l and F(X + l.dX)
		double l = line_backtrack(Ea, PTa, PLa, Nparticles, Mass, Sign, Degeneracy, Baryon, X, dX, dX_abs, f, F);

		for(int i = 0; i < 3; i++)
	    {
	    	X[i] += (l * dX[i]);								// update solution X and convergence values
	    }

	    // printf("l = %lf\n\n", l);

	    double F_abs = sqrt(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]);	// update convergence values
	    //double F_abs = fabs(F[0]) + fabs(F[1]) + fabs(F[2]);	// not sure if l1norm is better

	    dX_abs *= l;

		if(X[0] < 0 || X[1] < 0 || X[2] < 0)				// check if any variable goes negative
		{
			aniso_variables variables;
			variables.lambda = lambda_0;
			variables.aT = aT_0;
			variables.aL = aL_0;
			variables.did_not_find_solution = 1;
			variables.number_of_iterations = n + 1;

			free_2D(J, 3);
			gsl_permutation_free(p);
       		gsl_vector_free(x);

			return variables;								// solution failed (unphysical)
		}
		else if(dX_abs <= tol_dX && F_abs <= tol_F)			// check for convergence
		{
			aniso_variables variables;
			variables.lambda = X[0];
			variables.aT = X[1];
			variables.aL = X[2];
			variables.did_not_find_solution = 0;
			variables.number_of_iterations = n + 1;

   			free_2D(J, 3);
			gsl_permutation_free(p);
   			gsl_vector_free(x);

			return variables;								// found solution
		}
	}	// newton iteration (n)

	aniso_variables variables;
	variables.lambda = lambda_0;							// try setting to T,1,1
	variables.aT = aT_0;
	variables.aL = aL_0;
	variables.did_not_find_solution = 1;
	variables.number_of_iterations = N_max;					// solution failed to converge (use previous guess)

	free_2D(J, 3);
	gsl_permutation_free(p);
  	gsl_vector_free(x);

	return variables;
}


famod_coefficient compute_famod_coefficient(double lambda, double aT, double aL, int Nparticles, double *Mass, double *Sign, double *Degeneracy, double * Baryon)
{
	famod_coefficient famod;

	double lambda2 = lambda * lambda;
	double aT2 = aT * aT;				// useful expressions
	double aL2 = aL * aL;
	double aT2_minus_aL2 = aT2 - aL2;
	double common_factor = aT2 * aL * lambda * lambda2 * lambda2 / four_pi2_hbarC3;

    double J_402m1 = 0;					// anisotropic integrals
    double J_421m1 = 0;

    for(int n = 0; n < Nparticles; n++)	// loop over hadrons in PDG
	{
		double mass = Mass[n];			// particle info from PDG arrays
		double sign = Sign[n];
		double degeneracy = Degeneracy[n];
		double baryon = Baryon[n];

		if(mass == 0)
		{
			continue;					// skip photons
		}

		double mbar = mass / lambda;
		double mbar2 = mbar * mbar;
		double chem = 0;				// replace w/ baryon * upsilonB (anisotropic version of baryon * alphaB)

    	double J_402m1_n = 0;			// individual hadron contribution to anisotropic integrals
    	double J_421m1_n = 0;

    	for(int i = 0; i < pbar_pts; i++)		// radial momentum integration loop
		{
			double pbar =     pbar_root_a3[i];	// pbar roots and weights for a = 3 (a = n + s)
			double weight = pbar_weight_a3[i];

			double pbar2 = pbar * pbar;
			double Ebar = sqrt(pbar2  +  mbar2);

			double w = sqrt(aL2  +  mbar2 / pbar2);
			double z = aT2_minus_aL2 / (w * w);
			double z2 =  z * z;

	  		double t_402;						// hypergeometric functions
	  		double t_421;

			if(z > delta)						// compute hypergeometric functions
			{
				double sqrtz = sqrt(z);
				double t = atan(sqrtz) / sqrtz;

				t_402 = (3. * (z - 1.)  +  (z * (3.*z - 2.) + 3.) * t) / (4. * z2);
				t_421 = (3.  +  z  +  (1. + z) * (z - 3.) * t) / (4. * z2);
			}
			else if(z < -delta && z > -1.)
			{
				double sqrtmz = sqrt(-z);
				double t = atanh(sqrtmz) / sqrtmz;

				t_402 = (3. * (z - 1.)  +  (z * (3.*z - 2.) + 3.) * t) / (4. * z2);
				t_421 = (3.  +  z  +  (1. + z) * (z - 3.) * t) / (4. * z2);
			}
			else if(fabs(z) <= delta)
			{
				double z3 = z2 * z;
				double z4 = z3 * z;
				double z5 = z4 * z;
				double z6 = z5 * z;

				t_402 = 1.0666666666666667 - 0.4571428571428572*z + 0.3047619047619048*z2 - 0.23088023088023088*z3 + 0.1864801864801865*z4 - 0.15664335664335666*z5 + 0.13514328808446457*z6;

				t_421 = 0.2666666666666666 - 0.0761904761904762*z + 0.0380952380952381*z2 - 0.023088023088023088*z3 + 0.015540015540015537*z4 - 0.011188811188811189*z5 + 0.00844645550527904*z6;
			}
			else
			{
			#ifdef FLAGS
				printf("compute_famod_coefficient flag: z = %lf is out of bounds\n", z);
			#endif
			}

			double quantum_stat = exp(Ebar + chem) + sign;
			double common_weight = weight * exp(pbar + Ebar) / (quantum_stat * quantum_stat);

			J_402m1_n += pbar2 / Ebar * common_weight * t_402 / w;
	    	J_421m1_n += pbar2 / Ebar * common_weight * t_421 / w;
		}

    	J_402m1_n *= degeneracy;		// multiply by degeneracy factor
    	J_421m1_n *= degeneracy;

    	J_402m1 += J_402m1_n;			// add individual hadron contribution
    	J_421m1 += J_421m1_n;
    }

	J_402m1 *= common_factor * aT2 * aT2 / 8.;	// multiply by common factor (and other things)
	J_421m1 *= common_factor * aT2 * aL2 / 2.;

	famod.betapiperp = J_402m1 / (aT2 * lambda);
	famod.betaWperp  = J_421m1 / (aT * aL * lambda);

	return famod;
}




