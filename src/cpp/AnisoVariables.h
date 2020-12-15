
#ifndef ANISOVARIABLES_H_
#define ANISOVARIABLES_H_

const int N_max = 30;				// max number of iterations
const int partial_backtracks = 20;	// max number of partial backtracks (default = 20)
const double tol_dX = 1.e-4;		// tolerance error for dX (1.e-4)
const double tol_F = 1.e-4;			// tolerance error for F

const double delta = 0.01;		// piecewise interval where hypergeometric functions are taylor expanded

const int pbar_pts = 16;


// gauss-laguerre tables for anisotropic particle density:

const double pbar_root_a1[pbar_pts] =  {0.21614030523945225517,
										0.72638824325180395399,
										1.5335931603735413179,
										2.6449709986119109578,
										4.0709781608801906464,
										5.8258555151056045459,
										7.928504185306667092,
										10.403808289951039264,
										13.284661070707038206,
										16.615173216866612565,
										20.456006020027216858,
										24.893847025351910848,
										30.059862920202576072,
										36.170694543679177983,
										43.64036518417683703,
										53.529151160268420405};

const double pbar_weight_a1[pbar_pts] = {0.063277332879539430926,
										0.23109046152071973177,
										0.31693354216399909707,
										0.23789421787524460737,
										0.11027274335924788215,
										0.033089588356396983661,
										0.0065257164077924663047,
										0.00084271143679204720664,
										0.000069961909555681874201,
										3.6120860421335547395e-6,
										1.1015412514742446856e-7,
										1.8356825127325042141e-9,
										1.4814888658850781417e-11,
										4.7348453590396131986e-14,
										4.0841460550022559944e-17,
										3.6239549624789140618e-21};

// gauss-laguerre tables for F calculation:

const double pbar_root_a2[pbar_pts] =  {0.37761350834474073436,
										1.0174919576025657044,
										1.947758020424243826,
										3.1769272448898686839,
										4.7162400697917956248,
										6.5805882657749125043,
										8.7894652706470741271,
										11.368323082833318905,
										14.350626727437044399,
										17.781095724841646049,
										21.721084796571308914,
										26.258138675111068012,
										31.524596004275818683,
										37.738921002528939074,
										45.318546110089842558,
										55.3325835388358122};

const double pbar_weight_a2[pbar_pts] = {0.04860640946707870253,
										0.29334739019044023149,
										0.58321936338355098941,
										0.5818741485961734789,
										0.33818053747379021731,
										0.12210596394497989094,
										0.028114625800663730181,
										0.004143149192482256704,
										0.00038564853376743830033,
										0.000022015800563109071234,
										7.3423624381565175075e-7,
										1.3264604420480413669e-8,
										1.1526664829084294661e-10,
										3.9470691512460869719e-13,
										3.6379782563605336032e-16,
										3.4545761231361240027e-20};

// gauss-laguerre tables for J calculation:

const double pbar_root_a3[pbar_pts] =   {0.56744345899157414477,
											1.3329077327598933409,
											2.3814824800700584858,
											3.7238266420934269684,
											5.3721239521618709961,
											7.3419366282613523442,
											9.6533321372612290461,
											12.332301407018237936,
											15.412850065407715512,
											18.940275575860272298,
											22.976595715601869244,
											27.610181447426074872,
											32.974509203258495494,
											39.289823253364147934,
											46.976896276710312577,
											57.113514023753468805};

const double pbar_weight_a3[pbar_pts] = {0.065098112100944931481,
											0.56527322423640238298,
											1.4898931385790650188,
											1.8612470448965348651,
											1.3004755484827235985,
											0.54781964685691509622,
											0.14384304320810570104,
											0.023749847242162275954,
											0.0024424446935061086306,
											0.00015234059955885076823,
											5.5012653050135610492e-6,
											1.0684262964359744438e-7,
											9.9252433589722221362e-10,
											3.6186334278007682335e-12,
											3.5437327807321532759e-15,
											3.5818035528777921894e-19};


typedef struct
{
	double lambda;
	double aT;
	double aL;
	bool did_not_find_solution;
	int number_of_iterations;

} aniso_variables;

typedef struct
{
	double betapiperp;
	double betaWperp;

} famod_coefficient;

aniso_variables find_anisotropic_variables(double E, double pl, double pt, double lambda_0, double aT_0, double aL_0, int Nparticles, double *Mass, double *Sign, double *Degeneracy, double *Baryon);

famod_coefficient compute_famod_coefficient(double lambda, double aT, double aL, int Nparticles, double *Mass, double *Sign, double *Degeneracy, double *Baryon);


#endif




