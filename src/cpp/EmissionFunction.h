#ifndef EMISSIONFUNCTION_H
#define EMISSIONFUNCTION_H

#include<string>
#include<vector>
#include<random>
#include "Table.h"
#include "iS3D.h"
#include "ParameterReader.h"
#include "DeltafData.h"
#include "SampledParticle.h"
#include "LocalRestFrame.h"

using namespace std;


// thermal particle density (just for crosschecking)
//double equilibrium_particle_density(double mass, double degeneracy, double sign, double T, double chem);

double compute_detA(Shear_Stress pimunu, double betapi, double bulk_mod);

bool is_linear_pion0_density_negative(double T, double neq_pion0, double J20_pion0, double bulkPi, double F, double betabulk);

bool does_feqmod_breakdown(double mass_pion0, double T, double F, double bulkPi, double betabulk, double detA, double detA_min, double z, Gauss_Laguerre * laguerre, int df_mode, int fast, double Tavg, double F_avg, double betabulk_avg);


class EmissionFunctionArray
{
private:
  ParameterReader* paraRdr;

  long CORES;   // number of cores for openmp

  int OPERATION; // calculate smooth spectra or sample distributions
  int MODE; //vh or vah , ...

  int DF_MODE;  // delta-f type
  string df_correction;

  int DIMENSION; // hydro d+1 dimensions (2+1 or 3+1)
  int INCLUDE_BULK_DELTAF;
  int INCLUDE_SHEAR_DELTAF;
  int INCLUDE_BARYONDIFF_DELTAF;

  int REGULATE_DELTAF;
  int OUTFLOW;

  int INCLUDE_BARYON;
  double DETA_MIN;
  int GROUP_PARTICLES;
  double PARTICLE_DIFF_TOLERANCE;

  double MASS_PION0;

  int LIGHTEST_PARTICLE; //mcid of lightest resonance to calculate in decay feed-down
  int DO_RESONANCE_DECAYS; // smooth resonance decays option

  int OVERSAMPLE; // whether or not to iteratively oversample surface
  int FAST;                 // switch to compute mean hadron number quickly using an averaged (T,muB)
  double MIN_NUM_HADRONS; //min number of particles summed over all samples
  double MAX_NUM_SAMPLES; // max number of events sampled
  long int SAMPLER_SEED; //the seed for the particle sampler. If chosen < 0, seed set with clocktime

  int TEST_SAMPLER;

  long Nevents = 1;                  // default number of sampled events

  // for binning sampled particles (for sampler tests)
  double PT_MIN;
  double PT_MAX;
  int PT_BINS;
  double PT_WIDTH;

  double Y_CUT;
  int Y_BINS;
  double Y_WIDTH;

  int PHIP_BINS;
  double PHIP_WIDTH;

  double ETA_CUT;
  int ETA_BINS;
  double ETA_WIDTH;

  double TAU_MIN;
  double TAU_MAX;
  int TAU_BINS;
  double TAU_WIDTH;

  double R_MIN;
  double R_MAX;
  int R_BINS;
  double R_WIDTH;

  // for sampler test (2+1d)
  double **dN_dy_count;          // event-averaged momentum distributions
  double **dN_2pipTdpTdy_count;
  double **dN_dphipdy_count;

  double ***vn_real_count;      // event-averaged Vn's
  double ***vn_imag_count;
  const int K_MAX = 7;          // {v1, ..., v7}
  double **pT_count;            // count in each pT bin

  double **dN_deta_count;       // event-averaged spacetime distribution
  double **dN_taudtaudy_count;
  double **dN_twopirdrdy_count;
  double **dN_dphisdy_count;


  Table *pT_tab, *phi_tab, *y_tab, *eta_tab;
  long pT_tab_length, phi_tab_length, y_tab_length, eta_tab_length;
  long FO_length;
  double *dN_pTdpTdphidy; //to hold smooth CF 3D spectra of all species
  double *logdN_PTdPTdPhidY; // hold log of smooth CF 3D spectra of parent (set in res decay for linear interpolation)

  double *St, *Sx, *Sy, *Sn; //to hold the polarization vector of all species
  double *Snorm; //the normalization of the polarization vector of all species

  std::vector<Sampled_Particle> particle_list;                        // to hold sampled particle list (inactive)
  std::vector< std::vector<Sampled_Particle> > particle_event_list;   // holds sampled particle list of all events

  int *chosen_particles_01_table;       // has length Nparticle, 0 means miss, 1 means include
  int *chosen_particles_sampling_table; // store particle index; the sampling process follows the order specified by this table
  int Nparticles;
  int number_of_chosen_particles;
  particle_info* particles;       // contains all the particle info from pdg.dat
  FO_surf* surf_ptr;
  Deltaf_Data * df_data;
  bool particles_are_the_same(int, int);

public:

  // constructor
  EmissionFunctionArray(ParameterReader* paraRdr_in, Table* chosen_particle, Table* pT_tab_in, Table* phi_tab_in, Table* y_tab_in, Table* eta_tab_in, particle_info* particles_in, int Nparticles, FO_surf* FOsurf_ptr_in, long FO_length_in, Deltaf_Data * df_data_in);

  ~EmissionFunctionArray();

  // main function
  void calculate_spectra(std::vector<std::vector<Sampled_Particle>> &particle_event_list_in);


  // continuous spectra routines:
  //:::::::::::::::::::::::::::::::::::::::::::::::::

  // continuous spectra with feq + df
  void calculate_dN_pTdpTdphidy(double *Mass, double *Sign, double *Degeneracy, double *Baryon, double *T_fo, double *P_fo, double *E_fo, double *tau_fo, double *eta_fo, double *ux_fo, double *uy_fo, double *un_fo, double *dat_fo, double *dax_fo, double *day_fo, double *dan_fo, double *pixx_fo, double *pixy_fo, double *pixn_fo, double *piyy_fo, double *piyn_fo, double *bulkPi_fo, double *muB_fo, double *nB_fo, double *Vx_fo, double *Vy_fo, double *Vn_fo, Deltaf_Data *df_data);

  // continuous spectra with feqmod
  void calculate_dN_pTdpTdphidy_feqmod(double *Mass, double *Sign, double *Degeneracy, double *Baryon, double *T_fo, double *P_fo, double *E_fo, double *tau_fo, double *eta_fo, double *ux_fo, double *uy_fo, double *un_fo, double *dat_fo, double *dax_fo, double *day_fo, double *dan_fo, double *pixx_fo, double *pixy_fo, double *pixn_fo, double *piyy_fo, double *piyn_fo, double *bulkPi_fo, double *muB_fo, double *nB_fo, double *Vx_fo, double *Vy_fo, double *Vn_fo, Gauss_Laguerre * laguerre, Deltaf_Data * df_data);

  void calculate_dN_pTdpTdphidy_famod(double *Mass, double *Sign, double *Degeneracy, double *Baryon, double *T_fo, double *P_fo, double *E_fo, double *tau_fo, double *eta_fo, double *ux_fo, double *uy_fo, double *un_fo, double *dat_fo, double *dax_fo, double *day_fo, double *dan_fo, double *pixx_fo, double *pixy_fo, double *pixn_fo, double *piyy_fo, double *piyn_fo, double *bulkPi_fo, double *muB_fo, double *nB_fo, double *Vx_fo, double *Vy_fo, double *Vn_fo, int Nparticles, double *Mass_PDG, double *Sign_PDG, double *Degeneracy_PDG, double *Baryon_PDG);


  void calculate_dN_dX(int *MCID, double *Mass, double *Sign, double *Degeneracy, double *Baryon,
  double *T_fo, double *P_fo, double *E_fo, double *tau_fo, double *x_fo, double *y_fo, double *eta_fo, double *ux_fo, double *uy_fo, double *un_fo,
  double *dat_fo, double *dax_fo, double *day_fo, double *dan_fo,
  double *pixx_fo, double *pixy_fo, double *pixn_fo, double *piyy_fo, double *piyn_fo, double *bulkPi_fo,
  double *muB_fo, double *nB_fo, double *Vx_fo, double *Vy_fo, double *Vn_fo, Deltaf_Data *df_data);

   void calculate_dN_dX_feqmod(int *MCID, double *Mass, double *Sign, double *Degeneracy, double *Baryon, double *T_fo, double *P_fo, double *E_fo, double *tau_fo, double *x_fo, double *y_fo, double *eta_fo, double *ux_fo, double *uy_fo, double *un_fo, double *dat_fo, double *dax_fo, double *day_fo, double *dan_fo, double *pixx_fo, double *pixy_fo, double *pixn_fo, double *piyy_fo, double *piyn_fo, double *bulkPi_fo, double *muB_fo, double *nB_fo, double *Vx_fo, double *Vy_fo, double *Vn_fo, Gauss_Laguerre * laguerre, Deltaf_Data * df_data);



  //:::::::::::::::::::::::::::::::::::::::::::::::::


  // sampling spectra routines:
  //:::::::::::::::::::::::::::::::::::::::::::::::::

  // calculate average total particle yield from freezeout surface to determine number of events to sample
  double calculate_total_yield(double * Equilibrium_Density, double * Bulk_Density, double * Diffusion_Density, double *T_fo, double *P_fo, double *E_fo, double *tau_fo, double *ux_fo, double *uy_fo, double *un_fo, double *dat_fo, double *dax_fo, double *day_fo, double *dan_fo, double *pixx_fo, double *pixy_fo, double *pixn_fo, double *piyy_fo, double *piyn_fo, double *bulkPi_fo, double *muB, double *nB, double *Vx_fo, double *Vy_fo, double *Vn_fo, Deltaf_Data * df_data, Gauss_Laguerre * laguerre);

  // sample particles with feq + df14, feq + dfCE, PTM feqmod or PTB feqmod
  void sample_dN_pTdpTdphidy(double *Mass, double *Sign, double *Degeneracy, double *Baryon, int *MCID, double *Equilibrium_Density, double *Bulk_Density, double *Diffusion_Density, double *T_fo, double *P_fo, double *E_fo, double *tau_fo, double *x_fo, double *y_fo, double *eta_fo, double *ux_fo, double *uy_fo, double *un_fo, double *dat_fo, double *dax_fo, double *day_fo, double *dan_fo, double *pixx_fo, double *pixy_fo, double *pixn_fo, double *piyy_fo, double *piyn_fo, double *bulkPi_fo, double *muB_fo, double *nB_fo, double *Vx_fo, double *Vy_fo, double *Vn_fo, Deltaf_Data *df_data, Gauss_Laguerre * laguerre, Gauss_Legendre * legendre);


  // sample particles with fa or PTM famod
  void sample_dN_pTdpTdphidy_famod(double *Mass, double *Sign, double *Degeneracy, double *Baryon, int *MCID, double *T_fo, double *P_fo, double *E_fo, double *tau_fo, double *x_fo, double *y_fo, double *eta_fo, double *ux_fo, double *uy_fo, double *un_fo, double *dat_fo, double *dax_fo, double *day_fo, double *dan_fo, double *pixx_fo, double *pixy_fo, double *pixn_fo, double *piyy_fo, double *piyn_fo, double *bulkPi_fo, double *muB_fo, double *nB_fo, double *Vx_fo, double *Vy_fo, double *Vn_fo, int Nparticles, double *Mass_PDG, double *Sign_PDG, double *Degeneracy_PDG, double *Baryon_PDG);


  // add counts for sampled distributions
  void sample_dN_dy(int chosen_index, double y);
  void sample_dN_deta(int chosen_index, double eta);
  void sample_dN_dphipdy(int chosen_index, double px, double py);
  void sample_dN_2pipTdpTdy(int chosen_index, double px, double py);
  void sample_vn(int chosen_index, double px, double py);
  void sample_dN_dX(int chosen_index, double tau, double x, double y);


  //:::::::::::::::::::::::::::::::::::::::::::::::::

  // spin polarization:
  void calculate_spin_polzn(double *Mass, double *Sign, double *Degeneracy,
  double *tau_fo, double *eta_fo, double *ux_fo, double *uy_fo, double *un_fo,
  double *dat_fo, double *dax_fo, double *day_fo, double *dan_fo,
  double *wtx_fo, double *wty_fo, double *wtn_fo, double *wxy_fo, double *wxn_fo, double *wyn_fo, Plasma * QGP);


  // write to file functions:
  //:::::::::::::::::::::::::::::::::::::::::::::::::

  void write_dN_pTdpTdphidy_toFile(int *MCID); // write invariant 3D spectra to file
  void write_dN_dphidy_toFile(int *MCID);
  void write_dN_twopipTdpTdy_toFile(int *MCID);
  void write_dN_dy_toFile(int *MCID);
  void write_continuous_vn_toFile(int *MCID);
  void write_polzn_vector_toFile(); //write components of spin polarization vector to file

  void write_particle_list_toFile();              // write sampled particle list
  void write_particle_list_OSC();                 // write sampled particle list in OSCAR format for UrQMD/SMASH

  // for sampler test
  void write_sampled_dN_dy_to_file_test(int * MCID);
  void write_sampled_dN_deta_to_file_test(int * MCID);
  void write_sampled_dN_2pipTdpTdy_to_file_test(int * MCID);
  void write_sampled_dN_dphipdy_to_file_test(int * MCID);
  void write_sampled_vn_to_file_test(int * MCID);
  void write_sampled_dN_dX_to_file_test(int * MCID);

};

#endif
