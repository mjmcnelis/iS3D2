#ifndef READINDATA_H
#define READINDATA_H

#include "iS3D.h"
#include "ParameterReader.h"
#include <fstream>

using namespace std;

class Gauss_Laguerre
{
  private:
    int alpha;   // # generalized powers (0 : alpha - 1)

  public:
    int points;  // # quadrature points
    double ** root;
    double ** weight;

    Gauss_Laguerre();
    void load_roots_and_weights(string file_name);
};

class Gauss_Legendre
{
  private:

  public:
    int points;  // # quadrature points
    double * root;
    double * weight;

    Gauss_Legendre();
    void load_roots_and_weights(string file_name);
};

class Plasma
{
  private:
    // I just wanted a cool name
  public:                               // average thermo quantities over freezeout surface
    double temperature;                 // GeV
    double energy_density;              // GeV / fm^3
    double pressure;                    // GeV / fm^3
    double baryon_chemical_potential;   // GeV
    double net_baryon_density;          // fm^-3

    Plasma();
    void load_thermodynamic_averages();
};

typedef struct
{
  long int mc_id; // Monte Carlo number according PDG
  string name;
  double mass;
  double width;
  int gspin; // spin degeneracy
  int baryon;
  int strange;
  int charm;
  int bottom;
  int gisospin; // isospin degeneracy
  int charge;
  int decays; // amount of decays listed for this resonance
  int stable; // defines whether this particle is considered as stable
  int decays_Npart[Maxdecaychannel];
  double decays_branchratio[Maxdecaychannel];
  int decays_part[Maxdecaychannel][Maxdecaypart];
  int sign; //Bose-Einstein or Dirac-Fermi statistics

  // ~ particle number / cell volume (for sampler routine)
  double equilibrium_density; // equilibrium density  (thermal number / u.dsigma)
  double bulk_density;        // bulk correction      (bulk number / u.dsigma / bulkPi)
  double diff_density;        // diffusion correction (diffusion number / V.dsigma)

} particle_info;

typedef struct
{
   double tau, x, y, eta;                 // contravariant spacetime position x^\mu
   double dat, dax, day, dan;             // covariant surface normal vector d\sigma_\mu
   double ux, uy, un;                     // contravariant fluid velocity u^\mu
   double E, T, P;                        // energy density E, temperature T and equilibrium pressure P
   double pixx, pixy, pixn, piyy, piyn;   // contravariant shear stress pi^\munu
   double bulkPi;                         // bulk viscous pressure Pi
   double muB, nB, Vx, Vy, Vn;            // net-baryon chemical potential muB, number density nB and contravariant diffusion current V^\mu
   double wtx, wty, wtn, wxy, wxn, wyn;   // contravariant thermal vorticity wbar^\mu\nu

   // double muE, muS; // electric and strange chemical potentials (might be needed in long run)
} FO_surf;

typedef struct
{
  // coefficients of Grad 14-moment approximation (vh)
  // df ~ ((c0-c2)m^2 + b.c1(u.p) + (4c2-c0)(u.p)^2).Pi + (b.c3 + c4(u.p))p_u.V^u + c5.p_u.p_v.pi^uv
  double c0;
  double c1;
  double c2;
  double c3;
  double c4;
  double shear14_coeff;

  // coefficients of RTA Chapman-Enskog expansion (vh)
  // df ~ ((c0-c2)m^2 + b.c1(u.p) + (4c2-c0)(u.p)^2).Pi + (b.c3 + c4(u.p))p_u.V^u + c5.p_u.p_v.pi^uv
  double F;
  double G;
  double betabulk;
  double betaV;
  double betapi;

  // PTB exclusive coefficients (vh)
  double lambda;
  double z;

  double delta_lambda;    // linearize lambda = 0 + delta_lambda
  double delta_z;         // linearize z = 1 + delta_z

} deltaf_coefficients;    // df coefficients for df corrections (vh)


const int max_digits = 10;  // max number of mcid digits (goes up to nuclei)

class read_mcid
{
  // determine the remaining particle properties based on mcid
  // (borrows some functionality from smash's pdgcode.hpp)

  // note: only have hadrons in mind (print errors if have leptons, etc)

  private:
    long int mcid;          // mcid of particle
  public:
    bool is_deuteron;       // label particle as deuteron (nothing yet for fake dibaryon d')
    bool is_hadron;         // label particle as hadron
    bool is_meson;          // label particle as meson
    bool is_baryon;         // label particle as baryon
    bool has_antiparticle;  // does particle have a distinct antiparticle?
    int baryon;             // baryon number
    int spin;               // spin x 2
    int gspin;              // spin degeneracy
    int sign;               // quantum statistics sign (BE, FD) = (-1, 1)

    uint32_t nJ  : 4;       // spin quantum number nJ = 2J + 1
    uint32_t nq3 : 4;       // third quark field
    uint32_t nq2 : 4;       // second quark field
    uint32_t nq1 : 4;       // first quark field, 0 for mesons
    uint32_t nL  : 4;       // "angular momentum"
    uint32_t nR  : 4;       // "radial excitation"
    uint32_t n   : 4, :3;   // first field: "counter"
    uint32_t n8;
    uint32_t n9;
    uint32_t n10;           // nuclei have 10-digits

    read_mcid(long int mcid_in);

    void is_particle_a_deuteron();                    // determine if particle is a deuteron
    void is_particle_a_hadron();                      // determine if particle is a hadron
    void is_particle_a_meson();                       // determine if particle is a meson
    void is_particle_a_baryon();                      // determine if the hadron is a baryon
    void get_baryon();                                // get the baryon number
    void get_spin();                                  // get the spin x 2
    void get_gspin();                                 // get the spin degeneracy
    void get_sign();                                  // get the quantum statistics sign
    void does_particle_have_distinct_antiparticle();  // determine if there's a distinct antiparticle

};


class FO_data_reader
{
    private:
        ParameterReader* paraRdr;
        int mode;                   // hydro code that constructed freezeout surface
        int dimension;              // dimension of freezeout surface
        int include_baryon;         // switch to include baryon chemical potential
        int number_of_cells;        // number of freezeout cells in freezeout surface file

    public:
        FO_data_reader(ParameterReader * paraRdr_in, string pathToInput);
        ~FO_data_reader();

        int get_number_cells();

        void read_freezeout_surface(FO_surf * surf_ptr);
        void read_surface_cpu_vh(FO_surf * surf_ptr);       // 1 (or 5 to include thermal vorticity)
        void read_surface_music(FO_surf* surf_ptr);         // 6
        void read_surface_hic_eventgen(FO_surf* surf_ptr);  // 7
};


class PDG_Data
{
  private:
    ParameterReader * paraRdr;
    int hrg_eos;
    string urqmd = "PDG/pdg-urqmd_v3.3+.dat"; // list of available pdg files
    string smash = "PDG/pdg_smash.dat";
    string smash_box = "PDG/pdg_box.dat";
  public:
    // read resonances from pdg file
    PDG_Data(ParameterReader * paraRdr_in);
    ~PDG_Data();

    int read_resonances_conventional(particle_info * particle, string pdg_filename);
    int read_resonances_smash_box(particle_info * particle, string pdg_filename);
    int read_resonances(particle_info * particle);
};

#endif
