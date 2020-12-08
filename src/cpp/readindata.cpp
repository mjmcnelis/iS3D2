
#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<stdlib.h>

#include "iS3D.h"
#include "Macros.h"
#include "readindata.h"
#include "arsenal.h"
#include "ParameterReader.h"
#include "Table.h"
#include "gaussThermal.h"

using namespace std;


Gauss_Laguerre::Gauss_Laguerre()
{

}


void Gauss_Laguerre::load_roots_and_weights(string file_name)
{
  stringstream file;
  file << file_name;

  FILE * gauss_file = fopen(file.str().c_str(), "r");
  if(gauss_file == NULL) printf("Error: couldn't open gauss laguerre file\n");
  // get the powers and number of gauss points
  fscanf(gauss_file, "%d\t%d", &alpha, &points);
  // allocate memory for the roots and weights
  root = (double **)calloc(alpha, sizeof(double));
  weight = (double **)calloc(alpha, sizeof(double));

  for(int i = 0; i < alpha; i++)
  {
    root[i] = (double *)calloc(points, sizeof(double));
    weight[i] = (double *)calloc(points, sizeof(double));
  }

  int dummy;  // dummy alpha index

  // load the arrays
  for(int i = 0; i < alpha; i++)
  {
    for(int j = 0; j < points; j++)
    {
      fscanf(gauss_file, "%d\t%lf\t%lf", &dummy, &root[i][j], &weight[i][j]);
    }
  }
  fclose(gauss_file);
}


Gauss_Legendre::Gauss_Legendre()
{

}


void Gauss_Legendre::load_roots_and_weights(string file_name)
{
  stringstream file;
  file << file_name;

  FILE * gauss_file = fopen(file.str().c_str(), "r");

  if(gauss_file == NULL) printf("Error: couldn't open gauss legendre file\n");

  fscanf(gauss_file, "%d", &points);

  // allocate memory for the roots and weights
  root = (double *)calloc(points, sizeof(double));
  weight = (double *)calloc(points, sizeof(double));

  // load the arrays
  for(int i = 0; i < points; i++)
  {
    fscanf(gauss_file, "%lf\t%lf", &root[i], &weight[i]);
  }

  fclose(gauss_file);
}


Plasma::Plasma()
{

}


void Plasma::load_thermodynamic_averages()
{
  // reads the averaged thermodynamic quantities of the freezeout surface
  // ** assumes the file has already been created in read_surf_switch

  FILE * thermodynamic_file = fopen("average_thermodynamic_quantities.dat", "r");

  if(thermodynamic_file == NULL) printf("Error opening average thermodynamic file\n");
  fscanf(thermodynamic_file, "%lf\n%lf\n%lf\n%lf\n%lf", &temperature, &energy_density, &pressure, &baryon_chemical_potential, &net_baryon_density);
  fclose(thermodynamic_file);
}


FO_data_reader::FO_data_reader(ParameterReader* paraRdr_in, string path_in)
{
  paraRdr = paraRdr_in;
  mode = paraRdr->getVal("mode");                         // change name to hydro_code
  dimension = paraRdr->getVal("dimension");
  include_baryon = paraRdr->getVal("include_baryon");
}


FO_data_reader::~FO_data_reader()
{

}


int FO_data_reader::get_number_cells()
{
  ostringstream surface_file;
  surface_file << "input/surface.dat";
  Table block_file(surface_file.str().c_str());

  number_of_cells = block_file.getNumberOfRows();

  return number_of_cells;
}


void FO_data_reader::read_freezeout_surface(FO_surf* surf_ptr)
{

  if(mode == 1 || mode == 5)
  {
    read_surface_cpu_vh(surf_ptr);                    // read 2+1d or 3+1d surface file from cpu vh (or cpu vah)
  }
  else if (mode == 6)
  {
    read_surface_music(surf_ptr);                     // read 2+1d or 3+1d surface file from MUSIC (public version)
  }
  else if (mode == 7)
  {
    read_surface_hic_eventgen(surf_ptr);              // read 2+1d surface file from HIC-EventGen
  }
}


void FO_data_reader::read_surface_cpu_vh(FO_surf* surf_ptr)
{
  printf("from input/surface.dat and undoing hbarc = 1 units...");
  if(mode == 5)
  {
    printf(" (includes thermal vorticity)");          // only Derek's version of cpu vh outputs thermal vorticity wbar^\munu
  }

  printf("\n\nHydrodynamic code = CPU VH (or CPU VAH)\n\n");
  printf("\thttps://github.com/derekeverett/cpu-vh\t(CPU VH)\n");
  printf("\thttps://github.com/mjmcnelis/cpu_vah\t(CPU VAH)\n\n");

  printf("Please check that input/surface.dat has the following format (and 1 blank line eof):\n\n\t");

  if(include_baryon)
  {
    if(mode == 5)
    {
      printf("[t x y n ds_t ds_x ds_y ds_n u^x u^y u^n E T P pi^xx pi^xy pi^xn pi^yy pi^yn Pi muB nB V^x V^y V^n wbar^tx wbar^ty wbar^tn wbar^xy wbar^xn wbar^yn]\n\n");
    }
    else
    {
      printf("[t x y n ds_t ds_x ds_y ds_n u^x u^y u^n E T P pi^xx pi^xy pi^xn pi^yy pi^yn Pi muB nB V^x V^y V^n]\n\n");
    }
  }
  else
  {
    if(mode == 5)
    {
      printf("[t x y n ds_t ds_x ds_y ds_n u^x u^y u^n E T P pi^xx pi^xy pi^xn pi^yy pi^yn Pi wbar^tx wbar^ty wbar^tn wbar^xy wbar^xn wbar^yn]\n\n");
    }
    else
    {
      printf("[t x y n ds_t ds_x ds_y ds_n u^x u^y u^n E T P pi^xx pi^xy pi^xn pi^yy pi^yn Pi]\n\n");
    }
  }


  ostringstream surfdat_stream;                       // prepare to read in surface.dat
  surfdat_stream << "input/surface.dat";
  ifstream surfdat(surfdat_stream.str().c_str());
  double dummy;


  double T_avg = 0;                                   // average thermodynamic variables across freezeout surface
  double E_avg = 0;
  double P_avg = 0;
  double muB_avg = 0;
  double nB_avg = 0;
  double max_volume = 0;                              // max volume of freezeout surface


  for(long i = 0; i < number_of_cells; i++)           // loop over freezeout cells
  {
    // contravariant spacetime position x^\mu
    surfdat >> surf_ptr[i].tau;                       // \tau [fm]
    surfdat >> surf_ptr[i].x;                         // x [fm]
    surfdat >> surf_ptr[i].y;                         // y [fm]
    surfdat >> surf_ptr[i].eta;                       // \eta_s [1]


    // covariant surface normal vector d\sigma_\mu
    surfdat >> surf_ptr[i].dat;                       // d\sigma_\tau [fm^-2]
    surfdat >> surf_ptr[i].dax;                       // d\sigma_x [fm^-2]
    surfdat >> surf_ptr[i].day;                       // d\sigma_y [fm^-2]
    surfdat >> surf_ptr[i].dan;                       // d\sigma_\eta [fm^-1]


    // contravariant fluid velocity u^\mu
    surfdat >> surf_ptr[i].ux;                        // u^x [1]
    surfdat >> surf_ptr[i].uy;                        // u^y [1]
    surfdat >> surf_ptr[i].un;                        // u^\eta [fm^-1]


    // thermodynamic variables
    surfdat >> dummy;                                 // energy density [fm^-4]
    double E = dummy * hbarC;                         // convert to [GeV/fm^3]
    surf_ptr[i].E = E;

    surfdat >> dummy;                                 // temperature [fm^-1]
    double T = dummy * hbarC;                         // convert to [GeV]
    surf_ptr[i].T = T;

    surfdat >> dummy;                                 // equilibrium pressure [fm^-4]
    double P = dummy * hbarC;                         // convert to [GeV/fm^3]
    surf_ptr[i].P = P;


    // contravariant shear stress pi^\munu
    surfdat >> dummy;                                 // pi^xx [fm^-4]
    surf_ptr[i].pixx = dummy * hbarC;                 // convert to [GeV/fm^3]

    surfdat >> dummy;                                 // pi^xy [fm^-4]
    surf_ptr[i].pixy = dummy * hbarC;                 // pi^x\eta [fm^-5]

    surfdat >> dummy;                                 // pi^x\eta [fm^-5]
    surf_ptr[i].pixn = dummy * hbarC;                 // convert to [GeV/fm^4]

    surfdat >> dummy;                                 // pi^yy [fm^-4]
    surf_ptr[i].piyy = dummy * hbarC;                 // convert to [GeV/fm^3]

    surfdat >> dummy;                                 // pi^y\eta [fm^-5]
    surf_ptr[i].piyn = dummy * hbarC;                 // convert to [GeV/fm^4]


    // bulk viscous pressure
    surfdat >> dummy;                                 // Pi [fm^-4]
    surf_ptr[i].bulkPi = dummy * hbarC;               // convert to [GeV/fm^3]


    double muB = 0;                                   // default value for net-baryon chemical potential
    double nB = 0;                                    // default value for net-baryon density

    if(include_baryon)
    {
      // net-baryon chemical potential and density
      surfdat >> dummy;                               // muB [fm^-1]
      muB = dummy * hbarC;                            // convert to [GeV]
      surf_ptr[i].muB = muB;

      surfdat >> nB;                                  // nB[fm^-3]
      surf_ptr[i].nB = nB;


      // contravariant baryon diffusion V^\mu
      surfdat >> surf_ptr[i].Vx;                      // V^x [fm^-3]
      surfdat >> surf_ptr[i].Vy;                      // V^y [fm^-3]
      surfdat >> surf_ptr[i].Vn;                      // V^\eta [fm^-4] (12/2/20 check this --> fixed units on 10/8/18)
    }


    // thermal vorticity wbar^\mu\nu
    if(mode == 5)                                     // contravariant, dimensionless? undo hbarc = 1?
    {
      surfdat >> surf_ptr[i].wtx;                     // ask Derek for definition and units (any conversion?)
      surfdat >> surf_ptr[i].wty;                     // all upper indices?
      surfdat >> surf_ptr[i].wtn;
      surfdat >> surf_ptr[i].wxy;
      surfdat >> surf_ptr[i].wxn;
      surfdat >> surf_ptr[i].wyn;
    }


    // check whether 2+1d freezeout cells are really boost-invariant
    if(dimension == 2)
    {
      if(surf_ptr[i].eta != 0)
      {
      #ifdef FLAGS
        printf("read_surface_cpu_vh flag: setting spacetime rapidity of boost-invariant freezeout cell to eta = 0\n");
      #endif

        surf_ptr[i].eta = 0;
      }
      if(surf_ptr[i].dan != 0 || surf_ptr[i].un != 0 || surf_ptr[i].pixn != 0 || surf_ptr[i].piyn != 0)
      {
      #ifdef FLAGS
        printf("read_surface_cpu_vh flag: dimension = 2 but freezeout cell %ld is not boost-invariant (please check format in surface.dat)\n", i);
      #endif
      }
    }


    // compute the averaged thermodynamic quantities (for fast df coefficients)
    double tau = surf_ptr[i].tau;
    double tau2 = tau * tau;
    double ux = surf_ptr[i].ux;
    double uy = surf_ptr[i].uy;
    double un = surf_ptr[i].un;
    double ut = sqrt(1.  +  ux * ux  +  uy * uy  +  tau2 * un * un);
    double dat = surf_ptr[i].dat;
    double dax = surf_ptr[i].dax;
    double day = surf_ptr[i].day;
    double dan = surf_ptr[i].dan;

    double uds = ut * dat  +  ux * dax  +  uy * day  +  un * dan;                 // u^\mu . d\sigma_\mu
    double ds_ds = dat * dat  -  dax * dax  -  day * day  -  dan * dan / tau2;    // d\sigma^\mu . d\sigma_\mu
    double ds_max = fabs(uds)  +  sqrt(fabs(uds * uds  -  ds_ds));                // max volume element |ds|

    max_volume += ds_max;         // append values
    E_avg += (E * ds_max);
    T_avg += (T * ds_max);
    P_avg += (P * ds_max);
    muB_avg += (muB * ds_max);
    nB_avg += (nB * ds_max);
  }

  surfdat.close();

  T_avg /= max_volume;            // divide by total max volume
  E_avg /= max_volume;
  P_avg /= max_volume;
  muB_avg /= max_volume;
  nB_avg /= max_volume;


  // write averaged thermodynamic variables to file (what happens if read from memory again?)
  ofstream thermal_average("average_thermodynamic_quantities.dat", ios_base::out);
  thermal_average << setprecision(15) << T_avg << "\n" << E_avg << "\n" << P_avg << "\n" << muB_avg << "\n" << nB_avg;
  thermal_average.close();
}




void FO_data_reader::read_surface_music(FO_surf* surf_ptr)
{
  printf("from input/surface.dat and undoing hbarc = 1 units and tau factors...\n");

  printf("\nHydrodynamic code = MUSIC (public version)\n\n");
  printf("\thttps://github.com/MUSIC-fluid/MUSIC\n\n");

  printf("Please check that input/surface.dat has the following format (and 1 blank line eof):\n\n\t");

  if(include_baryon)
  {
    printf("[t x y n ds_t/t ds_x/t ds_y/t ds_n/t u^t u^x u^y t.u^n E T muB muS muC (E+P)/T pi^tt pi^tx pi^ty t.pi^tn pi^xx pi^xy t.pi^xn pi^yy t.pi^yn t2.pi^nn Pi nB V^t V^x V^y t.V^n]\n\n");

    // muB, muS, muC = baryon, strange, charm chemical potentials
  }
  else
  {
    printf("[t x y n ds_t/t ds_x/t ds_y/t ds_n/t u^t u^x u^y t.u^n E T muB muS muC (E+P)/T pi^tt pi^tx pi^ty t.pi^tn pi^xx pi^xy t.pi^xn pi^yy t.pi^yn t2.pi^nn Pi]\n\n");
  }


  ostringstream surfdat_stream;                       // prepare to read in surface.dat
  surfdat_stream << "input/surface.dat";
  ifstream surfdat(surfdat_stream.str().c_str());
  double dummy;


  double T_avg = 0;                                   // average thermodynamic variables across freezeout surface
  double E_avg = 0;
  double P_avg = 0;
  double muB_avg = 0;
  double nB_avg = 0;
  double max_volume = 0;                              // max volume of freezeout surface


  for(long i = 0; i < number_of_cells; i++)           // loop over freezeout cells
  {
    // contravariant spacetime position x^\mu
    surfdat >> dummy;                                 // \tau [fm]
    double tau = dummy;
    surf_ptr[i].tau = tau;

    surfdat >> surf_ptr[i].x;                         // x [fm]
    surfdat >> surf_ptr[i].y;                         // y [fm]
    surfdat >> surf_ptr[i].eta;                       // \eta_s [1]


    // covariant surface normal vector d\sigma_\mu / \tau
    surfdat >> dummy;                                 // d\sigma_\tau / \tau [fm^-3]
    surf_ptr[i].dat = dummy * tau;                    // multiply by \tau

    surfdat >> dummy;                                 // d\sigma_x / \tau [fm^-3]
    surf_ptr[i].dax = dummy * tau;                    // multiply by \tau

    surfdat >> dummy;                                 // d\sigma_y / \tau [fm^-3]
    surf_ptr[i].day = dummy * tau;                    // multiply by \tau

    surfdat >> dummy;                                 // d\sigma_\eta / \tau [fm^-4]
    surf_ptr[i].dan = dummy * tau;                    // multiply by \tau


    // contravariant fluid velocity u^\mu
    surfdat >> dummy;                                 // u^\tau [1]
    surfdat >> surf_ptr[i].ux;                        // u^x [1]
    surfdat >> surf_ptr[i].uy;                        // u^y [1]

    surfdat >> dummy;                                 // \tau . u^\eta [1]
    surf_ptr[i].un = dummy / tau;                     // divide by \tau


    // thermodynamic variables
    surfdat >> dummy;                                 // energy density [fm^-4]
    double E = dummy * hbarC;                         // convert to [GeV/fm^3]
    surf_ptr[i].E = E;

    surfdat >> dummy;                                 // temperature [fm^-1]
    double T = dummy * hbarC;                         // convert to [GeV]
    surf_ptr[i].T = T;

    surfdat >> dummy;                                 // net-baryon chemical potential [fm^-1]
    double muB = dummy * hbarC;                       // convert to [GeV]
    surf_ptr[i].muB = muB;

    surfdat >> dummy;                                 // strange chemical potential (units?)
    surfdat >> dummy;                                 // charm chemical potential (these don't seem to be used here...)

    surfdat >> dummy;                                 // (E + P) / T [fm^-3]
    double P = dummy * T  -  E;                       // equilibrium pressure [GeV/fm^3]
    surf_ptr[i].P = P;


    // contravariant shear stress pi^\munu
    surfdat >> dummy;                                 // pi^\tau\tau [fm^-4]
    surfdat >> dummy;                                 // pi^\taux [fm^-4]
    surfdat >> dummy;                                 // pi^\tauy [fm^-4]
    surfdat >> dummy;                                 // tau . pi^\tau\eta [fm^-4]

    surfdat >> dummy;                                 // pi^xx [fm^-4]
    surf_ptr[i].pixx = dummy * hbarC;                 // convert to [GeV/fm^3]

    surfdat >> dummy;                                 // pi^xy [fm^-4]
    surf_ptr[i].pixy = dummy * hbarC;                 // convert to [GeV/fm^3]

    surfdat >> dummy;                                 // \tau . pi^x\eta [fm^-4]
    surf_ptr[i].pixn = dummy * hbarC / tau;           // convert to [GeV/fm^4] (divided by \tau)

    surfdat >> dummy;                                 // pi^yy [fm^-4]
    surf_ptr[i].piyy = dummy * hbarC;                 // convert to [GeV/fm^3]

    surfdat >> dummy;                                 // \tau . pi^y\eta [fm^-4]
    surf_ptr[i].piyn = dummy * hbarC / tau;           // convert to [GeV/fm^4] (divided by \tau)

    surfdat >> dummy;                                 // \tau^2 . pi^\eta\eta [fm^-4]


    // bulk viscous pressure
    surfdat >> dummy;                                 // Pi [fm^-4]
    surf_ptr[i].bulkPi = dummy * hbarC;               // convert to [GeV/fm^3]


    double nB = 0.0;                                  // default value for net-baryon density

    if(include_baryon)
    {
      // net-baryon density
      surfdat >> nB;                                  // nB [fm^-3]
      surf_ptr[i].nB = nB;


      // contravariant net-baryon diffusion V^\mu
      surfdat >> dummy;                               // V^\tau [fm^-3]
      surfdat >> surf_ptr[i].Vx;                      // V^x [fm^-3]
      surfdat >> surf_ptr[i].Vy;                      // V^y [fm^-3]

      surfdat >> dummy;                               // \tau . V^\eta [fm^-3] (need to check music)
      surf_ptr[i].Vn = dummy / tau;                   // divide by \tau
    }


    // check whether 2+1d freezeout cells are really boost-invariant
    if(dimension == 2)
    {
      if(surf_ptr[i].eta != 0)
      {
      #ifdef FLAGS
        printf("read_surface_music flag: setting spacetime rapidity of boost-invariant freezeout cell to eta = 0\n");
      #endif

        surf_ptr[i].eta = 0;
      }
      if(surf_ptr[i].dan != 0 || surf_ptr[i].un != 0 || surf_ptr[i].pixn != 0 || surf_ptr[i].piyn != 0)
      {
      #ifdef FLAGS
        printf("read_surface_music flag: dimension = 2 but freezeout cell %ld is not boost-invariant (please check format in surface.dat)\n", i);
      #endif
      }
    }


    // compute average thermodynamic quantities
    double tau2 = tau * tau;
    double ux = surf_ptr[i].ux;
    double uy = surf_ptr[i].uy;
    double un = surf_ptr[i].un;
    double ut = sqrt(1.  +  ux * ux  +  uy * uy  +  tau2 * un * un);
    double dat = surf_ptr[i].dat;
    double dax = surf_ptr[i].dax;
    double day = surf_ptr[i].day;
    double dan = surf_ptr[i].dan;

    double uds = ut * dat  +  ux * dax  +  uy * day  +  un * dan;                 // u^\mu . d\sigma_\mu
    double ds_ds = dat * dat  -  dax * dax  -  day * day  -  dan * dan / tau2;    // d\sigma^\mu . d\sigma_\mu
    double ds_max = fabs(uds)  +  sqrt(fabs(uds * uds  -  ds_ds));                // max volume element |ds|

    max_volume += ds_max;         // append values
    E_avg += (E * ds_max);
    T_avg += (T * ds_max);
    P_avg += (P * ds_max);
    muB_avg += (muB * ds_max);
    nB_avg += (nB * ds_max);
  }

  surfdat.close();

  T_avg /= max_volume;            // divide by total max volume
  E_avg /= max_volume;
  P_avg /= max_volume;
  muB_avg /= max_volume;
  nB_avg /= max_volume;


  // write averaged thermodynamic variables to file
  ofstream thermal_average("average_thermodynamic_quantities.dat", ios_base::out);
  thermal_average << setprecision(15) << T_avg << "\n" << E_avg << "\n" << P_avg << "\n" << muB_avg << "\n" << nB_avg;
  thermal_average.close();
}


void FO_data_reader::read_surface_hic_eventgen(FO_surf* surf_ptr)
{
  printf("from input/surface.dat and undoing tau factors...\n");

  printf("\nHydrodynamic code = HIC-EventGen\n\n");
  printf("\thttps://github.com/Duke-QCD/hic-eventgen\n\n");

  if(dimension != 2)
  {
    printf("read_surface_hic_eventgen error: HIC-EventGen is boost-invariant (need to set dimension = 2)\n");
    exit(-1);
  }
  else if(include_baryon)
  {
    printf("read_surface_hic_eventgen error: HIC-EventGen does not consider baryon chemical potential (need to set include_baryon = 0)\n");
    exit(-1);
  }

  printf("Please check that input/surface.dat has the following format (and 1 blank line eof):\n\n\t");
  printf("[t x y n ds_t/t ds_x/t ds_y/t ds_n/t v^x v^y t.v^n pi^tt pi^tx pi^ty t.pi^tn pi^xx pi^xy t.pi^xn pi^yy t.pi^yn t2.pi^nn Pi T E P muB]\n\n");


  ostringstream surfdat_stream;                       // prepare to read in surface.dat
  surfdat_stream << "input/surface.dat";
  ifstream surfdat(surfdat_stream.str().c_str());
  double dummy;


  double T_avg = 0;                                   // average thermodynamic variables across freezeout surface
  double E_avg = 0;
  double P_avg = 0;
  double muB_avg = 0;
  double nB_avg = 0;
  double max_volume = 0;                              // max volume of freezeout surface


  for(long i = 0; i < number_of_cells; i++)           // loop over freezeout cells
  {
    // contravariant spacetime position x^\mu
    surfdat >> dummy;                                  // \tau [fm]
    double tau = dummy;
    surf_ptr[i].tau = tau;

    surfdat >> surf_ptr[i].x;                         // x [fm]
    surfdat >> surf_ptr[i].y;                         // y [fm]

    surfdat >> dummy;                                 // \eta_s [1]
    surf_ptr[i].eta = 0;


    // covariant surface normal vector d\sigma_\mu / \tau
    surfdat >> dummy;                                 // d\sigma_\tau / \tau [fm^-3]
    surf_ptr[i].dat = dummy * tau;                    // multiply by \tau

    surfdat >> dummy;                                 // d\sigma_x / \tau [fm^-3]
    surf_ptr[i].dax = dummy * tau;                    // multiply by \tau

    surfdat >> dummy;                                 // d\sigma_y / \tau [fm^-3]
    surf_ptr[i].day = dummy * tau;                    // multiply by \tau

    surfdat >> dummy;                                 // d\sigma_\eta / \tau [fm^-4]
    surf_ptr[i].dan = 0;



    // puzzled about this...
    // covariant fluid velocity                       // ask Derek if covariant...
    double vx;
    double vy;

    surfdat >> vx;                                    // u^x / u^\tau [1]
    surfdat >> vy;                                    // u^y / u^\tau [1]
    surfdat >> dummy;                                 // \tau . u^\eta / u^\tau

    double ut =  1. / sqrt(fabs(1.  -  vx * vx  -  vy * vy));

    surf_ptr[i].ux = ut * vx;                         // if covariant, would need minus sign...
    surf_ptr[i].uy = ut * vy;
    surf_ptr[i].un = 0;



    // contravariant shear stress pi^\mu\nu
    surfdat >> dummy;                                 // pi^\tau\tau [GeV/fm^3]
    surfdat >> dummy;                                 // pi^\taux [GeV/fm^3]
    surfdat >> dummy;                                 // pi^\taux [GeV/fm^3]
    surfdat >> dummy;                                 // \tau . pi^\tau\eta [GeV/fm^3]
    surfdat >> surf_ptr[i].pixx;                      // pi^xx [GeV/fm^3]
    surfdat >> surf_ptr[i].pixy;                      // pi^xy [GeV/fm^3]

    surfdat >> dummy;                                 // \tau . pi^x\eta [GeV/fm^3]
    surf_ptr[i].pixn = 0;

    surfdat >> surf_ptr[i].piyy;                      // pi^yy [GeV/fm^3]

    surfdat >> dummy;                                 // \tau . pi^y\eta [GeV/fm^3]
    surf_ptr[i].piyn = 0;

    surfdat >> dummy;                                 // \tau^2 . pi^\eta\eta [GeV/fm^3]


    // bulk viscous pressure
    surfdat >> surf_ptr[i].bulkPi;                    // Pi [GeV/fm^3]


    // thermodynamic variables
    surfdat >> dummy;                                 // temperature [GeV]
    double T = dummy;
    surf_ptr[i].T = T;

    surfdat >> dummy;                                 // energy density [GeV/fm^3]
    double E = dummy;
    surf_ptr[i].E = E;

    surfdat >> dummy;                                 // equilibrium pressure [GeV/fm^3]
    double P = dummy;
    surf_ptr[i].P = P;

    surfdat >> dummy;                                 // baryon chemical potential [GeV]
    double muB = dummy;
    surf_ptr[i].muB = muB;

    double nB = 0.0;                                  // default value for net-baryon density

    double tau2 = tau * tau;
    double ux = surf_ptr[i].ux;
    double uy = surf_ptr[i].uy;
    double un = 0;
    double dat = surf_ptr[i].dat;
    double dax = surf_ptr[i].dax;
    double day = surf_ptr[i].day;
    double dan = 0;

    double uds = ut * dat  +  ux * dax  +  uy * day  +  un * dan;                 // u^\mu . d\sigma_\mu
    double ds_ds = dat * dat  -  dax * dax  -  day * day  -  dan * dan / tau2;    // d\sigma^\mu . d\sigma_\mu
    double ds_max = fabs(uds)  +  sqrt(fabs(uds * uds  -  ds_ds));                // max volume element |ds|

    max_volume += ds_max;         // append values
    E_avg += (E * ds_max);
    T_avg += (T * ds_max);
    P_avg += (P * ds_max);
    muB_avg += (muB * ds_max);
    nB_avg += (nB * ds_max);
  }

  surfdat.close();

  T_avg /= max_volume;
  E_avg /= max_volume;
  P_avg /= max_volume;
  muB_avg /= max_volume;
  nB_avg /= max_volume;

  // write averaged thermodynamic variables to file
  ofstream thermal_average("average_thermodynamic_quantities.dat", ios_base::out);
  thermal_average << setprecision(15) << T_avg << "\n" << E_avg << "\n" << P_avg << "\n" << muB_avg << "\n" << nB_avg;
  thermal_average.close();

  return;
}




read_mcid::read_mcid(long int mcid_in)
{
  mcid = mcid_in;

  char sign = '\0';

  if(mcid < 0)
  {
    sign = '-';
    printf("Error: should only be particles (not antiparticles) in pdg_test.dat\n");
  }

  int * mcid_holder = (int *)calloc(max_digits, sizeof(int));

  long int x = abs(mcid_in);
  int digits = 1;

  // get individual digits from right to left
  for(int i = 0; i < max_digits; i++)
  {
    mcid_holder[i] = (int)(x % (long int)10);

    x /= (long int)10;

    if(x > 0)
    {
      digits++;
      if(digits > max_digits) printf("Error: mcid %ld is > %d digits\n", mcid_in, max_digits);
    }
  }

  // set the quantum numbers
  nJ = mcid_holder[0];      // hadrons have 7-digit codes
  nq3 = mcid_holder[1];
  nq2 = mcid_holder[2];
  nq1 = mcid_holder[3];
  nL = mcid_holder[4];
  nR = mcid_holder[5];
  n = mcid_holder[6];
  n8 = mcid_holder[7];      // there are some hadrons with 8-digit codes
  n9 = mcid_holder[8];
  n10 = mcid_holder[9];     // nuclei have 10-digit codes

  nJ += n8;                 // I think n8 adds to nJ if spin > 9

  // test print the mcid
  //cout << sign << n << nR << nL << nq1 << nq2 << nq3 << nJ << endl;

  free(mcid_holder);  // free memory

  // get relevant particle properties
  is_particle_a_deuteron();
  is_particle_a_hadron();
  is_particle_a_meson();
  is_particle_a_baryon();
  get_spin();
  get_gspin();
  get_baryon();
  get_sign();
  does_particle_have_distinct_antiparticle();
}


void read_mcid::is_particle_a_deuteron()
{
  // check if particle is a deuteron (easy way)
  is_deuteron = (mcid == 1000010020);

  if(is_deuteron) printf("Error: there is a deuteron in HRG\n");
}

void read_mcid::is_particle_a_hadron()
{
  // check if particle is a hadron
  is_hadron = (!is_deuteron && nq3 != 0 && nq2 != 0);
}

void read_mcid::is_particle_a_meson()
{
  // check if particle is a meson
  is_meson = (is_hadron && nq1 == 0);
}

void read_mcid::is_particle_a_baryon()
{
  // check if particle is a baryon
  is_baryon = (is_hadron && nq1 != 0);
}

void read_mcid::get_spin()
{
  // get the spin x 2 of the particle
  if(is_hadron)
  {
    if(nJ == 0)
    {
      spin = 0;  // special cases: K0_L=0x130 & K0_S=0x310
      return;
    }
    else
    {
      spin = nJ - 1;
      return;
    }
  }
  else if(is_deuteron)
  {
    spin = 2;
    return;
  }
  else
  {
    printf("Error: particle is not a deuteron or hadron\n");
    // this assumes that we only have white particles (no single
    // quarks): Electroweak fermions have 11-17, so the
    // second-to-last-digit is the spin. The same for the Bosons: they
    // have 21-29 and 2spin = 2 (this fails for the Higgs).
    spin = nq3;
    return;
  }

}


void read_mcid::get_gspin()
{
  // get the spin degeneracy of the particle
  if(is_hadron && nJ > 0)
  {
    gspin = nJ;
    return;
  }
  else if(is_deuteron)
  {
    gspin = 3;  // isospin is 0 -> spin = 1 (generalize for any I later)
    return;
  }
  else
  {
    printf("Error: particle is not a deuteron or hadron\n");
    gspin = spin + 1; // lepton
    return;
  }

}

void read_mcid::get_baryon()
{
  // get the baryon number of the particle
  if(is_deuteron)
  {
    //baryon = nq3  +  10 * nq2  +  100 * nq1; // nucleus
    baryon = 2;
    return;
  }
  else if(is_hadron)
  {
    if(is_meson)
    {
      baryon = 0;
      return;
    }
    else if(is_baryon)
    {
      baryon = 1;
      return;
    }
  }
  else
  {
    printf("Error: particle is not a deuteron or hadron\n");
    baryon = 0;
    return;
  }
}

void read_mcid::get_sign()
{
  // get the quantum statistics sign of the particle
  if(is_deuteron)
  {
    sign = -1;   // deuteron is boson
    return;
  }
  else if(is_hadron)
  {
    if(is_meson)
    {
      sign = -1;
      return;
    }
    else if(is_baryon)
    {
      sign = 1;
      return;
    }
  }
  else
  {
    printf("Error: particle is not a deuteron or hadron\n");
    sign = spin % 2;
    return;
  }
}

void read_mcid::does_particle_have_distinct_antiparticle()
{
  // check if the particle has distinct antiparticle
  if(is_hadron) // hadron
  {
    has_antiparticle = ((baryon != 0) || (nq2 != nq3));
    return;
  }
  else if(is_deuteron)
  {
    has_antiparticle = true;
  }
  else
  {
    printf("Error: particle is not a deuteron or hadron\n");
    has_antiparticle = (nq3 == 1); // lepton
    return;
  }
}


PDG_Data::PDG_Data(ParameterReader * paraRdr_in)
{
  paraRdr = paraRdr_in;
  hrg_eos = paraRdr->getVal("hrg_eos");
}


PDG_Data::~PDG_Data()
{
  //////////////////////////////
}


int PDG_Data::read_resonances_conventional(particle_info * particle, string pdg_filename)
{
  double eps = 1.e-15;
  int Nparticle = 0;

  ifstream resofile(pdg_filename);
  int local_i = 0;
  int dummy_int;

  while(!resofile.eof())
  {
    resofile >> particle[local_i].mc_id;
    resofile >> particle[local_i].name;
    resofile >> particle[local_i].mass;
    resofile >> particle[local_i].width;
    resofile >> particle[local_i].gspin;	      //spin degeneracy
    resofile >> particle[local_i].baryon;
    resofile >> particle[local_i].strange;
    resofile >> particle[local_i].charm;
    resofile >> particle[local_i].bottom;
    resofile >> particle[local_i].gisospin;     //isospin degeneracy
    resofile >> particle[local_i].charge;
    resofile >> particle[local_i].decays;

    for(int j = 0; j < particle[local_i].decays; j++)
    {
      resofile >> dummy_int;
      resofile >> particle[local_i].decays_Npart[j];
      resofile >> particle[local_i].decays_branchratio[j];
      resofile >> particle[local_i].decays_part[j][0];
      resofile >> particle[local_i].decays_part[j][1];
      resofile >> particle[local_i].decays_part[j][2];
      resofile >> particle[local_i].decays_part[j][3];
      resofile >> particle[local_i].decays_part[j][4];
    }

    //decide whether particle is stable under strong interactions
    if (particle[local_i].decays_Npart[0] == 1) particle[local_i].stable = 1;
    else particle[local_i].stable = 0;

    //add anti-particle entry
    if (particle[local_i].baryon > 0) // changed on Feb. 2019
    {
      local_i++;
      particle[local_i].mc_id = -particle[local_i-1].mc_id;

      ostringstream antiname;
      antiname << "Anti-baryon-" << particle[local_i-1].name;

      particle[local_i].name = antiname.str();
      particle[local_i].mass = particle[local_i-1].mass;
      particle[local_i].width = particle[local_i-1].width;
      particle[local_i].gspin = particle[local_i-1].gspin;
      particle[local_i].baryon = -particle[local_i-1].baryon;
      particle[local_i].strange = -particle[local_i-1].strange;
      particle[local_i].charm = -particle[local_i-1].charm;
      particle[local_i].bottom = -particle[local_i-1].bottom;
      particle[local_i].gisospin = particle[local_i-1].gisospin;
      particle[local_i].charge = -particle[local_i-1].charge;
      particle[local_i].decays = particle[local_i-1].decays;
      particle[local_i].stable = particle[local_i-1].stable;

      for (int j = 0; j < particle[local_i].decays; j++)
      {
        particle[local_i].decays_Npart[j]=particle[local_i-1].decays_Npart[j];
        particle[local_i].decays_branchratio[j]=particle[local_i-1].decays_branchratio[j];

        for (int k=0; k< Maxdecaypart; k++)
        {
          if(particle[local_i-1].decays_part[j][k] == 0) particle[local_i].decays_part[j][k] = (particle[local_i-1].decays_part[j][k]);
          else
          {
            int idx;
            // find the index for decay particle
            for(idx = 0; idx < local_i; idx++)
            if (particle[idx].mc_id == particle[local_i-1].decays_part[j][k]) break;
            if(idx == local_i && particle[local_i-1].stable == 0 && particle[local_i-1].decays_branchratio[j] > eps)
            {
              cout << "Error: can not find decay particle index for anti-baryon!" << endl;
              cout << "particle mc_id : " << particle[local_i-1].decays_part[j][k] << endl;
              exit(1);
            }
            if (particle[idx].baryon == 0 && particle[idx].charge == 0 && particle[idx].strange == 0) particle[local_i].decays_part[j][k] = (particle[local_i-1].decays_part[j][k]);
            else particle[local_i].decays_part[j][k] = (- particle[local_i-1].decays_part[j][k]);
          }
        }
      }
    }
    local_i++;	// Add one to the counting variable "i" for the meson/baryon
  }
  resofile.close();
  Nparticle = local_i - 1; //take account the final fake one
  for (int i = 0; i < Nparticle; i++)
  {
    //if (particle[i].baryon == 0) particle[i].sign = -1; // this made deuteron a fermion
    if (particle[i].baryon % 2 == 0) particle[i].sign = -1;
    else particle[i].sign = 1;
  }

  // count the number of mesons, baryons and antibaryons
  int meson = 0;
  int baryon = 0;
  int antibaryon = 0;

  for(int i = 0; i < Nparticle; i++)
  {
    if(particle[i].baryon == 0) meson++;
    else if(particle[i].baryon > 0) baryon++;
    else if(particle[i].baryon < 0) antibaryon++;
  }
  if(baryon != antibaryon) printf("Error: (anti)baryons not paired correctly\n");

  printf("\nNumber of resonances = %d\n\n\t", Nparticle);
  printf("%d mesons\n\t", meson);
  printf("%d baryons\n\t", baryon);
  printf("%d antibaryons\n\n", antibaryon);

  particle_info last_particle = particle[Nparticle - 1];

  printf("Last particle has mcid = %ld, %s, m = %lf GeV (please check this) \n\n", last_particle.mc_id, last_particle.name.c_str(), last_particle.mass);

  return Nparticle;
}


int PDG_Data::read_resonances_smash_box(particle_info * particle, string pdg_filename)
{
  //******************************|
  //******************************|
  const int mcid_entries = 4;   //|   (increase if > 4 mcid entries per line)
  //******************************|
  //******************************|

  printf("\nNumber of mcid entries per row is set to mcid_entries = %d (increase if needed)\n", mcid_entries);

  long int * mc_id = (long int*)calloc(mcid_entries, sizeof(long int));

  string name;
  double mass;
  double width;
  char parity;

  string line;

  ifstream pdg_smash_box(pdg_filename);

  int i = 0;  // particle index

  while(getline(pdg_smash_box, line))
  {
      istringstream current_line(line);

      // add mc_id[..] if increase mcid_entries
      current_line >> name >> mass >> width >> parity >> mc_id[0] >> mc_id[1] >> mc_id[2] >> mc_id[3];

      if(line.empty() || line.at(0) == '#')
      {
        // skip lines that are blank or start with comments
        continue;
      }

      for(int k = 0; k < mcid_entries; k++)
      {
        // add particles with nonzero mc_id's to particle struct
        if(mc_id[k] != 0)
        {
          // get remaining particle info from the mcid
          read_mcid mcid_info(mc_id[k]);

          particle[i].name = name;
          particle[i].mass = mass;
          particle[i].width = width;
          particle[i].mc_id = mc_id[k];
          particle[i].gspin = mcid_info.gspin;
          particle[i].baryon = mcid_info.baryon;
          particle[i].sign = mcid_info.sign;

          i++;

          if(mcid_info.has_antiparticle)
          {
            // create antiparticle next to hadron
            ostringstream antiname;
              antiname << "Anti-" << name;

            particle[i].name = antiname.str();
            particle[i].mass = mass;
            particle[i].width = width;
            particle[i].mc_id = -mc_id[k];
            particle[i].gspin = mcid_info.gspin;
            particle[i].baryon = -mcid_info.baryon;
            particle[i].sign = mcid_info.sign;

            i++;

          } // add antiparticle

        } // check if mc_id if nonzero

        // reset mc_id's to zero
        mc_id[k] = 0;

      } // mcid index

      if(i > (Maxparticle -1))
      {
        printf("\nError: number of particles in file exceeds Maxparticle = %d. Exiting...\n\n", Maxparticle);
        exit(-1);
      }

  } // scanning file

  pdg_smash_box.close();

  free(mc_id);

  int Nparticle = i;  // total number of resonances

  // count the number of mesons, baryons and antibaryons
  int meson = 0;
  int baryon = 0;
  int antibaryon = 0;

  for(int i = 0; i < Nparticle; i++)
  {
    if(particle[i].baryon == 0) meson++;
    else if(particle[i].baryon > 0) baryon++;
    else if(particle[i].baryon < 0) antibaryon++;
  }
  if(baryon != antibaryon) printf("Error: (anti)baryons not paired correctly\n");

  printf("\nNumber of resonances = %d\n\n\t", Nparticle);
  printf("%d mesons\n\t", meson);
  printf("%d baryons\n\t", baryon);
  printf("%d antibaryons\n\n", antibaryon);

  particle_info last_particle = particle[Nparticle - 1];

  printf("Last particle is mcid = %ld, %s, m = %lf GeV (please check this) \n\n", last_particle.mc_id, last_particle.name.c_str(), last_particle.mass);

  return Nparticle;
}


int PDG_Data::read_resonances(particle_info * particle)
{
  int Nparticle;

  switch(hrg_eos)
  {
    case 1:
    {
      printf("PDG/pdg-urqmd_v3.3+.dat... (please check if 1 blank line eof)\n");

      Nparticle = read_resonances_conventional(particle, urqmd);      // read urqmd
      break;
    }
    case 2:
    {
      printf("PDG/pdg_smash.dat... (please check if 1 blank line eof)\n");

      Nparticle = read_resonances_conventional(particle, smash);      // read smash
      break;
    }
    case 3:
    {
      printf("PDG/pdg_box.dat...\n");

      Nparticle = read_resonances_smash_box(particle, smash_box);     // read smash box
      break;
    }
    default:
    {
      printf("\nread_resonances error: need to set hrg_eos = (1,2,3)\n");
      exit(-1);
    }
  }

  return Nparticle;
}

