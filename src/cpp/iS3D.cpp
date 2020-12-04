
#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include <string.h>
#include<sstream>
#include<cmath>
#include<vector>
#include<sys/time.h>
#include "iS3D.h"
#include "Macros.h"
#include "Table.h"
#include "readindata.h"
#include "emissionfunction.h"
#include "arsenal.h"
#include "ParameterReader.h"
#include "deltafReader.h"


IS3D::IS3D()
{

}


IS3D::~IS3D()
{

}


void IS3D::read_fo_surf_from_memory(
                                    std::vector<double> tau_in,
                                    std::vector<double> x_in,
                                    std::vector<double> y_in,
                                    std::vector<double> eta_in,
                                    std::vector<double> dsigma_tau_in,
                                    std::vector<double> dsigma_x_in,
                                    std::vector<double> dsigma_y_in,
                                    std::vector<double> dsigma_eta_in,
                                    std::vector<double> E_in,
                                    std::vector<double> T_in,
                                    std::vector<double> P_in,
                                    std::vector<double> ux_in,
                                    std::vector<double> uy_in,
                                    std::vector<double> un_in,
                                    std::vector<double> pixx_in,
                                    std::vector<double> pixy_in,
                                    std::vector<double> pixn_in,
                                    std::vector<double> piyy_in,
                                    std::vector<double> piyn_in,
                                    std::vector<double> pinn_in,
                                    std::vector<double> Pi_in
                                   )
{
  tau = tau_in;
  x = x_in;
  y = y_in;
  eta = eta_in;
  dsigma_tau = dsigma_tau_in;
  dsigma_x = dsigma_x_in;
  dsigma_y = dsigma_y_in;
  dsigma_eta = dsigma_eta_in;
  E = E_in;
  T = T_in;
  P = P_in;
  ux = ux_in;
  uy = uy_in;
  un = un_in;
  pixx = pixx_in;
  pixy = pixy_in;
  pixn = pixn_in;
  piyy = piyy_in;
  piyn = piyn_in;
  pinn = pinn_in;
  Pi = Pi_in;
}

void IS3D::run_particlization(int fo_from_file)
{
  printf("\n::::::::::::::::::::::::::::::::::::::::\n");
  printf("::                                    ::\n");
  printf("::    Starting iS3D particlization    ::\n");
  printf("::                                    ::\n");
  printf("::::::::::::::::::::::::::::::::::::::::\n\n");


#ifdef PRINT_PROCESSOR
  print_processor();                                    // print processor in use (for benchmarking)
#endif


  printf("\n\nReading in parameters...\n\n");
  ParameterReader *paraRdr = new ParameterReader;       // parameter reader class
  paraRdr->readFromFile("iS3D_parameters.dat");
  int include_baryon = paraRdr->getVal("include_baryon");

#ifdef PRINT_PARAMETERS
  paraRdr->echo();
#endif


  printf("\n\nReading in freezeout surface ");
  long FO_length = 0;
  FO_data_reader freeze_out_data(paraRdr, "input");     // freezeout reader class

  if(fo_from_file == 1)
  {
    FO_length = freeze_out_data.get_number_cells();     // get length of freezeout surface from file
  }
  else
  {
    FO_length = tau.size();                             // get length of freezeout surface from memory
  }


  FO_surf* surf_ptr = new FO_surf[FO_length];           // freezeout surface pointer

  if(fo_from_file == 1)
  {
    printf("file ");
    freeze_out_data.read_freezeout_surface(surf_ptr);   // load freezeout surface info from file
  }
  else
  {
    printf("from memory (please check that you've already undone hbarc = 1 units in hydro module)...\n\n");

    for(long is = 0; is < FO_length; is++)              // load freezeout surface info from memory
    {
      // contravariant spacetime position x^\mu
      surf_ptr[is].tau = tau[is];                       // \tau [fm]
      surf_ptr[is].x = x[is];                           // x [fm]
      surf_ptr[is].y = y[is];                           // y [fm]
      surf_ptr[is].eta = eta[is];                       // \eta_s [1]

      // covariant surface normal vector
      surf_ptr[is].dat = dsigma_tau[is];                // d\sigma_\tau [fm^-2]
      surf_ptr[is].dax = dsigma_x[is];                  // d\sigma_x [fm^-2]
      surf_ptr[is].day = dsigma_y[is];                  // d\sigma_y [fm^-2]
      surf_ptr[is].dan = dsigma_eta[is];                // d\sigma_\eta [fm^-3]

      // contravariant fluid velocity u^\mu
      surf_ptr[is].ux = ux[is];                         // u^x [1]
      surf_ptr[is].uy = uy[is];                         // u^y [1]
      surf_ptr[is].un = un[is];                         // u^\eta [fm^-1]

      // thermodynamic variables
      surf_ptr[is].E = E[is];                           // energy density [GeV/fm^3]
      surf_ptr[is].T = T[is];                           // temperature [GeV]
      surf_ptr[is].P = P[is];                           // equilibrium pressure [GeV/fm^3]

      // contravariant shear stress pi^\mu\nu
      surf_ptr[is].pixx = pixx[is];                     // pi^xx [GeV/fm^3]
      surf_ptr[is].pixy = pixy[is];                     // pi^xy [GeV/fm^3]
      surf_ptr[is].pixn = pixn[is];                     // pi^x\eta [GeV/fm^4]
      surf_ptr[is].piyy = piyy[is];                     // pi^yy [GeV/fm^3]
      surf_ptr[is].piyn = piyn[is];                     // pi^y\eta [GeV/fm^4]

      // bulk viscous pressure
      surf_ptr[is].bulkPi = Pi[is];                     // Pi [GeV/fm^3]

      // JETSCAPE does not consider (muB, nB, V^\mu) at the moment

    }
  }


  printf("\n\nReading in particle info from ");
  particle_info *particle_data = new particle_info[Maxparticle];  // particle info pointer
  PDG_Data pdg(paraRdr);                                          // PDG class
  int Nparticle = pdg.read_resonances(particle_data);             // get number of resonances in PDG file


  printf("Reading in chosen particles table... (please check that PDG/chosen_particles.dat has 1 blank line eof)\n\n");
  Table chosen_particles("PDG/chosen_particles.dat");             // chosen particles table


  // df coefficient data
  Deltaf_Data * df_data = new Deltaf_Data(paraRdr);               // df data pointer
  df_data->load_df_coefficient_data();                            // read in df coefficient tables

  if(!include_baryon)
  {
    df_data->construct_cubic_splines();                           // prepare cubic spline interpolation (muB = 0)
    df_data->compute_jonah_coefficients(particle_data, Nparticle);// compute PTB exclusive coefficients (muB = 0)
  }

  df_data->compute_particle_densities(particle_data, Nparticle);  // compute resonances' particle density (T = T_avg, muB = muB_avg)
  df_data->test_df_coefficients(-0.1);                            // test df coefficients for bulk pressure Pi = -Peq/10


  printf("Number of freezeout cells = %ld\n", FO_length);
  printf("Number of chosen particles = %d\n", chosen_particles.getNumberOfRows());


  // read in momentum and spacetime rapidity tables
  Table pT_tab("tables/momentum/pT_table.dat");             // pT table
  Table phi_tab("tables/momentum/phi_table.dat");           // phi table
  Table y_tab("tables/momentum/y_table.dat");               // y table (for 3+1d smooth CFF)
  Table eta_tab("tables/spacetime_rapidity/eta_table.dat"); // eta table (for 2+1d smooth CFF)

  EmissionFunctionArray efa(paraRdr, &chosen_particles, &pT_tab, &phi_tab, &y_tab, &eta_tab, particle_data, Nparticle, surf_ptr, FO_length, df_data);

  std::vector<Sampled_Particle> particle_event_list_in;
  efa.calculate_spectra(particle_event_list_in);

  //copy final particle list to memory to pass to JETSCAPE module
  int operation = paraRdr->getVal("operation");
  if (operation == 2)
  {
    cout << "Copying final particle list to memory" << endl;
    cout << "Particle list contains " << particle_event_list_in.size() << " particles" << endl;
    final_particles_ = particle_event_list_in;
  }

  delete [] surf_ptr;
  delete paraRdr;
  delete df_data;

  printf("\nFinished calculating particle spectra\n");

}
