#ifndef PARTICLE_H
#define PARTICLE_H

// class particle
// {
// public:
//   long int 	mcID;     //Montecarlo number according PDG
//   char	name[26];
//   double	mass;
//   double	width;
//   int	    gspin;      //spin degeneracy
//   int	    baryon;
//   int	    strange;
//   int	    charm;
//   int	    bottom;
//   int	    gisospin;  // isospin degeneracy
//   int	    charge;
//   int	    decays;    // amount of decays listed for this resonance
//   int	    stable;     // defines whether this particle is considered as stable
// };

// class particleDecay
// {
// public:
//   int	reso;       // Montecarlo number of decaying resonance
//   int	numpart;    // number of daughter particles after decay
//   double branch;  // branching ratio
//   int	part[5];    // array of daughter particles Montecarlo numbers
// };


class Sampled_Particle
{
  // info of a sampled particle in the particle list, which gets passed to an afterburner

public:

  int chosen_index = 0;     // chosen particle index
  int mcID = 0;             // Monte-Carlo ID number
  double mass = 0;          // mass [GeV]

  double tau = 0;           // spacetime position of sampled particle
  double x = 0;
  double y = 0;
  double eta = 0;

  double t = 0;             // cartesian coordinates
  double z = 0;

  double E = 0;             // cartesian momentum
  double px = 0;
  double py = 0;
  double pz = 0;
};

#endif