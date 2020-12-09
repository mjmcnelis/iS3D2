
#ifndef MOMENTUM_H
#define MOMENTUM_H

#include "LocalRestFrame.h"


typedef struct
{
  double E;     // u.p
  double px;    // -X.p
  double py;    // -Y.p
  double pz;    // -Z.p
  double feq; 	// thermal distribution
} LRF_Momentum;	// local rest frame momentum


class Lab_Momentum
{
  private:          // LRF momentum components
    double E_LRF;
    double px_LRF;
    double py_LRF;
    double pz_LRF;

  public:           // contravariant lab frame momentum p^mu (milne components)
    double ptau;    // p^tau
    double px;      // p^x
    double py;      // p^y
    double pn;      // p^eta

    Lab_Momentum(LRF_Momentum pLRF_in);
    void boost_pLRF_to_lab_frame(Milne_Basis basis_vectors, double ut, double ux, double uy, double un);
};

#endif