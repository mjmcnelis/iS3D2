
#include <cmath>

#include "iS3D.h"
#include "EmissionFunction.h"

using namespace std;

void EmissionFunctionArray::sample_dN_dy(int chosen_index, double y)
{
  // bin sampled dN/dy

  // bin index
  int iy = (int)floor((y + Y_CUT) / Y_WIDTH);

  // add counts
  if(iy >= 0 && iy < Y_BINS)
  {
    dN_dy_count[chosen_index][iy] += 1.0;
  }
}

void EmissionFunctionArray::sample_dN_deta(int chosen_index, double eta)
{
  // bin sampled dN/deta

  // bin index
  int ieta = (int)floor((eta + ETA_CUT) / ETA_WIDTH);

  if(ieta >= 0 && ieta < ETA_BINS)
  {
    dN_deta_count[chosen_index][ieta] += 1.0;
  }
}

void EmissionFunctionArray::sample_dN_dphipdy(int chosen_index, double px, double py)
{
  // bin sampled dN/dphipdy

  double phip = atan2(py, px);

  if(phip < 0.0)
  {
    phip += two_pi;
  }

  // bin index
  int iphip = (int)floor(phip / PHIP_WIDTH);

  // add counts
  if(iphip >= 0 && iphip < PHIP_BINS)
  {
    dN_dphipdy_count[chosen_index][iphip] += 1.0;
  }
}

void EmissionFunctionArray::sample_dN_2pipTdpTdy(int chosen_index, double px, double py)
{
  // bin sampled dN/2pipTdpTdy

  double pT = sqrt(px*px + py*py);

  // bin index
  int ipT = (int)floor((pT - PT_MIN) / PT_WIDTH);

  if(ipT >= 0 && ipT < PT_BINS)
  {
    dN_2pipTdpTdy_count[chosen_index][ipT] += 1.0;
  }
}

void EmissionFunctionArray::sample_vn(int chosen_index, double px, double py)
{
  // bin vn(pT)

  double pT = sqrt(px*px + py*py);
  double phi = atan2(py, px);

  if(phi < 0.0)
  {
    phi += two_pi;
  }

  // pT bin index
  int ipT = (int)floor((pT - PT_MIN) / PT_WIDTH);

  if(ipT >= 0 && ipT < PT_BINS)
  {
    // pT count
    pT_count[chosen_index][ipT] += 1.0;

    for(int k = 0; k < K_MAX; k++)
    {
      // Vn count
      vn_real_count[k][chosen_index][ipT] += cos(((double)k + 1.0) * phi);
      vn_imag_count[k][chosen_index][ipT] += sin(((double)k + 1.0) * phi);
    }
  }
}

void EmissionFunctionArray::sample_dN_dX(int chosen_index, double tau, double x, double y)
{
  // construct spacetime distributions
  // assume dN/deta = dN/dy (boost-invariance)

  double r = sqrt(x*x + y*y);
  double phi = atan2(y, x);

  if(phi < 0.0)
  {
    phi += two_pi;
  }

  // bin indices
  int itau = (int)floor((tau - TAU_MIN) / TAU_WIDTH);
  int ir = (int)floor((r - R_MIN) / R_WIDTH);
  int iphi = (int)floor(phi / PHIP_WIDTH);

  if(itau >= 0 && itau < TAU_BINS)
  {
    dN_taudtaudy_count[chosen_index][itau] += 1.0;
  }

  if(ir >= 0 && ir < R_BINS)
  {
    dN_twopirdrdy_count[chosen_index][ir] += 1.0;
  }

  if(iphi >= 0 && iphi < PHIP_BINS)
  {
    dN_dphisdy_count[chosen_index][iphi] += 1.0;
  }
}
