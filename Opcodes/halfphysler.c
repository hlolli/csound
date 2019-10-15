/*
  halfphysler.c:

  This is a C port of the C++ csound-plugin halfphysler
  Copyright (C) 2018 - Sebastian Schmutzhard, Alex Hofmann,
  Gokberk Erdogan and Vasileios Chatziioannou
  GitHub: https://github.com/ketchupok/half-physler
  Ported to C from C++ in October 2019 - Hlöðver Sigurðsson

  This file is part of Csound.

  The Csound Library is free software; you can redistribute it
  and/or modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  Csound is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with Csound; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02110-1301 USA
*/

#include <stdlib.h>
#include "interlocks.h"
#include "csoundCore.h"

#ifndef MIN
#define MIN(a,b) ((a>b)?(b):(a))
#endif

typedef struct {
  OPDS h;
  MYFLT *out_feedback,
    *out_sound,
    *in,
    *L,                         // Length of the tube
    *r,                         // Tube radius
    *slope,                     // Slope of the tube
    *mult_alpha,                // user_multiplier for radiation losses
    *mult_rho,                  // user_multiplier for density
    *pickup_pos;                // user defined pickup position in tube
  AUXCH pold;                   // Pressure value at (n)th time grid
  AUXCH vold;                   // Velocity value at (n)th time grid
  AUXCH pnew;                   // Pressure value at (n+1)th time grid
  AUXCH vnew;                   // velocity value at (n+1)th time grid
  AUXCH S;                      // Cross sectional area
  uint32_t M;                   // Number of steps
  uint32_t Mold;                // Previous number of steps
  MYFLT Lold;                   // Previous length of the tube
  MYFLT dt;                     // Time grid
  MYFLT fs;                     // Sampling rate
  MYFLT dx;                     // Step size
  MYFLT dxold;                  // Previous step size
  MYFLT rad_alphaS;             // Radiation parameter alpha
  MYFLT rad_betaS;              // Radiation parameter beta
  MYFLT r_old;                  // Prev. tube radius
  MYFLT slope_old;              // Prev. slope of the tube
} CONE_RADIATION_LOSSES;

MYFLT c = 3.4386e+02;      // speed of sound
MYFLT rho = 1.2000e+00;    // density
int32_t Mmax  = 400;       // Maximal number of grid points

// constants taken from Bilbao Numerical Sound Synthesis (2009) formula (9.11)
MYFLT rad_alpha  =  2.890027475795721e+00;
MYFLT rad_beta   =  6.646516378651400e-01;


void grid_init(MYFLT Len, MYFLT dt, MYFLT * dx, uint32_t * M, MYFLT * L) {
  // TODO(AH) - replace weird arrays with [subscript]
  // INPUTS: (total) length
  // calculate time step and grid size, taking stability into account
  // return length and grid size (M and dx)

  // In order to satisfy the stability conditions dx/dt>c and M*dx = L
  L[0] = Len;
  dx[0]  = c*dt;
  dx[0]  = dx[0]/0.9;  // Guarantees the first stability condition

  M[0]   = (uint32_t)floor(L[0] / dx[0]);  // Number of steps has to be an integer
  M[0]   = MIN(M[0], Mmax);
  dx[0]  = L[0] / (MYFLT) M[0];  // Guarantees both conditions
}

MYFLT cross_area(MYFLT radius, MYFLT x, MYFLT prelength, MYFLT slope) {
  /*x is the current position(m*dx), prelength is the length of the total tube
    up until the current segment (used for concatenation), radius is the
    input radius and slope is a ratio (dr/dx)*/
  MYFLT area;
  area = (radius+(x-prelength)*slope)*(radius+(x-prelength)*slope)*PI_F;
  return area;
}

void interpolation(uint32_t M, uint32_t Mold, MYFLT Lold, MYFLT dx,
                   MYFLT dxold, MYFLT *knew, MYFLT *kold) {

  MYFLT x, xl, xr;
  uint32_t m1, m2;
  for (int m = 0; m <= MIN(M, FLOOR(Lold / dx)); m++) {
    m1 = FLOOR(m*(dx / dxold));
    m2 = m1 + 1;  // Equivalent to m2 = ceil(m*dx/dxold)

    xl = m1 * dxold;
    xr = m2 * dxold;
    x  = m * dx;
    kold[m] = knew[m1]*((xr - x)/(xr - xl)) + knew[m2]*((x-xl)/(xr-xl));
  }

  // when L increases k is set to last old value (e.g. pressure at tube end)
  for (int m = MIN(M, (uint32_t) FLOOR(Lold/dx)) + 1; m <= M; m++ ) {
    kold[m] = knew[Mold];
  }
}

void update_vp(int M, MYFLT dt, MYFLT dx, \
               const MYFLT c, const MYFLT rho_user,    \
               const MYFLT *S,                         \
               const MYFLT *pold, const MYFLT *vold,   \
               MYFLT  *pnew, MYFLT *vnew) {
  //  calculate flow up to last grid point
  for (int m = 1; m <= M; m++) {  // Solve the first diff. eq. wrt vnew
    vnew[m] = -dt/rho_user * (pold[m] - pold[m-1]) / dx + vold[m];
  }

  //  calculate pressure up to one grid point before the last
  // (radiation in resonator cpp file)
  for (int m = 0; m <= M-1; m++) {  // Solve the second diff. eq. wrt pnew
    pnew[m] = - dt *rho_user *c*c / S[m] * (S[m+1]*vnew[m+1] - S[m]*vnew[m])/dx + pold[m];
  }
}

static int32_t cone_radiation_losses_init(CSOUND* csound, CONE_RADIATION_LOSSES* p) {

  // Allocate array elements
  uint32_t aux_size = sizeof(MYFLT) * (Mmax + 1);
  csound->AuxAlloc(csound, aux_size, &(p->pold));
  csound->AuxAlloc(csound, aux_size, &(p->vold));
  csound->AuxAlloc(csound, aux_size, &(p->pnew));
  csound->AuxAlloc(csound, aux_size, &(p->vnew));
  csound->AuxAlloc(csound, aux_size, &(p->S));

  MYFLT *pold = (MYFLT*)p->pold.auxp;
  MYFLT *vold = (MYFLT*)p->vold.auxp;
  MYFLT *pnew = (MYFLT*)p->pnew.auxp;
  MYFLT *vnew = (MYFLT*)p->vnew.auxp;
  MYFLT *S    = (MYFLT*)p->S.auxp;

  p->dt = 1.0/CS_ESR;

  grid_init(* p->L, p->dt, &p->dx, &p->M, p->L);

  // Set all array elements to 0
  for (int m = 0; m <= p->M; m++) {

    pold[m] = 0.0;
    vold[m] = 0.0;
    pnew[m] = 0.0;
    vnew[m] = 0.0;

    S[m] = cross_area(*p->r, m * p->dx, 0, *p->slope);  // Input radius of cone
  }

  p->rad_alphaS = (rad_alpha * *p->mult_alpha) / sqrt(S[p->M]);  // normalize radiation parameters
  p->rad_betaS = rad_beta / c;
  p->Lold      = *p->L;
  p->r_old     = *p->r;
  p->slope_old = *p->slope;
  p->Mold      = p->M;
  p->dxold     = p->dx;

  return OK;

}


static int32_t cone_radiation_losses_aperf(CSOUND* csound, CONE_RADIATION_LOSSES *p) {

  MYFLT *out_sound = p->out_sound;
  MYFLT *out_feedback = p->out_feedback;

  MYFLT *pold = (MYFLT*)p->pold.auxp;
  MYFLT *vold = (MYFLT*)p->vold.auxp;
  MYFLT *pnew = (MYFLT*)p->pnew.auxp;
  MYFLT *vnew = (MYFLT*)p->vnew.auxp;
  MYFLT *S    = (MYFLT*)p->S.auxp;

  uint32_t aux_size = sizeof(MYFLT) * (Mmax + 1);
  uint32_t offset = p->h.insdshead->ksmps_offset;
  uint32_t early  = p->h.insdshead->ksmps_no_end;
  uint32_t i, nsmps = CS_KSMPS;


  p->rad_alphaS = (rad_alpha * *p->mult_alpha) / sqrt(S[p->M]);  // normalization

  if (UNLIKELY(offset)) {
    memset(p->out_sound, '\0', offset*sizeof(MYFLT));
    memset(p->out_feedback, '\0', offset*sizeof(MYFLT));
  }
  if (UNLIKELY(early)) {
    nsmps -= early;
    memset(&p->out_sound[nsmps], '\0', early*sizeof(MYFLT));
    memset(&p->out_feedback[nsmps], '\0', early*sizeof(MYFLT));
  }

  // re-calculate the grid
  if (*p->L != p->Lold || *p->r != p->r_old || *p->slope != p->slope_old) {
    // Ensures that new calculations are made only when L,r or slope changed
    grid_init(* p->L, p->dt, &p->dx, &p->M, p->L); // make grid with new geometry

    // new cross-sectional area
    for (int m = 0; m <= p->M; m++) {
      S[m] = cross_area(*p->r, m * p->dx, 0, *p->slope);
    }

    // interpolate old grid status to new grid for each point
    interpolation(p->M, p->Mold, p->Lold, p->dx, p->dxold, &pnew[0], &pold[0]);
    interpolation(p->M, p->Mold, p->Lold, p->dx, p->dxold, &vnew[0], &vold[0]);
    p->rad_alphaS = (rad_alpha * *p->mult_alpha) / sqrt(S[p->M]);  // normalization
  }  // Ending bracket of changed geometry

  for (i = offset; i < nsmps; i++) {
    /* ---  Coupling reed and resonator ---
       - omits latency compensation as it is complicated on a PC
    */
    out_feedback[i] = pnew[0];

    // sound card delayed out_feeback
    // if latency based low-freq noise occures compensate in csound with hpf
    vnew[0] = p->in[i];
    update_vp(p->M, p->dt, p->dx, c, (*p->mult_rho * rho), &S[0], &pold[0], &vold[0], &pnew[0], &vnew[0]);

    // Boundary condition at tube end has radiation losses, damps traveling wave
    pnew[p->M] = (pold[p->M] * p->rad_betaS / (*p->mult_rho * rho) + vnew[p->M] - vold[p->M]) /
      (p->rad_betaS / (*p->mult_rho * rho) + p->rad_alphaS / (*p->mult_rho * rho)*p->dt);

    // Copying p(n+1) to p(n) and v(n+1) to v(n)
    // i.e. new becomes old grid for next call
    memcpy(&pold[0], &pnew[0], aux_size + 1);
    memcpy(&vold[0], &vnew[0], aux_size + 1);

    // sound is obtained at given position
    int pickup_idx = MIN((int) CEIL(*p->pickup_pos * (*p->L / p->dx)), p->M - 1);

    out_sound[i] = pnew[pickup_idx];  // Output the damped ending of the tube to csound
  }

  p->Lold      = *p->L;
  p->r_old     = *p->r;
  p->slope_old = *p->slope;
  p->Mold      = p->M;
  p->dxold     = p->dx;

  return OK;
}

static OENTRY halfphysler_localops[] =
  {
   {
    "halfphysler", (int32_t) sizeof(CONE_RADIATION_LOSSES), 0, 3, "aa", "akkkkkk",
    (SUBR)cone_radiation_losses_init,
    (SUBR)cone_radiation_losses_aperf
   }
  };

LINKAGE_BUILTIN(halfphysler_localops)
