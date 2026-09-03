#include <stddef.h>  // size_t
#include <stdio.h>
#include <stdlib.h>  // calloc, free, exit, EXIT_FAILURE
#include <omp.h>

/*--------------------------------------------------------------------
  compute_mass_cg_and_cg_sens

  Purpose
  -------
  Computes the structural mass, global center of gravity (CG),
  and sensitivities of the CG coordinates with respect to the
  physical element densities used in topology optimization.

  The routine accounts for both active and passive elements.
  Passive element IDs are supplied through passiveIDs and are
  used only for mass bookkeeping and reporting purposes.
  The total structural mass includes contributions from both
  active and passive regions.

  Method
  ------
  1. Construct a lookup table identifying passive elements.
  2. Compute:
         M  = Σ (ρ_mat ρ_i V_i)

         Mx = Σ (ρ_mat ρ_i V_i x_i)
         My = Σ (ρ_mat ρ_i V_i y_i)
         Mz = Σ (ρ_mat ρ_i V_i z_i)

     where:
         ρ_mat = material density
         ρ_i   = physical element density
         V_i   = element volume

  3. Compute global center of gravity:

         CGx = Mx / M
         CGy = My / M
         CGz = Mz / M

  4. Compute sensitivities:

         dCGx/dρ_k =
             (ρ_mat V_k / M)(x_k - CGx)

         dCGy/dρ_k =
             (ρ_mat V_k / M)(y_k - CGy)

         dCGz/dρ_k =
             (ρ_mat V_k / M)(z_k - CGz)

  Parallelization
  ---------------
  The mass and first-moment accumulations are parallelized
  using OpenMP reductions. Sensitivity evaluations are also
  parallelized since each element sensitivity is independent.

  Inputs
  ------
  ne          : Number of elements.
  eleVol      : Element volumes [m^3].
  rhoPhys     : Physical element densities.
  elCG        : Element centroids stored as
                [x1,y1,z1,x2,y2,z2,...].
  mat_dens    : Material density [kg/m^3].
  passiveIDs  : Array containing passive element IDs (1-based).
  numPassive  : Number of passive elements.

  Outputs
  -------
  M           : Total structural mass [kg].
  cgx,cgy,cgz : Global center of gravity coordinates.
  dCGx_dRho   : Sensitivity of CGx w.r.t. rhoPhys.
  dCGy_dRho   : Sensitivity of CGy w.r.t. rhoPhys.
  dCGz_dRho   : Sensitivity of CGz w.r.t. rhoPhys.

  Notes
  -----
  - If M <= 0, the CG is set to zero and all sensitivities
    are returned as zero.
  - Passive element mass is reported separately for
    bookkeeping purposes.
--------------------------------------------------------------------*/
void compute_mass_cg_and_cg_sens(
    size_t ne,
    const double *eleVol,
    const double *rhoPhys,
    const double *elCG,
    double *M,
    double *cgx, double *cgy, double *cgz,
    double *dCGx_dRho,
    double *dCGy_dRho,
    double *dCGz_dRho,
    const double mat_dens,
    const int *passiveIDs,
    const int numPassive)
  {
    // ---------- 1) Compute total mass and first moments ----------
    double mass  = 0.0;   // Σ m_i
    double Mx    = 0.0;   // Σ m_i * x_i
    double My    = 0.0;   // Σ m_i * y_i
    double Mz    = 0.0;   // Σ m_i * z_i

    double passiveMass = 0.0;

    /* Passive-element lookup table */
    int *is_skin = calloc(ne, sizeof(int));


    if (!is_skin)
    {
        perror("calloc failed for is_skin");
        exit(EXIT_FAILURE);
    }

    /* passiveIDs are 1-based */
    if (passiveIDs && numPassive > 0)
    {
      #pragma omp parallel for schedule(static)
      for (int k = 0; k < numPassive; ++k)
      {
        int id = passiveIDs[k] - 1;
        if (id >= 0 && id < (int)ne)
          {
            is_skin[id] = 1;
          }
      }
    }

    #pragma omp parallel for schedule(static) reduction(+:mass,Mx,My,Mz,passiveMass)
    for (size_t i = 0; i < ne; ++i) 
    {
      /* Mass depends on element density, material density and element volume */
      const double mi = mat_dens * rhoPhys[i] * eleVol[i];   // m_i = ρ_i V_i
      mass += mi;
      Mx   += mi * elCG[3*i + 0];
      My   += mi * elCG[3*i + 1];
      Mz   += mi * elCG[3*i + 2];

      if (is_skin[i])
      {
          passiveMass += mat_dens * eleVol[i];
      }
    }

    // Guard: no mass → zero CG and zero sensitivities
    if (mass <= 0.0) 
    {
        *M   = 0.0;
        *cgx = 0.0; *cgy = 0.0; *cgz = 0.0;

        /* Only fill sensitivities if non-NULL*/
        if (dCGx_dRho)
        {
          #pragma omp parallel for schedule(static)
          for (size_t k = 0; k < ne; ++k) 
          {
            dCGx_dRho[k] = 0.0;
            dCGy_dRho[k] = 0.0;
            dCGz_dRho[k] = 0.0;
          }
        }
        free(is_skin);
        return;
    }

    // ---------- 2) Compute global CG ----------
    const double invM = 1.0 / mass;
    const double CGx  = Mx * invM;
    const double CGy  = My * invM;
    const double CGz  = Mz * invM;

    /* Return the total mass */
    *M   = mass;
    *cgx = CGx; 
    *cgy = CGy; 
    *cgz = CGz;

    printf("Total mass:           %12.3f kg\n", mass);
    printf("Passive element mass: %12.3f kg\n", passiveMass);
    printf("Active element mass:  %12.3f kg\n", mass - passiveMass);
    fflush(stdout);


    // ---------- 3) Sensitivities of CG wrt rho_k ----------
    /* NOTE: CG is computed using the full structure, including passive elements
    However, passive/skin element sensitivities are set to zero */

  if (dCGx_dRho)
  {
    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < ne; ++k) 
    {
      if (is_skin[k])
      {
          dCGx_dRho[k] = 0.0;
          dCGy_dRho[k] = 0.0;
          dCGz_dRho[k] = 0.0;
      }
      else
      {
          const double scale = mat_dens * eleVol[k] * invM;

          dCGx_dRho[k] = scale * (elCG[3*k + 0] - CGx);
          dCGy_dRho[k] = scale * (elCG[3*k + 1] - CGy);
          dCGz_dRho[k] = scale * (elCG[3*k + 2] - CGz);
      }
    }
  }
  free(is_skin);
}

/*--------------------------------------------------------------------
  compute_mass

  Purpose
  -------
  Computes the total structural mass from the physical element
  densities and element volumes. The routine additionally
  computes and reports the mass contribution from passive
  elements.

  Method
  ------
  The total mass is computed as

      M = Σ (ρ_mat ρ_i V_i)

  where

      ρ_mat = material density
      ρ_i   = physical element density
      V_i   = element volume

  Passive elements are identified using passiveIDs and their
  mass contribution is accumulated separately as

      M_passive = Σ (ρ_mat V_i)

  assuming passive elements remain fully solid
  (ρ_i = 1.0).

  The active-domain mass is then

      M_active = M - M_passive

  Parallelization
  ---------------
  The mass accumulation is parallelized using OpenMP
  reductions. Passive-element identification is performed
  through a lookup table constructed once at the beginning
  of the routine.

  Inputs
  ------
  ne          : Number of elements.
  eleVol      : Element volumes [m^3].
  rhoPhys     : Physical element densities.
  mat_dens    : Material density [kg/m^3].
  passiveIDs  : Array containing passive element IDs (1-based).
  numPassive  : Number of passive elements.

  Outputs
  -------
  M           : Total structural mass [kg].

  Notes
  -----
  - The total mass includes both active and passive elements.
  - Passive and active masses are printed separately for
    verification and bookkeeping.
  - The routine is intended for topology optimization
    applications where passive regions remain fixed and
    fully solid throughout the optimization process.
--------------------------------------------------------------------*/


void compute_mass(
    size_t ne,
    const double *eleVol,
    const double *rhoPhys,
    const double mat_dens,
    const int *passiveIDs,
    const int numPassive)
{
    // ---------- 1) Compute total mass and first moments ----------
    double mass  = 0.0;   // Σ m_i

    double passiveMass = 0.0;

    /* Passive-element lookup table */
    int *is_skin = calloc(ne, sizeof(int));


    if (!is_skin)
    {
        perror("calloc failed for is_skin");
        exit(EXIT_FAILURE);
    }

    /* passiveIDs are 1-based */
    
    if (passiveIDs && numPassive > 0)
    {
      #pragma omp parallel for schedule(static)
      for (int k = 0; k < numPassive; ++k)
      {
        int id = passiveIDs[k] - 1;
        if (id >= 0 && id < (int)ne)
          {
            is_skin[id] = 1;
          }
      }
    }


    #pragma omp parallel for schedule(static) reduction(+:mass,passiveMass)
    for (size_t i = 0; i < ne; ++i) 
    {
      const double mi = mat_dens *rhoPhys[i]* eleVol[i];   // m_i = ρ_i V_i
      //const double mi = 2700;
      mass += mi;

      if (is_skin[i])
      {
        passiveMass += mat_dens * eleVol[i];
      }
    }

    printf("Total mass:           %12.3f kg\n", mass);
    printf("Passive element mass: %12.3f kg\n", passiveMass);
    printf("Active element mass:  %12.3f kg\n", mass - passiveMass);
    fflush(stdout);
    
    free(is_skin);
}
