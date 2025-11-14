#include <stddef.h>  // size_t

/*--------------------------------------------------------------
  compute_mass_cg_and_cg_sens

  Purpose:
    1) Compute total mass M and global center of gravity (cgx,cgy,cgz)
       from element volumes, densities, and centroids.
    2) Compute the sensitivities of each CG component w.r.t. each
       element's physical density rho_k:
         dCGx/d(rho_k) = (V_k / M) * (x_k - CGx)
         dCGy/d(rho_k) = (V_k / M) * (y_k - CGy)
         dCGz/d(rho_k) = (V_k / M) * (z_k - CGz)

  Inputs:
    ne        - number of elements
    eleVol    - [ne] element volumes V_i
    rhoPhys   - [ne] element densities rho_i
    elCG      - [3*ne] element centroids:
                elCG[3*i+0] = x_i, elCG[3*i+1] = y_i, elCG[3*i+2] = z_i

  Outputs:
    *M        - total mass (Σ rho_i * V_i)
    *cgx,*cgy,*cgz - global CG components
    dCGx_dRho - [ne] sensitivities of CGx w.r.t. rho_i
    dCGy_dRho - [ne] sensitivities of CGy w.r.t. rho_i
    dCGz_dRho - [ne] sensitivities of CGz w.r.t. rho_i

  Notes:
    - Units must be consistent (e.g., V in m^3, rho in kg/m^3 → mass kg).
    - If M <= 0, CG is set to (0,0,0) and sensitivities are zero.
----------------------------------------------------------------*/
void compute_mass_cg_and_cg_sens(
    size_t ne,
    const double *eleVol,
    const double *rhoPhys,
    const double *elCG,
    double *M,
    double *cgx, double *cgy, double *cgz,
    double *dCGx_dRho,
    double *dCGy_dRho,
    double *dCGz_dRho)
{
    // ---------- 1) Compute total mass and first moments ----------
    double mass  = 0.0;   // Σ m_i
    double Mx    = 0.0;   // Σ m_i * x_i
    double My    = 0.0;   // Σ m_i * y_i
    double Mz    = 0.0;   // Σ m_i * z_i

    for (size_t i = 0; i < ne; ++i) 
    {
        const double mi = rhoPhys[i] * eleVol[i];   // m_i = ρ_i V_i
        mass += mi;
        Mx   += mi * elCG[3*i + 0];
        My   += mi * elCG[3*i + 1];
        Mz   += mi * elCG[3*i + 2];
    }

    // Guard: no mass → zero CG and zero sensitivities
    if (mass <= 0.0) 
    {
        *M   = 0.0;
        *cgx = 0.0; *cgy = 0.0; *cgz = 0.0;

        /* Only fill sensitivities if non-NULL*/
        if (dCGx_dRho)
        {
          for (size_t k = 0; k < ne; ++k) 
          {
            dCGx_dRho[k] = 0.0;
            dCGy_dRho[k] = 0.0;
            dCGz_dRho[k] = 0.0;
          }
        }
        return;
    }

    // ---------- 2) Compute global CG ----------
    const double invM = 1.0 / mass;
    const double CGx  = Mx * invM;
    const double CGy  = My * invM;
    const double CGz  = Mz * invM;

    *M   = mass;
    *cgx = CGx; 
    *cgy = CGy; 
    *cgz = CGz;

    // ---------- 3) Sensitivities of CG wrt rho_k ----------
    // d(CG)/d(rho_k) = (V_k / M) * (r_k - CG)
    // Only fill if not NULL
    if (dCGx_dRho)
    {
      for (size_t k = 0; k < ne; ++k) 
      {
          const double scale = eleVol[k] * invM;        // V_k / M
          dCGx_dRho[k] = scale * (elCG[3*k + 0] - CGx); // x_k - CGx
          dCGy_dRho[k] = scale * (elCG[3*k + 1] - CGy); // y_k - CGy
          dCGz_dRho[k] = scale * (elCG[3*k + 2] - CGz); // z_k - CGz
      }
    }
}
