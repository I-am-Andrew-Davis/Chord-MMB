#ifdef CH_LANG_CC
/*
 *      _______               __
 *     / ___/ /  ___  __  ___/ /
 *    / /__/ _ \/ _ \/ _\/ _  /
 *    \___/_//_/\___/_/  \_._/
 *    Please refer to Copyright.txt, in Chord's root directory.
 */
#endif

/******************************************************************************/
/**
 * \file RiemannFunctions.H
 *
 * \brief Contains the functions for the Riemann solver and other
 *        misc functions
 *//*+*************************************************************************/

#ifndef _RIEMANNFUNCTIONS_H_
#define _RIEMANNFUNCTIONS_H_

//----- Standard Library -----//

#include <string>
#include <vector>

//----- Chombo Library -----//

#include "REAL.H"
#include "RootSolver.H"
#include "UsingNamespace.H"


/*==============================================================================
 * Definition of static variables
 *============================================================================*/

// Tolerance for exact Riemann solver
Real rTol = 3.E-15;
// Ratio test for adaptive Riemann solver
Real Quser = 2.;
// Tolerance test for adaptive Riemann solver
Real aTol = 1.E-4;


/*==============================================================================
 * Riemann functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Functions that are used in the Riemann solution
/** General nomenclature for these functions. Some might have a number
 *  to signify that two different values for one variable must be input
 *  \param[in]  a_a     Speed of sound
 *  \param[in]  a_u     Velocity
 *  \param[in]  a_gamma Specific heat ratio
 *  \param[in]  a_p     Pressure
 *  Note: the functions called RM mean the wave is right moving and LM 
 *        mean the wave is left moving
 *//*-----------------------------------------------------------------*/

// Equation for the pressure behind an expansion wave given the velocity
// spread, right moving wave
Real EWaveP2RM(const Real& a_a,
               const Real& a_p,
               const Real& a_u1,
               const Real& a_gamma)
{
  Real pstar = a_p*pow((1. - (a_gamma - 1.)*a_u1/(2.*a_a)),
                       (2.*a_gamma/(a_gamma - 1.)));
  return pstar;
}

// Equation for the pressure behind an expansion wave given the velocity
// spread, left moving wave
Real EWaveP2LM(const Real& a_a,
               const Real& a_p,
               const Real& a_u1,
               const Real& a_gamma)
{
  Real pstar = a_p*pow((1. - (a_gamma - 1.)*(-a_u1)/(2.*a_a)),
                       (2.*a_gamma/(a_gamma - 1.)));
  return pstar;
}

// Equation for the pressure behind a shock, right moving wave
Real PStarRM(const Real& a_a,
             const Real& a_p,
             const Real& a_u,
             const Real& a_gamma)
{
  Real ratio = a_u/a_a;
  Real pval = a_p*(1. + a_gamma*(a_gamma + 1.)*pow(ratio,2)/4. -
                   a_gamma*ratio*sqrt(1. + pow(((a_gamma + 1.)*ratio),2)/16.));
  return pval;
}

// Equation for the pressure behind a shock, left moving wave
Real PStarLM(const Real& a_a,
             const Real& a_p,
             const Real& a_u,
             const Real& a_gamma)
{
  Real ratio = a_u/a_a;
  Real pval = a_p*(1. + a_gamma*(a_gamma + 1.)*pow(ratio,2)/4. +
                   a_gamma*ratio*sqrt(1. + pow(((a_gamma + 1.)*ratio),2)/16.));
  return pval;
}

// Equation for the velocity induced behind an expansion wave, right moving wave
Real EWaveVel2RM(const Real& a_a,
                 const Real& a_p1,
                 const Real& a_u,
                 const Real& a_gamma,
                 const Real& a_p2)
{
  Real waveVel = a_u - 2.*a_a/(a_gamma - 1.)*
    (1. - pow(a_p2/a_p1,(a_gamma - 1.)/(2.*a_gamma)));
  return waveVel;
}

// Equation for the velocity induced behind an expansion wave, left moving wave
Real EWaveVel2LM(const Real& a_a,
                 const Real& a_p1,
                 const Real& a_u,
                 const Real& a_gamma,
                 const Real& a_p2)
{
  Real waveVel = a_u + 2.*a_a/(a_gamma - 1.)*
    (1. - pow(a_p2/a_p1,(a_gamma - 1.)/(2.*a_gamma)));
  return waveVel;
}

// Equation for the velocity induced behind a shock, right moving wave
Real ShockVel2RM(const Real& a_a,
                 const Real& a_p1,
                 const Real& a_u,
                 const Real& a_gamma,
                 const Real& a_p2)
{
  Real presRatio = a_p2/a_p1;
  Real waveVel = a_u + a_a/a_gamma*(presRatio - 1.)*
    sqrt((2.*a_gamma/(a_gamma + 1.))/
         (presRatio + (a_gamma - 1.)/(a_gamma + 1.)));
  return waveVel;
}

// Equation for the velocity induced behind a shock, left moving wave
Real ShockVel2LM(const Real& a_a,
                 const Real& a_p1,
                 const Real& a_u,
                 const Real& a_gamma,
                 const Real& a_p2)
{
  Real presRatio = a_p2/a_p1;
  Real waveVel = a_u - a_a/a_gamma*(presRatio - 1.)*
    sqrt((2.*a_gamma/(a_gamma + 1.))/
         (presRatio + (a_gamma - 1.)/(a_gamma + 1.)));
  return waveVel;
}

// Equation for the density behind a shock
Real ShockRho(const Real& a_rho,
              const Real& a_p1,
              const Real& a_gamma,
              const Real& a_p2)
{
  Real presRatio = a_p2/a_p1;
  Real gammaRatio = (a_gamma + 1.)/(a_gamma - 1.);
  Real rhoVal = a_rho*(1. + gammaRatio*presRatio)/(gammaRatio + presRatio);
  return rhoVal;
}

// Equation for the temperature behind a shock
Real ShockT(const Real& a_T,
            const Real& a_p1,
            const Real& a_gamma,
            const Real& a_p2)
{
  Real presRatio = a_p2/a_p1;
  Real gammaRatio = (a_gamma + 1.)/(a_gamma - 1.);
  Real TVal = a_T*(gammaRatio + presRatio)/(1. + gammaRatio*presRatio);
  return TVal;
}

// Equation for the density behind an expansion wave
Real EWaveRho(const Real& a_rho,
              const Real& a_p1,
              const Real& a_gamma,
              const Real& a_p2)
{
  Real presRatio = a_p2/a_p1;
  Real gammaRatio = 1./a_gamma;
  Real rhoVal = a_rho*pow(presRatio,gammaRatio);
  return rhoVal;
}

// Equation for the density behind an expansion wave
Real EWaveT(const Real& a_T,
            const Real& a_p1,
            const Real& a_gamma,
            const Real& a_p2)
{
  Real presRatio = a_p2/a_p1;
  Real gammaRatio = (a_gamma - 1.)/a_gamma;
  Real TVal = a_T*pow(presRatio,gammaRatio);
  return TVal;
}

// Equation for the speed of sound behind the shock
Real ShockSonic(const Real& a_a1,
                const Real& a_p1,
                const Real& a_gamma,
                const Real& a_p2)
{
  const Real presRatio = a_p2/a_p1;
  const Real gammaRatio = (a_gamma + 1)/(a_gamma - 1);
  Real aVal = a_a1*sqrt(presRatio*((gammaRatio + presRatio)/
                                   (1. + gammaRatio*presRatio)));
  return aVal;
}

// Equation for the speed of sound behind the expansion wave
Real EWaveSonic(const Real& a_a1,
                const Real& a_p1,
                const Real& a_gamma,
                const Real& a_p2)
{
  Real aVal = a_a1*pow(a_p2/a_p1,(a_gamma - 1.)/(2.*a_gamma));
  return aVal;
}

// Equation for the speed of the shock, right moving wave
Real ShockWVelRM(const Real& a_a,
                 const Real& a_p1,
                 const Real& a_u,
                 const Real& a_gamma,
                 const Real& a_p2)
{
  Real shockVel = a_u + a_a*sqrt((a_gamma + 1.)/(2.*a_gamma)*a_p2/a_p1 +
                                 (a_gamma - 1.)/(2.*a_gamma));
  return shockVel;
}

// Equation for the speed of the shock, left moving wave
Real ShockWVelLM(const Real& a_a,
                 const Real& a_p1,
                 const Real& a_u,
                 const Real& a_gamma,
                 const Real& a_p2)
{
  Real shockVel = a_u - a_a*sqrt((a_gamma + 1.)/(2.*a_gamma)*a_p2/a_p1 +
                                 (a_gamma - 1.)/(2.*a_gamma));
  return shockVel;
}

// Equation for state in expansion fan, right moving wave
// The values labeled Star are modified
void EWaveStateRM(Real&       a_rhoStar,
                  Real&       a_uStar,
                  Real&       a_pStar,
                  const Real& a_rho,
                  const Real& a_p,
                  const Real& a_a,
                  const Real& a_u,
                  const Real& a_gamma)
{
  Real gRatio1 = (a_gamma - 1.)/(a_gamma + 1.);
  Real gRatio2 = 2./(a_gamma + 1.);
  Real gRatio3 = 2./(a_gamma - 1.);
  a_rhoStar = a_rho*pow(gRatio2 - gRatio1/a_a*a_u, gRatio3);
  a_uStar = gRatio2*(-a_a + a_u/gRatio3);
  a_pStar = a_p*pow(gRatio2 - gRatio1/a_a*a_u, gRatio3*a_gamma);
}

// Equation for state in expansion fan, left moving wave
// The values labeled Star are modified
void EWaveStateLM(Real&       a_rhoStar,
                  Real&       a_uStar,
                  Real&       a_pStar,
                  const Real& a_rho,
                  const Real& a_p,
                  const Real& a_a,
                  const Real& a_u,
                  const Real& a_gamma)
{
  Real gRatio1 = (a_gamma - 1.)/(a_gamma + 1.);
  Real gRatio2 = 2./(a_gamma + 1.);
  Real gRatio3 = 2./(a_gamma - 1.);
  a_rhoStar = a_rho*pow(gRatio2 + gRatio1/a_a*a_u, gRatio3);
  a_uStar = gRatio2*(a_a + a_u/gRatio3);
  a_pStar = a_p*pow(gRatio2 + gRatio1/a_a*a_u, gRatio3*a_gamma);
}

/*==============================================================================
 * Misc functions and structs
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Non-linear function for consToPrim solution, called with the Brent solver
/**
 *  \param[in]  a_CNP   Reference to the ThermPhysics class
 *  \param[in]  a_E     Conservative energy divided by density minus 
 *                      kinetic energy and heat of formation
 *  \param[in]  a_Rgas  Gas constant based on mass fractions
 *  \param[in]  a_spec  Mass fractions
 *//*-----------------------------------------------------------------*/

struct C2PFunc
{
  C2PFunc(const ThermPhysics&                        a_CNP,
          const Real                                 a_E,
          const Real                                 a_Rgas,
          const int                                  a_numSpecies,
          const std::vector<Real, StackAlloc<Real>>& a_spec,
          std::vector<Real, StackAlloc<Real>>&       a_Pvals,
          std::vector<int, StackAlloc<int>>&         a_curL)
    :
    m_CNP(a_CNP),
    m_E(a_E),
    m_Rgas(a_Rgas),
    m_numSpecies(a_numSpecies),
    m_spec(a_spec),
    m_Pvals(a_Pvals),
    m_curL(a_curL)
    { }
  Real operator()(const stc::RVec<1>& a_T, const int& a_comp) const
  {
    CH_assert(false);  // Is this function used?
    Real T = a_T[a_comp];
    Real result = this->operator()(T);
    return result;
  }
  Real operator()(const Real& a_T) const
    {
      Real sumfn = 0.;
      m_CNP.lookupRefs(a_T, m_curL, m_Pvals);
      for (int sp = 0; sp != m_numSpecies; ++sp)
        {
          for (int nt = 0; nt != m_CNP.m_interpN; ++nt)
            {
              sumfn += m_spec[sp]*m_CNP.m_lookupH[sp][m_curL[nt]]*m_Pvals[nt];
            }
        }
      return (sumfn - m_Rgas*a_T) - m_E;
    }
  const ThermPhysics& m_CNP;
  const Real m_E;
  const Real m_Rgas;
  const int m_numSpecies;
  const std::vector<Real, StackAlloc<Real>>& m_spec;
  std::vector<Real, StackAlloc<Real>>& m_Pvals;
  std::vector<int, StackAlloc<int>>& m_curL;
};


/*--------------------------------------------------------------------*/
//  Derivative of non-linear function for consToPrim solution,
//  called with the Newton solver
/**
 *  \param[in]  a_CNP   Reference to the ThermPhysics class
 *  \param[in]  a_E     Conservative energy divided by density minus
 *                      kinetic energy and heat of formation
 *  \param[in]  a_Rgas  Gas constant based on mass fractions
 *  \param[in]  a_spec  Mass fractions
 *//*-----------------------------------------------------------------*/

struct C2PDeriv
{
  C2PDeriv(const ThermPhysics&      a_CNP,
           const Real&              a_Rgas,
           const int&               a_numSpecies,
           const std::vector<Real>& a_spec)
    :
    m_CNP(a_CNP),
    m_Rgas(a_Rgas),
    m_numSpecies(a_numSpecies),
    m_spec(a_spec)
    { }
  Real operator()(const stc::RVec<1>& a_T,
                  const int& a_comp,
                  const stc::RVec<1>& a_deriv) const
  {
    (void)a_deriv; // Suppress compiler warnings about unused var
    Real T = a_T[a_comp];
    Real result = this->operator()(T);
    return result;
  }
  Real operator()(const Real& a_T) const
    {
      Real sumfn = 0.;
      for (int sp = 0; sp != m_numSpecies; ++sp)
        {
          // 7-coefficient calculation
          if(m_CNP.m_thermFileFormat == 0)
            {
              CRD::msg << "RiemannFunctions::C2PDeriv: error "
                       << "derivative for 7-coefficient calculation not "
                       << " yet implemented" << CRD::error;
            }
          else if (m_CNP.m_thermFileFormat == 1)
            {
              if(a_T >= m_CNP.m_midLookupT)
                {
                  sumfn += m_spec[sp]*(m_CNP.m_Rn[sp]*
                                       (m_CNP.m_hnA1H[sp] +
                                        a_T*(m_CNP.m_hnA2H[sp] +
                                             a_T*(m_CNP.m_hnA3H[sp] +
                                                  a_T*(m_CNP.m_hnA4H[sp] +
                                                       a_T*m_CNP.m_hnA5H[sp])))));
                }
              else
                {
                  sumfn += m_spec[sp]*(m_CNP.m_Rn[sp]*
                                       (m_CNP.m_hnA1L[sp] +
                                        a_T*(m_CNP.m_hnA2L[sp] +
                                             a_T*(m_CNP.m_hnA3L[sp] +
                                                  a_T*(m_CNP.m_hnA4L[sp] +
                                                       a_T*m_CNP.m_hnA5L[sp])))));
                }
            }
        }
      return sumfn - m_Rgas;
    }
  const ThermPhysics& m_CNP;
  const Real& m_Rgas;
  const int& m_numSpecies;
  const std::vector<Real>& m_spec;
};

/*--------------------------------------------------------------------*/
//  Non-linear function for Riemann solution, called with the Brent solver
/**
 *  \param[in]  a_p1    Pressure on side 1
 *  \param[in]  a_p2    Pressure on side 2
 *  \param[in]  a_a1    Sound speed on side 1
 *  \param[in]  a_a2    Sound speed on side 2
 *  \param[in]  a_g1    Gamma on side 1
 *  \param[in]  a_g2    Gamma on side 2
 *  \param[in]  a_u1    Normal velocity on side 1
 *  \param[in]  a_u2    Normal velocity on side 2
 *//*-----------------------------------------------------------------*/

struct RFunc
{
  RFunc(const Real& a_p1,
        const Real& a_p2,
        const Real& a_a1,
        const Real& a_a2,
        const Real& a_g1,
        const Real& a_g2,
        const Real& a_u1,
        const Real& a_u2)
    :
    m_p1(a_p1),
    m_p2(a_p2),
    m_a1(a_a1),
    m_a2(a_a2),
    m_g1(a_g1),
    m_g2(a_g2),
    m_u1(a_u1),
    m_u2(a_u2)
    { }
  Real operator()(const Real& a_p3f) const
    {
      const Real u3f1 = (a_p3f <= m_p1) ?
        EWaveVel2RM(m_a1,m_p1,m_u1,m_g1,a_p3f) :
        ShockVel2RM(m_a1,m_p1,m_u1,m_g1,a_p3f);
      const Real u3f2 = (a_p3f <= m_p2) ?
        EWaveVel2LM(m_a2,m_p2,m_u2,m_g2,a_p3f) :
        ShockVel2LM(m_a2,m_p2,m_u2,m_g2,a_p3f);
      return u3f2 - u3f1;
    }
  Real m_p1;
  Real m_p2;
  Real m_a1;
  Real m_a2;
  Real m_g1;
  Real m_g2;
  Real m_u1;
  Real m_u2;
};

/*--------------------------------------------------------------------*/
//  Solve for the star state exactly or approximately
/** \param[in]  a_p3    Middle pressure approximated using PVRS
 *  \param[out] a_p3    Middle pressure
 *  \param[out] a_u3    Middle velocity
 *  \param[in]  a_CL    Left aL*rhoL
 *  \param[in]  a_CR    Right aR*rhoR
 *  \param[in]  a_aL    Left speed of sound
 *  \param[in]  a_aR    Left speed of sound
 *  \param[in]  a_pL    Left pressure
 *  \param[in]  a_pR    Right pressure
 *  \param[in]  a_uL    Left velocity
 *  \param[in]  a_uR    Right velocity
 *  \param[in]  a_gammaL
 *                      Left specific heat ratio
 *  \param[in]  a_gammaR
 *                      Right specific heat ratio
 *  \param[out] a_runCheck
 *                      1 if nonlinear solver was used, 0 otherwise
 *  \return     errorB  Flag if Brents method failed
 *//*-----------------------------------------------------------------*/

// Exact always solves the nonlinear solution for the middle pressure
struct ExactRStarSolve
{
  int operator()(Real&       a_p3,
                 Real&       a_u3,
                 const Real& a_CL,
                 const Real& a_CR,
                 const Real& a_aL,
                 const Real& a_aR,
                 const Real& a_pL,
                 const Real& a_pR,
                 const Real& a_uL,
                 const Real& a_uR,
                 const Real& a_gammaL,
                 const Real& a_gammaR,
                 int&        a_runCheck) const
    {
      // Solve for p_3, the pressure between waves
      Real pstarR = (a_uR >= 0.) ? EWaveP2RM(a_aR,a_pR,a_uR,a_gammaR) :
        PStarRM(a_aR,a_pR,a_uR,a_gammaR);
      Real pstarL = (a_uL <= 0.) ? EWaveP2LM(a_aL,a_pL,a_uL,a_gammaL) :
        PStarLM(a_aL,a_pL,a_uL,a_gammaL);
      Real pMin = std::min(pstarR,pstarL);
      Real pMax = std::max(pstarR,pstarL);
      // Flag for Brents method success or failing
      int errorB = 0;
      a_runCheck = 1;
      if (pMax - pMin > a_p3*rTol)
        {
          // Number of iterations within the Brent method
          int iterBrent = 0;
          const RFunc& f = RFunc(a_pR,a_pL,a_aR,a_aL,a_gammaR,a_gammaL,
                                 a_uR,a_uL);
          a_p3 = RootSolver::BrentER(iterBrent,errorB,f,pMax,pMin);
          // Solve for the velocity between waves
          a_u3 = (a_pR > a_p3) ? EWaveVel2RM(a_aR,a_pR,a_uR,a_gammaR,a_p3) :
            ShockVel2RM(a_aR,a_pR,a_uR,a_gammaR,a_p3);
        }
      else
        {
          a_runCheck = 0;
        }
      return errorB;
    }
};

// Approximate uses the PVRS method 
struct ApproxRStarSolve
{
  int operator()(Real&       a_p3,
                 Real&       a_u3,
                 const Real& a_CL,
                 const Real& a_CR,
                 const Real& a_aL,
                 const Real& a_aR,
                 const Real& a_pL,
                 const Real& a_pR,
                 const Real& a_uL,
                 const Real& a_uR,
                 const Real& a_gammaL,
                 const Real& a_gammaR,
                 int&        a_runCheck) const
    {
      a_runCheck = 0;
      return 0;
    }
};

// Adaptive method uses criteria to choose between approximate and exact
struct AdaptiveRStarSolve
{
  int operator()(Real&       a_p3,
                 Real&       a_u3,
                 const Real& a_CL,
                 const Real& a_CR,
                 const Real& a_aL,
                 const Real& a_aR,
                 const Real& a_pL,
                 const Real& a_pR,
                 const Real& a_uL,
                 const Real& a_uR,
                 const Real& a_gammaL,
                 const Real& a_gammaR,
                 int&        a_runCheck) const
    {
      // Solve for $p_M$ and $u_M$, the pressure and velocity between the waves
      Real pSMin = std::min(a_pL, a_pR);
      Real pSMax = std::max(a_pL, a_pR);
      Real test1 = pSMax/pSMin;
      Real test2 = std::abs(a_pL - a_pR)/pSMin;
      // Flag for Brents method success or failing
      int errorB = 0;
      a_runCheck = 1;
      if (test1 > Quser || test2 > aTol)
        {
          Real pstarR = (a_uR >= 0.) ? EWaveP2RM(a_aR,a_pR,a_uR,a_gammaR) :
            PStarRM(a_aR,a_pR,a_uR,a_gammaR);
          Real pstarL = (a_uL <= 0.) ? EWaveP2LM(a_aL,a_pL,a_uL,a_gammaL) :
            PStarLM(a_aL,a_pL,a_uL,a_gammaL);
          Real pMin = std::min(pstarR,pstarL);
          Real pMax = std::max(pstarR,pstarL);
          // Number of iterations within the Brent method
          int iterBrent = 0;
          if (pMax - pMin > a_p3*rTol)
            {
              const RFunc& f = RFunc(a_pR,a_pL,a_aR,a_aL,a_gammaR,a_gammaL,
                                     a_uR,a_uL);
              a_p3 = RootSolver::BrentER(iterBrent,errorB,f,pMax,pMin);
              a_u3 = (a_pR > a_p3) ? EWaveVel2RM(a_aR,a_pR,a_uR,a_gammaR,a_p3) :
                ShockVel2RM(a_aR,a_pR,a_uR,a_gammaR,a_p3);
            }
          else
            {
              a_runCheck = 0;
            }
        }
      else
        {
          a_runCheck = 0;
        }
      return errorB;
    }
};

void
approxFinalValues(std::vector<Real, StackAlloc<Real>>&       a_cnStar,
                  Real&                                      a_rhoOut,
                  Real&                                      a_uOut,
                  Real&                                      a_pOut,
                  const Real&                                a_p3,
                  const Real&                                a_u3,
                  const Real&                                a_rhoL,
                  const Real&                                a_rhoR,
                  const Real&                                a_aL,
                  const Real&                                a_aR,
                  const Real&                                a_pL,
                  const Real&                                a_pR,
                  const Real&                                a_uL,
                  const Real&                                a_uR,
                  const Real&                                a_gammaL,
                  const Real&                                a_gammaR,
                  const std::vector<Real, StackAlloc<Real>>& a_cnL,
                  const std::vector<Real, StackAlloc<Real>>& a_cnR)
{
  Real small = 1.E-6;
  Real gammaStar;
  Real rho0, u0, p0, a0, sig;
  Real uStar = a_u3;
  Real pStar = a_p3;
  if (uStar >= 0.)
    {
      rho0 = a_rhoL;
      p0 = a_pL;
      u0 = a_uL;
      a0 = a_aL;
      sig = 1.;
      gammaStar = a_gammaL;
    }
  else
    {
      rho0 = a_rhoR;
      p0 = a_pR;
      u0 = a_uR;
      a0 = a_aR;
      sig = -1.;
      gammaStar = a_gammaR;
    }
  Real rhoStar = rho0 + (pStar - p0)/(a0*a0);
  Real aStar = sqrt(std::abs(gammaStar*pStar/rhoStar));
  Real wStar = 0.5*(aStar*rhoStar + a0*rho0);
  Real spout = a0 - sig*u0;
  Real spin = aStar - sig*uStar;
  Real uShock = wStar/rhoStar - sig*uStar;
  if (pStar > p0)
    {
      spout = uShock;
      spin = uShock;
    }
  Real frac = (1. + (spout + spin)/std::max(spout - spin,small))*0.5;
  frac = std::max(0.,std::min(1., frac));
  Real finalRho = rho0 + frac*(rhoStar - rho0);
  Real finalU = u0 + frac*(uStar - u0);
  Real finalP = p0 + frac*(pStar - p0);

  if (spout <= 0.)
    {
      finalRho = rho0;
      finalU = u0;
      finalP = p0;
    }
  if (spin > 0.)
    {
      finalRho = rhoStar;
      finalU = uStar;
      finalP = pStar;
    }
  a_rhoOut = finalRho;
  a_uOut = finalU;
  a_pOut = finalP;
}

void
exactFinalValues(std::vector<Real, StackAlloc<Real>>&       a_cnStar,
                 Real&                                      a_rhoOut,
                 Real&                                      a_uOut,
                 Real&                                      a_pOut,
                 const Real&                                a_p3,
                 const Real&                                a_u3,
                 const Real&                                a_rhoL,
                 const Real&                                a_rhoR,
                 const Real&                                a_aL,
                 const Real&                                a_aR,
                 const Real&                                a_pL,
                 const Real&                                a_pR,
                 const Real&                                a_uL,
                 const Real&                                a_uR,
                 const Real&                                a_gammaL,
                 const Real&                                a_gammaR,
                 const std::vector<Real, StackAlloc<Real>>& a_cnL,
                 const std::vector<Real, StackAlloc<Real>>& a_cnR)
{
  // Final values
  Real rhoStar, uStar, pStar, aStar;
  Real sonicV; // Either speed of the shock or tail speed
  Real sonicH; // Speed of the head of expansion wave     
  if (a_u3 >= 0.)
    {
      // Check if the wave is a shock wave
      if (a_p3 > a_pL)
        {
          sonicV = ShockWVelLM(a_aL,a_pL,a_uL,a_gammaL,a_p3);
          // Check if point is before the contact surface
          if (sonicV > 0.)
            {
              rhoStar = a_rhoL;
              uStar = a_uL;
              pStar = a_pL;
              aStar = a_aL;
            }
          else
            {
              rhoStar = ShockRho(a_rhoL,a_pL,a_gammaL,a_p3);
              aStar = ShockSonic(a_aL,a_pL,a_gammaL,a_p3);
              uStar = a_u3;
              pStar = a_p3;
            }
        }
      // Check if wave is expansion fan
      else
        {
          aStar = EWaveSonic(a_aL,a_pL,a_gammaL,a_p3);
          sonicV = a_u3 - aStar;
          sonicH = a_uL - a_aL;
          // Check if point is behind the head
          if (sonicH > 0.)
            {
              rhoStar = a_rhoL;
              uStar = a_uL;
              pStar = a_pL;
              aStar = a_aL;
            }
          // Check if point is in expansion fan
          else if (sonicV > 0.)
            {
              EWaveStateLM(rhoStar,uStar,pStar,a_rhoL,a_pL,a_aL,a_uL,a_gammaL);
            }
          // Otherwise, the point is between fan and contact surface
          else
            {
              rhoStar = EWaveRho(a_rhoL,a_pL,a_gammaL,a_p3);
              uStar = a_u3;
              pStar = a_p3;
            }
        }
    }
  else
    {
      if (a_p3 > a_pR)
        {
          sonicV = ShockWVelRM(a_aR,a_pR,a_uR,a_gammaR,a_p3);
          if (sonicV < 0.)
            {
              rhoStar = a_rhoR;
              uStar = a_uR;
              pStar = a_pR;
              aStar = a_aR;
            }
          else
            {
              rhoStar = ShockRho(a_rhoR,a_pR,a_gammaR,a_p3);
              aStar = ShockSonic(a_aR,a_pR,a_gammaR,a_p3);
              uStar = a_u3;
              pStar = a_p3;
            }
        }
      else
        {
          aStar = EWaveSonic(a_aR,a_pR,a_gammaR,a_p3);
          sonicV = a_u3 + aStar;
          sonicH = a_uR + a_aR;
          if (sonicH < 0.)
            {
              rhoStar = a_rhoR;
              uStar = a_uR;
              pStar = a_pR;
              aStar = a_aR;
            }
          else if (sonicV < 0.)
            {
              EWaveStateRM(rhoStar,uStar,pStar,a_rhoR,a_pR,a_aR,a_uR,a_gammaR);
            }
          else
            {
              rhoStar = EWaveRho(a_rhoR,a_pR,a_gammaR,a_p3);
              uStar = a_u3;
              pStar = a_p3;
            }
        }
    }
  a_rhoOut = rhoStar;
  a_pOut = pStar;
  a_uOut = uStar;
}

#endif /* ! defined _RIEMANNFUNCTIONS_H_ */
