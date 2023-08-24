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
 * \file DCFlattening.cpp
 *
 * \brief Member functions for DCFlattening
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "LoHiCenter.H"
#include "MOLUtilities.H"
#include "GodunovUtilitiesF_F.H"
#include "DebugOut.H"

//----- Internal -----//

#include "DataTemp.H"
#include "CRDPhysics.H"
#include "CRDparam.H"
#include "CRDmsg.H"
#include "DCFlattening.H"
#include "CRDmsg.H"
#include "CNSIBC.H"
#include "CRDutil.H"

/*******************************************************************************
 *
 * Class DCFlattening: member definitions
 *
 ******************************************************************************/

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default Constructor
/**
 *//*-----------------------------------------------------------------*/

DCFlattening::DCFlattening(const Real& a_nlvTol,
                           const Real& a_fcorTol)
  :
  m_nlvTol(a_nlvTol),
  m_fcorTol(a_fcorTol)
{
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

DCFlattening::~DCFlattening()
{
}

/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Apply DC flattening
/** \param[in]  a_box   Box to limit values over
 *  \param[in]  a_HOFab High-order variables
 *  \param[out] a_HOFab Flattened value where necessary
 *  \param[in]  a_LOFab Low-order variables
 *  \param[in]  a_numComp
 *                      Number of components to correct
 *  \param[in]  a_useExtraChecks
 *                      Check to apply changes based on c_n and velocity
 *//*-----------------------------------------------------------------*/

void
DCFlattening::directDCFlatten(const Box&       a_box,
                              FArrayBox&       a_HOFab,
                              const FArrayBox& a_LOFab,
                              const int        a_numComp,
                              const bool       a_useExtraChecks) const
{
  // Leave intervals as empty unless extra checks are used
  Interval specIntv;
  Interval velIntv;
  if (a_useExtraChecks)
    {
      specIntv = CRDparam::g_CRDPhysics->speciesPrimInterval();
      velIntv = CRDparam::g_CRDPhysics->velocityInterval();
    }
  CH_assert(a_HOFab.contains(a_box));
  CH_assert(a_LOFab.contains(a_box));
  for (int comp = 0; comp != a_numComp; ++comp)
    {
      MD_BOXLOOP(a_box, i)
        {
          Real varHO = a_HOFab[MD_IX(i, comp)];
          Real varLO = a_LOFab[MD_IX(i, comp)];
          Real c1 = 1.E-16;
          if (velIntv.contains(comp))
            {
              c1 = 1.;
            }
          if (varLO > 1.E100 || varLO != varLO)
            {
              IntVect iervect(D_DECL(i0,i1,i2));
              CRD::msg << "DCFlattening::directDCFlatten: Error "
                       << "in non-linear solution at " << iervect
                       << " in box " << a_box << " for test variable "
                       << comp << " at value " << varLO << CRD::error;
            }
          if (!specIntv.contains(comp))
            {
              Real tol = m_nlvTol*(std::abs(varLO) + c1);
              Real ftest = std::abs(varHO - varLO);
              if (ftest > tol && ftest > 1.E-12)
                {
                  a_HOFab[MD_IX(i, comp)] = varLO;
                }
            }
          else
            {
              Real tol = m_nlvTol*(std::abs(varLO) + std::abs(varHO) + c1);
              Real ftest = std::abs(varHO - varLO);
              if (ftest > tol && ftest > 1.E-12)
                {
                  a_HOFab[MD_IX(i, comp)] = varLO;
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compare HO and LO flux values and flatten when necessary
/** \param[in]  a_Wminus
 *                      High-order limited minus side face state 
 *  \param[out] a_Wminus
 *                      Applied FCOR to minus side face state
 *  \param[in]  a_Wplus High-order limited plus side face state 
 *  \param[out] a_Wplus Applied FCOR to plus side face state
 *  \param[in]  a_W     High-order cell-averaged state
 *  \param[in]  a_WLO   Second-order cell-averaged state
 *  \param[in]  a_WminusLO
 *                      Second-order limited minus side face state
 *  \param[in]  a_WplusLO
 *                      Second-order limited plus side face state
 *  \param[in]  a_numSlopes
 *                      Number of slopes to apply FCOR
 *  \param[in]  a_dir   Direction of faces
 *  \param[in]  a_box   Box to apply FCOR over
 *  \param[in]  a_domain
 *                      Problem domain
 *//*-----------------------------------------------------------------*/

void
DCFlattening::ppmFCOR(FArrayBox&           a_Wminus,
                      FArrayBox&           a_Wplus,
                      const FArrayBox&     a_W,
                      const FArrayBox&     a_WLO,
                      const FArrayBox&     a_WminusLO,
                      const FArrayBox&     a_WplusLO,
                      const int&           a_numSlopes,
                      const int&           a_dir,
                      const Box&           a_box,
                      const ProblemDomain& a_domain) const
{
  const Real tol = m_fcorTol;
  Box box1cells = grow(a_box, BASISV(a_dir));
  int hasLo, hasHi;
  Box loBox, hiBox, centerBox, entireBox;
  loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox, box1cells,
             a_domain, a_dir);
  const int MD_ID(o, a_dir);
  MD_BOXLOOP(centerBox, i)
    {
      for (int c = 0; c != a_numSlopes; ++c)
        {
          Real W_i = a_W[MD_IX(i, c)];
          Real W_ip1 = a_W[MD_OFFSETIX(i,+,o,c)];
          Real W_im1 = a_W[MD_OFFSETIX(i,-,o,c)];
          Real W2_i = a_WLO[MD_IX(i, c)];
          Real W2_ip1 = a_WLO[MD_OFFSETIX(i,+,o,c)];
          Real W2_im1 = a_WLO[MD_OFFSETIX(i,-,o,c)];
          Real Wleft = a_Wminus[MD_IX(i, c)] + W_i;
          Real Wright = a_Wplus[MD_IX(i, c)] + W_i;
          Real WleftLO = a_WminusLO[MD_IX(i, c)] + W2_i;
          Real WrightLO = a_WplusLO[MD_IX(i, c)] + W2_i;
          // Find local maximum values
          Real maxEr = std::max({std::abs(W_i), std::abs(W_ip1),
                std::abs(W_im1)});
          {
            Real Er = std::abs(WleftLO - Wleft);
            Real D1 = (-9.*W_im1 - 8.*Wleft + 18.*W_i - W_ip1)/12.;
            Real D1LO = (-9.*W2_im1 - 8.*WleftLO + 18.*W2_i - W2_ip1)/12.;
            Real maxD1 = std::max(std::abs(D1), std::abs(D1LO));
            if (std::abs(D1 - D1LO) > maxD1*tol &&
                Er > maxEr*tol && maxD1 > 1.2E-14*maxEr + 1.E-12)
              {
                a_Wminus[MD_IX(i, c)] = a_WminusLO[MD_IX(i, c)];
              }
          }
          {
            Real Er = std::abs(WrightLO - Wright);
            Real D1 = (9.*W_ip1 + 8.*Wright - 18.*W_i + W_im1)/12.;
            Real D1LO = (9.*W2_ip1 + 8.*WrightLO - 18.*W2_i + W2_im1)/12.;
            Real maxD1 = std::max(std::abs(D1), std::abs(D1LO));
            if (std::abs(D1 - D1LO) > maxD1*tol &&
                Er > maxEr*tol && maxD1 > 1.2E-14*maxEr + 1.E-12)
              {
                a_Wplus[MD_IX(i, c)] = a_WplusLO[MD_IX(i, c)];
              }
          }
        }
    }
  if (hasLo)
    {
      MD_BOXLOOP(loBox, i)
        {
          for (int c = 0; c != a_numSlopes; ++c)
            {
              Real W_i = a_W[MD_IX(i, c)];
              Real W_ip1 = a_W[MD_OFFSETIX(i,+,o,c)];
              Real W_ip2 = a_W[MD_OFFSETIX(i,+,2*o,c)];
              Real W2_i = a_WLO[MD_IX(i, c)];
              Real W2_ip1 = a_WLO[MD_OFFSETIX(i,+,o,c)];
              Real W2_ip2 = a_WLO[MD_OFFSETIX(i,+,2*o,c)];
              Real Wleft = a_Wminus[MD_IX(i, c)] + W_i;
              Real Wright = a_Wplus[MD_IX(i, c)] + W_i;
              Real WleftLO = a_WminusLO[MD_IX(i, c)] + W2_i;
              Real WrightLO = a_WplusLO[MD_IX(i, c)] + W2_i;
              Real maxEr = std::max(std::abs(W_i), std::abs(W_ip1));
              {
                Real Er = std::abs(WleftLO - Wleft);
                Real D1 = (-184.*Wleft + 225.*W_i - 50.*W_ip1 + 9.*W_ip2)/60.;
                Real D1LO = (-184.*WleftLO + 225.*W2_i - 50.*W2_ip1
                             + 9.*W2_ip2)/60.;
                Real maxD1 = std::max(std::abs(D1), std::abs(D1LO));
                if (std::abs(D1 - D1LO) > maxD1*tol &&
                    Er > maxEr*tol && maxD1 > 1.2E-14*maxEr + 1.E-12)
                  {
                    a_Wminus[MD_IX(i, c)] = a_WminusLO[MD_IX(i, c)];
                  }
              }
              {
                Real Er = std::abs(WrightLO - Wright);
                Real D1 = (Wleft - 6.*W_i + 3.*Wright + 2.*W_ip1)/3.;
                Real D1LO = (WleftLO - 6.*W2_i + 3.*WrightLO + 2.*W2_ip1)/3.;
                Real maxD1 = std::max(std::abs(D1), std::abs(D1LO));
                if (std::abs(D1 - D1LO) > maxD1*tol &&
                    Er > maxEr*tol && maxD1 > 1.2E-14*maxEr + 1.E-12)
                  {
                    a_Wplus[MD_IX(i, c)] = a_WplusLO[MD_IX(i, c)];
                  }
              }
            }
        }
    }
  if (hasHi)
    {
      MD_BOXLOOP(hiBox, i)
        {
          for (int c = 0; c != a_numSlopes; ++c)
            {
              Real W_i = a_W[MD_IX(i, c)];
              Real W_im1 = a_W[MD_OFFSETIX(i,-,o,c)];
              Real W_im2 = a_W[MD_OFFSETIX(i,-,2*o,c)];
              Real W2_i = a_WLO[MD_IX(i, c)];
              Real W2_im1 = a_WLO[MD_OFFSETIX(i,-,o,c)];
              Real W2_im2 = a_WLO[MD_OFFSETIX(i,-,2*o,c)];
              Real Wleft = a_Wminus[MD_IX(i, c)] + W_i;
              Real Wright = a_Wplus[MD_IX(i, c)] + W2_i;
              Real WleftLO = a_WminusLO[MD_IX(i, c)] + W2_i;
              Real WrightLO = a_WplusLO[MD_IX(i, c)] + W2_i;
              Real maxEr = std::max(std::abs(W_i), std::abs(W_im1));
              {
                Real Er = std::abs(WleftLO - Wleft);
                Real D1 = (-Wright + 6.*W_i - 3.*Wleft - 2.*W_im1)/3.;
                Real D1LO = (-WrightLO + 6.*W2_i - 3.*WleftLO - 2.*W2_im1)/3.;
                Real maxD1 = std::max(std::abs(D1), std::abs(D1LO));
                if (std::abs(D1 - D1LO) > maxD1*tol &&
                    Er > maxEr*tol && maxD1 > 1.2E-14*maxEr + 1.E-12)
                  {
                    a_Wminus[MD_IX(i, c)] = a_WminusLO[MD_IX(i, c)];
                  }
              }
              {
                Real Er = std::abs(WrightLO - Wright);
                Real D1 = (184.*Wright - 225.*W_i + 50.*W_im1 - 9.*W_im2)/60.;
                Real D1LO = (184.*WrightLO - 225.*W2_i + 50.*W2_im1
                             - 9.*W2_im2)/60.;
                Real maxD1 = std::max(std::abs(D1), std::abs(D1LO));
                if (std::abs(D1 - D1LO) > maxD1*tol &&
                    Er > maxEr*tol && maxD1 > 1.2E-14*maxEr + 1.E-12)
                  {
                    a_Wplus[MD_IX(i, c)] = a_WplusLO[MD_IX(i, c)];
                  }
              }
            }
        }
    }
}

/*==============================================================================
 * Private member functions
 *============================================================================*/
