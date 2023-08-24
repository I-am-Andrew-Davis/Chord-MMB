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
 * \file CNSIBCSpatiallyEvolvingShear.cpp
 *
 * \brief Member functions for CNSIBCSpatiallyEvolvingShear
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelGridMetrics.H"
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"
#include "UnitNormalsF_F.H"

//----- Internal -----//

#include "CNSIBCSpatiallyEvolvingShear.H"
#include "ChordInput.H"
#include "CRDPhysics.H"
#include "CRDState.H"
#include "PatchMappedFunc.H"
#include "CRDparam.H"
#include "DataTemp.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"
#include "TagMethodValue.H"
#include "TagMethodVMS.H"
#include "CRDutil.H"
#include "ViscousTensor4thOrderOpF_F.H"
#include "LoHiCenter.H"

/*******************************************************************************
 *
 * Class CNSIBCSpatiallyEvolvingShear: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCSpatiallyEvolvingShear::CNSIBCSpatiallyEvolvingShear()
  :
  CNSIBCGeneralized(),
  m_layerShift(0.),
  m_momThick(0.),
  m_rhoThick(0.),
  m_specThick(0.),
  m_loRho(m_initRho),
  m_loU(-m_initVel),
  m_loT(m_initT),
  m_singleShear(true),
  m_reaction(false),
  m_O2comp(-1),
  m_CH4comp(-1),
  m_N2comp(-1)
{
  readBCInfo();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCSpatiallyEvolvingShear::~CNSIBCSpatiallyEvolvingShear()
{
}

/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Return a name describing the IBC
/** \return             Name of IBC
 *//*-----------------------------------------------------------------*/

const char *const
CNSIBCSpatiallyEvolvingShear::IBCName() const
{
  return "Mixing layer case";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCSpatiallyEvolvingShear::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Initialize \<U\>
/** Sets the initial state for a solution.  This routine must compute
 *  \<U\>.
 *  \param[out] a_U     State on the level to be initialized in this
 *                      routine
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_time  Initial time
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

void
CNSIBCSpatiallyEvolvingShear::initialize(
  LevelData<FArrayBox>&      a_U,
  LevelGridMetrics&          a_gridMetrics,
  const LayoutData<FluxBox>& a_unitNormals,
  const Real                 a_time,
  const int                  a_level) const
{
#if (CH_SPACEDIM != 1)
  const int numWVar  = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();

  // High and low states for the inlet condition
  const CRDState& lowState = CRDState::get(m_idxStateInletLow);
  const CRDState& highState = CRDState::get(m_idxStateInletHigh);

  // A stream-function is used to compute an initial divergence-free condition
  // psi = -xi*u_0*tanh(xi/(2*delta_theta)) + u_m
  //
  // xi = y + exp(-eta*\abs(y))*(\sum b_i*sin(2*pi*omega_i*x/L + phi_i))
  //
  // Next, we numerically differentiate psi with-respect-to y to get u-velocity
  //   and with-respect-to x to get v-velocity

  // IC constants: a) _1 refers to upper layer, _0 refers to lower layer
  //               b) c_1prime is mean y-offset for upper shear layer
  const RealVect domLen(CRDparam::g_domainLength);
  D_TERM(,Real c_1prime = (0.5 + m_layerShift)*domLen[1];
          Real c_0prime = 0.5*(domLen[1] - c_1prime);,);
  D_TERM(,c_1prime = 0.5*(domLen[1] + c_1prime);,);

  const int  numModesHi = m_perturbHi.size();
  const int  numModesLo = m_perturbLo.size();
  const Real gamma      = CRDparam::g_gamma;
  const Real vortThick  = m_momThick/2.;
  const Real c_exp      = 2.*PI;

  // State information for upper and lower portions of the inlet
  const Real rhoHi = highState.density();
  const RealVect velHi = highState.velocity();
  const Real pressHi = highState.pressure();
  const Real tempHi = highState.temperature();

  const RealVect velLo = lowState.velocity();

  const Real u_0        = 0.5*(velHi[0] - velLo[0]);
  const Real u_m        = 0.5*(velHi[0] + velLo[0]);
  const Real u_absMean  = 0.5*(std::abs(velHi[0]) + std::abs(velLo[0]));
  D_TERM(, Real d = 0.25*domLen[1];,);

  // Assume m_initP and m_initRho are always set (check this)
  Real mach = u_absMean/sqrt(gamma*pressHi/rhoHi);

  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box6Dom = grow(box, 6);
      box6Dom &= blockDomain;
      Box box5Dom = grow(box, 5);
      box5Dom &= blockDomain;
      Box box4Dom = grow(box, 4);
      box4Dom &= blockDomain;
      Box box3Dom = grow(box, 3);
      box3Dom &= blockDomain;
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));

      // Get physical coordinates
      FABSTACKTEMP(XiFab, box6Dom, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, box6Dom, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(box6Dom, XiFab, XFab, blockCoordSys);

      // Pointwise values of U
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(0.);
      CH_assert(UFab.box().contains(box2Dom));

      // Pointwise values of W
      FABSTACKTEMP(Wc, box2Dom, numWVar);
      Wc.setVal(0.);

      // Pointwise values of streamfunction psi -- this assumes two-dimensional
      //   flow with variability of x and y velocity in the z-direction
      //   Note: a discontinuity is present in the streamfunction where
      //   the upper and lower halves of the domain are joined. In order that
      //   the differentiation of this provide a smooth velocity, the separate
      //   stream definitions must be stored as separate variables
      FABSTACKTEMP(psiPntFab, box6Dom, 2);
      psiPntFab.setVal(0.);
      // Average values of streamfunction psi
      FABSTACKTEMP(psiAvgFab, box5Dom, 2);
      psiAvgFab.setVal(0.);
      // True shear locations
      FABSTACKTEMP(yLoc, box6Dom, 2);
      // Sine function summations
      FABSTACKTEMP(trigSum, box6Dom, 2);

      MD_BOXLOOP(box6Dom, i)
        {
          // Current, absolute physical-space location
          const RealVect loc(D_DECL(XFab[MD_IX(i, 0)], XFab[MD_IX(i, 1)],
                                    XFab[MD_IX(i, 2)]));
          // Define true upper shear location offset
          D_TERM(,Real y_1 = loc[1] - c_1prime;,);
          Real phase_xy1 = 0.; // Define hi x-offset for all z-modes
          for (int mode = 0; mode != numModesHi; ++mode)
            {
              D_TERM(,,
                     Real modeAmp0 = m_modeAmpHi[mode][2]*domLen[1]*0.25;
                     if (m_singleShear)
                       {
                         modeAmp0 *= 2.;
                       }
                     const Real trigArg0 = m_phaseHi[mode][2]
                     + 2.*PI*m_perturbHi[mode][2]*loc[2]/domLen[2];
                     y_1 -= modeAmp0*std::sin(trigArg0);
                     const Real modeAmp1 = m_modeAmpHi[mode][1]*domLen[0];
                     const Real trigArg1 = m_phaseHi[mode][1]
                     + 2.*PI*m_perturbHi[mode][1]*loc[2]/domLen[2];
                     phase_xy1 += modeAmp1*std::sin(trigArg1););
            }
          yLoc[MD_IX(i, 0)] = y_1;
          // Define true lower shear location offset
          D_TERM(,Real y_0 = loc[1] - c_0prime;,
                 Real phase_xy0 = 0.;); // Define low x-offset for all z-modes
          for (int mode = 0; mode != numModesLo; ++mode)
            {
              D_TERM(,,
                     Real modeAmp0 = m_modeAmpLo[mode][2]*domLen[1]*0.25;
                     if (m_singleShear)
                       {
                         modeAmp0 *= 2.;
                       }
                     const Real trigArg0 = m_phaseLo[mode][2]
                     + 2.*PI*m_perturbLo[mode][2]*loc[2]/domLen[2];
                     y_0 -= modeAmp0*std::sin(trigArg0);
                     const Real modeAmp1 = m_modeAmpLo[mode][1]*domLen[0];
                     const Real trigArg1 = m_phaseLo[mode][1]
                     + 2.*PI*m_perturbLo[mode][1]*loc[2]/domLen[2];
                     phase_xy0 += modeAmp1*std::sin(trigArg1););
            }
          yLoc[MD_IX(i, 1)] = y_0;
          Real trig_s1 = 0.; // Define upper sine-function summation
          for (int mode = 0; mode != numModesHi; ++mode)
            {
              D_TERM(,
                     const Real modeAmp = m_modeAmpHi[mode][0]*domLen[1]*0.25;
                     const Real trigArg = m_phaseHi[mode][0] + phase_xy1
                     + 2.*PI*m_perturbHi[mode][0]*loc[0]/domLen[0];
                     trig_s1 += modeAmp*std::sin(trigArg);,);
            }
          trigSum[MD_IX(i, 0)] = trig_s1;
          Real trig_s0 = 0.; // Define lower sine-function summation
          for (int mode = 0; mode != numModesLo; ++mode)
            {
              D_TERM(,
                     const Real modeAmp = m_modeAmpLo[mode][0]*domLen[1]*0.25;
                     const Real trigArg = m_phaseLo[mode][0] + phase_xy1
                     + 2.*PI*m_perturbLo[mode][0]*loc[0]/domLen[0];
                     trig_s0 += modeAmp*std::sin(trigArg);,);
            }
          trigSum[MD_IX(i, 1)] = trig_s0;

          const Real e_0  = std::exp(-c_exp*std::abs(y_0)/d);
          const Real e_1  = std::exp(-c_exp*std::abs(y_1)/d);
          const Real xi_0 = y_0 + e_0*trig_s0;
          const Real xi_1 = y_1 + e_1*trig_s1;
          const Real c_s0 = std::tanh(xi_0/vortThick);
          const Real c_s1 = std::tanh(xi_1/vortThick);

          // Always assume double-shear
          // Assign Psi point values
          // Upper definition
          psiPntFab[MD_IX(i, 0)] = xi_1*u_0*c_s1 + u_m*loc[1];
          // Lower definition
          psiPntFab[MD_IX(i, 1)] = -xi_0*u_0*c_s0 + u_m*loc[1];
        }

      // Convolve Psi to get average values
      psiAvgFab.copy(psiPntFab);
      Interval psiInterval(0,1);
      CRDutil::deconvolve(psiAvgFab, psiPntFab, box5Dom, blockDomain,
                          psiInterval, 4, -1, false);

      // Compute the gradient of Psi
      FABSTACKTEMP(psiGradAvgFab, box3Dom, 2*SpaceDim);
      psiGradAvgFab.setVal(0.);
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          Box gradBox = grow(box, 3);
          gradBox.grow(dir, 1); // pecularity of loHiCenter function
          const int MD_ID(o, dir);
          Box loBox, nextLoBox, hiBox, nextHiBox, centerBox, innerCenterBox,
            entireBox;
          int hasLo, hasHi;
          loHiCenter5(loBox, nextLoBox, hasLo, hiBox, nextHiBox, hasHi,
                      centerBox, innerCenterBox, entireBox, gradBox,
                      blockDomain, dir);
          for (int comp = 0; comp != 2; ++comp)
            {
              Real f_0 = 1./(12.*a_gridMetrics.dxVect()[dir]);
              MD_BOXLOOP(innerCenterBox, i)
                {
                  Real u_im2 = psiAvgFab[MD_OFFSETIX(i,-,2*o,comp)];
                  Real u_im1 = psiAvgFab[MD_OFFSETIX(i,-,1*o,comp)];
                  Real u_ip1 = psiAvgFab[MD_OFFSETIX(i,+,1*o,comp)];
                  Real u_ip2 = psiAvgFab[MD_OFFSETIX(i,+,2*o,comp)];
                  Real d1 = f_0*(u_im2 - 8.*u_im1 + 8.*u_ip1 - u_ip2);
                  psiGradAvgFab[MD_IX(i, comp*SpaceDim + dir)] = d1;
                }
              if (hasLo)
                {
                  MD_BOXLOOP(loBox, i)
                    {
                      Real u_i   = psiAvgFab[MD_IX(i,comp)];
                      Real u_ip1 = psiAvgFab[MD_OFFSETIX(i,+,1*o,comp)];
                      Real u_ip2 = psiAvgFab[MD_OFFSETIX(i,+,2*o,comp)];
                      Real u_ip3 = psiAvgFab[MD_OFFSETIX(i,+,3*o,comp)];
                      Real u_ip4 = psiAvgFab[MD_OFFSETIX(i,+,4*o,comp)];
                      Real d1 = f_0*(-25.*u_i + 48.*u_ip1 - 36.*u_ip2
                                     + 16.*u_ip3 - 3.*u_ip4);
                      psiGradAvgFab[MD_IX(i, comp*SpaceDim + dir)] = d1;
                    }
                  MD_BOXLOOP(nextLoBox, i)
                    {
                      Real u_im1 = psiAvgFab[MD_OFFSETIX(i,-,1*o,comp)];
                      Real u_i   = psiAvgFab[MD_IX(i,comp)];
                      Real u_ip1 = psiAvgFab[MD_OFFSETIX(i,+,1*o,comp)];
                      Real u_ip2 = psiAvgFab[MD_OFFSETIX(i,+,2*o,comp)];
                      Real u_ip3 = psiAvgFab[MD_OFFSETIX(i,+,3*o,comp)];
                      Real d1 = f_0*(-3.*u_im1 - 10.*u_i + 18.*u_ip1
                                     - 6.*u_ip2 + u_ip3);
                      psiGradAvgFab[MD_IX(i, comp*SpaceDim + dir)] = d1;
                    }
                }
              if (hasHi)
                {
                  MD_BOXLOOP(hiBox, i)
                    {
                      Real u_i   = psiAvgFab[MD_IX(i,comp)];
                      Real u_im1 = psiAvgFab[MD_OFFSETIX(i,-,1*o,comp)];
                      Real u_im2 = psiAvgFab[MD_OFFSETIX(i,-,2*o,comp)];
                      Real u_im3 = psiAvgFab[MD_OFFSETIX(i,-,3*o,comp)];
                      Real u_im4 = psiAvgFab[MD_OFFSETIX(i,-,4*o,comp)];
                      Real d1 = -f_0*(-25.*u_i + 48.*u_im1 - 36.*u_im2
                                      + 16.*u_im3 - 3.*u_im4);
                      psiGradAvgFab[MD_IX(i, comp*SpaceDim + dir)] = d1;
                    }
                  MD_BOXLOOP(nextHiBox, i)
                    {
                      Real u_ip1 = psiAvgFab[MD_OFFSETIX(i,+,1*o,comp)];
                      Real u_i   = psiAvgFab[MD_IX(i,comp)];
                      Real u_im1 = psiAvgFab[MD_OFFSETIX(i,-,1*o,comp)];
                      Real u_im2 = psiAvgFab[MD_OFFSETIX(i,-,2*o,comp)];
                      Real u_im3 = psiAvgFab[MD_OFFSETIX(i,-,3*o,comp)];
                      Real d1 = -f_0*(-3.*u_ip1 - 10.*u_i + 18.*u_im1
                                      - 6.*u_im2 + u_im3);
                      psiGradAvgFab[MD_IX(i, comp*SpaceDim + dir)] = d1;
                    }
                }
            }
        }

      // Deconvolve the gradient of psiAvgFab
      FABSTACKTEMP(psiGradPntFab, box2Dom, 2*SpaceDim);
      psiGradPntFab.setVal(0.);
      Interval psiGradInterval(0, 2*SpaceDim - 1);
      CRDutil::deconvolve(psiGradPntFab, psiGradAvgFab, box2Dom, blockDomain,
                          psiGradInterval, 4, -1, false);

      // Fill Wc with data
      MD_ARRAY_RESTRICT(arrayW, Wc);
      MD_ARRAY_RESTRICT(arrayX, XFab);
      MD_BOXLOOP(box2Dom, i)
        {
          // Current, absolute physical-space location
          const RealVect loc(D_DECL(arrayX[MD_IX(i, 0)], arrayX[MD_IX(i, 1)],
                                    arrayX[MD_IX(i, 2)]));
          // Define true upper shear location offset
          D_TERM(,Real y_1 = loc[1] - c_1prime;,);
          Real phase_xy1 = 0.; // Define hi x-offset for all z-modes
          for (int mode = 0; mode != numModesHi; ++mode)
            {
              D_TERM(,,
                     const Real modeAmp1 = m_modeAmpHi[mode][1]*domLen[0];
                     const Real trigArg1 = m_phaseHi[mode][1]
                     + 2.*PI*m_perturbHi[mode][1]*loc[2]/domLen[2];
                     phase_xy1 += modeAmp1*std::sin(trigArg1););
            }
          // Define true lower shear location offset
          D_TERM(,Real y_0 = loc[1] - c_0prime;,
                 Real phase_xy0 = 0.;); // Define low x-offset for all z-modes
          for (int mode = 0; mode != numModesLo; ++mode)
            {
              D_TERM(,,
                     const Real modeAmp1 = m_modeAmpLo[mode][1]*domLen[0];
                     const Real trigArg1 = m_phaseLo[mode][1]
                     + 2.*PI*m_perturbLo[mode][1]*loc[2]/domLen[2];
                     phase_xy0 += modeAmp1*std::sin(trigArg1););
            }
          y_1 = yLoc[MD_IX(i, 0)];
          y_0 = yLoc[MD_IX(i, 1)];
          Real trig_s1 = 0.; // Define upper sine-function summation
          Real trig_c1 = 0.; // Define upper cosine-function summation
          for (int mode = 0; mode != numModesHi; ++mode)
            {
              D_TERM(,
                     const Real modeAmp = m_modeAmpHi[mode][0]*domLen[1]*0.25;
                     const Real c_3     = 2.*m_perturbHi[mode][0]*PI/domLen[1];
                     const Real trigArg = m_phaseHi[mode][0] + phase_xy1
                     + 2.*PI*m_perturbHi[mode][0]*loc[0]/domLen[0];
                     trig_s1 += modeAmp*std::sin(trigArg);
                     trig_c1 += modeAmp*c_3*std::cos(trigArg);,);
            }
          Real trig_s0 = 0.; // Define lower sine-function summation
          Real trig_c0 = 0.; // Define lower cosine-function summation
          for (int mode = 0; mode != numModesLo; ++mode)
            {
              D_TERM(,
                     const Real modeAmp = m_modeAmpLo[mode][0]*domLen[1]*0.25;
                     const Real c_3     = 2.*m_perturbLo[mode][0]*PI/domLen[1];
                     const Real trigArg = m_phaseLo[mode][0] + phase_xy1
                     + 2.*PI*m_perturbLo[mode][0]*loc[0]/domLen[0];
                     trig_s0 += modeAmp*std::sin(trigArg);
                     trig_c0 += modeAmp*c_3*std::cos(trigArg);,);
            }

          const Real e_0  = std::exp(-c_exp*std::abs(y_0)/d);
          const Real e_1  = std::exp(-c_exp*std::abs(y_1)/d);
          const Real c_s0 = std::tanh((y_0 + e_0*trig_s0)/vortThick);
          const Real c_s1 = std::tanh((y_1 + e_1*trig_s1)/vortThick);

          Real xVel  = 0.;
          Real yVel  = 0.;
          Real rho   = rhoHi;
          Real press = pressHi;
          Real temp  = tempHi;

          // Always assume double-shear
          D_TERM(
            ,
            if (loc[1] >= 0.5*domLen[1])
              {
                xVel = psiGradPntFab[MD_IX(i,1)];
                yVel = -psiGradPntFab[MD_IX(i,0)];
                const Real a = velLo[0]/velHi[0];
                const Real b = (u_m + u_0*c_s1)/velHi[0];

                const Real d_1 = 1. - b;
                const Real d_3 = b - a;

                // Assume m_initP and m_initRho are always set (check this)
                Real c_4 = 1. + 0.5*(gamma - 1.)*mach*mach*d_1*d_3;
                press = pressHi/c_4;
                rho   = rhoHi/(c_4*c_4);
              }
            else
              {
                xVel = psiGradPntFab[MD_IX(i,SpaceDim + 1)];
                yVel = -psiGradPntFab[MD_IX(i,SpaceDim)];
                const Real a = velLo[0]/velHi[0];
                const Real b = (u_m - u_0*c_s0)/velHi[0];

                const Real d_1 = 1. - b;
                const Real d_3 = b - a;

                // Assume m_initP and m_initRho are always set (check this)
                Real c_4 = 1. + 0.5*(gamma - 1.)*mach*mach*d_1*d_3;
                press = pressHi/c_4;
                rho   = rhoHi/(c_4*c_4);
              }
            ,);

          D_TERM(arrayW[MD_IX(i,WvelIndx)]   = xVel;,
                 arrayW[MD_IX(i,WvelIndx+1)] = yVel;,
                 arrayW[MD_IX(i,WvelIndx+2)] = 0.;);
          arrayW[MD_IX(i, rhoIndx)]  = rho;
          arrayW[MD_IX(i, presIndx)] = press;
          arrayW[MD_IX(i, tempIndx)] = temp;
        }

      // Initialize the values in UFab
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box2Dom);
      fourthOrderAverageCell(UFab, blockDomain, box1Dom);
    }
#endif
}

/*--------------------------------------------------------------------*/
//  Add body force
/** \param[out] a_sourceFab
 *                      Cell-centered source terms
 *  \param[out] a_invDtFab
 *                      Inverse of the time step
 *  \param[in]  a_Wcell Cell-centered primitive variables
 *  \param[in]  a_UcellAvg
 *                      Cell-averaged conservative variables
 *  \param[in]  a_WfaceAvgFxb
 *                      Face-averaged primitive state
 *  \param[in]  a_domain
 *                      Problem domain
 *  \param[in]  a_gridMetrics
 *                      Level grid metrics
 *  \param[in]  a_time  Solution time
 *  \param[in]  a_stageWeight
 *                      Contribution of this stage to the final flux
 *                      over the time step
 *  \param[in]  a_level Current level
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[in]  a_solveBox
 *                      Box to solve over
 *  \param[in]  a_GlobalKE
 *                      Global mean kinetic energy (turbulent forcing)
 *  \param[in]  a_globalHelicity
 *                      Global mean helicity (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void
CNSIBCSpatiallyEvolvingShear::addSourceTerm(
  FArrayBox&           a_sourceFab,
  FArrayBox&           a_invDtFab,
  const FArrayBox&     a_Wcell,
  const FArrayBox&     a_UcellAvg,
  const FluxBox&       a_WfaceAvgFxb,
  const ProblemDomain& a_domain,
  LevelGridMetrics&    a_gridMetrics,
  const Real           a_time,
  const Real           a_stageWeight,
  const int            a_level,
  const Box&           a_disjointBox,
  const Box&           a_solveBox,
  const DataIndex&     a_didx,
  const Real           a_GlobalKE,
  const Real           a_globalHelicity) const
{
  return;
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCSpatiallyEvolvingShear::haveExactSol() const
{
  return false;
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set the imposed (exterior or farfield) primitive state at flow BCC
/** State is set according to reference conditions
 *  \param[in]  a_Wface Primitive state on wall from interior scheme
 *                      (face-centered)
 *  \param[out] a_Wface Primitive state at exterior to boundary
 *  \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
 *  \param[in]  a_Wcell Primitive state in cells at interior of
 *                      domain.  The cells have been shifted by half
 *                      so that the first interior layer of cells
 *                      overlaps a_Wface.
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces.  (This is the same basis used
 *                      for solving the Riemann problem on mapped
 *                      grids
 *  \param[in]  a_dir   Direction of the boundary
 *  \param[in]  a_side  LoHi side of the boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *  \param[in]  a_disjointBox
 *                      Complete current interior disjoint box
 *//*-----------------------------------------------------------------*/

void
CNSIBCSpatiallyEvolvingShear::setImposedBCprimState(
  FArrayBox&                    a_Wface,
  const Box&                    a_boundaryFaceBox,
  const FArrayBox&              a_Wcell,
  const FArrayBox&              a_unitNormalBasisFab,
  const BoundaryIndex&          a_bcIdx,
  const Box&                    a_disjointBox,
  LevelGridMetrics&             a_gridMetrics,
  const Real                    a_time,
  const int                     a_level,
  const BCInfo&                 a_bcInfo) const
{
  // Set all the component variables and intervals
  const int numSpecies       = CRDparam::g_numSpecies;
  const Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
  const int velIndx      = velIntv.begin();
  const int rhoIndx      = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx     = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx         = CRDparam::g_CRDPhysics->temperatureIndex();
  const int numWComp         = CRDparam::g_CRDPhysics->numPrimitive();
  std::vector<Real> cnBC(numSpecies);
  Real Rval = CRDparam::g_R;
  if (numSpecies > 0)
    {
      Rval = 0.;
      CRD::msg << "Spatially-evolving shear does not run with multispecies!"
               << CRD::abort;
    }

  const Side::LoHiSide a_side = a_bcIdx.m_side;
  const int a_dir = a_bcIdx.m_dir;

  // High and low states for the inlet condition
  const CRDState& meanState = CRDState::get(m_idxStateInit);
  const CRDState& lowState = CRDState::get(m_idxStateInletLow);
  const CRDState& highState = CRDState::get(m_idxStateInletHigh);

  // A stream-function is used to compute an initial divergence-free condition
  // psi = -xi*u_0*tanh(xi/(2*delta_theta)) + u_m
  //
  // xi = y + exp(-eta*\abs(y))*(\sum b_i*sin(2*pi*omega_i*x/L + phi_i))
  //
  // Next, we numerically differentiate psi with-respect-to y to get u-velocity
  //   and with-respect-to x to get v-velocity

  // IC constants: a) _1 refers to upper layer, _0 refers to lower layer
  //               b) c_1prime is mean y-offset for upper shear layer
  const RealVect domLen(CRDparam::g_domainLength);
  D_TERM(,Real c_1prime = (0.5 + m_layerShift)*domLen[1];
          Real c_0prime = 0.5*(domLen[1] - c_1prime);,);
  D_TERM(,c_1prime = 0.5*(domLen[1] + c_1prime);,);

  const int  numModesHi = m_perturbHi.size();
  const int  numModesLo = m_perturbLo.size();
  const Real gamma      = CRDparam::g_gamma;
  const Real vortThick  = m_momThick/2.;
  const Real c_exp      = 2.*PI;

  // State information for upper and lower portions of the inlet
  const Real rhoHi = highState.density();
  const RealVect velHi = highState.velocity();
  const Real pressHi = highState.pressure();

  const RealVect velLo = lowState.velocity();

  const Real u_0        = 0.5*(velHi[0] - velLo[0]);
  const Real u_m        = 0.5*(velHi[0] + velLo[0]);
  const Real u_absMean  = 0.5*(std::abs(velHi[0]) + std::abs(velLo[0]));
  D_TERM(, Real d = 0.25*domLen[1];,);

  // Assume m_initP and m_initRho are always set (check this)
  Real mach = u_absMean/sqrt(gamma*pressHi/rhoHi);

  Real Temp, Pres, Rho;

  FABSTACKTEMP(nttVel, a_boundaryFaceBox, velIntv.size());
  nttVel.copy(a_Wface, a_boundaryFaceBox, velIndx,
              a_boundaryFaceBox, 0, velIntv.size());
  // Transform velocity into normal and tangent components
  FORT_FORWARDTRANSFORMF(CHF_FRA(nttVel),
                         CHF_CONST_FRA(a_unitNormalBasisFab),
                         CHF_BOX(a_boundaryFaceBox));
  switch (a_bcInfo.m_type)
    {
    case CRDparam::DomainBCTypeDirichlet:
    case CRDparam::DomainBCTypeRelaxedCBCFar:
    case CRDparam::DomainBCTypeFarfield:
    {
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      for (int comp = 0, comp_end = CRDparam::g_CRDPhysics->numPrimitive();
           comp != comp_end; ++comp)
        {
          a_Wface.setVal(state(comp), a_boundaryFaceBox, comp);
        }
      break;
    }
    case CRDparam::DomainBCTypeRelaxedCBCOut:
    {
      // Let's just target the interior value and see what happens
      const Real alpha = a_bcInfo.m_relaxCBCStateParam;
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      for (int comp = 0, comp_end = CRDparam::g_CRDPhysics->numPrimitive();
           comp != comp_end; ++comp)
        {
          MD_BOXLOOP(a_boundaryFaceBox, i)
            {
              const Real interior = a_Wcell[MD_IX(i, comp)];
              const Real exterior = state(comp);
              a_Wface[MD_IX(i, comp)] =
                alpha*interior + (1. - alpha)*exterior;
            }
        }
      break;
    }
    case CRDparam::DomainBCTypeOutflow:
    {
//**Hack for relax test
      CH_assert(!a_boundaryFaceBox.isEmpty());
      bool haveBoRelax = true;
      auto ins = m_hackRelax.insert({ a_boundaryFaceBox, nullptr });
      if (ins.second)  // Insertion took place
        {
          haveBoRelax = false;
          ins.first->second = new FArrayBox(a_boundaryFaceBox, 1);
        }
      FArrayBox& boRelaxData = *(ins.first->second);
//**End hack for relax test
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          const RealVect Vref = meanState.velocity();
          const Real p0 = meanState.pressure();
          const Real T0 = meanState.temperature();
          // Interior estimate of velocity
          RealVect intVel;
          for (const int idxVel : EachDir)
            {
              intVel[idxVel] = a_Wface[MD_IX(i, velIndx + idxVel)];
            }
          const RealVect intVel1(intVel - Vref);
          const Real intVel1Sq = stc::dot(intVel1, intVel1);
          //**FIXME Probably need to replace use of cp with h
          const Real cp = CRDparam::g_CRDPhysics->cp();
          bool inflow;
          Real intT;
          // Temperature state based on characteristics
          if (nttVel[MD_IX(i, 0)]*Side::sign(a_side) < 0.)  // Inflow
            {
              inflow = true;
              intT = T0 - intVel1Sq/(2*cp);
            }
          else                                              // Outflow
            {
              inflow = false;
              intT = a_Wface[MD_IX(i, tempIndx)];
            }
          const Real intc1 = 1 + intVel1Sq/(2*cp*intT);
          // Whether inflow or outflow, this is the assigned pressure.  Here,
          // assume a frozen mixture with gamma from the interior.
          Real p;
#if 0
          if (inflow)
            {
              p = p0*std::pow(intc1,
                              -CRDparam::g_gamma/(CRDparam::g_gamma - 1.0));
            }
          else
            {
              p= p0;  // To account for variable velocity at exit
            }
#else
          bool debug = false;
          if (MD_GETIV(i) == IntVect{ 256, 33 } && CNSIBC::s_lastRKStage)
            {
              debug = true;
              CRD::msg.setPrecFloatSN(4);
            }
          if (inflow)
            // Hard reset
            {
              p = p0*std::pow(intc1,
                              -CRDparam::g_gamma/(CRDparam::g_gamma - 1.0));
              boRelaxData[MD_IX(i, 0)] = p0;
              if (debug)
                {
                  CRD::msg << "Reset : static " << p0 << CRD::end;
                }
            }
          else
            // Relax
            {
              Real pOld = p0;
              if (haveBoRelax)
                {
                  pOld = boRelaxData[MD_IX(i, 0)];
                }
              Real pNew = a_Wface[MD_IX(i, presIndx)];
              if (debug)
                {
                  CRD::msg << "Status: ref: " << p0 << ", old: " << pOld << ", new: " << pNew;
                }
              // p0new is between p0 and pOld: Use p0new or new relax
              //   (whicheve is closer to p0)
              // p0 is between pNew and pOld: Use pNew
              // pOld is between p0 and pNew: Use pNew
              if ((pNew - p0)*(pNew - pOld) <= 0)
                {
                  if (debug)
                    {
                      CRD::msg << " R";  // New p closer to p0 than pOld
                    }
                  constexpr Real rate = 0.05;
                  const Real pOldAdj = pOld + rate*(p0 - pOld);
                  if (rate == 1.0 || (pOldAdj - p0)*(pOldAdj - pNew) <= 0)
                    {
                      if (debug)
                        {
                          CRD::msg << 'R';  // Using relaxed p
                        }
                      pNew = pOldAdj;
                    }
                  else if (debug)
                    {
                      CRD::msg << '-';  // Just using new p
                    }
                }
              else if (debug)
                {
                  CRD::msg << " N-";  // Just using new p
                }
              if (debug)
                {
                  CRD::msg << ", adj: " << pNew << CRD::end;
                           // << CNSIBC::s_lastRKStage
                           // << ' ' << haveBoRelax << CRD::end;
                }
              // Save state
              if (CNSIBC::s_lastRKStage || !haveBoRelax)
                {
                  boRelaxData[MD_IX(i, 0)] = pNew;
                }
              p = pNew;
            }
#endif
          a_Wface[MD_IX(i, presIndx)] = p;
          // Temperature and density assigned only if inflow
          if (nttVel[MD_IX(i, 0)]*Side::sign(a_side) < 0.)  // Have inflow
            {
              const Real T = T0/intc1;
              a_Wface[MD_IX(i, rhoIndx)]  = p/(Rval*T);
              a_Wface[MD_IX(i, tempIndx)] = T;
            }
        }
      break;
    }
    case CRDparam::DomainBCTypeRelaxedCBCIn:
    case CRDparam::DomainBCTypeInflow:
    {
      // Set time-dependent inlet for spatially-evolving shear-layer

      Box WfaceAvgBox = a_boundaryFaceBox;
      Box WcellAvgBox = WfaceAvgBox;
      WcellAvgBox.grow(a_dir, 1); // Add an extra face to switch to a cell-box
      WcellAvgBox.enclosedCells(a_dir); // Change the box to a cell-box
      WcellAvgBox.grow(a_dir, 2); // Now grow the box by two, normal to the face
      Box PsiCellAvgBox = grow(WcellAvgBox, 2);
      Box PsiCellPntBox = grow(PsiCellAvgBox, 1);

      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys =
        *(a_gridMetrics.getCoordSys(a_disjointBox));

      // Get physical coordinates
      FABSTACKTEMP(XiFab, PsiCellPntBox, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, PsiCellPntBox, SpaceDim);   // Physical coordinates
      FABSTACKTEMP(XiFaceFab, WfaceAvgBox, SpaceDim);
      FABSTACKTEMP(XFaceFab, WfaceAvgBox, SpaceDim);
      this->CNSIBC::getCellCoordinates(PsiCellPntBox,
                                       XiFab,
                                       XFab,
                                       blockCoordSys);
      this->CNSIBC::getFaceCoordinates(WfaceAvgBox,
                                       XiFaceFab,
                                       XFaceFab,
                                       a_dir,
                                       blockCoordSys);

      // Pointwise values of streamfunction psi -- this assumes two-dimensional
      //   flow with variability of x and y velocity in the z-direction
      //   Note: a discontinuity is present in the streamfunction where
      //   the upper and lower halves of the domain are joined. In order that
      //   the differentiation of this provide a smooth velocity, the separate
      //   stream definitions must be stored as separate variables
      FABSTACKTEMP(psiPntFab, PsiCellPntBox, 2);
      psiPntFab.setVal(0.);
      // Average values of streamfunction psi
      FABSTACKTEMP(psiAvgFab, PsiCellAvgBox, 2);
      psiAvgFab.setVal(0.);
      // Cell-centered values of primitive states
      FABSTACKTEMP(WcellPntFab, PsiCellPntBox, numWComp);
      WcellPntFab.setVal(0.);
      // Cell-averaged values of primitive states
      FABSTACKTEMP(WcellAvgFab, PsiCellAvgBox, numWComp);
      WcellAvgFab.setVal(0.);

      // Set cell-centered values of psi and primitive states
      MD_BOXLOOP(PsiCellPntBox, i)
        {
          // The x-coordinate must be based on time
          Real t_0 = u_m*a_time;
          // Current, absolute physical-space location
          const RealVect loc(D_DECL(XFab[MD_IX(i, 0)] - t_0, XFab[MD_IX(i, 1)],
                                    XFab[MD_IX(i, 2)]));
          // Define true upper shear location offset
          D_TERM(,Real y_1 = loc[1] - c_1prime;,);
          Real phase_xy1 = 0.; // Define hi x-offset for all z-modes
          for (int mode = 0; mode != numModesHi; ++mode)
            {
              D_TERM(,,
                     Real modeAmp0 = m_modeAmpHi[mode][2]*domLen[1]*0.25;
                     if (m_singleShear)
                       {
                         modeAmp0 *= 2.;
                       }
                     const Real trigArg0 = m_phaseHi[mode][2]
                     + 2.*PI*m_perturbHi[mode][2]*loc[2]/domLen[2];
                     y_1 -= modeAmp0*std::sin(trigArg0);
                     const Real modeAmp1 = m_modeAmpHi[mode][1]*domLen[0];
                     const Real trigArg1 = m_phaseHi[mode][1]
                     + 2.*PI*m_perturbHi[mode][1]*loc[2]/domLen[2];
                     phase_xy1 += modeAmp1*std::sin(trigArg1););
            }
          // Define true lower shear location offset
          D_TERM(,Real y_0 = loc[1] - c_0prime;,
                 Real phase_xy0 = 0.;); // Define low x-offset for all z-modes
          for (int mode = 0; mode != numModesLo; ++mode)
            {
              D_TERM(,,
                     Real modeAmp0 = m_modeAmpLo[mode][2]*domLen[1]*0.25;
                     if (m_singleShear)
                       {
                         modeAmp0 *= 2.;
                       }
                     const Real trigArg0 = m_phaseLo[mode][2]
                     + 2.*PI*m_perturbLo[mode][2]*loc[2]/domLen[2];
                     y_0 -= modeAmp0*std::sin(trigArg0);
                     const Real modeAmp1 = m_modeAmpLo[mode][1]*domLen[0];
                     const Real trigArg1 = m_phaseLo[mode][1]
                     + 2.*PI*m_perturbLo[mode][1]*loc[2]/domLen[2];
                     phase_xy0 += modeAmp1*std::sin(trigArg1););
            }
          Real trig_s1 = 0.; // Define upper sine-function summation
          for (int mode = 0; mode != numModesHi; ++mode)
            {
              D_TERM(,
                     const Real modeAmp = m_modeAmpHi[mode][0]*domLen[1]*0.25;
                     const Real trigArg = m_phaseHi[mode][0] + phase_xy1
                     + 2.*PI*m_perturbHi[mode][0]*loc[0]/domLen[0];
                     trig_s1 += modeAmp*std::sin(trigArg);,);
            }
          Real trig_s0 = 0.; // Define lower sine-function summation
          for (int mode = 0; mode != numModesLo; ++mode)
            {
              D_TERM(,
                     const Real modeAmp = m_modeAmpLo[mode][0]*domLen[1]*0.25;
                     const Real trigArg = m_phaseLo[mode][0] + phase_xy1
                     + 2.*PI*m_perturbLo[mode][0]*loc[0]/domLen[0];
                     trig_s0 += modeAmp*std::sin(trigArg);,);
            }

          const Real e_0  = std::exp(-c_exp*std::abs(y_0)/d);
          const Real e_1  = std::exp(-c_exp*std::abs(y_1)/d);
          const Real xi_0 = y_0 + e_0*trig_s0;
          const Real xi_1 = y_1 + e_1*trig_s1;
          const Real c_s0 = std::tanh(xi_0/vortThick);
          const Real c_s1 = std::tanh(xi_1/vortThick);

          // Assign Psi point values
          // Upper definition
          psiPntFab[MD_IX(i, 0)] = xi_1*u_0*c_s1 + u_m*y_0;
          // Lower definition
          psiPntFab[MD_IX(i, 1)] = -xi_0*u_0*c_s0 + u_m*y_0;

          // Additional primitive states
          D_TERM(
            ,
            if (loc[1] >= 0.5*domLen[1])
              {
                const Real a = velLo[0]/velHi[0];
                const Real b = (u_m + u_0*c_s1)/velHi[0];
                Real c_4  = 1.;

                const Real d_1 = 1. - b;
                const Real d_3 = b - a;

                // Assume m_initP and m_initRho are always set (check this)
                c_4 = 1. + 0.5*(gamma - 1.)*mach*mach*d_1*d_3;
                Pres = pressHi/c_4;
                Rho  = rhoHi/(c_4*c_4);
                Temp = c_4*pressHi/(rhoHi*CRDparam::g_R);
              }
            else
              {
                const Real a = velLo[0]/velHi[0];
                const Real b = (u_m - u_0*c_s0)/velHi[0];
                Real c_4  = 1.;

                const Real d_1 = 1. - b;
                const Real d_3 = b - a;

                // Assume m_initP and m_initRho are always set (check this)
                c_4 = 1. + 0.5*(gamma - 1.)*mach*mach*d_1*d_3;
                Pres = pressHi/c_4;
                Rho  = rhoHi/(c_4*c_4);
                Temp = c_4*pressHi/(rhoHi*CRDparam::g_R);
              }
            ,);
          WcellPntFab[MD_IX(i, rhoIndx)]  = Rho;
          WcellPntFab[MD_IX(i, presIndx)] = Pres;
          WcellPntFab[MD_IX(i, tempIndx)] = Temp;
        }

      psiAvgFab.copy(psiPntFab);
      // Convolve to get cell-averaged values of psi
      for (int comp = 0; comp != 2; ++comp)
        {
          Real factor = 1./24.;
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              const int MD_ID(o, dir);
              MD_BOXLOOP(PsiCellAvgBox, i)
                {
                  Real u_i  = psiPntFab[MD_IX(i, comp)];
                  Real u_ip = psiPntFab[MD_OFFSETIX(i, +, o, comp)];
                  Real u_im = psiPntFab[MD_OFFSETIX(i, -, o, comp)];
                  Real d2 = u_ip - 2.*u_i + u_im;
                  psiAvgFab[MD_IX(i,comp)] += factor*d2;
                }
            }
        }

      // Compute gradient to get cell-averaged velocity
      FABSTACKTEMP(psiGradAvgFab, WcellAvgBox, 2*SpaceDim);
      psiGradAvgFab.setVal(0.);
      for (int comp = 0; comp != 2; ++comp)
        {
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              const int MD_ID(o, dir);
              Real f_0 = 1./(12.*a_gridMetrics.dxVect()[dir]);
              MD_BOXLOOP(WcellAvgBox, i)
                {
                  Real u_im2 = psiAvgFab[MD_OFFSETIX(i,-,2*o,comp)];
                  Real u_im1 = psiAvgFab[MD_OFFSETIX(i,-,1*o,comp)];
                  Real u_ip1 = psiAvgFab[MD_OFFSETIX(i,+,1*o,comp)];
                  Real u_ip2 = psiAvgFab[MD_OFFSETIX(i,+,2*o,comp)];
                  Real d1 = f_0*(u_im2 - 8.*u_im1 + 8.*u_ip1 - u_ip2);
                  psiGradAvgFab[MD_IX(i, comp*SpaceDim + dir)] = d1;
                }
            }
        }

      // Convolve to get cell-averaged values of additional primitive states
      WcellAvgFab.copy(WcellPntFab);
      for (int comp = 0; comp != numWComp; ++comp)
        {
          Real factor = 1./24.;
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              const int MD_ID(o, dir);
              MD_BOXLOOP(PsiCellAvgBox, i)
                {
                  Real u_i  = WcellPntFab[MD_IX(i, comp)];
                  Real u_ip = WcellPntFab[MD_OFFSETIX(i, +, o, comp)];
                  Real u_im = WcellPntFab[MD_OFFSETIX(i, -, o, comp)];
                  Real d2 = u_ip - 2.*u_i + u_im;
                  WcellAvgFab[MD_IX(i,comp)] += factor*d2;
                }
            }
        }

      // Interpolate the cell-averaged velocity to face-average velocity
      FABSTACKTEMP(psiGradAvgFaceFab, WfaceAvgBox, 2*SpaceDim);
      psiGradAvgFaceFab.setVal(0.);
      for (int comp = 0; comp != 2*SpaceDim; ++comp)
        {
          const int MD_ID(o, a_dir);
          Real f_0 = 1./12.;
          MD_BOXLOOP(WfaceAvgBox, i)
            {
              Real u_im2 = psiGradAvgFab[MD_OFFSETIX(i,-,2*o,comp)];
              Real u_im1 = psiGradAvgFab[MD_OFFSETIX(i,-,1*o,comp)];
              Real u_i   = psiGradAvgFab[MD_IX(i,comp)];
              Real u_ip1 = psiGradAvgFab[MD_OFFSETIX(i,+,1*o,comp)];
              Real d0 = f_0*7.*(u_im1 + u_i) - f_0*(u_im2 + u_ip1);
              psiGradAvgFaceFab[MD_IX(i, comp)] = d0;
            }
        }

      // Interpolate the cell-averaged primitive states to face-averages
      FABSTACKTEMP(WfaceAvgFab, WfaceAvgBox, numWComp);
      WfaceAvgFab.setVal(0.);
      for (int comp = 0; comp != numWComp; ++comp)
        {
          const int MD_ID(o, a_dir);
          Real f_0 = 1./12.;
          MD_BOXLOOP(WfaceAvgBox, i)
            {
              Real u_im2 = WcellAvgFab[MD_OFFSETIX(i,-,2*o,comp)];
              Real u_im1 = WcellAvgFab[MD_OFFSETIX(i,-,1*o,comp)];
              Real u_i   = WcellAvgFab[MD_IX(i,comp)];
              Real u_ip1 = WcellAvgFab[MD_OFFSETIX(i,+,1*o,comp)];
              Real d0 = f_0*7.*(u_im1 + u_i) - f_0*(u_im2 + u_ip1);
              WfaceAvgFab[MD_IX(i, comp)] = d0;
            }
        }

      // Assign the primitive state to a_Wface
      MD_BOXLOOP(WfaceAvgBox, i)
        {
          RealVect extVel;
          if (XFaceFab[MD_IX(i, 1)] >= 0.5*domLen[1])
            {
              extVel = RealVect{  psiGradAvgFaceFab[MD_IX(i, 1)],
                                 -psiGradAvgFaceFab[MD_IX(i, 0)],
                                  0. };
            }
          else
            {
              extVel = RealVect{  psiGradAvgFaceFab[MD_IX(i, SpaceDim + 1)],
                                 -psiGradAvgFaceFab[MD_IX(i, SpaceDim)],
                                  0. };
            }
          if (CRDparam::DomainBCTypeRelaxedCBCIn & a_bcInfo.m_type)
            {
              for (const int idxVel : EachDir)
                {
                  a_Wface[MD_IX(i, velIndx + idxVel)] = extVel[idxVel];
                }
              a_Wface[MD_IX(i, presIndx)] = WfaceAvgFab[MD_IX(i, presIndx)];
              a_Wface[MD_IX(i, rhoIndx)] = WfaceAvgFab[MD_IX(i, rhoIndx)];
              a_Wface[MD_IX(i, tempIndx)] = a_Wface[MD_IX(i, presIndx)]/
                (a_Wface[MD_IX(i, rhoIndx)]*Rval);
            }
          else
            {
              // Save interior estimate of velocity and set imposed velocity
              RealVect intVel;
              for (const int idxVel : EachDir)
                {
                  intVel[idxVel] = a_Wface[MD_IX(i, velIndx + idxVel)];
                  a_Wface[MD_IX(i, velIndx + idxVel)] = extVel[idxVel];
                }
              const RealVect Vref = meanState.velocity();
              // Pressure is determined from the interior.  First, find a
              // stagnation pressure for the given reference frame.  Assume
              // gas is frozen at upstream conditions.
              const RealVect intVel1(intVel - Vref);
              const Real intVel1Sq = stc::dot(intVel1, intVel1);
              const Real cp = CRDparam::g_CRDPhysics->cp();
              const Real intT = a_Wface[MD_IX(i, tempIndx)];
              const Real intc1 = 1 + intVel1Sq/(2*cp*intT);
              const Real p0 = a_Wface[MD_IX(i, presIndx)]*
                std::pow(intc1, CRDparam::g_gamma/(CRDparam::g_gamma - 1.));
              // Find imposed pressure from stagnation pressure
              const RealVect extVel1(extVel - Vref);
              const Real extVel1Sq = stc::dot(extVel1, extVel1);
              const Real extT =  WfaceAvgFab[MD_IX(i, tempIndx)];
              const Real extc1 = 1 + extVel1Sq/(2*cp*extT);
              a_Wface[MD_IX(i, presIndx)] =
                p0*std::pow(extc1, -gamma/(gamma - 1.0));
              // Inflow cannot become outflow, so always set density/temperature
              a_Wface[MD_IX(i, rhoIndx)] = WfaceAvgFab[MD_IX(i, rhoIndx)];
              //**FIXME need to be more careful about Rval for mixtures (it
              //**      changes per cell)
              a_Wface[MD_IX(i, tempIndx)] = a_Wface[MD_IX(i, presIndx)]/
                (a_Wface[MD_IX(i, rhoIndx)]*Rval);
            }
        }
      break;
    }
    default:
      CH_assert(false);
    }
}

/*==============================================================================
 * Private member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Read any information related to the IBC from input
/** \note
 *  <ul>
 *    <li> No output should be printed aside from errors
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
CNSIBCSpatiallyEvolvingShear::readBCInfo()
{
  ParmParse ppIBC("ibc");
  ppIBC.query("lowLayer_density", m_loRho);
  std::vector<Real> inputInitVel(SpaceDim, 0.);
  ppIBC.queryarr("lowLayer_velocity", inputInitVel, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_loU.dataPtr(),
                                           &inputInitVel.front());
  ppIBC.query("lowLayer_temperature", m_loT);
  ppIBC.query("shearLayer_position", m_layerShift);
  ppIBC.query("momentum_thickness", m_momThick);
  ppIBC.query("density_transition_thickness", m_rhoThick);
  ppIBC.query("species_transition_thickness", m_specThick);
  ppIBC.query("single_shear_layer", m_singleShear);
  ppIBC.query("combustion", m_reaction);
  
  const int numModesHi = std::max(ppIBC.countval("modes_xy_hi"),
                                  std::max(ppIBC.countval("modes_x_hi"),
                                           ppIBC.countval("modes_z_hi")));
  const int numModesLo = std::max(ppIBC.countval("modes_xy_lo"),
                                  std::max(ppIBC.countval("modes_x_lo"),
                                           ppIBC.countval("modes_z_lo")));
  m_perturbHi.resize(numModesHi);
  m_perturbHi.assign(numModesHi,IntVect::Zero);
  m_perturbLo.resize(numModesLo);
  m_perturbLo.assign(numModesLo,IntVect::Zero);
  m_phaseHi.resize(numModesHi);
  m_phaseHi.assign(numModesHi,RealVect::Zero);
  m_phaseLo.resize(numModesLo);
  m_phaseLo.assign(numModesLo,RealVect::Zero);
  m_modeAmpHi.resize(numModesHi);
  m_modeAmpHi.assign(numModesHi,RealVect::Zero);
  m_modeAmpLo.resize(numModesLo);
  m_modeAmpLo.assign(numModesLo,RealVect::Zero);
  for (int mode = 0; mode != numModesHi; ++mode)
    {
      IntVect modeWavelength(IntVect::Zero);
      D_TERM(ppIBC.query("modes_x_hi", modeWavelength[0],mode);,
             ppIBC.query("modes_xy_hi",modeWavelength[1],mode);,
             ppIBC.query("modes_z_hi", modeWavelength[2],mode););
      RealVect modePhaseShift(RealVect::Zero);
      D_TERM(ppIBC.query("mode_phaseShift_x_hi", modePhaseShift[0],mode);,
             ppIBC.query("mode_phaseShift_xy_hi",modePhaseShift[1],mode);,
             ppIBC.query("mode_phaseShift_z_hi", modePhaseShift[2],mode););
      RealVect modeAmplitude(RealVect::Zero);
      D_TERM(ppIBC.query("mode_amplitude_x_hi", modeAmplitude[0],mode);,
             ppIBC.query("mode_amplitude_xy_hi",modeAmplitude[1],mode);,
             ppIBC.query("mode_amplitude_z_hi", modeAmplitude[2],mode););
      m_perturbHi[mode] = modeWavelength;
      m_phaseHi[mode]   = modePhaseShift;
      m_modeAmpHi[mode] = modeAmplitude;
    }
  for (int mode = 0; mode != numModesLo; ++mode)
    {
      IntVect modeWavelength(IntVect::Zero);
      D_TERM(ppIBC.query("modes_x_lo", modeWavelength[0],mode);,
             ppIBC.query("modes_xy_lo",modeWavelength[1],mode);,
             ppIBC.query("modes_z_lo", modeWavelength[2],mode););
      RealVect modePhaseShift(RealVect::Zero);
      D_TERM(ppIBC.query("mode_phaseShift_x_lo", modePhaseShift[0],mode);,
             ppIBC.query("mode_phaseShift_xy_lo",modePhaseShift[1],mode);,
             ppIBC.query("mode_phaseShift_z_lo", modePhaseShift[2],mode););
      RealVect modeAmplitude(RealVect::Zero);
      D_TERM(ppIBC.query("mode_amplitude_x_lo", modeAmplitude[0],mode);,
             ppIBC.query("mode_amplitude_xy_lo",modeAmplitude[1],mode);,
             ppIBC.query("mode_amplitude_z_lo", modeAmplitude[2],mode););
      m_perturbLo[mode] = modeWavelength;
      m_phaseLo[mode]   = modePhaseShift;
      m_modeAmpLo[mode] = modeAmplitude;
    }
  
  // Retrieve specific components (FIXME)
  if (m_reaction)
    {
      int numSpecies = CRDparam::g_numSpecies;
      for(int i = 0; i != numSpecies; ++i)
        {
          if(CRDparam::g_speciesNames[i] == "CH4")
            {
              m_CH4comp = i;
            }
          else if(CRDparam::g_speciesNames[i] == "O2")
            {
              m_O2comp = i;
            }
          else if(CRDparam::g_speciesNames[i] == "N2")
            {
              m_N2comp = i;
            }
        }
      if (m_O2comp < 0 || m_CH4comp < 0)
        {
          CRD::msg << "Input (MixingLayerReaction IBC): Must have O2 and CH4!"
                   << CRD::error;
        }

      // If the solution uses multispecies physics
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          // Specify the initial mass fractions in the flow field
          m_hiMassFractions.resize(numSpecies);
          m_hiMassFractions.assign(numSpecies,0.);
          // Call function to assign mass fraction values
          int massCheck = assignMassFractions(m_hiMassFractions,
                                              "hi_specs",
                                              "hi_mfs");
          if (massCheck == 1)
            {
              CRD::msg << "Input (MixingLayerReaction IBC): 'hi_mfs' "
                       << "must be equal to 1." << CRD::error;
            }
        }
      // If the solution uses multispecies physics
      if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
        {
          // Specify the initial mass fractions in the flow field
          m_loMassFractions.resize(numSpecies);
          m_loMassFractions.assign(numSpecies,0.);
          // Call function to assign mass fraction values
          int massCheck = assignMassFractions(m_loMassFractions,
                                              "lo_specs",
                                              "lo_mfs");
          if (massCheck == 1)
            {
              CRD::msg << "Input (MixingLayerReaction IBC): 'lo_mfs' "
                       << "must be equal to 1." << CRD::error;
            }
        }
    }

  // States for high and low sections of inlet condition
  const char* initStateNameHdrHigh = "initial_state_inlet_high";
  const char* initStateNameHdrLow  = "initial_state_inlet_low";
  if (ppIBC.contains(initStateNameHdrHigh))
    {
      std::string initStateName;
      ppIBC.get(initStateNameHdrHigh, initStateName);
      m_idxStateInletHigh = CRDState::nameIndex(initStateName);
    }
  else
    {
      CRD::msg << "Input (spatially-evolving shear): " << initStateNameHdrHigh
               << " must be given!" << CRD::error;
    }
  if (ppIBC.contains(initStateNameHdrLow))
    {
      std::string initStateName;
      ppIBC.get(initStateNameHdrLow, initStateName);
      m_idxStateInletLow = CRDState::nameIndex(initStateName);
    }
  else
    {
      CRD::msg << "Input (spatially-evolving shear): " << initStateNameHdrLow
               << " must be given!" << CRD::error;
    }

  m_readInput = true;
}
