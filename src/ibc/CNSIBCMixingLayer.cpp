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
 * \file CNSIBCMixingLayer.cpp
 *
 * \brief Member functions for CNSIBCMixingLayer
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

//----- Internal -----//

#include "CNSIBCMixingLayer.H"
#include "ChordInput.H"
#include "CRDPhysics.H"
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

/*******************************************************************************
 *
 * Class CNSIBCMixingLayer: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCMixingLayer::CNSIBCMixingLayer()
  :
  CNSIBCGeneralizedSingleBlock(),
  m_layerShift(0.),
  m_momThick(0.),
  m_rhoThick(0.),
  m_specThick(0.),
  m_loRho(m_initRho),
  m_loU(-m_initVel),
  m_loT(m_initT),
  m_singleShear(true)
{
  readBCInfo();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCMixingLayer::~CNSIBCMixingLayer()
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
CNSIBCMixingLayer::IBCName() const
{
  return "Mixing layer case";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCMixingLayer::writeIBCInfo() const
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
CNSIBCMixingLayer::initialize(LevelData<FArrayBox>&      a_U,
                              LevelGridMetrics&          a_gridMetrics,
                              const LayoutData<FluxBox>& a_unitNormals,
                              const Real                 a_time,
                              const int                  a_level) const
{
#if (CH_SPACEDIM != 1)
  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();

  // IC constants: a) _1 refers to upper layer, _0 refers to lower layer
  //               b) c_1prime is mean y-offset for upper shear layer
  const RealVect domLen(CRDparam::g_domainLength);
  D_TERM(,Real c_1prime = (0.5 + m_layerShift)*domLen[1];
          Real c_0prime = 0.5*(domLen[1] - c_1prime);,);
  if (!m_singleShear)
    {
      D_TERM(,c_1prime = 0.5*(domLen[1] + c_1prime);,);
    }
  const int  numModesHi = m_perturbHi.size();
  const int  numModesLo = m_perturbLo.size();
  const Real gamma      = CRDparam::g_gamma;
  Real       rgas       = CRDparam::g_R;
  const Real rhoRatio   = m_initRho/m_loRho;
  const Real tempRatio  = m_initT/m_loT;
  const Real vortThick  = m_momThick/2.;
  const Real c_exp      = 2.*PI;
  const Real u_0        = 0.5*(m_initVel[0] - m_loU[0]);
  const Real u_m        = 0.5*(m_initVel[0] + m_loU[0]);
  D_TERM(, Real d = 0.25*domLen[1];,);
  if (m_singleShear)
    {
      D_TERM(,d *= 2.;,);
    }
  Real mach = 0.;
  if ((m_initP > 0.) && (m_initRho > 0.))
    {
      mach = m_initVel[0]/sqrt(gamma*m_initP/m_initRho);
    }
  else
    {
      mach = m_initVel[0]/sqrt(gamma*rgas*m_initT);
    }

  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get the domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));

      // Get physical coordinates
      FABSTACKTEMP(XiFab, box2Dom, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, box2Dom, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(box2Dom, XiFab, XFab, blockCoordSys);

      // Pointwise values of U
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(0.);
      CH_assert(UFab.box().contains(box2Dom));

      // Pointwise values of W
      FABSTACKTEMP(Wc, box2Dom, numWVar);
      Wc.setVal(0.);

      // Fill Wc with data
      MD_ARRAY_RESTRICT(arrayW, Wc);
      MD_ARRAY_RESTRICT(arrayX, XFab);
      MD_BOXLOOP(box2Dom, i)
        {
          // Current physical space location
          const RealVect loc(D_DECL(arrayX[MD_IX(i, 0)],
                                    arrayX[MD_IX(i, 1)],
                                    arrayX[MD_IX(i, 2)]));
          D_TERM(,Real y_1 = loc[1];,);
          y_1 -= c_1prime;     // Define true upper shear location offset
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
          D_TERM(,Real y_0 = loc[1];,);
          y_0 -= c_0prime;     // Define true lower shear location offset
          // Real phase_xy0 = 0.; // Define low x-offset for all z-modes
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
                     // const Real modeAmp1 = m_modeAmpLo[mode][1]*domLen[0];
                     // const Real trigArg1 = m_phaseLo[mode][1]
                     // + 2.*PI*m_perturbLo[mode][1]*loc[2]/domLen[2];
                     // phase_xy0 += modeAmp1*std::sin(trigArg1);
                );
            }
          Real sgn_y0 = 1.;
          Real sgn_y1 = 1.;
          if (y_0 < 0)
            {
              sgn_y0 = -1.;
            }
          if (y_1 < 0)
            {
              sgn_y1 = -1.;
            }
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
          const Real c_d0 = 1. - (c_exp*e_0*sgn_y0*trig_s0)/d;
          const Real c_d1 = 1. - (c_exp*e_1*sgn_y1*trig_s1)/d;
          const Real c_e0 = u_0*y_0 + u_0*e_0*trig_s0;
          const Real c_e1 = u_0*y_1 + u_0*e_1*trig_s1;

          const Real u_h = u_m + u_0*c_s1*c_d1
            - (c_d1*c_e1*(c_s1*c_s1 - 1.))/vortThick;
          const Real u_l = u_m - u_0*c_s0*c_d0
            + (c_d0*c_e0*(c_s0*c_s0 - 1.))/vortThick;
          const Real v_h = -u_0*e_1*c_s1*trig_c1
            + (e_1*(c_s1*c_s1 - 1.)*trig_c1*c_e1)/vortThick;
          const Real v_l =  u_0*e_0*c_s0*trig_c0
            - (e_0*(c_s0*c_s0 - 1.)*trig_c0*c_e0)/vortThick;

          Real xVel      = 0.;
          Real yVel      = 0.;
          Real rho       = m_initRho;
          Real press     = m_initP;
          Real temp      = m_initT;

          // Velocity and density/temperature profiles
          if (m_singleShear)
            {
              xVel = u_h;
              const Real u_p  = u_0*c_s1;
              const Real u_pp = u_h - u_p;
              const Real u_c  = u_p*u_pp;
              const Real a = m_loU[0]/m_initVel[0];
              const Real b = c_s1;
              Real c_4  = 1.;

              const Real d_1 = 1. - b;
              const Real d_2 = 1. - a;
              const Real d_3 = b - a;
              if ((m_initP > 0.) && (m_initRho > 0.))
                {
                  c_4 = (1./rhoRatio)*d_1/d_2 + d_3/d_2
                    + 0.5*(gamma - 1.)*mach*mach*d_1*d_3;
                  rho   = m_initRho*(1./c_4);
                  press = m_initP - 0.5*m_initRho*(u_pp*u_pp + 2.*u_c
                                                   + v_h*v_h)*(gamma - 1.);
                }
              else if ((m_initP < 0.) && (m_initRho > 0.))
                {
                  // Not tested
                  c_4 = (1./rhoRatio)*d_1/d_2 + d_3/d_2
                    + 0.5*(gamma - 1.)*mach*mach*d_1*d_3;
                  rho  = m_initRho*(1./c_4);
                  const Real p_diff = 0.5*m_initRho*(u_pp*u_pp + 2.*u_c
                                                     + v_h*v_h)*(gamma - 1.);
                  temp = m_initT*c_4 - p_diff/(m_initRho*rgas);
                }
              else
                {
                  // Not tested
                  press = m_initP - 0.5*m_initRho*(u_pp*u_pp + 2.*u_c
                                                   + v_h*v_h)*(gamma - 1.);
                  c_4 = tempRatio*d_1/d_2 + d_3/d_2
                    + 0.5*(gamma - 1.)*mach*mach*d_1*d_3;
                  temp = m_initT*c_4;
                }
            }
          else
            {
              D_TERM(
                ,
                if (loc[1] >= 0.5*domLen[1])
                  {
                    xVel = u_h;
                    yVel = v_h;
                    const Real u_p  = u_0*c_s1;
                    const Real u_pp = u_h - u_p;
                    const Real u_c  = u_p*u_pp;
                    const Real a = m_loU[0]/m_initVel[0];
                    const Real b = c_s1;
                    Real c_4  = 1.;

                    const Real d_1 = 1. - b;
                    const Real d_2 = 1. - a;
                    const Real d_3 = b - a;
                    if ((m_initP > 0.) && (m_initRho > 0.))
                      {
                        c_4 = (1./rhoRatio)*d_1/d_2 + d_3/d_2
                          + 0.5*(gamma - 1.)*mach*mach*d_1*d_3;
                        rho   = m_initRho*(1./c_4);
                        press = m_initP - 0.5*m_initRho
                          *(u_pp*u_pp + 2.*u_c + v_h*v_h)*(gamma - 1.);
                      }
                    else if ((m_initP < 0.) && (m_initRho > 0.))
                      {
                        // Not tested
                        c_4 = (1./rhoRatio)*d_1/d_2 + d_3/d_2
                          + 0.5*(gamma - 1.)*mach*mach*d_1*d_3;
                        rho  = m_initRho*(1./c_4);
                        const Real p_diff = 0.5*m_initRho
                          *(u_pp*u_pp + 2.*u_c + v_h*v_h)*(gamma - 1.);
                        temp = m_initT*c_4 - p_diff/(m_initRho*rgas);
                      }
                    else
                      {
                        // Not tested
                        c_4 = tempRatio*d_1/d_2 + d_3/d_2
                          + 0.5*(gamma - 1.)*mach*mach*d_1*d_3;
                        press = m_initP - 0.5*m_initRho
                          *(u_pp*u_pp + 2.*u_c + v_h*v_h)*(gamma - 1.);
                        temp = m_initT*c_4;
                      }
                  }
                else
                  {
                    xVel = u_l;
                    yVel = v_l;
                    const Real u_p  = u_0*c_s0;
                    const Real u_pp = u_l - u_p;
                    const Real u_c  = u_p*u_pp;
                    const Real a = m_loU[0]/m_initVel[0];
                    const Real b = c_s0;
                    Real c_4  = 1.;

                    const Real d_1 = 1. - b;
                    const Real d_2 = 1. - a;
                    const Real d_3 = b - a;
                    if ((m_initP > 0.) && (m_initRho > 0.))
                      {
                        c_4 = (1./rhoRatio)*d_1/d_2 + d_3/d_2
                          + 0.5*(gamma - 1.)*mach*mach*d_1*d_3;
                        rho   = m_initRho*(1./c_4);
                        press = m_initP - 0.5*m_initRho
                          *(u_pp*u_pp + 2.*u_c + v_l*v_l)*(gamma - 1.);
                      }
                    else if ((m_initP < 0.) && (m_initRho > 0.))
                      {
                        // Not tested
                        c_4 = (1./rhoRatio)*d_1/d_2 + d_3/d_2
                          + 0.5*(gamma - 1.)*mach*mach*d_1*d_3;
                        rho  = m_initRho*(1./c_4);
                        const Real p_diff = 0.5*m_initRho
                          *(u_pp*u_pp + 2.*u_c + v_l*v_l)*(gamma - 1.);
                        temp = m_initT*c_4 - p_diff/(m_initRho*rgas);
                      }
                    else
                      {
                        // Not tested
                        c_4 = tempRatio*d_1/d_2 + d_3/d_2
                          + 0.5*(gamma - 1.)*mach*mach*d_1*d_3;
                        temp = m_initT*c_4;
                      }
                  }
                ,);
            }

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
 *  \param[in]  a_globalKE
 *                      Global mean kinetic energy (turbulent forcing)
 *  \param[in]  a_globalHelicity
 *                      Global mean helicity (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void
CNSIBCMixingLayer::addSourceTerm(FArrayBox&           a_sourceFab,
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
                                 const Real           a_globalKE,
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
CNSIBCMixingLayer::haveExactSol() const
{
  return false;
}

/*--------------------------------------------------------------------*/
//  Compute the exact solution state \<U\> in the cells
/**
 *  \param[out] a_Ux    Exact solution
 *  \param[in]  a_box   Cell locations to update
 *  \param[in]  a_problemDomain
 *                      The problem domain on the level
 *  \param[in]  a_time  Current solution time
 *  \param[in]  a_level Grid level
 *  \return             0 - Successfully computed exact solution
 *                      1 - Exact solution is not known
 *//*-----------------------------------------------------------------*/

int
CNSIBCMixingLayer::exactSol(FArrayBox&        a_Ux, 
                            const Box&        a_box, 
                            const Box&        a_disjointBox,
                            LevelGridMetrics& a_gridMetrics,
                            const Real        a_time,
                            const int         a_level) const
{
  return 1;
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
CNSIBCMixingLayer::readBCInfo()
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
  m_readInput = true;
}
