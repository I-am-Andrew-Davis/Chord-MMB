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
 * \file CNSIBCMMS.cpp
 *
 * \brief Member functions for CNSIBCMMS
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelGridMetrics.H"
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"
#include "LoHiCenter.H"

//----- Internal -----//

#include "CNSIBCMMS.H"
#include "ChordInput.H"
#include "CRDPhysics.H"
#include "CRDState.H"
#include "CRDutil.H"
#include "PatchMappedFunc.H"
#include "DataTemp.H"

/*******************************************************************************
 *
 * Class CNSIBCMMS: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCMMS::CNSIBCMMS()
  :
  CNSIBCGeneralized(),
  m_rA(0.),
  m_uA(RealVect::Zero),
  m_pA(0.),
  m_TA(0.),
  m_sOmegaRho(RealVect::Zero),
  m_tOmegaRho(0.),
  m_sPhiRho(RealVect::Zero),
  m_tPhiRho(0.),
  m_sOmegaU(RealVect::Zero),
  m_tOmegaU(0.),
  m_sPhiU(RealVect::Zero),
  m_tPhiU(0.),
  m_sOmegaV(RealVect::Zero),
  m_tOmegaV(0.),
  m_sPhiV(RealVect::Zero),
  m_tPhiV(0.),
  m_sOmegaW(RealVect::Zero),
  m_tOmegaW(0.),
  m_sPhiW(RealVect::Zero),
  m_tPhiW(0.),
  m_sOmegaP(RealVect::Zero),
  m_tOmegaP(0.),
  m_sPhiP(RealVect::Zero),
  m_tPhiP(0.),
  m_sOmegaT(RealVect::Zero),
  m_tOmegaT(0.),
  m_sPhiT(RealVect::Zero),
  m_tPhiT(0.),
  m_charLength(RealVect::Unit),
  m_charTime(1.),
  m_bcTestCase(0),
  m_exactInit(0),
  m_sixthOrderInit(0)
{
  CNSIBCMMS::readBCInfo();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCMMS::~CNSIBCMMS()
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
CNSIBCMMS::IBCName() const
{
  return "Method of manufactured solutions";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCMMS::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeTagInfo();
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
CNSIBCMMS::initialize(LevelData<FArrayBox>&      a_U,
                      LevelGridMetrics&          a_gridMetrics,
                      const LayoutData<FluxBox>& a_unitNormals,
                      const Real                 a_time,
                      const int                  a_level) const
{
  // Global indices and variables
  const int numWVar  = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  // Obtain the mean values and the wave-amplitude values
  const CRDState& mean = CRDState::get(m_idxStateInit);
  const Real rM = mean.density();
  const RealVect velM = mean.velocity();
  const Real pM = mean.pressure();
  const Real TM = mean.temperature();

  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get the domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box3Dom = grow(box, 3);
      box3Dom &= blockDomain;
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));
      // Get physical coordinates
      FABSTACKTEMP(XiFab, box3Dom, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, box3Dom, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(box3Dom, XiFab, XFab, blockCoordSys);

      // Pointwise values of U
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(0.);
      CH_assert(UFab.box().contains(box2Dom));
      // Pointwise values of W
      FABSTACKTEMP(Wc, box3Dom, numWVar);
      Wc.setVal(-1.);

      // Fill Wc with data
      MD_BOXLOOP(box3Dom, i)
        {
          // Physical location
          const RealVect loc(
            D_DECL(XFab[MD_IX(i, 0)], XFab[MD_IX(i, 1)], XFab[MD_IX(i, 2)]));
          // Spatial location normalized by characteristic length
          const RealVect sLoc = loc/m_charLength;
          // Time normalized by characteristic time
          const Real tLoc = a_time/m_charTime;
          // Spatial argument for density field
          const RealVect rhoArg = 2.*PI*m_sOmegaRho*sLoc + m_sPhiRho;
          // Temporal argument for density field
          const Real rhoTimeArg = 2.*PI*m_tOmegaRho*tLoc + m_tPhiRho;
          // Spatial + temporal arguments for velocity fields
          D_TERM(const RealVect uArg = 2.*PI*m_sOmegaU*sLoc + m_sPhiU;
                 const Real uTimeArg = 2.*PI*m_tOmegaU*tLoc + m_tPhiU;,
                 const RealVect vArg = 2.*PI*m_sOmegaV*sLoc + m_sPhiV;
                 const Real vTimeArg = 2.*PI*m_tOmegaV*tLoc + m_tPhiV;,
                 const RealVect wArg = 2.*PI*m_sOmegaW*sLoc + m_sPhiW;
                 const Real wTimeArg = 2.*PI*m_tOmegaW*tLoc + m_tPhiW;);
          // Spatial argument for pressure field
          const RealVect pArg = 2.*PI*m_sOmegaP*sLoc + m_sPhiP;
          // Temporal argument for pressure field
          const Real pTimeArg = 2.*PI*m_tOmegaP*tLoc + m_tPhiP;
          // Spatial argument for temperature field
          const RealVect TArg = 2.*PI*m_sOmegaT*sLoc + m_sPhiT;
          // Temporal argument for temperature field
          const Real TTimeArg = 2.*PI*m_tOmegaT*tLoc + m_tPhiT;
          // Define the density field
          Wc[MD_IX(i, rhoIndx)] = rM + m_rA*(std::sin(rhoTimeArg))*D_TERM(
            std::sin(rhoArg[0]),*std::sin(rhoArg[1]),*std::sin(rhoArg[2]));
          // Define the velocity fields -- first, the X-velocity field
          D_TERM(Wc[MD_IX(i, WvelIndx)]     = velM[0] +
                 m_uA[0]*(std::sin(uTimeArg))*D_TERM(
                   std::sin(uArg[0]),*std::sin(uArg[1]),*std::sin(uArg[2]));,
                 // Y-velocity field
                 Wc[MD_IX(i, WvelIndx + 1)] = velM[1] +
                 m_uA[1]*(std::sin(vTimeArg))*D_TERM(
                   std::sin(vArg[0]),*std::sin(vArg[1]),*std::sin(vArg[2]));,
                 // Z-velocity field
                 Wc[MD_IX(i, WvelIndx + 2)] = velM[2] +
                 m_uA[2]*(std::sin(wTimeArg))*D_TERM(
                   std::sin(wArg[0]),*std::sin(wArg[1]),*std::sin(wArg[2])););
          if (!m_bcTestCase)
            {
              // Define the pressure field
              Wc[MD_IX(i, presIndx)] = pM + m_pA*(std::sin(pTimeArg))*D_TERM(
                std::sin(pArg[0]),*std::sin(pArg[1]),*std::sin(pArg[2]));
            }
          else
            {
              // Define the temperature field
              Wc[MD_IX(i, tempIndx)] = TM + m_TA*(std::sin(TTimeArg))*D_TERM(
                std::sin(TArg[0]),*std::sin(TArg[1]),*std::sin(TArg[2]));
            }
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

      if (m_sixthOrderInit)
        {
          // Testing another approach
          FABSTACKTEMP(UFabTemp, box3Dom, UFab.nComp());
          UFabTemp.setVal(0.);
          CRDparam::g_CRDPhysics->initialize(UFabTemp,
                                             Wc,
                                             a_gridMetrics,
                                             a_unitNormals[dit],
                                             dit(),
                                             box,
                                             box3Dom);
          //**NOTE: Adding 6th-order convolution function just for this
          const Real d2xFactor = 1./24.;
          const Real d4xFactor = 1./1920.;
          for (int comp = 0; comp != UFab.nComp(); ++comp)
            {
              // Compute the 2nd-derivative and the 4th-derivative using 5 cells
              FABSTACKTEMP(D2Fab, box1Dom, 1);
              FABSTACKTEMP(D4Fab, box1Dom, 1);
              D2Fab.setVal(0.);
              D4Fab.setVal(0.);
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  Box loBox, nextLoBox, hiBox, nextHiBox, centerBox,
                    innerCenterBox, entireBox;
                  int hasLo, hasHi;
                  Box inBox(box1Dom);
                  inBox.grow(dir, 1);
                  loHiCenter5(loBox, nextLoBox, hasLo, hiBox, nextHiBox, hasHi,
                              centerBox, innerCenterBox, entireBox, inBox,
                              blockDomain, dir);
                  const int MD_ID(o, dir);
                  MD_BOXLOOP(innerCenterBox, i)
                    {
                      Real u_m2 = UFabTemp[MD_OFFSETIX(i,-,2*o,comp)];
                      Real u_m1 = UFabTemp[MD_OFFSETIX(i,-,o,comp)];
                      Real u    = UFabTemp[MD_IX(i,comp)];
                      Real u_p1 = UFabTemp[MD_OFFSETIX(i,+,o,comp)];
                      Real u_p2 = UFabTemp[MD_OFFSETIX(i,+,2*o,comp)];
                      D2Fab[MD_IX(i, 0)] +=
                        (-u_m2 + 16.*u_m1 - 30.*u + 16.*u_p1 - u_p2)/12.;
                      D4Fab[MD_IX(i, 0)] +=
                        u_m2 - 4.*u_m1 + 6.*u - 4*u_p1 + u_p2;
                    }
                  if (hasLo)
                    {
                      MD_BOXLOOP(loBox, i)
                        {
                          Real u    = UFabTemp[MD_IX(i,comp)];
                          Real u_p1 = UFabTemp[MD_OFFSETIX(i,+,o,comp)];
                          Real u_p2 = UFabTemp[MD_OFFSETIX(i,+,2*o,comp)];
                          Real u_p3 = UFabTemp[MD_OFFSETIX(i,+,3*o,comp)];
                          Real u_p4 = UFabTemp[MD_OFFSETIX(i,+,4*o,comp)];
                          D2Fab[MD_IX(i, 0)] += (35.*u - 104.*u_p1 + 114.*u_p2
                                                 - 56.*u_p3 + 11.*u_p4)/12.;
                          D4Fab[MD_IX(i, 0)] += u - 4.*u_p1 + 6.*u_p2 - 4.*u_p3
                            + u_p4;
                        }
                      MD_BOXLOOP(nextLoBox, i)
                        {
                          Real u_m1 = UFabTemp[MD_OFFSETIX(i,-,o,comp)];
                          Real u    = UFabTemp[MD_IX(i,comp)];
                          Real u_p1 = UFabTemp[MD_OFFSETIX(i,+,o,comp)];
                          Real u_p2 = UFabTemp[MD_OFFSETIX(i,+,2*o,comp)];
                          Real u_p3 = UFabTemp[MD_OFFSETIX(i,+,3*o,comp)];
                          D2Fab[MD_IX(i, 0)] +=
                            (11.*u_m1 - 20.*u + 6.*u_p1 + 4.*u_p2 - u_p3)/12.;
                          D4Fab[MD_IX(i, 0)] += u_m1 - 4.*u + 6.*u_p1 - 4.*u_p2
                            + u_p3;
                        }
                    }
                  if (hasHi)
                    {
                      MD_BOXLOOP(hiBox, i)
                        {
                          Real u_m4 = UFabTemp[MD_OFFSETIX(i,-,4*o,comp)];
                          Real u_m3 = UFabTemp[MD_OFFSETIX(i,-,3*o,comp)];
                          Real u_m2 = UFabTemp[MD_OFFSETIX(i,-,2*o,comp)];
                          Real u_m1 = UFabTemp[MD_OFFSETIX(i,-,o,comp)];
                          Real u    = UFabTemp[MD_IX(i,comp)];
                          D2Fab[MD_IX(i, 0)] += (11.*u_m4 - 56.*u_m3 + 114.*u_m2
                                                 - 104.*u_m1 + 35.*u)/12.;
                          D4Fab[MD_IX(i, 0)] +=
                            u_m4 - 4.*u_m3 + 6.*u_m2 - 4.*u_m1
                            + u;
                        }
                      MD_BOXLOOP(nextHiBox, i)
                        {
                          Real u_m3 = UFabTemp[MD_OFFSETIX(i,-,3*o,comp)];
                          Real u_m2 = UFabTemp[MD_OFFSETIX(i,-,2*o,comp)];
                          Real u_m1 = UFabTemp[MD_OFFSETIX(i,-,o,comp)];
                          Real u    = UFabTemp[MD_IX(i,comp)];
                          Real u_p1 = UFabTemp[MD_OFFSETIX(i,+,o,comp)];
                          D2Fab[MD_IX(i, 0)] +=
                            (-u_m3 + 4.*u_m2 + 6.*u_m1 - 20.*u + 11.*u_p1)/12.;
                          D4Fab[MD_IX(i, 0)] += u_m3 - 4.*u_m2 + 6.*u_m1 - 4.*u
                            + u_p1;
                        }
                    }
                }
              MD_BOXLOOP(box1Dom, i)
                {
                  UFab[MD_IX(i, comp)] = UFabTemp[MD_IX(i, comp)]
                    + d2xFactor*D2Fab[MD_IX(i, 0)]
                    + d4xFactor*D4Fab[MD_IX(i, 0)];
                }
            }
        }
      if (CRDparam::g_cartesian && m_exactInit)
        {
          // First, get node locations of the cells
          //**NOTE: a special box is necessary where the high-end has been
          //        grown by 1 so that the high nodes can be accounted for
          Box box2DomNodes = box2Dom;
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              box2DomNodes.growHi(dir, 1);
            }
          FABSTACKTEMP(XiFabNodes, box2DomNodes, SpaceDim);
          FABSTACKTEMP(XFabNodes, box2DomNodes, SpaceDim);
          this->CNSIBC::getNodeCoordinates(box2DomNodes,
                                           XiFabNodes,
                                           XFabNodes,
                                           blockCoordSys);
          Real R_cell = CRDparam::g_R;
          Real gamma = CRDparam::g_gamma;
          MD_BOXLOOP(box1Dom, i)
            {
              XFabNodes[MD_OFFSETIV(i,+,IntVect::Unit, 0)];
              // Physical location of high side of cell
              RealVect highLoc(
                D_DECL(XFabNodes[MD_OFFSETIV(i,+,IntVect::Unit, 0)],
                       XFabNodes[MD_OFFSETIV(i,+,IntVect::Unit, 1)],
                       XFabNodes[MD_OFFSETIV(i,+,IntVect::Unit, 2)]));
              // Physical location of low side of cell
              RealVect lowLoc(D_DECL(XFabNodes[MD_IX(i, 0)],
                                     XFabNodes[MD_IX(i, 1)],
                                     XFabNodes[MD_IX(i, 2)]));
              // Length of the cell in each direction
              RealVect length = highLoc - lowLoc;
              // Cell volume
              Real vol = length.product();
              // Spatial location normalized by characteristic length
              const RealVect sLocHigh = highLoc/m_charLength;
              const RealVect sLocLow  = lowLoc/m_charLength;
              // Time normalized by characteristic time
              const Real tLoc = a_time/m_charTime;

              // Spatial argument for density field -- high side evaluation
              const RealVect rhoArgH =
                2.*PI*m_sOmegaRho*sLocHigh + m_sPhiRho;
              // Spatial argument for density field -- low side evaluation
              const RealVect rhoArgL  =
                2.*PI*m_sOmegaRho*sLocLow + m_sPhiRho;
              // Derivative of spatial argument for density field
              const RealVect rhoArgD = 2.*PI*m_sOmegaRho;
              // Temporal argument for density field
              const Real rhoTimeArg = 2.*PI*m_tOmegaRho*tLoc + m_tPhiRho;
              // Cell-averaged density
              // Note there are 2 different cases for this integral
                // 1) m_sOmegaRho[0] = 0
                //    -- this requires just multiplying the sinusoid by length
                // 2) m_sOmegaRho[0] != 0
                //    -- this requires just integrating the sinusoid
              RealVect rho_integralDir = RealVect::Zero;
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  if (m_sOmegaRho[dir] == 0)
                    {
                      rho_integralDir[dir] =
                        std::sin(rhoArgL[dir])*length[dir];
                    }
                  else
                    {
                      rho_integralDir[dir] =
                        (std::cos(rhoArgL[dir])
                       - std::cos(rhoArgH[dir]))/rhoArgD[dir];
                    }
                }
              Real rho_cellAvg = rM
                + (m_rA/vol)*(std::sin(rhoTimeArg))*rho_integralDir.product();

              // Set the cell-averaged density to the exact solution
              UFab[MD_IX(i, rhoIndx)] = rho_cellAvg;

              // Spatial + temporal arguments for velocity fields
              D_TERM(
                const RealVect uArgH = 2.*PI*m_sOmegaU*sLocHigh + m_sPhiU;
                const RealVect uArgL = 2.*PI*m_sOmegaU*sLocLow + m_sPhiU;
                const RealVect uArgD = 2.*PI*m_sOmegaU;
                const Real uTimeArg  = 2.*PI*m_tOmegaU*tLoc + m_tPhiU;,
                const RealVect vArgH = 2.*PI*m_sOmegaV*sLocHigh + m_sPhiV;
                const RealVect vArgL = 2.*PI*m_sOmegaV*sLocLow + m_sPhiV;
                const RealVect vArgD = 2.*PI*m_sOmegaV;
                const Real vTimeArg  = 2.*PI*m_tOmegaV*tLoc + m_tPhiV;,
                const RealVect wArgH = 2.*PI*m_sOmegaW*sLocHigh + m_sPhiW;
                const RealVect wArgL = 2.*PI*m_sOmegaW*sLocLow + m_sPhiW;
                const RealVect wArgD = 2.*PI*m_sOmegaW;
                const Real wTimeArg  = 2.*PI*m_tOmegaW*tLoc + m_tPhiW;);
              // Define the spatial integrals of the velocity sinusoids
              // Note there are 2 different cases for these integrals
                // 1) m_sOmegaU[0] = 0
                //    -- this requires just multiplying the sinusoid by length
                // 2) m_sOmegaU[0] != 0
                //    -- this requires just integrating the sinusoid
              D_TERM(
                RealVect u_integralDir = RealVect::Zero;
                for (int dir = 0; dir != SpaceDim; ++dir)
                  {
                    if (m_sOmegaU[dir] == 0)
                      {
                        u_integralDir[dir] = std::sin(uArgL[dir])*length[dir];
                      }
                    else
                      {
                        u_integralDir[dir] = (std::cos(uArgL[dir])
                                            - std::cos(uArgH[dir]))/uArgD[dir];
                      }
                  },
                RealVect v_integralDir = RealVect::Zero;
                for (int dir = 0; dir != SpaceDim; ++dir)
                  {
                    if (m_sOmegaV[dir] == 0)
                      {
                        v_integralDir[dir] = std::sin(vArgL[dir])*length[dir];
                      }
                    else
                      {
                        v_integralDir[dir] = (std::cos(vArgL[dir])
                                            - std::cos(vArgH[dir]))/vArgD[dir];
                      }
                  },
                RealVect w_integralDir = RealVect::Zero;
                for (int dir = 0; dir != SpaceDim; ++dir)
                  {
                    if (m_sOmegaW[dir] == 0)
                      {
                        w_integralDir[dir] = std::sin(wArgL[dir])*length[dir];
                      }
                    else
                      {
                        w_integralDir[dir] = (std::cos(wArgL[dir])
                                            - std::cos(wArgH[dir]))/wArgD[dir];
                      }
                  });
              // Define the momentum fields -- first, the X-momentum field
              D_TERM(
                Real u_0 = rM*velM[0]; // mean rho*u component
                // Mean u*rho_sinusoid
                Real u_1 = velM[0]*(m_rA/vol)*(std::sin(rhoTimeArg))
                *rho_integralDir.product();
                // Mean rho*u_sinusoid
                Real u_2 = rM*(m_uA[0]/vol)*(std::sin(uTimeArg))
                *u_integralDir.product();
                // rho_sinusoid*u_sinusoid
                // Note there are 5 different cases for this integral
                // 1) m_sOmegaU[0] = 0 and m_sOmegaRho = 0
                //    -- this requires just multiplying the two sinusoids
                // 2) m_sOmegaU[0] = 0 and m_sOmegaRho != 0
                //    -- this requires integrating just the rho sinusoid
                // 3) m_sOmegaU[0] != 0 and m_sOmegaRho = 0
                //    -- this requires integrating just the u sinusoid
                // 4) m_sOmegaU[0] != 0 and m_sOmegaRho != 0 and
                //    m_sOmegaU[0] = m_sOmegaRho[0]
                //    -- this requires integration using same frequencies
                // 5) m_sOmegaU[0] != 0 and m_sOmegaRho != 0 and
                //    m_sOmegaU[0] != m_sOmegaRho[0]
                //    -- this requires integration using different frequencies
                RealVect rhoU_integralDir = RealVect::Zero;
                for (int dir = 0; dir != SpaceDim; ++dir)
                  {
                    if ((uArgD[dir] == 0) && (rhoArgD[dir] == 0))
                      {
                        rhoU_integralDir[dir] = std::sin(m_sPhiU[dir])
                          *std::sin(m_sPhiRho[dir])*length[dir];
                      }
                    else if ((uArgD[dir] == 0) && (rhoArgD[dir] != 0))
                      {
                        rhoU_integralDir[dir] =
                          std::sin(m_sPhiU[dir])*rho_integralDir[dir];
                      }
                    else if ((uArgD[dir] != 0) && (rhoArgD[dir] == 0))
                      {
                        rhoU_integralDir[dir] =
                          u_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                      }
                    else if (uArgD[dir] == rhoArgD[dir])
                      {
                        rhoU_integralDir[dir] = (1./(4.*rhoArgD[dir]))*(
                          2.*rhoArgD[dir]*sLocHigh[dir]
                          *std::cos(m_sPhiRho[dir] - m_sPhiU[dir])
                          - std::sin(rhoArgH[dir] + uArgH[dir])
                          - 2.*rhoArgD[dir]*sLocLow[dir]
                          *std::cos(m_sPhiRho[dir] - m_sPhiU[dir])
                          + std::sin(rhoArgL[dir] + uArgL[dir]));
                      }
                    else
                      {
                        rhoU_integralDir[dir] = 0.5*(
                        std::sin(rhoArgH[dir]-uArgH[dir])/(rhoArgD[dir]
                                                          -uArgD[dir])
                      - std::sin(rhoArgH[dir]+uArgH[dir])/(rhoArgD[dir]
                                                          +uArgD[dir])
                      - std::sin(rhoArgL[dir]-uArgL[dir])/(rhoArgD[dir]
                                                          -uArgD[dir])
                      + std::sin(rhoArgL[dir]+uArgL[dir])/(rhoArgD[dir]
                                                          +uArgD[dir]));
                      }
                  }
                Real u_4 = m_rA*(m_uA[0]/vol)
                *(std::sin(rhoTimeArg))*(std::sin(uTimeArg))
                *(rhoU_integralDir.product());
                Real rhoU_cellAvg = u_0 + u_1 + u_2 + u_4;,
                // Y-momentum field
                Real v_0 = rM*velM[1]; // mean rho*v component
                // Mean v*rho_sinusoid
                Real v_1 = velM[1]*(m_rA/vol)*(std::sin(rhoTimeArg))
                *rho_integralDir.product();
                // Mean rho*v_sinusoid
                Real v_2 = rM*(m_uA[1]/vol)*(std::sin(vTimeArg))
                *v_integralDir.product();
                // rho_sinusoid*v_sinusoid
                // Note there are 5 different cases for this integral
                // 1) m_sOmegaV[0] = 0 and m_sOmegaRho = 0
                //    -- this requires just multiplying the two sinusoids
                // 2) m_sOmegaV[0] = 0 and m_sOmegaRho != 0
                //    -- this requires integrating just the rho sinusoid
                // 3) m_sOmegaV[0] != 0 and m_sOmegaRho = 0
                //    -- this requires integrating just the v sinusoid
                // 4) m_sOmegaV[0] != 0 and m_sOmegaRho != 0 and
                //    m_sOmegaV[0] = m_sOmegaRho[0]
                //    -- this requires integration using same frequencies
                // 5) m_sOmegaV[0] != 0 and m_sOmegaRho != 0 and
                //    m_sOmegaV[0] != m_sOmegaRho[0]
                //    -- this requires integration using different frequencies
                RealVect rhoV_integralDir = RealVect::Zero;
                for (int dir = 0; dir != SpaceDim; ++dir)
                  {
                    if ((vArgD[dir] == 0) && (rhoArgD[dir] == 0))
                      {
                        rhoV_integralDir[dir] = std::sin(m_sPhiV[dir])
                          *std::sin(m_sPhiRho[dir])*length[dir];
                      }
                    else if ((vArgD[dir] == 0) && (rhoArgD[dir] != 0))
                      {
                        rhoV_integralDir[dir] =
                          std::sin(m_sPhiV[dir])*rho_integralDir[dir];
                      }
                    else if ((vArgD[dir] != 0) && (rhoArgD[dir] == 0))
                      {
                        rhoV_integralDir[dir] =
                          v_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                      }
                    else if (vArgD[dir] == rhoArgD[dir])
                      {
                        rhoV_integralDir[dir] = (1./(4.*rhoArgD[dir]))*(
                          2.*rhoArgD[dir]*sLocHigh[dir]
                          *std::cos(m_sPhiRho[dir] - m_sPhiV[dir])
                          - std::sin(rhoArgH[dir] + vArgH[dir])
                          - 2.*rhoArgD[dir]*sLocLow[dir]
                          *std::cos(m_sPhiRho[dir] - m_sPhiV[dir])
                          + std::sin(rhoArgL[dir] + vArgL[dir]));
                      }
                    else
                      {
                        rhoV_integralDir[dir] = 0.5*(
                        std::sin(rhoArgH[dir]-vArgH[dir])/(rhoArgD[dir]
                                                          -vArgD[dir])
                      - std::sin(rhoArgH[dir]+vArgH[dir])/(rhoArgD[dir]
                                                          +vArgD[dir])
                      - std::sin(rhoArgL[dir]-vArgL[dir])/(rhoArgD[dir]
                                                          -vArgD[dir])
                      + std::sin(rhoArgL[dir]+vArgL[dir])/(rhoArgD[dir]
                                                          +vArgD[dir]));
                      }
                  }
                Real v_4 = m_rA*(m_uA[1]/vol)
                *(std::sin(rhoTimeArg))*(std::sin(vTimeArg))
                *(rhoV_integralDir.product());
                Real rhoV_cellAvg = v_0 + v_1 + v_2 + v_4;,
                // Z-momentum field
                Real w_0 = rM*velM[2]; // mean rho*w component
                // Mean w*rho_sinusoid
                Real w_1 = velM[2]*(m_rA/vol)*(std::sin(rhoTimeArg))
                *rho_integralDir.product();
                // Mean rho*w_sinusoid
                Real w_2 = rM*(m_uA[2]/vol)*(std::sin(wTimeArg))
                *w_integralDir.product();
                // rho_sinusoid*w_sinusoid
                // Note there are 5 different cases for this integral
                // 1) m_sOmegaW[0] = 0 and m_sOmegaRho = 0
                //    -- this requires just multiplying the two sinusoids
                // 2) m_sOmegaW[0] = 0 and m_sOmegaRho != 0
                //    -- this requires integrating just the rho sinusoid
                // 3) m_sOmegaW[0] != 0 and m_sOmegaRho = 0
                //    -- this requires integrating just the w sinusoid
                // 4) m_sOmegaW[0] != 0 and m_sOmegaRho != 0 and
                //    m_sOmegaW[0] = m_sOmegaRho[0]
                //    -- this requires integration using same frequencies
                // 5) m_sOmegaW[0] != 0 and m_sOmegaRho != 0 and
                //    m_sOmegaW[0] != m_sOmegaRho[0]
                //    -- this requires integration using different frequencies
                RealVect rhoW_integralDir = RealVect::Zero;
                for (int dir = 0; dir != SpaceDim; ++dir)
                  {
                    if ((wArgD[dir] == 0) && (rhoArgD[dir] == 0))
                      {
                        rhoW_integralDir[dir] = std::sin(m_sPhiW[dir])
                          *std::sin(m_sPhiRho[dir])*length[dir];
                      }
                    else if ((wArgD[dir] == 0) && (rhoArgD[dir] != 0))
                      {
                        rhoW_integralDir[dir] =
                          std::sin(m_sPhiW[dir])*rho_integralDir[dir];
                      }
                    else if ((wArgD[dir] != 0) && (rhoArgD[dir] == 0))
                      {
                        rhoW_integralDir[dir] =
                          w_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                      }
                    else if (wArgD[dir] == rhoArgD[dir])
                      {
                        rhoW_integralDir[dir] = (1./(4.*rhoArgD[dir]))*(
                          2.*rhoArgD[dir]*sLocHigh[dir]
                          *std::cos(m_sPhiRho[dir] - m_sPhiW[dir])
                          - std::sin(rhoArgH[dir] + wArgH[dir])
                          - 2.*rhoArgD[dir]*sLocLow[dir]
                          *std::cos(m_sPhiRho[dir] - m_sPhiW[dir])
                          + std::sin(rhoArgL[dir] + wArgL[dir]));
                      }
                    else
                      {
                        rhoW_integralDir[dir] = 0.5*(
                        std::sin(rhoArgH[dir]-wArgH[dir])/(rhoArgD[dir]
                                                          -wArgD[dir])
                      - std::sin(rhoArgH[dir]+wArgH[dir])/(rhoArgD[dir]
                                                          +wArgD[dir])
                      - std::sin(rhoArgL[dir]-wArgL[dir])/(rhoArgD[dir]
                                                          -wArgD[dir])
                      + std::sin(rhoArgL[dir]+wArgL[dir])/(rhoArgD[dir]
                                                          +wArgD[dir]));
                      }
                  }
                Real w_4 = m_rA*(m_uA[2]/vol)
                *(std::sin(rhoTimeArg))*(std::sin(wTimeArg))
                *(rhoW_integralDir.product());
                Real rhoW_cellAvg = w_0 + w_1 + w_2 + w_4;);

              // Set the cell-averaged momentum to the exact solution
              D_TERM(UFab[MD_IX(i, WvelIndx)]   = rhoU_cellAvg;,
                     UFab[MD_IX(i, WvelIndx+1)] = rhoV_cellAvg;,
                     UFab[MD_IX(i, WvelIndx+2)] = rhoW_cellAvg;);

              // Spatial argument for pressure field -- high side evaluation
              const RealVect pArgH = 2.*PI*m_sOmegaP*sLocHigh + m_sPhiP;
              // Spatial argument for pressure field -- low side evaluation
              const RealVect pArgL  = 2.*PI*m_sOmegaP*sLocLow + m_sPhiP;
              // Derivative of spatial argument for pressure field
              const RealVect pArgD = 2.*PI*m_sOmegaP;
              // Temporal argument for pressure field
              const Real pTimeArg = 2.*PI*m_tOmegaP*tLoc + m_tPhiP;
              Real p_cellAvg = 0.;

              // Spatial argument for temperature field -- high side evaluation
              const RealVect TArgH = 2.*PI*m_sOmegaT*sLocHigh + m_sPhiT;
              // Spatial argument for temperature field -- low side evaluation
              const RealVect TArgL  = 2.*PI*m_sOmegaT*sLocLow + m_sPhiT;
              // Derivative of spatial argument for temperature field
              const RealVect TArgD = 2.*PI*m_sOmegaT;
              // Temporal argument for temperature field
              const Real TTimeArg = 2.*PI*m_tOmegaT*tLoc + m_tPhiT;

              if (!m_bcTestCase)
                {
                  // Note there are 2 different cases for this integral
                  // 1) m_sOmegaP[0] = 0
                  //    -- this requires multiplying the sinusoid by length
                  // 2) m_sOmegaP[0] != 0
                  //    -- this requires integrating the sinusoid
                  RealVect p_integralDir = RealVect::Zero;
                  for (int dir = 0; dir != SpaceDim; ++dir)
                    {
                      if (m_sOmegaP[dir] == 0)
                        {
                          p_integralDir[dir] =
                            std::sin(pArgL[dir])*length[dir];
                        }
                      else
                        {
                          p_integralDir[dir] =
                            (std::cos(pArgL[dir])
                           - std::cos(pArgH[dir]))/pArgD[dir];
                        }
                    }
                  p_cellAvg = pM
                    + (m_pA/vol)*(std::sin(pTimeArg))*p_integralDir.product();
                }
              else
                {
                  // Note there are 2 different cases for this integral
                  // 1) m_sOmegaT[0] = 0
                  //    -- this requires multiplying the sinusoid by length
                  // 2) m_sOmegaT[0] != 0
                  //    -- this requires integrating the sinusoid
                  RealVect T_integralDir = RealVect::Zero;
                  for (int dir = 0; dir != SpaceDim; ++dir)
                    {
                      if (m_sOmegaT[dir] == 0)
                        {
                          T_integralDir[dir] =
                            std::sin(TArgL[dir])*length[dir];
                        }
                      else
                        {
                          T_integralDir[dir] =
                            (std::cos(TArgL[dir])
                           - std::cos(TArgH[dir]))/TArgD[dir];
                        }
                    }
                  // Temperature*R*rho
                  Real T_0 = rM*TM; // mean rho*T component
                  // Mean T*rho_sinusoid
                  Real T_1 = TM*(m_rA/vol)*(std::sin(rhoTimeArg))
                    *rho_integralDir.product();
                  // Mean rho*T_sinusoid
                  Real T_2 = rM*(m_TA/vol)*(std::sin(TTimeArg))
                    *T_integralDir.product();
                  // rho_sinusoid*T_sinusoid
                  // Note there are 5 different cases for this integral
                  // 1) m_sOmegaT[0] = 0 and m_sOmegaRho = 0
                  //    -- this requires just multiplying the two sinusoids
                  // 2) m_sOmegaT[0] = 0 and m_sOmegaRho != 0
                  //    -- this requires integrating just the rho sinusoid
                  // 3) m_sOmegaT[0] != 0 and m_sOmegaRho = 0
                  //    -- this requires integrating just the u sinusoid
                  // 4) m_sOmegaT[0] != 0 and m_sOmegaRho != 0 and
                  //    m_sOmegaT[0] = m_sOmegaRho[0]
                  //    -- this requires integration using same frequencies
                  // 5) m_sOmegaT[0] != 0 and m_sOmegaRho != 0 and
                  //    m_sOmegaT[0] != m_sOmegaRho[0]
                  //    -- this requires integration using different frequencies
                  RealVect rhoT_integralDir = RealVect::Zero;
                  for (int dir = 0; dir != SpaceDim; ++dir)
                    {
                      if ((TArgD[dir] == 0) && (rhoArgD[dir] == 0))
                        {
                          rhoT_integralDir[dir] = std::sin(m_sPhiT[dir])
                            *std::sin(m_sPhiRho[dir])*length[dir];
                        }
                      else if ((TArgD[dir] == 0) && (rhoArgD[dir] != 0))
                        {
                          rhoT_integralDir[dir] =
                            std::sin(m_sPhiT[dir])*rho_integralDir[dir];
                        }
                      else if ((TArgD[dir] != 0) && (rhoArgD[dir] == 0))
                        {
                          rhoT_integralDir[dir] =
                            T_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                        }
                      else if (TArgD[dir] == rhoArgD[dir])
                        {
                          rhoT_integralDir[dir] = (1./(4.*rhoArgD[dir]))*(
                            2.*rhoArgD[dir]*sLocHigh[dir]
                            *std::cos(m_sPhiRho[dir] - m_sPhiT[dir])
                            - std::sin(rhoArgH[dir] + TArgH[dir])
                            - 2.*rhoArgD[dir]*sLocLow[dir]
                            *std::cos(m_sPhiRho[dir] - m_sPhiT[dir])
                            + std::sin(rhoArgL[dir] + TArgL[dir]));
                        }
                      else
                        {
                          rhoT_integralDir[dir] = 0.5*(
                            std::sin(rhoArgH[dir]-TArgH[dir])/(rhoArgD[dir]
                                                              -TArgD[dir])
                          - std::sin(rhoArgH[dir]+TArgH[dir])/(rhoArgD[dir]
                                                              +TArgD[dir])
                          - std::sin(rhoArgL[dir]-TArgL[dir])/(rhoArgD[dir]
                                                              -TArgD[dir])
                          + std::sin(rhoArgL[dir]+TArgL[dir])/(rhoArgD[dir]
                                                              +TArgD[dir]));
                        }
                    }
                  Real T_4 = m_rA*(m_TA/vol)
                    *(std::sin(rhoTimeArg))*(std::sin(TTimeArg))
                    *(rhoT_integralDir.product());
                  Real rhoT_cellAvg = T_0 + T_1 + T_2 + T_4;
                  p_cellAvg = rhoT_cellAvg*R_cell;
                }

              // Now for the fun part -- cell-averaged kinetic energy
              D_TERM(
                Real uu_0 = rM*velM[0]*velM[0]; // mean rho*u*u component
                // Mean u*u*rho_sinusoid
                Real uu_1 = velM[0]*velM[0]*(m_rA/vol)
                *(std::sin(rhoTimeArg))*rho_integralDir.product();
                // Mean rho*u*u_sinusoid
                Real uu_2 = 2.*velM[0]*rM*(m_uA[0]/vol)*(std::sin(uTimeArg))
                *u_integralDir.product();
                // Mean u*rho_sinusoid*u_sinusoid
                Real uu_3 = 2.*velM[0]*m_rA*(m_uA[0]/vol)
                *(std::sin(rhoTimeArg))*(std::sin(uTimeArg))
                *(rhoU_integralDir.product());
                // Mean rho*u_sinusoid*u_sinusoid
                RealVect UU_integralDir = RealVect::Zero;
                // Note: because this is u*u there are only two options
                // 1) m_sOmegaU[0] = 0
                //    -- this requires just multiplying the two sinusoids
                // 2) m_sOmegaU[0] != 0
                //    -- this requires integration using same frequencies
                for (int dir = 0; dir != SpaceDim; ++dir)
                  {
                    if (uArgD[dir] == 0)
                      {
                        UU_integralDir[dir] = std::sin(m_sPhiU[dir])
                          *std::sin(m_sPhiU[dir])*length[dir];
                      }
                    else
                      {
                        UU_integralDir[dir] = (1./(4.*uArgD[dir]))*(
                          2.*uArgH[dir] - std::sin(2.*uArgH[dir])
                        - 2.*uArgL[dir] + std::sin(2.*uArgL[dir]));
                      }
                  }
                Real uu_4 = rM*m_uA[0]*(m_uA[0]/vol)
                *(std::sin(uTimeArg))*(std::sin(uTimeArg))
                *(UU_integralDir.product());
                // Rho_sinusoid*u_sinusoid*u_sinusoid
                RealVect rhoUU_integralDir = RealVect::Zero;
                // Note: there are 15 cases here, 10 of them are unnecessary in
                //       our case (r*u*u will always have u = u). Additionally
                //       case 5 has two variants
                //       a) m_sOmegaRho = 2.*m_sOmegaU[0]
                //       b) m_sOmegaRho != 2.*m_sOmegaU[0]
                // 1) m_sOmegaU[0] = 0 and m_sOmegaRho = 0
                //    -- this requires just multiplying the two sinusoids
                // 2) m_sOmegaU[0] = 0 and m_sOmegaRho != 0
                //    -- this requires integrating just the rho sinusoid
                // 3) m_sOmegaU[0] != 0 and m_sOmegaRho = 0
                //    -- this requires integrating just the u sinusoid
                // 4) m_sOmegaU[0] != 0 and m_sOmegaRho != 0 and
                //    m_sOmegaU[0] = m_sOmegaRho[0]
                //    -- this requires integration using same frequencies
                // 5) m_sOmegaU[0] != 0 and m_sOmegaRho != 0 and
                //    m_sOmegaU[0] != m_sOmegaRho[0]
                //    -- this requires integration using different frequencies
                for (int dir = 0; dir != SpaceDim; ++dir)
                  {
                    if ((uArgD[dir] == 0) && (rhoArgD[dir] == 0))
                      {
                        rhoUU_integralDir[dir] =
                          std::sin(m_sPhiU[dir])*std::sin(m_sPhiU[dir])
                         *std::sin(m_sPhiRho[dir])*length[dir];
                      }
                    else if ((uArgD[dir] == 0) && (rhoArgD[dir] != 0))
                      {
                        rhoUU_integralDir[dir] =
                          std::sin(m_sPhiU[dir])*std::sin(m_sPhiU[dir])
                         *rho_integralDir[dir];
                      }
                    else if ((uArgD[dir] != 0) && (rhoArgD[dir] == 0))
                      {
                        rhoUU_integralDir[dir] =
                          UU_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                      }
                    else if (uArgD[dir] == rhoArgD[dir])
                      {
                        rhoUU_integralDir[dir] = (1./(12.*uArgD[dir]))
                          *(std::cos(rhoArgH[dir] + uArgH[dir] + uArgH[dir])
                       - 3.*std::cos(rhoArgH[dir] - uArgH[dir] + uArgH[dir])
                       - 3.*std::cos(rhoArgH[dir] + uArgH[dir] - uArgH[dir])
                       - 3.*std::cos(rhoArgH[dir] - uArgH[dir] - uArgH[dir])
                          - std::cos(rhoArgL[dir] + uArgL[dir] + uArgL[dir])
                       + 3.*std::cos(rhoArgL[dir] - uArgL[dir] + uArgL[dir])
                       + 3.*std::cos(rhoArgL[dir] + uArgL[dir] - uArgL[dir])
                       + 3.*std::cos(rhoArgL[dir] - uArgL[dir] - uArgL[dir]));
                      }
                    else
                      {
                        if (rhoArgD[dir] == (2.*uArgD[dir]))
                          {
                            rhoUU_integralDir[dir] =
                              (std::cos(rhoArgH[dir] + 2.*uArgH[dir])
                               - 4.*std::cos(rhoArgH[dir])
                               - 2.*rhoArgD[dir]*sLocHigh[dir]
                               *std::sin(rhoArgH[dir] - 2.*uArgH[dir])
                               - std::cos(rhoArgL[dir] + 2.*uArgL[dir])
                               + 4.*std::cos(rhoArgL[dir])
                               + 2.*rhoArgD[dir]*sLocLow[dir]
                               *std::sin(rhoArgL[dir] - 2.*uArgL[dir])
                                )/(16.*uArgD[dir]);
                          }
                        else
                          {
                            rhoUU_integralDir[dir] =
                              (std::cos(rhoArgH[dir] - 2.*uArgH[dir])/(
                                rhoArgD[dir] - 2.*uArgD[dir])
                               + std::cos(rhoArgH[dir] + 2.*uArgH[dir])/(
                                 rhoArgD[dir] + 2.*uArgD[dir])
                               + 2.*std::sin(rhoArgD[dir]*sLocHigh[dir])
                               *std::sin(m_sPhiRho[dir])/rhoArgD[dir]
                               - 2.*std::cos(rhoArgD[dir]*sLocHigh[dir])
                               *std::cos(m_sPhiRho[dir])/rhoArgD[dir]
                               - std::cos(rhoArgL[dir] - 2.*uArgL[dir])/(
                                 rhoArgD[dir] - 2.*uArgD[dir])
                               - std::cos(rhoArgL[dir] + 2.*uArgL[dir])/(
                                 rhoArgD[dir] + 2.*uArgD[dir])
                               - 2.*std::sin(rhoArgD[dir]*sLocLow[dir])
                               *std::sin(m_sPhiRho[dir])/rhoArgD[dir]
                               + 2.*std::cos(rhoArgD[dir]*sLocLow[dir])
                               *std::cos(m_sPhiRho[dir])/rhoArgD[dir])/4;
                          }
                      }
                  }
                Real uu_5 = m_rA*m_uA[0]*(m_uA[0]/vol)
                *std::sin(uTimeArg)*std::sin(uTimeArg)*std::sin(rhoTimeArg)
                *(rhoUU_integralDir.product());
                Real rhoUU_cellAvg = uu_0 + uu_1 + uu_2 + uu_3 + uu_4 + uu_5;,

                // Rho*v*v term
                Real vv_0 = rM*velM[1]*velM[1]; // mean rho*v*v component
                // Mean v*v*rho_sinusoid
                Real vv_1 = velM[1]*velM[1]*(m_rA/vol)
                *(std::sin(rhoTimeArg))*rho_integralDir.product();
                // Mean rho*v*v_sinusoid
                Real vv_2 = 2.*velM[1]*rM*(m_uA[1]/vol)*(std::sin(vTimeArg))
                *v_integralDir.product();
                // Mean v*rho_sinusoid*v_sinusoid
                Real vv_3 = 2.*velM[1]*m_rA*(m_uA[1]/vol)
                *(std::sin(rhoTimeArg))*(std::sin(vTimeArg))
                *(rhoV_integralDir.product());
                // Mean rho*v_sinusoid*v_sinusoid
                RealVect VV_integralDir = RealVect::Zero;
                // Note: because this is v*v there are only two options
                // 1) m_sOmegaV[0] = 0
                //    -- this requires just multiplying the two sinusoids
                // 2) m_sOmegaV[0] != 0
                //    -- this requires integration using same frequencies
                for (int dir = 0; dir != SpaceDim; ++dir)
                  {
                    if (vArgD[dir] == 0)
                      {
                        VV_integralDir[dir] = std::sin(m_sPhiV[dir])
                          *std::sin(m_sPhiV[dir])*length[dir];
                      }
                    else
                      {
                        VV_integralDir[dir] = (1./(4.*vArgD[dir]))*(
                          2.*vArgH[dir] - std::sin(2.*vArgH[dir])
                        - 2.*vArgL[dir] + std::sin(2.*vArgL[dir]));
                      }
                  }
                Real vv_4 = rM*m_uA[1]*(m_uA[1]/vol)
                *(std::sin(vTimeArg))*(std::sin(vTimeArg))
                *(VV_integralDir.product());
                // Rho_sinusoid*v_sinusoid*v_sinusoid
                RealVect rhoVV_integralDir = RealVect::Zero;
                // Note: there are 15 cases here, 10 of them are unnecessary in
                //       our case (r*v*v will always have v = v). Additionally
                //       case 5 has two variants
                //       a) m_sOmegaRho = 2.*m_sOmegaV[0]
                //       b) m_sOmegaRho != 2.*m_sOmegaV[0]
                // 1) m_sOmegaV[0] = 0 and m_sOmegaRho = 0
                //    -- this requires just multiplying the two sinusoids
                // 2) m_sOmegaV[0] = 0 and m_sOmegaRho != 0
                //    -- this requires integrating just the rho sinusoid
                // 3) m_sOmegaV[0] != 0 and m_sOmegaRho = 0
                //    -- this requires integrating just the u sinusoid
                // 4) m_sOmegaV[0] != 0 and m_sOmegaRho != 0 and
                //    m_sOmegaV[0] = m_sOmegaRho[0]
                //    -- this requires integration using same frequencies
                // 5) m_sOmegaV[0] != 0 and m_sOmegaRho != 0 and
                //    m_sOmegaV[0] != m_sOmegaRho[0]
                //    -- this requires integration using different frequencies
                for (int dir = 0; dir != SpaceDim; ++dir)
                  {
                    if ((vArgD[dir] == 0) && (rhoArgD[dir] == 0))
                      {
                        rhoVV_integralDir[dir] =
                          std::sin(m_sPhiV[dir])*std::sin(m_sPhiV[dir])
                         *std::sin(m_sPhiRho[dir])*length[dir];
                      }
                    else if ((vArgD[dir] == 0) && (rhoArgD[dir] != 0))
                      {
                        rhoVV_integralDir[dir] =
                          std::sin(m_sPhiV[dir])*std::sin(m_sPhiV[dir])
                         *rho_integralDir[dir];
                      }
                    else if ((vArgD[dir] != 0) && (rhoArgD[dir] == 0))
                      {
                        rhoVV_integralDir[dir] =
                          VV_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                      }
                    else if (vArgD[dir] == rhoArgD[dir])
                      {
                        rhoVV_integralDir[dir] = (1./(12.*vArgD[dir]))
                          *(std::cos(rhoArgH[dir] + vArgH[dir] + vArgH[dir])
                       - 3.*std::cos(rhoArgH[dir] - vArgH[dir] + vArgH[dir])
                       - 3.*std::cos(rhoArgH[dir] + vArgH[dir] - vArgH[dir])
                       - 3.*std::cos(rhoArgH[dir] - vArgH[dir] - vArgH[dir])
                          - std::cos(rhoArgL[dir] + vArgL[dir] + vArgL[dir])
                       + 3.*std::cos(rhoArgL[dir] - vArgL[dir] + vArgL[dir])
                       + 3.*std::cos(rhoArgL[dir] + vArgL[dir] - vArgL[dir])
                       + 3.*std::cos(rhoArgL[dir] - vArgL[dir] - vArgL[dir]));
                      }
                    else
                      {
                        if (rhoArgD[dir] == (2.*vArgD[dir]))
                          {
                            rhoVV_integralDir[dir] =
                              (std::cos(rhoArgH[dir] + 2.*vArgH[dir])
                               - 4.*std::cos(rhoArgH[dir])
                               - 2.*rhoArgD[dir]*sLocHigh[dir]
                               *std::sin(rhoArgH[dir] - 2.*vArgH[dir])
                               - std::cos(rhoArgL[dir] + 2.*vArgL[dir])
                               + 4.*std::cos(rhoArgL[dir])
                               + 2.*rhoArgD[dir]*sLocLow[dir]
                               *std::sin(rhoArgL[dir] - 2.*vArgL[dir])
                                )/(16.*vArgD[dir]);
                          }
                        else
                          {
                            rhoVV_integralDir[dir] =
                              (std::cos(rhoArgH[dir] - 2.*vArgH[dir])/(
                                rhoArgD[dir] - 2.*vArgD[dir])
                               + std::cos(rhoArgH[dir] + 2.*vArgH[dir])/(
                                 rhoArgD[dir] + 2.*vArgD[dir])
                               + 2.*std::sin(rhoArgD[dir]*sLocHigh[dir])
                               *std::sin(m_sPhiRho[dir])/rhoArgD[dir]
                               - 2.*std::cos(rhoArgD[dir]*sLocHigh[dir])
                               *std::cos(m_sPhiRho[dir])/rhoArgD[dir]
                               - std::cos(rhoArgL[dir] - 2.*vArgL[dir])/(
                                 rhoArgD[dir] - 2.*vArgD[dir])
                               - std::cos(rhoArgL[dir] + 2.*vArgL[dir])/(
                                 rhoArgD[dir] + 2.*vArgD[dir])
                               - 2.*std::sin(rhoArgD[dir]*sLocLow[dir])
                               *std::sin(m_sPhiRho[dir])/rhoArgD[dir]
                               + 2.*std::cos(rhoArgD[dir]*sLocLow[dir])
                               *std::cos(m_sPhiRho[dir])/rhoArgD[dir])/4;
                          }
                      }
                  }
                Real vv_5 = m_rA*m_uA[1]*(m_uA[1]/vol)
                *std::sin(vTimeArg)*std::sin(vTimeArg)*std::sin(rhoTimeArg)
                *(rhoVV_integralDir.product());
                Real rhoVV_cellAvg = vv_0 + vv_1 + vv_2 + vv_3 + vv_4 + vv_5;,

                // Rho*w*w terms
                Real ww_0 = rM*velM[2]*velM[2]; // mean rho*w*w component
                // Mean w*w*rho_sinusoid
                Real ww_1 = velM[2]*velM[2]*(m_rA/vol)
                *(std::sin(rhoTimeArg))*rho_integralDir.product();
                // Mean rho*w*w_sinusoid
                Real ww_2 = 2.*velM[2]*rM*(m_uA[2]/vol)*(std::sin(wTimeArg))
                *w_integralDir.product();
                // Mean w*rho_sinusoid*w_sinusoid
                Real ww_3 = 2.*velM[2]*m_rA*(m_uA[2]/vol)
                *(std::sin(rhoTimeArg))*(std::sin(wTimeArg))
                *(rhoW_integralDir.product());
                // Mean rho*w_sinusoid*w_sinusoid
                RealVect WW_integralDir = RealVect::Zero;
                // Note: because this is w*w there are only two options
                // 1) m_sOmegaW[0] = 0
                //    -- this requires just multiplying the two sinusoids
                // 2) m_sOmegaW[0] != 0
                //    -- this requires integration using same frequencies
                for (int dir = 0; dir != SpaceDim; ++dir)
                  {
                    if (wArgD[dir] == 0)
                      {
                        WW_integralDir[dir] = std::sin(m_sPhiW[dir])
                          *std::sin(m_sPhiW[dir])*length[dir];
                      }
                    else
                      {
                        WW_integralDir[dir] = (1./(4.*wArgD[dir]))*(
                          2.*wArgH[dir] - std::sin(2.*wArgH[dir])
                        - 2.*wArgL[dir] + std::sin(2.*wArgL[dir]));
                      }
                  }
                Real ww_4 = rM*m_uA[2]*(m_uA[2]/vol)
                *(std::sin(wTimeArg))*(std::sin(wTimeArg))
                *(WW_integralDir.product());
                // Rho_sinusoid*w_sinusoid*w_sinusoid
                RealVect rhoWW_integralDir = RealVect::Zero;
                // Note: there are 15 cases here, 10 of them are unnecessary in
                //       our case (r*w*w will always have w = w). Additionally
                //       case 5 has two variants
                //       a) m_sOmegaRho = 2.*m_sOmegaW[0]
                //       b) m_sOmegaRho != 2.*m_sOmegaW[0]
                // 1) m_sOmegaW[0] = 0 and m_sOmegaRho = 0
                //    -- this requires just multiplying the two sinusoids
                // 2) m_sOmegaW[0] = 0 and m_sOmegaRho != 0
                //    -- this requires integrating just the rho sinusoid
                // 3) m_sOmegaW[0] != 0 and m_sOmegaRho = 0
                //    -- this requires integrating just the u sinusoid
                // 4) m_sOmegaW[0] != 0 and m_sOmegaRho != 0 and
                //    m_sOmegaW[0] = m_sOmegaRho[0]
                //    -- this requires integration using same frequencies
                // 5) m_sOmegaW[0] != 0 and m_sOmegaRho != 0 and
                //    m_sOmegaW[0] != m_sOmegaRho[0]
                //    -- this requires integration using different frequencies
                for (int dir = 0; dir != SpaceDim; ++dir)
                  {
                    if ((wArgD[dir] == 0) && (rhoArgD[dir] == 0))
                      {
                        rhoWW_integralDir[dir] =
                          std::sin(m_sPhiW[dir])*std::sin(m_sPhiW[dir])
                         *std::sin(m_sPhiRho[dir])*length[dir];
                      }
                    else if ((wArgD[dir] == 0) && (rhoArgD[dir] != 0))
                      {
                        rhoWW_integralDir[dir] =
                          std::sin(m_sPhiW[dir])*std::sin(m_sPhiW[dir])
                         *rho_integralDir[dir];
                      }
                    else if ((wArgD[dir] != 0) && (rhoArgD[dir] == 0))
                      {
                        rhoWW_integralDir[dir] =
                          WW_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                      }
                    else if (wArgD[dir] == rhoArgD[dir])
                      {
                        rhoWW_integralDir[dir] = (1./(12.*wArgD[dir]))
                          *(std::cos(rhoArgH[dir] + wArgH[dir] + wArgH[dir])
                       - 3.*std::cos(rhoArgH[dir] - wArgH[dir] + wArgH[dir])
                       - 3.*std::cos(rhoArgH[dir] + wArgH[dir] - wArgH[dir])
                       - 3.*std::cos(rhoArgH[dir] - wArgH[dir] - wArgH[dir])
                          - std::cos(rhoArgL[dir] + wArgL[dir] + wArgL[dir])
                       + 3.*std::cos(rhoArgL[dir] - wArgL[dir] + wArgL[dir])
                       + 3.*std::cos(rhoArgL[dir] + wArgL[dir] - wArgL[dir])
                       + 3.*std::cos(rhoArgL[dir] - wArgL[dir] - wArgL[dir]));
                      }
                    else
                      {
                        if (rhoArgD[dir] == (2.*wArgD[dir]))
                          {
                            rhoWW_integralDir[dir] =
                              (std::cos(rhoArgH[dir] + 2.*wArgH[dir])
                               - 4.*std::cos(rhoArgH[dir])
                               - 2.*rhoArgD[dir]*sLocHigh[dir]
                               *std::sin(rhoArgH[dir] - 2.*wArgH[dir])
                               - std::cos(rhoArgL[dir] + 2.*wArgL[dir])
                               + 4.*std::cos(rhoArgL[dir])
                               + 2.*rhoArgD[dir]*sLocLow[dir]
                               *std::sin(rhoArgL[dir] - 2.*wArgL[dir])
                                )/(16.*wArgD[dir]);
                          }
                        else
                          {
                            rhoWW_integralDir[dir] =
                              (std::cos(rhoArgH[dir] - 2.*wArgH[dir])/(
                                rhoArgD[dir] - 2.*wArgD[dir])
                               + std::cos(rhoArgH[dir] + 2.*wArgH[dir])/(
                                 rhoArgD[dir] + 2.*wArgD[dir])
                               + 2.*std::sin(rhoArgD[dir]*sLocHigh[dir])
                               *std::sin(m_sPhiRho[dir])/rhoArgD[dir]
                               - 2.*std::cos(rhoArgD[dir]*sLocHigh[dir])
                               *std::cos(m_sPhiRho[dir])/rhoArgD[dir]
                               - std::cos(rhoArgL[dir] - 2.*wArgL[dir])/(
                                 rhoArgD[dir] - 2.*wArgD[dir])
                               - std::cos(rhoArgL[dir] + 2.*wArgL[dir])/(
                                 rhoArgD[dir] + 2.*wArgD[dir])
                               - 2.*std::sin(rhoArgD[dir]*sLocLow[dir])
                               *std::sin(m_sPhiRho[dir])/rhoArgD[dir]
                               + 2.*std::cos(rhoArgD[dir]*sLocLow[dir])
                               *std::cos(m_sPhiRho[dir])/rhoArgD[dir])/4;
                          }
                      }
                  }
                Real ww_5 = m_rA*m_uA[2]*(m_uA[2]/vol)
                *std::sin(wTimeArg)*std::sin(wTimeArg)*std::sin(rhoTimeArg)
                *(rhoWW_integralDir.product());
                Real rhoWW_cellAvg = ww_0 + ww_1 + ww_2 + ww_3 + ww_4 + ww_5;
                );
              Real ke = 0.5*(D_TERM(rhoUU_cellAvg,
                                  + rhoVV_cellAvg,
                                  + rhoWW_cellAvg));
              Real energy = p_cellAvg/(gamma - 1.) + ke;
              UFab[MD_IX(i, presIndx)] = energy;
            }
        }
    }
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
 *  \param[in]  a_didx  Index of current box in layout
 *  \param[in]  a_globalKE
 *                      Global mean kinetic energy (turbulent forcing)
 *  \param[in]  a_globalHelicity
 *                      Global mean helicity (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void
CNSIBCMMS::addSourceTerm(FArrayBox&           a_sourceFab,
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
  // NOTE: we need some switches here
  // 1) ability to test multispecies
  // 2) ability to test additional models (e.g. SV SGS LES model)

  // NOTE: this all assumes calorically perfect physics
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  // Get physical coordinates
  FABSTACKTEMP(XiFab, a_solveBox, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, a_solveBox, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(a_solveBox, XiFab, XFab, blockCoordSys);
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const Real gamma  = CRDparam::g_gamma;
  const Real mu     = CRDparam::g_mu;
  const Real lambda = CRDparam::g_lambda;
  const Real R      = CRDparam::g_R;
  const Real kappa  = CRDparam::g_K;
  // Obtain the mean values and the wave-amplitude values
  const CRDState& mean = CRDState::get(m_idxStateInit);
  const Real rM = mean.density();
  const RealVect velM = mean.velocity();
  const Real pM = mean.pressure();
  const Real TM = mean.temperature();

  MD_BOXLOOP(a_solveBox, i)
    {
      // Physical location
      const RealVect loc(
        D_DECL(XFab[MD_IX(i, 0)], XFab[MD_IX(i, 1)], XFab[MD_IX(i, 2)]));
      // Spatial location normalized by characteristic length
      const RealVect sLoc = loc/m_charLength;
      // Time normalized by characteristic time
      const Real tLoc = a_time/m_charTime;

      // Spatial argument for density field
      const RealVect rA = 2.*PI*m_sOmegaRho*sLoc + m_sPhiRho;
      // Derivative with respect to space of spatial argument for density field
      const RealVect rB = 2.*PI*m_sOmegaRho/m_charLength;
      // Temporal argument for density field
      const Real rTA = 2.*PI*m_tOmegaRho*tLoc + m_tPhiRho;
      // Derivative with respect to time of temporal argument for density field
      const Real rTB = 2.*PI*m_tOmegaRho/m_charTime;

      // Spatial + temporal arguments for velocity fields
      D_TERM(const RealVect uA = 2.*PI*m_sOmegaU*sLoc + m_sPhiU;
             // Temporal argument for x-velocity field
             const Real uTA = 2.*PI*m_tOmegaU*tLoc + m_tPhiU;
             // Derivative with respect to space of uA
             const RealVect uB = 2.*PI*m_sOmegaU/m_charLength;
             // Derivative with respect to time of uTA
             const Real uTB = 2.*PI*m_tOmegaU/m_charTime;,
             // Spatial argument for y-velocity field
             const RealVect vA = 2.*PI*m_sOmegaV*sLoc + m_sPhiV;
             // Temporal argument for y-velocity field
             const Real vTA = 2.*PI*m_tOmegaV*tLoc + m_tPhiV;
             // Derivative with respect to space of vA
             const RealVect vB = 2.*PI*m_sOmegaV/m_charLength;
             // Derivative with respect to time of vTA
             const Real vTB = 2.*PI*m_tOmegaV/m_charTime;,
             // Spatial argument for z-velocity field
             const RealVect wA = 2.*PI*m_sOmegaW*sLoc + m_sPhiW;
             // Temporal argument for z-velocity field
             const Real wTA = 2.*PI*m_tOmegaW*tLoc + m_tPhiW;
             // Derivative with respect to space of wA
             const RealVect wB = 2.*PI*m_sOmegaW/m_charLength;
             // Derivative with respect to time of wTA
             const Real wTB = 2.*PI*m_tOmegaW/m_charTime;);

      // Spatial argument for pressure field
      const RealVect pA = 2.*PI*m_sOmegaP*sLoc + m_sPhiP;
      // Derivative with respect to space of spatial argument for pressure
      const RealVect pB = 2.*PI*m_sOmegaP/m_charLength;
      // Temporal argument for pressure field
      const Real pTA = 2.*PI*m_tOmegaP*tLoc + m_tPhiP;
      // Derivative with respect to time of temporal argument for pressure
      const Real pTB = 2.*PI*m_tOmegaP/m_charTime;

      // Spatial argument for temperature field
      const RealVect TA = 2.*PI*m_sOmegaT*sLoc + m_sPhiT;
      // Derivative with respect to space of spatial argument for temperature
      const RealVect TB = 2.*PI*m_sOmegaT/m_charLength;
      // Temporal argument for temperature field
      const Real TTA = 2.*PI*m_tOmegaT*tLoc + m_tPhiT;
      // Derivative with respect to time of temporal argument for temperature
      const Real TTB = 2.*PI*m_tOmegaT/m_charTime;

      // Spatial sine terms for density field
      const RealVect rS(
        D_DECL(std::sin(rA[0]),std::sin(rA[1]),std::sin(rA[2])));
      // Temporal sine term for density field
      const Real rTS = std::sin(rTA);
      // Spatial cosine terms for density field
      const RealVect rC(
        D_DECL(std::cos(rA[0]),std::cos(rA[1]),std::cos(rA[2])));
      // Temporal cosine term for density field
      const Real rTC = std::cos(rTA);

      // Density sine-series
      const Real Rho = rM + m_rA*rTS*D_TERM(rS[0],*rS[1],*rS[2]);
      // Spatial derivatives of density sine-series
      D_TERM(
        const Real rAdX = m_rA*rTS*D_TERM(rB[0]*rC[0],*rS[1],*rS[2]);,
        const Real rAdY = m_rA*rTS*D_TERM(rS[0],*rB[1]*rC[1],*rS[2]);,
        const Real rAdZ = m_rA*rTS*D_TERM(rS[0],*rS[1],*rB[2]*rC[2]););
      // Temporal derivative of density sine-series
      const Real rAdT = m_rA*rTC*rTB*D_TERM(rS[0],*rS[1],*rS[2]);
      // Second-derivatives of density sine-series
      D_TERM(
        const Real rAdXX =
        -m_rA*rTS*D_TERM(rB[0]*rB[0]*rS[0],*rS[1],*rS[2]);,
        const Real rAdYY =
        -m_rA*rTS*D_TERM(rS[0],*rB[1]*rB[1]*rS[1],*rS[2]);,
        const Real rAdZZ =
        -m_rA*rTS*D_TERM(rS[0],*rS[1],*rB[2]*rB[2]*rS[2]););

      // Spatial gradient of density
      const RealVect GradRho(D_DECL(rAdX,rAdY,rAdZ));
      // Diagonal terms of the density Hessian
      const RealVect SecondDerivRho(D_DECL(rAdXX, rAdYY, rAdZZ));
      // Square of the gradient of density
      const RealVect GradRhoSqrd(D_DECL(rAdX*rAdX,rAdY*rAdY,rAdZ*rAdZ));

      D_TERM(
        // Spatial sine terms for the x-velocity field
        const RealVect uS(
          D_DECL(std::sin(uA[0]),std::sin(uA[1]),std::sin(uA[2])));
        // Temporal sine term for x-velocity field
        const Real uTS = std::sin(uTA);
        // Spatial cosine term for x-velocity field
        const RealVect uC(
          D_DECL(std::cos(uA[0]),std::cos(uA[1]),std::cos(uA[2])));
        // Temporal cosine term for x-velocity field
        const Real uTC = std::cos(uTA);

        // X-velocity
        const Real UVel = velM[0] + m_uA[0]*uTS*D_TERM(uS[0],*uS[1],*uS[2]);
        // Spatial derivatives of x-velocity sine-series
        D_TERM(
          const Real uAdX = m_uA[0]*uTS*D_TERM(uB[0]*uC[0],*uS[1],*uS[2]);,
          const Real uAdY = m_uA[0]*uTS*D_TERM(uS[0],*uB[1]*uC[1],*uS[2]);,
          const Real uAdZ = m_uA[0]*uTS*D_TERM(uS[0],*uS[1],*uB[2]*uC[2]););
        // Temporal derivative of x-velocity sine-series
        const Real uAdT = m_uA[0]*uTB*uTC*D_TERM(uS[0],*uS[1],*uS[2]);
        // Spatial second-derivatives of x-velocity sine-series
        D_TERM(
          const Real uAdXX =
            -m_uA[0]*uTS*D_TERM(uB[0]*uB[0]*uS[0],*uS[1],*uS[2]);,
          const Real uAdXY =
            m_uA[0]*uTS*D_TERM(uB[0]*uC[0],*uB[1]*uC[1],*uS[2]);
          const Real uAdYY =
            -m_uA[0]*uTS*D_TERM(uS[0],*uB[1]*uB[1]*uS[1],*uS[2]);,
          const Real uAdXZ =
            m_uA[0]*uTS*D_TERM(uB[0]*uC[0],*uS[1],*uB[2]*uC[2]);
          const Real uAdZZ =
            -m_uA[0]*uTS*D_TERM(uS[0],*uS[1],*uB[2]*uB[2]*uS[2]););
        // Gradient of x-velocity
        const RealVect GradUVel(D_DECL(uAdX,uAdY,uAdZ));
        // Gradient of x-velocity squared
        const RealVect GradUVelSqrd(D_DECL(uAdX*uAdX,uAdY*uAdY,uAdZ*uAdZ));
        // Gradient of x-derivative of x-velocity
        const RealVect DXGradUVel(D_DECL(uAdXX,uAdXY,uAdXZ));
        // Laplacian of x-velocity
        const Real LapUVel = D_TERM(uAdXX, + uAdYY, + uAdZZ);,

        // Spatial sine terms for the y-velocity field
        const RealVect vS(
          D_DECL(std::sin(vA[0]),std::sin(vA[1]),std::sin(vA[2])));
        // Temporal sine term for y-velocity field
        const Real vTS = std::sin(vTA);
        // Spatial cosine term for y-velocity field
        const RealVect vC(
          D_DECL(std::cos(vA[0]),std::cos(vA[1]),std::cos(vA[2])));
        // Temporal cosine term for y-velocity field
        const Real vTC = std::cos(vTA);

        // Y-velocity
        const Real VVel = velM[1] + m_uA[1]*vTS*D_TERM(vS[0],*vS[1],*vS[2]);
        // Spatial derivatives of y-velocity sine-series
        D_TERM(
          const Real vAdX = m_uA[1]*vTS*D_TERM(vB[0]*vC[0],*vS[1],*vS[2]);,
          const Real vAdY = m_uA[1]*vTS*D_TERM(vS[0],*vB[1]*vC[1],*vS[2]);,
          const Real vAdZ = m_uA[1]*vTS*D_TERM(vS[0],*vS[1],*vB[2]*vC[2]););
        // Temporal derivative of y-velocity sine-series
        const Real vAdT = m_uA[1]*vTB*vTC*D_TERM(vS[0],*vS[1],*vS[2]);
        // Spatial second-derivatives of y-velocity sine-series
        D_TERM(
          const Real vAdXX =
            -m_uA[1]*vTS*D_TERM(vB[0]*vB[0]*vS[0],*vS[1],*vS[2]);,
          const Real vAdXY =
            m_uA[1]*vTS*D_TERM(vB[0]*vC[0],*vB[1]*vC[1],*vS[2]);
          const Real vAdYY =
            -m_uA[1]*vTS*D_TERM(vS[0],*vB[1]*vB[1]*vS[1],*vS[2]);,
          const Real vAdYZ =
            m_uA[1]*vTS*D_TERM(vS[0],*vB[1]*vC[1],*vB[2]*vC[2]);
          const Real vAdZZ =
            -m_uA[1]*vTS*D_TERM(vS[0],*vS[1],*vB[2]*vB[2]*vS[2]););
        // Gradient of y-velocity
        const RealVect GradVVel(D_DECL(vAdX,vAdY,vAdZ));
        // Gradient of y-velocity squared
        const RealVect GradVVelSqrd(D_DECL(vAdX*vAdX,vAdY*vAdY,vAdZ*vAdZ));
        // Gradient of y-derivative of y-velocity
        const RealVect DYGradVVel(D_DECL(vAdXY,vAdYY,vAdYZ));
        // Laplacian of y-velocity
        const Real LapVVel = D_TERM(vAdXX, + vAdYY, + vAdZZ);,

        // Spatial sine terms for the z-velocity field
        const RealVect wS(
          D_DECL(std::sin(wA[0]),std::sin(wA[1]),std::sin(wA[2])));
        // Temporal sine term for z-velocity field
        const Real wTS = std::sin(wTA);
        // Spatial cosine term for z-velocity field
        const RealVect wC(
          D_DECL(std::cos(wA[0]),std::cos(wA[1]),std::cos(wA[2])));
        // Temporal cosine term for z-velocity field
        const Real wTC = std::cos(wTA);

        // Z-velocity
        const Real WVel = velM[2] + m_uA[2]*wTS*D_TERM(wS[0],*wS[1],*wS[2]);
        // Spatial derivatives of z-velocity sine-series
        D_TERM(
          const Real wAdX = m_uA[2]*wTS*D_TERM(wB[0]*wC[0],*wS[1],*wS[2]);,
          const Real wAdY = m_uA[2]*wTS*D_TERM(wS[0],*wB[1]*wC[1],*wS[2]);,
          const Real wAdZ = m_uA[2]*wTS*D_TERM(wS[0],*wS[1],*wB[2]*wC[2]););
        // Temporal derivative of z-velocity sine-series
        const Real wAdT = m_uA[2]*wTB*wTC*D_TERM(wS[0],*wS[1],*wS[2]);
        // Spatial second-derivatives of z-velocity sine-series
        D_TERM(
          const Real wAdXX =
            -m_uA[2]*wTS*D_TERM(wB[0]*wB[0]*wS[0],*wS[1],*wS[2]);,
          const Real wAdYY =
            -m_uA[2]*wTS*D_TERM(wS[0],*wB[1]*wB[1]*wS[1],*wS[2]);,
          const Real wAdXZ =
            m_uA[2]*wTS*D_TERM(wB[0]*wC[0],*wS[1],*wB[2]*wC[2]);
          const Real wAdYZ =
            m_uA[2]*wTS*D_TERM(wS[0],*wB[1]*wC[1],*wB[2]*wC[2]);
          const Real wAdZZ =
            -m_uA[2]*wTS*D_TERM(wS[0],*wS[1],*wB[2]*wB[2]*wS[2]););
        // Gradient of z-velocity
        const RealVect GradWVel(D_DECL(wAdX,wAdY,wAdZ));
        // Gradient of z-velocity squared
        const RealVect GradWVelSqrd(D_DECL(wAdX*wAdX,wAdY*wAdY,wAdZ*wAdZ));
        // Gradient of z-derivative of z-velocity
        const RealVect DZGradWVel(D_DECL(wAdXZ,wAdYZ,wAdZZ));
        // Laplacian of z-velocity
        const Real LapWVel = D_TERM(wAdXX, + wAdYY, + wAdZZ););

      // Thermodynamic variables must be computed based on BCs
      Real Press = 0.;
      Real T = 0.;
      RealVect GradP(RealVect::Zero);
      RealVect GradT(RealVect::Zero);
      Real LapT = 0.; // Laplacian of temperature
      Real pAdT = 0.;
      Real TAdT = 0.;

      if (!m_bcTestCase)
        {
          // Spatial sine terms for pressure field
          const RealVect pS(
            D_DECL(std::sin(pA[0]),std::sin(pA[1]),std::sin(pA[2])));
          // Temporal sine term for pressure field
          const Real pTS = std::sin(pTA);
          // Spatial cosine terms for pressure field
          const RealVect pC(
            D_DECL(std::cos(pA[0]),std::cos(pA[1]),std::cos(pA[2])));
          // Temporal cosine term for pressure field
          const Real pTC = std::cos(pTA);

          // Pressure sine-series
          Press = pM + m_pA*pTS*D_TERM(pS[0],*pS[1],*pS[2]);
          // Spatial derivatives of pressure sine-series
          D_TERM(
            const Real pAdX = m_pA*pTS*D_TERM(pB[0]*pC[0],*pS[1],*pS[2]);,
            const Real pAdY = m_pA*pTS*D_TERM(pS[0],*pB[1]*pC[1],*pS[2]);,
            const Real pAdZ = m_pA*pTS*D_TERM(pS[0],*pS[1],*pB[2]*pC[2]););
          // Temporal derivative of pressure sine-series
          pAdT = m_pA*pTC*pTB*D_TERM(pS[0],*pS[1],*pS[2]);
          // Second-derivatives of pressure sine-series
          D_TERM(
            const Real pAdXX =
            -m_pA*pTS*D_TERM(pB[0]*pB[0]*pS[0],*pS[1],*pS[2]);,
            const Real pAdYY =
            -m_pA*pTS*D_TERM(pS[0],*pB[1]*pB[1]*pS[1],*pS[2]);,
            const Real pAdZZ =
            -m_pA*pTS*D_TERM(pS[0],*pS[1],*pB[2]*pB[2]*pS[2]););
          // Spatial gradient of pressure
          GradP = RealVect(D_DECL(pAdX,pAdY,pAdZ));
          // Diagonal terms of the density Hessian
          const RealVect SecondDerivP(D_DECL(pAdXX, pAdYY, pAdZZ));

          // Pressure gradient times density gradient
          const RealVect GradPTimesGradRho(
            D_DECL(pAdX*rAdX, pAdY*rAdY, pAdZ*rAdZ));
          // Inverse density
          const Real invRho = 1./Rho;
          // Inverse R
          const Real invR = 1./R;
          // Temperature
          T = Press/(Rho*R);
          // Temperature Laplacian
          RealVect LapTVect =
            invRho*(invR*(SecondDerivP - 2.*invRho*GradPTimesGradRho)
                    + T*(2.*invRho*GradRhoSqrd - SecondDerivRho));
          LapT = LapTVect.sum();
        }
      else
        {
          // Spatial sine terms for temperature field
          const RealVect TS(
            D_DECL(std::sin(TA[0]),std::sin(TA[1]),std::sin(TA[2])));
          // Temporal sine term for temperature field
          const Real TTS = std::sin(TTA);
          // Spatial cosine terms for temperature field
          const RealVect TC(
            D_DECL(std::cos(TA[0]),std::cos(TA[1]),std::cos(TA[2])));
          // Temporal cosine term for temperature field
          const Real TTC = std::cos(TTA);

          // Temperature sine-series
          T = TM + m_TA*TTS*D_TERM(TS[0],*TS[1],*TS[2]);
          // Spatial derivatives of temperature sine-series
          D_TERM(
            const Real TAdX = m_TA*TTS*D_TERM(TB[0]*TC[0],*TS[1],*TS[2]);,
            const Real TAdY = m_TA*TTS*D_TERM(TS[0],*TB[1]*TC[1],*TS[2]);,
            const Real TAdZ = m_TA*TTS*D_TERM(TS[0],*TS[1],*TB[2]*TC[2]););
          // Spatial gradient of temperature
          GradT = RealVect(D_DECL(TAdX,TAdY,TAdZ));
          // Second-derivatives of temperature sine-series
          D_TERM(
            const Real TAdXX =
            -m_TA*TTS*D_TERM(TB[0]*TB[0]*TS[0],*TS[1],*TS[2]);,
            const Real TAdYY =
            -m_TA*TTS*D_TERM(TS[0],*TB[1]*TB[1]*TS[1],*TS[2]);,
            const Real TAdZZ =
            -m_TA*TTS*D_TERM(TS[0],*TS[1],*TB[2]*TB[2]*TS[2]););
          // Laplacian of temperature
          LapT = D_TERM(TAdXX, + TAdYY, + TAdZZ);
          // Temporal derivative of temperature sine-series
          TAdT = m_TA*TTC*TTB*D_TERM(TS[0],*TS[1],*TS[2]);

          // Pressure
          Press = Rho*R*T;
          // Pressure gradient
          GradP = GradRho*R*T + Rho*R*GradT;
          pAdT = rAdT*R*T + Rho*R*TAdT;
        }

      // Velocity vector
      const RealVect Vel(D_DECL(UVel, VVel, WVel));
      // Divergence of the velocity vector
      const Real DivVel = D_TERM(GradUVel[0], + GradVVel[1], + GradWVel[2]);
      // Gradient of velocity multiplied by velocity
      const RealVect VelGradTimesVel(D_DECL(GradUVel.dotProduct(Vel),
                                            GradVVel.dotProduct(Vel),
                                            GradWVel.dotProduct(Vel)));
      // Kinetic Energy
      const Real KE = D_TERM(UVel*UVel, + VVel*VVel, + WVel*WVel);
      // Time derivative of velocity
      const RealVect VelDT(D_DECL(uAdT, vAdT, wAdT));
      // Gradient of the kinetic energy
      const RealVect GradKE(D_DECL(D_TERM(
        GradUVel[0]*Vel[0], + GradVVel[0]*Vel[1], + GradWVel[0]*Vel[2]),D_TERM(
        GradUVel[1]*Vel[0], + GradVVel[1]*Vel[1], + GradWVel[1]*Vel[2]),D_TERM(
        GradUVel[2]*Vel[0], + GradVVel[2]*Vel[1], + GradWVel[2]*Vel[2])));
      // Standard form of the kinetic energy (0.5*rho*(u\dot u))
      const Real HalfRhoKE = 0.5*Rho*KE;

      // Finally, the source terms that we need for the inviscid terms
      Real rhoInertialS = 0.;
      RealVect rhoUInertialS(RealVect::Zero);
      Real rhoEInertialS = 0.;
      if (CRDparam::g_physicsModels & CRDparam::PhysicsInertial)
        {
          rhoInertialS = GradRho.dotProduct(Vel) + Rho*DivVel;

          rhoUInertialS = GradP + Rho*DivVel*Vel
            + Vel*GradRho.dotProduct(Vel) + Rho*VelGradTimesVel;

          rhoEInertialS =
            (GradP.dotProduct(Vel))/(gamma - 1.) + Rho*Vel.dotProduct(GradKE)
            + 0.5*KE*GradRho.dotProduct(Vel) + GradP.dotProduct(Vel)
            + Press*DivVel + DivVel*(HalfRhoKE + Press/(gamma - 1.));
        }

      // Source terms that we need for the viscous terms
      RealVect rhoUViscS(RealVect::Zero);
      Real rhoEViscousS = 0.;
      if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
        {
          // Laplacian of the velocity vector
          const RealVect LapVelVec(D_DECL(LapUVel, LapVVel, LapWVel));
          // Gradient of the divergence of velocity
          const RealVect GradDivVel =
            D_TERM(DXGradUVel,+DYGradVVel,+DZGradWVel);
          // Sum of the squares of all velocity derivatives
          const Real GradVelTensorProd =
            D_TERM(GradUVelSqrd.sum(),+GradVVelSqrd.sum(),+GradWVelSqrd.sum());
          // Sum of all comps of (Grad(Velocity))*(Grad(Velocity))^transpose
          // In other notation, Grad(Velocity)_{ij}*Grad(Velocity)_{ji}
          const Real GradVelTransposeTensorProd = D_TERM(
            GradUVelSqrd[0], + GradVVelSqrd[1] + 2.*GradUVel[1]*GradVVel[0],
          + GradWVelSqrd[2] + 2.*GradUVel[2]*GradWVel[0]
          + 2.*GradVVel[2]*GradWVel[1]);

          const Real lFact = (1. - 2.*lambda);
          rhoUViscS = mu*(LapVelVec + lFact*GradDivVel);
          rhoEViscousS =
            mu*(LapVelVec.dotProduct(Vel) + lFact*GradDivVel.dotProduct(Vel)
                + GradVelTensorProd + GradVelTransposeTensorProd
                - 2.*lambda*DivVel*DivVel) + kappa*LapT;
        }

      // Source terms that we need for temporal derivative terms
      Real rhoTimeS = rAdT;
      RealVect rhoUTimeS = Rho*VelDT + rAdT*Vel;
      Real rhoETimeS =
        pAdT/(gamma - 1.) + 0.5*rAdT*KE + Rho*Vel.dotProduct(VelDT);

      // Final, summed source terms for each conservation equation
      a_sourceFab[MD_IX(i, rhoIndx)]      = rhoTimeS + rhoInertialS;
      D_TERM(
        a_sourceFab[MD_IX(i, WvelIndx)]   =
          rhoUTimeS[0] + rhoUInertialS[0] - rhoUViscS[0];,
        a_sourceFab[MD_IX(i, WvelIndx+1)] =
          rhoUTimeS[1] + rhoUInertialS[1] - rhoUViscS[1];,
        a_sourceFab[MD_IX(i, WvelIndx+2)] =
          rhoUTimeS[2] + rhoUInertialS[2] - rhoUViscS[2];);
      a_sourceFab[MD_IX(i, presIndx)]     =
        rhoETimeS + rhoEInertialS - rhoEViscousS;
    }
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCMMS::haveExactSol() const
{
  return true;
}

/*--------------------------------------------------------------------*/
//  Compute the exact solution state \<U\> in the cells
/**
 *  \param[out] a_Ux    Exact solution
 *  \param[in]  a_box   Cell locations to update
 *  \param[in]  a_disjointBox
 *                      Current disjoint box
 *  \param[in]  a_gridMetrics
 *                      Level grid metrics
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_didx  Current DataIndex on the current disjoint box
 *  \param[in]  a_time  Current solution time
 *  \param[in]  a_level Grid level
 *  \return             0 - Successfully computed exact solution
 *                      1 - Exact solution is not known
 *//*-----------------------------------------------------------------*/

int
CNSIBCMMS::exactSol(FArrayBox&              a_Ux,
                    const Box&              a_box,
                    const Box&              a_disjointBox,
                    const LevelGridMetrics& a_gridMetrics,
                    const FluxBox&          a_unitNormals,
                    const DataIndex&        a_didx,
                    const Real              a_time,
                    const int               a_level) const
{
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_disjointBox);
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  Box box1Dom = grow(a_box, 1);
  box1Dom &= blockDomain;
  Box box2Dom = grow(a_box, 2);
  box2Dom &= blockDomain;
  Box initBox = box2Dom;
  // Get physical coordinates
  FABSTACKTEMP(XiFab, initBox, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, initBox, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(initBox, XiFab, XFab, blockCoordSys);
  const int numWVar  = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  // Obtain the mean values and the wave-amplitude values
  const CRDState& mean = CRDState::get(m_idxStateInit);
  const Real rM = mean.density();
  const RealVect velM = mean.velocity();
  const Real pM = mean.pressure();
  const Real TM = mean.temperature();

  // Pointwise values of W
  FABSTACKTEMP(Wc, initBox, numWVar);
  Wc.setVal(-1.);

  MD_BOXLOOP(initBox, i)
    {
      // Physical location
      const RealVect loc(
        D_DECL(XFab[MD_IX(i, 0)], XFab[MD_IX(i, 1)], XFab[MD_IX(i, 2)]));
      // Spatial location normalized by characteristic length
      const RealVect sLoc = loc/m_charLength;
      // Time normalized by characteristic time
      const Real tLoc = a_time/m_charTime;
      // Spatial argument for density field
      const RealVect rhoArg = 2.*PI*m_sOmegaRho*sLoc + m_sPhiRho;
      // Temporal argument for density field
      const Real rhoTimeArg = 2.*PI*m_tOmegaRho*tLoc + m_tPhiRho;
      // Spatial + temporal arguments for velocity fields
      D_TERM(const RealVect uArg = 2.*PI*m_sOmegaU*sLoc + m_sPhiU;
             const Real uTimeArg = 2.*PI*m_tOmegaU*tLoc + m_tPhiU;,
             const RealVect vArg = 2.*PI*m_sOmegaV*sLoc + m_sPhiV;
             const Real vTimeArg = 2.*PI*m_tOmegaV*tLoc + m_tPhiV;,
             const RealVect wArg = 2.*PI*m_sOmegaW*sLoc + m_sPhiW;
             const Real wTimeArg = 2.*PI*m_tOmegaW*tLoc + m_tPhiW;);
      // Spatial argument for pressure field
      const RealVect pArg = 2.*PI*m_sOmegaP*sLoc + m_sPhiP;
      // Temporal argument for pressure field
      const Real pTimeArg = 2.*PI*m_tOmegaP*tLoc + m_tPhiP;
      // Spatial argument for temperature field
      const RealVect TArg = 2.*PI*m_sOmegaT*sLoc + m_sPhiT;
      // Temporal argument for temperature field
      const Real TTimeArg = 2.*PI*m_tOmegaT*tLoc + m_tPhiT;
      // Define the density field
      Wc[MD_IX(i, rhoIndx)] = rM + m_rA*(std::sin(rhoTimeArg))*D_TERM(
        std::sin(rhoArg[0]),*std::sin(rhoArg[1]),*std::sin(rhoArg[2]));
      // Define the velocity fields -- first, the X-velocity field
      D_TERM(Wc[MD_IX(i, WvelIndx)]     = velM[0] +
             m_uA[0]*(std::sin(uTimeArg))*D_TERM(
               std::sin(uArg[0]),*std::sin(uArg[1]),*std::sin(uArg[2]));,
             // Y-velocity field
             Wc[MD_IX(i, WvelIndx + 1)] = velM[1] +
             m_uA[1]*(std::sin(vTimeArg))*D_TERM(
               std::sin(vArg[0]),*std::sin(vArg[1]),*std::sin(vArg[2]));,
             // Z-velocity field
             Wc[MD_IX(i, WvelIndx + 2)] = velM[2] +
             m_uA[2]*(std::sin(wTimeArg))*D_TERM(
               std::sin(wArg[0]),*std::sin(wArg[1]),*std::sin(wArg[2])););
      if (!m_bcTestCase)
        {
          // Define the pressure field
          Wc[MD_IX(i, presIndx)] = pM + m_pA*(std::sin(pTimeArg))*D_TERM(
            std::sin(pArg[0]),*std::sin(pArg[1]),*std::sin(pArg[2]));
        }
      else
        {
          // Define the temperature field
          Wc[MD_IX(i, tempIndx)] = TM + m_TA*(std::sin(TTimeArg))*D_TERM(
            std::sin(TArg[0]),*std::sin(TArg[1]),*std::sin(TArg[2]));
        }
    }
  // Initialize the values in UFab
  CRDparam::g_CRDPhysics->initialize(a_Ux,
                                     Wc,
                                     a_gridMetrics,
                                     a_unitNormals,
                                     a_didx,
                                     a_disjointBox,
                                     box1Dom);
  fourthOrderAverageCell(a_Ux, blockDomain, a_box);

  if (m_sixthOrderInit)
    {
      // Testing another approach
      FABSTACKTEMP(UFab, initBox, a_Ux.nComp());
      UFab.setVal(0.);
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals,
                                         a_didx,
                                         a_disjointBox,
                                         box2Dom);
      //**NOTE: Adding 6th-order convolution function just for this
      const Real d2xFactor = 1./24.;
      const Real d4xFactor = 1./1920.;
      for (int comp = 0; comp != a_Ux.nComp(); ++comp)
        {
          // Compute the 2nd-derivative and the 4th-derivative using 5 cells
          FABSTACKTEMP(D2Fab, a_box, 1);
          FABSTACKTEMP(D4Fab, a_box, 1);
          D2Fab.setVal(0.);
          D4Fab.setVal(0.);
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              Box loBox, nextLoBox, hiBox, nextHiBox, centerBox, innerCenterBox,
                entireBox;
              int hasLo, hasHi;
              Box inBox(a_box);
              inBox.grow(dir, 1);
              loHiCenter5(loBox, nextLoBox, hasLo, hiBox, nextHiBox, hasHi,
                          centerBox, innerCenterBox, entireBox, inBox,
                          blockDomain, dir);
              const int MD_ID(o, dir);
              MD_BOXLOOP(innerCenterBox, i)
                {
                  Real u_m2 = UFab[MD_OFFSETIX(i,-,2*o,comp)];
                  Real u_m1 = UFab[MD_OFFSETIX(i,-,o,comp)];
                  Real u    = UFab[MD_IX(i,comp)];
                  Real u_p1 = UFab[MD_OFFSETIX(i,+,o,comp)];
                  Real u_p2 = UFab[MD_OFFSETIX(i,+,2*o,comp)];
                  D2Fab[MD_IX(i, 0)] +=
                    (-u_m2 + 16.*u_m1 - 30.*u + 16.*u_p1 - u_p2)/12.;
                  D4Fab[MD_IX(i, 0)] += u_m2 - 4.*u_m1 + 6.*u - 4*u_p1 + u_p2;
                }
              if (hasLo)
                {
                  MD_BOXLOOP(loBox, i)
                    {
                      Real u    = UFab[MD_IX(i,comp)];
                      Real u_p1 = UFab[MD_OFFSETIX(i,+,o,comp)];
                      Real u_p2 = UFab[MD_OFFSETIX(i,+,2*o,comp)];
                      Real u_p3 = UFab[MD_OFFSETIX(i,+,3*o,comp)];
                      Real u_p4 = UFab[MD_OFFSETIX(i,+,4*o,comp)];
                      D2Fab[MD_IX(i, 0)] +=
                        (35.*u - 104.*u_p1 + 114.*u_p2 - 56.*u_p3
                         + 11.*u_p4)/12.;
                      D4Fab[MD_IX(i, 0)] +=
                        u - 4.*u_p1 + 6.*u_p2 - 4.*u_p3 + u_p4;
                    }
                  MD_BOXLOOP(nextLoBox, i)
                    {
                      Real u_m1 = UFab[MD_OFFSETIX(i,-,o,comp)];
                      Real u    = UFab[MD_IX(i,comp)];
                      Real u_p1 = UFab[MD_OFFSETIX(i,+,o,comp)];
                      Real u_p2 = UFab[MD_OFFSETIX(i,+,2*o,comp)];
                      Real u_p3 = UFab[MD_OFFSETIX(i,+,3*o,comp)];
                      D2Fab[MD_IX(i, 0)] +=
                        (11.*u_m1 - 20.*u + 6.*u_p1 + 4.*u_p2 - u_p3)/12.;
                      D4Fab[MD_IX(i, 0)] +=
                        u_m1 - 4.*u + 6.*u_p1 - 4.*u_p2 + u_p3;
                    }
                }
              if (hasHi)
                {
                  MD_BOXLOOP(hiBox, i)
                    {
                      Real u_m4 = UFab[MD_OFFSETIX(i,-,4*o,comp)];
                      Real u_m3 = UFab[MD_OFFSETIX(i,-,3*o,comp)];
                      Real u_m2 = UFab[MD_OFFSETIX(i,-,2*o,comp)];
                      Real u_m1 = UFab[MD_OFFSETIX(i,-,o,comp)];
                      Real u    = UFab[MD_IX(i,comp)];
                      D2Fab[MD_IX(i, 0)] +=
                        (11.*u_m4 - 56.*u_m3 + 114.*u_m2 - 104.*u_m1
                         + 35.*u)/12.;
                      D4Fab[MD_IX(i, 0)] +=
                        u_m4 - 4.*u_m3 + 6.*u_m2 - 4.*u_m1 + u;
                    }
                  MD_BOXLOOP(nextHiBox, i)
                    {
                      Real u_m3 = UFab[MD_OFFSETIX(i,-,3*o,comp)];
                      Real u_m2 = UFab[MD_OFFSETIX(i,-,2*o,comp)];
                      Real u_m1 = UFab[MD_OFFSETIX(i,-,o,comp)];
                      Real u    = UFab[MD_IX(i,comp)];
                      Real u_p1 = UFab[MD_OFFSETIX(i,+,o,comp)];
                      D2Fab[MD_IX(i, 0)] +=
                        (-u_m3 + 4.*u_m2 + 6.*u_m1 - 20.*u + 11.*u_p1)/12.;
                      D4Fab[MD_IX(i, 0)] +=
                        u_m3 - 4.*u_m2 + 6.*u_m1 - 4.*u + u_p1;
                    }
                }
            }
          MD_BOXLOOP(a_box, i)
            {
              a_Ux[MD_IX(i, comp)] = UFab[MD_IX(i, comp)]
                + d2xFactor*D2Fab[MD_IX(i, 0)] + d4xFactor*D4Fab[MD_IX(i, 0)];
            }
        }
    }

  if (CRDparam::g_cartesian && m_exactInit)
    {
      // First, get node locations of the cells
      //**NOTE: a special box is necessary where the high-end has been
      //        grown by 1 so that the high nodes can be accounted for
      Box boxDomNodes = a_box;
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          boxDomNodes.growHi(dir, 1);
        }
      FABSTACKTEMP(XiFabNodes, boxDomNodes, SpaceDim);
      FABSTACKTEMP(XFabNodes, boxDomNodes, SpaceDim);
      this->CNSIBC::getNodeCoordinates(boxDomNodes,
                                       XiFabNodes,
                                       XFabNodes,
                                       blockCoordSys);
      Real R_cell = CRDparam::g_R;
      Real gamma = CRDparam::g_gamma;
      MD_BOXLOOP(a_box, i)
        {
          XFabNodes[MD_OFFSETIV(i,+,IntVect::Unit, 0)];
          // Physical location of high side of cell
          RealVect highLoc(
            D_DECL(XFabNodes[MD_OFFSETIV(i,+,IntVect::Unit, 0)],
                   XFabNodes[MD_OFFSETIV(i,+,IntVect::Unit, 1)],
                   XFabNodes[MD_OFFSETIV(i,+,IntVect::Unit, 2)]));
          // Physical location of low side of cell
          RealVect lowLoc(D_DECL(XFabNodes[MD_IX(i, 0)],
                                 XFabNodes[MD_IX(i, 1)],
                                 XFabNodes[MD_IX(i, 2)]));
          // Length of the cell in each direction
          RealVect length = highLoc - lowLoc;
          // Cell volume
          Real vol = length.product();
          // Spatial location normalized by characteristic length
          const RealVect sLocHigh = highLoc/m_charLength;
          const RealVect sLocLow  = lowLoc/m_charLength;
          // Time normalized by characteristic time
          const Real tLoc = a_time/m_charTime;

          // Spatial argument for density field -- high side evaluation
          const RealVect rhoArgH =
            2.*PI*m_sOmegaRho*sLocHigh + m_sPhiRho;
          // Spatial argument for density field -- low side evaluation
          const RealVect rhoArgL  =
            2.*PI*m_sOmegaRho*sLocLow + m_sPhiRho;
          // Derivative of spatial argument for density field
          const RealVect rhoArgD = 2.*PI*m_sOmegaRho;
          // Temporal argument for density field
          const Real rhoTimeArg = 2.*PI*m_tOmegaRho*tLoc + m_tPhiRho;
          // Cell-averaged density
          // Note there are 2 different cases for this integral
          // 1) m_sOmegaRho[0] = 0
          //    -- this requires just multiplying the sinusoid by length
          // 2) m_sOmegaRho[0] != 0
          //    -- this requires just integrating the sinusoid
          RealVect rho_integralDir = RealVect::Zero;
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              if (m_sOmegaRho[dir] == 0)
                {
                  rho_integralDir[dir] =
                    std::sin(rhoArgL[dir])*length[dir];
                }
              else
                {
                  rho_integralDir[dir] =
                    (std::cos(rhoArgL[dir])
                     - std::cos(rhoArgH[dir]))/rhoArgD[dir];
                }
            }
          Real rho_cellAvg = rM
            + (m_rA/vol)*(std::sin(rhoTimeArg))*rho_integralDir.product();

          // Set the cell-averaged density to the exact solution
          a_Ux[MD_IX(i, rhoIndx)] = rho_cellAvg;

          // Spatial + temporal arguments for velocity fields
          D_TERM(
            const RealVect uArgH = 2.*PI*m_sOmegaU*sLocHigh + m_sPhiU;
            const RealVect uArgL = 2.*PI*m_sOmegaU*sLocLow + m_sPhiU;
            const RealVect uArgD = 2.*PI*m_sOmegaU;
            const Real uTimeArg  = 2.*PI*m_tOmegaU*tLoc + m_tPhiU;,
            const RealVect vArgH = 2.*PI*m_sOmegaV*sLocHigh + m_sPhiV;
            const RealVect vArgL = 2.*PI*m_sOmegaV*sLocLow + m_sPhiV;
            const RealVect vArgD = 2.*PI*m_sOmegaV;
            const Real vTimeArg  = 2.*PI*m_tOmegaV*tLoc + m_tPhiV;,
            const RealVect wArgH = 2.*PI*m_sOmegaW*sLocHigh + m_sPhiW;
            const RealVect wArgL = 2.*PI*m_sOmegaW*sLocLow + m_sPhiW;
            const RealVect wArgD = 2.*PI*m_sOmegaW;
            const Real wTimeArg  = 2.*PI*m_tOmegaW*tLoc + m_tPhiW;);
          // Define the spatial integrals of the velocity sinusoids
          // Note there are 2 different cases for these integrals
          // 1) m_sOmegaU[0] = 0
          //    -- this requires just multiplying the sinusoid by length
          // 2) m_sOmegaU[0] != 0
          //    -- this requires just integrating the sinusoid
          D_TERM(
            RealVect u_integralDir = RealVect::Zero;
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if (m_sOmegaU[dir] == 0)
                  {
                    u_integralDir[dir] = std::sin(uArgL[dir])*length[dir];
                  }
                else
                  {
                    u_integralDir[dir] = (std::cos(uArgL[dir])
                                          - std::cos(uArgH[dir]))/uArgD[dir];
                  }
              },
            RealVect v_integralDir = RealVect::Zero;
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if (m_sOmegaV[dir] == 0)
                  {
                    v_integralDir[dir] = std::sin(vArgL[dir])*length[dir];
                  }
                else
                  {
                    v_integralDir[dir] = (std::cos(vArgL[dir])
                                          - std::cos(vArgH[dir]))/vArgD[dir];
                  }
              },
            RealVect w_integralDir = RealVect::Zero;
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if (m_sOmegaW[dir] == 0)
                  {
                    w_integralDir[dir] = std::sin(wArgL[dir])*length[dir];
                  }
                else
                  {
                    w_integralDir[dir] = (std::cos(wArgL[dir])
                                          - std::cos(wArgH[dir]))/wArgD[dir];
                  }
              });
          // Define the momentum fields -- first, the X-momentum field
          D_TERM(
            Real u_0 = rM*velM[0]; // mean rho*u component
            // Mean u*rho_sinusoid
            Real u_1 = velM[0]*(m_rA/vol)*(std::sin(rhoTimeArg))
            *rho_integralDir.product();
            // Mean rho*u_sinusoid
            Real u_2 = rM*(m_uA[0]/vol)*(std::sin(uTimeArg))
            *u_integralDir.product();
            // rho_sinusoid*u_sinusoid
            // Note there are 5 different cases for this integral
            // 1) m_sOmegaU[0] = 0 and m_sOmegaRho = 0
            //    -- this requires just multiplying the two sinusoids
            // 2) m_sOmegaU[0] = 0 and m_sOmegaRho != 0
            //    -- this requires integrating just the rho sinusoid
            // 3) m_sOmegaU[0] != 0 and m_sOmegaRho = 0
            //    -- this requires integrating just the u sinusoid
            // 4) m_sOmegaU[0] != 0 and m_sOmegaRho != 0 and
            //    m_sOmegaU[0] = m_sOmegaRho[0]
            //    -- this requires integration using same frequencies
            // 5) m_sOmegaU[0] != 0 and m_sOmegaRho != 0 and
            //    m_sOmegaU[0] != m_sOmegaRho[0]
            //    -- this requires integration using different frequencies
            RealVect rhoU_integralDir = RealVect::Zero;
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if ((uArgD[dir] == 0) && (rhoArgD[dir] == 0))
                  {
                    rhoU_integralDir[dir] = std::sin(m_sPhiU[dir])
                      *std::sin(m_sPhiRho[dir])*length[dir];
                  }
                else if ((uArgD[dir] == 0) && (rhoArgD[dir] != 0))
                  {
                    rhoU_integralDir[dir] =
                      std::sin(m_sPhiU[dir])*rho_integralDir[dir];
                  }
                else if ((uArgD[dir] != 0) && (rhoArgD[dir] == 0))
                  {
                    rhoU_integralDir[dir] =
                      u_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                  }
                else if (uArgD[dir] == rhoArgD[dir])
                  {
                    rhoU_integralDir[dir] = (1./(4.*rhoArgD[dir]))*(
                      2.*rhoArgD[dir]*sLocHigh[dir]
                      *std::cos(m_sPhiRho[dir] - m_sPhiU[dir])
                      - std::sin(rhoArgH[dir] + uArgH[dir])
                      - 2.*rhoArgD[dir]*sLocLow[dir]
                      *std::cos(m_sPhiRho[dir] - m_sPhiU[dir])
                      + std::sin(rhoArgL[dir] + uArgL[dir]));
                  }
                else
                  {
                    rhoU_integralDir[dir] = 0.5*(
                      std::sin(rhoArgH[dir]-uArgH[dir])/(rhoArgD[dir]
                                                         -uArgD[dir])
                      - std::sin(rhoArgH[dir]+uArgH[dir])/(rhoArgD[dir]
                                                           +uArgD[dir])
                      - std::sin(rhoArgL[dir]-uArgL[dir])/(rhoArgD[dir]
                                                           -uArgD[dir])
                      + std::sin(rhoArgL[dir]+uArgL[dir])/(rhoArgD[dir]
                                                           +uArgD[dir]));
                  }
              }
            Real u_4 = m_rA*(m_uA[0]/vol)
            *(std::sin(rhoTimeArg))*(std::sin(uTimeArg))
            *(rhoU_integralDir.product());
            Real rhoU_cellAvg = u_0 + u_1 + u_2 + u_4;,
            // Y-momentum field
            Real v_0 = rM*velM[1]; // mean rho*v component
            // Mean v*rho_sinusoid
            Real v_1 = velM[1]*(m_rA/vol)*(std::sin(rhoTimeArg))
            *rho_integralDir.product();
            // Mean rho*v_sinusoid
            Real v_2 = rM*(m_uA[1]/vol)*(std::sin(vTimeArg))
            *v_integralDir.product();
            // rho_sinusoid*v_sinusoid
            // Note there are 5 different cases for this integral
            // 1) m_sOmegaV[0] = 0 and m_sOmegaRho = 0
            //    -- this requires just multiplying the two sinusoids
            // 2) m_sOmegaV[0] = 0 and m_sOmegaRho != 0
            //    -- this requires integrating just the rho sinusoid
            // 3) m_sOmegaV[0] != 0 and m_sOmegaRho = 0
            //    -- this requires integrating just the v sinusoid
            // 4) m_sOmegaV[0] != 0 and m_sOmegaRho != 0 and
            //    m_sOmegaV[0] = m_sOmegaRho[0]
            //    -- this requires integration using same frequencies
            // 5) m_sOmegaV[0] != 0 and m_sOmegaRho != 0 and
            //    m_sOmegaV[0] != m_sOmegaRho[0]
            //    -- this requires integration using different frequencies
            RealVect rhoV_integralDir = RealVect::Zero;
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if ((vArgD[dir] == 0) && (rhoArgD[dir] == 0))
                  {
                    rhoV_integralDir[dir] = std::sin(m_sPhiV[dir])
                      *std::sin(m_sPhiRho[dir])*length[dir];
                  }
                else if ((vArgD[dir] == 0) && (rhoArgD[dir] != 0))
                  {
                    rhoV_integralDir[dir] =
                      std::sin(m_sPhiV[dir])*rho_integralDir[dir];
                  }
                else if ((vArgD[dir] != 0) && (rhoArgD[dir] == 0))
                  {
                    rhoV_integralDir[dir] =
                      v_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                  }
                else if (vArgD[dir] == rhoArgD[dir])
                  {
                    rhoV_integralDir[dir] = (1./(4.*rhoArgD[dir]))*(
                      2.*rhoArgD[dir]*sLocHigh[dir]
                      *std::cos(m_sPhiRho[dir] - m_sPhiV[dir])
                      - std::sin(rhoArgH[dir] + vArgH[dir])
                      - 2.*rhoArgD[dir]*sLocLow[dir]
                      *std::cos(m_sPhiRho[dir] - m_sPhiV[dir])
                      + std::sin(rhoArgL[dir] + vArgL[dir]));
                  }
                else
                  {
                    rhoV_integralDir[dir] = 0.5*(
                      std::sin(rhoArgH[dir]-vArgH[dir])/(rhoArgD[dir]
                                                         -vArgD[dir])
                      - std::sin(rhoArgH[dir]+vArgH[dir])/(rhoArgD[dir]
                                                           +vArgD[dir])
                      - std::sin(rhoArgL[dir]-vArgL[dir])/(rhoArgD[dir]
                                                           -vArgD[dir])
                      + std::sin(rhoArgL[dir]+vArgL[dir])/(rhoArgD[dir]
                                                           +vArgD[dir]));
                  }
              }
            Real v_4 = m_rA*(m_uA[1]/vol)
            *(std::sin(rhoTimeArg))*(std::sin(vTimeArg))
            *(rhoV_integralDir.product());
            Real rhoV_cellAvg = v_0 + v_1 + v_2 + v_4;,
            // Z-momentum field
            Real w_0 = rM*velM[2]; // mean rho*w component
            // Mean w*rho_sinusoid
            Real w_1 = velM[2]*(m_rA/vol)*(std::sin(rhoTimeArg))
            *rho_integralDir.product();
            // Mean rho*w_sinusoid
            Real w_2 = rM*(m_uA[2]/vol)*(std::sin(wTimeArg))
            *w_integralDir.product();
            // rho_sinusoid*w_sinusoid
            // Note there are 5 different cases for this integral
            // 1) m_sOmegaW[0] = 0 and m_sOmegaRho = 0
            //    -- this requires just multiplying the two sinusoids
            // 2) m_sOmegaW[0] = 0 and m_sOmegaRho != 0
            //    -- this requires integrating just the rho sinusoid
            // 3) m_sOmegaW[0] != 0 and m_sOmegaRho = 0
            //    -- this requires integrating just the w sinusoid
            // 4) m_sOmegaW[0] != 0 and m_sOmegaRho != 0 and
            //    m_sOmegaW[0] = m_sOmegaRho[0]
            //    -- this requires integration using same frequencies
            // 5) m_sOmegaW[0] != 0 and m_sOmegaRho != 0 and
            //    m_sOmegaW[0] != m_sOmegaRho[0]
            //    -- this requires integration using different frequencies
            RealVect rhoW_integralDir = RealVect::Zero;
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if ((wArgD[dir] == 0) && (rhoArgD[dir] == 0))
                  {
                    rhoW_integralDir[dir] = std::sin(m_sPhiW[dir])
                      *std::sin(m_sPhiRho[dir])*length[dir];
                  }
                else if ((wArgD[dir] == 0) && (rhoArgD[dir] != 0))
                  {
                    rhoW_integralDir[dir] =
                      std::sin(m_sPhiW[dir])*rho_integralDir[dir];
                  }
                else if ((wArgD[dir] != 0) && (rhoArgD[dir] == 0))
                  {
                    rhoW_integralDir[dir] =
                      w_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                  }
                else if (wArgD[dir] == rhoArgD[dir])
                  {
                    rhoW_integralDir[dir] = (1./(4.*rhoArgD[dir]))*(
                      2.*rhoArgD[dir]*sLocHigh[dir]
                      *std::cos(m_sPhiRho[dir] - m_sPhiW[dir])
                      - std::sin(rhoArgH[dir] + wArgH[dir])
                      - 2.*rhoArgD[dir]*sLocLow[dir]
                      *std::cos(m_sPhiRho[dir] - m_sPhiW[dir])
                      + std::sin(rhoArgL[dir] + wArgL[dir]));
                  }
                else
                  {
                    rhoW_integralDir[dir] = 0.5*(
                      std::sin(rhoArgH[dir]-wArgH[dir])/(rhoArgD[dir]
                                                         -wArgD[dir])
                      - std::sin(rhoArgH[dir]+wArgH[dir])/(rhoArgD[dir]
                                                           +wArgD[dir])
                      - std::sin(rhoArgL[dir]-wArgL[dir])/(rhoArgD[dir]
                                                           -wArgD[dir])
                      + std::sin(rhoArgL[dir]+wArgL[dir])/(rhoArgD[dir]
                                                           +wArgD[dir]));
                  }
              }
            Real w_4 = m_rA*(m_uA[2]/vol)
            *(std::sin(rhoTimeArg))*(std::sin(wTimeArg))
            *(rhoW_integralDir.product());
            Real rhoW_cellAvg = w_0 + w_1 + w_2 + w_4;);

          // Set the cell-averaged momentum to the exact solution
          D_TERM(a_Ux[MD_IX(i, WvelIndx)]   = rhoU_cellAvg;,
                 a_Ux[MD_IX(i, WvelIndx+1)] = rhoV_cellAvg;,
                 a_Ux[MD_IX(i, WvelIndx+2)] = rhoW_cellAvg;);

          // Spatial argument for pressure field -- high side evaluation
          const RealVect pArgH = 2.*PI*m_sOmegaP*sLocHigh + m_sPhiP;
          // Spatial argument for pressure field -- low side evaluation
          const RealVect pArgL  = 2.*PI*m_sOmegaP*sLocLow + m_sPhiP;
          // Derivative of spatial argument for pressure field
          const RealVect pArgD = 2.*PI*m_sOmegaP;
          // Temporal argument for pressure field
          const Real pTimeArg = 2.*PI*m_tOmegaP*tLoc + m_tPhiP;
          Real p_cellAvg = 0.;

          // Spatial argument for temperature field -- high side evaluation
          const RealVect TArgH = 2.*PI*m_sOmegaT*sLocHigh + m_sPhiT;
          // Spatial argument for temperature field -- low side evaluation
          const RealVect TArgL  = 2.*PI*m_sOmegaT*sLocLow + m_sPhiT;
          // Derivative of spatial argument for temperature field
          const RealVect TArgD = 2.*PI*m_sOmegaT;
          // Temporal argument for temperature field
          const Real TTimeArg = 2.*PI*m_tOmegaT*tLoc + m_tPhiT;

          if (!m_bcTestCase)
            {
              // Note there are 2 different cases for this integral
              // 1) m_sOmegaP[0] = 0
              //    -- this requires multiplying the sinusoid by length
              // 2) m_sOmegaP[0] != 0
              //    -- this requires integrating the sinusoid
              RealVect p_integralDir = RealVect::Zero;
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  if (m_sOmegaP[dir] == 0)
                    {
                      p_integralDir[dir] =
                        std::sin(pArgL[dir])*length[dir];
                    }
                  else
                    {
                      p_integralDir[dir] =
                        (std::cos(pArgL[dir])
                         - std::cos(pArgH[dir]))/pArgD[dir];
                    }
                }
              p_cellAvg = pM
                + (m_pA/vol)*(std::sin(pTimeArg))*p_integralDir.product();
            }
          else
            {
              // Note there are 2 different cases for this integral
              // 1) m_sOmegaT[0] = 0
              //    -- this requires multiplying the sinusoid by length
              // 2) m_sOmegaT[0] != 0
              //    -- this requires integrating the sinusoid
              RealVect T_integralDir = RealVect::Zero;
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  if (m_sOmegaT[dir] == 0)
                    {
                      T_integralDir[dir] =
                        std::sin(TArgL[dir])*length[dir];
                    }
                  else
                    {
                      T_integralDir[dir] =
                        (std::cos(TArgL[dir])
                         - std::cos(TArgH[dir]))/TArgD[dir];
                    }
                }
              // Temperature*R*rho
              Real T_0 = rM*TM; // mean rho*T component
              // Mean T*rho_sinusoid
              Real T_1 = TM*(m_rA/vol)*(std::sin(rhoTimeArg))
                *rho_integralDir.product();
              // Mean rho*T_sinusoid
              Real T_2 = rM*(m_TA/vol)*(std::sin(TTimeArg))
                *T_integralDir.product();
              // rho_sinusoid*T_sinusoid
              // Note there are 5 different cases for this integral
              // 1) m_sOmegaT[0] = 0 and m_sOmegaRho = 0
              //    -- this requires just multiplying the two sinusoids
              // 2) m_sOmegaT[0] = 0 and m_sOmegaRho != 0
              //    -- this requires integrating just the rho sinusoid
              // 3) m_sOmegaT[0] != 0 and m_sOmegaRho = 0
              //    -- this requires integrating just the u sinusoid
              // 4) m_sOmegaT[0] != 0 and m_sOmegaRho != 0 and
              //    m_sOmegaT[0] = m_sOmegaRho[0]
              //    -- this requires integration using same frequencies
              // 5) m_sOmegaT[0] != 0 and m_sOmegaRho != 0 and
              //    m_sOmegaT[0] != m_sOmegaRho[0]
              //    -- this requires integration using different frequencies
              RealVect rhoT_integralDir = RealVect::Zero;
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  if ((TArgD[dir] == 0) && (rhoArgD[dir] == 0))
                    {
                      rhoT_integralDir[dir] = std::sin(m_sPhiT[dir])
                        *std::sin(m_sPhiRho[dir])*length[dir];
                    }
                  else if ((TArgD[dir] == 0) && (rhoArgD[dir] != 0))
                    {
                      rhoT_integralDir[dir] =
                        std::sin(m_sPhiT[dir])*rho_integralDir[dir];
                    }
                  else if ((TArgD[dir] != 0) && (rhoArgD[dir] == 0))
                    {
                      rhoT_integralDir[dir] =
                        T_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                    }
                  else if (TArgD[dir] == rhoArgD[dir])
                    {
                      rhoT_integralDir[dir] = (1./(4.*rhoArgD[dir]))*(
                        2.*rhoArgD[dir]*sLocHigh[dir]
                        *std::cos(m_sPhiRho[dir] - m_sPhiT[dir])
                        - std::sin(rhoArgH[dir] + TArgH[dir])
                        - 2.*rhoArgD[dir]*sLocLow[dir]
                        *std::cos(m_sPhiRho[dir] - m_sPhiT[dir])
                        + std::sin(rhoArgL[dir] + TArgL[dir]));
                    }
                  else
                    {
                      rhoT_integralDir[dir] = 0.5*(
                        std::sin(rhoArgH[dir]-TArgH[dir])/(rhoArgD[dir]
                                                           -TArgD[dir])
                        - std::sin(rhoArgH[dir]+TArgH[dir])/(rhoArgD[dir]
                                                             +TArgD[dir])
                        - std::sin(rhoArgL[dir]-TArgL[dir])/(rhoArgD[dir]
                                                             -TArgD[dir])
                        + std::sin(rhoArgL[dir]+TArgL[dir])/(rhoArgD[dir]
                                                             +TArgD[dir]));
                    }
                }
              Real T_4 = m_rA*(m_TA/vol)
                *(std::sin(rhoTimeArg))*(std::sin(TTimeArg))
                *(rhoT_integralDir.product());
              Real rhoT_cellAvg = T_0 + T_1 + T_2 + T_4;
              p_cellAvg = rhoT_cellAvg*R_cell;
            }

          // Now for the fun part -- cell-averaged kinetic energy
          D_TERM(
            Real uu_0 = rM*velM[0]*velM[0]; // mean rho*u*u component
            // Mean u*u*rho_sinusoid
            Real uu_1 = velM[0]*velM[0]*(m_rA/vol)
            *(std::sin(rhoTimeArg))*rho_integralDir.product();
            // Mean rho*u*u_sinusoid
            Real uu_2 = 2.*velM[0]*rM*(m_uA[0]/vol)*(std::sin(uTimeArg))
            *u_integralDir.product();
            // Mean u*rho_sinusoid*u_sinusoid
            Real uu_3 = 2.*velM[0]*m_rA*(m_uA[0]/vol)
            *(std::sin(rhoTimeArg))*(std::sin(uTimeArg))
            *(rhoU_integralDir.product());
            // Mean rho*u_sinusoid*u_sinusoid
            RealVect UU_integralDir = RealVect::Zero;
            // Note: because this is u*u there are only two options
            // 1) m_sOmegaU[0] = 0
            //    -- this requires just multiplying the two sinusoids
            // 2) m_sOmegaU[0] != 0
            //    -- this requires integration using same frequencies
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if (uArgD[dir] == 0)
                  {
                    UU_integralDir[dir] = std::sin(m_sPhiU[dir])
                      *std::sin(m_sPhiU[dir])*length[dir];
                  }
                else
                  {
                    UU_integralDir[dir] = (1./(4.*uArgD[dir]))*(
                      2.*uArgH[dir] - std::sin(2.*uArgH[dir])
                      - 2.*uArgL[dir] + std::sin(2.*uArgL[dir]));
                  }
              }
            Real uu_4 = rM*m_uA[0]*(m_uA[0]/vol)
            *(std::sin(uTimeArg))*(std::sin(uTimeArg))
            *(UU_integralDir.product());
            // Rho_sinusoid*u_sinusoid*u_sinusoid
            RealVect rhoUU_integralDir = RealVect::Zero;
            // Note: there are 15 cases here, 10 of them are unnecessary in
            //       our case (r*u*u will always have u = u). Additionally
            //       case 5 has two variants
            //       a) m_sOmegaRho = 2.*m_sOmegaU[0]
            //       b) m_sOmegaRho != 2.*m_sOmegaU[0]
            // 1) m_sOmegaU[0] = 0 and m_sOmegaRho = 0
            //    -- this requires just multiplying the two sinusoids
            // 2) m_sOmegaU[0] = 0 and m_sOmegaRho != 0
            //    -- this requires integrating just the rho sinusoid
            // 3) m_sOmegaU[0] != 0 and m_sOmegaRho = 0
            //    -- this requires integrating just the u sinusoid
            // 4) m_sOmegaU[0] != 0 and m_sOmegaRho != 0 and
            //    m_sOmegaU[0] = m_sOmegaRho[0]
            //    -- this requires integration using same frequencies
            // 5) m_sOmegaU[0] != 0 and m_sOmegaRho != 0 and
            //    m_sOmegaU[0] != m_sOmegaRho[0]
            //    -- this requires integration using different frequencies
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if ((uArgD[dir] == 0) && (rhoArgD[dir] == 0))
                  {
                    rhoUU_integralDir[dir] =
                      std::sin(m_sPhiU[dir])*std::sin(m_sPhiU[dir])
                      *std::sin(m_sPhiRho[dir])*length[dir];
                  }
                else if ((uArgD[dir] == 0) && (rhoArgD[dir] != 0))
                  {
                    rhoUU_integralDir[dir] =
                      std::sin(m_sPhiU[dir])*std::sin(m_sPhiU[dir])
                      *rho_integralDir[dir];
                  }
                else if ((uArgD[dir] != 0) && (rhoArgD[dir] == 0))
                  {
                    rhoUU_integralDir[dir] =
                      UU_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                  }
                else if (uArgD[dir] == rhoArgD[dir])
                  {
                    rhoUU_integralDir[dir] = (1./(12.*uArgD[dir]))
                      *(std::cos(rhoArgH[dir] + uArgH[dir] + uArgH[dir])
                        - 3.*std::cos(rhoArgH[dir] - uArgH[dir] + uArgH[dir])
                        - 3.*std::cos(rhoArgH[dir] + uArgH[dir] - uArgH[dir])
                        - 3.*std::cos(rhoArgH[dir] - uArgH[dir] - uArgH[dir])
                        - std::cos(rhoArgL[dir] + uArgL[dir] + uArgL[dir])
                        + 3.*std::cos(rhoArgL[dir] - uArgL[dir] + uArgL[dir])
                        + 3.*std::cos(rhoArgL[dir] + uArgL[dir] - uArgL[dir])
                        + 3.*std::cos(rhoArgL[dir] - uArgL[dir] - uArgL[dir]));
                  }
                else
                  {
                    if (rhoArgD[dir] == (2.*uArgD[dir]))
                      {
                        rhoUU_integralDir[dir] =
                          (std::cos(rhoArgH[dir] + 2.*uArgH[dir])
                           - 4.*std::cos(rhoArgH[dir])
                           - 2.*rhoArgD[dir]*sLocHigh[dir]
                           *std::sin(rhoArgH[dir] - 2.*uArgH[dir])
                           - std::cos(rhoArgL[dir] + 2.*uArgL[dir])
                           + 4.*std::cos(rhoArgL[dir])
                           + 2.*rhoArgD[dir]*sLocLow[dir]
                           *std::sin(rhoArgL[dir] - 2.*uArgL[dir])
                            )/(16.*uArgD[dir]);
                      }
                    else
                      {
                        rhoUU_integralDir[dir] =
                          (std::cos(rhoArgH[dir] - 2.*uArgH[dir])/(
                            rhoArgD[dir] - 2.*uArgD[dir])
                           + std::cos(rhoArgH[dir] + 2.*uArgH[dir])/(
                             rhoArgD[dir] + 2.*uArgD[dir])
                           + 2.*std::sin(rhoArgD[dir]*sLocHigh[dir])
                           *std::sin(m_sPhiRho[dir])/rhoArgD[dir]
                           - 2.*std::cos(rhoArgD[dir]*sLocHigh[dir])
                           *std::cos(m_sPhiRho[dir])/rhoArgD[dir]
                           - std::cos(rhoArgL[dir] - 2.*uArgL[dir])/(
                             rhoArgD[dir] - 2.*uArgD[dir])
                           - std::cos(rhoArgL[dir] + 2.*uArgL[dir])/(
                             rhoArgD[dir] + 2.*uArgD[dir])
                           - 2.*std::sin(rhoArgD[dir]*sLocLow[dir])
                           *std::sin(m_sPhiRho[dir])/rhoArgD[dir]
                           + 2.*std::cos(rhoArgD[dir]*sLocLow[dir])
                           *std::cos(m_sPhiRho[dir])/rhoArgD[dir])/4;
                      }
                  }
              }
            Real uu_5 = m_rA*m_uA[0]*(m_uA[0]/vol)
            *std::sin(uTimeArg)*std::sin(uTimeArg)*std::sin(rhoTimeArg)
            *(rhoUU_integralDir.product());
            Real rhoUU_cellAvg = uu_0 + uu_1 + uu_2 + uu_3 + uu_4 + uu_5;,

            // Rho*v*v term
            Real vv_0 = rM*velM[1]*velM[1]; // mean rho*v*v component
            // Mean v*v*rho_sinusoid
            Real vv_1 = velM[1]*velM[1]*(m_rA/vol)
            *(std::sin(rhoTimeArg))*rho_integralDir.product();
            // Mean rho*v*v_sinusoid
            Real vv_2 = 2.*velM[1]*rM*(m_uA[1]/vol)*(std::sin(vTimeArg))
            *v_integralDir.product();
            // Mean v*rho_sinusoid*v_sinusoid
            Real vv_3 = 2.*velM[1]*m_rA*(m_uA[1]/vol)
            *(std::sin(rhoTimeArg))*(std::sin(vTimeArg))
            *(rhoV_integralDir.product());
            // Mean rho*v_sinusoid*v_sinusoid
            RealVect VV_integralDir = RealVect::Zero;
            // Note: because this is v*v there are only two options
            // 1) m_sOmegaV[0] = 0
            //    -- this requires just multiplying the two sinusoids
            // 2) m_sOmegaV[0] != 0
            //    -- this requires integration using same frequencies
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if (vArgD[dir] == 0)
                  {
                    VV_integralDir[dir] = std::sin(m_sPhiV[dir])
                      *std::sin(m_sPhiV[dir])*length[dir];
                  }
                else
                  {
                    VV_integralDir[dir] = (1./(4.*vArgD[dir]))*(
                      2.*vArgH[dir] - std::sin(2.*vArgH[dir])
                      - 2.*vArgL[dir] + std::sin(2.*vArgL[dir]));
                  }
              }
            Real vv_4 = rM*m_uA[1]*(m_uA[1]/vol)
            *(std::sin(vTimeArg))*(std::sin(vTimeArg))
            *(VV_integralDir.product());
            // Rho_sinusoid*v_sinusoid*v_sinusoid
            RealVect rhoVV_integralDir = RealVect::Zero;
            // Note: there are 15 cases here, 10 of them are unnecessary in
            //       our case (r*v*v will always have v = v). Additionally
            //       case 5 has two variants
            //       a) m_sOmegaRho = 2.*m_sOmegaV[0]
            //       b) m_sOmegaRho != 2.*m_sOmegaV[0]
            // 1) m_sOmegaV[0] = 0 and m_sOmegaRho = 0
            //    -- this requires just multiplying the two sinusoids
            // 2) m_sOmegaV[0] = 0 and m_sOmegaRho != 0
            //    -- this requires integrating just the rho sinusoid
            // 3) m_sOmegaV[0] != 0 and m_sOmegaRho = 0
            //    -- this requires integrating just the u sinusoid
            // 4) m_sOmegaV[0] != 0 and m_sOmegaRho != 0 and
            //    m_sOmegaV[0] = m_sOmegaRho[0]
            //    -- this requires integration using same frequencies
            // 5) m_sOmegaV[0] != 0 and m_sOmegaRho != 0 and
            //    m_sOmegaV[0] != m_sOmegaRho[0]
            //    -- this requires integration using different frequencies
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if ((vArgD[dir] == 0) && (rhoArgD[dir] == 0))
                  {
                    rhoVV_integralDir[dir] =
                      std::sin(m_sPhiV[dir])*std::sin(m_sPhiV[dir])
                      *std::sin(m_sPhiRho[dir])*length[dir];
                  }
                else if ((vArgD[dir] == 0) && (rhoArgD[dir] != 0))
                  {
                    rhoVV_integralDir[dir] =
                      std::sin(m_sPhiV[dir])*std::sin(m_sPhiV[dir])
                      *rho_integralDir[dir];
                  }
                else if ((vArgD[dir] != 0) && (rhoArgD[dir] == 0))
                  {
                    rhoVV_integralDir[dir] =
                      VV_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                  }
                else if (vArgD[dir] == rhoArgD[dir])
                  {
                    rhoVV_integralDir[dir] = (1./(12.*vArgD[dir]))
                      *(std::cos(rhoArgH[dir] + vArgH[dir] + vArgH[dir])
                        - 3.*std::cos(rhoArgH[dir] - vArgH[dir] + vArgH[dir])
                        - 3.*std::cos(rhoArgH[dir] + vArgH[dir] - vArgH[dir])
                        - 3.*std::cos(rhoArgH[dir] - vArgH[dir] - vArgH[dir])
                        - std::cos(rhoArgL[dir] + vArgL[dir] + vArgL[dir])
                        + 3.*std::cos(rhoArgL[dir] - vArgL[dir] + vArgL[dir])
                        + 3.*std::cos(rhoArgL[dir] + vArgL[dir] - vArgL[dir])
                        + 3.*std::cos(rhoArgL[dir] - vArgL[dir] - vArgL[dir]));
                  }
                else
                  {
                    if (rhoArgD[dir] == (2.*vArgD[dir]))
                      {
                        rhoVV_integralDir[dir] =
                          (std::cos(rhoArgH[dir] + 2.*vArgH[dir])
                           - 4.*std::cos(rhoArgH[dir])
                           - 2.*rhoArgD[dir]*sLocHigh[dir]
                           *std::sin(rhoArgH[dir] - 2.*vArgH[dir])
                           - std::cos(rhoArgL[dir] + 2.*vArgL[dir])
                           + 4.*std::cos(rhoArgL[dir])
                           + 2.*rhoArgD[dir]*sLocLow[dir]
                           *std::sin(rhoArgL[dir] - 2.*vArgL[dir])
                            )/(16.*vArgD[dir]);
                      }
                    else
                      {
                        rhoVV_integralDir[dir] =
                          (std::cos(rhoArgH[dir] - 2.*vArgH[dir])/(
                            rhoArgD[dir] - 2.*vArgD[dir])
                           + std::cos(rhoArgH[dir] + 2.*vArgH[dir])/(
                             rhoArgD[dir] + 2.*vArgD[dir])
                           + 2.*std::sin(rhoArgD[dir]*sLocHigh[dir])
                           *std::sin(m_sPhiRho[dir])/rhoArgD[dir]
                           - 2.*std::cos(rhoArgD[dir]*sLocHigh[dir])
                           *std::cos(m_sPhiRho[dir])/rhoArgD[dir]
                           - std::cos(rhoArgL[dir] - 2.*vArgL[dir])/(
                             rhoArgD[dir] - 2.*vArgD[dir])
                           - std::cos(rhoArgL[dir] + 2.*vArgL[dir])/(
                             rhoArgD[dir] + 2.*vArgD[dir])
                           - 2.*std::sin(rhoArgD[dir]*sLocLow[dir])
                           *std::sin(m_sPhiRho[dir])/rhoArgD[dir]
                           + 2.*std::cos(rhoArgD[dir]*sLocLow[dir])
                           *std::cos(m_sPhiRho[dir])/rhoArgD[dir])/4;
                      }
                  }
              }
            Real vv_5 = m_rA*m_uA[1]*(m_uA[1]/vol)
            *std::sin(vTimeArg)*std::sin(vTimeArg)*std::sin(rhoTimeArg)
            *(rhoVV_integralDir.product());
            Real rhoVV_cellAvg = vv_0 + vv_1 + vv_2 + vv_3 + vv_4 + vv_5;,

            // Rho*w*w terms
            Real ww_0 = rM*velM[2]*velM[2]; // mean rho*w*w component
            // Mean w*w*rho_sinusoid
            Real ww_1 = velM[2]*velM[2]*(m_rA/vol)
            *(std::sin(rhoTimeArg))*rho_integralDir.product();
            // Mean rho*w*w_sinusoid
            Real ww_2 = 2.*velM[2]*rM*(m_uA[2]/vol)*(std::sin(wTimeArg))
            *w_integralDir.product();
            // Mean w*rho_sinusoid*w_sinusoid
            Real ww_3 = 2.*velM[2]*m_rA*(m_uA[2]/vol)
            *(std::sin(rhoTimeArg))*(std::sin(wTimeArg))
            *(rhoW_integralDir.product());
            // Mean rho*w_sinusoid*w_sinusoid
            RealVect WW_integralDir = RealVect::Zero;
            // Note: because this is w*w there are only two options
            // 1) m_sOmegaW[0] = 0
            //    -- this requires just multiplying the two sinusoids
            // 2) m_sOmegaW[0] != 0
            //    -- this requires integration using same frequencies
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if (wArgD[dir] == 0)
                  {
                    WW_integralDir[dir] = std::sin(m_sPhiW[dir])
                      *std::sin(m_sPhiW[dir])*length[dir];
                  }
                else
                  {
                    WW_integralDir[dir] = (1./(4.*wArgD[dir]))*(
                      2.*wArgH[dir] - std::sin(2.*wArgH[dir])
                      - 2.*wArgL[dir] + std::sin(2.*wArgL[dir]));
                  }
              }
            Real ww_4 = rM*m_uA[2]*(m_uA[2]/vol)
            *(std::sin(wTimeArg))*(std::sin(wTimeArg))
            *(WW_integralDir.product());
            // Rho_sinusoid*w_sinusoid*w_sinusoid
            RealVect rhoWW_integralDir = RealVect::Zero;
            // Note: there are 15 cases here, 10 of them are unnecessary in
            //       our case (r*w*w will always have w = w). Additionally
            //       case 5 has two variants
            //       a) m_sOmegaRho = 2.*m_sOmegaW[0]
            //       b) m_sOmegaRho != 2.*m_sOmegaW[0]
            // 1) m_sOmegaW[0] = 0 and m_sOmegaRho = 0
            //    -- this requires just multiplying the two sinusoids
            // 2) m_sOmegaW[0] = 0 and m_sOmegaRho != 0
            //    -- this requires integrating just the rho sinusoid
            // 3) m_sOmegaW[0] != 0 and m_sOmegaRho = 0
            //    -- this requires integrating just the u sinusoid
            // 4) m_sOmegaW[0] != 0 and m_sOmegaRho != 0 and
            //    m_sOmegaW[0] = m_sOmegaRho[0]
            //    -- this requires integration using same frequencies
            // 5) m_sOmegaW[0] != 0 and m_sOmegaRho != 0 and
            //    m_sOmegaW[0] != m_sOmegaRho[0]
            //    -- this requires integration using different frequencies
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if ((wArgD[dir] == 0) && (rhoArgD[dir] == 0))
                  {
                    rhoWW_integralDir[dir] =
                      std::sin(m_sPhiW[dir])*std::sin(m_sPhiW[dir])
                      *std::sin(m_sPhiRho[dir])*length[dir];
                  }
                else if ((wArgD[dir] == 0) && (rhoArgD[dir] != 0))
                  {
                    rhoWW_integralDir[dir] =
                      std::sin(m_sPhiW[dir])*std::sin(m_sPhiW[dir])
                      *rho_integralDir[dir];
                  }
                else if ((wArgD[dir] != 0) && (rhoArgD[dir] == 0))
                  {
                    rhoWW_integralDir[dir] =
                      WW_integralDir[dir]*std::sin(m_sPhiRho[dir]);
                  }
                else if (wArgD[dir] == rhoArgD[dir])
                  {
                    rhoWW_integralDir[dir] = (1./(12.*wArgD[dir]))
                      *(std::cos(rhoArgH[dir] + wArgH[dir] + wArgH[dir])
                        - 3.*std::cos(rhoArgH[dir] - wArgH[dir] + wArgH[dir])
                        - 3.*std::cos(rhoArgH[dir] + wArgH[dir] - wArgH[dir])
                        - 3.*std::cos(rhoArgH[dir] - wArgH[dir] - wArgH[dir])
                        - std::cos(rhoArgL[dir] + wArgL[dir] + wArgL[dir])
                        + 3.*std::cos(rhoArgL[dir] - wArgL[dir] + wArgL[dir])
                        + 3.*std::cos(rhoArgL[dir] + wArgL[dir] - wArgL[dir])
                        + 3.*std::cos(rhoArgL[dir] - wArgL[dir] - wArgL[dir]));
                  }
                else
                  {
                    if (rhoArgD[dir] == (2.*wArgD[dir]))
                      {
                        rhoWW_integralDir[dir] =
                          (std::cos(rhoArgH[dir] + 2.*wArgH[dir])
                           - 4.*std::cos(rhoArgH[dir])
                           - 2.*rhoArgD[dir]*sLocHigh[dir]
                           *std::sin(rhoArgH[dir] - 2.*wArgH[dir])
                           - std::cos(rhoArgL[dir] + 2.*wArgL[dir])
                           + 4.*std::cos(rhoArgL[dir])
                           + 2.*rhoArgD[dir]*sLocLow[dir]
                           *std::sin(rhoArgL[dir] - 2.*wArgL[dir])
                            )/(16.*wArgD[dir]);
                      }
                    else
                      {
                        rhoWW_integralDir[dir] =
                          (std::cos(rhoArgH[dir] - 2.*wArgH[dir])/(
                            rhoArgD[dir] - 2.*wArgD[dir])
                           + std::cos(rhoArgH[dir] + 2.*wArgH[dir])/(
                             rhoArgD[dir] + 2.*wArgD[dir])
                           + 2.*std::sin(rhoArgD[dir]*sLocHigh[dir])
                           *std::sin(m_sPhiRho[dir])/rhoArgD[dir]
                           - 2.*std::cos(rhoArgD[dir]*sLocHigh[dir])
                           *std::cos(m_sPhiRho[dir])/rhoArgD[dir]
                           - std::cos(rhoArgL[dir] - 2.*wArgL[dir])/(
                             rhoArgD[dir] - 2.*wArgD[dir])
                           - std::cos(rhoArgL[dir] + 2.*wArgL[dir])/(
                             rhoArgD[dir] + 2.*wArgD[dir])
                           - 2.*std::sin(rhoArgD[dir]*sLocLow[dir])
                           *std::sin(m_sPhiRho[dir])/rhoArgD[dir]
                           + 2.*std::cos(rhoArgD[dir]*sLocLow[dir])
                           *std::cos(m_sPhiRho[dir])/rhoArgD[dir])/4;
                      }
                  }
              }
            Real ww_5 = m_rA*m_uA[2]*(m_uA[2]/vol)
            *std::sin(wTimeArg)*std::sin(wTimeArg)*std::sin(rhoTimeArg)
            *(rhoWW_integralDir.product());
            Real rhoWW_cellAvg = ww_0 + ww_1 + ww_2 + ww_3 + ww_4 + ww_5;
            );
          Real ke = 0.5*(D_TERM(rhoUU_cellAvg,
                                + rhoVV_cellAvg,
                                + rhoWW_cellAvg));
          Real energy = p_cellAvg/(gamma - 1.) + ke;
          a_Ux[MD_IX(i, presIndx)] = energy;
        }
    }
  return 0;
}

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
 *//*-----------------------------------------------------------------*/

void
CNSIBCMMS::setImposedBCprimState(
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
  const Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
  const int cVel   = velIntv.begin();
  const int cRho   = CRDparam::g_CRDPhysics->densityIndex();
  const int cPres  = CRDparam::g_CRDPhysics->pressureIndex();
  const int cTemp   = CRDparam::g_CRDPhysics->temperatureIndex();
  const int numPrim = CRDparam::g_CRDPhysics->numPrimitive();

  switch (a_bcInfo.m_type)
    {
    case CRDparam::DomainBCTypeDirichlet:
    {
      // This allows us to test almost any grid by disregarding the exact
      // form of the boundaries

      // Global indices and variables
      const int numWVar  = CRDparam::g_CRDPhysics->numPrimitive();
      // Obtain the mean values and the wave-amplitude values
      const CRDState& mean = CRDState::get(m_idxStateInit);
      const Real rM = mean.density();
      const RealVect velM = mean.velocity();
      const Real pM = mean.pressure();
      const Real TM = mean.temperature();
      const Real R      = CRDparam::g_R;

      // Create a box to fill fill with face-centered data
      const int faceDir = a_bcIdx.m_dir;
      Box facePntBox = a_boundaryFaceBox;
      facePntBox.grow(1); // Grow by 1 in all directions
      facePntBox.grow(faceDir, -1); // Shrink in the face-normal direction

      // Get the physical-space locations at the face centers
      const BlockCoordSys& blockCoordSys =
        *(a_gridMetrics.getCoordSys(a_disjointBox));
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(a_disjointBox);
      FABSTACKTEMP(faceXFab, facePntBox, SpaceDim);
      FABSTACKTEMP(faceXiFab, facePntBox, SpaceDim);
      CRDparam::g_CNSIBC->getFaceCoordinates(
        facePntBox, faceXiFab, faceXFab, faceDir, blockCoordSys);

      // Create Wc to store the primitive-state data
      FABSTACKTEMP(Wc, facePntBox, numWVar);
      // Initialize to a constant state
      for (int comp = 0; comp != numPrim; ++comp)
        {
          Wc.setVal(mean(comp), facePntBox, comp);
        }

      // Fill Wc with location-specific data
      MD_BOXLOOP(facePntBox, i)
        {
          // Physical location
          const RealVect loc(D_DECL(faceXFab[MD_IX(i, 0)],
                                    faceXFab[MD_IX(i, 1)],
                                    faceXFab[MD_IX(i, 2)]));
          // Spatial location normalized by characteristic length
          const RealVect sLoc = loc/m_charLength;
          // Time normalized by characteristic time
          const Real tLoc = a_time/m_charTime;
          // Spatial argument for density field
          const RealVect rhoArg = 2.*PI*m_sOmegaRho*sLoc + m_sPhiRho;
          // Temporal argument for density field
          const Real rhoTimeArg = 2.*PI*m_tOmegaRho*tLoc + m_tPhiRho;
          // Spatial + temporal arguments for velocity fields
          D_TERM(const RealVect uArg = 2.*PI*m_sOmegaU*sLoc + m_sPhiU;
                 const Real uTimeArg = 2.*PI*m_tOmegaU*tLoc + m_tPhiU;,
                 const RealVect vArg = 2.*PI*m_sOmegaV*sLoc + m_sPhiV;
                 const Real vTimeArg = 2.*PI*m_tOmegaV*tLoc + m_tPhiV;,
                 const RealVect wArg = 2.*PI*m_sOmegaW*sLoc + m_sPhiW;
                 const Real wTimeArg = 2.*PI*m_tOmegaW*tLoc + m_tPhiW;);
          // Spatial argument for pressure field
          const RealVect pArg = 2.*PI*m_sOmegaP*sLoc + m_sPhiP;
          // Temporal argument for pressure field
          const Real pTimeArg = 2.*PI*m_tOmegaP*tLoc + m_tPhiP;
          // Spatial argument for temperature field
          const RealVect TArg = 2.*PI*m_sOmegaT*sLoc + m_sPhiT;
          // Temporal argument for temperature field
          const Real TTimeArg = 2.*PI*m_tOmegaT*tLoc + m_tPhiT;
          // Define the density field
          Wc[MD_IX(i, cRho)] = rM + m_rA*(std::sin(rhoTimeArg))*D_TERM(
            std::sin(rhoArg[0]),*std::sin(rhoArg[1]),*std::sin(rhoArg[2]));
          // Define the velocity fields -- first, the X-velocity field
          D_TERM(Wc[MD_IX(i, cVel)]     = velM[0] +
                 m_uA[0]*(std::sin(uTimeArg))*D_TERM(
                   std::sin(uArg[0]),*std::sin(uArg[1]),*std::sin(uArg[2]));,
                 // Y-velocity field
                 Wc[MD_IX(i, cVel + 1)] = velM[1] +
                 m_uA[1]*(std::sin(vTimeArg))*D_TERM(
                   std::sin(vArg[0]),*std::sin(vArg[1]),*std::sin(vArg[2]));,
                 // Z-velocity field
                 Wc[MD_IX(i, cVel + 2)] = velM[2] +
                 m_uA[2]*(std::sin(wTimeArg))*D_TERM(
                   std::sin(wArg[0]),*std::sin(wArg[1]),*std::sin(wArg[2])););
          if (!m_bcTestCase)
            {
              // Define the pressure field
              Wc[MD_IX(i, cPres)] = pM + m_pA*(std::sin(pTimeArg))*D_TERM(
                std::sin(pArg[0]),*std::sin(pArg[1]),*std::sin(pArg[2]));

              // Define the temperature field from the rest of the data
              Wc[MD_IX(i, cTemp)] = Wc[MD_IX(i, cPres)]/(R*Wc[MD_IX(i, cRho)]);
            }
          else
            {
              // Define the temperature field
              Wc[MD_IX(i, cTemp)] = TM + m_TA*(std::sin(TTimeArg))*D_TERM(
                std::sin(TArg[0]),*std::sin(TArg[1]),*std::sin(TArg[2]));

              // Define the pressure field from the rest of the data
              Wc[MD_IX(i, cPres)] = R*Wc[MD_IX(i, cRho)]*Wc[MD_IX(i, cTemp)];
            }
        }
      // Convolve the primitive state to get the face-averaged values
      int order = 4;
      FABSTACKTEMP(WcAvg, a_boundaryFaceBox, numWVar);
      CRDutil::convolveFace(WcAvg, Wc, a_boundaryFaceBox, blockDomain,
                            Interval(0,numWVar-1), faceDir, order,
                            false, false, false);
      // Fill a_Wface with WcAvg
      for (int comp = 0; comp != numWVar; ++comp)
        {
          MD_BOXLOOP(a_boundaryFaceBox, i)
            {
              a_Wface[MD_IX(i, comp)] = WcAvg[MD_IX(i, comp)];
            }
        }
      break;
    }
    case CRDparam::DomainBCTypeRelaxedCBCFar:
    case CRDparam::DomainBCTypeFarfield:
    {
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      for (int comp = 0; comp != numPrim; ++comp)
        {
          a_Wface.setVal(state(comp), a_boundaryFaceBox, comp);
        }
      break;
    }
    case CRDparam::DomainBCTypeRelaxedCBCOut:
    {
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      const Real alpha = a_bcInfo.m_relaxCBCStateParam; // 1/0 = non/reflecting
      for (int c = 0; c != numPrim; ++c)
        {
          MD_BOXLOOP(a_boundaryFaceBox, i)
            {
              const Real interior = a_Wcell[MD_IX(i, c)];
              const Real exterior = state(c);
              a_Wface[MD_IX(i, c)] = alpha*interior + (1. - alpha)*exterior;
            }
        }
      break;
    }
    case  CRDparam::DomainBCTypeOutflow:
    {
      // Transform velocity into normal and tangent components
      FABSTACKTEMP(nttVel, a_boundaryFaceBox, velIntv.size());
      nttVel.copy(a_Wface, a_boundaryFaceBox, cVel,
                  a_boundaryFaceBox, 0, velIntv.size());
      PatchMappedFunc::forwardTransform(nttVel,
                                        a_unitNormalBasisFab,
                                        a_boundaryFaceBox);
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          //**FIXME multispecies
          const Real gamma = CRDparam::g_CRDPhysics->gamma();
          //**FIXME Probably need to replace use of cp with h
          const Real cp    = CRDparam::g_CRDPhysics->cp();
          const Real Rgas  = CRDparam::g_CRDPhysics->Rgas();
          // Interior estimate of velocity
          RealVect intVel;
          for (const int idxVel : EachDir)
            {
              intVel[idxVel] = a_Wface[MD_IX(i, cVel + idxVel)];
            }
          const RealVect intVel1(intVel - state.m_frameVelocity);
          const Real intVel1Sq = stc::dot(intVel1, intVel1);
          Real intT;
          // Temperature state based on characteristics
          if (nttVel[MD_IX(i, 0)]*Side::sign(a_bcIdx.m_side) < 0.)  // Inflow
            {
              intT = state.temperature() - intVel1Sq/(2*cp);
            }
          else                                                      // Outflow
            {
              intT = a_Wface[MD_IX(i, cTemp)];
            }
          const Real intc1 = 1 + intVel1Sq/(2*cp*intT);
          // Whether inflow or outflow, this is the assigned pressure.  Here,
          // assume a frozen mixture with gamma from the interior.
          const Real p = state.pressure()*std::pow(intc1, -gamma/(gamma - 1.0));
          a_Wface[MD_IX(i, cPres)] = p;
          // Temperature and density assigned only if inflow
          if (nttVel[MD_IX(i, 0)]*Side::sign(a_bcIdx.m_side) < 0.)  // Inflow
            {
              const Real T = state.temperature()/intc1;
              a_Wface[MD_IX(i, cRho)]  = p/(Rgas*T);
              a_Wface[MD_IX(i, cTemp)] = T;
            }
        }
      break;
    }
    case  CRDparam::DomainBCTypeRelaxedCBCIn:
    case  CRDparam::DomainBCTypeInflow:
    {
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          //**FIXME multispecies
          const Real gamma = CRDparam::g_CRDPhysics->gamma();
          //**FIXME Probably need to replace use of cp with h
          const Real cp    = CRDparam::g_CRDPhysics->cp();
          const Real Rgas  = CRDparam::g_CRDPhysics->Rgas();
          const RealVect extVel = state.velocity();
          RealVect intVel;
          for (const int idxVel : EachDir)
            {
              intVel[idxVel] = a_Wface[MD_IX(i, cVel + idxVel)];
              a_Wface[MD_IX(i, cVel + idxVel)] = extVel[idxVel];
            }
          if (CRDparam::DomainBCTypeRelaxedCBCIn & a_bcInfo.m_type)
            {
              for (int comp = 0; comp != numPrim; ++comp)
                {
                  a_Wface[MD_IX(i, comp)] = state(comp);
                }
            }
          else
            {
              // Pressure is determined from the interior.  First, find a
              // stagnation pressure for the given reference frame.  Assume
              // gas is frozen at upstream conditions.
              const RealVect intVel1(intVel - state.m_frameVelocity);
              const Real intVel1Sq = stc::dot(intVel1, intVel1);
              const Real intT = a_Wface[MD_IX(i, cTemp)];
              const Real intc1 = 1 + intVel1Sq/(2*cp*intT);
              const Real p0 = a_Wface[MD_IX(i, cPres)]*
                std::pow(intc1, gamma/(gamma - 1.));
              // Find imposed pressure from stagnation pressure
              const RealVect extVel1(extVel - state.m_frameVelocity);
              const Real extVel1Sq = stc::dot(extVel1, extVel1);
              const Real extc1 = 1 + extVel1Sq/(2*cp*state.temperature());
              a_Wface[MD_IX(i, cPres)] =
                p0*std::pow(extc1, -gamma/(gamma - 1.0));
              // Inflow cannot become outflow, so we always set density and
              // temperature
              a_Wface[MD_IX(i, cRho)] = state.density();
              a_Wface[MD_IX(i, cTemp)] = a_Wface[MD_IX(i, cPres)]/
                (state.density()*Rgas);
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
CNSIBCMMS::readBCInfo()
{
  ParmParse ppIBC("ibc");

  // Read density-wave amplitude
  ppIBC.query("wave_amp_density", m_rA);
  // Read velocity-wave amplitude
  std::vector<Real> velA(SpaceDim);
  ppIBC.queryarr("wave_amp_velocity", velA, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_uA.dataPtr(), &velA.front());
  // Read pressure-wave amplitude
  ppIBC.query("wave_amp_pressure", m_pA);
  // Read temperature-wave amplitude
  ppIBC.query("wave_amp_temperature", m_TA);

  // Read density spatial-frequency
  std::vector<Real> sOmegaRho(SpaceDim);
  ppIBC.queryarr("space_omega_density", sOmegaRho, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_sOmegaRho.dataPtr(),
                                           &sOmegaRho.front());
  // Read density temporal-frequency
  ppIBC.query("time_omega_density", m_tOmegaRho);
  // Read density spatial phase-shift
  std::vector<Real> sPhiRho(SpaceDim);
  ppIBC.queryarr("space_phi_density", sPhiRho, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_sPhiRho.dataPtr(),
                                           &sPhiRho.front());
  // Read density temporal phase-shift
  ppIBC.query("time_phi_density", m_tPhiRho);

  // Read x-velocity spatial-frequency
  std::vector<Real> sOmegaU(SpaceDim);
  ppIBC.queryarr("space_omega_x-velocity", sOmegaU, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_sOmegaU.dataPtr(),
                                           &sOmegaU.front());
  // Read x-velocity temporal-frequency
  ppIBC.query("time_omega_x-velocity", m_tOmegaU);
  // Read x-velocity spatial phase-shift
  std::vector<Real> sPhiU(SpaceDim);
  ppIBC.queryarr("space_phi_x-velocity", sPhiU, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_sPhiU.dataPtr(), &sPhiU.front());
  // Read x-velocity temporal phase-shift
  ppIBC.query("time_phi_x-velocity", m_tPhiU);

  // Read y-velocity spatial-frequency
  std::vector<Real> sOmegaV(SpaceDim);
  ppIBC.queryarr("space_omega_y-velocity", sOmegaV, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_sOmegaV.dataPtr(),
                                           &sOmegaV.front());
  // Read y-velocity temporal-frequency
  ppIBC.query("time_omega_y-velocity", m_tOmegaV);
  // Read y-velocity spatial phase-shift
  std::vector<Real> sPhiV(SpaceDim);
  ppIBC.queryarr("space_phi_y-velocity", sPhiV, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_sPhiV.dataPtr(), &sPhiV.front());
  // Read y-velocity temporal phase-shift
  ppIBC.query("time_phi_y-velocity", m_tPhiV);

  // Read z-velocity spatial-frequency
  std::vector<Real> sOmegaW(SpaceDim);
  ppIBC.queryarr("space_omega_z-velocity", sOmegaW, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_sOmegaW.dataPtr(),
                                           &sOmegaW.front());
  // Read z-velocity temporal-frequency
  ppIBC.query("time_omega_z-velocity", m_tOmegaW);
  // Read z-velocity spatial phase-shift
  std::vector<Real> sPhiW(SpaceDim);
  ppIBC.queryarr("space_phi_z-velocity", sPhiW, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_sPhiW.dataPtr(), &sPhiW.front());
  // Read z-velocity temporal phase-shift
  ppIBC.query("time_phi_z-velocity", m_tPhiW);

  // Read pressure spatial-frequency
  std::vector<Real> sOmegaP(SpaceDim);
  ppIBC.queryarr("space_omega_pressure", sOmegaP, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_sOmegaP.dataPtr(),
                                           &sOmegaP.front());
  // Read pressure temporal-frequency
  ppIBC.query("time_omega_pressure", m_tOmegaP);
  // Read pressure spatial phase-shift
  std::vector<Real> sPhiP(SpaceDim);
  ppIBC.queryarr("space_phi_pressure", sPhiP, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_sPhiP.dataPtr(), &sPhiP.front());
  // Read pressure temporal phase-shift
  ppIBC.query("time_phi_pressure", m_tPhiP);

  // Read temperature spatial-frequency
  std::vector<Real> sOmegaT(SpaceDim);
  ppIBC.queryarr("space_omega_temperature", sOmegaT, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_sOmegaT.dataPtr(),
                                           &sOmegaT.front());
  // Read temperature temporal-frequency
  ppIBC.query("time_omega_temperature", m_tOmegaT);
  // Read temperature spatial phase-shift
  std::vector<Real> sPhiT(SpaceDim);
  ppIBC.queryarr("space_phi_temperature", sPhiT, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_sPhiT.dataPtr(), &sPhiT.front());
  // Read temperature temporal phase-shift
  ppIBC.query("time_phi_temperature", m_tPhiT);

  // Read characteristic length
  std::vector<Real> charL(SpaceDim);
  ppIBC.queryarr("characteristic_length", charL, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_charLength.dataPtr(),
                                           &charL.front());
  // Read characteristic time
  ppIBC.query("characteristic_time", m_charTime);
  // Read case number
  ppIBC.query("bc_test_case", m_bcTestCase);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  //                                                      //
  //   Set up some variables based on specific BC cases   //
  //                                                      //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++ //

  const CRDState& mean = CRDState::get(m_idxStateInit);
  const RealVect velM = mean.velocity();
  // BC case 1:
  // Left boundary: periodic
  // Right boundary: periodic
  // Top boundary: adiabatic no-slip wall
  //   a) all velocities must be zero
  //   b) wall-normal temperature gradient must be zero
  // Bottom boundary: adiabatic no-slip wall
  //   a) all velocities must be zero
  //   b) wall-normal temperature gradient must be zero
  // Mean x-velocity: 0
  // Mean y-velocity: 0
  // Mean z-velocity: 0
  if (m_bcTestCase == 1)
    {
      // Enforce zero x-velocity at the top and bottom walls
      m_sPhiU[1] = 0.;
      // Enforce integer wave-number in the y-direction
      Real fraction = std::modf(m_sOmegaU[1], &m_sOmegaU[1]);
      // Enforce zero y-velocity at the top and bottom walls
      m_sPhiV[1] = 0.;
      // Enforce integer wave-number in the y-direction
      fraction = std::modf(m_sOmegaV[1], &m_sOmegaV[1]);
      // Enforce zero temperature gradient at the top and bottom walls
      m_sPhiT[1] = 0.;
      // Enforce integer wave-number in the y-direction
      fraction = std::modf(m_sOmegaT[1], &m_sOmegaT[1]);
      // Remove annoying warning
      fraction = fraction;
    }
  // BC case 2:
  // Left boundary: inlet
  //   a) wall-normal velocity must be into domain
  // Right boundary: outlet
  //   a) anything is allowed
  // Top boundary: adiabatic no-slip wall
  //   a) all velocities must be zero
  //   b) wall-normal temperature gradient must be zero
  // Bottom boundary: isothermal no-slip wall
  //   a) all velocities must be zero
  //   b) wall-normal temperature gradient can be anything
  else if (m_bcTestCase == 2)
    {
      // Inlet-normal velocity must have a non-zero mean
    }
  ppIBC.query("exact_initialization", m_exactInit);
  ppIBC.query("sixth_order_initialization", m_sixthOrderInit);

  m_readInput = true;
}
