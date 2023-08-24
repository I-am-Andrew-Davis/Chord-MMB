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
 * \file CNSIBCBluffBodyCombustion.cpp
 *
 * \brief Member functions for CNSIBCBluffBodyCombustion.H
 *
 *//*+*************************************************************************/

#include <algorithm>

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
#include "LoHiCenter.H"
 
//----- Internal -----//

#include "CNSIBCBluffBodyCombustion.H"
#include "CNSIBCCombustionTestF_F.H"
#include "ViscousTensor4thOrderOpF_F.H"
#include "CNSIBCF_F.H"
#include "ChordInput.H"
#include "CRDPhysics.H"
#include "CRDState.H"
#include "PatchMappedFunc.H"

/*******************************************************************************
 *
 * Class CNSIBCBluffBodyCombustion: member definitions
 *
 ******************************************************************************/

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCBluffBodyCombustion::CNSIBCBluffBodyCombustion(const bool a_readInputs)
  :
  CNSIBCGeneralized(),
  m_idxStateRegion2(-1),
  m_region2_lo(RealVect::Zero),
  m_region2_hi(RealVect::Zero)
{
  // read in all other info
  readBCInfo();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCBluffBodyCombustion::~CNSIBCBluffBodyCombustion()
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
CNSIBCBluffBodyCombustion::IBCName() const
{
  return "Bluff-body combustion case";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/
void
CNSIBCBluffBodyCombustion::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
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
CNSIBCBluffBodyCombustion::initialize(LevelData<FArrayBox>&      a_U,
                                      LevelGridMetrics&          a_gridMetrics,
                                      const LayoutData<FluxBox>& a_unitNormals,
                                      const Real                 a_time,
                                      const int                  a_level) const
{
  // Set the initial values
  const int numWVar      = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx     = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx  = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx     = CRDparam::g_CRDPhysics->temperatureIndex();
  const int speciesComp  = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();

  const int numSpecies = CRDparam::g_numSpecies;
  const CRDState& state = CRDState::get(m_idxStateInit);
  const CRDState& stateRegion2 = CRDState::get(m_idxStateRegion2);
  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));
      // Get the domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Get physical coordinates
      FABSTACKTEMP(XiFab, box2Dom, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, box2Dom, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(box2Dom, XiFab, XFab, blockCoordSys);

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(0.);
      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box2Dom));
      FABSTACKTEMP(Wc, box2Dom, numWVar);
      Wc.setVal(0.);
      
      //Wc.setVal(state.density(), rhoIndx);
      const Real tempRho = -1.0;
      Wc.setVal(tempRho, rhoIndx);
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          Wc.setVal(state.velocity()[dir], box2Dom, WvelIndx + dir, 1);
        }
      Wc.setVal(state.pressure(), presIndx);
      Wc.setVal(state.temperature(), tempIndx);
      for (int i = 0; i != numSpecies; ++i)
        {
          int comp = speciesComp + i;
          Wc.setVal(state(comp), comp);
        }

      //**FIXME: NEED TO ADD THE SMOOTHING ACROSSING THE INITIAL INTERFACE
      MD_BOXLOOP(box2Dom, i)
        {
          // Region-2 spcification
          const RealVect loc(D_DECL(XFab[MD_IX(i, 0)], XFab[MD_IX(i, 1)],
                                    XFab[MD_IX(i, 2)]));
          
          if (D_TERM(loc[0] >= m_region2_lo[0] && loc[0] <= m_region2_hi[0],
                     && loc[1] >= m_region2_lo[1] && loc[1] <= m_region2_hi[1],
                     && loc[2] >= m_region2_lo[2] && loc[2] <= m_region2_hi[2])
                )
            {
              Wc[MD_IX(i, tempIndx)] = stateRegion2.temperature();
              for (int s = 0; s != numSpecies; ++s)
                {
                  int comp = speciesComp + s;
                  Wc[MD_IX(i, comp)] = stateRegion2(comp);
                }
            }
        }
      
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box2Dom);
      fourthOrderAverageCell(UFab, blockDomain, box1Dom);
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
 *  \param[in]  a_globalKE
 *                      Global mean kinetic energy (turbulent forcing)
 *  \param[in]  a_globalHelicity
 *                      Global mean helicity (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void
CNSIBCBluffBodyCombustion::addSourceTerm(
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
  const Real           a_globalKE,
  const Real           a_globalHelicity) const
{
  const int cTemp = CRDparam::g_CRDPhysics->temperatureIndex();
  const int cEnergy = CRDparam::g_CRDPhysics->energyFluxIndex();
  // Add a heat-sink sponge-region near the walls to control temperature
  // fluctuations due to near-wall numerical error

  // We need to set a temperature threshold for the source term
  // Above the threshold, the source term will lower the temperature

  const Real tempThreshold = 1900.; // Needs to be user-defined

  // Source terms are added to the solution, so that requires
  //
  // dS/dt = ((rho*e)_target - (rho*e)_current)/dt
  //
  // The key is, we don't want this update to perfectly match
  // the target energy value in a single step as this would likely
  // be numerically unstable. So, we have to provide a fractional
  // scaling value.

  const Real f_0 = 0.05; // Reach the desired target in 20 time-steps
  const Real dt = 1.e-6; // Update this function to provide current time-step
                         // Until then, it's safer to use a larger number

  // Get the domain for the block for this level
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_disjointBox);

  int isYLoWall = 0; // Default to not being a no-slip wall boundary
  int isYHiWall = 0;

  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  // Get physical coordinates
  FABSTACKTEMP(XiFab, a_solveBox, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, a_solveBox, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(a_solveBox, XiFab, XFab, blockCoordSys);
  Box boundaryFaceBoxLo;
  CRDparam::g_CNSIBC->getBoundaryFaces(
        boundaryFaceBoxLo, a_solveBox, blockDomain, 1, Side::Lo);
  Box boundaryFaceBoxHi;
  CRDparam::g_CNSIBC->getBoundaryFaces(
        boundaryFaceBoxHi, a_solveBox, blockDomain, 1, Side::Hi);

  if (!boundaryFaceBoxHi.isEmpty())
    {
      // Check to make sure it's not the bluff body
      Real maxY = -999;
      MD_BOXLOOP(a_solveBox, i)
        {
          RealVect xloc(D_DECL(XFab[MD_IX(i,0)],
                               XFab[MD_IX(i,1)],
                               XFab[MD_IX(i,2)]));
          maxY=std::max(xloc[1],maxY);
        }
      if (maxY >= 0.01)
        {
          isYHiWall = 1;
        }
    }
  if (!boundaryFaceBoxLo.isEmpty())
    {
      // Check to make sure it's not the bluff body
      Real minY = 999;
      MD_BOXLOOP(a_solveBox, i)
        {
          RealVect xloc(D_DECL(XFab[MD_IX(i,0)],
                               XFab[MD_IX(i,1)],
                               XFab[MD_IX(i,2)]));
          minY=std::min(xloc[1],minY);

        }
      if (minY <= -0.01)
        {
          isYLoWall = 1;
        }
    }
  // Get the number of cells in the spanwise and Y direction
  // Only want to do this in 3D 
  D_TERM(,
         const int numCellsY = blockDomain.domainBox().size()[1];
         const int hiYIndx = numCellsY - 1;,
         const int numCellsZ = blockDomain.domainBox().size()[2];
         const int hiZIndx = numCellsZ - 1;);

  // Create the target primitive state
  FABSTACKTEMP(Wtarget, a_solveBox, a_Wcell.nComp());
  Wtarget.copy(a_Wcell);
  Wtarget.setVal(tempThreshold, cTemp);
  // Compute the target conservative state
  FABSTACKTEMP(Utarget, a_solveBox, CRDparam::g_CRDPhysics->numConservative());
  FABSTACKTEMP(Ucurrent, a_solveBox, CRDparam::g_CRDPhysics->numConservative());
  // Yes, we're actually using this
  CRDparam::g_CRDPhysics->primToCons(Utarget, Wtarget, a_solveBox);
  CRDparam::g_CRDPhysics->primToCons(Ucurrent, a_Wcell, a_solveBox);

  // Compute the difference and scale the source term
  MD_BOXLOOP(a_solveBox, i)
    {
      if (a_Wcell[MD_IX(i, cTemp)] >= tempThreshold)
        {
          // We also need to restrict this to right near the z and y-walls
          Real zc_0 = 0.;
          Real yc_0 = 0.;
          D_TERM(,
                 // Y Walls
                 const int yIndx = MD_GETIV(i)[1];
                 if (((yIndx == 0) && (isYLoWall == 1)) ||
                     ((yIndx == hiYIndx) && (isYHiWall == 1)))
                   {
                     yc_0 = 1.;
                   }
                 else if (((yIndx == 1) && (isYLoWall == 1)) ||
                          ((yIndx == (hiYIndx - 1)) && (isYHiWall == 1)))
                   {
                     yc_0 = 0.5;
                   }
                 else if (((yIndx == 2) && (isYLoWall == 1)) ||
                          ((yIndx == (hiYIndx - 2)) && (isYHiWall == 1)))
                   {
                     yc_0 = 0.25;
                   },
                 // Z Walls
                 const int zIndx = MD_GETIV(i)[2];
                 if ((zIndx == 0) || (zIndx == hiZIndx))
                   {
                     zc_0 = 1.;
                   }
                 else if ((zIndx == 1) || (zIndx == (hiZIndx - 1)))
                   {
                     zc_0 = 0.5;
                   }
                 else if ((zIndx == 2) || (zIndx == (hiZIndx - 2)))
                   {
                     zc_0 = 0.25;
                   }
            );
          const Real rhoETarget = Utarget[MD_IX(i, cEnergy)];
          const Real rhoECurrent = Ucurrent[MD_IX(i, cEnergy)];
          // Take the larger component
          Real c_0 = std::max(zc_0, yc_0);
          a_sourceFab[MD_IX(i, cEnergy)] = 
            c_0 * (f_0/dt) * (rhoETarget - rhoECurrent);
        }
    }
  return;
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
 *//*-----------------------------------------------------------------*/

void
CNSIBCBluffBodyCombustion::setImposedBCprimState(
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
  const int cTemp      = CRDparam::g_CRDPhysics->temperatureIndex();
  const int cSpecies   = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;

  switch (a_bcInfo.m_type)
    {
    case CRDparam::DomainBCTypeDirichlet:
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
      std::vector<Real> cn(CRDparam::g_numSpecies, 0.);
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          for (int comp = 0; comp != CRDparam::g_numSpecies; ++comp)
            {
              const int wComp = comp + cSpecies;
              cn[comp] = a_Wcell[MD_IX(i, wComp)];
            }
          const Real cp = CRDparam::g_CRDPhysics->cp(a_Wface[MD_IX(i, cTemp)],cn.data());
          const Real gamma = CRDparam::g_CRDPhysics->gamma(a_Wface[MD_IX(i, cTemp)],cn.data());
          const Real Rgas  = CRDparam::g_CRDPhysics->Rgas(cn.data());
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
              for (int j = 0; j != numSpecies; ++j)
                {
                  a_Wface[MD_IX(i, j+cSpecies)]  = state(j+cSpecies);
                }
              a_Wface[MD_IX(i, cRho)]  = p/(Rgas*T);
              a_Wface[MD_IX(i, cTemp)] = T;
            }
        }
      break;
    }
    case  CRDparam::DomainBCTypeInflow:
    {
      const CRDState& state = CRDState::get(a_bcInfo.m_idxState);
      std::vector<Real> cn(CRDparam::g_numSpecies, 0.);   
      MD_BOXLOOP(a_boundaryFaceBox, i)
        {
          //**FIXME Probably need to replace use of cp with h
          // const Real cp    = CRDparam::g_CRDPhysics->cp();
          for (int comp = 0; comp != CRDparam::g_numSpecies; ++comp)
            {
              const int wComp = comp + cSpecies;
              cn[comp] = state(wComp);
            }
          const Real cp = CRDparam::g_CRDPhysics->cp(state.temperature(),
                                                     cn.data());
          const Real gamma = CRDparam::g_CRDPhysics->gamma(state.temperature(),
                                                           cn.data());
          const Real Rgas  = CRDparam::g_CRDPhysics->Rgas(cn.data());
          const RealVect extVel = state.velocity();
          RealVect intVel;
          for (const int idxVel : EachDir)
            {
              intVel[idxVel] = a_Wface[MD_IX(i, cVel + idxVel)];
              a_Wface[MD_IX(i, cVel + idxVel)] = extVel[idxVel];
            }
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
          a_Wface[MD_IX(i, cPres)] = p0*std::pow(extc1, -gamma/(gamma - 1.0));
          // Inflow cannot become outflow, so we always set density and
          // temperature
          for (int j = 0; j != numSpecies; ++j)
            {
              a_Wface[MD_IX(i, j+cSpecies)]  = state(j+cSpecies);
            }
          a_Wface[MD_IX(i, cRho)] = state.density();
          a_Wface[MD_IX(i, cTemp)] = a_Wface[MD_IX(i, cPres)]/
            (state.density()*Rgas);
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
CNSIBCBluffBodyCombustion::readBCInfo()
{

// Most IBC info is already read in from generalized
//-- Read the initial state
  ParmParse ppIBC("ibc");

  // State for intitialization
  const char* initStateNameHdr = "initial_state_region2";
  if (ppIBC.contains(initStateNameHdr))
    {
      std::string initStateName;
      ppIBC.get(initStateNameHdr, initStateName);
      m_idxStateRegion2 = CRDState::nameIndex(initStateName);
    }
  else
    {
      CRD::msg << "Input (Bluff-body combustion IBC): " << initStateNameHdr
               << " must be given" << CRD::error;
    }

  std::vector<Real> inputRegion2Lo(SpaceDim, 0.);
  ppIBC.queryarr("init_region2_lo", inputRegion2Lo, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_region2_lo.dataPtr(),
                                           &inputRegion2Lo.front());
  std::vector<Real> inputRegion2Hi(SpaceDim, 0.);
  ppIBC.queryarr("init_region2_hi", inputRegion2Hi, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_region2_hi.dataPtr(),
                                           &inputRegion2Hi.front());
}








