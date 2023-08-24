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
 * \file CNSIBCVortex.cpp
 *
 * \brief Member functions for CNSIBCTaylorGreen
 * 
 *        Here's what this IBC class does. First, you can run the standard
 *        Taylor-Green vortex case (simple and straightforward). Or, you can
 *        run a forced turbulence case. The case can then be restarted
 *        without forcing and the case becomes decaying turbulence. To run
 *        the forcing case, there are a few things to consider. First, what
 *        type of forcing do you want? Right now, there's linear, velocity
 *        forcing. Soon, there will be linear, vorticity forcing. There may
 *        also be some nonlinear varieties soon. Next, choose two parameters
 *        for the simulation.
 *
 *        (1) First, note that the integral length scale of the simulation
 *            is essentially fixed at 1/5 the domain size. Most researchers
 *            report that this is the consistent integral length scale for
 *            linearly forced turbulence simulations such as these.
 *
 *        (2) Choose a target Taylor-microscale Reynolds number.
 *
 *        (3) Choose a target turbulent velocity to match a target turbulent
 *            Mach number. This velocity is the variance of the velocity field.
 *            As such, turbulent kinetic energy is a priori estimated to be 
 *
 *            ke_turb = (3/2)*u_turb^2
 *
 *        (4) Compute the viscosity.
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

#include "CNSIBCTaylorGreen.H"
#include "ChordInput.H"
#include "CRDPhysics.H"
#include "PatchMappedFunc.H"
#include "CRDparam.H"
#include "CRDState.H"
#include "DataTemp.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"
#include "TagMethodValue.H"
#include "TagMethodVMS.H"
#include "CRDutil.H"

/*******************************************************************************
 *
 * Class CNSIBCTaylorGreen: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCTaylorGreen::CNSIBCTaylorGreen()
  :
  CNSIBCGeneralized(),
  m_turbForcing(0),
  m_numVort(2*RealVect::Unit),
  m_zVelIC(0),
  m_coeffG(0.)
{
  CNSIBCTaylorGreen::readBCInfo();
  if (SpaceDim < 3)
    {
      CRD::msg << "Problem not defined for less than 3D!" << CRD::abort;
    }
  
  // BC are all periodic by default  
  // setAllDomainBC(CRDparam::DomainBCTypePeriodic);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCTaylorGreen::~CNSIBCTaylorGreen()
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
CNSIBCTaylorGreen::IBCName() const
{
  return "Taylor-Green Vortex case";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCTaylorGreen::writeIBCInfo() const
{
  const CRDState& turbState = CRDState::get(m_turbState);
  const CRDState& meanState = CRDState::get(m_idxStateInit);
  const Real gamma = CRDparam::g_gamma;
  const Real domLen = CRDparam::g_domainLength[0];
  const Real intLen = domLen/5.;
  const Real tgLen = domLen/(2.*PI);
  const Real mu = CRDparam::g_mu;
  const RealVect meanVel = meanState.velocity();
  const RealVect turbVel = turbState.velocity();
  const Real meanRho = meanState.density();
  const Real meanPress = meanState.pressure();
  const Real M_t = (turbVel[0])/(std::sqrt(gamma*meanPress/meanRho));

  const Real Re_taylor = std::sqrt(15.*meanRho*turbVel[0]*intLen/mu);
  const Real Re_t = meanRho*turbVel[0]*intLen/mu;
  const Real Re_s = meanRho*turbVel[0]*tgLen/mu;
  const Real epsilon = turbVel[0]*turbVel[0]*turbVel[0]/intLen;
  const Real B = turbVel[0]/(intLen*2.924);
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Domain-mean density\n"        << meanRho   << CRD::var;
  CRD::msg << "Domain-mean velocity\n"       << meanVel   << CRD::var;
  CRD::msg << "Turbulent velocity (sqrt(u*u*1/3), RMS velocity)\n"
                                             << turbVel   << CRD::var;
  CRD::msg << "Domain-mean pressure\n"       << meanPress << CRD::var;
  CRD::msg << "Turbulent Mach number\n"      << M_t       << CRD::var;
  CRD::msg << "Taylor Reynolds number\n"     << Re_taylor << CRD::var;
  CRD::msg << "Turbulent Reynolds number\n"  << Re_t      << CRD::var;
  CRD::msg << "Standard Reynolds number\n"   << Re_s      << CRD::var;
  CRD::msg << "Number of vortices\n"         << m_numVort << CRD::var;
  if (CRDparam::g_useTurbForce)
    {
      CRD::msg << "Turbulent dissipation rate\n" << epsilon << CRD::var;
      CRD::msg << "Forcing coefficient (turb dissipation/KE)\n"
                                                 << B << CRD::var;
    }
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
CNSIBCTaylorGreen::initialize(LevelData<FArrayBox>&      a_U,
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

  // Turbulent and mean initial states
  const CRDState& turbState = CRDState::get(m_turbState);
  const CRDState& meanState = CRDState::get(m_idxStateInit);

  // IC constants
  const Real Len = CRDparam::g_domainLength[0];

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
      MD_BOXLOOP(box2Dom, i)
        {
          RealVect loc(D_DECL(XFab[MD_IX(i,0)],
                              XFab[MD_IX(i,1)],
                              XFab[MD_IX(i,2)]));
          RealVect intLoc = m_numVort*loc*PI/Len;
          RealVect turbVel = turbState.velocity();
          if (m_turbForcing)
            {
              turbVel *= std::sqrt(12.);
            }
          const RealVect meanVel = meanState.velocity();
          Real OrigVel = turbVel[0]*D_TERM(std::sin(intLoc[0]),
                                          *std::cos(intLoc[1]),
                                          *std::sin(intLoc[2]));
          Real yVel = -turbVel[0]*D_TERM(std::cos(intLoc[0]),
                                        *std::sin(intLoc[1]),
                                        *std::sin(intLoc[2]));
          if (m_zVelIC)
            {
              OrigVel = turbVel[0]*D_TERM(std::sin(intLoc[0]),
                                         *std::cos(intLoc[1]),
                                         *std::cos(intLoc[2]));
              yVel = -turbVel[0]*D_TERM(std::cos(intLoc[0]),
                                       *std::sin(intLoc[1]),
                                       *std::cos(intLoc[2]));
            }

          D_TERM(Wc[MD_IX(i,WvelIndx)] = meanVel[0] + OrigVel;,
                 Wc[MD_IX(i,WvelIndx+1)] = yVel;,
                 Wc[MD_IX(i,WvelIndx+2)] = 0.;);
          const Real meanPress = meanState.pressure();
          const Real meanRho = meanState.density();
          Real numer = D_SELECT(1,1,std::cos(2.*intLoc[2]) + 2.);
          Real pres = meanPress + meanRho*turbVel[0]*turbVel[0]*
            (D_TERM(std::cos(2.*intLoc[0]),+std::cos(2.*intLoc[1]),))*numer/16.;
          Real rho = pres*meanRho/meanPress;
          if (m_uniformRhoIC)
            {
              rho = meanRho;
            }
          Wc[MD_IX(i, presIndx)] = pres;
          Wc[MD_IX(i, rhoIndx)] = rho;
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
CNSIBCTaylorGreen::addSourceTerm(FArrayBox&           a_sourceFab,
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
  // --------------------------------------------------------------- //
  // Forcing parameter from Rosales 2005                             //
  // --------------------------------------------------------------- //
  // Description: The forcing term is formulated as                  //
  // B = (eps + PD)/(rho*u_t*u_t), eps is turbulent dissipation rate //
  // PD is pressure dilatation (approximate to be very small at low  //
  // Mach number), specify eps according to Rosales (u_t^3/intLen)   //
  // Dynamically compute KE, change velocity initialization to       //
  // balance KE properly                                             //
  // --------------------------------------------------------------- //

  // Turbulent and mean initial states
  const CRDState& turbState = CRDState::get(m_turbState);
  const CRDState& meanState = CRDState::get(m_idxStateInit);

  const Real meanRho = meanState.density();
  const RealVect turbVel = turbState.velocity();

  const int cVel =
    CRDparam::g_CRDPhysics->vectorFluxInterval().begin();

  // Rosales (2005) state intLngthScale to be 1/5 of the domain
  Real intLngthScale = (1./(5.))*CRDparam::g_domainLength[0];
  const Real u_t = turbVel[0];
  // Turbulent dissipation is u_t^3/intLength based on dimensional arguments
  const Real eps = u_t*u_t*u_t/intLngthScale;
  // Target kinetic energy
  const Real ke_target = (3./2.)*u_t*u_t*meanRho;
  // Proportional controller coefficient (Bassenne et al. 2016)
  const Real timeScale = intLngthScale/u_t;
  const Real controller = m_coeffG*(a_globalKE - ke_target)/timeScale;

  // Compute forcing coefficient
  Real B = (eps - controller)/(2.*a_globalKE);
  if (a_globalKE < 1e-12)
    {
      // If we have little global KE, we force so that we do have some
      B = turbVel[0]/(intLngthScale*2.924);;
    }

  // Compute forcing coefficient for divergence-free forcing
  if (m_turbForcing == 2)
    {
      Real eps2 = 1.;
      Real c_0 = a_globalHelicity;
      if (a_globalHelicity <= eps2)
        {
          c_0 = eps2;
        }
      B = (eps - controller)/(c_0);
    }

  if (m_turbForcing == 1)
    {
      for (int momComp = 0; momComp != SpaceDim; ++momComp)
        {
          MD_BOXLOOP(a_solveBox,i)
            {
              Real momentum = a_UcellAvg[MD_IX(i,cVel + momComp)];
              a_sourceFab[MD_IX(i,cVel + momComp)] = B*momentum;
            }
        }
    }
  else if (m_turbForcing == 2)
    {
      Box box1Dom = grow(a_solveBox, 1);
      box1Dom &= a_domain;
      // Compute the cell-centered vorticity
      int order = ((CRDparam::g_cellDeconvolveFlatten < 2) &&
                   (CRDparam::g_cellConvolveFlatten < 2)) ? 4 : 2;
      const bool fourthOrder = (order == 4) ? true : false;
      FABSTACKTEMP(cellPntVelGrad, box1Dom, SpaceDim*SpaceDim);
      for (int row = 0; row != SpaceDim; ++row)
        {
          for (int col = 0; col != SpaceDim; ++col)
            {
              const int valComp = cVel + row;
              const int derivComp = row + col*SpaceDim;
              PatchMappedFunc::cellAvgDerivFromCellAvgCS(
                cellPntVelGrad, a_Wcell, valComp, derivComp, a_domain,
                box1Dom, a_gridMetrics.dxVect(), col, fourthOrder, false);
            }
        }

      // Map the cell-centered velocity gradient
      FABSTACKTEMP(cellPntPhysVelGrad, box1Dom, SpaceDim*SpaceDim);
      PatchMappedFunc::gradientCStoPS(
        box1Dom, cellPntPhysVelGrad, cellPntVelGrad,
        a_gridMetrics.m_cellNtJ[a_didx]);

      // Compute the cell-centered vorticity
      const int numVorticityComp = PatchMappedFunc::m_numCurlComps;
      FABSTACKTEMP(cellPntVorticity, box1Dom, numVorticityComp);
      PatchMappedFunc::curlPS(box1Dom, cellPntVorticity, cellPntPhysVelGrad);

      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          MD_BOXLOOP(a_solveBox,i)
            {
              a_sourceFab[MD_IX(i,cVel + dir)] =
                B*cellPntVorticity[MD_IX(i,dir)];
            }
        }
    }
  return;
}

/*--------------------------------------------------------------------*/
//  Output pnt-values and domain-sums specific to simulation
/** \param[in] a_WfaceAvgFxb
 *                      Face-averaged primitive state
 *  \param[in] a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in] a_WcellPntFab
 *                      Cell-centered primitive state
 *  \param[in] a_unitNormalsFxb
 *                      Unit-normal basis for faces
 *  \param[in] a_gridMetrics
 *                      Grid metrics for current level
 *  \param[in] a_time   Current time on current level
 *  \param[in] a_level  Current level
 *//*-----------------------------------------------------------------*/

void
CNSIBCTaylorGreen::inSituSumPntProcessing(
  const LevelData<FluxBox>&   a_faceAvgPlotData,
  const LevelData<FluxBox>&   a_WfaceAvgFxb,
  const LevelData<FArrayBox>& a_WcellAvgFab,
  const LevelData<FArrayBox>& a_WcellPntFab,
  const LayoutData<FluxBox>&  a_unitNormalsFxb,
  const LevelGridMetrics&     a_gridMetrics,
  const Real                  a_time,
  const int                   a_stage,
  const int                   a_level) const
{
  CH_TIME("CNSIBCTaylorGreen::inSituSumPntProcessing");

  const int cRho = CRDparam::g_CRDPhysics->densityIndex();
  const int cVel = CRDparam::g_CRDPhysics->vectorFluxInterval().begin();
  const int cPres = CRDparam::g_CRDPhysics->pressureIndex();
  const Real vol =
    a_gridMetrics.dxVect().product()/(CRDparam::g_domainLength.product());

  Real locMeanKE = 0.;
  Real locMeanEnstrophy = 0.;
  Real locMeanDivVel = 0.;
  Real locMeanDilatation = 0.;
  Real locMeanHelicity = 0.;
  Real locMeanEps1 = 0.;
  Real locMeanEps2 = 0.;
  Real locMeanEps3 = 0.;

  for (DataIterator dit = a_WcellAvgFab.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = a_WcellAvgFab.disjointBoxLayout()[dit];
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      const FluxBox& unitNormalsFxb = a_unitNormalsFxb[dit];
      const FluxBox& WfaceAvgFxb = a_WfaceAvgFxb[dit];
      const FArrayBox& WcellPntFab = a_WcellPntFab[dit];

      // (1) Compute and volume-weighted-sum the cell-averaged kinetic energy
      FABSTACKTEMP(cellPntKE, box1Dom, 1);
      MD_BOXLOOP(box1Dom, i)
        {
          const Real rho = WcellPntFab[MD_IX(i, cRho)];
          const RealVect vel(D_DECL(WcellPntFab[MD_IX(i, cVel)],
                                    WcellPntFab[MD_IX(i, cVel+1)],
                                    WcellPntFab[MD_IX(i, cVel+2)]));
          cellPntKE[MD_IX(i, 0)] = 0.5*rho*(vel.radSquared());
        }
      FABSTACKTEMP(cellAvgKE, box, 1);
      // If higher-order, convolve the kinetic energy
      int order = ((CRDparam::g_cellDeconvolveFlatten < 2) &&
                   (CRDparam::g_cellConvolveFlatten < 2)) ? 4 : 2;
      int convolveSign = -1;
      CRDutil::deconvolve(cellAvgKE, cellPntKE, box, blockDomain,
                          Interval(0,0), order, convolveSign, false);
      // Sum the cell-averaged kinetic energy
      MD_BOXLOOP(box, i)
        {
          locMeanKE += vol*cellAvgKE[MD_IX(i, 0)];
        }

      // (2) Compute the cell-centered velocity gradient
      const bool fourthOrder = (order == 4) ? true : false;
      FABSTACKTEMP(cellPntVelGrad, box1Dom, SpaceDim*SpaceDim);
      for (int row = 0; row != SpaceDim; ++row)
        {
          for (int col = 0; col != SpaceDim; ++col)
            {
              const int valComp = cVel + row;
              const int derivComp = row + col*SpaceDim;
              PatchMappedFunc::cellAvgDerivFromCellAvgCS(
                cellPntVelGrad, WcellPntFab, valComp, derivComp, blockDomain,
                box1Dom, a_gridMetrics.dxVect(), col, fourthOrder, false);
            }
        }

      // Map the cell-centered velocity gradient
      FABSTACKTEMP(cellPntPhysVelGrad, box1Dom, SpaceDim*SpaceDim);
      PatchMappedFunc::gradientCStoPS(
        box1Dom,cellPntPhysVelGrad,cellPntVelGrad,a_gridMetrics.m_cellNtJ[dit]);

      // Compute the cell-centered vorticity
      const int numVorticityComp = PatchMappedFunc::m_numCurlComps;
      FABSTACKTEMP(cellPntVorticity, box1Dom, numVorticityComp);
      PatchMappedFunc::curlPS(box1Dom, cellPntVorticity, cellPntPhysVelGrad);

      // (3) Compute the cell-averaged enstrophy --- only works in 3D
      FABSTACKTEMP(cellPntEnstrophy, box1Dom, 1);
      MD_BOXLOOP(box1Dom, i)
        {
          const RealVect vorticity(D_DECL(cellPntVorticity[MD_IX(i, 0)],
                                          cellPntVorticity[MD_IX(i, 1)],
                                          cellPntVorticity[MD_IX(i, 2)]));
          cellPntEnstrophy[MD_IX(i, 0)] = vorticity.radSquared();
        }
      // Convolve the cell-centered enstrophy
      FABSTACKTEMP(cellAvgEnstrophy, box1Dom, 1);
      CRDutil::deconvolve(cellAvgEnstrophy, cellPntEnstrophy, box, blockDomain,
                          Interval(0,0), order, convolveSign, false);
      // Sum the cell-averaged enstrophy
      MD_BOXLOOP(box, i)
        {
          locMeanEnstrophy += vol*cellAvgEnstrophy[MD_IX(i, 0)];
        }

      // (4) Compute the cell-averaged helicity
      FABSTACKTEMP(cellPntHelicity, box1Dom, 1);
      MD_BOXLOOP(box1Dom, i)
        {
          const RealVect vorticity(D_DECL(cellPntVorticity[MD_IX(i, 0)],
                                          cellPntVorticity[MD_IX(i, 1)],
                                          cellPntVorticity[MD_IX(i, 2)]));
          const RealVect vel(D_DECL(WcellPntFab[MD_IX(i, cVel)],
                                    WcellPntFab[MD_IX(i, cVel+1)],
                                    WcellPntFab[MD_IX(i, cVel+2)]));
          cellPntHelicity[MD_IX(i, 0)] = vorticity.dotProduct(vel);
        }
      // Convolve the cell-centered helicity
      FABSTACKTEMP(cellAvgHelicity, box1Dom, 1);
      CRDutil::deconvolve(cellAvgHelicity, cellPntHelicity, box, blockDomain,
                          Interval(0,0), order, convolveSign, false);
      // Sum the cell-averaged helicity
      MD_BOXLOOP(box, i)
        {
          locMeanHelicity += vol*cellAvgHelicity[MD_IX(i, 0)];
        }

      // (5) Compute the cell-averaged velocity divergence
      // Map the face-averaged velocity into face-normal components
      FLUXBOXSTACKTEMP(faceAvgVelFxb, box, SpaceDim);
      faceAvgVelFxb.copy(WfaceAvgFxb, cVel, 0, SpaceDim);
      FABSTACKTEMP(cellAvgDivVel, box, 1);
      cellAvgDivVel.setVal(0.);
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          Box faceBox = box;
          faceBox.surroundingNodes(dir);
          const int MD_ID(ii, dir);
          const Real invDx = 1./(a_gridMetrics.dxVect()[dir]);
          FArrayBox& faceAvgVelFab = faceAvgVelFxb[dir];
          const FArrayBox& unitNormalsFab = unitNormalsFxb[dir];
          PatchMappedFunc::forwardTransform(
            faceAvgVelFab, unitNormalsFab, faceBox);
          // Compute the difference in face-values in dir
          MD_BOXLOOP(box, i)
            {
              const Real hiFace = faceAvgVelFab[MD_OFFSETIX(i, +, ii, dir)];
              const Real loFace = faceAvgVelFab[MD_IX(i, dir)];
              cellAvgDivVel[MD_IX(i, 0)] += invDx*(hiFace - loFace);
            }
        }
      // Sum the cell-averaged absolute value of velocity divergence
      MD_BOXLOOP(box, i)
        {
          locMeanDivVel += vol*std::abs(cellAvgDivVel[MD_IX(i, 0)]);
          locMeanDilatation +=
            vol*(cellAvgDivVel[MD_IX(i, 0)])*(cellAvgDivVel[MD_IX(i, 0)]);
        }

      // (6) Compute epsilon1 (viscous dissipation of kinetic energy)
      // Compute cell-centered divergence of velocity
      FABSTACKTEMP(cellPntDivVel, box1Dom, 1);
      cellPntDivVel.setVal(0.);
      for (int row = 0; row != SpaceDim; ++row)
        {
          const int divComp = row + row*SpaceDim;
          MD_BOXLOOP(box1Dom, i)
            {
              cellPntDivVel[MD_IX(i, 0)] +=
                cellPntPhysVelGrad[MD_IX(i, divComp)];
            }
        }
      // Compute cell-centered viscous stress tensor
      const Real mu = CRDparam::g_mu;
      FABSTACKTEMP(cellPntTau, box1Dom, SpaceDim*SpaceDim);
      for (int row = 0; row != SpaceDim; ++row)
        {
          for (int col = 0; col != SpaceDim; ++col)
            {
              const int tauComp = row + col*SpaceDim;
              const int invComp = col + row*SpaceDim;
              MD_BOXLOOP(box1Dom, i)
                {
                  cellPntTau[MD_IX(i, tauComp)] =
                    mu*(cellPntPhysVelGrad[MD_IX(i, tauComp)]
                      + cellPntPhysVelGrad[MD_IX(i, invComp)]);
                }
              if (row == col)
                {
                  MD_BOXLOOP(box1Dom, i)
                    {
                      cellPntTau[MD_IX(i, tauComp)] -=
                        2.*mu*CRDparam::g_lambda*(cellPntDivVel[MD_IX(i, 0)]);
                    }
                }
            }
        }
      // Project cell-centered viscous stress tensor onto velocity gradient
      FABSTACKTEMP(cellPntEps1, box1Dom, 1);
      cellPntEps1.setVal(0.);
      for (int comp = 0; comp != (SpaceDim*SpaceDim); ++comp)
        {
          MD_BOXLOOP(box1Dom, i)
            {
              cellPntEps1[MD_IX(i, 0)] +=
                cellPntTau[MD_IX(i, comp)]*cellPntPhysVelGrad[MD_IX(i, comp)];
            }
        }
      // Convolve cell-centered epsilon1
      FABSTACKTEMP(cellAvgEps1, box, 1);
      CRDutil::deconvolve(cellAvgEps1, cellPntEps1, box, blockDomain,
                          Interval(0,0), order, convolveSign, false);
      // Sum cell-averaged epsilon1
      MD_BOXLOOP(box, i)
        {
          locMeanEps1 += vol*cellAvgEps1[MD_IX(i, 0)];
        }

      // (7) Compute epsilon2 (viscous compressive dissipation of ke)
      FABSTACKTEMP(cellPntEps2, box1Dom, 1);
      MD_BOXLOOP(box1Dom, i)
        {
          const Real lambda_b = 2.*CRDparam::g_lambda;
          cellPntEps2[MD_IX(i, 0)] =
            mu*lambda_b*cellPntDivVel[MD_IX(i, 0)]*cellPntDivVel[MD_IX(i, 0)];
        }
      // Convolve cell-centered epsilon2
      FABSTACKTEMP(cellAvgEps2, box, 1);
      CRDutil::deconvolve(cellAvgEps2, cellPntEps2, box, blockDomain,
                          Interval(0,0), order, convolveSign, false);
      // Sum cell-averaged epsilon2
      MD_BOXLOOP(box, i)
        {
          locMeanEps2 += vol*cellAvgEps2[MD_IX(i, 0)];
        }

      // (8) Compute epsilon3 (pressure dissipation of kinetic energy)
      FABSTACKTEMP(cellPntEps3, box1Dom, 1);
      MD_BOXLOOP(box1Dom, i)
        {
          const Real pressure = WcellPntFab[MD_IX(i, cPres)];
          cellPntEps3[MD_IX(i, 0)] =
            pressure*cellPntDivVel[MD_IX(i, 0)];
        }
      // Convolve cell-centered epsilon3
      FABSTACKTEMP(cellAvgEps3, box, 1);
      CRDutil::deconvolve(cellAvgEps3, cellPntEps3, box, blockDomain,
                          Interval(0,0), order, convolveSign, false);
      // Sum cell-averaged epsilon3
      MD_BOXLOOP(box, i)
        {
          locMeanEps3 += vol*cellAvgEps3[MD_IX(i, 0)];
        }
    }

  Real meanGlobalKE = locMeanKE;
  Real meanGlobalEnstrophy = locMeanEnstrophy;
  Real meanGlobalHelicity = locMeanHelicity;
  Real meanGlobalDivVel = locMeanDivVel;
  Real meanGlobalDilatation = locMeanDilatation;
  Real meanGlobalEps1 = locMeanEps1;
  Real meanGlobalEps2 = locMeanEps2;
  Real meanGlobalEps3 = locMeanEps3;
#ifdef CH_MPI
  MPI_Allreduce(&locMeanKE, &meanGlobalKE, 1,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

  MPI_Allreduce(&locMeanEnstrophy, &meanGlobalEnstrophy, 1,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

  MPI_Allreduce(&locMeanHelicity, &meanGlobalHelicity, 1,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

  MPI_Allreduce(&locMeanDivVel, &meanGlobalDivVel, 1,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

  MPI_Allreduce(&locMeanDilatation, &meanGlobalDilatation, 1,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

  MPI_Allreduce(&locMeanEps1, &meanGlobalEps1, 1,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

  MPI_Allreduce(&locMeanEps2, &meanGlobalEps2, 1,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

  MPI_Allreduce(&locMeanEps3, &meanGlobalEps3, 1,
                MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
#endif
  // We use pout so the text doesn't line-break from CRDmsg
  pout() << "t = " << a_time << "; Level = " << a_level
         << "; KE = " << meanGlobalKE << "; Ens = " << meanGlobalEnstrophy
         << "; Hel = " << meanGlobalHelicity << "; DivU = "
         << meanGlobalDivVel << "; DivU_Squared = " << meanGlobalDilatation
         << "; Eps1 = " << meanGlobalEps1 << "; Eps2 = " << meanGlobalEps2
         << "; Eps3 = " << meanGlobalEps3 << std::endl;
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCTaylorGreen::haveExactSol() const
{
  return false;
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
CNSIBCTaylorGreen::readBCInfo()
{
  ParmParse ppIBC("ibc");

  // Turbulent initialization state (turbulent velocity magnitude)
  const char* initStateNameHdr = "turb_state";
  if (ppIBC.contains(initStateNameHdr))
    {
      std::string initStateName;
      ppIBC.get(initStateNameHdr, initStateName);
      m_turbState = CRDState::nameIndex(initStateName);
    }
  else
    {
      CRD::msg << "Input (Taylor-Green Vortex IBC): " << initStateNameHdr
               << " must be given" << CRD::error;
    }
  ppIBC.query("z_velocity_initialization", m_zVelIC);
  ppIBC.query("initialize_uniform_density", m_uniformRhoIC);
  m_numVort = RealVect::Zero;
  std::vector<Real> tnumVort(SpaceDim);
  ppIBC.queryarr("vortex_number", tnumVort, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_numVort.dataPtr(),
                                           &tnumVort.front());

  ParmParse ppTURB("turb");
  ppTURB.query("turb_forcing", m_turbForcing);
  ppTURB.query("mean_energy_controller_amplitude", m_coeffG);

  m_readInput = true;
}
