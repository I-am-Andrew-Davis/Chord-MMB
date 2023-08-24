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
 * \file CNSIBCMixedCouette.cpp
 *
 * \brief Member functions for CNSIBCMixedCouette
 *
 *//*+*************************************************************************/
#include <fstream>
//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
//**FIXME
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"
 
//----- Internal -----//

#include "DataTemp.H"
#include "CNSIBCMixedCouette.H"
#include "CNSIBCTransientCouetteF_F.H"
#include "CNSIBCF_F.H"
#include "CRDPhysics.H"
#include "ChordInput.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"


/*******************************************************************************
 *
 * Class CNSIBCMixedCouette: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_wallVelocity
 *                      Velocity of the wall.  Default is 
 *                      \f[
 *                        \sqrt{\frac{V^2}{\mathrm{SpaceDim}-1}}
 *                      \f]
 *                      where \f$V=\frac{\mathrm{Re}\nu}
 *                      {\mathrm{Domain length}}\f$.  The normal
 *                      component is forced to zero.
 *  \param[in]  a_wallNormalDir
 *                      Direction normal to the wall.  Default is
 *                      SpaceDim-1
 *//*-----------------------------------------------------------------*/

CNSIBCMixedCouette::CNSIBCMixedCouette(
  const RealVect a_wallVelocity,
  const int      a_wallNormalDir)
  :
  CNSIBCReferenceCubeBC(),
  m_wallVelocity(a_wallVelocity),
  m_wallNormalDir(a_wallNormalDir),
  m_movWall(0),
  m_initVel(-1.),
  m_sonicTunnel(false)
{

  readBCInfo();

  // Set BC Type
  setAllDomainBC(CRDparam::DomainBCTypePeriodic);
  // Set stationary adiabatic walls
  if(m_sonicTunnel)
    { 
      setDomainBC(m_wallNormalDir, 0, CRDparam::DomainBCTypeAdiabaticWall);
      setDomainBC(m_wallNormalDir, 1, CRDparam::DomainBCTypeAdiabaticWall);
    }
  else // Otherwise, set to Couette flow walls
    {
      // NOTE: We set both walls to be adiabatic walls because we can set
      //  the moving wall part within this file
      int movSide = 0;
      int statSide = 1;
      if(m_movWall == 1)
        {
          movSide = 1;
          statSide = 0;
        }
      // Low side is moving
      setDomainBC(m_wallNormalDir, movSide,
                  CRDparam::DomainBCTypeAdiabaticWall);
      // Hi side is stationary
      setDomainBC(m_wallNormalDir, statSide,
                  CRDparam::DomainBCTypeAdiabaticWall);
    }

  // Set conditions on the moving wall
  m_wallVelocity[m_wallNormalDir] = 0.;
  D_TERM(
    setReferenceBCState(m_wallNormalDir, 0, WVELX, m_wallVelocity[0]);,
    setReferenceBCState(m_wallNormalDir, 0, WVELY, m_wallVelocity[1]);,
    setReferenceBCState(m_wallNormalDir, 0, WVELZ, m_wallVelocity[2]);)
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCMixedCouette::~CNSIBCMixedCouette()
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
CNSIBCMixedCouette::IBCName() const
{
  return "Cold mixed Couette flow";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCMixedCouette::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Initial mass fractions\n(";
  CRD::msg << m_initMassFraction[0];
  for(int species = 1; species != CRDparam::g_numSpecies; ++species)
    {
      CRD::msg << ", " << m_initMassFraction[species];
    }
  CRD::msg << ")" << CRD::var;
  CRD::msg << "Wall velocity\n(";
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      if (dir != 0)
        {
          CRD::msg << ',';
        }
      CRD::msg << m_wallVelocity[dir];
    }
  CRD::msg << ')' << CRD::var;
  CRD::msg << "Wall normal direction\n" << m_wallNormalDir << CRD::var;
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
CNSIBCMixedCouette::initialize(LevelData<FArrayBox>&      a_U,
                               LevelGridMetrics&          a_gridMetrics,
                               const LayoutData<FluxBox>& a_unitNormals,
                               const Real                 a_time,
                               const int                  a_level) const
{
  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box2Dom = grow(box, 2);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, 1);
      box1Dom &= blockDomain;

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      FArrayBox& UFab = a_U[dit];
      this->exactSol(UFab,
                     box1Dom,
                     box,
                     a_gridMetrics,
                     a_unitNormals[dit],
                     dit(),
                     a_time,
                     a_level);
    }
}

/*--------------------------------------------------------------------*/
//  Set the tagging method if one can be associated with IBC
/** The default implementation just sets the buffer and otherwise does
 *  no tagging.
 *  \param[in]  a_tagBufferSize
 *                      Requested tag buffer size (should be
 *                      respected).
 *
 *  \note
 *  <ul>
 *    <li> Allocate levels and methods with 'new' and do not delete
 *  </ul>
 *//*-----------------------------------------------------------------*/

TagLevelFactory*
CNSIBCMixedCouette::setTagMethod(const int a_tagBufferSize)
{
  const int maxAMRLevel = CRDparam::maxAMRLevel();
  TagLevel* tagLevel;
  std::vector<TagLevel*> tagLevelVec;
  std::vector<int> levelMapVec;
  const RealVect dxVect = CRDparam::g_domainLength/
    CRDparam::g_domainBaseSize;

  switch (maxAMRLevel)
    {
    case 0:
    case 1:
    {
      // One TagLevel
      tagLevelVec.resize(1);
      levelMapVec.resize(1);
      // Apply this TagLevel to all levels
      levelMapVec[0] = 0;
      // Tag a box grown by the max displacement of the BL
      tagLevel = new TagLevel;
      tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
      IntVect lo;
      IntVect hi;
      for(int dir = 0; dir != SpaceDim; ++dir)
        {
          lo[dir] = std::ceil(CRDparam::g_domainLength[dir]/dxVect[dir]*2/8);
          hi[dir] = std::floor(CRDparam::g_domainLength[dir]/dxVect[dir]*6/8-1);
        }
      Box tagBox(lo,hi);
      tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      // Add in tag buffer
      tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
      tagLevelVec[0] = tagLevel;
      break;
    }
    case 2:
    {
      // Two TagLevel
      tagLevelVec.resize(2);
      levelMapVec.resize(2);
      // Apply the first to level 0 and the second to level 1
      levelMapVec[0] = 0;
      levelMapVec[1] = 1;

      // First TagLevel is a box grown by 4 times the max displacement of the BL
      // This is supposed to do a better resolution of the inviscid flow
      tagLevel = new TagLevel;
      tagLevel->FIXMEsetTagBufferSize(0);
      Box tagBox;
      IntVect lo;
      IntVect hi;
      for(int dir = 0; dir != SpaceDim; ++dir)
        {
          lo[dir] = std::ceil(CRDparam::g_domainLength[dir]/dxVect[dir]*2/8);
          hi[dir] = std::floor(CRDparam::g_domainLength[dir]/dxVect[dir]*6/8-1);
        }
      tagBox.define(lo,hi);
      tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      // No tag buffer
      tagLevelVec[0] = tagLevel;
      for(int dir = 0; dir != SpaceDim; ++dir)
        {
          lo[dir] = std::ceil(CRDparam::g_domainLength[dir]/dxVect[dir]*2/8);
          hi[dir] = std::floor(CRDparam::g_domainLength[dir]/dxVect[dir]*6/8-1);
        }
      tagBox.define(lo,hi);
      // Second TagLevel is a box grown by the max displacement of the BL
      tagLevel = new TagLevel;
      tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);

      tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      // Add in tag buffer
      tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
      tagLevelVec[1] = tagLevel;
      break;
    }
    }

    // Return the factory
    return new TagLevelFactory(tagLevelVec, levelMapVec);
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCMixedCouette::haveExactSol() const
{
  if(m_sonicTunnel)
    {
      return false;
    }
  else
    {
      return true;
    }
}

/*--------------------------------------------------------------------*/
//  Compute the exact solution state \<U\> in the cells
/** The exact solution is for \f$Re = 200\f$, air at 300K
 *  \f$\nu = 1.56e-5 m^2/s\f$ \f$h = 0.1m\f$, solution time = 500s,
 *  the number n in the series is 500.  The exact solution form
 *  (moving wall is on the low side of y direction) is
 *  \f[
 *     U(y,t) = U_0\left(1 - \frac{y}{h}\right) -
 *       2 U_0\sum\limits_{n=1}^{N}\left(\frac{\sin(n\pi \frac{y}{h})}
 *       {n\pi}e^{-\frac{\displaystyle n^2\pi^2\nu t}
 *       {\displaystyle h^2}}\right)
 *  \f]
 *  \param[out] a_Ux    Exact solution
 *  \param[in]  a_box   Cell locations to update
 *  \param[in]  a_problemDomain
 *                      The problem domain on the level
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
CNSIBCMixedCouette::exactSol(FArrayBox&              a_Ux,
                             const Box&              a_box,
                             const Box&              a_disjointBox,
                             const LevelGridMetrics& a_gridMetrics,
                             const FluxBox&          a_unitNormals,
                             const DataIndex&        a_didx,
                             const Real              a_time,
                             const int               a_level) const
{
  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rComp = CRDparam::g_CRDPhysics->densityIndex();
  const int pComp = CRDparam::g_CRDPhysics->pressureIndex();
  const int tComp = CRDparam::g_CRDPhysics->temperatureIndex();
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  const int nTerms = 5000;  // 40K at y = 0.0001*h
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_disjointBox);
  const Real rho = CRDparam::g_rho;
  const Real tval = CRDparam::g_T;
  // Working set boxes
  Box box1Dom = grow(a_box, 1);
  box1Dom &= blockDomain;
  Box initBox(box1Dom);
  // Create FAB containing the primitive variables
  FABSTACKTEMP(Wc, initBox, numWVar);
  Wc.setVal(rho, rComp);
  Wc.setVal(0., initBox, WvelIndx, SpaceDim);
  Wc.setVal(-1., pComp);
  Wc.setVal(tval, tComp);
  for(int i = 0; i != numSpecies; ++i)
    {
      int comp = wCompStart + i;
      Wc.setVal(m_initMassFraction[i], comp);
    }
  if(m_sonicTunnel)
    {
      Wc.setVal(m_initVel, initBox, WvelIndx);
    }
  // Get physical coordinates
  FABSTACKTEMP(XiFab, initBox, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, initBox, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(initBox, XiFab, XFab, blockCoordSys);
  
  // Set all the components to state with zero velocity before computing the
  // exact solution unless m_initVel is specified in input file
  const Real nu = CRDparam::g_mu/CRDparam::g_rho;
  const Real height = CRDparam::g_domainLength[m_wallNormalDir];
  CRDparam::g_CRDPhysics->initialize(a_Ux,
                                     Wc,
                                     a_gridMetrics,
                                     a_unitNormals,
                                     a_didx,
                                     a_disjointBox,
                                     initBox);
  if(m_sonicTunnel)
    {
      return 0;
    }
  const Real tstar = nu*a_time/(height*height);
  // At time 0, the velocity everywhere is 0
  if (tstar <= 0.)
    {
      fourthOrderAverageCell(a_Ux, blockDomain, a_box);
      return 0;
    }

  // Between 0 and 0.01, we have not validated the exact solution
  if (tstar < 0.01)
    {
      Real finalTime = height*height*0.01/nu;
      CRD::msg << "Analytic solution time must be > " << finalTime
               << CRD::error;
    }

  Real wallSpeed = m_wallVelocity.vectorLength()/rho;
  if (wallSpeed == 0.)
    {
      return 0;
    }
  // comp is the component to solve for in a_Ux
  int comp = -1;
  D_TERM(
  if(m_wallVelocity[0] > 0.)
    {
      comp = UMOMX;
    },
  else if(m_wallVelocity[1] > 0.)
    {
      comp = UMOMY;
    },
  else if(m_wallVelocity[2] > 0.)
    {
      comp = UMOMZ;
    });
  FORT_CNSIBCTRANSIENTCOUETTEEXACTSOLPRIM(
    CHF_FRA(Wc),
    CHF_BOX(initBox),
    CHF_CONST_FRA1(XFab,m_wallNormalDir),
    CHF_CONST_REAL(a_time),
    CHF_CONST_REAL(height),
    CHF_CONST_REAL(nu),
    CHF_CONST_REAL(wallSpeed),
    CHF_CONST_INT(m_movWall),
    CHF_CONST_INT(nTerms),
    CHF_CONST_INT(comp));

  CRDparam::g_CRDPhysics->initialize(a_Ux,
                                     Wc,
                                     a_gridMetrics,
                                     a_unitNormals,
                                     a_didx,
                                     a_disjointBox,
                                     initBox);
  fourthOrderAverageCell(a_Ux, blockDomain, a_box);
  if(!m_outputFileName.empty()
#ifdef CH_MPI
     && CHprocID() == 0
#endif
    )
    {
      const int numIx = 101;
      std::vector<Real> yvec(numIx);
      std::vector<Real> usol(numIx);
      const Real datDX = height/(numIx-1);
      for(int i = 0; i != numIx; ++i)
        {
          yvec[i] = i*datDX;
        }
      FORT_CNSIBCEXACTSOLOUTPUT(
        CHF_VR(usol),
        CHF_CONST_VR(yvec),
        CHF_CONST_REAL(a_time),
        CHF_CONST_REAL(height),
        CHF_CONST_REAL(nu),
        CHF_CONST_REAL(wallSpeed),
        CHF_CONST_INT(m_movWall),
        CHF_CONST_INT(nTerms));
      std::ofstream udat;
      std::string instring(m_outputFileName);
      char* conv_cstr = new char[instring.length()+1];
      strcpy(conv_cstr,instring.c_str());
      udat.open(conv_cstr);
      udat << "y u" << endl;
      for(int i = 0; i != numIx; ++i)
        {
          udat << yvec[i] << " " << usol[i] << endl;
        }
      udat.close();
    }
  return 0;
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set the primitive state at wall BC
/** In this override, the wallVelocity is not zero on low-side walls
 *  \param[in]  a_Wface Primitive state on wall from interior scheme
 *                      (face-centered)
 *  \param[out] a_Wface Primitive state at exterior to boundary
 *  \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
 *  \param[in]  a_Wcell Primitive state in cells at interior of
 *                      domain.  The cells have been shifted by half
 *                      so that the first interior layer of cells
 *                      overlaps a_Wface.
 *  \param[in]  a_boundarySlipVelocity
 *                      Slip velocities on the boundary provided by
 *                      turbulence model wall-model
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
 *  \param[in]  a_domT  Unused for this case
 *//*-----------------------------------------------------------------*/

void
CNSIBCMixedCouette::setWallBCprimState(
  FArrayBox&           a_Wface,
  const Box&           a_boundaryFaceBox,
  const FArrayBox&     a_Wcell,
  const FArrayBox&     a_boundarySlipVelocity,
  const FArrayBox&     a_unitNormalBasisFab,
  const BoundaryIndex& a_bcIdx,
  const Box&           a_disjointBox,
  LevelGridMetrics&    a_gridMetrics,
  const Real           a_time,
  const int            a_level,
  const BCInfo&        a_domT) const
{
  const int lohiSign = Side::sign(a_bcIdx.m_side);
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();

  int viscousSlip = 0;
  if (CRDparam::DomainBCTypeSlipWall & a_domT.m_type)
    {
      viscousSlip = 1;
    }

  RealVect wallVelocity = RealVect::Zero;
  if ((a_bcIdx.m_side == Side::Lo && m_movWall == 0) ||
      (a_bcIdx.m_side == Side::Hi && m_movWall == 1))
    {
      CH_assert(m_wallVelocity[a_bcIdx.m_dir] == 0.);
      wallVelocity = m_wallVelocity;
    }
  // Gamma values in cells
  FABSTACKTEMP(gammaCell, a_boundaryFaceBox, 1);
  // Gamma values on faces
  FABSTACKTEMP(gammaFace, a_boundaryFaceBox, 1);
  CRDparam::g_CRDPhysics->calcGamma(a_boundaryFaceBox,
                                    gammaCell,
                                    a_Wcell);
  CRDparam::g_CRDPhysics->calcGamma(a_boundaryFaceBox,
                                    gammaFace,
                                    a_Wface);
  FORT_CNSIBCPRIMSTATEWALL(CHF_FRA(a_Wface),
                           CHF_BOX(a_boundaryFaceBox),
                           CHF_CONST_FRA(a_Wcell),
                           CHF_CONST_FRA1(gammaFace,0),
                           CHF_CONST_FRA1(gammaCell,0),
                           CHF_CONST_REALVECT(wallVelocity),
                           CHF_CONST_INT(rhoIndx),
                           CHF_CONST_INT(presIndx),
                           CHF_CONST_INT(velIndx),
                           CHF_CONST_INT(viscousSlip),
                           CHF_CONST_INT(a_bcIdx.m_dir),
                           CHF_CONST_INT(lohiSign));
  CRDparam::g_CRDPhysics->temperature(a_Wface, a_boundaryFaceBox);
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
CNSIBCMixedCouette::readBCInfo()
{
  CH_assert(!m_readInput);
  int numSpecies = CRDparam::g_numSpecies;
  ParmParse ppIBC("ibc");
  ppIBC.query("wall_normal_dir", m_wallNormalDir);
  if (m_wallNormalDir < 0 || m_wallNormalDir >= SpaceDim)
    {
      CRD::msg << "Input (MixedCouette IBC): 'wall_normal_dir' must "
        "be >= 0 and < " << SpaceDim << '!' << CRD::error;
    }
  // Determine if low wall is moving or high wall
  ppIBC.query("moving_wall", m_movWall);
  // Solve wall velocity using Reynolds number if "wall_velocity" isn't
  // specified in input file
  const Real height = CRDparam::g_domainLength[m_wallNormalDir];
  const Real nu = CRDparam::g_mu/CRDparam::g_rho;
  int velDir = m_wallNormalDir - 1;
  if(velDir < 0)
    {
      velDir = SpaceDim - 1;
    }
  const Real wallSpeed = CRDparam::g_Re*nu/height;
  m_wallVelocity = RealVect::Zero;
  m_wallVelocity[velDir] = wallSpeed;
  // Wall velocity will be set to input file value
  std::vector<Real> IBCwallVelocity(SpaceDim);
  if(ppIBC.contains("wall_velocity"))
    { 
      ppIBC.queryarr("wall_velocity", IBCwallVelocity, 0, SpaceDim);
      SpaceDimArray<Real, Real>::loadFromArray(m_wallVelocity.dataPtr(),
                                               &IBCwallVelocity.front());
    }
  //**FIXME This needs to be update to the new BC format (address other fixmes
  //**      in this file as well)
  CH_assert(false);
  // if (m_movWall == 1)
  //   {
  //     m_hiBCVel[m_wallNormalDir][velDir] = wallSpeed;
  //   }
  // else
  //   {
  //     m_loBCVel[m_wallNormalDir][velDir] = wallSpeed;
  //   }
  // Specify initial velocity value of the flow
  m_initVel = 0.;
  ppIBC.query("initial_velocity", m_initVel);
  // Specify the initial mass fractions in the flow field
  m_initMassFraction.resize(numSpecies);
  m_initMassFraction.assign(numSpecies,0.);
  // Call function to assign mass fraction values
  int massCheck = assignMassFractions(m_initMassFraction,
                                      "init_specs",
                                      "init_mfs");
  if(massCheck == 1)
    {
      CRD::msg << "Input (MixedCouette IBC): 'mass_fractions' must "
        "be equal to 1!" << CRD::error;
    }
  ppIBC.query("dat_file",m_outputFileName);
  m_sonicTunnel = false;
  ppIBC.query("sonic_tunnel", m_sonicTunnel);
  // FIXME: This is a quick fix necessary for artificial viscosity
  
  m_readInput = true;
}
