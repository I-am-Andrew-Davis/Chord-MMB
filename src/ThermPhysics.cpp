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
 * \file ThermPhysics.cpp
 *
 * \brief Member functions for ThermPhysics
 *
 *//*+*************************************************************************/

#include <iostream>
#include <string>
#include <fstream>
#include <array>

//----- Chombo -----//

#include "LoHiCenter.H"
#include "CellToEdge.H"
#include "CH_System.H"
#include "EdgeToCellF_F.H"
#include "MaxScaledAcousticSpeedF_F.H"
#include "FourthOrderUtil.H"
#include "LevelGridMetrics.H"
#include "ParmParse.H"
#include "RootSolver.H"
#include "CHMatrixOps.H"
#include "InsertionSort.H"
#include "Misc.H"

//----- Internal -----//

#include "ThermPhysics.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "PatchMappedFunc.H"
#include "DataTemp.H"
#include "FileThermParser.H"
#include "FileThermCKParser.H"
#include "RiemannFunctions.H" // Also contains RFunc, C2PFunc


/*==============================================================================
 * Define how is a fortran function named
 *============================================================================*/

#ifdef CH_FORT_UPPERCASE
  #ifdef CH_FORT_UNDERSCORE
    #define FORTRAN_NAME( NAME ,name ) NAME ## _
  #else
    #define FORTRAN_NAME( NAME ,name ) NAME
  #endif
#else
  #ifdef CH_FORT_UNDERSCORE
    #define FORTRAN_NAME( NAME ,name ) name ## _
  #else
    #define FORTRAN_NAME( NAME ,name ) name
  #endif
#endif

#define LBFGSBINTERFACE FORTRAN_NAME(LBFGSB_INTERFACE, lbfgsb_interface)

extern "C" void LBFGSBINTERFACE(int*,  int*, Real*, Real*, Real*, int*, Real*,
                                Real*, Real*, Real*, Real*, int*, char*, int*,
                                char*, int*, int*, Real*);


/*******************************************************************************
 *
 * Class ThermPhysics: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default Constructor
/**
 *//*-----------------------------------------------------------------*/

ThermPhysics::ThermPhysics()
  :
  m_Rmol(8.3144598), // Universal gas constant in J/mol-K
  m_molMass(CRDparam::g_numSpecies,0.),
  m_Rn(CRDparam::g_numSpecies,1.),
  m_riemannSolver(RiemannSolver::Adaptive),
  m_reactionFile(""),
  m_thermoFile(""),
  m_transFile(""),
  m_thermFileFormat(0),
  m_lewisNum(1.),
  m_schmidtNum(-1.),
  m_lookupLoT(200.),
  m_lookupHiT(4000.),
  m_lookupDelT(1.),
  m_midLookupT(1000.),
  m_lookupSize(1),
  m_lookupH(CRDparam::g_numSpecies),
  m_lookupMu(CRDparam::g_numSpecies),
  m_lookupKappa(CRDparam::g_numSpecies),
  m_interpN(1),
  m_diffCL(1,0),
  m_interpDenoms(1,0.),
  m_H0(CRDparam::g_numSpecies,0.),
  m_hnA1L(CRDparam::g_numSpecies,-1.),
  m_hnA2L(CRDparam::g_numSpecies,0.),
  m_hnA3L(CRDparam::g_numSpecies,0.),
  m_hnA4L(CRDparam::g_numSpecies,0.),
  m_hnA5L(CRDparam::g_numSpecies,0.),
  m_hnA6L(CRDparam::g_numSpecies,0.),
  m_hnA7L(CRDparam::g_numSpecies,0.),
  m_hnB1L(CRDparam::g_numSpecies,0.),
  m_hnB2L(CRDparam::g_numSpecies,0.),
  m_hnA1H(CRDparam::g_numSpecies,-1.),
  m_hnA2H(CRDparam::g_numSpecies,0.),
  m_hnA3H(CRDparam::g_numSpecies,0.),
  m_hnA4H(CRDparam::g_numSpecies,0.),
  m_hnA5H(CRDparam::g_numSpecies,0.),
  m_hnA6H(CRDparam::g_numSpecies,0.),
  m_hnA7H(CRDparam::g_numSpecies,0.),
  m_hnB1H(CRDparam::g_numSpecies,0.),
  m_hnB2H(CRDparam::g_numSpecies,0.),
  m_altSpecNames(CRDparam::g_numSpecies, "U"),
  m_N2idx(-1),
  m_intSeq(CRDparam::g_numSpecies),
  m_preAmodUnits(1.E-6), // Default converts m^3 to cm^3
  m_numReactions(0), // Number of reactions
  m_numTBReact(0),   // Number of third body reactions
  m_numPRReact(0),
  m_numREVReact(0),
  m_numArbReact(0),
  m_numRects(1,0),
  m_numProds(1,0),
  m_Tcutoff(300.),
  m_lookupKfwd(CRDparam::g_numReactions),
  m_lookupKbkwd(CRDparam::g_numReactions),
  m_lookupKeq(CRDparam::g_numReactions),
  m_maxSpecPerReact(3), 
  m_Fc(1.),
  m_betai(CRDparam::g_numReactions),
  m_EAR(CRDparam::g_numReactions)
{
  CRD::msg << CRD::fv3 << "ThermPhysics::ThermPhysics " << CRD::end;
  const int numSpecies = CRDparam::g_numSpecies;
  // Set the maximum number of reactants and products allowed per reaction
  // FIXME: This will need to be increased if any reactions have more
  // than 3 reactants or products
  m_maxSpecPerReact = 3;
  // Set the universal gas constant using Boltzmann constant and Avo's number
  const Real K = 1.38064852E-23;
  const Real Na = 6.0221409E23;
  m_Rmol = K*Na;
  ParmParse ppTHERM("therm");
  // Apply alternate names for species, like Ar = AR
  // Must be hard-coded in for current implementation
  // FIXME: List is nowhere near complete
  // FIXME: Hydrocarbon species names are complicated and are listed multiple
  // times. Must determine a fix for this.
  for (int cn = 0; cn != numSpecies; ++cn)
    {
      std::string spec = CRDparam::g_speciesNames[cn];
      m_altSpecNames[cn] = spec;
      if (spec == "Ar" || spec == "He" || spec == "Be")
        {
          m_altSpecNames[cn][1] = toupper(m_altSpecNames[cn][1]);
        }
      else if (spec == "OH")
        {
          m_altSpecNames[cn] = "HO";
        }
    }
  ppTHERM.query("thermo_file_format", m_thermFileFormat);
  FileThermParser* thermParser = nullptr;
  if (m_thermFileFormat == 0)
    {
      thermParser = new FileThermParser(*this);
    }
  else if (m_thermFileFormat == 1)
    {
      thermParser = new FileThermCKParser(*this);
    }
  else
    {
      CRD::msg << "Input: 'thermo_file_format' not recognized" << CRD::error;
    }
  thermParser->setLookupTableVals();
  // Read file name containing thermodynamic data
  m_thermoFile = "Mechanism_files/thermo.inp";
  ppTHERM.query("thermo_file", m_thermoFile);
  thermParser->readThermFile();
  thermParser->fillThermoTables();
  if (ppTHERM.contains("reaction_file") && CRDparam::g_numReactions > 0)
    {
      ppTHERM.get("reaction_file", m_reactionFile);
      if (!CH_System::fileExists(m_reactionFile.c_str()))
        {
          CRD::msg << "Input: 'reaction_file' " << m_reactionFile
                   << " does not exist!" << CRD::error;
        }
      thermParser->readReactionFile();
      // Any temperature below this value, reactions are not solved
      ppTHERM.query("cutoff_temperature", m_Tcutoff);
      if (m_Tcutoff < 0)
        {
          CRD::msg << "Input: 'cutoff_temperature' must be >= 0!"
                   << CRD::error;
        }
      if (m_numPRReact > 0)
        {
          // Change the Lindemann fit number if necessary
          ppTHERM.query("lindemann_fit_number",m_Fc);
        }
    }
  if ((CRDparam::g_physicsModels & CRDparam::PhysicsViscous) &&
      (CRDparam::g_mu < 0.))
    {
      // Make vector of arrays to for mu and kappa coefficients
      std::vector<std::array<Real, 4> > muVecL(numSpecies);
      std::vector<std::array<Real, 4> > muVecH(numSpecies);
      std::vector<std::array<Real, 4> > kVecL(numSpecies);
      std::vector<std::array<Real, 4> > kVecH(numSpecies);
      // Read thermal transport file name
      m_transFile = "Mechanism_files/trans.inp";
      ppTHERM.query("transport_file", m_transFile);
      // Read trans input file and set muVec and kVec
      thermParser->readTransportFile(muVecL, muVecH, kVecL, kVecH);
      // Fill lookup tables for transport properties
      thermParser->fillTransportTables(muVecL, muVecH, kVecL, kVecH);
    }
  // Use species Lewis numbers by default unless Schmidt number is specified
  if (ppTHERM.contains("schmidt_number") &&
      (CRDparam::g_physicsModels & CRDparam::PhysicsViscous))
    {
      ppTHERM.get("schmidt_number", m_schmidtNum);
      if (m_schmidtNum <= 0.)
        {
          CRD::msg << "Input: 'schmidt_number' must be > 0!" << CRD::error;
        }
    }
  else if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
    {
      ppTHERM.query("lewis_number", m_lewisNum);
      if (m_lewisNum <= 0.)
        {
          CRD::msg << "Input: 'lewis_number' must be > 0!" << CRD::error;
        }
    }
  std::string riemannSolverType = "adaptive";
  ppTHERM.query("riemann_solver", riemannSolverType);
  if (riemannSolverType == "adaptive")
    {
      m_riemannSolver = RiemannSolver::Adaptive;
    }
  else if (riemannSolverType == "exact")
    {
      m_riemannSolver = RiemannSolver::Exact;
    }
  else if (riemannSolverType == "approx")
    {
      m_riemannSolver= RiemannSolver::Approx;
    }
  else
    {
      CRD::msg << "ThermPhysics::ThermPhysics: 'riemann_solver' not recognized."
               << " Must be 'exact', 'approx', or 'adaptive'!" << CRD::error;
    }
  delete thermParser;

  // Determine which component is N2 for doing species correction
  for (int sp = 0; sp != numSpecies; ++sp)
    {
      auto ins = m_mapSpecNameToIndex.insert(
        { CRDparam::g_speciesNames[sp], sp });
      CH_assert(ins.second); // Insertion must have succeeded
      if (CRDparam::g_speciesNames[sp] == "N2")
        {
          m_N2idx = sp;
        }
    }

  // Set the indirection array for a non-sorted array to an integer sequence
  for (int j = 0; j != numSpecies; ++j)
    {
      m_intSeq[j] = j;
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

ThermPhysics::~ThermPhysics()
{
}

/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Zero-based index of a species
/** \param[in]  a_specName
 *                      Species name
 *  \return             Index counting from zero
 *//*-----------------------------------------------------------------*/

int
ThermPhysics::speciesZeroIndex(const std::string& a_specName) const
{
  auto iter = m_mapSpecNameToIndex.find(a_specName);
  return (iter == m_mapSpecNameToIndex.end()) ? -1 : iter->second;
}

/*--------------------------------------------------------------------*/
//  Component index of species named 'a_specName' in primitive
//  variables
/** \param[in]  a_specName
 *                      Species name
 *  \return             Component index in range of primitive
 *                      variables
 *//*-----------------------------------------------------------------*/

int
ThermPhysics::speciesPrimIndex(const std::string& a_specName) const
{
  int idx = speciesZeroIndex(a_specName);
  if (idx != -1) idx += speciesPrimInterval().begin();
  return idx;
}

/*--------------------------------------------------------------------*/
//  Component index of species named 'a_specName' in conservative
//  variables
/** \param[in]  a_specName
 *                      Species name
 *  \return             Component index in range of conservative
 *                      variables
 *//*-----------------------------------------------------------------*/

int
ThermPhysics::speciesConsIndex(const std::string& a_specName) const
{
  int idx = speciesZeroIndex(a_specName);
  if (idx != -1) idx += speciesConsInterval().begin();
  return idx;
}

/*--------------------------------------------------------------------*/
//  Compute the speed of sound
/** \param[out] a_speed Sound speed in cells
 *  \param[in]  a_W     Native primitive state
 *  \param[in]  a_box   Box to compute sound speed in
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::soundSpeed(FArrayBox&       a_speed,
                         const FArrayBox& a_W,
                         const Box&       a_box) const
{
  CH_TIME("ThermPhysics::soundSpeed");
  CH_assert(a_W.contains(a_box));
  CH_assert(a_speed.contains(a_box));
  const int primSpecStart = speciesPrimInterval().begin();
  const int rhoIndx = densityIndex();
  const int presIndx = pressureIndex();
  const int numSpecies = CRDparam::g_numSpecies;
  VECSTACKTEMPSIZE(spec, Real, numSpecies);
  MD_ARRAY_RESTRICT(arrW, a_W);
  MD_ARRAY_RESTRICT(arrC, a_speed);
  MD_BOXLOOP(a_box, i)
    {
      Real pres = arrW[MD_IX(i, presIndx)];
      Real rho = arrW[MD_IX(i, rhoIndx)];
      Real Rgas = 0.;
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          int specIndx = primSpecStart + sp;
          spec[sp] = arrW[MD_IX(i, specIndx)];
          Rgas += spec[sp]*m_Rn[sp];
        }
      const Real T = pres/(rho*Rgas);
      const Real gammaVal = gamma(T, spec.data());
      arrC[MD_IX(i, 0)] = std::sqrt(gammaVal*pres/rho);
    }
}

/*--------------------------------------------------------------------*/
//  Compute the intertial flux from primitive variable values on a
//  face
/** \param[out] a_flux  Flux on the faces
 *  \param[in]  a_WFace Primitive state on the face
 *  \param[in]  a_dir   Direction of the face
 *  \param[in]  a_box   Box on which to compute flux
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::getFlux(FArrayBox&       a_flux,
                      const FArrayBox& a_WFace,
                      const int&       a_dir,
                      const Box&       a_box) const
{
  CH_TIME("ThermPhysics::getFlux");
  CH_assert(a_flux.contains(a_box));
  CH_assert(a_WFace.contains(a_box));

  const int numSpecies = CRDparam::g_numSpecies;
  const int primSpecStart = speciesPrimInterval().begin();
  const int consSpecStart = speciesConsInterval().begin();
  const int presIndx = pressureIndex();
  const int rhoIndx = densityIndex();
  const int tempIndx = temperatureIndex();
  const int engIndx = energyFluxIndex();
  const int velWIndx = velocityInterval().begin();
  const int velUIndx = vectorFluxInterval().begin();
  D_TERM(
    const int inorm = velWIndx + a_dir;,
    const int itan1 = velWIndx + ((a_dir + 1)%SpaceDim);,
    const int itan2 = velWIndx + ((a_dir + 2)%SpaceDim););

  D_TERM(
    const int inormc = velUIndx + a_dir;,
    const int itanc1 = velUIndx + ((a_dir + 1)%SpaceDim);,
    const int itanc2 = velUIndx + ((a_dir + 2)%SpaceDim););

  MD_ARRAY_RESTRICT(arrFlux, a_flux);
  MD_ARRAY_RESTRICT(arrWFace, a_WFace);
  // Vector of locations
  VECSTACKTEMPASSIGN(curL, int, m_interpN, 0);
  // Vector of coefficients
  VECSTACKTEMPASSIGN(Pvals, Real, m_interpN, 0.);
  MD_BOXLOOP(a_box, i)
    {
      Real rho = arrWFace[MD_IX(i, rhoIndx)];
      Real T = arrWFace[MD_IX(i, tempIndx)];
      // Assign proper indices and weighting to curL and Pvals
      lookupRefs(T, curL, Pvals);
      Real pres = arrWFace[MD_IX(i, presIndx)];
      D_TERM(
        Real uu = arrWFace[MD_IX(i, inorm)];,
        Real vv = arrWFace[MD_IX(i, itan1)];,
        Real ww = arrWFace[MD_IX(i, itan2)];);
      Real ke = 0.5*(D_TERM(uu*uu,+ vv*vv,+ ww*ww));
      // Calculate enthalpy
      Real sumfn = 0.;
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          int specIndxW = primSpecStart + sp;
          int specIndxU = consSpecStart + sp;
          Real spec = arrWFace[MD_IX(i, specIndxW)];
          arrFlux[MD_IX(i, specIndxU)] = rho*spec*uu;
          Real hval = m_H0[sp];
          for (int nt = 0; nt != m_interpN; ++nt)
            {
              hval += m_lookupH[sp][curL[nt]]*Pvals[nt];
            }
          sumfn += spec*hval;
        }
      // Calculate the continuity, momentum, and energy inertial flux
      arrFlux[MD_IX(i, rhoIndx)] = rho*uu;
      D_TERM(
        arrFlux[MD_IX(i, inormc)] = rho*uu*uu + pres;,
        arrFlux[MD_IX(i, itanc1)] = rho*uu*vv;,
        arrFlux[MD_IX(i, itanc2)] = rho*uu*ww;);
      arrFlux[MD_IX(i, engIndx)] = rho*uu*(ke + sumfn);
    }
  if (CRDparam::g_turbModelType)
    {
      m_turbModel->addTurbConvectiveFlux(a_flux, a_WFace, a_dir, a_box);
    }
  if (CRDparam::g_numTransportScalars > 0)
    {
      const int cVel = vectorFluxInterval().begin() + a_dir;
      MD_BOXLOOP(a_box, i)
        {
          for (int cTr = transportConsInterval().begin(),
                 cTrEnd = transportConsInterval().end() + 1;
               cTr != cTrEnd; ++cTr)
            {
              arrFlux[MD_IX(i, cTr)] =
                arrWFace[MD_IX(i, cVel)]*arrWFace[MD_IX(i, cTr)];
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute the solution to the Riemann problem
/** The CRD version of the Riemann problem does NOT apply boundary
 *  conditions. This is different from most Chombo versions.
 *  \param[out] a_WStar Face-centered solution to Riemann problem
 *  \param[in]  a_WLeft Left primitive state, on cells to left of each
 *                      face
 *  \param[in]  a_WRight
 *                      Right primitive state, on cells to right of
 *                      each face
 *  \param[in]  a_dir   Direction of the faces
 *  \param[in]  a_box   Face-centered box on which to compute a_WStar
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::riemann(FArrayBox&       a_WStar,
                      const FArrayBox& a_WLeft,
                      const FArrayBox& a_WRight,
                      const int&       a_dir,
                      const Box&       a_box) const
{
  CH_TIME("ThermPhysics::riemann");
  CRD::msg << CRD::fv4 << "ThermPhysics::riemann " << CRD::end;
  CH_assert(a_WStar.contains(a_box));
  // Get the numbers of relevant variables
  // Cast away "const" inputs so their boxes can be shifted left or right
  // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WLeft;
  FArrayBox& shiftWRight = (FArrayBox&)a_WRight;

  // Shift the left and right primitive variable boxes 1/2 cell so they are
  // face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);

  if (m_riemannSolver == RiemannSolver::Adaptive)
    { 
      AdaptiveRStarSolve func{};
      riemannSolution(func,
                      a_WStar,
                      shiftWLeft,
                      shiftWRight,
                      a_dir,
                      a_box);
    }
  else if (m_riemannSolver == RiemannSolver::Exact)
    { 
      ExactRStarSolve func{};
      riemannSolution(func,
                      a_WStar,
                      shiftWLeft,
                      shiftWRight,
                      a_dir,
                      a_box);
    }
  else
    {
      ApproxRStarSolve func{};
      riemannSolution(func,
                      a_WStar,
                      shiftWLeft,
                      shiftWRight,
                      a_dir,
                      a_box);
    }

  // Shift the left and right primitive variable boxes back to their original
  // position.
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}

/*--------------------------------------------------------------------*/
//  Compute the solution to the Riemann problem on the boundaries
/** 
 * \param[out]  a_WavgFace
 *                      Face-centered solution to Riemann problem
 *  \param[in]  a_WLeft Left primitive state, on cells to left of each
 *                      face
 *  \param[in]  a_WRight
 *                      Right primitive state, on cells to right of
 *                      each face
 *  \param[in]  a_unitNormalBasisFab
 *                      Unit normal basis Fab
 *  \param[in]  a_dir   Direction of the faces
 *  \param[in]  a_boundaryFaceBox
 *                      Face-centered box on which to compute a_WavgFace
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::riemannBC(FArrayBox&            a_WavgFace,
                        FArrayBox&            a_WLeft,
                        FArrayBox&            a_WRight,
                        const FArrayBox&      a_unitNormalBasisFab,
                        const int             a_dir,
                        const Side::LoHiSide& a_side,
                        const Box&            a_boundaryFaceBox) const
{
  CH_TIME("ThermPhysics::riemannBC");
  // Ensure we want to do limiting
  // Forward transform, the following is an alias
  FArrayBox velLeftFab(velocityInterval(), a_WLeft);
  PatchMappedFunc::forwardTransform(velLeftFab,
                                    a_unitNormalBasisFab,
                                    a_boundaryFaceBox);
  // Forward transform, the following is an alias
  FArrayBox velRightFab(velocityInterval(), a_WRight);
  PatchMappedFunc::forwardTransform(velRightFab,
                                    a_unitNormalBasisFab,
                                    a_boundaryFaceBox);
  if (m_riemannSolver == RiemannSolver::Adaptive)
    { 
      AdaptiveRStarSolve func{};
      riemannSolution(func,
                      a_WavgFace,
                      a_WLeft,
                      a_WRight,
                      a_dir,
                      a_boundaryFaceBox);
    }
  else if (m_riemannSolver == RiemannSolver::Exact)
    { 
      ExactRStarSolve func{};
      riemannSolution(func,
                      a_WavgFace,
                      a_WLeft,
                      a_WRight,
                      a_dir,
                      a_boundaryFaceBox);
    }
  else
    {
      ApproxRStarSolve func{};
      riemannSolution(func,
                      a_WavgFace,
                      a_WLeft,
                      a_WRight,
                      a_dir,
                      a_boundaryFaceBox);
    }
  // Reverse transform, the following is an alias
  FArrayBox velStarFab(velocityInterval(),
                       a_WavgFace);
  PatchMappedFunc::reverseTransform(velStarFab,
                                    a_unitNormalBasisFab,
                                    a_boundaryFaceBox);
  temperature(a_WavgFace, a_boundaryFaceBox);
}

/*--------------------------------------------------------------------*/
//  Compute primitive variables from conserved variables
/** \param[out] a_W     Computed primitive state
 *  \param[in]  a_U     Conservative state
 *  \param[in]  a_box   Where to compute
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::consToPrim(FArrayBox&       a_W,
                         const FArrayBox& a_U,
                         const Box&       a_box,
                         const FArrayBox& a_WOld) const
{
  CH_TIME("ThermPhysics::consToPrim");
  CRD::msg << CRD::fv4 << "ThermPhysics::consToPrim" << CRD::end;
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  const int consSpecStart = speciesConsInterval().begin();
  const int primSpecStart = speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  const int presIndx = pressureIndex();
  const int rhoIndx = densityIndex();
  const int engIndx = energyFluxIndex();
  const int velUIndx = vectorFluxInterval().begin();
  const int velWIndx = velocityInterval().begin();
  // Vector of locations
  VECSTACKTEMPASSIGN(curL, int, m_interpN, 0);
  // Vector of coefficients
  VECSTACKTEMPASSIGN(Pvals, Real, m_interpN, 0.);
  // Vector of mass fractions
  VECSTACKTEMPSIZE(spec, Real, numSpecies);
  // Indices into sorted mass fractions
  VECSTACKTEMPSIZE(sortcnj, int, numSpecies);
  MD_ARRAY_RESTRICT(arrW, a_W);
  MD_ARRAY_RESTRICT(arrU, a_U);
  MD_ARRAY_RESTRICT(arrWOld, a_WOld);
  MD_BOXLOOP(a_box, i)
    {
      // Sum the species to get rho used for bounding.  The idea here is to
      // only apply the lower bound so cn >= 0.
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int cSp = consSpecStart + sp;
          spec[sp] = std::max(0., arrU[MD_IX(i, cSp)]);
        }
      // Sort for a precise sum
      normalizePrimSpecies(NormalizeTypeNone,  // no normalization
                           false,              // already bounded
                           spec.data(),        // masses
                           1,                  // stride
                           sortcnj.data());    // indirection from sorting
      // Now sum to get rho for finding total mass per unit volume.  This is not
      // a density we trust, just a good one for finding cn from rhocn.
      Real rho = 0.;
      for (int j = 0; j != numSpecies; ++j)
        {
          const int sp = sortcnj[j];
          rho += spec[sp];
        }
      // Up to here, the only point was to get a very precise rho for finding
      // cn.
      // Now get the bounded species mass fractions (the bounds here should not
      // be necessary but are added just for certainty)
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          spec[sp] = std::max(0., std::min(1., spec[sp]/rho));
        }
      // Normalize the mass fractions so they sum to 1
      normalizePrimSpecies(NormalizeTypeRedistribute,
                           false,            // already bounded
                           spec.data(),      // mass fractions
                           1,                // stride
                           sortcnj.data());  // indirection from sorting
      // Compute properties by summing species in order from low cn to high
      // using indirection array
      Real Rval = 0.;
      Real heatofform = 0.;
      for (int j = 0; j != numSpecies; ++j)
        {
          const int sp = sortcnj[j];
          Rval += spec[sp]*m_Rn[sp];
          heatofform += spec[sp]*m_H0[sp];
        }

      // This is the actual density we trust.  The point density of each species
      // will then be spec[sp]*rho.
      rho = arrU[MD_IX(i, rhoIndx)];
      arrW[MD_IX(i, rhoIndx)] = rho;
      // Check to make sure uninitialized values are not being solved
      CH_assert(rho < 1.0E300);
      CH_assert(rho != 0.);

      // Velocities
      D_TERM(
        Real uval = arrU[MD_IX(i, velUIndx)]/rho;,
        Real vval = arrU[MD_IX(i, velUIndx+1)]/rho;,
        Real wval = arrU[MD_IX(i, velUIndx+2)]/rho;);
      D_TERM(
        arrW[MD_IX(i, velWIndx)] = uval;,
        arrW[MD_IX(i, velWIndx+1)] = vval;,
        arrW[MD_IX(i, velWIndx+2)] = wval;);

      // Compute temperature
      Real ke = 0.5*(D_TERM(uval*uval,+ vval*vval,+ wval*wval));
      Real energy = arrU[MD_IX(i, engIndx)]/rho - (heatofform + ke);
      Real Tmin = m_lookupLoT + 10.;
      Real Tmax = m_lookupHiT - 10.;
      const C2PFunc& f = C2PFunc(*this,energy,Rval,numSpecies,spec,Pvals,curL);
      int iterBrent = 0;
      int errorBrent = 0;
      Real Tval = RootSolver::BrentER(iterBrent,errorBrent,f,Tmin,Tmax);

      // Handle thermodynamic state that failed to find a temperature from
      // Brent solver
      if (errorBrent != 0 || Tval != Tval || Tval > Tmax || Tval < Tmin)
        {
          const int tempIndx = temperatureIndex();
          Tval = arrWOld[MD_IX(i, tempIndx)];
          // if (CRDparam::g_verbosity >= 3)
            {
              CRD::msg << CRD::fv1 << "-->: Brent failed with old temperature "
                "guess: " << Tval << CRD::end;
              CRD::msg << "On level " << CRDparam::g_level << " in cell "
                       << MD_GETIV(i) << CRD::end;
              CRD::msg << "Internal Energy: " << energy << " Rgas: " << Rval
                       << CRD::end;
              //           01234567-01234567890123-01234567
              CRD::msg << "  Name   Mass fraction    Rgas" << CRD::end;
              constexpr int c_bufsize = 33;
              char strbuf[c_bufsize];
              for (int j = 0; j != numSpecies; ++j)
                {
                  snprintf(strbuf, c_bufsize, "%-8s %14.7g %8.3f",
                           CRDparam::g_speciesNames[j].c_str(), spec[j],
                           m_Rn[j]);
                  CRD::msg << strbuf << CRD::end;
                }
              CRD::msg << "Brent Solver Error Value: " << errorBrent
                       << " Brent Iters: " << iterBrent
                       << " Tmin: " << Tmin << " Tmax: " << Tmax << CRD::end;
              CRD::msg << "Bad temperature in ThermPhysics::consToPrim temp: "
                       << Tval << CRD::end;
              CRD::msg << "Falling back to optimization solver to try and find "
                "a valid thermo state.  Calling L-BFGS-B solver "
                "ThermPhysics::findTempSpeciesFromFailedBrentSolve."
                       << CRD::warn;
            }

          findTempSpeciesFromFailedBrentSolve(Tval, spec, sortcnj, energy);
          CRD::msg << CRD::fv1 << "<--: Temperature found from L-BFGS-B solver "
            "is: " << Tval << CRD::end;

          // Compute properties by summing species in order from low cn to high
          // using indirection array
          Real Rval = 0.;
          for (int j = 0; j != numSpecies; ++j)
            {
              const int sp = sortcnj[j];
              Rval += spec[sp]*m_Rn[sp];
            }
        }

      for(int sp = 0; sp != numSpecies; ++sp)
        {
          int specIndxW = primSpecStart + sp;
          arrW[MD_IX(i, specIndxW)] = spec[sp];
        }

      const Real pres = Tval*Rval*rho;
      if (pres < CRDparam::g_smallp || std::isnan(pres) || std::isinf(pres))
        {
          CRD::msg << "Failure in ThermPhysics::ConsToPrim: invalid pressure: "
                   << pres << CRD::error;
        }
      arrW[MD_IX(i, presIndx)] = pres;
    } // Close BoxLoop

  // Wall model
  if (CRDparam::g_turbModelType)
    {
      m_turbModel->turbConsToPrim(a_W, a_U, a_box);
    }
  if (CRDparam::g_numTransportScalars > 0)
    {
      a_W.copy(a_box,
               transportPrimInterval(),
               a_box,
               a_U,
               transportConsInterval());
    }
}

/*--------------------------------------------------------------------*/
//  Compute conservative variables from primitive variables
/** \param[out] a_U     Computed conservative state
 *  \param[in]  a_W     Primitive state
 *  \param[in]  a_box   Where to compute
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::primToCons(FArrayBox&       a_U,
                         const FArrayBox& a_W,
                         const Box&       a_box) const
{
  CH_TIME("ThermPhysics::primToCons");
  // CH_assert(false);  // I don't think this is ever used.  If so, I would like
  //                    // to know when
  // We're using this now in CNSIBCBluffBodyCombustion in addSourceTerm
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  const int consSpecStart = speciesConsInterval().begin();
  const int primSpecStart = speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  const int presIndx = pressureIndex();
  const int rhoIndx = densityIndex();
  const int tempIndx = temperatureIndex();
  const int engIndx = energyFluxIndex();
  const int velIndx = velocityInterval().begin();
  const int momIndx = vectorFluxInterval().begin();

  VECSTACKTEMPSIZE(spec, Real, numSpecies);
  // Vector of locations
  VECSTACKTEMPASSIGN(curL, int, m_interpN, 0);
  // Vector of coefficients
  VECSTACKTEMPASSIGN(Pvals, Real, m_interpN, 0.);
  MD_BOXLOOP(a_box, i)
    {
      Real T = a_W[MD_IX(i, tempIndx)];
      // Assign proper indices and weighting to curL and Pvals
      lookupRefs(T, curL, Pvals);
      Real Rgas = 0.;
      // Calculate enthalpy
      Real sumfn = 0.;
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int specIndx = primSpecStart + sp;
          spec[sp] = a_W[MD_IX(i, specIndx)];
          Real hval = m_H0[sp];
          for (int nt = 0; nt != m_interpN; ++nt)
            {
              hval += m_lookupH[sp][curL[nt]]*Pvals[nt];
            }
          sumfn += spec[sp]*hval;
          Rgas += spec[sp]*m_Rn[sp];
        }
      Real rho = a_W[MD_IX(i, rhoIndx)];
      Real pres = a_W[MD_IX(i, presIndx)];
      a_U[MD_IX(i, rhoIndx)] = rho;
      D_TERM(
        Real uval = a_W[MD_IX(i, velIndx)];,
        Real vval = a_W[MD_IX(i, velIndx+1)];,
        Real wval = a_W[MD_IX(i, velIndx+2)];);
      D_TERM(
        a_U[MD_IX(i, momIndx)] = rho*uval;,
        a_U[MD_IX(i, momIndx+1)] = rho*vval;,
        a_U[MD_IX(i, momIndx+2)] = rho*wval;);
      Real ke = 0.5*(D_TERM(uval*uval, + vval*vval, + wval*wval));
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int specIndx = consSpecStart + sp;
          a_U[MD_IX(i, specIndx)] = rho*spec[sp]; 
        }
      a_U[MD_IX(i, engIndx)] = rho*ke + rho*sumfn - pres;
    }
}

/*--------------------------------------------------------------------*/
//  Set full average primitive state from native primitive state
/** The native state is the minimal state (e.g., pressure and density)
 *  from which other values are derived (.e.g, T = p/(rho*R)).
 *  This routine assumes that species mass fractions have already been
 *  bounded and normalized as required.
 *  \param[out] a_Wx    Location to store extra primitive state
 *  \param[in]  a_compWxBeg
 *                      Location in 'a_Wx' where extra components (T)
 *                      begin
 *  \param[in]  a_W     Native primitive state (rho, vel, pres)
 *  \param[in]  a_box   Where to compute the extra primitive state
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::extraPrimitiveState(FArrayBox&       a_Wx,
                                  const int        a_compWxBeg,
                                  const FArrayBox& a_Wp,
                                  const Box&       a_box) const
{
  CH_TIME("ThermPhysics::extraPrimitiveState");
  CH_assert(a_Wx.box().contains(a_box));
  CH_assert(a_Wp.box().contains(a_box));
  const int rhoIndx = densityIndex();
  const int presIndx = pressureIndex();
  const int primSpecStart = speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  VECSTACKTEMPSIZE(spec, Real, numSpecies);
  VECSTACKTEMPSIZE(sortcnj, int, numSpecies);
  MD_BOXLOOP(a_box, i)
    {
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int specIndx = primSpecStart + sp;
          spec[sp] = a_Wp[MD_IX(i, specIndx)];
        }
      normalizePrimSpecies(NormalizeTypeNone,  // no normalization
                           false,              // do not bound
                           spec.data(),        // mass fractions
                           1,                  // stride
                           sortcnj.data());    // indirection from sorting
      // Compute properties by summing species in order from low cn to high
      Real Rval = 0.;
      for (int j = 0; j != numSpecies; ++j)
        {
          const int sp = sortcnj[j];
          Rval += spec[sp]*m_Rn[sp];
        }
      Real rho = a_Wp[MD_IX(i, rhoIndx)];
      Real pres = a_Wp[MD_IX(i, presIndx)];
      Real T = pres/(rho*Rval);
      if (T < m_lookupLoT)
        {
          CRD::msg << "ThermPhysics::extraPrimitiveState: index " << MD_GETIV(i)
                   << " on level " << CRDparam::g_level << ", temperature "
            "adjusted from " << T << " to " << m_lookupLoT << " to permit "
            "table lookup (but violates perfect gas law!)" << CRD::warn;
          T = m_lookupLoT;
        }
      a_Wx[MD_IX(i, a_compWxBeg)] = T;
    }
}

/*--------------------------------------------------------------------*/
//  Set full average primitive state from native primitive state
//  in the same FAB
/** The native state is the minimal state (e.g., pressure and density)
 *  from which other values are derived (.e.g, T = p/(rho*R)).  For
 *  CNS, the full state adds temperature.
 *  \param[out] a_W     Primitive state variables to solve for extra
 *  \param[in]  a_box   Where to compute the extra primitive state
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::extraPrimitiveState(FArrayBox& a_W,
                                  const Box& a_box) const
{
  extraPrimitiveState(a_W,
                      temperatureIndex(),
                      a_W,
                      a_box);
}

/*--------------------------------------------------------------------*/
//  Compute the temperature from primary primitive variables
/** \param[in]  a_W     Computed primary primitive state
 *  \param[out] a_W     Updated with temperature
 *  \param[in]  a_U     Conservative variables
 *  \param[in]  a_box   Where to compute the temperature
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::temperature(FArrayBox& a_W,
                          const Box& a_box) const
{
  CH_TIME("ThermPhysics::temperature");
  CH_assert(a_W.box().contains(a_box));
  const int rhoIndx = densityIndex();
  const int presIndx = pressureIndex();
  const int tempIndx = temperatureIndex();
  const int primSpecStart = speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  MD_BOXLOOP(a_box, i)
    {
      Real Rval = 0.;
      Real rho = a_W[MD_IX(i, rhoIndx)];
      Real pres = a_W[MD_IX(i, presIndx)];
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int specIndx = primSpecStart + sp;
          Rval += m_Rn[sp]*a_W[MD_IX(i, specIndx)];
        }
      a_W[MD_IX(i, tempIndx)] = pres/(rho*Rval);
    }
}

/*--------------------------------------------------------------------*/
//  Compute the pressure from temperature and density
/** \param[in]  a_W     Computed primary primitive state
 *  \param[out] a_W     Updated with temperature
 *  \param[in]  a_U     Conservative variables
 *  \param[in]  a_box   Where to compute the pressure
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::pressure(FArrayBox& a_W,
                       const Box& a_box) const
{
  CH_TIME("ThermPhysics::pressure");
  CH_assert(a_W.box().contains(a_box));
  const int rhoIndx = densityIndex();
  const int presIndx = pressureIndex();
  const int tempIndx = temperatureIndex();
  const int primSpecStart = speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  MD_BOXLOOP(a_box, i)
    {
      Real Rval = 0.;
      Real rho = a_W[MD_IX(i, rhoIndx)];
      Real temp = a_W[MD_IX(i, tempIndx)];
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int specIndx = primSpecStart + sp;
          Rval += m_Rn[sp]*a_W[MD_IX(i, specIndx)];
        }
      a_W[MD_IX(i, presIndx)] = rho*Rval*temp;
    }
}

/*--------------------------------------------------------------------*/
//  Compute the density from temperature and pressure
/** \param[in]  a_W     Computed primary primitive state
 *  \param[out] a_W     Updated with temperature
 *  \param[in]  a_U     Conservative variables
 *  \param[in]  a_box   Where to compute density
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::density(FArrayBox& a_W,
                      const Box& a_box) const
{
  CH_TIME("ThermPhysics::density");
  CH_assert(a_W.box().contains(a_box));
  const int rhoIndx = densityIndex();
  const int presIndx = pressureIndex();
  const int tempIndx = temperatureIndex();
  const int primSpecStart = speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  MD_BOXLOOP(a_box, i)
    {
      Real Rval = 0.;
      Real pres = a_W[MD_IX(i, presIndx)];
      Real temp = a_W[MD_IX(i, tempIndx)];
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int specIndx = primSpecStart + sp;
          Rval += m_Rn[sp]*a_W[MD_IX(i, specIndx)];
        }
      a_W[MD_IX(i, rhoIndx)] = pres/(temp*Rval);
    }
}

/*--------------------------------------------------------------------*/
//  Solves for the linearly related primitive variables
/** \param[out] a_W    FAB of high-order primitive variables
 *  \param[in]  a_U    FAB of conservative variables
 *  \param[in]  a_box  Box to solve over
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::intermediateConsToPrim(FArrayBox&       a_W,
                                     const FArrayBox& a_U,
                                     const Box&       a_box,
                                     const FArrayBox& a_WOld) const
{
  CH_TIME("ThermPhysics::intermediateConsToPrim");
  constexpr int cRho   = densityIndex();
  constexpr int cMom0  = CRDPhysics::velocityInterval().begin();

  MD_BOXLOOP(a_box, i)
    {
      const Real rhoInv = 1./a_U[MD_IX(i, cRho)];
      a_W[MD_IX(i, cRho)] = a_U[MD_IX(i, cRho)];
      D_TERM(a_W[MD_IX(i, cMom0 + 0)] = a_U[MD_IX(i, cMom0 + 0)]*rhoInv;,
             a_W[MD_IX(i, cMom0 + 1)] = a_U[MD_IX(i, cMom0 + 1)]*rhoInv;,
             a_W[MD_IX(i, cMom0 + 2)] = a_U[MD_IX(i, cMom0 + 2)]*rhoInv;)
    }

  constexpr Real minWSp = 0.;
  constexpr Real maxWSp = 1.;
  const int numSp = CRDparam::g_numSpecies;
  VECSTACKTEMPSIZE(species, Real, numSp);

  MD_BOXLOOP(a_box, i)
    {
      const Real rhoInv = 1./a_U[MD_IX(i, cRho)];
      Real sumWSp = 0.;
      for (int idxSp = numSp, cUSp = speciesConsInterval().end(); idxSp--;
           --cUSp)
        {
          // Species concentration with physical bounds imposed
          const Real WSp =
            std::max(minWSp, std::min(maxWSp, a_U[MD_IX(i, cUSp)]*rhoInv));
          species[idxSp] = WSp;
          sumWSp += WSp;
        }
      for (int idxSp = numSp, cWSp = speciesPrimInterval().end(); idxSp--;
           --cWSp)
        {
          a_W[MD_IX(i, cWSp)] = species[idxSp]/sumWSp;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Calculate the intermediate temperature and the derivative
//  of f(T) = <e> 
/** \param[out] a_Tbar FAB of the intermediate temperature
 *  \param[out] a_dFdTbar
 *                     FAB of d <e>/d Tbar
 *  \param[in]  a_WcellPntFab
 *                     FAB of cell-centered primitive variables rho, u, cn
 *  \param[in]  a_box  Box to solve over
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::intermediateTandDiffF(FArrayBox&       a_Tbar,
                                    FArrayBox&       a_dFdTbar,
                                    const FArrayBox& a_WcellPntFab,
                                    const FArrayBox& a_engCellAvgFab,
                                    const Box&       a_box) const
{
  // const int primSpecStart = speciesPrimInterval().begin();
  // const int numSpecies = CRDparam::g_numSpecies;
  // std::vector<Real> spec(numSpecies,0.);
  // MD_ARRAY_RESTRICT(arrW, a_WcellPntFab);
  // MD_ARRAY_RESTRICT(arrTBar, a_Tbar);
  // MD_ARRAY_RESTRICT(arrDFDT, a_dFdTbar);
  // MD_ARRAY_RESTRICT(arrEng, a_engCellAvgFab);
  // MD_BOXLOOP(a_box, i)
  //   {
  //     Real Rgas = 0.;
  //     Real heatofform = 0.;
  //     for (int sp = 0; sp != numSpecies; ++sp)
  //       {
  //         int wSpecIndx = primSpecStart + sp;
  //         spec[sp] = arrW[MD_IX(i, wSpecIndx)];
  //         Rgas += spec[sp]*m_Rn[sp];
  //         heatofform += spec[sp]*m_H0[sp];
  //       }
  //     Real energy = arrEng[MD_IX(i, 0)] - heatofform;
  //     Real Tmin = m_lookupLoT + 10.;
  //     Real Tmax = m_lookupHiT - 10.;
  //     const C2PFunc& f = C2PFunc(*this,energy,Rgas,spec);
  //     int iterBrent = 0;
  //     int errorBrent = 0;
  //     Real Tbar = RootSolver::BrentER(iterBrent,errorBrent,f,Tmin,Tmax);
  //     Real Cp = calcCp(Tbar, spec);
  //     Real dFdTbar = Cp - Rgas;
  //     arrTBar[MD_IX(i, 0)] = Tbar;
  //     arrDFDT[MD_IX(i, 0)] = dFdTbar;
  //   }
}

/*--------------------------------------------------------------------*/
//  Solve for the thermal conductivity and dynamic viscosity
/** \param[in]  a_box   Cell box.  Flux needs to be computed on
 *                      surrounding faces.
 *  \param[out] a_muFab
 *                      FAB of dynamic viscosity
 *  \param[out] a_kappaFab
 *                      FAB of thermal conductivity
 *  \param[in]  a_WfacePntFab
 *                      FAB containing the face-centered primitive variables
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::calcCoeffKappaMu(const Box&       a_box,
                               FArrayBox&       a_muFab,
                               FArrayBox&       a_kappaFab,
                               const FArrayBox& a_WfacePntFab) const
{
  CH_TIME("ThermPhysics::calcCoeffKappaMu");
  CRD::msg << CRD::fv4 << "ThermPhysics::calcCoeffKappaMu" << CRD::end;
  CH_assert(a_muFab.contains(a_box));
  CH_assert(a_WfacePntFab.contains(a_box));
  // If mu and kappa are set to constants, we use the constants
  if (CRDparam::g_mu >= 0.)
    {
      a_kappaFab.setVal(CRDparam::g_K);
      a_muFab.setVal(CRDparam::g_mu);
      return;
    }
  // Otherwise, solve for the values using temperature
  const int numSpecies = CRDparam::g_numSpecies;
  const int primSpecStart = speciesPrimInterval().begin();
  const int tempIndx = temperatureIndex();
  // Vector of species mass fractions
  VECSTACKTEMPSIZE(spec, Real, numSpecies);
  // Vector of mu for each species
  VECSTACKTEMPASSIGN(mui, Real, numSpecies, 0.);
  // Vector of kappa for each species
  VECSTACKTEMPASSIGN(kappai, Real, numSpecies, 0.);
  MD_ARRAY_RESTRICT(arrW, a_WfacePntFab);
  MD_ARRAY_RESTRICT(arrMu, a_muFab);
  MD_ARRAY_RESTRICT(arrKappa, a_kappaFab);
  // Vector of locations
  VECSTACKTEMPASSIGN(curL, int, m_interpN, 0);
  // Vector of coefficients
  VECSTACKTEMPASSIGN(Pvals, Real, m_interpN, 0.);
  MD_BOXLOOP(a_box, i)
    {
      const Real T = std::max(m_lookupLoT + 150., arrW[MD_IX(i, tempIndx)]);
      // Molecular weight of mixture
      Real mwMix = 0.;
      // Fill vector of locations and coefficients for lookup table
      lookupRefs(T, curL, Pvals);
      for (int comp = 0; comp != numSpecies; ++comp)
        {
          const int wComp = comp + primSpecStart;
          // Extract mass fractions
          spec[comp] = arrW[MD_IX(i, wComp)];
          mwMix += spec[comp]/m_molMass[comp];
          Real muval = 0.;
          Real kappaval = 0.;
          for (int nt = 0; nt != m_interpN; ++nt)
            {
              muval = muval + m_lookupMu[comp][curL[nt]]*Pvals[nt];
              kappaval = kappaval + m_lookupKappa[comp][curL[nt]]*Pvals[nt];
            }
          mui[comp] = muval;
          kappai[comp] = kappaval;
        }
      mwMix = 1./mwMix;
      Real kappaValD = 0.;
      Real kappaValN = 0.;
      Real muValD = 0.;
      Real muValN = 0.;
      for (int comp = 0; comp != numSpecies; ++comp)
        {
          // Molar fraction
          const Real xi = spec[comp]*mwMix/m_molMass[comp];
          kappaValD += kappai[comp]*xi;
          muValD += mui[comp]*xi;
          kappaValN += xi/kappai[comp];
          muValN += xi/mui[comp];
        }
      Real muValF = 0.5*(muValD + 1./muValN);
      arrMu[MD_IX(i, 0)] = muValF;
      arrKappa[MD_IX(i, 0)] = 0.5*(kappaValD + 1./kappaValN);
    }
  return;
}

/*--------------------------------------------------------------------*/
//  Solves for the mass diffusivity
/** \param[in]  a_box   Box of faces to solve over
 *  \param[in]  a_JnfacePntFab
 *                      FAB of the face-centered mapped gradients
 *  \param[out] a_JnfacePntFab
 *                      FAB of gradients multiplied by D_n
 *  \param[in]  a_muFab
 *                      FAB of the dynamic viscosity on each face
 *  \param[in]  a_kappaFab
 *                      FAB of the thermal conductivity values
 *  \param[in]  a_WfacePntFab
 *                      FAB containing face-centered primitive variables
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::speciesDiffHeatFlux(const Box&       a_box,
                                  FArrayBox&       a_JnfacePntFab,
                                  FArrayBox&       a_energyFab,
                                  const FArrayBox& a_muFab,
                                  const FArrayBox& a_kappaFab,
                                  const FArrayBox& a_WfacePntFab) const
{
  CH_TIME("ThermPhysics::speciesDiffHeatFlux");
  CRD::msg << CRD::fv4 << "ThermPhysics::speciesDiffusivity" << CRD::end;
  const int numSpecies = CRDparam::g_numSpecies;
  CH_assert(a_JnfacePntFab.nComp() == SpaceDim*numSpecies);
  const int tempIndx = temperatureIndex();
  const int primSpecStart = speciesPrimInterval().begin();
  // Vector of locations
  VECSTACKTEMPASSIGN(curL, int, m_interpN, 0);
  // Vector of coefficients
  VECSTACKTEMPASSIGN(Pvals, Real, m_interpN, 0.);
  // If Lewis # is not specified in the input file, then
  // rho*D_n = mu/Sc
  if (m_schmidtNum > 0.)
    {
      MD_BOXLOOP(a_box, i)
        {
          Real mu = a_muFab[MD_IX(i, 0)];
          Real T = a_WfacePntFab[MD_IX(i, tempIndx)];
          Real Dn = mu/m_schmidtNum;
          lookupRefs(T, curL, Pvals);
          for (int sp = 0; sp != numSpecies; ++sp)
            {
              Real hval = m_H0[sp];
              for (int nt = 0; nt != m_interpN; ++nt)
                {
                  hval += m_lookupH[sp][curL[nt]]*Pvals[nt];
                }
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  const int jnComp = SpaceDim*sp + dir;
                  Real Jn = Dn*a_JnfacePntFab[MD_IX(i, jnComp)];
                  a_JnfacePntFab[MD_IX(i, jnComp)] = Jn;
                  a_energyFab[MD_IX(i, dir)] += hval*Jn;
                }
            }
        }
    }
  // If Lewis # is specified in input file, then
  // rho*D_n = kappa/(c_p*Le_i)
  else
    {
      // Vector of enthalpies
      VECSTACKTEMPSIZE(hn, Real, numSpecies);
      // Vector of mass fractions
      VECSTACKTEMPSIZE(spec, Real, numSpecies);
      MD_BOXLOOP(a_box, i)
        {
          Real kappa = a_kappaFab[MD_IX(i, 0)];
          Real T = a_WfacePntFab[MD_IX(i, tempIndx)];
          lookupRefs(T, curL, Pvals);
          for (int sp = 0; sp != numSpecies; ++sp)
            {
              const int specIndx = primSpecStart + sp;
              Real hval = m_H0[sp];
              for (int nt = 0; nt != m_interpN; ++nt)
                {
                  hval += m_lookupH[sp][curL[nt]]*Pvals[nt];
                }
              hn[sp] = hval;
              spec[sp] = a_WfacePntFab[MD_IX(i, specIndx)];
            }
          Real cpVal = cp(T, spec.data());
          Real Dn = kappa/(m_lewisNum*cpVal);
          for (int sp = 0; sp != numSpecies; ++sp)
            {
              for (int dir = 0; dir != SpaceDim; ++dir)
                {
                  const int jnComp = SpaceDim*sp + dir;
                  Real Jn = Dn*a_JnfacePntFab[MD_IX(i, jnComp)];
                  a_JnfacePntFab[MD_IX(i, jnComp)] = Jn;
                  a_energyFab[MD_IX(i, dir)] += hn[sp]*Jn;
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Normalize the species mass fractions (primitive state)
/** This general implementation supports the approaches of
 *  redistribution and lumping in an inert species.  One can also
 *  first bound the fractions and sort the fractions before summing
 *  for increased accuracy.
 *  \param[in]  a_normalizeType
 *                      Type of normalization:
 *                        NormalizeTypeRedistribute: Redistribute all
 *                          species so they sum to 1
 *                        NormalizeTypeLumpInert: Lump any errors into
 *                          the inert species so al sum to 1
 *  \param[in]  a_bound Bound all mass fractions between 0 and 1
 *                      before normalizing
 *  \param[in]  a_cn    Species mass fractions
 *  \param[out] a_cn    Normalized species mass fractions
 *  \param[in]  a_cnStride
 *                      Stride through species mass fractions
 *  \param[out] a_sortcnj
 *                      An indirection array of size >=
 *                      CRDparam::g_numSpecies that permits sorting
 *                      of species before summing.  Will be filled
 *                      with an indirections from lowest to highest
 *                      mass fraction from _before_ normalization.
 *                      The stride through this array must be 1.
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::normalizePrimSpecies(const SpeciesNormalizeType a_normalizeType,
                                   const bool                 a_bound,
                                   Real *const                a_cn,
                                   const std::ptrdiff_t       a_cnStride,
                                   int*                       a_sortcnj) const
{
  const int numSpecies = CRDparam::g_numSpecies;

  // Apply bounds if requested
  if (a_bound)
    {
      Real* p_cn = a_cn;
      for (int n = numSpecies; n--;)
        {
          *p_cn = std::max(0., std::min(1., *p_cn));
          p_cn += a_cnStride;
        }
    }

  // Sort if indirection array is present
  if (a_sortcnj != const_cast<int*>(m_intSeq.data()))
    {
      if (a_sortcnj)
        {
          for (int j = 0; j != numSpecies; ++j)
            {
              a_sortcnj[j] = j;
            }
          Sort::insertion(
            numSpecies,
            a_sortcnj,
            Sort::CmpLessAbsIndexStride<Real>(a_cn, a_cnStride));
        }
      else
        {
          a_sortcnj = const_cast<int*>(m_intSeq.data());
        }
    }

  if (a_normalizeType == NormalizeTypeNone)
    {
      return;
    }

  // Normalize
  Real sumcj = 0.;
  for (int j = 0; j != numSpecies; ++j)
    {
      const int sp = a_sortcnj[j];
      sumcj += a_cn[sp*a_cnStride];
    }
  switch (a_normalizeType)
    {
    case NormalizeTypeRedistribute:
      for (int j = 0; j != numSpecies; ++j)
        {
          const int sp = a_sortcnj[j];
          a_cn[sp*a_cnStride] /= sumcj;
        }
      break;
    case NormalizeTypeLumpInert:
    default:
      a_cn[m_N2idx*a_cnStride] += 1. - sumcj;
      break;
    }

  // The differences should not be large enough to necessitate a resort
}

/*--------------------------------------------------------------------*/
//  Normalize the species mass fractions (primitive state)
/** This general implementation supports the approaches of
 *  redistribution and lumping in an inert species.  One can also
 *  first bound the fractions and sort the fractions before summing
 *  for increased accuracy.
 *  \param[in]  a_normalizeType
 *                      Type of normalization:
 *                        NormalizeTypeRedistribute: Redistribute all
 *                          species so they sum to 1
 *                        NormalizeTypeLumpInert: Lump any errors into
 *                          the inert species so al sum to 1
 *  \param[in]  a_bound Bound all mass fractions between 0 and 1
 *                      before normalizing
 *  \param[in]  a_sortcnj
 *                      Whether the mass fractions should be sorted
 *                      from low to high before summing
 *  \param[in]  a_box   Cells to operate on
 *  \param[in]  a_W     Primitive state with species mass fractions
 *  \param[out] a_W     Normalized species mass fractions
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::normalizePrimSpecies(const SpeciesNormalizeType a_normalizeType,
                                   const bool                 a_bound,
                                   const bool                 a_sortcnj,
                                   const Box&                 a_box,
                                   FArrayBox&                 a_W) const
{
  CH_assert(a_W.contains(a_box));
  CH_assert(speciesPrimInterval().end() < a_W.nComp());

  const int compPrimSpecBeg = speciesPrimInterval().begin();
  // Indices into sorted mass fractions
  VECSTACKTEMPSIZE(sortcnjVec, int, CRDparam::g_numSpecies);
  int *const sortcnj = (a_sortcnj) ? sortcnjVec.data() : nullptr;
  const std::ptrdiff_t stride = a_W.strideComp();
  MD_BOXLOOP(a_box, i)
    {
      normalizePrimSpecies(a_normalizeType,
                           a_bound,
                           &a_W[MD_IX(i, compPrimSpecBeg)],
                           stride,
                           sortcnj);
    }
}

/*--------------------------------------------------------------------*/
//  Perform the species correction
/** \param[in]  a_box   Cell box.  Flux needs to be computed on
 *                      surrounding faces.
 *  \param[in]  a_U     FAB of the cell-averaged conservative variables
 *  \param[out] a_U     Species corrected cell-averaged conservative variables
 *  \param[in]  a_J     FAB of the values of mapping variable J
 *  \param[in]  a_domain
 *                      Problem domain
 *  \param[in]  a_dx    Grid spacing
 *  \param[in]  a_tolNeg
 *                      This is the maximum of Jrho that is allowed
 *                      to be negative.  I.e.,
 *                      tolJrhoNeg = -a_tolNeg*Jrho.  Note that this
 *                      must be >= 0.  The default is 0.01 or 1 %.
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::speciesCorrection(const Box&           a_box,
                                FArrayBox&           a_JU,
                                const ProblemDomain& a_domain,
                                const RealVect&      a_dx,
                                const Real           a_tolNeg) const
{
  CH_TIME("ThermPhysics::speciesCorrection");
  CRD::msg << CRD::fv4 << "ThermPhysics::speciesCorrection" << CRD::end;
  CH_assert(a_tolNeg >= 0.);
  const int numSpecies = CRDparam::g_numSpecies;
  // Density component
  const int rhoIndx = densityIndex();
  // First component for the species conservative variables
  const int consSpecStart = speciesConsInterval().begin();
#if 1
  // Vector of mass fractions
  VECSTACKTEMPSIZE(spec, Real, numSpecies);

/* rho is our primary trusted density and we think of rhocn as designating
   fractions of rho.  In this correction, we want \sum_rhocn = rho.  We can
   achieve that without alterning the definition of rhocn as a fraction.
   However, we also want to restrict or bound the amount of negative rhocn.
   Some is permitted because these values are averages of near zero quantities.
   We restrict the negative fraction to 1 % of rho.  Applying bounds can alter
   the remaining fractions.

   The main thought process here is that we don't know what else prevents
   some species value from going to a ridiculous negative value.
*/

  const int zeroScalingCoef = (a_tolNeg == (Real)0) ? -1 : 0;
  const int zeroScalingMult = !zeroScalingCoef;  // Note: 0 or 1
  MD_BOXLOOP(a_box, i)
    {
      // The trusted Jrho
      const Real Jrho = a_JU[MD_IX(i, rhoIndx)];

      Real JrhoPos = 0.;
      Real JrhoNeg = 0.;
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          spec[sp] = a_JU[MD_IX(i, consSpecStart + sp)];
          JrhoPos += std::max(0., spec[sp]);
          JrhoNeg += std::min(0., spec[sp]);
        }
      const Real tolJrhoNeg = -a_tolNeg*Jrho;
      const Real allowedJrhoNeg = std::max(JrhoNeg, tolJrhoNeg);
      // Always < 0.  Either:
      // -eps            if (a_tolNeg == 0.)
      // eps*tolJrhoNeg  otherwise
      const Real zeroJrhoNeg = std::numeric_limits<Real>::epsilon()*
        (zeroScalingCoef + zeroScalingMult*tolJrhoNeg);

      const Real JrhoSum = JrhoPos + allowedJrhoNeg;
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          if (spec[sp] < 0.)
            {
              spec[sp] = allowedJrhoNeg*(spec[sp]/(JrhoNeg + zeroJrhoNeg));
            }
          spec[sp] /= JrhoSum;
        }
      // Sum over spec should now be 1.

      Real sumJrhocn = 0.;
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const Real Jrhocn = spec[sp]*Jrho;
          sumJrhocn += Jrhocn;
          a_JU[MD_IX(i, consSpecStart + sp)] = Jrhocn;
        }
      // Final safety correction
      a_JU[MD_IX(i, consSpecStart + m_N2idx)] += Jrho - sumJrhocn;
    }
#else
  const int N2ConsIndx = consSpecStart+m_N2idx;
  VECSTACKTEMPSIZE(rspec, Real, numSpecies);
  MD_BOXLOOP(a_box, i)
    {
      Real rho = a_JU[MD_IX(i, rhoIndx)];
      Real sumWithoutN2 = 0.;
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int specIndx = consSpecStart + sp;
          rspec[sp] = std::min(rho, std::max(0., a_JU[MD_IX(i, specIndx)]));
          if(sp != m_N2idx)
            {
              sumWithoutN2 += rspec[sp];
              a_JU[MD_IX(i, specIndx)] = rspec[sp];
            }
        }
      a_JU[MD_IX(i,N2ConsIndx)] = rho - sumWithoutN2;
    }
#endif
}

/*--------------------------------------------------------------------*/
//  Solves for rho*omega, reaction rate and finds the max Jacobian of
//  the reaction source terms with respect to the conservative variables
/** \param[in]  a_box   Cell box to solve the source term
 *  \param[out] a_RCTcellPntFab
 *                      FAB to solve cell-centered reaction rate source
 *  \param[out] a_invDtFab
 *                      Inverse of the sum of the time step sizes
 *  \param[in]  a_WcellPntFab
 *                      FAB containing cell-centered primitive variables
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[out] a_minChemDt
 *                      Mininum chemical time step
 *  \param[out] a_minChemDtCell
 *                      Cell with minimum chemical time step
 *  \return     cnlim   Species associated with the smallest time step
 *//*-----------------------------------------------------------------*/

int
ThermPhysics::addReactionSource(const Box&       a_box,
                                FArrayBox&       a_RCTcellPntFab,
                                FArrayBox&       a_invDtFab,
                                const FArrayBox& a_WcellPntFab,
                                const Real       a_time,
                                const int        a_level,
                                Real&            a_minChemDt,
                                IntVect&         a_minChemDtCell) const
{
  CH_TIME("ThermPhysics::addReactionSource");
  CRD::msg << CRD::fv4 << "ThermPhysics::addReactionSource" << CRD::end;
  CH_assert(a_invDtFab.contains(a_box));
  // Limiting species for reaction time step
  int cnlim = -1;
  // First primitive index for species 
  const int primSpecStart = speciesPrimInterval().begin();
  // Number of species
  const int numSpecies = CRDparam::g_numSpecies;
  // Number of reactions
  const int numReactions = CRDparam::g_numReactions;
  // Get density index
  const int rhoIndx = densityIndex();
  // Get temperature index
  const int tempIndx = temperatureIndex();
  // Cutoff temperature
  const Real Tcutoff = m_Tcutoff;
  Real refFromBase = CRDparam::g_refFromBase[a_level];
  Real maxARKLevel = CRDparam::g_ARKmaxLevel;
  Real chemDtScale = CRDparam::g_chemicalDtScale;
  Real useChemDt = 1.;
  if(a_level <= maxARKLevel && CRDparam::g_additiveRK)
    {
      useChemDt = 0.;
    }
  a_RCTcellPntFab.setVal(0.);
  VECSTACKTEMPSIZE(mmass1, Real, numSpecies);
  VECSTACKTEMPSIZE(mmass2, Real, numSpecies);
  VECSTACKTEMPSIZE(spec, Real, numSpecies);
  VECSTACKTEMPSIZE(destrate, Real, numSpecies);

  VECSTACKTEMPASSIGN(Pvals, Real, m_interpN, 0.);
  VECSTACKTEMPASSIGN(curL, int, m_interpN, 0);
  for (int sp = 0; sp != numSpecies; ++sp)
    {
      mmass1[sp] = m_preAmodUnits/m_molMass[sp];
      mmass2[sp] = m_molMass[sp]/m_preAmodUnits;
    }
  MD_ARRAY_RESTRICT(arrRCT, a_RCTcellPntFab);
  MD_ARRAY_RESTRICT(arrW, a_WcellPntFab);
  MD_ARRAY_RESTRICT(arrInvDt, a_invDtFab);
  MD_BOXLOOP(a_box, i)
    {
      Real rho = arrW[MD_IX(i, rhoIndx)];
      Real T = arrW[MD_IX(i, tempIndx)];
      if (T < Tcutoff)
        {
          continue;
        }
      const Real Tinv = 1./T;
      lookupRefs(T, curL, Pvals);
      // Find molar concentration and reset destrate
      Real maxConc = std::numeric_limits<Real>::min();
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          int specIndx = primSpecStart + sp;
          const Real conc = std::min(1.,std::max(0.,arrW[MD_IX(i, specIndx)]))
            *mmass1[sp]*rho;
          maxConc = std::max(maxConc, conc);
          spec[sp] = conc;
          destrate[sp] = 0.;
        }
      // An effective zero values is set based on maxConc
      const Real zeroConc = maxConc*std::numeric_limits<Real>::epsilon();
      // Current third-body reaction index
      int curTB = 0;
      // Current pressure-dependent reaction index
      int curPR = 0;
      // Current arbitrary order reaction index
      int curARB = 0;
      // Loop over reactions
      for (int rctcomp = 0; rctcomp != numReactions; ++rctcomp)
        {
          // Number of reactants for this reaction
          const int cnumr = m_numRects[rctcomp];
          // Number of products for this reaction
          const int cnump = m_numProds[rctcomp];
          // Current index for nup and nupp, m_maxSpecPerReact
          // is the stride
          const int nucomp = m_maxSpecPerReact*rctcomp;
          // Solve for $k_{f,r}$
          Real kfwd = 0.;
          Real kbkwd = 0.;
          Real keq = 0.;
          for (int nt = 0; nt != m_interpN; ++nt)
            {
              kfwd += m_lookupKfwd[rctcomp][curL[nt]]*Pvals[nt];
              kbkwd += m_lookupKbkwd[rctcomp][curL[nt]]*Pvals[nt];
              keq += m_lookupKeq[rctcomp][curL[nt]]*Pvals[nt];
            }
          Real sumtbval = 1.;
          // If a three-body reaction
          if (m_thirdBody[rctcomp] == 1)
            {
              sumtbval = 0.;
              for (int sp = 0; sp != numSpecies; ++sp)
                {
                  const int tbcomp = sp*m_numTBReact + curTB;
                  sumtbval += m_tbAlpha[tbcomp]*spec[sp];
                }
              curTB++;
            }
          else if (m_thirdBody[rctcomp] == -1)
            {
              curTB++;
            }
          // If pressure dependent reaction
          if (m_numPRReact > 0)
            {
              int presDep = m_presDep[rctcomp];
              if (presDep != 0)
                {
                  Real k0, kinfty;
                  Real k1 = kfwd;
                  Real k2 = m_PRpreAF[curPR]*std::pow(T, m_PRbetai[curPR])*
                    std::exp(-m_PREAR[curPR]*Tinv);
                  if (presDep > 0)
                    {
                      k0 = k1;
                      kinfty = k2;
                    }
                  else
                    {
                      k0 = k2;
                      kinfty = k1;
                    }
                  Real Pr = k0*sumtbval/kinfty;
                  Real F = m_Fc;
                  if (presDep == 2)
                    {
                      F = TroeSolution(Pr, T, curPR);
                    }
                  else if (presDep == 3)
                    {
                      F = SRISolution(Pr, T, curPR);
                    }
                  if (presDep > 0)
                    {
                      kfwd = k0*(1./(Pr + 1.))*F;
                    }
                  else
                    {
                      kfwd = kinfty*(Pr/(Pr + 1.))*F;
                    }
                  curPR++;
                  kbkwd = 0;
                  // For the irreversible, third-body, and pressure-dependent
                  // reaction, the backward reaction rate should only be set to
                  // be zero. Otherwise, one would get the inf value on kbkwd
                  // because keq is zero. This was observed in Z66 propane-air
                  // kinetics.
                  if (m_revReact[rctcomp] >= 0)
                    {
                      kbkwd = kfwd/keq;
                    }
                }
            } 
          // ffwd is the product of [x_n]^nup
          Real ffwd = 1.;
          // fback is the pruduct of [x_n]^nupp
          Real fback = 1.;
          if (std::abs(m_revReact[rctcomp]) == 2)
            {
              for (int curs = 0; curs != cnumr; ++curs)
                {
                  const int offcomp = nucomp + curs;
                  const int spcomp1 = m_refNup[offcomp];
                  const int arbcomp = curARB*m_maxSpecPerReact + curs;
                  ffwd *= std::pow(spec[spcomp1],m_ford[arbcomp]);
                }
              // If reaction is reversible
              if (m_revReact[rctcomp] >= 0)
                {
                  for (int curs = 0; curs != cnump; ++curs)
                    {
                      const int offcomp = nucomp + curs;
                      const int spcomp1 = m_refNupp[offcomp];
                      const int arbcomp = curARB*m_maxSpecPerReact + curs;
                      fback *= std::pow(spec[spcomp1],m_rord[arbcomp]);
                    }
                  kbkwd *= sumtbval*fback;
                }
              ++curARB;
            }
          else
            {
              for (int offcomp = nucomp; offcomp != nucomp+cnumr; ++offcomp)
                {
                  const int spcomp1 = m_refNup[offcomp];
                  ffwd *= std::pow(spec[spcomp1],m_nup[offcomp]);
                }
              // If reaction is reversible
              if (m_revReact[rctcomp] >= 0)
                {
                  for (int offcomp = nucomp; offcomp != nucomp+cnump; ++offcomp)
                    {
                      const int spcomp1 = m_refNupp[offcomp];
                      fback *= std::pow(spec[spcomp1],m_nupp[offcomp]);
                    }
                  kbkwd *= sumtbval*fback;
                }
            }
          kfwd *= sumtbval*ffwd;
          Real diffK = (kfwd - kbkwd);
          for (int offcomp = nucomp; offcomp != nucomp+cnumr; ++offcomp)
            {
              const int spcomp1 = m_refNup[offcomp];
              Real ratechange = mmass2[spcomp1]*m_nup[offcomp]
                *diffK;
              arrRCT[MD_IX(i, spcomp1)] -= ratechange;
              destrate[spcomp1] += m_nup[offcomp]*kfwd;
            }
          for (int offcomp = nucomp; offcomp != nucomp+cnump; ++offcomp)
            {
              const int spcomp1 = m_refNupp[offcomp];
              Real ratechange = mmass2[spcomp1]*m_nupp[offcomp]
                *diffK;
              arrRCT[MD_IX(i, spcomp1)] += ratechange;
              destrate[spcomp1] += m_nupp[offcomp]*kbkwd;
            }
        }  // Loop over reactions
      Real minTau = std::numeric_limits<Real>::max();
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          Real tempdt =
            (spec[sp] + zeroConc)/(destrate[sp] + zeroConc*zeroConc);
          if (tempdt < minTau && tempdt > 0.)
            {
              minTau = tempdt;
              cnlim = sp;
            }
        }
      Real chemDt = refFromBase*minTau;
      // If ARK, useChemDt == 0, otherwise useChemDt == 1
      Real invChemDt = useChemDt/(chemDtScale*chemDt);
      arrInvDt[MD_IX(i, 0)] += invChemDt;
      if (chemDt < a_minChemDt)
        {
          a_minChemDt = chemDt;
          a_minChemDtCell = MD_GETIV(i);
        }
    }  // Loop over cells
  return cnlim;
}

/*--------------------------------------------------------------------*/
//  Solves for rho*omega, reaction rate and finds the max Jacobian of
//  the reaction source terms with respect to the conservative variables
/** \param[in]  a_box   Cell box to solve the source term
 *  \param[out] a_RCTcellPntFab
 *                      FAB to solve cell-centered reaction rate source
 *  \param[out] a_invDtFab
 *                      Inverse of the sum of the time step sizes
 *  \param[in]  a_WcellPntFab
 *                      FAB containing cell-centered primitive variables
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[out] a_minChemDt
 *                      Mininum chemical time step
 *  \param[out] a_minChemDtCell
 *                      Cell with minimum chemical time step
 *  \return     cnlim   Species associated with the smallest time step
 *//*-----------------------------------------------------------------*/

int
ThermPhysics::ARSwithDiagonalFluxCorrection(
  const Box&       a_box,
  FArrayBox&       a_RCTcellPntFab,
  FArrayBox&       a_invDtFab,
  const FArrayBox& a_WcellPntFab,
  const FArrayBox& a_UcellPntFab,
  const Real       a_stageDt,
  const int        a_level,
  Real&            a_minChemDt,
  IntVect&         a_minChemDtCell) const
{
  CH_TIME("ThermPhysics::ARSwithDiagonalFluxCorrection");
  CRD::msg << CRD::fv4 << "ThermPhysics::ARSwithDiagonalFluxCorrection"
           << CRD::end;
  const int numSpecies = CRDparam::g_numSpecies;
  const int cRho  = densityIndex();
  const int cPres = pressureIndex();
  const int cTemp = temperatureIndex();
  const int cEng  = energyFluxIndex();
  const int cVelBeg = velocityInterval().begin();
  const int cPrimSpecBeg = speciesPrimInterval().begin();
  FABSTACKTEMP(WcellPntFabAdv, a_box, numPrimitive());
  // Vector of locations
  VECSTACKTEMPASSIGN(curL, int, m_interpN, 0);
  // Vector of coefficients
  VECSTACKTEMPASSIGN(Pvals, Real, m_interpN, 0.);
  // Vector of mass fractions
  VECSTACKTEMPSIZE(spec, Real, numSpecies);
  // Indices into sorted mass fractions
  VECSTACKTEMPSIZE(sortcnj, int, numSpecies);
  WcellPntFabAdv.copy(a_WcellPntFab, a_box, 0, a_box, 0, numPrimitive());

//--Save the internal energy in a velocity slot

  MD_BOXLOOP(a_box, i)
    {
      Real rho = WcellPntFabAdv[MD_IX(i, cRho)];
      // Get the bounded species mass fractions
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int cSp = cPrimSpecBeg + sp;
          spec[sp] = std::max(0., std::min(1., WcellPntFabAdv[MD_IX(i, cSp)]));
        }

      // Normalize the mass fractions by lumping error into an inert species
      normalizePrimSpecies(NormalizeTypeRedistribute,
                           false,            // already bounded
                           spec.data(),      // mass fractions
                           1,                // stride
                           sortcnj.data());  // indirection from sorting

      // Compute properties by summing species in order from low cn to high
      // using indirection array
      Real heatofform = 0.;
      for (int j = 0; j != numSpecies; ++j)
        {
          const int sp = sortcnj[j];
          heatofform += spec[sp]*m_H0[sp];
        }

      // Velocities
      RealVect vel;
      for (const int dir : EachDir)
        {
          vel[dir] = WcellPntFabAdv[MD_IX(i, cVelBeg + dir)];
        }

      // Compute energy and stash in the first velocity component
      const Real ke = 0.5*stc::dot(vel, vel);
      WcellPntFabAdv[MD_IX(i, cVelBeg)] =
        a_UcellPntFab[MD_IX(i, cEng)]/rho - (heatofform + ke);
    }

//--Iterate using restricted dt

  int iter = 0;
  Real time = 0.;
  Real dt;
  CH_assert(time < a_stageDt);
  while (time < a_stageDt)
    {
      ++iter;
      // Need per cell
      if (iter == 1 && a_box.contains(IntVect{ 284, 11 }))
        {
          CRD::msg << "Intro frac: " << WcellPntFabAdv(IntVect{ 284, 11 }, cPrimSpecBeg + 4) << CRD::end;
          CRD::msg << "      pres: " << WcellPntFabAdv(IntVect{ 284, 11 }, cPres) << CRD::end;
          CRD::msg << "      temp: " << WcellPntFabAdv(IntVect{ 284, 11 }, cTemp) << CRD::end;
        }
      addReactionSource(a_box,
                        a_RCTcellPntFab,
                        a_invDtFab,
                        WcellPntFabAdv,
                        time,
                        a_level,
                        a_minChemDt,
                        a_minChemDtCell);
      if (a_stageDt - time < a_minChemDt)
        {
          dt = a_stageDt - time;
          time = a_stageDt;
        }
      else
        {
          dt = a_minChemDt;
          time += dt;
        }
      //Update ODE: normalize, recalc temp, get p from perfect gas law
      MD_BOXLOOP(a_box, i)
        {
          Real rho = WcellPntFabAdv[MD_IX(i, cRho)];
          for (int sp = 0; sp != numSpecies; ++sp)
            {
              const int cSp = cPrimSpecBeg + sp;
              spec[sp] = WcellPntFabAdv[MD_IX(i, cSp)] +
                dt*(a_RCTcellPntFab[MD_IX(i, sp)]/rho);
            }
          normalizePrimSpecies(NormalizeTypeRedistribute,
                               true,             // bound
                               spec.data(),      // mass fractions
                               1,                // stride
                               sortcnj.data());  // indirection from sorting
          Real Rval = 0.;
          for (int j = 0; j != numSpecies; ++j)
            {
              // Unsorted ops
              const int cSp = cPrimSpecBeg + j;
              WcellPntFabAdv[MD_IX(i, cSp)] = spec[j];
              // Sorted ops
              const int sp = sortcnj[j];
              Rval += spec[sp]*m_Rn[sp];
            }

          const Real energy = WcellPntFabAdv[MD_IX(i, cVelBeg)];
          Real Tmin = m_lookupLoT + 10.;
          Real Tmax = m_lookupHiT - 10.;
          const C2PFunc& f =
            C2PFunc(*this, energy, Rval, numSpecies, spec, Pvals, curL);
          int iterBrent = 0;
          int errorBrent = 0;
          Real Tval = RootSolver::BrentER(iterBrent, errorBrent, f, Tmin, Tmax);
          if (errorBrent != 0 || Tval != Tval || Tval > Tmax || Tval < Tmin)
            {
              CRD::msg << "Brent Solver Error Value: " << errorBrent
                       << " Brent Iters: " << iterBrent
                       << " Tmin: " << Tmin << " Tmax: " << Tmax << CRD::error;
            }
          WcellPntFabAdv[MD_IX(i, cPres)] = Tval*Rval*rho;
          WcellPntFabAdv[MD_IX(i, cTemp)] = Tval;
          if (iter == 1 && MD_GETIV(i) == IntVect{ 284, 11 })
            {
              CRD::msg << "Beg iter rate: " << a_RCTcellPntFab[MD_IX(i, 4)]
                       << "  -- dt  : " << dt << CRD::end;
              CRD::msg << "  frac: " << WcellPntFabAdv[MD_IX(i, cPrimSpecBeg + 4)] << CRD::end;
              CRD::msg << "  pres: " << WcellPntFabAdv[MD_IX(i, cPres)] << CRD::end;
              CRD::msg << "  temp: " << WcellPntFabAdv[MD_IX(i, cTemp)] << CRD::end;
            }
        }
    }
  if (a_box.contains(IntVect{ 284, 11 }))
    {
      CRD::msg << "End iter rate: " << a_RCTcellPntFab(IntVect{ 284, 11 }, 4)
               << "  -- iter: " << iter << CRD::end;
      CRD::msg << "  frac: " << WcellPntFabAdv(IntVect{ 284, 11 }, cPrimSpecBeg + 4) << CRD::end;
      CRD::msg << "  pres: " << WcellPntFabAdv(IntVect{ 284, 11 }, cPres) << CRD::end;
      CRD::msg << "  temp: " << WcellPntFabAdv(IntVect{ 284, 11 }, cTemp) << CRD::end;
    }

  // if (a_box.contains(IntVect{284,11}))
  //   {
  //     CRD::msg << "num iter: " << iter << " dt: " << dt << CRD::end;
  //   }

//--Backout dU/dt

  // Compute dU/dt for species as rho*(WcellPntFabAdv - a_WcellPntFab)/a_stageDt
  const Real factor = 1./a_stageDt;
  for (int sp = 0; sp != numSpecies; ++sp)
    {
      const int cSp = cPrimSpecBeg + sp;
      MD_BOXLOOP(a_box, i)
        {
          const Real rho = WcellPntFabAdv[MD_IX(i, cRho)];
          a_RCTcellPntFab[MD_IX(i, sp)] = rho*factor*
            (WcellPntFabAdv[MD_IX(i, cSp)] -
             a_WcellPntFab[MD_IX(i, cSp)]);
        }
    }
  if (a_box.contains(IntVect{ 284, 11 }))
    {
      CRD::msg << "End rate: " << a_RCTcellPntFab(IntVect{ 284, 11 }, 4)
               << CRD::end;
    }

  return 0;
}

/*--------------------------------------------------------------------*/
//  Compute the reaction source term Jacobian, dW/dU where only the
//  derivatives of temperature and c_n are taken with respect to the
//  conservative quantities
/** \param[in]  a_box   Cell box to compute the source jacobian
 *  \param[out] a_dW_dU
 *                      Matrix containing the dW_dU matrix
 *  \param[in]  a_WcellPntFab
 *                      FAB containing the cell-centered primitive variables
 *//*-----------------------------------------------------------------*/
void
ThermPhysics::computeReactionJacobian_dW_dU(CHMatrix&                                  a_dW_dU,
                                            const Real&                                a_rhoInv,
                                            const Real&                                a_temperature,
                                            const Real&                                a_gamma,
                                            const Real&                                a_RInv,
                                            const std::vector<Real, StackAlloc<Real>>& a_hvals,
                                            const unsigned long&                       a_numSpecies,
                                            const int&                                 a_JacTindx) const
{
  CH_TIME("ThermPhysics::computeReactionJacobian_dW_dU");
  CRD::msg << CRD::fv4 << "ThermPhysics::computeReactionJacobian_dW_dU" << CRD::end;

  // Calculate the derivatives of d(W)/d(U)
  a_dW_dU = 0; // Initialize to zero
  for(int cn = 0; cn != a_numSpecies; ++cn)
    {
      int cnIndx = cn + 1;
      // Set the temperature derivative for each species
      a_dW_dU(a_JacTindx, cn) = (a_gamma-1)*a_RInv*a_rhoInv*(m_Rn[cn]*a_temperature-a_hvals[cn]);
      // Set the derivative of c_n with respect to rho*c_n
      a_dW_dU(cnIndx, cn) = a_rhoInv;
    }
}

/*--------------------------------------------------------------------*/
//  Compute the reaction source term Jacobian, dS/dW where only the
//  derivatives of the reacting source term are taken with respect to
//  temperature and mass fraction.
/** \param[in]  a_box   Cell box to compute the source jacobian
 *  \param[out] a_dS_dW
 *                      Matrix containing the Jacobian dS/dW
 *  \param[in]  a_WcellPntFab
 *                      FAB containing the cell-centered primitive variables
 *//*-----------------------------------------------------------------*/
void
ThermPhysics::computeReactionJacobian_dS_dW(CHMatrix&                                  a_dS_dW,
                                            const Real&                                a_T,
                                            const Real&                                a_rho,
                                            const std::vector<Real, StackAlloc<Real>>& a_spec,
                                            const std::vector<Real, StackAlloc<Real>>& a_mmass1,
                                            const std::vector<Real, StackAlloc<Real>>& a_mmass2,
                                            const std::vector<Real, StackAlloc<Real>>& a_Pvals,
                                            const std::vector<int, StackAlloc<int>>&   a_curL,
                                            const unsigned long&                       a_numSpecies,
                                            const unsigned long&                       a_numReactions,
                                            const int&                                 a_JacTindx) const
{
  CH_TIME("ThermPhysics::computeReactionJacobian_dS_dW");
  CRD::msg << CRD::fv4 << "ThermPhysics::computeReactionJacobian_dS_dW" << CRD::end;

  a_dS_dW = 0; // Initialize to zero
  Real Tinv = 1./a_T;
  Real RhoInv = 1./a_rho;

  // Current third-body reaction index
  unsigned long curTB = 0;
  // Current pressure-dependent reaction index
  unsigned long curPR = 0;
  // Current arbitrary order reaction index
  unsigned long curARB = 0;
  // Loop over reactions
  for (unsigned long rctcomp = 0; rctcomp != static_cast<unsigned long>(a_numReactions); ++rctcomp)
    {
      // Number of reactants for this reaction
      const int cnumr = m_numRects[rctcomp];
      // Number of products for this reaction
      const int cnump = m_numProds[rctcomp];
      // Current index for nup and nupp, m_maxSpecPerReact
      // is the stride
      const unsigned long nucomp = static_cast<unsigned long>(m_maxSpecPerReact)*rctcomp;
      // Solve for $k_{f,r}$
      Real kfwd = 0.;
      Real kbkwd = 0.;
      Real keq = 0.;
      for (int nt = 0; nt != m_interpN; ++nt)
        {
          kfwd += m_lookupKfwd[rctcomp][static_cast<unsigned int>(a_curL[nt])]*a_Pvals[nt];
          kbkwd += m_lookupKbkwd[rctcomp][static_cast<unsigned int>(a_curL[nt])]*a_Pvals[nt];
          keq += m_lookupKeq[rctcomp][static_cast<unsigned int>(a_curL[nt])]*a_Pvals[nt];
        }
      Real sumtbval = 1.;
      unsigned long tbDerivComp = 0;
      Real rhoTBDeriv = 0;
      // If a three-body reaction
      if (m_thirdBody[rctcomp] == 1)
        {
          sumtbval = 0.;
          tbDerivComp = curTB;
          for (unsigned long sp = 0; sp != a_numSpecies; ++sp)
            {
              const unsigned long tbcomp = sp*m_numTBReact+curTB;
              sumtbval += m_tbAlpha[tbcomp]*a_spec[sp];
              rhoTBDeriv += m_tbAlpha[tbcomp]*a_spec[sp]*RhoInv;
            }
          curTB++;
        }
      else if (m_thirdBody[rctcomp] == -1)
        {
          curTB++;
        }
      // If pressure dependent reaction
      if (m_numPRReact > 0)
        {
          int presDep = m_presDep[rctcomp];
          if (presDep != 0)
            {
              Real k0, kinfty;
              Real k1 = kfwd;
              Real k2 = m_PRpreAF[curPR]*std::pow(a_T, m_PRbetai[curPR])*
                std::exp(-m_PREAR[curPR]*Tinv);
              if (presDep > 0)
                {
                  k0 = k1;
                  kinfty = k2;
                }
              else
                {
                  k0 = k2;
                  kinfty = k1;
                }
              Real Pr = k0*sumtbval/kinfty;
              Real F = m_Fc;
              if (presDep == 2)
                {
                  F = TroeSolution(Pr, a_T, static_cast<int>(curPR));
                }
              else if (presDep == 3)
                {
                  F = SRISolution(Pr, a_T, static_cast<int>(curPR));
                }
              if (presDep > 0)
                {
                  kfwd = k0*(1./(Pr + 1.))*F;
                }
              else
                {
                  kfwd = kinfty*(Pr/(Pr + 1.))*F;
                }
              curPR++;
              kbkwd = 0;
              if (m_revReact[rctcomp] >= 0)
                {
                  kbkwd = kfwd/keq;
                }
            }
        }
      // ffwd is the product of [x_n]^nup
      Real ffwd = 1.;
      // fback is the pruduct of [x_n]^nupp
      Real fback = 1.;
      Real sumNup = 0;
      Real sumNupp = 0;
      if (std::abs(m_revReact[rctcomp]) == 2)
        {
          for (unsigned long curs = 0; curs != static_cast<unsigned long>(cnumr); ++curs)
            {
              const unsigned long offcomp = static_cast<unsigned long>(nucomp) + curs;
              const unsigned long spcomp1 = static_cast<unsigned long>(m_refNup[offcomp]);
              const unsigned long arbcomp = curARB*static_cast<unsigned long>(m_maxSpecPerReact) + curs;
              ffwd *= std::pow(a_spec[spcomp1],m_ford[arbcomp]);
              sumNup += m_ford[arbcomp];
            }
          // If reaction is reversible
          if (m_revReact[rctcomp] >= 0)
            {
              for (unsigned long curs = 0; curs != static_cast<unsigned int>(cnump); ++curs)
                {
                  const unsigned long offcomp = nucomp + curs;
                  const unsigned long spcomp1 = static_cast<unsigned int>(m_refNupp[offcomp]);
                  const unsigned long arbcomp = curARB*static_cast<unsigned int>(m_maxSpecPerReact) + curs;
                  fback *= std::pow(a_spec[spcomp1],m_rord[arbcomp]);
                  sumNupp += m_rord[arbcomp];
                }
            }
          ++curARB;
        }
      else
        {
          for (unsigned long offcomp = nucomp; offcomp != nucomp+static_cast<unsigned int>(cnumr); ++offcomp)
            {
              const unsigned long spcomp1 = static_cast<unsigned int>(m_refNup[offcomp]);
              ffwd *= std::pow(a_spec[spcomp1],m_nup[offcomp]);
              sumNup += m_nup[offcomp];
            }
          // If reaction is reversible
          if (m_revReact[rctcomp] >= 0)
            {
              for (unsigned long offcomp = nucomp; offcomp != nucomp+static_cast<unsigned int>(cnump); ++offcomp)
                {
                  const unsigned long spcomp1 = static_cast<unsigned int>(m_refNupp[offcomp]);
                  fback *= std::pow(a_spec[spcomp1],m_nupp[offcomp]);
                  sumNupp += m_nupp[offcomp];
                }
            }
        }

      // Now we construct the jacobian. The jacobian is a measure of the reaction rate
      // of one species with respect to another. We only need to compute the jacobian
      // when both species are involved.
      //
      // For each species in the reaction (reaction rate species, index i)
      //   Loop over every other species in the reaction (species with respect of, index n)
      //     Add in this reaction's contribution to J_{in}
      Real kfwdDeriv = Tinv*kfwd*(m_betai[rctcomp] + m_EAR[rctcomp]*Tinv);
      Real kbkwdDeriv = 0.;
      if(keq > 1.e-16 && kbkwd > 1.e-16) // If we have a reversible reaction
        {
          kbkwdDeriv = Tinv*kfwdDeriv*(1./keq)*(m_betai[rctcomp] +
                                                m_EAR[rctcomp]*Tinv);
        }

      // Loop over reactants in this reaction
      for(unsigned long ratecomp = nucomp; ratecomp != nucomp+static_cast<unsigned int>(cnumr); ++ratecomp)
        {
          const unsigned long rcomp = static_cast<unsigned int>(m_refNup[ratecomp]);
          Real rnup = m_nup[ratecomp];

          // Compute derivatives of reactant with respect to temperature
          a_dS_dW(rcomp, a_JacTindx) += a_mmass2[rcomp]*sumtbval*rnup
            *(-kfwdDeriv*ffwd + kbkwdDeriv*fback);

          // Take the derivative of the reactant species reaction rate with
          // respect to this species (a.k.a. the derivative species)
          for(unsigned long dcomp = 0; dcomp != a_numSpecies; ++dcomp)
            {
              Real dnup = 0.;
              Real dnupp = 0.;
              const int derivSpecComp = dcomp+1;

              Real tbDeriv = 0.;
              if(m_thirdBody[rctcomp] > 0)
                {
                  const unsigned long tbcomp = dcomp*m_numTBReact+tbDerivComp;
                  tbDeriv = m_tbAlpha[tbcomp]*a_rho*a_mmass1[dcomp];
                }

              // If the derivative species is a reactant in this reaction, set dnup
              for(unsigned long derivcomp = nucomp; derivcomp != nucomp+static_cast<unsigned int>(cnumr); ++derivcomp)
                {
                  if(dcomp == m_refNup[derivcomp])
                    {
                      if(std::abs(m_revReact[rctcomp]) == 2)
                        {
                          const unsigned long arbcomp = curARB*static_cast<unsigned long>(m_maxSpecPerReact)
                            +(derivcomp-nucomp);
                          dnup = m_ford[arbcomp];
                        }
                      else
                        {
                          dnup = m_nup[derivcomp];
                        }
                    }
                }

              // If the derivative species is a product of this reaction, set dnupp
              for(unsigned long derivcomp = nucomp; derivcomp != nucomp+static_cast<unsigned int>(cnump); ++derivcomp)
                {
                  if(dcomp == m_refNupp[derivcomp])
                    {
                      if(std::abs(m_revReact[rctcomp]) == 2)
                        {
                          const unsigned long arbcomp = curARB*static_cast<unsigned long>(m_maxSpecPerReact)
                            +(derivcomp-nucomp);
                          dnupp = m_rord[arbcomp];
                        }
                      else
                        {
                          dnupp = m_nupp[derivcomp];
                        }
                    }
                }

              //  Finally compute the d(rho omega_rcomp)/d(c_dcomp) with rcomp the reactant
              //  species and dcomp the derivative species we are looping over
              Real contrib = tbDeriv*(kfwd*ffwd - kbkwd*fback);
              if(a_spec[dcomp] >= 1.e-14)
                {
                  contrib += sumtbval*a_rho*a_mmass1[dcomp]/a_spec[dcomp]*(dnup*kfwd*ffwd
                                                                     - dnupp*kbkwd*fback);
                }
              contrib *= -a_mmass2[rcomp]*rnup;

              a_dS_dW(rcomp, derivSpecComp) += contrib;
            }
        }

      // Loop over products of this reaction
      for(unsigned long ratecomp = nucomp; ratecomp != nucomp+static_cast<unsigned int>(cnump); ++ratecomp)
        {
          const unsigned long rcomp = static_cast<unsigned int>(m_refNupp[ratecomp]);
          Real rnupp = m_nupp[ratecomp];
          // Derivative of product species with respect to temperature
          a_dS_dW(rcomp, a_JacTindx) += a_mmass2[rcomp]*sumtbval*rnupp
            *(kfwdDeriv*ffwd - kbkwdDeriv*fback);

          // Take the derivative of the product species reaction rate with
          // respect to this species (a.k.a. the derivative species)
          for(unsigned long dcomp = 0; dcomp != a_numSpecies; ++dcomp)
            {
              Real dnup = 0.;
              Real dnupp = 0.;
              // const unsigned long dcomp = m_refNup[derivcomp];
              const int derivSpecComp = dcomp+1;

              Real tbDeriv = 0.;
              if(m_thirdBody[rctcomp] > 0)
                {
                  const unsigned long tbcomp = dcomp*m_numTBReact+tbDerivComp;
                  tbDeriv = m_tbAlpha[tbcomp]*a_rho*a_mmass1[dcomp];
                }

              // If the derivative species is a reactant in this reaction, set dnup
              for(unsigned long derivcomp = nucomp; derivcomp != nucomp+static_cast<unsigned int>(cnumr); ++derivcomp)
                {
                  if(dcomp == m_refNup[derivcomp])
                    {
                      if(std::abs(m_revReact[rctcomp]) == 2)
                        {
                          const unsigned long arbcomp = curARB*static_cast<unsigned long>(m_maxSpecPerReact)
                            +(derivcomp-nucomp);
                          dnup = m_ford[arbcomp];
                        }
                      else
                        {
                          dnup = m_nup[derivcomp];
                        }
                    }
                }

              // If the derivative species is a product of this reaction, set dnupp
              for(unsigned long derivcomp = nucomp; derivcomp != nucomp+static_cast<unsigned int>(cnump); ++derivcomp)
                {
                  if(dcomp == m_refNup[derivcomp])
                    {
                      if(std::abs(m_revReact[rctcomp]) == 2)
                        {
                          const unsigned long arbcomp = curARB*static_cast<unsigned long>(m_maxSpecPerReact)
                            +(derivcomp-nucomp);
                          dnupp = m_rord[arbcomp];
                        }
                      else
                        {
                          dnupp = m_nupp[derivcomp];
                        }
                    }
                }

              //  Finally compute the d(rho omega_rcomp)/d(c_sp) with rcomp the product
              //  species and sp the derivative species we are looping over
              Real contrib = tbDeriv*(kfwd*ffwd - kbkwd*fback);
              if(a_spec[dcomp] >= 1.e-14)
                {
                  contrib += sumtbval*a_rho*a_mmass1[dcomp]/a_spec[dcomp]*(dnup*kfwd*ffwd
                                                                     - dnupp*kbkwd*fback);
                }
              contrib *= a_mmass2[rcomp]*rnupp;

              a_dS_dW(rcomp,derivSpecComp) += contrib;
            }
        }
    }
}


/*--------------------------------------------------------------------*/
//  Compute the reaction source term Jacobian, dS/dU=dS/dW*dW/dU,
//  except the W computed here has temperature as a variable in place
//   of pressure. That is, W=[rho,u,v,w,T,c_i]^T instead of
//   W=[rho,u,v,w,p,c_i]^T
/** \param[in]  a_box   Cell box to compute the source jacobian
 *  \param[out] a_rxnJacobianFab
 *                      FAB to contain the reaction source Jacobian
 *  \param[in]  a_WcellPntFab
 *                      FAB containing the cell-centered primitive variables
 *//*-----------------------------------------------------------------*/
void
ThermPhysics::computeReactionJacobian(const Box&       a_box,
                                      FArrayBox&       a_rxnJacobianFab,
                                      const FArrayBox& a_WcellPntFab,
                                      const Real       a_dt) const
{
  CH_TIME("ThermPhysics::computeReactionJacobian");
  CRD::msg << CRD::fv4 << "ThermPhysics::computeReactionJacobian" << CRD::end;

  //  First primitive index for species
  const unsigned long primSpecStart = static_cast<unsigned long>(speciesPrimInterval().begin());
  // Number of species
  const unsigned long numSpecies = static_cast<unsigned long>(CRDparam::g_numSpecies);
  // Number of reactions
  const int numReactions = CRDparam::g_numReactions;
  // Get density index
  const int rhoIndx = densityIndex();
  // Get temperature index
  const int tempIndx = temperatureIndex();
  // Get the pressure index
  // const int JacTIndx = 0; // Unused
  // Cutoff temperature
  const Real Tcutoff = m_Tcutoff;

  VECSTACKTEMPSIZE(mmass1, Real, numSpecies);
  VECSTACKTEMPSIZE(mmass2, Real, numSpecies);
  VECSTACKTEMPSIZE(hval, Real, numSpecies);
  VECSTACKTEMPSIZE(spec, Real, numSpecies);
  VECSTACKTEMPSIZE(massFracs, Real, numSpecies);

  VECSTACKTEMPASSIGN(Pvals, Real, m_interpN, 0.);
  VECSTACKTEMPASSIGN(curL, int, m_interpN, 0);
  for (unsigned long sp = 0; sp != numSpecies; ++sp)
    {
      mmass1[sp] = m_preAmodUnits/m_molMass[sp];
      mmass2[sp] = m_molMass[sp]/m_preAmodUnits;
    }

  // Need to calculate gamma
  FABSTACKTEMP(gammaFab, a_box, 1); // Only need one component
  calcGamma(a_box, gammaFab, a_WcellPntFab);

  CHMatrix dW_dU(numSpecies+1, numSpecies);  // ([T,cn], rho*cn)
  CHMatrix dS_dW(numSpecies, numSpecies+1);  // (ratecn, [T,cn])
  CHMatrix dS_dU(numSpecies, numSpecies);    // (ratecn, rho*cn)

  MD_ARRAY_RESTRICT(arrGamma, gammaFab);
  MD_ARRAY_RESTRICT(arrJac, a_rxnJacobianFab);
  MD_ARRAY_RESTRICT(arrW, a_WcellPntFab);
  MD_BOXLOOP(a_box, i)
    {
      Real rho = arrW[MD_IX(i, rhoIndx)];
      Real T = arrW[MD_IX(i, tempIndx)];
      if (T < Tcutoff)
        {
          continue;
        }
      const Real gamma = arrGamma[MD_IX(i, 0)];
      // Precompute some cell constants
      Real RhoInv = 1./rho;
      // const Real Tinv = 1./T; // Unused
      Real Rval = 0.;
      lookupRefs(T, curL, Pvals);
      // Find molar concentration
      for (unsigned long sp = 0; sp != numSpecies; ++sp)
        {
          unsigned long specIndx = primSpecStart + sp;
          Real massFrac = std::min(1.,std::max(0.,arrW[MD_IX(i, static_cast<int>(specIndx))]));
          massFracs[sp] = massFrac;
          spec[sp] = massFrac*mmass1[sp]*rho;
          Rval += massFrac*m_Rn[sp];
        }

      Real RInv = 1./Rval;
      this->enthalpy(T, massFracs.data(), hval.data());

      //    Get the pressure index
      const int JacTIndx = 0;

      computeReactionJacobian_dW_dU(dW_dU,
                                    RhoInv,
                                    T,
                                    gamma,
                                    RInv,
                                    hval,
                                    numSpecies,
                                    JacTIndx);

      computeReactionJacobian_dS_dW(dS_dW,
                                    T,
                                    rho,
                                    spec,
                                    mmass1,
                                    mmass2,
                                    Pvals,
                                    curL,
                                    numSpecies,
                                    numReactions,
                                    JacTIndx);

      //  We have dS/dW and dW/dU so now do matrix mult to get dS/dU
      CHgemm(dS_dW, dW_dU, dS_dU);

      // Now put the matrix result into the jacobian FAB
      // The jacobian is row major
      //**FIXME use a transpose instead in the gemm routine.  Or better yet,
      //**FIXME use CHMatrix throughout.
      for(int row = 0; row != numSpecies; ++row)
      {
        for(int col = 0; col != numSpecies; ++col)
        {
          int linIndx = row*numSpecies + col;
          arrJac[MD_IX(i,linIndx)] = dS_dU(row,col);
        }
      }

    }

}

/*--------------------------------------------------------------------*/
//  Initialize the flow field
//  To use this function, you must specify all the primitive variables in a_W
//  Within a_W,the variables density, pressure, and temperature,
//  only two of the three
//  NOTE: the inputs for velocity, density, and pressure are unused, only
//  these values must be set in a_W
//  can be specified and the third variable MUST be set to -1
/** \param[out] a_U     Cell-centered conservative variables
 *  \param[in]  a_W     Cell-centered primitive variables
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_didx Current DataIndex on current disjoint box
 *  \param[in]  a_disjointBox
 *                     Current disjointBox
 *  \param[in]  a_box   Box to initialize over
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::initialize(FArrayBox&              a_U,
                         const FArrayBox&        a_W,
                         const LevelGridMetrics& a_gridMetrics,
                         const FluxBox&          a_unitNormals,
                         const DataIndex&        a_didx,
                         const Box&              a_disjointBox,
                         const Box&              a_box) const
{
  CH_TIME("ThermPhysics::initialize");
  CH_assert(a_W.box().contains(a_box));
  const int consSpecStart = speciesConsInterval().begin();
  const int primSpecStart = speciesPrimInterval().begin();
  const int uMomComp = vectorFluxInterval().begin();
  const int wVelComp = velocityInterval().begin();
  const int engIndx = energyFluxIndex();
  const int rhoIndx = densityIndex();
  const int presIndx = pressureIndex();
  const int tempIndx = temperatureIndex();
  const int numSpecies = CRDparam::g_numSpecies;
  int ier = 0;
  VECSTACKTEMPSIZE(spec, Real, numSpecies);
  MD_ARRAY_RESTRICT(arrW, a_W);
  MD_ARRAY_RESTRICT(arrU, a_U);
  MD_BOXLOOP(a_box, i)
    {
      Real Rgas = 0.;
      Real H0 = 0.;
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int specIndx = primSpecStart + sp;
          spec[sp] = arrW[MD_IX(i, specIndx)];
          Rgas += spec[sp]*m_Rn[sp];
          H0 += spec[sp]*m_H0[sp];
        }
      Real rho = arrW[MD_IX(i, rhoIndx)];
      Real T = arrW[MD_IX(i, tempIndx)];
      Real pres = arrW[MD_IX(i, presIndx)];
      if (rho < 0.)
        {
          rho = pres/(Rgas*T);
        }
      else if (T < 0.)
        {
          T = pres/(Rgas*rho);
        }
      else if (pres < 0.)
        {
          pres = rho*Rgas*T;
        }
      else
        {
          ier = 1;
        }
      arrU[MD_IX(i, rhoIndx)] = rho;
      D_TERM(
        Real uval = arrW[MD_IX(i, wVelComp)];,
        Real vval = arrW[MD_IX(i, wVelComp+1)];,
        Real wval = arrW[MD_IX(i, wVelComp+2)];);
      D_TERM(
        arrU[MD_IX(i, uMomComp)] = rho*uval;,
        arrU[MD_IX(i, uMomComp+1)] = rho*vval;,
        arrU[MD_IX(i, uMomComp+2)] = rho*wval;);
      Real ke = 0.5*(D_TERM(uval*uval, + vval*vval, + wval*wval));
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int specIndx = consSpecStart + sp;
          arrU[MD_IX(i, specIndx)] = rho*spec[sp]; 
        }
      Real sumfn = enthalpy(T, spec.data()) + H0;
      arrU[MD_IX(i, engIndx)] = rho*ke + rho*sumfn - pres;
    }
  if (ier > 0)
    {
      CRD::msg << "ThermPhysics::initialize: Rho, P, or T must be < 0 so the "
               << "problem is not overconstrained" << CRD::error;
    }
  // Initialize the turbulent variables
  if (CRDparam::g_turbModelType)
    {
      m_turbModel->turbInitialize(a_U,
                                  a_W,
                                  a_gridMetrics,
                                  a_unitNormals,
                                  a_didx,
                                  a_disjointBox,
                                  a_box);
    }
  CRD::msg.setTerminateOnError(true);
}

/*--------------------------------------------------------------------*/
//  Return the specific gas constant
/** \param[in]  a_cn    Beginning of pointer to species mass fractions
 *  \param[in]  a_stride
 *                      Stride to next species component.  Commonly
 *                      fab.compStride() for a BaseFab.
 *  \return             Specific gas constant for the mixture
 *
 *  \note
 *  <ul>
 *    <li> CRDparam::g_numSpecies is assumed
 *  </ul>
 *//*-----------------------------------------------------------------*/

Real
ThermPhysics::Rgas(const Real*          a_cn,
                   const std::ptrdiff_t a_cnStride) const
{
  //CH_TIME("ThermPhysics::Rgas");
  Real R = 0.;
  for (int idxSp = 0, idxSp_end = CRDparam::g_numSpecies; idxSp != idxSp_end;
       ++idxSp)
    {
      R += (*a_cn)*m_Rn[idxSp];
      a_cn += a_cnStride;
    }
  return R;
}

/*--------------------------------------------------------------------*/
//  Return the ratio of specific heats
/** \param[in]  a_T     Temperature of the mixture
 *  \param[in]  a_cn    Beginning of pointer to species mass fractions
 *  \param[in]  a_cnStride
 *                      Stride to next species component.  Commonly
 *                      fab.compStride() for a BaseFab.
 *  \return             Specific heat ratio for the mixture
 *
 *  \note
 *  <ul>
 *    <li> CRDparam::g_numSpecies is assumed
 *  </ul>
 *//*-----------------------------------------------------------------*/

Real
ThermPhysics::gamma(const Real           a_T,
                    const Real*          a_cn,
                    const std::ptrdiff_t a_cnStride) const
{
  //CH_TIME("ThermPhysics::gamma");
  // Create the specific heat value
  Real cp = this->cp(a_T, a_cn, a_cnStride);
  //**FIXME Cv = Cp - a_rgas is only valid for a single species.  You have to
  //**      individually sum up Cv and Cp.
  const Real Rgas = this->Rgas(a_cn, a_cnStride);
  return cp/(cp - Rgas);
}

/*--------------------------------------------------------------------*/
//  Solve for gamma using 'cp = h/T' and 'gamma = cp/(cp-Rgas)'
/** This assumes that the reacting cp is far different from the frozen
 *  cp (which is most likely true).  However, h and T are trusted.
 *  \param[in]  a_T     Temperature of the mixture
 *  \param[in]  a_Rgas  Specific gas constant for the mixture
 *  \param[in]  a_cn    Beginning of pointer to species mass fractions
 *  \param[in]  a_cnStride
 *                      Stride to next species component.  Commonly
 *                      fab.compStride() for a BaseFab.
 *  \return             Frozen specific heat ratio for the mixture
 *
 *  \note
 *  <ul>
 *    <li> CRDparam::g_numSpecies is assumed
 *  </ul>
 *//*-----------------------------------------------------------------*/

Real
ThermPhysics::frozenGamma(const Real           a_T,
                          const Real           a_Rgas,
                          const Real*          a_cn,
                          const std::ptrdiff_t a_cnStride) const
{
  CH_assert(a_cnStride == 1);  // Otherwise upgrade calcEnthalpy
  const Real h = enthalpy(a_T, a_cn, a_cnStride);
  const Real cp = h/a_T;
  return cp/(cp - a_Rgas);
}

/*--------------------------------------------------------------------*/
//  Solves for the gamma values
/** \param[in]  a_box   Box to put values for gamma in
 *  \param[out] a_gamma FAB containing gamma values
 *  \param[in]  a_W     Primitive variables
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::calcGamma(const Box&       a_box,
                        FArrayBox&       a_gamma,
                        const FArrayBox& a_W) const
{
  //CH_TIME("ThermPhysics::calcGamma");
  const int numSpecies = CRDparam::g_numSpecies;
  const int primSpecStart = speciesPrimInterval().begin();
  const int presIndx = pressureIndex();
  const int rhoIndx = densityIndex();
  VECSTACKTEMPSIZE(cn, Real, numSpecies);
  MD_ARRAY_RESTRICT(arrGamma, a_gamma);
  MD_ARRAY_RESTRICT(arrW, a_W);
  MD_BOXLOOP(a_box, i)
    {
      Real rho = arrW[MD_IX(i, rhoIndx)];
      Real pres = arrW[MD_IX(i, presIndx)];
      Real rgas = 0.;
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int specIndx = sp + primSpecStart;
          cn[sp] = arrW[MD_IX(i, specIndx)];
          rgas += cn[sp]*m_Rn[sp];
        }
      Real T = pres/(rgas*rho);
      arrGamma[MD_IX(i, 0)] = gamma(T, cn.data());
    }
}

/*--------------------------------------------------------------------*/
//  Solves for specific heat at constant pressure
/** \param[in]  a_T     Temperature of the mixture
 *  \param[in]  a_cn    Beginning of pointer to species mass fractions
 *  \param[in]  a_cnStride
 *                      Stride to next species component.  Commonly
 *                      fab.compStride() for a BaseFab.
 *  \return             Specific heat at constant pressure
 *//*-----------------------------------------------------------------*/

Real
ThermPhysics::cp(const Real           a_T,
                 const Real*          a_cn,
                 const std::ptrdiff_t a_cnStride) const
{
  //CH_TIME("ThermPhysics::cp");
  const int numSpecies = CRDparam::g_numSpecies;
  Real Tinv = 1./a_T;
  Real Tinv2 = Tinv/a_T;
  // Create the specific heat value
  Real cp = 0.;
  // 7-coefficient calculation
  if (m_thermFileFormat == 0)
    {
      if (a_T >= m_midLookupT)
        {
          for (int sp = 0; sp != numSpecies; ++sp)
            {
              cp += (*(a_cn + sp*a_cnStride))*m_Rn[sp]*
                (m_hnA1H[sp]*Tinv2+m_hnA2H[sp]*Tinv+m_hnA3H[sp]+a_T*
                 (m_hnA4H[sp]+a_T*(m_hnA5H[sp]+a_T*(m_hnA6H[sp]+
                                                    a_T*m_hnA7H[sp]))));
            }
        }
      else
        {
          for (int sp = 0; sp != numSpecies; ++sp)
            {
              cp += (*(a_cn + sp*a_cnStride))*m_Rn[sp]*
                (m_hnA1L[sp]*Tinv2+m_hnA2L[sp]*Tinv+m_hnA3L[sp]+a_T*
                 (m_hnA4L[sp]+a_T*(m_hnA5L[sp]+a_T*(m_hnA6L[sp]+
                                                    a_T*m_hnA7L[sp]))));
            }
        }
    }
  // 5-coefficient calculation
  else if (m_thermFileFormat == 1)
    {
      if (a_T >= m_midLookupT)
        {
          for (int sp = 0; sp != numSpecies; ++sp)
            {
              cp += (*(a_cn + sp*a_cnStride))*m_Rn[sp]*
                (m_hnA1H[sp] + a_T*(m_hnA2H[sp] +
                                    a_T*(m_hnA3H[sp] + a_T*
                                         (m_hnA4H[sp] + a_T*m_hnA5H[sp]))));
            }
        }
      else
        {
          for (int sp = 0; sp != numSpecies; ++sp)
            {
              cp += (*(a_cn + sp*a_cnStride))*m_Rn[sp]*
                (m_hnA1L[sp] + a_T*(m_hnA2L[sp] +
                                    a_T*(m_hnA3L[sp] + a_T*
                                         (m_hnA4L[sp] + a_T*m_hnA5L[sp]))));
            }
        }
    }
  return cp;
}

/*--------------------------------------------------------------------*/
//  Return the sensible enthalpy of the mixture
/** \param[in]  a_T     Temperature value
 *  \param[in]  a_cn    Beginning of pointer to species mass fractions
 *  \param[in]  a_cnStride
 *                      Stride to next species component.  Commonly
 *                      fab.compStride() for a BaseFab.
 *  \return             Enthalpy (units of J/kg)
 *
 *  \note
 *  <ul>
 *    <li> CRDparam::g_numSpecies is assumed
 *  </ul>
 *//*-----------------------------------------------------------------*/

Real
ThermPhysics::enthalpy(const Real           a_T,
                       const Real*          a_cn,
                       const std::ptrdiff_t a_cnStride) const
{
  VECSTACKTEMPSIZE(dummy, Real, CRDparam::g_numSpecies);
  return this->enthalpy(a_T, a_cn, dummy.data(), a_cnStride, 1);
}

/*--------------------------------------------------------------------*/
//  Solves for the sensible enthalpy with the heat of formation
/** \param[in]  a_T     Temperature value
 *  \param[in]  a_cn    Beginning of pointer to species mass fractions
 *  \param[out] a_hi    Vector of each species' enthalpy
 *  \param[in]  a_cnStride
 *                      Stride to next species component.  Commonly
 *                      fab.compStride() for a BaseFab.
 *  \param[in]  a_hiStride
 *                      Stride through heat of formation array
 *  \return             Enthalpy (units of J/kg)
 *
 *  \note
 *  <ul>
 *    <li> CRDparam::g_numSpecies is assumed
 *  </ul>
 *//*-----------------------------------------------------------------*/

Real
ThermPhysics::enthalpy(const Real           a_T,
                       const Real*          a_cn,
                       Real*                a_hi,
                       const std::ptrdiff_t a_cnStride,
                       const std::ptrdiff_t a_hiStride) const
{
  // Vector of locations
  VECSTACKTEMPASSIGN(curL, int, m_interpN, 0);
  // Vector of coefficients
  VECSTACKTEMPASSIGN(Pvals, Real, m_interpN, 0.);
  lookupRefs(a_T, curL, Pvals);
  Real sumfn = 0.;
  for (int sp = 0; sp != CRDparam::g_numSpecies; ++sp)
    {
      Real hval = 0.;
      for (int nt = 0; nt != m_interpN; ++nt)
        {
          hval += m_lookupH[sp][curL[nt]]*Pvals[nt];
        }
      (*a_hi) = hval;
      a_hi += a_hiStride;
      sumfn += (*a_cn)*hval;
      a_cn += a_cnStride;
    }
  return sumfn;
}

/*--------------------------------------------------------------------*/
//  Return the name of the primitive variable in the argument
//  component
/**
 *//*-----------------------------------------------------------------*/

const char *const
ThermPhysics::primStateName(const int a_iComp) const
{
  if (a_iComp < SpaceDim + 2)
    {
      static constexpr const char* NSname[] =
        {
          "density",
          D_DECL(
            "x-velocity",
            "y-velocity",
            "z-velocity"),
          "pressure"
        };
      return NSname[a_iComp];
    }
  if (speciesConsInterval().contains(a_iComp))
    {
      const int nameloc = a_iComp - speciesConsInterval().begin();
      return CRDparam::g_speciesNames[nameloc].c_str();
    }
  if (turbConsInterval().contains(a_iComp))
    {
      const int relativeComp = a_iComp - turbConsInterval().begin();
      return m_turbModel->turbStateName(relativeComp);
    }
  if (transportConsInterval().contains(a_iComp))
    {
      return "transport";  //**FIXME at some point when transport is used
    }
  if (a_iComp == temperatureIndex())
    {
      return "temperature";
    }
  return nullptr;
}

/*--------------------------------------------------------------------*/
//  Return the name of the conservative variable in the argument
//  component
/**
 *//*-----------------------------------------------------------------*/

const char *const
ThermPhysics::consvStateName(const int a_iComp) const
{
  if (a_iComp < SpaceDim + 2)
    {
      static constexpr const char* NSname[] =
        {
          "density",
          D_DECL(
            "x-momentum",
            "y-momentum",
            "z-momentum"),
          "energy-density"
        };
      return NSname[a_iComp];
    }
  if (speciesConsInterval().contains(a_iComp))
    {
      const int nameloc = a_iComp - speciesConsInterval().begin();
      return CRDparam::g_speciesNames[nameloc].c_str();
    }
  if (turbConsInterval().contains(a_iComp))
    {
      const int relativeComp = a_iComp - turbConsInterval().begin();
      return m_turbModel->turbStateName(relativeComp);
    }
  if (transportConsInterval().contains(a_iComp))
    {
      return "transport";  //**FIXME at some point when transport is used
    }
  if (a_iComp < numConservative() + 6)
    {
      static constexpr const char* name[] =
        {
          "temperature",
          "pressure",
          "cp",
          "enthalpy",
          "dynamic_viscosity",
          "thermal_conductivity"
        };
      return name[a_iComp - numConservative()];
    }
  return nullptr;
}

/*--------------------------------------------------------------------*/
//  Return name of the physics described by the class
/** \return             Name of physics
 *//*-----------------------------------------------------------------*/

const char *const
ThermPhysics::physicsName() const
{
  return "compressible Navier-Stokes with combustion";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the physics to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::writePhysicsInfo() const
{
  CRD::msg.setPrecFloatSN(5);
  CRD::msg << "Gas law\nthermally perfect ideal" << CRD::var;
  CRD::msg << "Specific heat ratio\nfunction of temperature" << CRD::var;
  if (m_numReactions > 0)
    {
      if (m_numTBReact > 0)
        {
          CRD::msg << "Number of third-body reactions\n" << m_numTBReact
                   << CRD::var;
        }
      if (m_numPRReact > 0)
        {
          CRD::msg << "Number of pressure dependent reactions\n" << m_numPRReact
                   << CRD::var;
        }
      if (m_numREVReact > 0)
        {
          CRD::msg << "Number of REV reactions\n" << m_numREVReact
                   << CRD::var;
        }
      if (m_numArbReact > 0)
        {
          CRD::msg << "Number of arbitrary order reactions\n" << m_numArbReact
                   << CRD::var;
        }
      CRD::msg << "Cutoff temperature\n" << m_Tcutoff << CRD::var;
    }
  if (CRDparam::g_useSpeciesCorrection)
    {
      CRD::msg << "Species correction\ntrue" << CRD::var;
    }
  else
    {
      CRD::msg << "Species correction\nfalse" << CRD::var;
    }
  CRD::msg << "Riemann solver\n";
    if (m_riemannSolver == RiemannSolver::Adaptive)
    {
      CRD::msg << "adaptive" << CRD::var;
    }
  else if (m_riemannSolver == RiemannSolver::Exact)
    {
      CRD::msg << "exact" << CRD::var;
    }
  else
    {
      CRD::msg << "approximate" << CRD::var;
    }
  CRD::msg.newline();
  CRD::msg.setFloatDefault();
  CRD::msg << "Parameters for mu and kappa lookup table" << CRD::var;
  CRD::msg << "Temperature bounds\n(" << m_lookupLoT << "," << m_lookupHiT
           << ")" << CRD::var;
  CRD::msg << "Number of interpolation points\n" << m_interpN << CRD::var;
  CRD::msg << "Temperature steps in lookup table\n" << m_lookupDelT << CRD::var;
  CRD::msg.newline();
  CRD::msg << "Combustion input files" << CRD::h2;
  CRD::msg << "Thermodynamic data file\n" << m_thermoFile << CRD::var;
  CRD::msg << "Thermodynamic file format\n";
  if (m_thermFileFormat == 0)
    {
      CRD::msg << "NASA - 7 coefficients" << CRD::var;
    }
  else if (m_thermFileFormat == 1)
    {
      CRD::msg << "Chemkin - 5 coefficients" << CRD::var;
    }
  if (CRDparam::g_physicsModels & CRDparam::PhysicsViscous)
    {
      if (CRDparam::g_mu < 0.)
        {
          CRD::msg << "Transport data file\n" << m_transFile << CRD::var;
        }
      if (m_schmidtNum > 0.)
        {
          CRD::msg << "Schmidt number\n" << m_schmidtNum << CRD::var;
        }
      else
        {
          CRD::msg << "Lewis number\n" << m_lewisNum << CRD::var;
        }
    }
  if (CRDparam::g_numReactions > 0)
    {
      CRD::msg << "Reaction file\n" << m_reactionFile << CRD::var;
    }
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Expressions for VisIt
/**
 *//*-----------------------------------------------------------------*/

#ifdef CH_USE_HDF5

void
ThermPhysics::expressions(HDF5HeaderData& a_holder) const
{
  const int numSpecies = CRDparam::g_numSpecies;
  // Add the species concentrations
  {
    for (int row = 0; row != numSpecies; ++row)
      {
        std::string name = "scalar rhoinv" + CRDparam::g_speciesNames[row];
        std::string speceq = "<" + CRDparam::g_speciesNames[row] +
          ">/<density>";
        a_holder.m_string[name] = speceq;
      }
  }
  a_holder.m_string["vector velocity"] = "momentum/density";
  a_holder.m_string["scalar kinetic_energy"] = "dot(velocity,velocity)/2";
  {
    // Set the solution of pressure to temperature*(sum(rho*c_n*R_n))
    // First find the R for the mixture
    std::string finalExpErCn("<density> - (");
    std::string finalExp("(");
    for (int row = 0; row != numSpecies; ++row)
      {
        finalExp += "zonal_constant(Mesh," + std::to_string(m_Rn[row]) + ")*<";
        std::string name = CRDparam::g_speciesNames[row];
        finalExp += "rhoinv" + name;
        finalExpErCn += "<"+name;
        if (row != numSpecies - 1)
          {
            finalExp += ">+";
            finalExpErCn += ">+";
          }
      }
    finalExp += ">)";
    // a_holder.m_string["scalar R_mix"] = finalExp;
    a_holder.m_string["scalar errorhoi"] = finalExpErCn + ">)";
  }
//  a_holder.m_string["scalar pressure"] = "<temperature>*<density>*<R_mix>";
  a_holder.m_string["scalar R_mix"] = "<pressure>/(<temperature>*<density>)";
  a_holder.m_string["scalar gamma"] = "<cp>/(<cp> - <R_mix>)";
  a_holder.m_string["scalar soundspeed"] =
    "sqrt(<temperature>*<R_mix>*<gamma>)";
  a_holder.m_string["scalar log10entropy"] =
    "log10(pressure) - <gamma>*log10(<density>)";
  D_TERM(
    a_holder.m_string["scalar x-velocity"] = "<x-momentum>/density";,
    a_holder.m_string["scalar y-velocity"] = "<y-momentum>/density";,
    a_holder.m_string["scalar z-velocity"] = "<z-momentum>/density";)
    a_holder.m_string["scalar machnumber"] =
    "magnitude(<velocity>)/<soundspeed>";
  // visit only calculates cartesian vorticity
  if (SpaceDim == 2)
    {
      a_holder.m_string["scalar cart_vorticity"] = "curl(velocity)";
    }
  else if (SpaceDim == 3)
    {
      a_holder.m_string["vector cart_vorticity"] = "curl(velocity)";
      a_holder.m_string["scalar cart_vorticity_mag"] =
        "magnitude(<cart_vorticity>)";
    }
  if (CRDparam::g_plotExtraVars)
    {
      if (CRDparam::g_K >= 0.)
        {
          std::string Kval = "zonal_constant(Mesh," +
            std::to_string(CRDparam::g_K) + ")";
          a_holder.m_string["scalar thermal_conductivity"] = Kval;
          std::string Muval = "zonal_constant(Mesh," +
            std::to_string(CRDparam::g_mu) + ")";
          a_holder.m_string["scalar dynamic_viscosity"] = Muval;
        }
      a_holder.m_string["scalar Prandtl_num"] =
        "<cp>*<dynamic_viscosity>/<thermal_conductivity>";
      if (m_schmidtNum > 0.)
        {
          std::string scNum = "zonal_constant(Mesh," +
            std::to_string(m_schmidtNum) + ")";
          a_holder.m_string["scalar Schmidt_num"] = scNum;
          a_holder.m_string["scalar Lewis_num"] =
            "<Schmidt_num>/<Prandtl_num>";
          a_holder.m_string["scalar Diffusion_coef"] =
            "<dynamic_viscosity>/(<density>*<Schmidt_num>)";
        }
      else
        {
          std::string lewNum = "zonal_constant(Mesh," +
            std::to_string(m_lewisNum) + ")";
          a_holder.m_string["scalar Lewis_num"] = lewNum;
          a_holder.m_string["scalar Schmidt_num"] =
            "<Prandtl_num>*<Lewis_num>";
          a_holder.m_string["scalar Diffusion_coef"] =
            "<thermal_conductivity>/(<Lewis_num>*<density>*<cp>)";
        }
    }
}

/*--------------------------------------------------------------------*/
//  Solve for the the output LevelData. This means we must solve for
//  temperature and speed of sound
/** \param[out] a_outputLD 
 *                      The level data to be outputted
 *  \param[in]  a_U     Level data containing the conservative variables
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::outputLevelData(
  LevelData<FArrayBox>&       a_outputLD,
  const LevelData<FArrayBox>& a_U,
  const LevelData<FArrayBox>& a_WOld,
  const LevelGridMetrics&     a_levelGridMetrics) const
{
  CRD::msg << CRD::fv4 << "ThermPhysics::outputLevelData" << CRD::end;
  const Interval consInterval(0, numConservative() - 1);
  a_U.copyTo(consInterval,
             a_outputLD,
             consInterval);
  int c = numConservative();
  // const int tempOutIndx  = c++;
  // const int cpOutIndx    = c++;
  // const int hIndx        = c++;
  // const int muOutIndx    = c++;
  // const int kappaOutIndx = c++;
  const int cOut_T        = c++;
  const int cOut_pressure = c++;
  const int cOut_cp       = c++;
  const int cOut_h        = c++;
  const int cOut_mu       = c++;
  const int cOut_kappa    = c++;
  // const int rhoIndx = densityIndex();
  // const int consSpecStart = speciesConsInterval().begin();
  // const int primSpecStart = speciesPrimInterval().begin();
  // const int numSpecies = CRDparam::g_numSpecies;
  // const int engIndx = energyFluxIndex();
  // const int uIndx = velocityInterval().begin();
  // const int tempIndx = temperatureIndex();
  const int numSpecies      = CRDparam::g_numSpecies;
  // const int c_rho           = densityIndex();
  // const int c_consSpecStart = speciesConsInterval().begin();
  const int c_primSpecStart = speciesPrimInterval().begin();
  // const int c_eng           = energyFluxIndex();
  // const int c_velStart      = velocityInterval().begin();
  const int c_pressure      = pressureIndex();
  const int c_T             = temperatureIndex();

  VECSTACKTEMPSIZE(spec, Real, numSpecies);
  VECSTACKTEMPASSIGN(Pvals, Real, m_interpN, 0.);
  VECSTACKTEMPASSIGN(curL, int, m_interpN, 0);

  for (DataIterator dit = a_outputLD.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& outputFab = a_outputLD[dit];

      Box boxOut = outputFab.box();
      const int idxBlk = a_outputLD.getBoxes().blockIndex(dit);
      CH_assert(idxBlk != -1);
      const BlockDomain& domain =
        a_levelGridMetrics.getCoordSys().problemDomain(idxBlk);
      boxOut &= domain;

      const FArrayBox& Ufab = a_U[dit];
      FABSTACKTEMP(Wfab, boxOut, numPrimitive());
      // This may not be the best idea and perhaps we should factor out the
      // temperature solve instead.
      consToPrim(Wfab, Ufab, boxOut, a_WOld[dit]);
      extraPrimitiveState(Wfab, boxOut);
      outputFab.copy(Wfab, boxOut, c_T, boxOut, cOut_T, 1);
      outputFab.copy(Wfab, boxOut, c_pressure, boxOut, cOut_pressure, 1);
      // if (!CRDparam::g_plotDACFDCheck)
      //   {
      //     outputFab.copy(Wfab, boxOut, c_pressure, boxOut, cOut_pressure, 1);
      //   }
      if ((!CRDparam::g_plotDACFDCheck) || CRDparam::g_plotExtraVars)
        {
          MD_BOXLOOP(boxOut, i)  // Per cell action!
            {
              // Get species at unit stride
              for (int sp = 0; sp != numSpecies; ++sp)
                {
                  spec[sp] = Wfab[MD_IX(i, c_primSpecStart + sp)];
                }
              const Real Tval = Wfab[MD_IX(i, c_T)];
              const Real cpVal = cp(Tval, spec.data());
              if (!CRDparam::g_plotDACFDCheck)
                {
                  outputFab[MD_IX(i, cOut_cp)] = cpVal;
                }
              if (CRDparam::g_plotExtraVars)
                {
                  CH_assert(!CRDparam::g_plotDACFDCheck);
                  lookupRefs(Tval, curL, Pvals);
                  if (CRDparam::g_K < 0.)
                    {
                      VECSTACKTEMPSIZE(mui, Real, numSpecies);
                      VECSTACKTEMPSIZE(kappai, Real, numSpecies);
                      Real mwMix = 0.;
                      Real enthalpy = 0.;
                      for (int sp = 0; sp != numSpecies; ++sp)
                        {
                          mwMix += spec[sp]/m_molMass[sp];
                          Real muval = 0.;
                          Real kappaval = 0.;
                          Real hval = m_H0[sp];
                          for (int nt = 0; nt != m_interpN; ++nt)
                            {
                              muval += m_lookupMu[sp][curL[nt]]*Pvals[nt];
                              kappaval += m_lookupKappa[sp][curL[nt]]*Pvals[nt];
                              hval += m_lookupH[sp][curL[nt]]*Pvals[nt];
                            }
                          mui[sp] = muval;
                          kappai[sp] = kappaval;
                          enthalpy += spec[sp]*hval;
                        }
                      mwMix = 1./mwMix;
                      Real kappaValD = 0.;
                      Real kappaValN = 0.;
                      Real muValD = 0.;
                      Real muValN = 0.;
                      for (int sp = 0; sp != numSpecies; ++sp)
                        {
                          // Molar fraction
                          const Real xi = spec[sp]*mwMix/m_molMass[sp];
                          kappaValD += kappai[sp]*xi;
                          muValD += mui[sp]*xi;
                          kappaValN += xi/kappai[sp];
                          muValN += xi/mui[sp];
                        }
                      Real muValF = 0.5*(muValD + 1./muValN);
                      outputFab[MD_IX(i, cOut_h)] = enthalpy;
                      outputFab[MD_IX(i, cOut_mu)] = muValF;
                      outputFab[MD_IX(i, cOut_kappa)] =
                        0.5*(kappaValD + 1./kappaValN);
                    }
                  else
                    {
                      Real enthalpy = 0.;
                      for (int sp = 0; sp != numSpecies; ++sp)
                        {
                          Real hval = m_H0[sp];
                          for (int nt = 0; nt != m_interpN; ++nt)
                            {
                              hval += m_lookupH[sp][curL[nt]]*Pvals[nt];
                            }
                          enthalpy += spec[sp]*hval;
                        }
                      outputFab[MD_IX(i, cOut_h)] = enthalpy;
                    }
                }
            }  // Loop over cells
        }  // If need to visit individual cells
    }  // Loop over boxes

#if 0
  std::vector<Real> spec(numSpecies,0.);
  std::vector<Real> Pvals(m_interpN,0.);
  std::vector<int> curL(m_interpN,0);
  if (CRDparam::g_plotExtraVars)
    {
      // Vector of mu for each species
      std::vector<Real> mui(numSpecies,0.);
      // Vector of kappa for each species
      std::vector<Real> kappai(numSpecies,0.);
      for (DataIterator dit = a_outputLD.dataIterator(); dit.ok(); ++dit)
        {
          FArrayBox& outputFab = a_outputLD[dit];
          const FArrayBox& WOldFab = a_WOld[dit];
          Box box = a_outputLD.getBoxes()[dit];
          MD_ARRAY_RESTRICT(arrU, outputFab);
          MD_BOXLOOP(box, i)
            {
              const Real rho = arrU[MD_IX(i, rhoIndx)];
              Real Rgas = 0.;
              Real heatofform = 0.;
              Real mwMix = 0.;
              for (int sp = 0; sp != numSpecies; ++sp)
                {
                  int specIndx = consSpecStart + sp;
                  spec[sp] = arrU[MD_IX(i, specIndx)]/rho;
                  mwMix += spec[sp]/m_molMass[sp];
                  Rgas += spec[sp]*m_Rn[sp];
                  heatofform += spec[sp]*m_H0[sp];
                }
              mwMix = 1./mwMix;
              D_TERM(
                Real uval = arrU[MD_IX(i, uIndx)]/rho;,
                Real vval = arrU[MD_IX(i, uIndx+1)]/rho;,
                Real wval = arrU[MD_IX(i, uIndx+2)]/rho;);
              Real ke = 0.5*(D_TERM(uval*uval,+ vval*vval,+ wval*wval));
              Real energy = arrU[MD_IX(i, engIndx)]/rho - (heatofform + ke);
              Real Tmin = m_lookupLoT + 10.;
              Real Tmax = m_lookupHiT - 10.;
              const C2PFunc& f = C2PFunc(*this,energy,Rgas,numSpecies,spec,
                                         Pvals,curL);
              int iterBrent = 0;
              int errorBrent = 0;
              Real Tval = RootSolver::BrentER(iterBrent,errorBrent,f,Tmin,Tmax);
              if (errorBrent != 0 || Tval != Tval)
                {
                  Tval = WOldFab[MD_IX(i, tempIndx)];
#if 0
                  if(Tval > Tmax || Tval < Tmin)
                    {
                      CRD::msg << "In cell: (" << i0 << "," << i1 << ")" << CRD::end;
                      CRD::msg << "Energy: " << energy << " Rgas: " << Rgas
                               << CRD::end;
                      CRD::msg << "SpeciesMassFracs (species corrected): " << CRD::end;
                      for(int i = 0; i != numSpecies; ++i)
                        {
                          CRD::msg << spec[i] << ", ";
                        }
                      CRD::msg << CRD::end;
                      CRD::msg << "SpeciesGasConstants: " << CRD::end;
                      for(int i = 0; i != numSpecies; ++i)
                        {
                          CRD::msg << m_Rn[i] << ", ";
                        }
                      CRD::msg << CRD::end;
                      CRD::msg << "SpeciesNames: " << CRD::end;
                      for(int i = 0; i != numSpecies; ++i)
                        {
                          CRD::msg << "'" << CRDparam::g_speciesNames[i] << "', ";
                        }
                      CRD::msg << CRD::end;
                      CRD::msg << "Brent Solver Error Value: " << errorBrent
                               << " Brent Iters: " << iterBrent
                               << " Tmin: " << Tmin << " Tmax: " << Tmax << CRD::end;
                      CRD::msg << "Bad temperature ThermPhysics::outputLevelData 1 temp: " << Tval << CRD::error;
                    }
#endif
                }
              arrU[MD_IX(i, tempOutIndx)] = Tval;
              const Real cp = calcCp(Tval, spec);
              arrU[MD_IX(i, cpOutIndx)] = cp;
              Real sumfn = 0.;
              if (CRDparam::g_K < 0.)
                {
                  lookupRefs(Tval, curL, Pvals);
                  for (int sp = 0; sp != numSpecies; ++sp)
                    {
                      Real muval = 0.;
                      Real kappaval = 0.;
                      Real hval = m_H0[sp];
                      for (int nt = 0; nt != m_interpN; ++nt)
                        {
                          muval += m_lookupMu[sp][curL[nt]]*Pvals[nt];
                          kappaval += m_lookupKappa[sp][curL[nt]]*Pvals[nt];
                          hval += m_lookupH[sp][curL[nt]]*Pvals[nt];
                        }
                      mui[sp] = muval;
                      kappai[sp] = kappaval;
                      sumfn += spec[sp]*hval;
                    }
                  Real kappaValD = 0.;
                  Real kappaValN = 0.;
                  Real muValD = 0.;
                  Real muValN = 0.;
                  for (int sp = 0; sp != numSpecies; ++sp)
                    {
                      // Molar fraction
                      const Real xi = spec[sp]*mwMix/m_molMass[sp];
                      kappaValD += kappai[sp]*xi;
                      muValD += mui[sp]*xi;
                      kappaValN += xi/kappai[sp];
                      muValN += xi/mui[sp];
                    }
                  Real muValF = 0.5*(muValD + 1./muValN);
                  arrU[MD_IX(i, muOutIndx)] = muValF;
                  arrU[MD_IX(i, kappaOutIndx)] = 0.5*(kappaValD + 1./kappaValN);
                }
              else
                {
                  sumfn = calcEnthalpy(Tval, spec) + heatofform;
                }
              arrU[MD_IX(i, hIndx)] = sumfn;
            }
        }
    }
  else if (CRDparam::g_plotDACFDCheck)
    {
      for (DataIterator dit = a_outputLD.dataIterator(); dit.ok(); ++dit)
        {
          FArrayBox& outputFab = a_outputLD[dit];
          Box box = a_outputLD.getBoxes()[dit];
          // Make a primitive fab that has the rho, c_n, and T
          MD_ARRAY_RESTRICT(arrU, outputFab);
          MD_BOXLOOP(box, i)
            {
              const Real rho = arrU[MD_IX(i, rhoIndx)];
              Real Rgas = 0.;
              Real heatofform = 0.;
              for (int sp = 0; sp != numSpecies; ++sp)
                {
                  int specIndx = consSpecStart + sp;
                  spec[sp] = arrU[MD_IX(i, specIndx)]/rho;
                  Rgas += spec[sp]*m_Rn[sp];
                  heatofform += spec[sp]*m_H0[sp];
                }
              D_TERM(
                Real uval = arrU[MD_IX(i, uIndx)]/rho;,
                Real vval = arrU[MD_IX(i, uIndx+1)]/rho;,
                Real wval = arrU[MD_IX(i, uIndx+2)]/rho;);
              Real ke = 0.5*(D_TERM(uval*uval,+ vval*vval,+ wval*wval));
              Real energy = arrU[MD_IX(i, engIndx)]/rho - (heatofform + ke);
              Real Tmin = m_lookupLoT + 10.;
              Real Tmax = m_lookupHiT - 10.;
              const C2PFunc& f = C2PFunc(*this,energy,Rgas,numSpecies,spec,
                                         Pvals,curL);
              int iterBrent = 0;
              int errorBrent = 0;
              Real Tval = RootSolver::BrentER(iterBrent,errorBrent,f,Tmin,Tmax);
              if (errorBrent != 0 || Tval != Tval)
                {
                  Tval = 298.0;
#if 0
                  // Use old value of temperature
                  if(Tval > Tmax || Tval < Tmin)
                    {
                      CRD::msg << "In cell: (" << i0 << "," << i1 << ")" << CRD::end;
                      CRD::msg << "Energy: " << energy << " Rgas: " << Rgas
                               << CRD::end;
                      CRD::msg << "SpeciesMassFracs (species corrected): " << CRD::end;
                      for(int i = 0; i != numSpecies; ++i)
                        {
                          CRD::msg << spec[i] << ", ";
                        }
                      CRD::msg << CRD::end;
                      CRD::msg << "SpeciesGasConstants: " << CRD::end;
                      for(int i = 0; i != numSpecies; ++i)
                        {
                          CRD::msg << m_Rn[i] << ", ";
                        }
                      CRD::msg << CRD::end;
                      CRD::msg << "SpeciesNames: " << CRD::end;
                      for(int i = 0; i != numSpecies; ++i)
                        {
                          CRD::msg << "'" << CRDparam::g_speciesNames[i] << "', ";
                        }
                      CRD::msg << CRD::end;
                      CRD::msg << "Brent Solver Error Value: " << errorBrent
                               << " Brent Iters: " << iterBrent
                               << " Tmin: " << Tmin << " Tmax: " << Tmax << CRD::end;
                      CRD::msg << "Bad temperature ThermPhysics::outputLevelData 2 temp: " << Tval << CRD::error;
                    }
#endif
                }
              arrU[MD_IX(i, tempOutIndx)] = Tval;
            }
        }
    
    }
  else
    {
      for (DataIterator dit = a_outputLD.dataIterator(); dit.ok(); ++dit)
        {
          FArrayBox& outputFab = a_outputLD[dit];
          Box box = a_outputLD.getBoxes()[dit];
          // Make a primitive fab that has the rho, c_n, and T
          MD_ARRAY_RESTRICT(arrU, outputFab);
          MD_BOXLOOP(box, i)
            {
              const Real rho = arrU[MD_IX(i, rhoIndx)];
              Real Rgas = 0.;
              Real heatofform = 0.;
              for (int sp = 0; sp != numSpecies; ++sp)
                {
                  int specIndx = consSpecStart + sp;
                  spec[sp] = arrU[MD_IX(i, specIndx)]/rho;
                  Rgas += spec[sp]*m_Rn[sp];
                  heatofform += spec[sp]*m_H0[sp];
                }
              D_TERM(
                Real uval = arrU[MD_IX(i, uIndx)]/rho;,
                Real vval = arrU[MD_IX(i, uIndx+1)]/rho;,
                Real wval = arrU[MD_IX(i, uIndx+2)]/rho;);
              Real ke = 0.5*(D_TERM(uval*uval,+ vval*vval,+ wval*wval));
              Real energy = arrU[MD_IX(i, engIndx)]/rho - (heatofform + ke);
              Real Tmin = m_lookupLoT + 10.;
              Real Tmax = m_lookupHiT - 10.;
              const C2PFunc& f = C2PFunc(*this,energy,Rgas,numSpecies,spec,
                                         Pvals,curL);
              int iterBrent = 0;
              int errorBrent = 0;
              Real Tval = RootSolver::BrentER(iterBrent,errorBrent,f,Tmin,Tmax);
              if (errorBrent != 0 || Tval != Tval)
                {
                  Tval = 298.0;
#if 0
                  // Use old value of temperature
                  if(Tval > Tmax || Tval < Tmin)
                    {
                      CRD::msg << "In cell: (" << i0 << "," << i1 << ")" << CRD::end;
                      CRD::msg << "Energy: " << energy << " Rgas: " << Rgas
                               << CRD::end;
                      CRD::msg << "SpeciesMassFracs (species corrected): " << CRD::end;
                      for(int i = 0; i != numSpecies; ++i)
                        {
                          CRD::msg << spec[i] << ", ";
                        }
                      CRD::msg << CRD::end;
                      CRD::msg << "SpeciesGasConstants: " << CRD::end;
                      for(int i = 0; i != numSpecies; ++i)
                        {
                          CRD::msg << m_Rn[i] << ", ";
                        }
                      CRD::msg << CRD::end;
                      CRD::msg << "SpeciesNames: " << CRD::end;
                      for(int i = 0; i != numSpecies; ++i)
                        {
                          CRD::msg << "'" << CRDparam::g_speciesNames[i] << "', ";
                        }
                      CRD::msg << CRD::end;
                      CRD::msg << "Brent Solver Error Value: " << errorBrent
                               << " Brent Iters: " << iterBrent
                               << " Tmin: " << Tmin << " Tmax: " << Tmax << CRD::end;
                      CRD::msg << "Bad temperature ThermPhysics::outputLevelData 3 temp: " << Tval << CRD::error;
                    }
#endif
                }
              arrU[MD_IX(i, tempOutIndx)] = Tval;
              const Real cp = calcCp(Tval, spec);
              arrU[MD_IX(i, cpOutIndx)] = cp;
            }
        }
    }
#endif
}
#endif

/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Solve the Riemann problem
/** This will probably be a temporary function until a better method to
 *  handle different Riemann solvers is implemented
 *  \param[in]  a_RStarSolve
 *                      Type of solver to use, either exact or approx
 *  \param[out] a_WStar Face-centered solution to Riemann problem
 *  \param[in]  a_WLeft Left primitive state, on cells to left of each
 *                      face
 *  \param[in]  a_WRight
 *                      Right primitive state, on cells to right of
 *                      each face
 *  \param[in]  a_dir   Direction of the faces
 *  \param[in]  a_box   Face-centered box on which to compute a_WStar    
 *//*-----------------------------------------------------------------*/

template <typename Func>
void
ThermPhysics::riemannSolution(const Func&      a_RStarSolve,
                              FArrayBox&       a_WStar,
                              const FArrayBox& a_WLeft,
                              const FArrayBox& a_WRight,
                              const int&       a_dir,
                              const Box&       a_box) const
{
  CH_TIME("ThermPhysics::riemannSolution");
  const Real smallp = 7.E-7;
  const Real smallr = CRDparam::g_smallr;
  const int numSpecies = CRDparam::g_numSpecies;
  const int primSpecStart = speciesPrimInterval().begin();
  const int presIndx = pressureIndex();
  const int rhoIndx = densityIndex();
  const int tempIndx = temperatureIndex();
  const int velIndx = velocityInterval().begin();
  // Left and right mass fractions
  VECSTACKTEMPSIZE(cnL, Real, numSpecies);
  VECSTACKTEMPSIZE(cnR, Real, numSpecies);
  // Indices into sorted mass fractions
  VECSTACKTEMPSIZE(sortcnjL, int, numSpecies);
  VECSTACKTEMPSIZE(sortcnjR, int, numSpecies);
  // Star state mass fractions
  VECSTACKTEMPSIZE(cnStar, Real, numSpecies);
  D_TERM(
    const int uNormIndx = velIndx + a_dir;,
    const int uTan1 = velIndx + (a_dir+1)%SpaceDim;,
    const int uTan2 = velIndx + (a_dir+2)%SpaceDim;);

  bool haveScalars = false;
  int scalarBeg = 0;
  int scalarEnd = 0;
  if (numTurbVar() + numTransportVar() > 0)
    {
      haveScalars = true;
      scalarBeg = turbPrimInterval().begin();
      scalarEnd = transportPrimInterval().end() + 1;
    }
  
  // VLA method
  MD_ARRAY_RESTRICT(arrWL, a_WLeft);
  MD_ARRAY_RESTRICT(arrWR, a_WRight);
  MD_ARRAY_RESTRICT(arrWStar, a_WStar);
  MD_BOXLOOP(a_box, i)
    {
      // Extract left and right normal velocity
      Real uL = arrWL[MD_IX(i, uNormIndx)];
      Real uR = arrWR[MD_IX(i, uNormIndx)];

      // Extract species mass fractions
      for (int sp = 0; sp != numSpecies; ++sp)
        {
          const int compSp = primSpecStart + sp;
          // cnL[sp] = arrWL[MD_IX(i, compSp)];
          // cnR[sp] = arrWR[MD_IX(i, compSp)];
          cnL[sp] = std::min(std::max(arrWL[MD_IX(i, compSp)], 0.), 1.);
          cnR[sp] = std::min(std::max(arrWR[MD_IX(i, compSp)], 0.), 1.);
        }
      if (CRDparam::g_useSpeciesCorrection)
        {
          normalizePrimSpecies(NormalizeTypeRedistribute,
                               false,             // already bounded
                               cnL.data(),        // mass fractions
                               1,                 // stride
                               sortcnjL.data());  // indirection from sorting
          normalizePrimSpecies(NormalizeTypeRedistribute,
                               false,             // already bounded
                               cnR.data(),        // mass fractions
                               1,                 // stride
                               sortcnjR.data());  // indirection from sorting
        }
      else
        {
          CH_assert(false);
          for (int j = 0; j != numSpecies; ++j)
            {
              sortcnjL[j] = j;
              sortcnjR[j] = j;
            }
        }

      // Gas constants
      Real rgasL = 0.;
      Real rgasR = 0.;
      for (int j = 0; j != numSpecies; ++j)
        {
          const int spL = sortcnjL[j];
          rgasL += cnL[spL]*m_Rn[spL];
          const int spR = sortcnjR[j];
          rgasR += cnR[spR]*m_Rn[spR];
        }

      Real rhoL = std::max(smallr, arrWL[MD_IX(i, rhoIndx)]);
      Real rhoR = std::max(smallr, arrWR[MD_IX(i, rhoIndx)]);
      Real pL = std::max(smallp, arrWL[MD_IX(i, presIndx)]);
      Real pR = std::max(smallp, arrWR[MD_IX(i, presIndx)]);
      Real TL = arrWL[MD_IX(i, tempIndx)];
      Real TR = arrWR[MD_IX(i, tempIndx)];
      // Real rhoL = std::max(smallr, pL/(rgasL*TL));
      // Real rhoR = std::max(smallr, pR/(rgasR*TR));
      // Real pL = std::max(smallp, rhoL*rgasL*TL);
      // Real pR = std::max(smallp, rhoR*rgasR*TR);
      rgasL = pL/(rhoL*TL);
      rgasR = pR/(rhoR*TR);
      Real gammaL = gamma(TL, cnL.data());
      Real gammaR = gamma(TR, cnR.data());
      Real aL = std::sqrt(gammaL*TL*rgasL);
      Real aR = std::sqrt(gammaR*TR*rgasR);
      Real CL = aL*rhoL;
      Real CR = aR*rhoR;
      // Middle pressure and velocity
      Real p3 = (CL*pR + CR*pL + CL*CR*(uL - uR))/(CL + CR);
      Real u3 = (CL*uL + CR*uR + (pL - pR))/(CL + CR);
      // Variable for if nonlinear solution was run:
      // 1 - nonlinear solution was run
      // 0 - approximate solution was used
      int runCheck = 1;
      int errorBrent = a_RStarSolve(p3, u3, CL, CR, aL, aR, pL, pR,
                                    uL, uR, gammaL, gammaR, runCheck);
      // Check if solution failed
      if (errorBrent != 0)
        {
          IntVect iervect(D_DECL(i0,i1,i2));
          CRD::msg << "ThermPhysics::riemann: error " << errorBrent
                   << " in non-linear solution at " << iervect << CRD::error;
        }

      // Final values
      Real rhoStar;
      Real uStar;
      Real pStar;
      Real RgasStar;
      // Set passive scalar values
      if (u3 >= 0.)
        {
          D_TERM(,arrWStar[MD_IX(i, uTan1)] = arrWL[MD_IX(i, uTan1)];,
                 arrWStar[MD_IX(i, uTan2)] = arrWL[MD_IX(i, uTan2)];);
          for (int spec = 0; spec != numSpecies; ++spec)
            {
              // cnStar[spec] = cnL[spec];
              const int compSp = primSpecStart + spec;
              arrWStar[MD_IX(i, compSp)] = cnL[spec];
            }
          if (haveScalars)
            {
              for (int sComp = scalarBeg; sComp != scalarEnd; ++sComp)
                {
                  arrWStar[MD_IX(i, sComp)] = arrWL[MD_IX(i, sComp)];
                }
            }
          RgasStar = rgasL;
        }
      else
        {
          D_TERM(,arrWStar[MD_IX(i, uTan1)] = arrWR[MD_IX(i, uTan1)];,
                 arrWStar[MD_IX(i, uTan2)] = arrWR[MD_IX(i, uTan2)];);
          for (int spec = 0; spec != numSpecies; ++spec)
            {
              // cnStar[spec] = cnR[spec];
              const int compSp = primSpecStart + spec;
              arrWStar[MD_IX(i, compSp)] = cnR[spec];
            }
          if (haveScalars)
            {
              for (int sComp = scalarBeg; sComp != scalarEnd; ++sComp)
                {
                  arrWStar[MD_IX(i, sComp)] = arrWR[MD_IX(i, sComp)];
                }
            }
          RgasStar = rgasR;
        }
      if (runCheck)
        {
          exactFinalValues(cnStar,rhoStar,uStar,pStar,p3,u3,rhoL,rhoR,aL,aR,
                           pL,pR,uL,uR,gammaL,gammaR,cnL,cnR);
        }
      else
        {
          approxFinalValues(cnStar,rhoStar,uStar,pStar,p3,u3,rhoL,rhoR,aL,aR,
                            pL,pR,uL,uR,gammaL,gammaR,cnL,cnR);
        }
      arrWStar[MD_IX(i, rhoIndx)] = rhoStar;
      arrWStar[MD_IX(i, uNormIndx)] = uStar;
      arrWStar[MD_IX(i, presIndx)] = pStar;
      arrWStar[MD_IX(i, tempIndx)] = pStar/(RgasStar*rhoStar);
      // if (a_dir == 0 && MD_GETIV(i) == IntVect(1, -2))
      //   {
      //     std::cerr << "State: " << rhoStar << ' ' << uStar << ' ' << pStar
      //              << ' ' << arrWStar[MD_IX(i, tempIndx)] << std::endl;
      //     std::cerr << "  P3: " << p3 << ", U3: " << u3 << ' ' << rhoL << ' ' << rhoR << ' ' << pL << ' ' << pR << ' ' << TL << ' ' << TR << ' ' << std::endl;
      //     std::cerr << "  " << rgasL << ' ' << rgasR << ' ' << gammaL << ' ' << gammaR << std::endl;
      //     CRD::msg << "State: " << rhoStar << ' ' << uStar << ' ' << pStar
      //              << ' ' << arrWStar[MD_IX(i, tempIndx)] << CRD::end;
      //   }
      // if (CRDparam::g_useSpeciesCorrection)
      //   {
      //     Real sumWithoutN2 = 0.;
      //     for (int sp = 0; sp != numSpecies; ++sp)
      //       {
      //         const int scomp = primSpecStart + sp;
      //         cnStar[sp] = std::min(std::max(cnStar[sp],0.),1.);

      //         if(sp != m_N2idx)
      //           {
      //             sumWithoutN2 += cnStar[sp];
      //             arrWStar[MD_IX(i, scomp)] = cnStar[sp];
      //           }
      //       }
      //     const int N2comp = primSpecStart + m_N2idx;
      //     arrWStar[MD_IX(i, N2comp)] = 1. - sumWithoutN2;
      //   }
      // else
      //   {
      //     for (int sp = 0; sp != numSpecies; ++sp)
      //       {
      //         const int scomp = primSpecStart + sp;
      //         arrWStar[MD_IX(i, scomp)] = cnStar[sp];
      //       }
      //   }
    }
}

/*--------------------------------------------------------------------*/
//  Return reference locations and values for lookup tables
/** \param[in]  a_T     Temperature value to lookup in table
 *  \param[out] a_locs  Vector of locations needed in lookup table
 *  \param[out] a_P     Vector of coefficients needed for lookup table
 *//*-----------------------------------------------------------------*/

void
ThermPhysics::lookupRefs(const Real&        a_T,
                         std::vector<int>&  a_locs,
                         std::vector<Real>& a_P) const
{
  if (a_T > m_lookupHiT)
    {
      CRD::msg << "Temperature " << a_T << " not in range of lookup table ("
               << m_lookupLoT << ':' << m_lookupHiT << ')' << CRD::error;
    }
  Real boundedT = std::max(a_T, m_lookupLoT); // Silently enforce lower bound
  int loc1 = int((boundedT - m_lookupLoT)/m_lookupDelT);
  Real delT = m_lookupLoT + loc1*m_lookupDelT - boundedT;
  for (int nt = 0; nt != m_interpN; ++nt)
    {
      a_locs[nt] = loc1 + m_diffCL[nt];
      a_P[nt] = m_interpDenoms[nt];
      for (int k = 0; k != m_interpN; ++k)
        {
          if (nt != k)
            {
              a_P[nt] *= -(delT + m_diffCL[k]*m_lookupDelT);
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Return reference locations and values for lookup tables
/** \param[in]  a_T     Temperature value to lookup in table
 *  \param[out] a_locs  Vector of locations needed in lookup table
 *  \param[out] a_P     Vector of coefficients needed for lookup table
 *//*-----------------------------------------------------------------*/

#ifdef USE_STACK
void
ThermPhysics::lookupRefs(const Real&                          a_T,
                         std::vector<int, StackAlloc<int>>&   a_locs,
                         std::vector<Real, StackAlloc<Real>>& a_P) const
{
  if (a_T > m_lookupHiT)
    {
      CRD::msg << "Temperature " << a_T << " not in range of lookup table ("
               << m_lookupLoT << ':' << m_lookupHiT << ')' << CRD::error;
    }
  Real boundedT = std::max(a_T, m_lookupLoT); // Silently enforce lower bound
  int loc1 = int((boundedT - m_lookupLoT)/m_lookupDelT);
  Real delT = m_lookupLoT + loc1*m_lookupDelT - boundedT;
  for (int nt = 0; nt != m_interpN; ++nt)
    {
      a_locs[nt] = loc1 + m_diffCL[nt];
      a_P[nt] = m_interpDenoms[nt];
      for (int k = 0; k != m_interpN; ++k)
        {
          if (nt != k)
            {
              a_P[nt] *= -(delT + m_diffCL[k]*m_lookupDelT);
            }
        }
    }
}
#endif

/*--------------------------------------------------------------------*/
//  Solve the F-value for a Troe pressure dependent reaction
/** \param[in]  a_Pr    Pressure dependent value
 *  \param[in]  a_T     Temperature of fluid
 *  \param[in]  a_curPR Current pressure dependent reaction
 *  \return     F       F-value for pressure dependent reaction
 *//*-----------------------------------------------------------------*/

Real
ThermPhysics::TroeSolution(const Real& a_Pr,
                           const Real& a_T,
                           const int   a_curPR) const
{
  const Real c1 = -0.4;
  const Real c2 = -0.67;
  const Real n1 = 0.75;
  const Real n2 = -1.27;
  const Real d = 0.14;
  const Real TroeAlpha = m_TroeAlpha[a_curPR];
  const Real TroeTS2 = m_TroeTS2[a_curPR];
  Real Fcent = (1. - TroeAlpha)*std::exp(-a_T/m_TroeTS3[a_curPR]) +
    TroeAlpha*std::exp(-a_T/m_TroeTS[a_curPR]);
  if (TroeTS2 > 0.)
    {
      Fcent += std::exp(-TroeTS2/a_T);
    }
  Real cval = c1 + c2*Fcent;
  Real nval = n1 + n2*Fcent;
  Real logPr = std::log(a_Pr);
  Real Ftemp = 1. + std::pow((logPr + cval)/(nval - d*(logPr + cval)), 2);
  return std::exp(std::log(Fcent)/Ftemp);
}

/*--------------------------------------------------------------------*/
//  Solve the F-value for a SRI pressure dependent reaction
/** \param[in]  a_Pr    Pressure dependent value
 *  \param[in]  a_T     Temperature of fluid
 *  \param[in]  a_curPR Current pressure dependent reaction
 *  \return     F       F-value for pressure dependent reaction
 *//*-----------------------------------------------------------------*/

Real
ThermPhysics::SRISolution(const Real& a_Pr,
                          const Real& a_T,
                          const int   a_curPR) const
{
  Real a = m_SRIa[a_curPR];
  Real b = m_SRIa[a_curPR];
  Real c = m_SRIa[a_curPR];
  Real d = m_SRIa[a_curPR];
  Real ev = m_SRIe[a_curPR];
  Real logPr = std::log(a_Pr);
  Real xval = 1./(1. + logPr*logPr);
  return d*std::pow(a*std::exp(-b/a_T) + std::exp(-a_T/c), xval)*
    std::pow(a_T, ev);
}

/**
 *  \param[out] a_sortcnj
 *                      An indirection array of mass fractions sorted
 *                      from low to high
 */
void
ThermPhysics::findTempSpeciesFromFailedBrentSolve(
  Real&                                a_temperature,
  std::vector<Real, StackAlloc<Real>>& a_massFracs,
  std::vector<int, StackAlloc<int>>&   a_sortcnj,
  const Real&                          a_energy) const
{
  int numComponents = CRDparam::g_numSpecies + 1; // One extra for temperature
  int numCorrections = 10;

  // Component ordering is T, c_n (for n = 1, ..., numSpecies)
  Vector<Real> x(numComponents);
  Vector<Real> lowerBounds(numComponents);
  Vector<Real> upperBounds(numComponents);
  Vector<int> boundsFlags(numComponents);

  Real func_val = 0;
  Vector<Real> gradient_vals(numComponents);

  x[0] = a_temperature;
  lowerBounds[0] = 290.;
  upperBounds[0] = 5590.;
  boundsFlags[0] = 2;
  gradient_vals[0] = 0.;
  for (int i = 1; i != numComponents; ++i)
    {
      x[i] = a_massFracs[i-1];
      lowerBounds[i] = 0.;
      upperBounds[i] = 1.;
      boundsFlags[i] = 2;
      gradient_vals[i] = 0.;
    }

  Real factr = 450359962737.0496;
  Real pgtol = 1.e-5;

  Vector<Real> wa(2*numComponents*numCorrections + 5*numComponents
                       + 11*numCorrections*numCorrections + 8*numCorrections);
  Vector<int> iwa(3*numComponents);
  char task[60] = "START";
  int iprint = -1;
  char csave[60];
  Vector<int> lsave(4);
  Vector<int> isave(44);
  Vector<Real> dsave(29);

  int numIters = 0;
  int maxIter = 15;
  while(true)
    {
      LBFGSBINTERFACE(&numComponents,
                      &numCorrections,
                      x.stdVector().data(),
                      lowerBounds.stdVector().data(),
                      upperBounds.stdVector().data(),
                      boundsFlags.stdVector().data(),
                      &func_val,
                      gradient_vals.stdVector().data(),
                      &factr,
                      &pgtol,
                      wa.stdVector().data(),
                      iwa.stdVector().data(),
                      task,
                      &iprint,
                      csave,
                      lsave.stdVector().data(),
                      isave.stdVector().data(),
                      dsave.stdVector().data());

      if (strncmp(task,"FG",2) == 0)
        {
          // Compute the function and gradient
          func_val = LBFGSB_Function_Eval(x.stdVector(), a_energy);
          LBFGSB_Gradient_Eval(gradient_vals.stdVector(), x.stdVector(), a_energy);
        }
      else if (strncmp(task,"NEW_X", 5) == 0)
        {
          numIters += 1;
          if (numIters >= maxIter)
            {
              strcpy(task,"STOP");
            }
        }
      else
        {
            break;
        }
    }

  // Normalize the species as they are in the objective function
  for (int c = 1; c != numComponents; ++c)
    {
      a_massFracs[c-1] = x[c];
    }
  normalizePrimSpecies(NormalizeTypeRedistribute,
                       false,               // solve should apply bounds
                       a_massFracs.data(),  // mass fractions
                       1,                   // stride
                       a_sortcnj.data());   // indirection from sorting
  a_temperature = x[0];
  if (a_temperature == 0.)
    {
      CRD::msg << "ThermPhysics::findTempSpeciesFromFailedBrentSolve "
        "temperature is equal to zero." << CRD::error;
    }

  // Check the sum
  Real sumcn = 0.;
  for (int j = 0, j_end = numComponents - 1; j != j_end; ++j)
    {
      const int sp = a_sortcnj[j];
      sumcn += a_massFracs[sp];
    }
  if (sumcn == 0.)
    {
      CRD::msg << "ThermPhysics::findTempSpeciesFromFailedBrentSolve "
        "species sum is equal to zero." << CRD::error;
    }
}

Real
ThermPhysics::LBFGSB_Function_Eval(
  const std::vector<Real>& a_currentSolution,
  Real                     a_internalEnergy) const
{
  VECSTACKTEMPASSIGN(Pvals, Real, m_interpN, 0.);
  VECSTACKTEMPASSIGN(curL, int, m_interpN, 0);

  const int numSpecies = CRDparam::g_numSpecies;
  const Real temperature = a_currentSolution[0];
  // Vector of mass fractions
  VECSTACKTEMPSIZE(massFracs, Real, numSpecies);
  // Indices into sorted mass fractions
  VECSTACKTEMPSIZE(sortcnj, int, numSpecies);

  // Compute normalization first
  for (int sp = 0; sp != numSpecies; ++sp)
    {
      massFracs[sp] = a_currentSolution[sp+1];
    }
  normalizePrimSpecies(NormalizeTypeRedistribute,
                       false,             // solver should apply bounds
                       massFracs.data(),  // mass fractions
                       1,                 // stride
                       sortcnj.data());   // indirection from sorting

  // Compute properties by summing species in order from low cn to high
  // using indirection array
  Real gasR = 0.;
  for (int j = 0; j != numSpecies; ++j)
    {
      const int sp = sortcnj[j];
      gasR += massFracs[sp]*m_Rn[sp];
    }

  auto func = C2PFunc(*this, a_internalEnergy, gasR, numSpecies,
                      massFracs, Pvals, curL);

  Real result = func(temperature);

  return Abs(result);
}

void
ThermPhysics::LBFGSB_Gradient_Eval(
  std::vector<Real>&       a_gradVals,
  const std::vector<Real>& a_currentSolution,
  Real                     a_internalEnergy) const
{
  VECSTACKTEMPASSIGN(Pvals, Real, m_interpN, 0.);
  VECSTACKTEMPASSIGN(curL, int, m_interpN, 0);

  const int numSpecies = CRDparam::g_numSpecies;
  const Real temperature = a_currentSolution[0];
  // Vector of mass fractions
  VECSTACKTEMPSIZE(massFracs, Real, numSpecies);
  // Indices into sorted mass fractions
  VECSTACKTEMPSIZE(sortcnj, int, numSpecies);

  // Compute normalization first
  for (int sp = 0; sp != numSpecies; ++sp)
    {
      massFracs[sp] = a_currentSolution[sp+1];
    }
  normalizePrimSpecies(NormalizeTypeRedistribute,
                       false,             // solve should apply bounds
                       massFracs.data(),  // mass fractions
                       1,                 // stride
                       sortcnj.data());   // indirection from sorting

  // Compute properties by summing species in order from low cn to high
  // using indirection array
  Real gasR = 0.;
  for (int j = 0; j != numSpecies; ++j)
    {
      const int sp = sortcnj[j];
      gasR += massFracs[sp]*m_Rn[sp];
    }

  auto func = C2PFunc(*this, a_internalEnergy, gasR, numSpecies,
                      massFracs, Pvals, curL);

  Real funcVal = func(temperature);
  int sign;
  if (funcVal > 0)
    {
      sign = 1;
    }
  else if (funcVal < 0)
    {
      sign = -1;
    }
  else
    {
      sign = 0;
    }

  LBFGSB_Gradient_Eval_wrt_T(a_gradVals, massFracs, temperature,
                             a_internalEnergy, sign);
  LBFGSB_Gradient_Eval_wrt_cn(a_gradVals, massFracs, temperature,
                              a_internalEnergy, sign);
}

void
ThermPhysics::LBFGSB_Gradient_Eval_wrt_T(std::vector<Real>&                          a_gradVals,
                                         const std::vector<Real, StackAlloc<Real> >& a_currentMassFracs,
                                         Real                                        a_temperature,
                                         Real                                        a_internalEnergy,
                                         int                                         a_sign) const
{
  const Real hprime = this->cp(a_temperature,  a_currentMassFracs.data());
  const Real Rgas = this->Rgas(a_currentMassFracs.data());

  a_gradVals[0] = (hprime-Rgas)*a_sign;
}

void
ThermPhysics::LBFGSB_Gradient_Eval_wrt_cn(std::vector<Real>&                          a_gradVals,
                                          const std::vector<Real, StackAlloc<Real> >& a_currentMassFracs,
                                          Real                                        a_temperature,
                                          Real                                        a_internalEnergy,
                                          int                                         a_sign) const
{
  const int numSpecies = CRDparam::g_numSpecies;

  std::vector<Real> enthalpies(numSpecies);
  enthalpy(a_temperature, a_currentMassFracs.data(), enthalpies.data());

  for(int sp = 0; sp != numSpecies; ++sp)
    {
      a_gradVals[sp+1] = (enthalpies[sp]-m_Rn[sp]*a_temperature)*a_sign;
    }
}
