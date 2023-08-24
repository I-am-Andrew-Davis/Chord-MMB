#ifdef CH_LANG_CC
/*
 *      _______               __
 *     / ___/ /  ___  __  ___/ /
 *    / /__/ _ \/ _ \/ _\/ _  /
 *    \___/_//_/\___/_/  \_._/
 *    Please refer to Copyright.txt, in Chord's root directory.
 */

/******************************************************************************
 ******************************************************************************
 ***********                                                         **********
 ******                              Chord                               ******
 ****                                                                      ****
 ****    Unsteady Solutions of the Compressible Reacting Navier-Stokes     ****
 ****          Equations Based on the Chombo Parallel AMR Library          ****
 ****                                                                      ****
 ****              Copyright 2014-2022 CFD & Propulsion Group              ****
 ****                      Colorado State University                       ****
 ****                                                                      ****
 ****    This file is part of Chord.                                       ****
 ****                                                                      ****
 ****    Chord is free software: you can redistribute it and/or modify     ****
 ****    it under the terms of the GNU General Public License as           ****
 ****    published by the Free Software Foundation, either version 3 of    ****
 ****    the License, or (at your option) any later version.               ****
 ****                                                                      ****
 ****    Chord is distributed in the hope that it will be useful,          ****
 ****    but WITHOUT ANY WARRANTY; without even the implied warranty of    ****
 ****    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     ****
 ****    GNU General Public License for more details.                      ****
 ****                                                                      ****
 ****    You should have received a copy of the GNU General Public         ****
 ****    License along with Chord.  If not, see                            ****
 ****    <http://www.gnu.org/licenses/>.                                   ****
 ****                                                                      ****
 ****    NOTE: A copy of the GPL is available in the 'doc' subdirectory    ****
 ****    of Chord                                                          ****
 ****                                                                      ****
 ****    NOTE: The Chombo library is released under a BSD-style            ****
 ****    license.  At the sole discretion of the CFD & Propulsion          ****
 ****    Group, software that is part of Chord may be moved to the         ****
 ****    Chombo library and dual-licensed under that license.              ****
 ******                                                                  ******
 **********                                                          **********
 ******************************************************************************
 ******************************************************************************/
#endif


/******************************************************************************/
/**
 * \mainpage Unsteady Solutions of the Compressible Reacting Navier-Stokes
 *
 * Chord is an application based on the Chombo Parallel AMR Library.  It is used
 * to solve the compressible reacting Navier-Stokes Equations.  Important
 * features include:
 * <ul>
 *   <li> Fourth-order accuracy in space and time
 *   <li> Adaptive mesh refinement
 *   <li> Highly parallel framework
 *   <li> Solutions on curvilinear grids by mapping
 *   <li> Monotone solutions of shock, flames and other discontinuities
 * </ul>
 *
 * \file chord.cpp
 *
 * \brief Main driver for Chord application
 *
 *//*+*************************************************************************/

/*------------------------------------------------------------------------------

  DEBUGGING WITH MPI - search for "attaching gdb" in code below

*//*--------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------

  CHECKING MEMORY WITH VALGRIND - in serial and parallel

  Turn off creation of plots and checkpoints in the input file of interest.  You
  must disable the timer and memory tracker.  Note that these  do not affect the
  build name so you probably want to fully clean Chord and Chombo before
  proceeding:

    cd Chombo/lib
    make realclean
    cd Chord/exe
    make NODEP=TRUE realclean

  For more robust reporting, you may need OPT=FALSE.  For a serial test:

    make -j16 USE_MT=FALSE USE_TIMER=FALSE example
    valgrind --leak-check=full ./chord2d.Linux.64.g++.gfortran.DEBUG.OPT.ex AMRShockBox.inputs

  There should be _NO_ leaks.

  In parallel

    make -j16 USE_MT=FALSE USE_TIMER=FALSE MPI=TRUE OPT=FALSE example
    mpirun -np 2 valgrind --leak-check=full --suppressions=/usr/share/openmpi/openmpi-valgrind.supp ./chord2d.Linux.64.mpicxx.gfortran.DEBUG.MPI.ex AMRShockBox.inputs

  You will see leaks but they should only be related to mpi or dl

  Don't forget to clean again before reverting back to standard builds with
  timers and memory tracking

*//*--------------------------------------------------------------------------*/


//----- Standard Library -----//

//----- System -----//

#include <unistd.h>

//----- Chombo Library -----//

#include "ParmParse.H"
#include "AMR.H"
#include "AMRLevel.H"
#include "CONSTANTS.H"
#include "CartesianCS.H"
#include "SingleBlockCSAdaptor.H"
#include "FABView.H"
#include "UsingNamespace.H"

//----- Internal -----//

#include "CRDparam.H"
#include "CRDmsg.H"
#include "DataTemp.H"
#include "ChordInput.H"
#include "TagMethodBuffer.H"
#include "TagLevel.H"
#include "TagLevelFactory.H"
#include "AMRLevelCNSFactory.H"
#include "DCFlattening.H"
#include "CRDMeshRefine.H"

// IBC
#include "CNSIBCFactory.H"
#include "CNSIBCGeneralized.H"
#include "CNSIBCGeneralizedSingleBlock.H"
#include "CNSIBCEulerAdvectionCube.H"
// AdvectionCylinder?
#include "CNSIBCEulerShockBox.H"
#include "CNSIBCMachReflection.H"
#include "CNSIBCShockTube.H"
#include "CNSIBCTransientCouette.H"
#include "CNSIBCTransientFlatplate.H"
#include "CNSIBCTransientPoiseuille.H"
#include "CNSIBCCombustionTest.H"
#include "CNSIBCLidDrivenCavity.H"
#include "CNSIBCMixedCouette.H"
#include "CNSIBCJetInlet.H"
#include "CNSIBCSpecShock.H"
#include "CNSIBCSpecShockBox.H"
#include "CNSIBCFlame.H"
#include "CNSIBCDetonation.H"
#include "CNSIBCBurnerTest.H"
#include "CNSIBCVortex.H"
#include "CNSIBCReactionAdvection.H"
#include "CNSIBCShockBubble.H"
#include "CNSIBCRMI.H"
#include "CNSIBCTemperatureDiffusion.H"
#include "CNSIBCSpecMachReflection.H"
// Channel ?
#include "CNSIBCTurbulentCouette.H"
#include "CNSIBCWedgeDetonation.H"
#include "CNSIBCShuOsher.H"
#include "CNSIBCCylinderFlow.H"
#include "CNSIBCTaylorGreen.H"
#include "CNSIBCIsotropicTurbulence.H"
#include "CNSIBCMixingLayer.H"
#include "CNSIBCGaussPulse.H"
#include "CNSIBCShear.H"
#include "CNSIBCRiemannCube.H"
#include "CNSIBCSpatiallyEvolvingShear.H"
#include "CNSIBCBluffBodyCombustion.H"
#include "CNSIBCMMS.H"
#include "CNSIBCRecirculatingInletTFP.H"
#include "CNSIBCSphere.H"
#include "CNSIBCObliqueWave.H"
#include "CNSIBCChannelSeparation.H"
#include "CNSIBCShockShock.H"
// Physics
#include "CNSPhysics.H"
#include "ThermPhysics.H"

// Turbulence Modeling
#include "TurbModeling.H"
#include "LES.H"
#include "SA.H"
#include "Smagorinsky.H"
#include "SSV.H"

#ifdef CH_CTHR
// Threads
#include "ThreadTeamArchitectChombo.H"
#endif

/*--------------------------------------------------------------------*
 * Definition of class global parameters declared in CRDparam.H and
 * CRDmsg.H
 *--------------------------------------------------------------------*/

#ifndef CH_MAKE_LIB
CRDparam::CRDparamVar CRDparam::CRDP;
const CRDPhysics* const& CRDparam::g_CRDPhysics = CRDparam::CRDP.CRDPhysicsOp();
const CNSIBC* const& CRDparam::g_CNSIBC = CRDparam::CRDP.CNSIBCOp();
const DCFlattening* const& CRDparam::g_DCF = CRDparam::CRDP.DCFOp();
#ifdef CH_CTHR
ThreadTools::ThreadTeamArchitectChombo_t* const & CRDparam::g_threads =
  CRDparam::CRDP.threads();
#endif
CRD::Msg& CRD::msg = CRD::Msg::getInstance();

//--Statics/globals in DataTemp

// No error checking here, but the static page size is compared against that
// used in the StackMem::alloc function where error checking is performed.
std::size_t StackMem::s_pageSize =
  static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
// Be aware that this may not allocate until used.
// This is 16 MB in 2D and 256 MB in 3D
thread_local StackMem tls_stack(16777216ull + (SpaceDim==3)*251658240);
#endif

/*--------------------------------------------------------------------*
 * Prototypes
 *--------------------------------------------------------------------*/

// From memtrack.H but used here even if memory tracking is disabled
// from the Chombo namespace
#include "NamespaceHeader.H"
void dumpmemoryatexit();
#include "NamespaceFooter.H"

// Precision to print floating point numbers
const int l_prec = 4;

const char* writeBitsInInt(const unsigned i)
{
  const int numBit = 8*sizeof(unsigned);
  static char buffer[numBit+1];
  static const char charBit[2] = { '0', '1' };
  for (int n = numBit; n--;)
    {
      buffer[numBit-n-1] = charBit[((i & (1 << n)) != 0)];
    }
  buffer[numBit] = '\0';
  return buffer;
}

struct Problem
{
  /// List of IBC types for each problem (must be same order as in CRDparam.H)
  using IBCs_tuple_type = std::tuple<
    CNSIBC,                           // 0  ProblemUndefined
    CNSIBC,                           // 1  ProblemExternal
    CNSIBCGeneralized,                // 2  ProblemGeneralized
    CNSIBCGeneralizedSingleBlock,     // 3  ProblemGeneralizedSingleBlock
    CNSIBCEulerAdvectionCube,         // 4  ProblemEulerAdvectionCube
    CNSIBC,                           // 5  ProblemEulerAdvectionCylinder
    CNSIBCEulerShockBox,              // 6  ProblemEulerShockBox
    CNSIBCMachReflection,             // 7  ProblemMachReflection
    CNSIBCShockTube,                  // 8  ProblemShockTube
    CNSIBCTransientCouette,           // 9  ProblemNavierStokesTransientCouette
    CNSIBCTransientFlatplate,         // 10 ProblemNSTransientFlatplate
    CNSIBCTransientPoiseuille,        // 11 ProblemNSTransientPoiseuille
    CNSIBCCombustionTest,             // 12 ProblemCombustion
    CNSIBCLidDrivenCavity,            // 13 ProblemLidDrivenCavity
    CNSIBCMixedCouette,               // 14 ProblemMixedCouette
    CNSIBCJetInlet,                   // 15 ProblemJetInlet
    CNSIBCSpecShock,                  // 16 ProblemSpecShock
    CNSIBCSpecShockBox,               // 17 ProblemSpecShockBox
    CNSIBCFlame,                      // 18 ProblemFlame
    CNSIBCDetonation,                 // 19 ProblemDetonation
    CNSIBCBurnerTest,                 // 20 ProblemBurnerTest
    CNSIBCVortex,                     // 21 ProblemVortex
    CNSIBCReactionAdvection,          // 22 ProblemReactionAdvection
    CNSIBCShockBubble,                // 23 ProblemShockBubble
    CNSIBCRMI,                        // 24 ProblemRMI
    CNSIBCTemperatureDiffusion,       // 25 ProblemTemperatureDiffusion
    CNSIBCSpecMachReflection,         // 26 ProblemSpecMachReflection
    CNSIBC,                           // 27 ProblemChannelFlow,
    CNSIBCTurbulentCouette,           // 28 ProblemTurbulentCouette
    CNSIBCWedgeDetonation,            // 29 ProblemWedgeDetonation
    CNSIBCShuOsher,                   // 30 ProblemShuOsher
    CNSIBCCylinderFlow,               // 31 ProblemCylinderFlow
    CNSIBCTaylorGreen,                // 32 ProblemTaylorGreen
    CNSIBCIsotropicTurbulence,        // 33 ProblemIsotropicTurbulence
    CNSIBCMixingLayer,                // 34 ProblemMixingLayer
    CNSIBCGaussPulse,                 // 35 ProblemGaussPulse
    CNSIBCShear,                      // 36 ProblemShear
    CNSIBCRiemannCube,                // 37 ProblemRiemannCube
    CNSIBCSpatiallyEvolvingShear,     // 38 ProblemSpatiallyEvolvingShear
    CNSIBCBluffBodyCombustion,        // 39 ProblemBluffBodyCombustion
    CNSIBCMMS,                        // 40 ProblemMMS
    CNSIBCRecirculatingInletTFP,      // 41 ProblemRecirculatingInletTFP
    CNSIBCSphere,                     // 42 ProblemSphere
    CNSIBCObliqueWave,                // 43 ProblemObliqueWave
    CNSIBCChannelSeparation,          // 44 ProblemChannelSeparation
    CNSIBCShockShock>;                // 45 ProblemShockShock
  /// Get IBC type for a problem type
  template <CRDparam::ProblemType PT>
  struct IBC
  {
    using type = typename std::tuple_element<static_cast<std::size_t>(PT),
                                             IBCs_tuple_type>::type;
  };
  /// Switch giving compile-time setup for run-time problem type
  template <template<typename> class SetupOp>
  static auto setup()
    {
      switch (CRDparam::g_problemType)
        {
        case CRDparam::ProblemUndefined:
          CRD::msg << "Problem type is not defined" << CRD::error;
          return SetupOp<int>::nullOp();
        case CRDparam::ProblemExternal:
          CRD::msg << "ProblemExternal is not supported" << CRD::error;
          return SetupOp<int>::nullOp();
        case CRDparam::ProblemGeneralized:
          return SetupOp<
            IBC<CRDparam::ProblemGeneralized>::type>::allocate();
        case CRDparam::ProblemGeneralizedSingleBlock:
          return SetupOp<
            IBC<CRDparam::ProblemGeneralizedSingleBlock>::type>::allocate();
        case CRDparam::ProblemEulerAdvectionCube:
          return SetupOp<
            IBC<CRDparam::ProblemEulerAdvectionCube>::type>::allocate();
        case CRDparam::ProblemEulerAdvectionCylinder:
          CRD::msg << "ProblemEulerAdvectionCylinder not supported"
                   << CRD::error;
          return SetupOp<int>::nullOp();
        case CRDparam::ProblemEulerShockBox:
          return SetupOp<
            IBC<CRDparam::ProblemEulerShockBox>::type>::allocate();
        case CRDparam::ProblemMachReflection:
          return SetupOp<
            IBC<CRDparam::ProblemMachReflection>::type>::allocate();
        case CRDparam::ProblemShockTube:
          return SetupOp<
            IBC<CRDparam::ProblemShockTube>::type>::allocate();
        case CRDparam::ProblemNavierStokesTransientCouette:
          return SetupOp<IBC<
            CRDparam::ProblemNavierStokesTransientCouette>::type>::allocate();
        case CRDparam::ProblemNavierStokesTransientFlatplate:
          return SetupOp<IBC<
            CRDparam::ProblemNavierStokesTransientFlatplate>::type>::allocate();
        case CRDparam::ProblemNavierStokesTransientPoiseuille:
          return SetupOp<IBC<CRDparam::ProblemNavierStokesTransientPoiseuille
                           >::type>::allocate();
        case CRDparam::ProblemCombustion:
          return SetupOp<
            IBC<CRDparam::ProblemCombustion>::type>::allocate();
        case CRDparam::ProblemLidDrivenCavity:
          return SetupOp<
            IBC<CRDparam::ProblemLidDrivenCavity>::type>::allocate();
        case CRDparam::ProblemMixedCouette:
          return SetupOp<
            IBC<CRDparam::ProblemMixedCouette>::type>::allocate();
        case CRDparam::ProblemJetInlet:
          return SetupOp<
            IBC<CRDparam::ProblemJetInlet>::type>::allocate();
        case CRDparam::ProblemSpecShock:
          return SetupOp<
            IBC<CRDparam::ProblemSpecShock>::type>::allocate();
        case CRDparam::ProblemSpecShockBox:
          return SetupOp<
            IBC<CRDparam::ProblemSpecShockBox>::type>::allocate();
        case CRDparam::ProblemFlame:
          return SetupOp<
            IBC<CRDparam::ProblemFlame>::type>::allocate();
        case CRDparam::ProblemDetonation:
          return SetupOp<
            IBC<CRDparam::ProblemDetonation>::type>::allocate();
        case CRDparam::ProblemBurnerTest:
          return SetupOp<
            IBC<CRDparam::ProblemBurnerTest>::type>::allocate();
        case CRDparam::ProblemVortex:
          return SetupOp<
            IBC<CRDparam::ProblemVortex>::type>::allocate();
        case CRDparam::ProblemReactionAdvection:
          return SetupOp<
            IBC<CRDparam::ProblemReactionAdvection>::type>::allocate();
        case CRDparam::ProblemShockBubble:
          return SetupOp<
            IBC<CRDparam::ProblemShockBubble>::type>::allocate();
        case CRDparam::ProblemRMI:
          return SetupOp<
            IBC<CRDparam::ProblemRMI>::type>::allocate();
        case CRDparam::ProblemTemperatureDiffusion:
          return SetupOp<
            IBC<CRDparam::ProblemTemperatureDiffusion>::type>::allocate();
        case CRDparam::ProblemSpecMachReflection:
          return SetupOp<
            IBC<CRDparam::ProblemSpecMachReflection>::type>::allocate();
        // case CRDparam::ProblemChannelFlow:
        //   return SetupOp<
        //     IBC<CRDparam::ProblemChannelFlow>::type>::allocate();
        case CRDparam::ProblemTurbulentCouette:
          return SetupOp<
            IBC<CRDparam::ProblemTurbulentCouette>::type>::allocate();
        case CRDparam::ProblemWedgeDetonation:
          return SetupOp<
            IBC<CRDparam::ProblemWedgeDetonation>::type>::allocate();
        case CRDparam::ProblemShuOsher:
          return SetupOp<
            IBC<CRDparam::ProblemShuOsher>::type>::allocate();
        case CRDparam::ProblemCylinderFlow:
          return SetupOp<
            IBC<CRDparam::ProblemCylinderFlow>::type>::allocate();
        case CRDparam::ProblemTaylorGreen:
          return SetupOp<
            IBC<CRDparam::ProblemTaylorGreen>::type>::allocate();
        case CRDparam::ProblemIsotropicTurbulence:
          return SetupOp<
            IBC<CRDparam::ProblemIsotropicTurbulence>::type>::allocate();
        case CRDparam::ProblemMixingLayer:
          return SetupOp<
            IBC<CRDparam::ProblemMixingLayer>::type>::allocate();
        case CRDparam::ProblemGaussPulse:
          return SetupOp<
            IBC<CRDparam::ProblemGaussPulse>::type>::allocate();
        case CRDparam::ProblemShear:
          return SetupOp<
            IBC<CRDparam::ProblemShear>::type>::allocate();
        case CRDparam::ProblemRiemannCube:
          return SetupOp<
            IBC<CRDparam::ProblemRiemannCube>::type>::allocate();
        case CRDparam::ProblemSpatiallyEvolvingShear:
          return SetupOp<
            IBC<CRDparam::ProblemSpatiallyEvolvingShear>::type>::allocate();
        case CRDparam::ProblemBluffBodyCombustion:
          return SetupOp<
            IBC<CRDparam::ProblemBluffBodyCombustion>::type>::allocate();
        case CRDparam::ProblemMMS:
          return SetupOp<
            IBC<CRDparam::ProblemMMS>::type>::allocate();
        case CRDparam::ProblemRecirculatingInletTFP:
          return SetupOp<
            IBC<CRDparam::ProblemRecirculatingInletTFP>::type>::allocate();
        case CRDparam::ProblemSphere:
          return SetupOp<
            IBC<CRDparam::ProblemSphere>::type>::allocate();
        case CRDparam::ProblemObliqueWave:
          return SetupOp<
            IBC<CRDparam::ProblemObliqueWave>::type>::allocate();
        case CRDparam::ProblemChannelSeparation:
          return SetupOp<
            IBC<CRDparam::ProblemChannelSeparation>::type>::allocate();
        case CRDparam::ProblemShockShock:
          return SetupOp<
            IBC<CRDparam::ProblemShockShock>::type>::allocate();
        default:
          CH_assert(false);
          return SetupOp<int>::nullOp();
        }
    }
};

// Setup operator that defines the physics.  Defers to IBC::definePhysics which
// must invoke CRDparam::definePhysics().
template <typename IBC>
struct IBCAllocatePhysics
{
  using return_type = std::vector<CRDPhysics*>;
  static return_type allocate()
    {
      return IBC::allocatePhysics();
    }
  static return_type nullOp()
    {
      return std::vector<CRDPhysics*>{};
    }
};

// Setup operator that constructs the new IBC object.
template <typename IBC>
struct IBCAllocateObject
{
  using return_type = CNSIBC*;
  static return_type allocate()
    {
      return new IBC;
    }
  static return_type nullOp()
    {
      return nullptr;
    }
};


/*******************************************************************************
 *
 * Unsteady Reacting Navier Stokes
 *
 ******************************************************************************/

#ifdef CH_MAKE_LIB
/*--------------------------------------------------------------------*/
//  Interface to the Chord program
/** Normally this is the main function.  However Chord can also be
 *  built as a library with the interface named 'chord'.  If so, there
 *  are a few extra optional arguments
 *  \param[in]  a_argc  Mimics number of command line arguments
 *  \param[in]  a_argv  Mimics command line arguments
 *  \param[in]  a_ibc   A custom IBC class to use
 *  \param[in]  a_physics
 *                      A custom physics class to use
 *  \param[in]  a_getFixedHierarchy
 *                      An optional function that produces a fixed
 *                      grid hierarchy if defined.  Often used for
 *                      tests.
 *//*-----------------------------------------------------------------*/
int chord(int a_argc, const char* a_argv[],
          CNSIBCFactory *const a_ibcFactory,
          CRDPhysics *const    a_physics,
          std::function<Vector<Vector<Box>>(const AMR&)> *const
                               a_getFixedHierarchy)
{
#else
int main(int a_argc, char* a_argv[])
{
  // Not used if this 'main' (only used in tests)
  CNSIBCFactory *const a_ibcFactory = nullptr;
  CRDPhysics *const a_physics = nullptr;
  std::function<Vector<Vector<Box>>(const AMR&)> *const a_getFixedHierarchy =
    nullptr;
#endif
#ifdef CH_MPI
  MPI_Init(&a_argc, const_cast<char***>(&a_argv));

//--For attaching gdb to a running MPI process, uncomment code below

  // 1) If using openmpi, start program with "mpirun -continuous -np ..."
  // 2) Start gdb like normal with program name
  // 3) Instead of 'r', attach to the process listed below using 'attach pid'
  // 4) Set breakpoints, etc
  // 5) Move up to this file and 'set var i = 1'
  // 6) 'c'
  // if (procID() == 0)
  //   {
  //     int i = 0;
  //     char hostname[256];
  //     gethostname(hostname, sizeof(hostname));
  //     pout() << "PID " << getpid() << " on " << hostname
  //            << " ready for attach" << std::endl;
  //     while (0 == i)
  //       sleep(5);
  //   }
  // MPI_Barrier(MPI_COMM_WORLD);
#endif

//--Set stream for output

  CRD::msg.setOutputStream(pout());


/*==============================================================================
 * Input
 *============================================================================*/

  // Check for an input file
  const char* inFile = nullptr;

  if (a_argc > 1)
    {
      inFile = a_argv[1];
    }
  else
    {
      CRD::msg << CRD::fbody << CRD::verb(-1)  // Make sure this is printed
               << "Usage: <executable name> <inputfile>" << CRD::end;
      CRD::msg << "No input file specified" << CRD::error;
      return 1;
    }

  ChomboParam chomboParam;
  ProblemDomain problemDomain;
  DCFlattening* DCF = NULL;

  // The factory for building levels. This is here because coordinate
  // system is set up poorly right now.  It is defined just after reading
  // input
  AMRLevelCNSFactory amrLevelFactory;
  
  {  // Input scope

//--The primary parser object -- keep in scope until input finished

    // Parse the command line and the input file (if any)
    ParmParse pp(a_argc-2, const_cast<char**>(a_argv+2), NULL, inFile);
    {
      // This determines the amount of diagnositic output generated
      int verbosity = 0;
      pp.query("verbosity", verbosity);
      if (verbosity < 0)
        {
          CRD::msg << "Input: 'verbosity' must be >= 0!" << CRD::error;
          verbosity = 0;
        }
      // Whether detailed dt info should be printed
      bool verbose_dt = false;
      if (pp.contains("verbose_dt"))
        {
          verbose_dt = true;
        }
      CRDparam::CRDP.defineInfo(verbosity, verbose_dt);
      CRD::msg.setVerbosity(verbosity);

      // Allow modification of linesize for I/O
      int terminalLineSize = CRD::defaultLineSize;
      pp.query("terminal_line_size", terminalLineSize);
      if (terminalLineSize <= 0)
        {
          CRD::msg << "Input: 'terminal_line_size' must be > 0!" << CRD::error;
          terminalLineSize = CRD::defaultLineSize;
        }
      CRD::msg.setLineSize(terminalLineSize);
    }

//--Print header

    {
      const char *const headerNames[2][3] =
        {
          { "*** chord ***",
            "Copyright (c) 2014-2022 CFD & Propulsion Group",
            "Colorado State University"
          },
          { "--- Chombo Library ---",
            "Copyright (c) 2000-2022, The Regents of the University of "
            "California",
            "Applied Numerical Algorithms Group, LBNL"
          }
        };
      for (int i = 0; i != 2; ++i)
        {
          for (int j = 0; j != 3; ++j)
            {
              const char *const name = headerNames[i][j];
              const int len = std::strlen(name);
              if (j == 0)
                {
                  CRD::msg.newline(-1);
                }
              if (len < CRD::msg.lineSize())
                {
                  for (int n = (CRD::msg.lineSize() - len)/2; n--;)
                    {
                      CRD::msg.pout().put(' ');
                    }
                  CRD::msg.pout() << name << std::endl;
                }
              else
                {
                  CRD::linePrint(name,
                                 CRD::msg.pout(),
                                 "",
                                 0,
                                 CRD::msg.lineSize());
                }
              if (j == 0)
                {
                  CRD::msg.newline(-1);
                }
            }
        }
      CRD::msg.newline(-1);
    }

    CRD::msg << "Starting " << SpaceDim << "D chord with input file \""
             << inFile << '\"' << CRD::h1;
    CRD::msg << "Verbosity\n" << CRDparam::g_verbosity << CRD::var;
    CRD::msg << "Verbose time step data\n"
             << (CRDparam::g_verboseDt ? "true" : "false") << CRD::var;
    CRD::msg << "Terminal Line Size\n" << CRD::msg.lineSize() << CRD::var;
    CRD::msg.newline();

#ifdef CH_MPI
    CRD::msg << "Parallel mode" << CRD::h1;
    CRD::msg << "MPI: this is processor " << CHprocID() << " of " << numProc()-1
             << " (" << numProc() << " total)" << CRD::h2;
    {
      char hostname[256];
      gethostname(hostname, sizeof(hostname));
      CRD::msg << "Hostname\n" << hostname << CRD::var;
      CRD::msg << "PID\n" << getpid() << CRD::var;
    }
#else
    CRD::msg << "Serial mode" << CRD::h1;
#endif

//--Read the input

    CRD::msg.newline();
    chordInput(chomboParam);
    CRD::msg.newline();
    if (CRD::msg.numErrors() > 0)
      {
        CRD::msg << "There " << CRD::wasNumError << " reading the input file.  "
          "Aborting now." << CRD::error;
      }

//--Define the factory now that we know the number of AMR levels

    amrLevelFactory.define();

//**FIXME Defaults are fine here, but there should be methods in IBC class
//**for customizing the CS and the tag strategy.  Maybe, as a default,
//**TagLevelFactory should itself read the input to construct a tag strategy.
//**IBC can do it's own thing otherwise.
//**UPDATE IBC can create a tag strategy.  Need to do the same still for CS.

//--Coordinate system construction.  If the pointer is non-zero, it was
//--constructed while reading input.

    // FIXME: Should be unused now
    if (chomboParam.COORDSYSfactory == NULL)
      {
        CH_assert(false);
        CRD::msg << "Coordinate System" << CRD::h2;
        switch (chomboParam.COORDSYStype)
          {
          case CRDparam::CoordSysSingleBlockExternal:
            CRD::msg << "External single block coordinate system is not yet"
              "supported.  Aborting now." << CRD::error;
          case CRDparam::CoordSysSingleBlockCartesian:
            chomboParam.COORDSYSfactory =
              new SingleBlockCSAdaptorFactory(
                new CartesianCSFactory(CRDparam::g_domainOrigin,
                                       RealVect::Unit),
                CRDparam::g_domainLength);
            CRD::msg << "Default coordinate system\nCartesian" << CRD::var;
            CRD::msg.setPrecFloatSN(l_prec);
            CRD::msg << "  Origin\n(";
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if (dir != 0)
                  {
                    CRD::msg << ',';
                  }
                CRD::msg << CRDparam::g_domainOrigin[dir];
              }
            CRD::msg << ')' << CRD::var;
            CRD::msg << "  Stretch\n(";
            for (int dir = 0; dir != SpaceDim; ++dir)
              {
                if (dir != 0)
                  {
                    CRD::msg << ',';
                  }
                CRD::msg << 1.;
              }
            CRD::msg << ')' << CRD::var;
            CRD::msg.setFloatDefault();
            CRD::msg.newline();
            break;
          case CRDparam::CoordSysSingleBlockWarped:
            CRD::msg << "Warped coordinate system must be defined in input "
              "file using keyword 'coordsys.type = warped'.  Aborting now."
                     << CRD::error;
          case CRDparam::CoordSysSingleBlockTwisted:
            CRD::msg << "Twisted coordinate system must be defined in input "
              "file using keyword 'coordsys.type = twisted'.  Aborting now."
                     << CRD::error;
          case CRDparam::CoordSysSingleBlockLogStretch:
            CRD::msg << "Log coordinate system must be defined in input "
              "file using keyword 'coordsys.type = logstretch'.  Aborting now."
                     << CRD::error;
          case CRDparam::CoordSysSingleBlockSkewed:
            CRD::msg << "Skewed coordinate system must be defined in input "
              "file using keyword 'coordsys.type = skewed'.  Aborting now."
                     << CRD::error;
          case CRDparam::CoordSysSingleBlockAnnulus:
            CRD::msg << "Annulus coordinate system must be defined in input "
              "file using keyword 'coordsys.type = annulus'.  Aborting now."
                     << CRD::error;
          // case CRDparam::CoordSysSingleBlockSchwarzChristoffelRamp:
          // {
          //   Real alpha, leadLength, rampLength;
          //   rampParameters(alpha, leadLength, rampLength);
          //   chomboParam.COORDSYSfactory =
          //     new SingleBlockCSAdaptorFactory(
          //       new SchwarzChristoffelRampCSFactory(
          //         CRDparam::g_domainBaseSize[0],
          //         alpha*Pi/180.,
          //         leadLength,
          //         rampLength));
          //   CRD::msg << "Default coordinate system\nSchwarz-Christoffel ramp"
          //            << CRD::var;
          //   CRD::msg << "  Angle (deg)\n" << alpha << CRD::var;
          //   CRD::msg << "  Lead length\n" << leadLength << CRD::var;
          //   CRD::msg << "  Ramp length\n" << rampLength << CRD::var;
          //   CRD::msg.newline();
          //   break;
          // }
          // case CRDparam::CoordSysSingleBlockLogSchwarzChristoffelRamp:
          // {
          //   Real alpha, leadLength, rampLength;
          //   rampParameters(alpha, leadLength, rampLength);
          //   RealVect shift;
          //   logshiftParameters(shift);
          //   chomboParam.COORDSYSfactory =
          //     new SingleBlockCSAdaptorFactory(
          //       new LogSchwarzChristoffelRampCSFactory(
          //         CRDparam::g_domainBaseSize[0],
          //         alpha*Pi/180.,
          //         leadLength,
          //         rampLength,
          //         shift));
          //   CRD::msg << "Default coordinate system\nLog-stretch "
          //            << "Schwarz-Christoffel ramp" << CRD::var;
          //   CRD::msg << "  Angle (deg)\n" << alpha << CRD::var;
          //   CRD::msg << "  Lead length\n" << leadLength << CRD::var;
          //   CRD::msg << "  Ramp length\n" << rampLength << CRD::var;
          //   CRD::msg.newline();
          //   break;
          // }
          case CRDparam::CoordSysMultiBlockExternal:
            CRD::msg << "External multi-block coordinate system is not yet"
              "supported.  Aborting now." << CRD::error;
          case CRDparam::CoordSysMultiBlockCartesian:
          default:
            break;
          }
      }
    if  (chomboParam.COORDSYSfactory == NULL)
      {
        CRD::msg << "Coordinate system was not provided." << CRD::error;
      }

//--Transfer some vectors to required datatype

    IntVect domainLo = IntVect_zero;
    IntVect domainHi = stc::max(IntVect_zero,
                                CRDparam::g_domainBaseSize - IntVect_unit);

    bool periodic[SpaceDim];
    SpaceDimArray<bool, int, ModConvertBool<int> >::loadFromArray(
      periodic,
      &chomboParam.GRIDperiodic.front(),
      ModConvertBool<int>(0));

//-- Build the coordinate system and make it globally accessible

    // !!FIXME!! There should be some policy here as data moves into the final
    //           data structures.  E.g., Should chomboParam.COORDSYSfactory
    //           ever be accessed after placing it in amrLevelFactory?
    amrLevelFactory.setCoordinateSystemFactory(chomboParam.COORDSYSfactory);
    // !!FIXME!! This setup is bad. We need a "domain" to set up a coordinate
    //           system, but really this describes the valid range of xi.
    //           periodicity *should* not matter, and is reset later
    problemDomain.define(domainLo, domainHi, periodic);
    // This is temporarily set here because it is needed to define the CS
    amrLevelFactory.setBaseProblemDomain(problemDomain);
    // !!FIXME!! This builds the base level, but coordinate system should be
    //           independent of the grid
    CRDparam::CRDP.defineCoordSys(amrLevelFactory.new_MultiBlockCoordSys(0));
    
//--Problem specific setup

    // Note: it does not matter if the problem domain sent to multi-block
    // coordinate system is periodic or not.  Only the number of cells in each
    // direction "may" be of interest.

    // Speaking of cells, the geometry is confused for ExternalMultiBlockCS.
    // Here, report the levelDomain as a temporary fix

    if (chomboParam.COORDSYStype == CRDparam::CoordSysMultiBlockExternal ||
        chomboParam.COORDSYStype == CRDparam::CoordSysExtrudedMultiBlock)
      {
        problemDomain = CRDparam::g_coordSys->levelDomain();
        std::vector<int> numCells(SpaceDim);
        std::vector<Real> domainLength(SpaceDim);
        std::vector<Real> domainOrigin(SpaceDim, 0.);
        for (const int dir : EachDir)
          {
            numCells[dir]     = problemDomain.size(dir);
            domainLength[dir] = problemDomain.size(dir);
          }
        CRDparam::CRDP.defineGeometry(numCells, domainLength, domainOrigin);
      }

    // First define the physics
    std::vector<CRDPhysics*> physics;
    if (a_physics)
      {
        physics.assign(1, a_physics);
      }
    else
      {
        physics = Problem::setup<IBCAllocatePhysics>();
      }
    if (physics.size() > 1)
      {
        // This invalidates CRDparam::g_CRDPhysics.  You must everywhere
        // replace CRDparam::g_CRDPhysics with CRDparam::CRDPhysicsSet()
        CH_assert(false);  // So you consider the above
        CRDparam::CRDP.defineCRDPhysics(physics.size(), physics.data());
      }
    else
      {
        CH_assert(physics[0] != nullptr);
        CRDparam::CRDP.defineCRDPhysics(physics[0]);
        CRD::msg << "Physics: " << CRDparam::g_CRDPhysics->physicsName()
                 << CRD::h2;
        CRDparam::g_CRDPhysics->writePhysicsInfo();
      }
    // Turbulence modeling
    TurbModeling* turbModel = nullptr;
    if (CRDparam::g_turbModelType)
      {
        if (CRDparam::g_turbModelType & CRDparam::TurbModelLES)
          {
            if (CRDparam::g_sgsModelType & CRDparam::SGSModelSmagorinsky)
              {
                turbModel = new Smagorinsky;
              }
            else if (CRDparam::g_sgsModelType &
                     CRDparam::SGSModelStretchedVortex)
              {
                turbModel = new SSV;
              }
          }
        else if (CRDparam::g_turbModelType & CRDparam::TurbModelSA)
          {
            turbModel = new SA;
          }
      }
    else
      {
        // Default is no model
        turbModel = new TurbModeling;
      }
    CH_assert(physics.size() == 1 && physics[0] != nullptr);
    physics[0]->defineTurbModeling(turbModel);
    if (CRDparam::g_turbModelType)
      {
        CRD::msg << "Turbulence: " << turbModel->modelName() << CRD::h2;
        turbModel->writeTurbModelInfo();
      }

    // Now define the IBC class object
    CNSIBC* ibc;
    if (a_ibcFactory)
      {
        // A custom IBC file is requested (usually for testing purposes when
        // Chord is configured as a library).
        ibc = a_ibcFactory->new_CNSIBC();
      }
    else
      {
        ibc = Problem::setup<IBCAllocateObject>();
      }
    CH_assert(ibc != nullptr);
    CRDparam::CRDP.defineIBC(ibc);
    CRD::msg << "Initial and boundary conditions: "
             << CRDparam::g_CNSIBC->IBCName() << CRD::h2;
    CRDparam::g_CNSIBC->writeIBCInfo();
    CRD::msg.newline();
    // Boundary conditions
    CRD::msg << "Boundary states and types: " << CRD::h2;
    CRDparam::g_CNSIBC->writeBoundaryConditions();

    //Create the tag level
    if (pp.countOfPrefix("tag") > 0)
      // User specified input
      {
        chomboParam.TAGfactory = new TagLevelFactory(
          chomboParam.AMRtagBufferSize);
      }
    else
      // User did not specify tagging so let the IBC classes do it.
      {
        chomboParam.TAGfactory =
          ibc->setTagMethod(chomboParam.AMRtagBufferSize);
      }
    CH_assert(chomboParam.TAGfactory != nullptr);

  }  // Input scope (cannot use PP after)

    // switch (CRDparam::g_problemType)
    //   {
    //   case CRDparam::ProblemUndefined:
    //     CRD::msg << "Problem type is note defined" << CRD::error;
    //     break;
    //   case CRDparam::ProblemExternal:
    //     CRD::msg << "ProblemExternal is not supported" << CRD::error;
    //     break;
    //   case CRDparam::ProblemEulerAdvectionCube:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCEulerAdvectionCube;
    //     physics = new CNSPhysics;
    //     break;
    //   case CRDparam::ProblemEulerAdvectionCylinder:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     CRD::msg << "ProblemEulerAdvectionCylinder not supported" << CRD::error;
    //     break;
    //   case CRDparam::ProblemEulerShockBox:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCEulerShockBox;
    //     physics = new CNSPhysics;
    //     break;
    //   case CRDparam::ProblemMachReflection:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCMachReflection;
    //     physics = new CNSPhysics;
    //     break;
    //   case CRDparam::ProblemSpecMachReflection:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCSpecMachReflection;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemShockTube:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCShockTube;
    //     physics = new CNSPhysics;
    //     break;
    //   case CRDparam::ProblemNavierStokesTransientCouette:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCTransientCouette;
    //     physics = new CNSPhysics;
    //     break;
    //   case CRDparam::ProblemTurbulentCouette:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCTransientCouette;
    //     if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    //       {
    //         physics = new ThermPhysics;
    //       }
    //     else
    //       {
    //         physics = new CNSPhysics;
    //       }
    //     break;
    //   case CRDparam::ProblemMixedCouette:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCTransientCouette;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemNavierStokesTransientFlatplate:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCTransientFlatplate;
    //     physics = new CNSPhysics;
    //     break;
    //   case CRDparam::ProblemNavierStokesTransientPoiseuille:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCTransientPoiseuille;
    //     physics = new CNSPhysics;
    //     break;
    //   case CRDparam::ProblemCombustion:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCCombustionTest;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemTemperatureDiffusion:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCTemperatureDiffusion;
    //     if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    //       {
    //         physics = new ThermPhysics;
    //       }
    //     else
    //       {
    //         physics = new CNSPhysics;
    //       }
    //     break;
    //   case CRDparam::ProblemLidDrivenCavity:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCLidDrivenCavity;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemJetInlet:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCJetInlet;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemSpecShock:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCSpecShock;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemSpecShockBox:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCSpecShockBox;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemFlame:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCFlame;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemDetonation:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCDetonation;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemShockBubble:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCShockBubble;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemRMI:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCRMI;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemBurnerTest:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCBurnerTest;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemVortex:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCVortex;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemReactionAdvection:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCReactionAdvection;
    //     physics = new ThermPhysics;
    //     break;
    //   case CRDparam::ProblemWedgeDetonation:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCWedgeDetonation;
    //     if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    //       {
    //         physics = new ThermPhysics;
    //       }
    //     else
    //       {
    //         physics = new CNSPhysics;
    //       }
    //     break;
    //   case CRDparam::ProblemCylinderFlow:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCCylinderFlow;
    //     physics = new CNSPhysics;
    //     break;
    //   case CRDparam::ProblemTaylorGreen:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCTaylorGreen;
    //     if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    //       {
    //         physics = new ThermPhysics;
    //       }
    //     else
    //       {
    //         physics = new CNSPhysics;
    //       }
    //     break;
    //   case CRDparam::ProblemMixingLayer:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCMixingLayer;
    //     if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    //       {
    //         physics = new ThermPhysics;
    //       }
    //     else
    //       {
    //         physics = new CNSPhysics;
    //       }
    //     break;
    //   case CRDparam::ProblemRiemannCube:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCRiemannCube;
    //     if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    //       {
    //         physics = new ThermPhysics;
    //       }
    //     else
    //       {
    //         physics = new CNSPhysics;
    //       }
    //     break;
    //   case CRDparam::ProblemGeneralized:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCGeneralized;
    //     if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    //       {
    //         physics = new ThermPhysics;
    //       }
    //     else
    //       {
    //         physics = new CNSPhysics;
    //       }
    //     break;
    //   case CRDparam::ProblemGeneralizedSingleBlock:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCGeneralizedSingleBlock;
    //     if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    //       {
    //         physics = new ThermPhysics;
    //       }
    //     else
    //       {
    //         physics = new CNSPhysics;
    //       }
    //     break;
    //   case CRDparam::ProblemGaussPulse:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCGaussPulse;
    //     if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    //       {
    //         physics = new ThermPhysics;
    //       }
    //     else
    //       {
    //         physics = new CNSPhysics;
    //       }
    //     break;
    //   case CRDparam::ProblemShear:
    //     problemDomain.define(domainLo, domainHi, periodic);
    //     ibc = new CNSIBCShear;
    //     if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    //       {
    //         physics = new ThermPhysics;
    //       }
    //     else
    //       {
    //         physics = new CNSPhysics;
    //       }
    //     break;
    //   }

    // Additional reporting of input
  CRD::msg << CRD::fh2 << CRD::verb(2) << "Number of ghosts" << CRD::end;
  CRD::msg << CRD::fvar << CRD::verb(2) << "<U> in cells\n"
           << CRDparam::queryGhosts(CRDparam::NumGhostUcellAvg) << CRD::end;
  CRD::msg << CRD::fvar << CRD::verb(2) << "Point U in cells\n"
           << CRDparam::queryGhosts(CRDparam::NumGhostUcellPnt) << CRD::end;
  CRD::msg << CRD::fvar << CRD::verb(2) << "Point W in cells\n"
           << CRDparam::queryGhosts(CRDparam::NumGhostWcellPnt) << CRD::end;
  CRD::msg << CRD::fvar << CRD::verb(2) << "<W> in cells\n"
           << CRDparam::queryGhosts(CRDparam::NumGhostWcellAvg) << CRD::end;
  CRD::msg << CRD::fvar << CRD::verb(2) << "<Wp> on tangent faces\n"
           << CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgTan)
           << CRD::end;
  CRD::msg << CRD::fvar << CRD::verb(2) << "Point W on tangent faces\n"
           << CRDparam::queryGhosts(CRDparam::NumGhostWfacePntTan)
           << CRD::end;
  CRD::msg << CRD::fvar << CRD::verb(2) << "<W> on tangent faces\n"
           << CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTan)
           << CRD::end;
  CRD::msg << CRD::fvar << CRD::verb(2) << "<W> on normal faces\n"
           << CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgNrm)
           << CRD::end;
  CRD::msg << CRD::fvar << CRD::verb(2) << "Point inertial F on faces\n"
           << CRDparam::queryGhosts(CRDparam::NumGhostInertialFfacePnt)
           << CRD::end;
  CRD::msg << CRD::fvar << CRD::verb(2) << "Inertial <F> on faces\n"
           << CRDparam::queryGhosts(CRDparam::NumGhostInertialFfaceAvg)
           << CRD::end;
  CRD::msg.newline(2);

//--Assign cell limiting method

  DCF = new DCFlattening(chomboParam.LIMITDCFlatTol,
                         chomboParam.LIMITFCORTol);
  CRDparam::CRDP.defineDCF(DCF);

//--Build the factory for creating AMR levels

  // AMRLevelCNSFactory amrLevelFactory;
  //**FIXME: Is this set twice?
  amrLevelFactory.setBaseProblemDomain(problemDomain);
  // amrLevelFactory.setCoordinateSystemFactory(chomboParam.COORDSYSfactory);
  amrLevelFactory.setTagLevelFactory(chomboParam.TAGfactory);

#ifdef CH_CTHR
  ThreadTools::ThreadTeamArchitectChombo_t* threads = nullptr;

  if (chomboParam.THREADuseThreads)
    {
      ThreadTools::setVerbosity(ThreadTools::Verb::v4);
      threads = new ThreadTools::ThreadTeamArchitectChombo_t();

      threads->setupTopology(chomboParam.THREADnumThreads,
                             chomboParam.THREADthreadsPerTeam,
                             0,
                             0,
                             chomboParam.THREADextraStackSize_MiB,
                             chomboParam.THREADuseHyperThreading,
                             true);

      for (unsigned teamCnt = 0; teamCnt < threads->numThreadTeams(); teamCnt++)
        {
          std::shared_ptr<ThreadTools::ChomboSubTasks_t> teamSubTasks =
            std::make_shared<ThreadTools::ChomboSubTasks_t>(
              threads->numThreadsPerTeam());
          ThreadTools::ThreadTeamSubTasksArg<DataIndex, Box, Box>
            arg{teamSubTasks};
          for (unsigned threadCnt = 0; threadCnt < threads->numThreadsPerTeam();
               threadCnt++)
            {
              threads->buildNextWorker(
                ThreadTools::ChomboSubTasks_t::threadTeamFunction, &arg);
            }
          threads->addTeamSubTasks(teamSubTasks);
        }

    CRDparam::CRDP.defineThreads(threads);
  }
#endif

/*==============================================================================
 * AMR scope - generate AMR class and run
 *============================================================================*/

  {  // Begin AMR scope

    AMR amr;

    if (CRDparam::g_coordSys->isMultiBlock())
      { 
        auto MBCSfactory =
          [&factory = amrLevelFactory]
          (const int a_level)
            {
              return factory.new_MultiBlockCoordSys(a_level).operator->();
            };
        RefCountedPtr<CRDMeshRefine<decltype(MBCSfactory)>> mbmr =
          RefCountedPtr<CRDMeshRefine<decltype(MBCSfactory)>>(
            new CRDMeshRefine<decltype(MBCSfactory)>(
              CRDparam::g_coordSys->levelDomain(),
              chomboParam.AMRrefRatios,
              chomboParam.AMRfillRatio,
              chomboParam.AMRblockFactor,
              chomboParam.AMRgridBufferSize,
              chomboParam.AMRmaxGridSize,
              chomboParam.AMRbaseLevel,
              std::move(MBCSfactory)));
        amr.setMeshRefine(mbmr);
        amr.define(chomboParam.AMRmaxLevel,
                   chomboParam.AMRrefRatios,
                   ProblemDomain{},    // AMR/AMRLevel has an empty problem
                   &amrLevelFactory);  // domain
      }
    else
      {
        amr.define(chomboParam.AMRmaxLevel,
                   chomboParam.AMRrefRatios,
                   problemDomain,
                   &amrLevelFactory);
      }

/*--ALL PROBLEM DOMAINS THAT ARE USED IN CHORD MUST BE OBTAINED FROM:
 *  a) MultiBlockCoordSys.problemDomain(Box) - for single block and multiblock
 *  b) MultiBlockCoordSys.problemDomain(0)   - allowed for single block
 *
 *  The only meaningful top-level ProblemDomain is stored in AMRLevelFactory and
 *  used solely to provide grid size information when constructing the
 *  coordinate system classes.  For single-block grids, that problem domain also
 *  informs on periodicity.
 *
 *  When building DisjointBoxLayouts, use
 *    MultiBlockCoordSys.levelDomain()
 *  to get a domain spanning all multi-block domains
 */

    amr.initialTime(chomboParam.SOLinitialTime);
    if (chomboParam.SOLfixedDt > 0.)
      {
        amr.fixedDt(chomboParam.SOLfixedDt);
      }

    // Set grid generation parameters
    amr.maxGridSize(chomboParam.AMRmaxGridSize);
    amr.blockFactor(chomboParam.AMRblockFactor);
    amr.fillRatio(chomboParam.AMRfillRatio);
    amr.regridIntervals(chomboParam.AMRregridIntervals);
    amr.gridBufferSize(chomboParam.AMRgridBufferSize);

    // Set output parameters
    amr.checkpointInterval(chomboParam.FILEcheckpointInterval);
    if (chomboParam.FILEplotPeriod > 0.)
      {
        amr.plotPeriod(chomboParam.FILEplotPeriod);
      }
    else
      {
        amr.plotInterval(chomboParam.FILEplotInterval);
      }
    amr.plotPrefix(chomboParam.FILEplotPrefix);
    amr.checkpointPrefix(chomboParam.FILEcheckpointPrefix);

    // Set time step parameters
    amr.maxDtGrow(chomboParam.SOLmaxDtGrowth);
    amr.dtToleranceFactor(chomboParam.SOLdtToleranceFactor);
    amr.useSubcyclingInTime(chomboParam.AMRuseSubcycling);

    amr.verbosity(CRDparam::g_verbosity);

    // Set up input files
    if (chomboParam.FILEdoRestart)
      {
#ifdef CH_USE_HDF5
        HDF5Handle handle(chomboParam.FILErestartFile, HDF5Handle::OPEN_RDONLY);
        // Read from checkpoint file
        amr.setupForRestart(handle, chomboParam.AMRregridOnRestart);
        handle.close();
#else
        CRD::msg << "Restart only defined with hdf5 enabled" << CRD::error;
#endif
      }
    else if (a_getFixedHierarchy)
      {
        Vector<Vector<Box>> levelBoxes = (*a_getFixedHierarchy)(amr);
        amr.setupForFixedHierarchyRun(levelBoxes, 0);  // Second argument unused
      }
    else
      {
        amr.setupForNewAMRRun();
      }

    // Run the computation
    CRD::msg << "Begin computation" << CRD::h1;
    amr.run(chomboParam.SOLmaxTime, chomboParam.SOLmaxStep);
    CRD::msg.newline();

    // Output the last plot file and statistics.
    CRD::msg << "End computation" << CRD::h1;
    amr.conclude();

#ifdef CH_CTHR
    if (threads)
      {
        threads->waitToFinish();
      }
#endif
  }  // End AMR scope

#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
  return 0;
}
