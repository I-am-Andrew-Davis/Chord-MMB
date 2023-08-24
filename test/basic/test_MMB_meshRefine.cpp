#ifdef CH_LANG_CC
/*
 *      _______               __
 *     / ___/ /  ___  __  ___/ /
 *    / /__/ _ \/ _ \/ _\/ _  /
 *    \___/_//_/\___/_/  \_._/
 *    Please refer to Copyright.txt, in Chord's root directory.
 */
#endif

/*******************************************************************************
 *
 * This test verifies the functionality of mesh refine on mapped multiblock
 * geometries.
 * Primary tested components are:
 *   - Chombo:BRMeshRefine.{H|cpp}, MeshRefine.{H.cpp}
 *   - Chord:CRDMeshRefine.{H|cpp}
 * Primary tested concepts are:
 *   - Proper nesting across multiblock boundaries
 *   - Multiblock mesh refine
 *
 ******************************************************************************/

#define TESTFUNC test_MMB_meshRefine

#define _STR_(x) #x                       /* Convert TESTFUNC name to string */
#define _NAMESTR_(x) _STR_(x)
#define TESTFUNCNAME _NAMESTR_(TESTFUNC)  /* String name */

//----- Standard Library -----//

#include <iostream>
#include <cstring>
#include <limits>

//----- System -----//

#ifdef CH_MPI
#include <mpi.h>
#endif

//----- Chombo Library -----//

#include "RealVect.H"
#include "CONSTANTS.H"
#include "parstream.H"
#include "parseTestOptions.H"
#include "FABView.H"
#include "DataTemp.H"
#include "ProblemDomain.H"
#include "LoHiCenter.H"
#include "GodunovUtilitiesF_F.H"
#include "BaseFabMacros.H"
#include "DoubleCartesianCS.H"
// #include "DoubleCartesianRotateCS.H"
#include "LoadBalance.H"
#include "LevelGridMetrics.H"
#include "MeshRefine.H"

#include "UsingNamespace.H"

//----- Internal -----//

#include "CRDparam.H"
#include "CRDmsg.H"
#include "AMRLevelCNSFactory.H"
#include "AMRLevelCNS.H"
#include "CRDMeshRefine.H"


/// Prototypes:
int
TESTFUNC();

/// Global variables for handling output:
static const char *indent = "   ", *indent2 = "      ";
static bool verbose = true;

int
main(int argc, const char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, const_cast<char***>(&argv));
#endif
  parseTestOptions(argc, argv, verbose);

  if (verbose)
    {
      pout() << indent2 << "Beginning " << TESTFUNCNAME << " ..." << std::endl;
    }

//--Run the tests

  int ret = TESTFUNC();
  if (ret == 0)
    {
      pout() << indent << TESTFUNCNAME << " passed all tests" << std::endl;
    }
  else
    {
      pout() << indent << TESTFUNCNAME << " failed " << ret << " test(s)"
             << std::endl;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}


/*--------------------------------------------------------------------*
 * Definition of class global parameters declared in CRDparam.H and
 * CRDmsg.H
 *--------------------------------------------------------------------*/

CRDparam::CRDparamVar CRDparam::CRDP;
const CRDPhysics* const& CRDparam::g_CRDPhysics = CRDparam::CRDP.CRDPhysicsOp();
const CNSIBC* const& CRDparam::g_CNSIBC = CRDparam::CRDP.CNSIBCOp();
const DCFlattening* const& CRDparam::g_DCF = CRDparam::CRDP.DCFOp();
#ifdef CH_CTHR
ThreadTools::ThreadTeamArchitectChombo_t* const & CRDparam::g_threads =
  CRDparam::CRDP.threads();
#endif
CRD::Msg& CRD::msg = CRD::Msg::getInstance();


/*******************************************************************************
 */
/// Routine test_MMB_meshRefine
/*
 ******************************************************************************/

/*
 *   Block layout
 *
 *   ^ y (canonical j)
 *   |
 *   o-------o-------o
 *   |       |       |
 *   |   2   |   3   |
 *   |       |       |
 *   o-------o-------o
 *   |       |       |
 *   |   0   |   1   |
 *   |       |       |
 *   o-------o-------o--> x (canonical i)
 */

int
TESTFUNC()
{
  constexpr int prec = std::numeric_limits<Real>::digits10 - 1;
  int status = 0;

//--Set up the CS

  // CRDparam defined in ChordInput
  std::vector<int> AMRrefRatios (4, 4);
  CRDparam::CRDP.defineAMR(8, 16, AMRrefRatios, true);

  std::vector<int>  GRIDnumCells    (SpaceDim, 32);
  std::vector<Real> GRIDdomainLength(SpaceDim, 1.);
  std::vector<Real> GRIDdomainOrigin(SpaceDim, 0.);
  CRDparam::CRDP.defineGeometry(GRIDnumCells,
                                GRIDdomainLength,
                                GRIDdomainOrigin);

  // Factory use requires CRDP.defineAMR() and CRDP.defineGeometry
  AMRLevelCNSFactory amrLevelFactory;
  // Needs CRDP.defineAMR()
  amrLevelFactory.define();
  amrLevelFactory.setCoordinateSystemFactory(
    new DoubleCartesianCSFactory(true));
  IntVect domainLo = IntVect_zero;
  IntVect domainHi = stc::max(IntVect_zero,
                              CRDparam::g_domainBaseSize - IntVect_unit);
  bool periodic[SpaceDim];
  for (const int dir : EachDir)
    {
      periodic[dir] = true;
    }
  ProblemDomain problemDomain;
  problemDomain.define(domainLo, domainHi, periodic);
  amrLevelFactory.setBaseProblemDomain(problemDomain);
  CRDparam::CRDP.defineCoordSys(amrLevelFactory.new_MultiBlockCoordSys(0));

//--Set up the layout

  Vector<Box> boxes;
  for (const Box& box : CRDparam::g_coordSys->mappingBlocks())
    {
      boxes.push_back(box);
    }
  Vector<int> procIDs;
  LoadBalance(procIDs, boxes);
  DisjointBoxLayout dbl(boxes,
                        procIDs,
                        CRDparam::g_coordSys->levelDomain());
  dbl.defineLocalBlocks(CRDparam::g_coordSys->mappingBlocks(), true);
  // Commented out becuase repetitive with next test
  if (verbose)
    {
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          pout() << "Box " << dbl[dit] << " on proc " << dbl.procID(dit())
                 << " on block " << dbl.blockIndex(dit) << std::endl;
        }
    }

  auto MBCSfactory =
    [&factory = amrLevelFactory]
    (const int a_level)
      {
        return factory.new_MultiBlockCoordSys(a_level).operator->();
      };
  int gridBufferSize = LevelGridMetrics::bufferSize4thO(
    AMRrefRatios,
    1,
    1);
  std::cout << "Grid buffer size: " << gridBufferSize << std::endl;
  RefCountedPtr<CRDMeshRefine<decltype(MBCSfactory)>> mbmr =
    RefCountedPtr<CRDMeshRefine<decltype(MBCSfactory)>>(
      new CRDMeshRefine<decltype(MBCSfactory)>(
        CRDparam::g_coordSys->levelDomain(),
        AMRrefRatios,
        0.75,
        8,
        gridBufferSize,
        16,
        std::move(MBCSfactory)));

//--Refine to create level 1 (single tag at 15*IntVect_unit)

  int baseLevel = 0;
  int topLevel = 0;
  Vector<Vector<Box>> oldGrids(topLevel + 1);
  Vector<Vector<Box>> newGrids;
  Vector<IntVectSet> tags(topLevel + 1);
  
  {
    Vector<Box>& oldGrids0 = oldGrids[0];
    oldGrids0.reserve(dbl.size());
    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit)
      {
        oldGrids0.push_back(dbl[lit]);
      }
  }
  tags[0] |= IntVect(15*IntVect_unit);

  // Execute tested elements
  topLevel = mbmr->regrid(newGrids,
                          tags,
                          baseLevel,
                          topLevel,
                          oldGrids);
  std::cout << "Top level: " << topLevel << std::endl;

  // Check results
  {
    int lstat = 0;
    std::cout << "Size: " << newGrids.size() << std::endl;
    for (int level = 0; level != newGrids.size(); ++level)
      {
        std::cout << "On level " << level << std::endl;
        for (const Box& box : newGrids[level])
          {
            std::cout << "  " << box << std::endl;
          }
      }
    if (lstat) ++status;
  }
  std::cout << "====================================\n";

//--Refine to create level 2 (single tag at 64*IntVect_unit and
//--192*IntVect_unit)

  oldGrids = newGrids;
  newGrids.clear();
  tags.resize(2);
  tags[0].makeEmpty();
  tags[1] |= IntVect(63*IntVect_unit);

  // Execute tested elements
  topLevel = mbmr->regrid(newGrids,
                          tags,
                          baseLevel,
                          topLevel,
                          oldGrids);
  std::cout << "Top level: " << topLevel << std::endl;

  // Check results
  {
    int lstat = 0;
    std::cout << "Size: " << newGrids.size() << std::endl;
    for (int level = 0; level != newGrids.size(); ++level)
      {
        std::cout << "On level " << level << std::endl;
        for (const Box& box : newGrids[level])
          {
            std::cout << "  " << box << std::endl;
          }
      }
    if (lstat) ++status;
  }



  // oldGrids.resize(2);
  // tags.resize(2);
  // {
  //   Vector<Box>& oldGrids0 = oldGrids[0];
  //   oldGrids0.reserve(dbl.size());
  //   for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit)
  //     {
  //       oldGrids0.push_back(dbl[lit]);
  //     }
  //   Vector<Box>& oldGrids1 = oldGrids[1];
  //   for (const Box& box : newGrids[1])
  //     {
  //       oldGrids1 
  //     }
    
    
  // }

  return status;
}
