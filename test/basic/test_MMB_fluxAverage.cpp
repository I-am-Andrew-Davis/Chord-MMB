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
 * This test verifies the functionality of averaging the flux across mapped
 * multiblock boundaries.
 * Primary tested components are:
 *   - Chombo:MultiBlockRegion.{H|cpp}
 *   - Chord:MMBSingleLevel.{H|cpp}
 * Primary tested concepts are:
 *   - Exchanges and averaging of the flux
 *   - Block orientations
 *   - The use of a compact layout for caching the flux for exchanges
 *
 ******************************************************************************/

#define TESTFUNC test_MMB_fluxAverage

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
#include "DoubleCartesianRotateCS.H"
#include "LoadBalance.H"
#include "LevelGridMetrics.H"

#include "UsingNamespace.H"

//----- Internal -----//

#include "CRDparam.H"
#include "CRDmsg.H"
#include "AMRLevelCNSFactory.H"
#include "AMRLevelCNS.H"


/// Prototypes:
int
TESTFUNC();

int
testUnalignedMappings();

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
  ret += testUnalignedMappings();
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

//--Statics/globals in DataTemp

// No error checking here, but the static page size is compared against that
// used in the StackMem::alloc function where error checking is performed.
std::size_t StackMem::s_pageSize =
  static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
// Be aware that this may not allocate until used.
// This is 16 MB in 2D and 256 MB in 3D
thread_local StackMem tls_stack(16777216ull + (SpaceDim==3)*251658240);


/*******************************************************************************
 */
/// Routine test_MMB_fluxAverage
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
  std::vector<int> AMRrefRatios (1, 1);
  CRDparam::CRDP.defineAMR(8, 32, AMRrefRatios, true);

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
  // if (verbose)
  //   {
  //     for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
  //       {
  //         pout() << "Box " << dbl[dit] << " on proc " << dbl.procID(dit())
  //                << " on block " << dbl.blockIndex(dit) << std::endl;
  //       }
  //   }

//--Some data (two components, 1 ghost)

  constexpr int numComp = 2;

//--Set up the regions

  LevelGridMetrics lgm(2, 4);
  lgm.define(nullptr,
             CRDparam::g_coordSys,
             nullptr,
             RealVect_unit,
             IntVect_zero);

  MultiBlockRegions& mbRegions =
    const_cast<MultiBlockRegions&>(lgm.getMBRegions());
  mbRegions.define(*CRDparam::g_coordSys, dbl);

  MMBSingleLevel mbLevelOp(lgm);
  mbLevelOp.define(numComp);

/*
  First experiment.  U = 0.  F(0) = 1 everywhere.  F(1) is 0 on blocks 0 and 3
  and 1 on blocks 1 and 2

  Where [] is a cell and the numbers are fluxes, at a boundary 0[]01[]1, after
  averaging this becomes 0[]0.5[]1 and both left and right cells see an increase
  in the flux difference:
  left : 0[]0 -> 0[]0.5 or diff 0 -> diff 0.5
  right: 1[]1 -> 0.5[]1 or diff 1 -> diff 0.5
  So the state in both adjacent cells will decrease.  In other words
  01 -> leads to a decrease in state in both adacent cells
  10 -> leads to an increase in state in both adjacent cells
*/

  // Set the data and fluxes for the test
  LevelData<FArrayBox> data(dbl, numComp, IntVect_unit);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      data[dit].setVal(0.);
    }
  LevelData<FluxBox> flux(dbl, numComp, IntVect_zero);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      for (const int dir : EachDir)
        {
          flux[dit][dir].setVal(1., 0);
        }
      const int idxBlk = dbl.blockIndex(dit);
      CH_assert(idxBlk >= 0);
      Real val = 0.;
      if (idxBlk == 1 || idxBlk == 2)
        {
          val = 1.;
        }
      for (const int dir : EachDir)
        {
          flux[dit][dir].setVal(val, 1);
        }
    }

  // // Do we have boundaries?
  // auto& cs = *CRDparam::g_coordSys;
  // for (int idxBlk = 0, idxBlk_end = cs.numBlocks(); idxBlk != idxBlk_end;
  //      ++idxBlk)
  //   {
  //     pout() << "For block " << idxBlk << ":\n";
  //     for (auto side : EachSide)
  //       {
  //         for (const int dir : EachDir)
  //           {
  //             const BlockBoundary& bb = cs.boundary(idxBlk, dir, side);
  //             pout() << "  neighbor(" << dir << ',' << (int)side << "): "
  //                       << bb.neighbor() << std::endl;
  //           }
  //       }
  //   }

  // Execute tested elements
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box disjointBox = dbl[dit];
      mbLevelOp.setFlux(disjointBox,
                        flux[dit],
                        dit(),
                        Interval(0, numComp - 1));
    }
  mbLevelOp.exchange();
  mbLevelOp.refluxAverage(data, RealVect_unit, Interval(0, numComp - 1));

  // Check results
  {
    int lstat = 0;
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
      {
        const Box disjointBox = dbl[dit];
        const FArrayBox& dataFab = data[dit];
        // Check component 0 (0. everywhere)
        MD_BOXLOOP(disjointBox, i)
          {
            if (Misc::compare(dataFab[MD_IX(i, 0)], 0., prec))
              {
                ++lstat;
              }
          }
        // Check component 1 (0. at interior, +-0.5 or more on block boundaries)
        const int idxBlk = dbl.blockIndex(dit);
        CH_assert(idxBlk >= 0);
        const Box blkBox = CRDparam::g_coordSys->mappingBlocks()[idxBlk];
        Box testBlkBox = grow(blkBox, -1);
        // Interior
        const Box interiorBox = disjointBox & testBlkBox;
        MD_BOXLOOP(interiorBox, i)
          {
            if (Misc::compare(dataFab[MD_IX(i, 1)], 0., prec))
              {
                ++lstat;
              }
          }
        // Block faces and corners (in i and j directions)
        for (const int dir : EachDir)
          {
            if (dir > 1) break;
            testBlkBox.grow(dir, 1);
            for (const auto side : EachSide)
              {
                if (disjointBox.sideEnd(side)[dir] == blkBox.sideEnd(side)[dir])
                  {
                    // Face
                    Box faceBox = adjCellBox(disjointBox, dir, side, -1);
                    faceBox &= testBlkBox;  // Crop corners
                    Real val;
                    if (idxBlk == 1 || idxBlk == 2)
                      {
                        val = Side::sign(side)*0.5;
                      }
                    else
                      {
                        val = Side::sign(Side::flip(side))*0.5;
                      }
                    MD_BOXLOOP(faceBox, i)
                      {
                        if (Misc::compare(dataFab[MD_IX(i, 1)], val, prec))
                          {
                            ++lstat;
                          }
                      }
                    // Corners (may be visited more than once)
                    for (const int dir2 : EachDir)
                      {
                        if (dir2 > 1) break;
                        if (dir2 == dir) continue;
                        for (const auto side2 : EachSide)
                          {
                            if (disjointBox.sideEnd(side2)[dir2] ==
                                blkBox.sideEnd(side2)[dir2])
                              {
                                // Corners were cropped from the faceBox, so
                                // here we need the outer cells
                                Box cornerBox =
                                  adjCellBox(faceBox, dir2, side2, 1);
                                // Same rules for value, just add up
                                Real val;
                                if (idxBlk == 1 || idxBlk == 2)
                                  {
                                    val = (Side::sign(side) +
                                           Side::sign(side2))*0.5;
                                  }
                                else
                                  {
                                    val = (Side::sign(Side::flip(side)) +
                                           Side::sign(Side::flip(side2)))*0.5;
                                  }
                                MD_BOXLOOP(cornerBox, i)
                                  {
                                    if (Misc::compare(dataFab[MD_IX(i, 1)],
                                                      val, prec))
                                      {
                                        ++lstat;
                                      }
                                  }
                              }
                          }
                      }  // Corner directions
                  }
              }  // Loop over sides
            testBlkBox.grow(dir, -1);
          }  // Loop over directions
      }  // Loop over boxes
    if (lstat) ++status;
  }

/*
  This experiment is similar to the previous for component 1.  Instead, we use
  a more refined checkerboard flux pattern (4 patches per block) to make sure
  the flux data is received in buffers in the correct order.  This test is only
  on component 1.
*/

  // Sub-boxes for each mapping block
  CH_assert(CRDparam::g_coordSys->mappingBlocks().size() == 4);
  Box subBlk[4][4];  // [mapping block][sub box]
  for (int idxBlk = 0; idxBlk != 4; ++idxBlk)
    {
      const Box& blkBox = CRDparam::g_coordSys->mappingBlocks()[idxBlk];
      const int i0HalfSize = blkBox.size(0)/2;
      CH_assert(2*i0HalfSize == blkBox.size(0));
      const int i1HalfSize = blkBox.size(1)/2;
      CH_assert(2*i1HalfSize == blkBox.size(1));
      const IntVect lo = blkBox.smallEnd();
      subBlk[idxBlk][0] = blkBox;  // Bottom left
      subBlk[idxBlk][0].setBig  (0, lo[0] + i0HalfSize - 1);
      subBlk[idxBlk][0].setBig  (1, lo[1] + i1HalfSize - 1);
      subBlk[idxBlk][1] = blkBox;  // Bottom right
      subBlk[idxBlk][1].setSmall(0, lo[0] + i0HalfSize);
      subBlk[idxBlk][1].setBig  (1, lo[1] + i1HalfSize - 1);
      subBlk[idxBlk][2] = blkBox;  // Top left
      subBlk[idxBlk][2].setBig  (0, lo[0] + i0HalfSize - 1);
      subBlk[idxBlk][2].setSmall(1, lo[1] + i1HalfSize);
      subBlk[idxBlk][3] = blkBox;  // Top right
      subBlk[idxBlk][3].setSmall(0, lo[0] + i0HalfSize);
      subBlk[idxBlk][3].setSmall(1, lo[1] + i1HalfSize);
    }

  // Set the data and fluxes for the test
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      data[dit].setVal(0.);
    }
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box disjointBox = dbl[dit];
      for (const int dir : EachDir)
        {
          flux[dit][dir].setVal(1., 0);
        }
      const int idxBlk = dbl.blockIndex(dit);
      for (int idxSubBlk = 0; idxSubBlk != 4; ++idxSubBlk)
        {
          Box subBox = subBlk[idxBlk][idxSubBlk];
          subBox &= disjointBox;
          Real val = 0.;
          if (idxSubBlk == 1 || idxSubBlk == 2)
            {
              val = 1.;
            }
          for (const int dir : EachDir)
            {
              Box subBoxDir = subBox;
              subBoxDir.surroundingNodes(dir);
              flux[dit][dir].setVal(val, subBoxDir, 1);
            }
        }
    }

  // Execute tested elements
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box disjointBox = dbl[dit];
      mbLevelOp.setFlux(disjointBox,
                        flux[dit],
                        dit(),
                        Interval(0, numComp - 1));
    }
  mbLevelOp.exchange();
  mbLevelOp.refluxAverage(data, RealVect_unit, Interval(0, numComp - 1));

  // Check results
  {
    int lstat = 0;
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
      {
        const Box disjointBox = dbl[dit];
        const FArrayBox& dataFab = data[dit];
        // Check component 0 (0. everywhere)
        MD_BOXLOOP(disjointBox, i)
          {
            if (Misc::compare(dataFab[MD_IX(i, 0)], 0., prec))
              {
                ++lstat;
              }
          }
        // Check component 1 (0. at interior, +-0.5 or more on block boundaries)
        const int idxBlk = dbl.blockIndex(dit);
        CH_assert(idxBlk >= 0);
        const Box blkBox = CRDparam::g_coordSys->mappingBlocks()[idxBlk];
        Box testBlkBox = grow(blkBox, -1);
        // Interior
        const Box interiorBox = disjointBox & testBlkBox;
        MD_BOXLOOP(interiorBox, i)
          {
            if (Misc::compare(dataFab[MD_IX(i, 1)], 0., prec))
              {
                ++lstat;
              }
          }
        // Block faces and corners (in i and j directions)
        for (const int dir : EachDir)
          {
            if (dir > 1) break;
            testBlkBox.grow(dir, 1);
            for (const auto side : EachSide)
              {
                if (disjointBox.sideEnd(side)[dir] == blkBox.sideEnd(side)[dir])
                  {
                    // Face
                    Box faceBox = adjCellBox(disjointBox, dir, side, -1);
                    faceBox &= testBlkBox;  // Crop corners
                    // Run through all sub blocks
                    for (int idxSubBlk = 0; idxSubBlk != 4; ++idxSubBlk)
                      {
                        Box subBox = subBlk[idxBlk][idxSubBlk];
                        if (subBox.sideEnd(side)[dir] !=
                            blkBox.sideEnd(side)[dir])
                          {
                            // Skip sub blocks that do not share the edge under
                            // consideration in the primary block
                            continue;
                          }
                        Box subFaceBox = faceBox & subBox;
                        Real val;
                        if (idxSubBlk == 1 || idxSubBlk == 2)
                          {
                            val = Side::sign(side)*0.5;
                          }
                        else
                          {
                            val = Side::sign(Side::flip(side))*0.5;
                          }
                        MD_BOXLOOP(subFaceBox, i)
                          {
                            if (Misc::compare(dataFab[MD_IX(i, 1)], val, prec))
                              {
                                ++lstat;
                              }
                          }
                        // Corners (may be visited more than once)
                        for (const int dir2 : EachDir)
                          {
                            if (dir2 > 1) break;
                            if (dir2 == dir) continue;
                            for (const auto side2 : EachSide)
                              {
                                if (
                                  // Edge of block
                                  disjointBox.sideEnd(side2)[dir2] ==
                                  blkBox.sideEnd(side2)[dir2] &&
                                  // Edge of sub block
                                  disjointBox.sideEnd(side2)[dir2] ==
                                  subBox.sideEnd(side2)[dir2] )
                                  {
                                    // Corners were cropped from the faceBox, so
                                    // here we need the outer cells
                                    Box cornerBox =
                                      adjCellBox(faceBox, dir2, side2, 1);
                                    // Same rules for value, just add up
                                    Real val;
                                    if (idxSubBlk == 1 || idxSubBlk == 2)
                                      {
                                        val = (Side::sign(side) +
                                               Side::sign(side2))*0.5;
                                      }
                                    else
                                      {
                                        val =
                                          (Side::sign(Side::flip(side)) +
                                           Side::sign(Side::flip(side2)))*0.5;
                                      }
                                    // lstat = 0;  // DEBUG
                                    MD_BOXLOOP(cornerBox, i)
                                      {
                                        if (Misc::compare(dataFab[MD_IX(i, 1)],
                                                          val, prec))
                                          {
                                            ++lstat;
                                          }
                                      }
                                  }
                              }
                          }  // Corner directions
                      }  // Loop over sub blocks
                  }
              }  // Loop over sides
            testBlkBox.grow(dir, -1);
          }  // Loop over directions
      }  // Loop over boxes
    if (lstat) ++status;
  }

  return status;
}

/*
 *   Block layout visualized with physical orientations (computation-space
 *   indices have been rotated as shown)
 *
 *   ^ y (physical space j)
 *   |       .      .        (47,47)   
 *   o-------o       o-------o
 *   | +-->j |       | i<--+ |
 *   |i| 2   |       |   3 |j|
 *   | v     |       |     v |
 *   o-------o       o-------o
 *                   (32,32)
 *
 *           (15,15)
 *   o-------o       o-------o
 *   | ^     |       |     ^ |
 *   |j| 0   |       |   1 |i|
 *   | +-->i |       | j<--+ |
 *   o-------o       o-------o--> x (physical space i)
 *   (0,0)
 * 
 *   In Xi and X coordinates, boxes have unit size.  16 cells are used per box.
 *   They are separated by a unit jump (16 cells).
 *
 *   Some coordinates in physical orientation are given above.  For each block,
 *   the coordinates in computation space have the same range but are reoriented
 *   according to the drawn origin.  I.e., (15,47) in block 2 points to the
 *   center of the domain which has physical coordinates (1.0,1.0) (roughly,
 *   need either cell center +0.5 conversions or advance to the node).
 *
 *   The dot at top represents the point (-1,47) in block 2 which becomes
 *   (48,31) in block 3.  Note that these locations are with respect to the
 *   orientations of the individual blocks.  Both have X=(0.96875,2.03125).
 *   This information is printed for block 2 if verbose.  This info is mostly
 *   for a sanity check because it is difficult to visualize these rotations.
 */

int
testUnalignedMappings()
{
  constexpr int prec = std::numeric_limits<Real>::digits10 - 1;
  int status = 0;

//--Set up the CS

  // CRDparam defined in ChordInput
  std::vector<int> AMRrefRatios (1, 1);
  CRDparam::CRDP.defineAMR(8, 32, AMRrefRatios, true);

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
    new DoubleCartesianRotateCSFactory(true));
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
  if (verbose)
    {
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          pout() << "Box " << dbl[dit] << " on proc " << dbl.procID(dit())
                 << " on block " << dbl.blockIndex(dit) << std::endl;
        }
    }

//--Some data (two components, 1 ghost)

  constexpr int numComp = 2;

//--Set up the regions

  LevelGridMetrics lgm(2, 4);
  lgm.define(nullptr,
             CRDparam::g_coordSys,
             nullptr,
             RealVect_unit,
             IntVect_zero);

  MultiBlockRegions& mbRegions =
    const_cast<MultiBlockRegions&>(lgm.getMBRegions());
  mbRegions.define(*CRDparam::g_coordSys, dbl);

  MMBSingleLevel mbLevelOp(lgm);
  mbLevelOp.define(numComp);

  // Print some info
  if (verbose)
    {
      const RealVect levelDx = RealVect_unit/16;
      //amrLevelFactory.getLevelDx(0);
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          const Box& disjointBox = dbl[dit];
          const int idxBlk = dbl.blockIndex(dit);
          RealVect loXi = levelDx*disjointBox.smallEnd();
          RealVect loX = CRDparam::g_coordSys->getCoordSys(idxBlk)
            ->realCoord(loXi);
          pout() << "Blk: " << idxBlk << ", lo: " << loXi << " --> " << loX
                 << std::endl;
          RealVect hiXi = levelDx*(disjointBox.bigEnd() + IntVect_unit);
          RealVect hiX = CRDparam::g_coordSys->getCoordSys(idxBlk)
            ->realCoord(hiXi);
          pout() << "Blk: " << idxBlk << ", hi: " << hiXi << " --> " << hiX
                 << std::endl;
          for (const int dir : EachDir)
            {
              for (const auto side : EachSide)
                {
                  pout() << "Dir: " << dir << ", side; " << side << std::endl;
                  BlockBoundary bb =
                    CRDparam::g_coordSys->boundary(idxBlk, dir, side);
                  pout() << "  nbridx : " << bb.neighbor() << std::endl;
                  pout() << "  nbrdir0: " << bb.dirOther(0) << std::endl;
                  pout() << "  nbrsgn0: " << bb.reorientFace(0) << std::endl;
                  pout() << "  nbrdir1: " << bb.dirOther(1) << std::endl;
                  pout() << "  nbrsgn1: " << bb.reorientFace(1) << std::endl;
                  if (idxBlk == 2 && dir == 1 && side == Side::Hi)
                    {
                      // In realspace (0.96875 2.03125)
                      IntVect i2{-1, 47, 0};
                      // In block 3, this is (48, 31)
                      IntVect i3 = bb.getTransformation().transformFwd(i2);
                      pout() << "Blk2: "
                             << CRDparam::g_coordSys->getCoordSys(idxBlk)
                        ->realCoord(levelDx*(RealVect(i2) + 0.5))
                             << std::endl;
                      pout() << "Blk3: "
                             << CRDparam::g_coordSys->getCoordSys(3)
                        ->realCoord(levelDx*(RealVect(i3) + 0.5)) << std::endl;
                    }
                }
            }
        }
    }

  // Test the above and also the rigid transforms
  {
    int lstat = 0;
    const RealVect levelDx = RealVect_unit/16;
    RealVect answersX[4][2] = {
      //                 x   y
      { RealVect{D_DECL6(0., 0., 0., 0., 0., 0.)},    // Blk 0, X_lo
        RealVect{D_DECL6(1., 1., 0., 0., 0., 0.)} },  // Blk 0, X_hi
      { RealVect{D_DECL6(2., 0., 0., 0., 0., 0.)},    // Blk 1, X_lo
        RealVect{D_DECL6(1., 1., 0., 0., 0., 0.)} },  // Blk 1, X_hi
      { RealVect{D_DECL6(0., 2., 0., 0., 0., 0.)},    // Blk 2, X_lo
        RealVect{D_DECL6(1., 1., 0., 0., 0., 0.)} },  // Blk 2, X_hi
      { RealVect{D_DECL6(2., 2., 0., 0., 0., 0.)},    // Blk 3, X_lo
        RealVect{D_DECL6(1., 1., 0., 0., 0., 0.)} }   // Blk 3, X_hi
    };
    int answersNbrIdx[4][2] = {
      // dir 0
      // |  dir 1
      {  1, 2 },  // Blk 0
      {  3, 0 },  // Blk 1
      {  0, 3 },  // Blk 2
      {  2, 1 }   // Blk 3
    };
    int answersNbrCS[2][4] = {  // Same for all blocks
      // nbrdir0
      // |   nbrsgn0
      // |   |  nbrdir1
      // |   |  |   nbrsgn1
      {  1, -1, 0,  1 },  // dir 0
      {  1,  1, 0, -1 },  // dir 1
    };
    RealVect answersLoSideRigTrm[4][2] = {
      // Blk 0
      { RealVect{D_DECL6( 1.,  0., 0., 0., 0., 0.)},  // dir 0
        RealVect{D_DECL6( 0.,  1., 0., 0., 0., 0.)}   // dir 1
      },
      // Blk 1
      { RealVect{D_DECL6( 0.,  1., 0., 0., 0., 0.)},  // dir 0
        RealVect{D_DECL6(-1.,  0., 0., 0., 0., 0.)}   // dir 1
      },
      // Blk 2
      { RealVect{D_DECL6( 0., -1., 0., 0., 0., 0.)},  // dir 0
        RealVect{D_DECL6( 1.,  0., 0., 0., 0., 0.)}   // dir 1
      },
      // Blk 3
      { RealVect{D_DECL6(-1.,  0., 0., 0., 0., 0.)},  // dir 0
        RealVect{D_DECL6( 0., -1., 0., 0., 0., 0.)}   // dir 1
      },
    };

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
      {
        const Box& disjointBox = dbl[dit];
        const int idxBlk = dbl.blockIndex(dit);
        RealVect loXi = levelDx*disjointBox.smallEnd();
        RealVect loX = CRDparam::g_coordSys->getCoordSys(idxBlk)
          ->realCoord(loXi);
        if (Misc::compare(loX[0], answersX[idxBlk][0][0], prec)) ++lstat;
        if (Misc::compare(loX[1], answersX[idxBlk][0][1], prec)) ++lstat;
        RealVect hiXi = levelDx*(disjointBox.bigEnd() + IntVect_unit);
        RealVect hiX = CRDparam::g_coordSys->getCoordSys(idxBlk)
          ->realCoord(hiXi);
        if (Misc::compare(hiX[0], answersX[idxBlk][1][0], prec)) ++lstat;
        if (Misc::compare(hiX[1], answersX[idxBlk][1][1], prec)) ++lstat;
        for (const int dir : EachDir)
          {
            for (const auto side : EachSide)
              {
                BlockBoundary bb =
                  CRDparam::g_coordSys->boundary(idxBlk, dir, side);
                if (bb.neighbor() != answersNbrIdx[idxBlk][dir]) ++lstat;
                if (bb.dirOther(0)     != answersNbrCS[dir][0]) ++lstat;
                if (bb.reorientFace(0) != answersNbrCS[dir][1]) ++lstat;
                if (bb.dirOther(1)     != answersNbrCS[dir][2]) ++lstat;
                if (bb.reorientFace(1) != answersNbrCS[dir][3]) ++lstat;
                switch (side)
                  {
                  case Side::Lo:
                    if (Misc::compare(
                          bb.getPhysTransformation().getTranslation()[0],
                          answersLoSideRigTrm[idxBlk][dir][0], prec)) ++lstat;
                    if (Misc::compare(
                          bb.getPhysTransformation().getTranslation()[1],
                          answersLoSideRigTrm[idxBlk][dir][1], prec)) ++lstat;
                    break;
                  case Side::Hi:
                    if (Misc::compare(
                          bb.getPhysTransformation().getTranslation()[0],
                          0., prec)) ++lstat;
                    if (Misc::compare(
                          bb.getPhysTransformation().getTranslation()[1],
                          0., prec)) ++lstat;
                    break;
                  default:
                    MayDay::Error("Invalid side");
                  }
              }
          }
      }
    if (lstat) ++status;
  }

/*
  First experiment.  U = 0.  F(0) = +-1, adjusted so that flux is in same
  direction depending on block orientation.  F(1) is 1 everywhere.

  In the first case, all fluxes should cancel.  In the second, flux accumulates
  at +1 on high-sides of blocks and -1 on low sides.

  For component 1, the fluxes look something like the following.  Flux
  directions are shown in both spaces at boundaries

  Computational space:
  ^ y
  |
  b-------d       e-------d
  |>^^^^^>|       |>^^^^^>|
  |>  1  >|       |>  1  >|
  |>^^^^^>|       |>^^^^^>|
  a-------c       a-------b



  c-------d       b-------d
  |>^^^^^>|       |>^^^^^>|
  |>  1  >|       |>  1  >|
  |>^^^^^>|       |>^^^^^>|
  a-------b       a-------e--> x

  Physical space:
  ^ y
  |
  a-------b       b-------a
  |>vvvvv>|       |<vvvvv<|
  |>  1  >|       |<  1  <|
  |>vvvvv>|       |<vvvvv<|
  c-------d       d-------e



  c-------d       d-------e
  |>^^^^^>|       |<^^^^^<|
  |>  1  >|       |<  1  <|
  |>^^^^^>|       |<^^^^^<|
  a-------b       b-------a--> x
 
*/

  // Set the data and fluxes for the test
  LevelData<FArrayBox> data(dbl, numComp, IntVect_unit);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      data[dit].setVal(0.);
    }
  LevelData<FluxBox> flux(dbl, numComp, IntVect_zero);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const int idxBlk = dbl.blockIndex(dit);
      CH_assert(idxBlk >= 0);
      const Box disjointBox = dbl[dit];
      for (const int dir : EachDir)
        {
          flux[dit][dir].setVal(0.);
        }
      // Counter-balance fluxes with directions
      switch (idxBlk)
        {
        case 0:  // Bottom left
          flux[dit][0].setVal( 1., 0);
          flux[dit][1].setVal( 1., 0);
          break;
        case 1:  // Bottom right
          flux[dit][0].setVal( 1., 0);
          flux[dit][1].setVal(-1., 0);
          break;
        case 2:  // Top left
          flux[dit][0].setVal(-1., 0);
          flux[dit][1].setVal( 1., 0);
          break;
        case 3:  // Top right
          flux[dit][0].setVal(-1., 0);
          flux[dit][1].setVal(-1., 0);
          break;
        }
      // All flux set to 1 (squeezes flux to physical interior)
      flux[dit][0].setVal(1., 1);
      flux[dit][1].setVal(1., 1);
    }

  // Execute tested elements
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box disjointBox = dbl[dit];
      mbLevelOp.setFlux(disjointBox,
                        flux[dit],
                        dit(),
                        Interval(0, numComp - 1));
    }
  mbLevelOp.exchange();
  mbLevelOp.refluxAverage(data, RealVect_unit, Interval(0, numComp - 1));

  // Check results
  {
    int lstat = 0;
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
      {
        const Box disjointBox = dbl[dit];
        const FArrayBox& dataFab = data[dit];
        // Check component 0 (0. everywhere)
        MD_BOXLOOP(disjointBox, i)
          {
            if (Misc::compare(dataFab[MD_IX(i, 0)], 0., prec))
              {
                ++lstat;
              }
          }
        // Check component 1 (-1 on low side, 1 on high side block boundaries)
        const int idxBlk = dbl.blockIndex(dit);
        CH_assert(idxBlk >= 0);
        const Box blkBox = CRDparam::g_coordSys->mappingBlocks()[idxBlk];
        Box testBlkBox = grow(blkBox, -1);
        // Interior
        const Box interiorBox = disjointBox & testBlkBox;
        MD_BOXLOOP(interiorBox, i)
          {
            if (Misc::compare(dataFab[MD_IX(i, 1)], 0., prec))
              {
                ++lstat;
              }
          }
        // Block faces and corners (in i and j directions)
        for (const int dir : EachDir)
          {
            if (dir > 1) break;
            testBlkBox.grow(dir, 1);
            for (const auto side : EachSide)
              {
                if (disjointBox.sideEnd(side)[dir] == blkBox.sideEnd(side)[dir])
                  {
                    // Face
                    Box faceBox = adjCellBox(disjointBox, dir, side, -1);
                    faceBox &= testBlkBox;  // Crop corners
                    Real val = Side::sign(side);
                    MD_BOXLOOP(faceBox, i)
                      {
                        if (Misc::compare(dataFab[MD_IX(i, 1)], val, prec))
                          {
                            ++lstat;
                          }
                      }
                    // Corners (may be visited more than once)
                    for (const int dir2 : EachDir)
                      {
                        if (dir2 > 1) break;
                        if (dir2 == dir) continue;
                        for (const auto side2 : EachSide)
                          {
                            if (disjointBox.sideEnd(side2)[dir2] ==
                                blkBox.sideEnd(side2)[dir2])
                              {
                                // Corners were cropped from the faceBox, so
                                // here we need the outer cells
                                Box cornerBox =
                                  adjCellBox(faceBox, dir2, side2, 1);
                                // Same rules for value, just add up
                                Real val = Side::sign(side) + Side::sign(side2);
                                MD_BOXLOOP(cornerBox, i)
                                  {
                                    if (Misc::compare(dataFab[MD_IX(i, 1)],
                                                      val, prec))
                                      {
                                        ++lstat;
                                      }
                                  }
                              }
                          }
                      }  // Corner directions
                  }
              }  // Loop over sides
            testBlkBox.grow(dir, -1);
          }  // Loop over directions
      }  // Loop over boxes
    if (lstat) ++status;
  }

/*
  Second experiment.  This is similar to the previous for component 1.  Instead,  we use a more refined checkerboard flux pattern (4 patches per block) to make
  sure the flux data is received in buffers in the correct order.  This test is
  only on component 1.

  For component 1, the fluxes look something like the following.  Flux
  directions are shown in both spaces

  Computational space:
  ^ y
  |
  o-------o       o-------o
  |>1^|   |       |>1^|   |
  |---+---|       |---+---|
  |   |>1^|       |   |>1^|
  o-------o       o-------o



  o-------o       o-------o
  |>1^|   |       |>1^|   |
  |---+---|       |---+---|
  |   |>1^|       |   |>1^|
  o-------o       o-------o--> x

  Physical space:
  ^ y
  |
  o-------o       o-------o
  |   |>1v|       |<1v|   |
  |---+---|       |---+---|
  |>1v|   |       |   |<1v|0.5
  o-------o       o-------o



  o-------o       o-------o
  |>1^|   |       |   |<1^|
  |---+---|       |---+---|
  |   |>1^|       |<1^|   |
  o-------o       o-------o--> x
*/

  // Sub-boxes for each mapping block
  CH_assert(CRDparam::g_coordSys->mappingBlocks().size() == 4);
  Box subBlk[4][4];  // [mapping block][sub box]
  for (int idxBlk = 0; idxBlk != 4; ++idxBlk)
    {
      const Box& blkBox = CRDparam::g_coordSys->mappingBlocks()[idxBlk];
      const int i0HalfSize = blkBox.size(0)/2;
      CH_assert(2*i0HalfSize == blkBox.size(0));
      const int i1HalfSize = blkBox.size(1)/2;
      CH_assert(2*i1HalfSize == blkBox.size(1));
      const IntVect lo = blkBox.smallEnd();
      subBlk[idxBlk][0] = blkBox;  // Bottom left
      subBlk[idxBlk][0].setBig  (0, lo[0] + i0HalfSize - 1);
      subBlk[idxBlk][0].setBig  (1, lo[1] + i1HalfSize - 1);
      subBlk[idxBlk][1] = blkBox;  // Bottom right
      subBlk[idxBlk][1].setSmall(0, lo[0] + i0HalfSize);
      subBlk[idxBlk][1].setBig  (1, lo[1] + i1HalfSize - 1);
      subBlk[idxBlk][2] = blkBox;  // Top left
      subBlk[idxBlk][2].setBig  (0, lo[0] + i0HalfSize - 1);
      subBlk[idxBlk][2].setSmall(1, lo[1] + i1HalfSize);
      subBlk[idxBlk][3] = blkBox;  // Top right
      subBlk[idxBlk][3].setSmall(0, lo[0] + i0HalfSize);
      subBlk[idxBlk][3].setSmall(1, lo[1] + i1HalfSize);
    }

  // Set the data and fluxes for the test
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      data[dit].setVal(0.);
    }
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const int idxBlk = dbl.blockIndex(dit);
      CH_assert(idxBlk >= 0);
      const Box disjointBox = dbl[dit];
      for (const int dir : EachDir)
        {
          flux[dit][dir].setVal(0.);
        }
      for (int idxSubBlk = 0; idxSubBlk != 4; ++idxSubBlk)
        {
          Box subBox = subBlk[idxBlk][idxSubBlk];
          subBox &= disjointBox;
          Real val = 0.;
          if (idxSubBlk == 1 || idxSubBlk == 2)
            {
              val = 1.;
            }
          for (const int dir : EachDir)
            {
              Box subBoxDir = subBox;
              subBoxDir.surroundingNodes(dir);
              flux[dit][dir].setVal(val, subBoxDir, 1);
            }
        }
    }

  // Execute tested elements
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box disjointBox = dbl[dit];
      mbLevelOp.setFlux(disjointBox,
                        flux[dit],
                        dit(),
                        Interval(0, numComp - 1));
    }
  mbLevelOp.exchange();
  mbLevelOp.refluxAverage(data, RealVect_unit, Interval(0, numComp - 1));

  // Check results
  {
    int lstat = 0;
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
      {
        const Box disjointBox = dbl[dit];
        const FArrayBox& dataFab = data[dit];
        // Check component 0 (0. everywhere)
        MD_BOXLOOP(disjointBox, i)
          {
            if (Misc::compare(dataFab[MD_IX(i, 0)], 0., prec))
              {
                ++lstat;
              }
          }
        // Check component 1 (-1 on low side, 1 on high side block boundaries)
        const int idxBlk = dbl.blockIndex(dit);
        CH_assert(idxBlk >= 0);
        const Box blkBox = CRDparam::g_coordSys->mappingBlocks()[idxBlk];
        Box testBlkBox = grow(blkBox, -1);
        // Interior
        const Box interiorBox = disjointBox & testBlkBox;
        MD_BOXLOOP(interiorBox, i)
          {
            if (Misc::compare(dataFab[MD_IX(i, 1)], 0., prec))
              {
                ++lstat;
              }
          }
        // Block faces and corners (in i and j directions)
        for (const int dir : EachDir)
          {
            if (dir > 1) break;
            testBlkBox.grow(dir, 1);
            for (const auto side : EachSide)
              {
                if (disjointBox.sideEnd(side)[dir] == blkBox.sideEnd(side)[dir])
                  {
                    // Face
                    Box faceBox = adjCellBox(disjointBox, dir, side, -1);
                    faceBox &= testBlkBox;  // Crop corners
                    // Run through all sub blocks
                    for (int idxSubBlk = 0; idxSubBlk != 4; ++idxSubBlk)
                      {
                        Box subBox = subBlk[idxBlk][idxSubBlk];
                        if (subBox.sideEnd(side)[dir] !=
                            blkBox.sideEnd(side)[dir])
                          {
                            // Skip sub blocks that do not share the edge under
                            // consideration in the primary block
                            continue;
                          }
                        Box subFaceBox = faceBox & subBox;
                        Real val = 0.;
                        if (idxSubBlk == 1 || idxSubBlk == 2)
                          {
                            val = Side::sign(side);
                          }
                        MD_BOXLOOP(subFaceBox, i)
                          {
                            if (Misc::compare(dataFab[MD_IX(i, 1)], val, prec))
                              {
                                ++lstat;
                              }
                          }
                        // Corners (may be visited more than once).  These are
                        // always 0 because the ones adjacent to changing fluxes
                        // are always adjacent to both low and high faces.
                        for (const int dir2 : EachDir)
                          {
                            if (dir2 > 1) break;
                            if (dir2 == dir) continue;
                            for (const auto side2 : EachSide)
                              {
                                if (
                                  // Edge of block
                                  disjointBox.sideEnd(side2)[dir2] ==
                                  blkBox.sideEnd(side2)[dir2] &&
                                  // Edge of sub block
                                  disjointBox.sideEnd(side2)[dir2] ==
                                  subBox.sideEnd(side2)[dir2] )
                                  {
                                    // Corners were cropped from the faceBox, so
                                    // here we need the outer cells
                                    Box cornerBox =
                                      adjCellBox(faceBox, dir2, side2, 1);
                                    // lstat = 0;  // DEBUG
                                    MD_BOXLOOP(cornerBox, i)
                                      {
                                        if (Misc::compare(dataFab[MD_IX(i, 1)],
                                                          0., prec))
                                          {
                                            ++lstat;
                                          }
                                      }
                                  }
                              }
                          }  // Corner directions
                      }  // Loop over sub blocks
                  }
              }  // Loop over sides
            testBlkBox.grow(dir, -1);
          }  // Loop over directions
      }  // Loop over boxes
    if (lstat) ++status;
  }

/*
  Second experiment.  Same as first, but add some crazy rotations
*/

  return status;
}
