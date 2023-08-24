#ifdef CH_LANG_CC
/*
 *      _______               __
 *     / ___/ /  ___  __  ___/ /
 *    / /__/ _ \/ _ \/ _\/ _  /
 *    \___/_//_/\___/_/  \_._/
 *    Please refer to Copyright.txt, in Chord's root directory.
 */
#endif

#include <iostream>
#include <cstring>
#include <limits>

#ifdef CH_MPI
#include <mpi.h>
#endif

#include "RealVect.H"
#include "CONSTANTS.H"
#include "parstream.H"
#include "parseTestOptions.H"
#include "DataTemp.H"
#include "ProblemDomain.H"
#include "LoHiCenter.H"
#include "GodunovUtilitiesF_F.H"
#include "BaseFabMacros.H"
#include "PatchMappedFunc.H"
#include "CRDutil.H"
#include "CRDparam.H"
#include "CRDmsg.H"

#include "UsingNamespace.H"

CRDparam::CRDparamVar CRDparam::CRDP;
const CRDPhysics* const& CRDparam::g_CRDPhysics = CRDparam::CRDP.CRDPhysicsOp();
const CNSIBC* const& CRDparam::g_CNSIBC = CRDparam::CRDP.CNSIBCOp();
const DCFlattening* const& CRDparam::g_DCF = CRDparam::CRDP.DCFOp();
#ifdef CH_CTHR
ThreadTools::ThreadTeamArchitectChombo_t* const & CRDparam::g_threads =
  CRDparam::CRDP.threads();
#endif
CRD::Msg& CRD::msg = CRD::Msg::getInstance();

/// Prototypes:
int
test_PatchMappedFunc_FaceAvgNormDerivFromCellAvg();

/// Global variables for handling output:
static const char *pgmname = "test_PatchMappedFunc_FaceAvgNormDerivFromCellAvg";
static const char *indent = "   ", *indent2 = "      ";
static bool verbose = true;

int
main(int argc, const char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions(argc, argv, verbose);

  if (verbose)
    {
      pout() << indent2 << "Beginning " << pgmname << " ..." << std::endl;
    }

//--Run the tests

  int ret = test_PatchMappedFunc_FaceAvgNormDerivFromCellAvg();
  if (ret == 0)
    {
      pout() << indent << pgmname << " passed all tests" << std::endl;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)"
             << std::endl;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}


/*******************************************************************************
 */
/// Routine test_CRDutil_Laplacian
/*
 ******************************************************************************/

int
test_PatchMappedFunc_FaceAvgNormDerivFromCellAvg()
{
  CH_TIMERS("test_PatchMappedFunc_FaceAvgNormDerivFromCellAvg");
  const int prec = std::numeric_limits<Real>::digits10 - 1;
  int status = 0;
//--Initial set up
//  We need to test several pieces of this function
//  1) periodic/non-periodic domain
//  2) order-of-accuracy

//--Initial conditions
//  The first IC is a linear profile. Using 2nd or higher order of accuracy
//  the derivative of this IC should be captured exactly. Additionally, it
//  should be a constant value, so we can easily test for issues
//
//  The second IC is a product of sine waves. This is set to
//  f(x,y,z) = sin(2*PI*w_x*x)*sin(2*PI*w_y*y)*sin(2*PI*w_z*z)
//  This way, the derivatives are
//  f_x = 2*PI*w_x*cos(2*PI*w_x*x)*sin(2*PI*w_y*y)*sin(2*PI*w_z*z)
//  f_y = 2*PI*w_y*sin(2*PI*w_x*x)*cos(2*PI*w_y*y)*sin(2*PI*w_z*z)
//  f_z = 2*PI*w_z*sin(2*PI*w_x*x)*sin(2*PI*w_y*y)*cos(2*PI*w_z*z)
//  Using this IC, the order of accuracy of the operator can be checked.

  constexpr int n = 16;
  const Box box(IntVect_unit*(-n/2), IntVect_unit*((n/2) - 1));
  const Box box1 = grow(box, 1);
  const Box box2 = grow(box, 2);

  FArrayBox linearFab(box, 1);
  FArrayBox linearFab1(box1, 1);
  FArrayBox linearFab2(box2, 1);
  const bool periodic[3]    = {true, true, true};
  const bool notPeriodic[3] = {false, false, false};
  const ProblemDomain pDomain(box, periodic);
  const ProblemDomain npDomain(box, notPeriodic);

  const RealVect dx = RealVect_unit/n;
  MD_BOXLOOP(box, i)
    {
      const RealVect loX = dx*RealVect{MD_EXPANDIX((Real)i)};
      const RealVect hiX = dx*(RealVect{MD_EXPANDIX((Real)i)} + RealVect::Unit);
      linearFab[MD_IX(i, 0)] = stc::sum(RealVect{1., 1.5, 2.}*(loX + hiX));
    }
  MD_BOXLOOP(box1, i)
    {
      const RealVect loX = dx*RealVect{MD_EXPANDIX((Real)i)};
      const RealVect hiX = dx*(RealVect{MD_EXPANDIX((Real)i)} + RealVect::Unit);
      linearFab1[MD_IX(i, 0)] = stc::sum(RealVect{1., 1.5, 2.}*(loX + hiX));
    }
  MD_BOXLOOP(box2, i)
    {
      const RealVect loX = dx*RealVect{MD_EXPANDIX((Real)i)};
      const RealVect hiX = dx*(RealVect{MD_EXPANDIX((Real)i)} + RealVect::Unit);
      linearFab2[MD_IX(i, 0)] = stc::sum(RealVect{1., 1.5, 2.}*(loX + hiX));
    }

/*--------------------------------------------------------------------*
 * Test normal derivative of linear data field
 *--------------------------------------------------------------------*/

  FLUXBOXSTACKTEMP(faceAvgNormDerivFxb, box, 1);
  const RealVect exactDerivVec(D_DECL(2., 3., 4.));
  // second-order interior computation
  {
    CH_TIME("test::secondOrderInteriorScheme");
    for (int dir = 0; dir != SpaceDim; ++dir)
      {
        FArrayBox& faceAvgNormDerivFab = faceAvgNormDerivFxb[dir];
        faceAvgNormDerivFab.setVal(0.);
        Box faceBox = surroundingNodes(box, dir);
        PatchMappedFunc::faceAvgNormDerivFromCellAvgCS(faceAvgNormDerivFab,
                                                       linearFab1,
                                                       0, // input comp
                                                       0, // output comp
                                                       pDomain,
                                                       faceBox,
                                                       dx,
                                                       dir,
                                                       false, // second-order
                                                       true); // interior
        int lstat = 0;
        const Real exactDeriv = exactDerivVec[dir];
        MD_BOXLOOP(faceBox, i)
          {
            const Real numericalDeriv = faceAvgNormDerivFab[MD_IX(i, 0)];
            if (Misc::compare1(numericalDeriv, exactDeriv, prec))
              {
                ++lstat;
              }
          }
        if (lstat) ++status;
      }
  }
  // fourth-order interior computation
  {
    CH_TIME("test::fourthOrderInteriorScheme");
    for (int dir = 0; dir != SpaceDim; ++dir)
      {
        FArrayBox& faceAvgNormDerivFab = faceAvgNormDerivFxb[dir];
        faceAvgNormDerivFab.setVal(0.);
        Box faceBox = surroundingNodes(box, dir);
        PatchMappedFunc::faceAvgNormDerivFromCellAvgCS(faceAvgNormDerivFab,
                                                       linearFab2,
                                                       0, // input comp
                                                       0, // output comp
                                                       pDomain,
                                                       faceBox,
                                                       dx,
                                                       dir,
                                                       true,  // fourth-order
                                                       true); // interior
        int lstat = 0;
        const Real exactDeriv = exactDerivVec[dir];
        MD_BOXLOOP(faceBox, i)
          {
            const Real numericalDeriv = faceAvgNormDerivFab[MD_IX(i, 0)];
            if (Misc::compare1(numericalDeriv, exactDeriv, prec))
              {
                ++lstat;
              }
          }
        if (lstat) ++status;
      }
  }
  // second-order bounded-domain computation
  {
    CH_TIME("test::secondOrderBoundaryScheme");
    for (int dir = 0; dir != SpaceDim; ++dir)
      {
        FArrayBox& faceAvgNormDerivFab = faceAvgNormDerivFxb[dir];
        faceAvgNormDerivFab.setVal(0.);
        Box faceBox = surroundingNodes(box, dir);
        PatchMappedFunc::faceAvgNormDerivFromCellAvgCS(faceAvgNormDerivFab,
                                                       linearFab,
                                                       0, // input comp
                                                       0, // output comp
                                                       npDomain,
                                                       faceBox,
                                                       dx,
                                                       dir,
                                                       false,  // second-order
                                                       false); // bounded domain
        int lstat = 0;
        const Real exactDeriv = exactDerivVec[dir];
        MD_BOXLOOP(faceBox, i)
          {
            const Real numericalDeriv = faceAvgNormDerivFab[MD_IX(i, 0)];
            if (Misc::compare1(numericalDeriv, exactDeriv, prec))
              {
                ++lstat;
              }
          }
        if (lstat) ++status;
      }
  }
  // fourth-order bounded-domain computation
  {
    CH_TIME("test::fourthOrderBoundaryScheme");
    for (int dir = 0; dir != SpaceDim; ++dir)
      {
        FArrayBox& faceAvgNormDerivFab = faceAvgNormDerivFxb[dir];
        faceAvgNormDerivFab.setVal(0.);
        Box faceBox = surroundingNodes(box, dir);
        PatchMappedFunc::faceAvgNormDerivFromCellAvgCS(faceAvgNormDerivFab,
                                                       linearFab,
                                                       0, // input comp
                                                       0, // output comp
                                                       npDomain,
                                                       faceBox,
                                                       dx,
                                                       dir,
                                                       true,   // fourth-order
                                                       false); // bounded domain
        int lstat = 0;
        const Real exactDeriv = exactDerivVec[dir];
        MD_BOXLOOP(faceBox, i)
          {
            const Real numericalDeriv = faceAvgNormDerivFab[MD_IX(i, 0)];
            if (Misc::compare1(numericalDeriv, exactDeriv, prec))
              {
                ++lstat;
              }
          }
        if (lstat) ++status;
      }
  }

/*--------------------------------------------------------------------*
 * Test normal derivative of sine-wave data field
 * -- testing order of accuracy
 *--------------------------------------------------------------------*/

  constexpr int numGrids = 4; // number of grids total (refinement is by 2)
  std::vector<Real> thisNorm(SpaceDim*3*4, 0.);
  std::vector<Real> prevNorm(SpaceDim*3*4, 0.);
  for (int grid = 0; grid != numGrids; ++grid)
    {
      const int m = std::pow(2, 4 + grid); // starts at 16, then 32, etc.
      const Box boxS(IntVect_unit*(-m/2), IntVect_unit*((m/2) - 1));
      const Box boxS1 = grow(boxS, 1);
      const Box boxS2 = grow(boxS, 2);

      FArrayBox sineFab(boxS, 1);
      FArrayBox sineFab1(boxS1, 1);
      FArrayBox sineFab2(boxS2, 1);
      const ProblemDomain pDomainS(boxS, periodic);
      const ProblemDomain npDomainS(boxS, notPeriodic);

      const RealVect dxS = RealVect_unit/m;
      const Real vol = dxS.product();
      const RealVect omega(D_DECL(1., 1., 1.)); // wavenumber in each direction
      const Real c_0 = 1./(std::pow(2.*PI, SpaceDim)*(omega.product())*vol);
      const RealVect c_1 = 2.*PI*omega;
      MD_BOXLOOP(boxS, i)
        {
          const RealVect loX = dxS*RealVect{MD_EXPANDIX((Real)i)};
          const RealVect hiX =
            dxS*(RealVect{MD_EXPANDIX((Real)i)} + RealVect::Unit);
          sineFab[MD_IX(i, 0)] = c_0*(
            D_TERM((std::cos(c_1[0]*loX[0]) - std::cos(c_1[0]*hiX[0])),
                  *(std::cos(c_1[1]*loX[1]) - std::cos(c_1[1]*hiX[1])),
                  *(std::cos(c_1[2]*loX[2]) - std::cos(c_1[2]*hiX[2]))));
        }
      MD_BOXLOOP(boxS1, i)
        {
          const RealVect loX = dxS*RealVect{MD_EXPANDIX((Real)i)};
          const RealVect hiX =
            dxS*(RealVect{MD_EXPANDIX((Real)i)} + RealVect::Unit);
          sineFab1[MD_IX(i, 0)] = c_0*(
            D_TERM((std::cos(c_1[0]*loX[0]) - std::cos(c_1[0]*hiX[0])),
                  *(std::cos(c_1[1]*loX[1]) - std::cos(c_1[1]*hiX[1])),
                  *(std::cos(c_1[2]*loX[2]) - std::cos(c_1[2]*hiX[2]))));
        }
      MD_BOXLOOP(boxS2, i)
        {
          const RealVect loX = dxS*RealVect{MD_EXPANDIX((Real)i)};
          const RealVect hiX =
            dxS*(RealVect{MD_EXPANDIX((Real)i)} + RealVect::Unit);
          sineFab2[MD_IX(i, 0)] = c_0*(
            D_TERM((std::cos(c_1[0]*loX[0]) - std::cos(c_1[0]*hiX[0])),
                  *(std::cos(c_1[1]*loX[1]) - std::cos(c_1[1]*hiX[1])),
                  *(std::cos(c_1[2]*loX[2]) - std::cos(c_1[2]*hiX[2]))));
        }

      // now test the derivatives and compare them to the exact values
      FLUXBOXSTACKTEMP(faceAvgNormDerivFxbS, boxS, 1);
      FLUXBOXSTACKTEMP(exactDerivFxb, boxS, 1);
      //  f_x = 2*PI*w_x*cos(2*PI*w_x*x)*sin(2*PI*w_y*y)*sin(2*PI*w_z*z)
      //  f_y = 2*PI*w_y*sin(2*PI*w_x*x)*cos(2*PI*w_y*y)*sin(2*PI*w_z*z)
      //  f_z = 2*PI*w_z*sin(2*PI*w_x*x)*sin(2*PI*w_y*y)*cos(2*PI*w_z*z)
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          Box faceBoxS = surroundingNodes(boxS, dir);
          FArrayBox& exactDerivFab = exactDerivFxb[dir];
          const IntVect dirVect(
            D_DECL(dir, ((dir + 1) % SpaceDim), ((dir + 2) % SpaceDim)));
          MD_BOXLOOP(faceBoxS, i)
            {
              RealVect xC = RealVect{MD_EXPANDIX((Real)i)};
              RealVect centerOffset(0.5*RealVect::Unit);
              centerOffset[dir] = 0.;
              xC += centerOffset;
              xC *= dxS;
              RealVect xL = RealVect{MD_EXPANDIX((Real)i)};
              xL *= dxS;
              RealVect xH = RealVect{MD_EXPANDIX((Real)i)};
              RealVect hiOffset(RealVect::Unit);
              hiOffset[dir] = 0.;
              xH += hiOffset;
              xH *= dxS;
              Real c_3 = D_TERM(
                2.*PI*omega[dirVect[0]],
                *(1./(2.*PI))*(1./(dxS[dirVect[1]]*omega[dirVect[1]])),
                *(1./(2.*PI))*(1./(dxS[dirVect[2]]*omega[dirVect[2]])));
              exactDerivFab[MD_IX(i, 0)] = c_3*(
                D_TERM(std::cos(2*PI*omega[dirVect[0]]*xC[dirVect[0]]),
                     *(std::cos(2*PI*omega[dirVect[1]]*xL[dirVect[1]])
                     - std::cos(2*PI*omega[dirVect[1]]*xH[dirVect[1]])),
                     *(std::cos(2*PI*omega[dirVect[2]]*xL[dirVect[2]])
                     - std::cos(2*PI*omega[dirVect[2]]*xH[dirVect[2]]))));
            }
        }
      // second-order interior computation
      {
        CH_TIME("test::secondOrderInteriorScheme");
        for (int dir = 0; dir != SpaceDim; ++dir)
          {
            FArrayBox& faceAvgNormDerivFab = faceAvgNormDerivFxbS[dir];
            faceAvgNormDerivFab.setVal(0.);
            Box faceBoxS = surroundingNodes(boxS, dir);
            PatchMappedFunc::faceAvgNormDerivFromCellAvgCS(faceAvgNormDerivFab,
                                                           sineFab1,
                                                           0, // input comp
                                                           0, // output comp
                                                           pDomainS,
                                                           faceBoxS,
                                                           dxS,
                                                           dir,
                                                           false, // 2nd-order
                                                           true); // interior
            int lstat = 0;
            FArrayBox& exactDerivFab = exactDerivFxb[dir];
            Real errorNormMax = 0.;
            Real errorNormL1 = 0.;
            Real errorNormL2 = 0.;
            MD_BOXLOOP(faceBoxS, i)
              {
                const Real numericalDeriv = faceAvgNormDerivFab[MD_IX(i, 0)];
                const Real exactDeriv = exactDerivFab[MD_IX(i, 0)];
                const Real error = numericalDeriv - exactDeriv;
                errorNormMax = std::fmax(errorNormMax, std::fabs(error));
                errorNormL1 += std::fabs(error);
                errorNormL2 += error*error;
              }
            const IntVect faceBoxSize = faceBoxS.size();
            const int numFaces = faceBoxSize.product();
            thisNorm[0 + dir*3] = errorNormL1/numFaces;
            thisNorm[1 + dir*3] = std::sqrt((1./numFaces)*errorNormL2);
            thisNorm[2 + dir*3] = errorNormMax;
            if (grid != 0)
              {
                // check order of accuracy here
                const Real slope1 = -(1./std::log(2.))*std::log(
                  thisNorm[0+dir*3]/prevNorm[0+dir*3]);
                const Real slope2 = -(1./std::log(2.))*std::log(
                  thisNorm[1+dir*3]/prevNorm[1+dir*3]);
                const Real slope3 = -(1./std::log(2.))*std::log(
                  thisNorm[2+dir*3]/prevNorm[2+dir*3]);
                if (slope1 < 1.85) ++lstat;
                if (slope2 < 1.85) ++lstat;
                if (slope3 < 1.85) ++lstat;
                if (lstat) ++status;
              }
            prevNorm[0+dir*3] = thisNorm[0+dir*3];
            prevNorm[1+dir*3] = thisNorm[1+dir*3];
            prevNorm[2+dir*3] = thisNorm[2+dir*3];
          }
      }
      // fourth-order interior computation
      {
        CH_TIME("test::fourthOrderInteriorScheme");
        for (int dir = 0; dir != SpaceDim; ++dir)
          {
            FArrayBox& faceAvgNormDerivFab = faceAvgNormDerivFxbS[dir];
            faceAvgNormDerivFab.setVal(0.);
            Box faceBoxS = surroundingNodes(boxS, dir);
            PatchMappedFunc::faceAvgNormDerivFromCellAvgCS(faceAvgNormDerivFab,
                                                           sineFab2,
                                                           0, // input comp
                                                           0, // output comp
                                                           pDomainS,
                                                           faceBoxS,
                                                           dxS,
                                                           dir,
                                                           true, // 4th-order
                                                           true); // interior
            int lstat = 0;
            FArrayBox& exactDerivFab = exactDerivFxb[dir];
            Real errorNormMax = 0.;
            Real errorNormL1 = 0.;
            Real errorNormL2 = 0.;
            MD_BOXLOOP(faceBoxS, i)
              {
                const Real numericalDeriv = faceAvgNormDerivFab[MD_IX(i, 0)];
                const Real exactDeriv = exactDerivFab[MD_IX(i, 0)];
                const Real error = numericalDeriv - exactDeriv;
                errorNormMax = std::fmax(errorNormMax, std::fabs(error));
                errorNormL1 += std::fabs(error);
                errorNormL2 += error*error;
              }
            const IntVect faceBoxSize = faceBoxS.size();
            const int numFaces = faceBoxSize.product();
            thisNorm[0+dir*3+SpaceDim*3] = errorNormL1/numFaces;
            thisNorm[1+dir*3+SpaceDim*3] = std::sqrt((1./numFaces)*errorNormL2);
            thisNorm[2+dir*3+SpaceDim*3] = errorNormMax;
            if (grid != 0)
              {
                // check order of accuracy here
                const Real slope1 = -(1./std::log(2.))*std::log(
                  thisNorm[0+dir*3+SpaceDim*3]/prevNorm[0+dir*3+SpaceDim*3]);
                const Real slope2 = -(1./std::log(2.))*std::log(
                  thisNorm[1+dir*3+SpaceDim*3]/prevNorm[1+dir*3+SpaceDim*3]);
                const Real slope3 = -(1./std::log(2.))*std::log(
                  thisNorm[2+dir*3+SpaceDim*3]/prevNorm[2+dir*3+SpaceDim*3]);
                if (slope1 < 3.85) ++lstat;
                if (slope2 < 3.85) ++lstat;
                if (slope3 < 3.85) ++lstat;
                if (lstat) ++status;
              }
            prevNorm[0+dir*3+SpaceDim*3] = thisNorm[0+dir*3+SpaceDim*3];
            prevNorm[1+dir*3+SpaceDim*3] = thisNorm[1+dir*3+SpaceDim*3];
            prevNorm[2+dir*3+SpaceDim*3] = thisNorm[2+dir*3+SpaceDim*3];
          }
      }
      // second-order bounded-domain computation
      {
        CH_TIME("test::secondOrderBoundaryScheme");
        for (int dir = 0; dir != SpaceDim; ++dir)
          {
            FArrayBox& faceAvgNormDerivFab = faceAvgNormDerivFxbS[dir];
            faceAvgNormDerivFab.setVal(0.);
            Box faceBoxS = surroundingNodes(boxS, dir);
            PatchMappedFunc::faceAvgNormDerivFromCellAvgCS(faceAvgNormDerivFab,
                                                           sineFab,
                                                           0, // input comp
                                                           0, // output comp
                                                           npDomainS,
                                                           faceBoxS,
                                                           dxS,
                                                           dir,
                                                           false,  // 2nd-order
                                                           false); // bounded
            int lstat = 0;
            FArrayBox& exactDerivFab = exactDerivFxb[dir];
            Real errorNormMax = 0.;
            Real errorNormL1 = 0.;
            Real errorNormL2 = 0.;
            MD_BOXLOOP(faceBoxS, i)
              {
                const Real numericalDeriv = faceAvgNormDerivFab[MD_IX(i, 0)];
                const Real exactDeriv = exactDerivFab[MD_IX(i, 0)];
                const Real error = numericalDeriv - exactDeriv;
                errorNormMax = std::fmax(errorNormMax, std::fabs(error));
                errorNormL1 += std::fabs(error);
                errorNormL2 += error*error;
              }
            const IntVect faceBoxSize = faceBoxS.size();
            const int numFaces = faceBoxSize.product();
            thisNorm[0+dir*3+SpaceDim*6] = errorNormL1/numFaces;
            thisNorm[1+dir*3+SpaceDim*6] = std::sqrt((1./numFaces)*errorNormL2);
            thisNorm[2+dir*3+SpaceDim*6] = errorNormMax;
            if (grid != 0)
              {
                // check order of accuracy here
                const Real slope1 = -(1./std::log(2.))*std::log(
                  thisNorm[0+dir*3+SpaceDim*6]/prevNorm[0+dir*3+SpaceDim*6]);
                const Real slope2 = -(1./std::log(2.))*std::log(
                  thisNorm[1+dir*3+SpaceDim*6]/prevNorm[1+dir*3+SpaceDim*6]);
                const Real slope3 = -(1./std::log(2.))*std::log(
                  thisNorm[2+dir*3+SpaceDim*6]/prevNorm[2+dir*3+SpaceDim*6]);
                if (slope1 < 1.85) ++lstat;
                if (slope2 < 1.85) ++lstat;
                if (slope3 < 1.85) ++lstat;
                if (lstat) ++status;
              }
            prevNorm[0+dir*3+SpaceDim*6] = thisNorm[0+dir*3+SpaceDim*6];
            prevNorm[1+dir*3+SpaceDim*6] = thisNorm[1+dir*3+SpaceDim*6];
            prevNorm[2+dir*3+SpaceDim*6] = thisNorm[2+dir*3+SpaceDim*6];
          }
      }
      // fourth-order bounded-domain computation
      {
        CH_TIME("test::fourthOrderBoundaryScheme");
        for (int dir = 0; dir != SpaceDim; ++dir)
          {
            FArrayBox& faceAvgNormDerivFab = faceAvgNormDerivFxbS[dir];
            faceAvgNormDerivFab.setVal(0.);
            Box faceBoxS = surroundingNodes(boxS, dir);
            PatchMappedFunc::faceAvgNormDerivFromCellAvgCS(faceAvgNormDerivFab,
                                                           sineFab,
                                                           0, // input comp
                                                           0, // output comp
                                                           npDomainS,
                                                           faceBoxS,
                                                           dxS,
                                                           dir,
                                                           true,   // 4th-order
                                                           false); // bounded
            int lstat = 0;
            FArrayBox& exactDerivFab = exactDerivFxb[dir];
            Real errorNormMax = 0.;
            Real errorNormL1 = 0.;
            Real errorNormL2 = 0.;
            MD_BOXLOOP(faceBoxS, i)
              {
                const Real numericalDeriv = faceAvgNormDerivFab[MD_IX(i, 0)];
                const Real exactDeriv = exactDerivFab[MD_IX(i, 0)];
                const Real error = numericalDeriv - exactDeriv;
                errorNormMax = std::fmax(errorNormMax, std::fabs(error));
                errorNormL1 += std::fabs(error);
                errorNormL2 += error*error;
              }
            const IntVect faceBoxSize = faceBoxS.size();
            const int numFaces = faceBoxSize.product();
            thisNorm[0+dir*3+SpaceDim*9] = errorNormL1/numFaces;
            thisNorm[1+dir*3+SpaceDim*9] = std::sqrt((1./numFaces)*errorNormL2);
            thisNorm[2+dir*3+SpaceDim*9] = errorNormMax;
            if (grid != 0)
              {
                // check order of accuracy here
                const Real slope1 = -(1./std::log(2.))*std::log(
                  thisNorm[0+dir*3+SpaceDim*9]/prevNorm[0+dir*3+SpaceDim*9]);
                const Real slope2 = -(1./std::log(2.))*std::log(
                  thisNorm[1+dir*3+SpaceDim*9]/prevNorm[1+dir*3+SpaceDim*9]);
                const Real slope3 = -(1./std::log(2.))*std::log(
                  thisNorm[2+dir*3+SpaceDim*9]/prevNorm[2+dir*3+SpaceDim*9]);
                if (slope1 < 3.85) ++lstat;
                if (slope2 < 3.85) ++lstat;
                if (slope3 < 3.85) ++lstat;
                if (lstat) ++status;
              }
            prevNorm[0+dir*3+SpaceDim*9] = thisNorm[0+dir*3+SpaceDim*9];
            prevNorm[1+dir*3+SpaceDim*9] = thisNorm[1+dir*3+SpaceDim*9];
            prevNorm[2+dir*3+SpaceDim*9] = thisNorm[2+dir*3+SpaceDim*9];
          }
      }

      // now take the norms of the differences and test order of accuracy
    }

  return status;
}
