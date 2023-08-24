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
test_PatchMappedFunc_CellAvgDerivFromFaceAvg();

/// Global variables for handling output:
static const char *pgmname = "test_PatchMappedFunc_CellAvgDerivFromFaceAvg";
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

  int ret = test_PatchMappedFunc_CellAvgDerivFromFaceAvg();
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
/// Routine test_PatchMappedFunc_CellAvgDerivFromFaceAvg
/*
 ******************************************************************************/

int
test_PatchMappedFunc_CellAvgDerivFromFaceAvg()
{
  CH_TIMERS("test_PatchMappedFunc_CellAvgDerivFromFaceAvg");
  const int prec = std::numeric_limits<Real>::digits10 - 1;
  int status = 0;
//--Initial set up
//  We need to test the exactness of this function

//  The first IC is a product of sine waves. This is set to
//  f(x,y,z) = sin(2*PI*w_x*x)*sin(2*PI*w_y*y)*sin(2*PI*w_z*z)
//  This way, the derivatives are
//  f_x = 2*PI*w_x*cos(2*PI*w_x*x)*sin(2*PI*w_y*y)*sin(2*PI*w_z*z)
//  f_y = 2*PI*w_y*sin(2*PI*w_x*x)*cos(2*PI*w_y*y)*sin(2*PI*w_z*z)
//  f_z = 2*PI*w_z*sin(2*PI*w_x*x)*sin(2*PI*w_y*y)*cos(2*PI*w_z*z)

/*--------------------------------------------------------------------*
 * Test normal derivative of sine-wave data field
 * -- testing order of accuracy
 *--------------------------------------------------------------------*/

  const int m = 32;
  const Box box(IntVect_unit*(-m/2), IntVect_unit*((m/2) - 1));
  const RealVect dx = RealVect_unit/m;
  const Real vol = dx.product();
  const RealVect omega(D_DECL(1., 1., 1.)); // wavenumber in each direction

  FLUXBOXSTACKTEMP(faceAvgFxb, box, 1);
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      Box faceBox = surroundingNodes(box, dir);
      FArrayBox exactCellGradFab(box, 1);
      FArrayBox numericalCellGradFab(box, 1);
      FArrayBox& faceAvgFab = faceAvgFxb[dir];

      const IntVect dirVect(
        D_DECL(dir, ((dir + 1) % SpaceDim), ((dir + 2) % SpaceDim)));

      // exact cell-averaged derivative
      MD_BOXLOOP(box, i)
        {
          RealVect xL = dx*RealVect{MD_EXPANDIX((Real)i)};
          RealVect xH = dx*(RealVect{MD_EXPANDIX((Real)i)} + RealVect::Unit);
          Real c_3 = D_TERM(
            (1./dx[dir]),
           *(1./(2.*PI))*(1./(dx[dirVect[1]]*omega[dirVect[1]])),
           *(1./(2.*PI))*(1./(dx[dirVect[2]]*omega[dirVect[2]])));
          exactCellGradFab[MD_IX(i, 0)] = c_3*(
            D_TERM((std::sin(2*PI*omega[dirVect[0]]*xH[dirVect[0]])
                  - std::sin(2*PI*omega[dirVect[0]]*xL[dirVect[0]])),
                  *(std::cos(2*PI*omega[dirVect[1]]*xL[dirVect[1]])
                  - std::cos(2*PI*omega[dirVect[1]]*xH[dirVect[1]])),
                  *(std::cos(2*PI*omega[dirVect[2]]*xL[dirVect[2]])
                  - std::cos(2*PI*omega[dirVect[2]]*xH[dirVect[2]]))));
        }
      // exact face-averaged values
      MD_BOXLOOP(faceBox, i)
        {
          RealVect xC = RealVect{MD_EXPANDIX((Real)i)};
          RealVect centerOffset(0.5*RealVect::Unit);
          centerOffset[dir] = 0.;
          xC += centerOffset;
          xC *= dx;
          RealVect xL = dx*RealVect{MD_EXPANDIX((Real)i)};
          RealVect xH = RealVect{MD_EXPANDIX((Real)i)};
          RealVect hiOffset(RealVect::Unit);
          hiOffset[dir] = 0.;
          xH += hiOffset;
          xH *= dx;
          Real c_3 =
            D_TERM(1.,
                 *(1./(2.*PI))*(1./(dx[dirVect[1]]*omega[dirVect[1]])),
                 *(1./(2.*PI))*(1./(dx[dirVect[2]]*omega[dirVect[2]])));
          faceAvgFab[MD_IX(i, 0)] = c_3*(
            D_TERM(std::sin(2*PI*omega[dirVect[0]]*xC[dirVect[0]]),
                 *(std::cos(2*PI*omega[dirVect[1]]*xL[dirVect[1]])
                 - std::cos(2*PI*omega[dirVect[1]]*xH[dirVect[1]])),
                 *(std::cos(2*PI*omega[dirVect[2]]*xL[dirVect[2]])
                 - std::cos(2*PI*omega[dirVect[2]]*xH[dirVect[2]]))));
        }
      // test the exactness of the scheme
      {
        CH_TIME("test::exactScheme");
        numericalCellGradFab.setVal(0.);
        PatchMappedFunc::cellAvgDerivFromFaceAvgCS(numericalCellGradFab,
                                                   faceAvgFab,
                                                   0, // input comp
                                                   0, // output comp
                                                   box,
                                                   dx,
                                                   dir); // interior
        int lstat = 0;
        MD_BOXLOOP(box, i)
          {
            const Real numericalVal = numericalCellGradFab[MD_IX(i, 0)];
            const Real exactVal = exactCellGradFab[MD_IX(i, 0)];
            if (Misc::compare1(numericalVal, exactVal, prec)) ++lstat;
          }
        if (lstat) ++status;
      }
    }

  return status;
}
