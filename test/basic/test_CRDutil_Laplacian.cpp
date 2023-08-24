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
test_CRDutil_Laplacian();

/// Global variables for handling output:
static const char *pgmname = "test_CRDutil_Laplacian";
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

  int ret = test_CRDutil_Laplacian();
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
test_CRDutil_Laplacian()
{
  CH_TIMERS("test_CRDutil_Laplacian");
  const int prec = std::numeric_limits<Real>::digits10 - 1;
  int status = 0;

//--Initial conditions

  constexpr int n = 16;
  const Box box(IntVect_unit*(1-n/2), IntVect_unit*(n/2));
  const Box box1 = grow(box, 1);

  FArrayBox fab(box1, 2);
  bool periodic[3]    = {true, true, true};
  bool notPeriodic[3] = {false, false, false};
  ProblemDomain pDomain(box, periodic);
  ProblemDomain npDomain(box, notPeriodic);

  RealVect dx = RealVect_unit/n;
  MD_BOXLOOP(box1, i)
    {
      RealVect x = dx*RealVect{MD_EXPANDIX((Real)i)} + 0.5;
      fab[MD_IX(i, 0)] = stc::sum(RealVect{2, 3, 4}*x*x);
      fab[MD_IX(i, 1)] = D_TERM(  std::sin(2*Pi*x[0]),
                                + std::sin(2*Pi*x[1]),
                                + std::sin(2*Pi*x[2]));
    }

/*--------------------------------------------------------------------*
 * Test Laplacian using Godunov Utilities (cache warm up)
 *--------------------------------------------------------------------*/

  FABSTACKTEMP(d2fabChF, box, 2);
  FABSTACKTEMP(d2fabChC, box, 2);
  FABSTACKTEMP(tfab, box, 2);
  d2fabChF.setVal(0.);
  d2fabChC.setVal(0.);
  {
    constexpr Real scale = 1.;
    for (int dir = 0; dir != SpaceDim; ++dir)
      {
        Box loBox, hiBox, centerBox, entireBox;
        int hasLo, hasHi;
        // Generate the domain boundary boxes, loBox and hiBox, if there are
        // domain boundaries there
        Box dirBox(box);
        dirBox.grow(dir, 1);
        loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                   dirBox, pDomain, dir);
        int nComp = 2;
        FORT_GETSECONDDIFF(CHF_FRA(tfab),
                           CHF_CONST_FRA(fab),
                           CHF_CONST_INT(nComp),
                           CHF_CONST_INT(dir),
                           CHF_BOX(loBox),
                           CHF_CONST_INT(hasLo),
                           CHF_BOX(hiBox),
                           CHF_CONST_INT(hasHi),
                           CHF_BOX(centerBox));
        d2fabChF.plus(tfab, scale);
        d2fabChC.plus(tfab, scale);
      }
  }

/*--------------------------------------------------------------------*
 * Test Laplacian using Godunov Utilities (periodic)
 *--------------------------------------------------------------------*/

  CH_TIMER("Godunov:Laplacian", t1);
  CH_START(t1);
  for (int n = 1000; n--;)
    {
      d2fabChF.setVal(0.);
      {
        constexpr Real scale = 1.;
        for (int dir = 0; dir != SpaceDim; ++dir)
          {
            Box loBox, hiBox, centerBox, entireBox;
            int hasLo, hasHi;
            // Generate the domain boundary boxes, loBox and hiBox, if there are
            // domain boundaries here
            Box box1dir(box);
            box1dir.grow(dir, 1);
            loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                       box1dir, pDomain, dir);
            int nComp = 2;
            FORT_GETSECONDDIFF(CHF_FRA(tfab),
                               CHF_CONST_FRA(fab),
                               CHF_CONST_INT(nComp),
                               CHF_CONST_INT(dir),
                               CHF_BOX(loBox),
                               CHF_CONST_INT(hasLo),
                               CHF_BOX(hiBox),
                               CHF_CONST_INT(hasHi),
                               CHF_BOX(centerBox));
            d2fabChF.plus(tfab, scale);
          }
      }
    }
  CH_STOP(t1);
  {
    const Real valid = stc::sum(2*RealVect{2, 3, 4});
    int lstat = 0;
    MD_BOXLOOP(box, i)
      {
        if (Misc::compare1(d2fabChF[MD_IX(i, 0)]/(dx[0]*dx[0]), valid, prec))
          {
            ++lstat;
          }
      }
    if (lstat) ++status;
  }

/*--------------------------------------------------------------------*
 * Test Laplacian using CRD utilities (periodic)
 *--------------------------------------------------------------------*/

  CH_TIMER("CRDutil:Laplacian", t2);
  CH_START(t2);
  for (int n = 1000; n--;)
    {
      CRDutil::Laplacian(d2fabChC,
                         fab,
                         box,
                         pDomain,
                         fab.interval());
    }
  CH_STOP(t2);
  {
    const Real valid = stc::sum(2*RealVect{2, 3, 4});
    int lstat = 0;
    MD_BOXLOOP(box, i)
      {
        if (Misc::compare1(d2fabChC[MD_IX(i, 0)]/(dx[0]*dx[0]), valid, prec))
          {
            ++lstat;
          }
      }
    if (lstat) ++status;
  }

//--Compare both approaches

  {
    int lstat = 0;
    MD_BOXLOOP(box, i)
      {
        if (Misc::compare1(d2fabChF[MD_IX(i, 0)], d2fabChC[MD_IX(i, 0)], prec))
          ++lstat;
        if (Misc::compare1(d2fabChF[MD_IX(i, 1)], d2fabChC[MD_IX(i, 1)], prec))
          ++lstat;
      }
    if (lstat) ++status;
  }

/*--------------------------------------------------------------------*
 * Test Laplacian using Godunov Utilities (not-periodic)
 *--------------------------------------------------------------------*/

  MD_BOXLOOP(box1, i)
    {
      if (!box.contains(IntVect{MD_EXPANDIX(i)}))
        {
          fab[MD_IX(i, 0)] = 0.;
          fab[MD_IX(i, 1)] = 0.;
        }
    }

  CH_TIMER("Godunov:Laplacian-np", t3);
  CH_START(t3);
  for (int n = 1000; n--;)
    {
      d2fabChF.setVal(0.);
      {
        constexpr Real scale = 1.;
        for (int dir = 0; dir != SpaceDim; ++dir)
          {
            Box loBox, hiBox, centerBox, entireBox;
            int hasLo, hasHi;
            // Generate the domain boundary boxes, loBox and hiBox, if there are
            // domain boundaries here
            Box box1dir(box);
            box1dir.grow(dir, 1);
            loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                       box1dir, npDomain, dir);
            int nComp = 2;
            FORT_GETSECONDDIFF(CHF_FRA(tfab),
                               CHF_CONST_FRA(fab),
                               CHF_CONST_INT(nComp),
                               CHF_CONST_INT(dir),
                               CHF_BOX(loBox),
                               CHF_CONST_INT(hasLo),
                               CHF_BOX(hiBox),
                               CHF_CONST_INT(hasHi),
                               CHF_BOX(centerBox));
            d2fabChF.plus(tfab, scale);
          }
      }
    }
  CH_STOP(t3);
  {
    const Real valid = stc::sum(2*RealVect{2, 3, 4});
    int lstat = 0;
    MD_BOXLOOP(box, i)
      {
        if (Misc::compare1(d2fabChF[MD_IX(i, 0)]/(dx[0]*dx[0]), valid, prec))
          {
            ++lstat;
          }
      }
    if (lstat) ++status;
  }

/*--------------------------------------------------------------------*
 * Test Laplacian using CRD utilities (not-periodic)
 *--------------------------------------------------------------------*/

  CH_TIMER("CRDutil:Laplacian-np", t4);
  CH_START(t4);
  for (int n = 1000; n--;)
    {
      CRDutil::Laplacian(d2fabChC,
                         fab,
                         box,
                         npDomain,
                         fab.interval());
    }
  CH_STOP(t4);
  {
    const Real valid = stc::sum(2*RealVect{2, 3, 4});
    int lstat = 0;
    MD_BOXLOOP(box, i)
      {
        if (Misc::compare1(d2fabChC[MD_IX(i, 0)]/(dx[0]*dx[0]), valid, prec))
          {
            ++lstat;
          }
      }
    if (lstat) ++status;
  }

//--Compare both approaches

  {
    int lstat = 0;
    MD_BOXLOOP(box, i)
      {
        if (Misc::compare1(d2fabChF[MD_IX(i, 0)], d2fabChC[MD_IX(i, 0)], prec))
          ++lstat;
        if (Misc::compare1(d2fabChF[MD_IX(i, 1)], d2fabChC[MD_IX(i, 1)], prec))
          ++lstat;
      }
    if (lstat) ++status;
  }

  return status;
}
