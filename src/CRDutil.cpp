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
 * \file CRDutil.cpp
 *
 * \brief Utilities for Chord
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "ProblemDomain.H"
#include "FArrayBox.H"
#include "BaseFabMacros.H"
#include "LoHiCenter.H"
#include "LoHiCenterApplyOp.H"
#include "Misc.H"
#include "GodunovUtilitiesF_F.H"
#include "UnitNormalsF_F.H"
#include "RootSolver.H"

//----- Internal -----//

#include "CRDparam.H"
#include "CRDutil.H"
#include "CRDPhysics.H"
#include "CNSIBC.H"
#include "DataTemp.H"
#include "PatchMappedFunc.H"

/*--------------------------------------------------------------------*/
//
/**
 *--------------------------------------------------------------------*/

/*==============================================================================
 * Definition of static variables
 *============================================================================*/

inline static Real
d1Func(const Real W_ip1,
       const Real W_i,
       Real       d2W_ip1h) noexcept
{
  // Use a quadratic interpolation to get derivatives on the
  // cells left and right of the face
  const Real d1W_L = W_ip1 - W_i - 0.5*d2W_ip1h;
  const Real d1W_R = W_ip1 - W_i + 0.5*d2W_ip1h;
  if (d1W_L*d1W_R < 0.)
    {
      d2W_ip1h = 0.;
    }
  else if (std::abs(d1W_L) > 3.*std::abs(d1W_R))
    {
      // Force d1W = 0 on right face
      d2W_ip1h = W_i - W_ip1;
    }
  else if (std::abs(d1W_R) > 3.*std::abs(d1W_L))
    {
      // Force d1W = 0 on left face
      d2W_ip1h = W_ip1 - W_i;
    }
  return 0.5*(W_i + W_ip1) - d2W_ip1h/6.;
}

void
CRDutil::PPMFaceValues(FArrayBox&           a_WFace,
                       const FArrayBox&     a_W,
                       const FArrayBox&     a_smoothTest,
                       const int&           a_numComp,
                       const bool&          a_useLimiting,
                       const int&           a_dir,
                       const Box&           a_faceBox,
                       const ProblemDomain& a_domain)
{
  CH_TIME("CRDutil::PPMFaceValues");
  // Tolerance for smoothness tests
  Box box1cells = grow(a_faceBox, BASISV(a_dir));
  box1cells.enclosedCells();

  const Real hiTol = 1e100; // high tolerance to check if data is valid
  const Real loTol = -1e100; // low tolerance to check if data is valid

  Box nearLoFaces, midLoFaces, farLoFaces;
  Box nearHiFaces, midHiFaces, farHiFaces;
  Box centerFaces, midCenterFaces, innerCenterFaces, entireFaces;
  int hasLoFaces, hasHiFaces;
  loHiCenterFace6(farLoFaces, midLoFaces, nearLoFaces, hasLoFaces,
                  farHiFaces, midHiFaces, nearHiFaces, hasHiFaces,
                  centerFaces, midCenterFaces, innerCenterFaces,
                  entireFaces, box1cells, a_domain, a_dir);

  CH_assert(a_WFace.box().contains(entireFaces));
  Box testBox(innerCenterFaces);
  testBox.enclosedCells();
  testBox.grow(a_dir, 3);
  testBox &= a_domain;
  CH_assert(a_W.box().contains(testBox));

  MD_ARRAY_RESTRICT(arrWFace, a_WFace);
  MD_ARRAY_RESTRICT(arrW, a_W);

  const int MD_ID(o, a_dir);

  // Shift boxes so that inner center index corresponds to i+1/2
  /*
    Let's say the inner center face has index 4.  When we write to i+1/2, we
    want to be writing to 4.  Set we set i = 3 so that the left face is i-1/2
    (also 3) and the right face is i+1/2 (index 4).  This way, we can always
    work with stencils centered around i+1/2.

    For nearLo and nearHi, the stencil is still centered around i+1/2
  */
  nearLoFaces.shift(a_dir, -1);  // Two cells from boundary
  innerCenterFaces.shift(a_dir, -1);
  nearHiFaces.shift(a_dir, -1);  // Two cells from boundary

  for (int c = 0; c != a_numComp; ++c)
    {
      if (hasLoFaces)
        {
          // Loop over faces near lower boundary
          MD_BOXLOOP(nearLoFaces, i)  // Cell index i or face index i-1/2
            {
              // Cell values
              const Real W_im1 = arrW[MD_OFFSETIX(i,-,1*o,c)];
              const Real W_i   = arrW[MD_IX(i,c)];
              const Real W_ip1 = arrW[MD_OFFSETIX(i,+,1*o,c)];
              const Real W_ip2 = arrW[MD_OFFSETIX(i,+,2*o,c)];
              const Real W_ip3 = arrW[MD_OFFSETIX(i,+,3*o,c)];

              CH_assert(W_im1 < hiTol && W_im1 > loTol);
              CH_assert(W_i   < hiTol && W_i   > loTol);
              CH_assert(W_ip1 < hiTol && W_ip1 > loTol);
              CH_assert(W_ip2 < hiTol && W_ip2 > loTol);
              CH_assert(W_ip3 < hiTol && W_ip3 > loTol);

              // Linear extension to predict W_im2
              const Real W_im2 = 2.*W_im1 - W_i;

              // 2nd derivatives on the faces at i-1/2 and i+1/2
              // If unconstrained (4 cells near low boundary)
              const Real d2W_im3h =
                0.5*(5.*W_im1 - 13.*W_i + 11.*W_ip1 - 3.*W_ip2);
              const Real d2W_im1h =
                0.5*(3.*W_im1 -  7.*W_i +  5.*W_ip1 -    W_ip2);
              // If constrained by d1W_im3h = 0 (3 cells near low boundary)
              // const Real d2W_im3h = (-25.*W_im1 + 32.*W_i - 7.*W_ip1)/11.;
              // const Real d2W_im1h = (    -W_im1 -  4.*W_i + 5.*W_ip1)/11.;
              // Centered d2W
              const Real d2W_ip1h = 0.5*((W_ip2 - W_ip1) - (W_i   - W_im1));
              const Real d2W_ip3h = 0.5*((W_ip3 - W_ip2) - (W_ip1 - W_i));

              // Compute fifth-order face values.  All use two cells and 2nd
              // derivatives
#if 1
              Real W_im3h = 0.5*(3.*W_im1 - W_i) -
                (-9.*d2W_im3h - 32.*d2W_im1h + d2W_ip1h)/120.;
              Real W_im1h =
                0.5*(W_im1 + W_i) - (d2W_im3h + 18.*d2W_im1h + d2W_ip1h)/120.;
              Real W_ip1h =
                0.5*(W_i + W_ip1) - (d2W_im1h + 18.*d2W_ip1h + d2W_ip3h)/120.;
#else
              // Compute the second/third/fourth-order face value (full fourth-
              // order is commented out)
              Real W_im3h = 0.5*(W_im1 + W_i);
              // Real W_im3h = (25.*W_im1 - 23.*W_i + 13.*W_ip1  - 3.*W_ip2)/12.;
              Real W_im1h = (2.*W_im1 + 5.*W_i - W_ip1)/6.;
              // Real W_im1h = ( 3.*W_im1 + 13.*W_i -  5.*W_ip1  +    W_ip2)/12.;
              Real W_ip1h = (   -W_im1 + 7.*(W_i +     W_ip1) -    W_ip2)/12.;
#endif

              // // Not bad
              // const int sd2W_im3h = 1 - 2*(int)std::signbit(d2W_im3h);
              // const int sd2W_im1h = 1 - 2*(int)std::signbit(d2W_im1h);
              // const int sd2W_ip1h = 1 - 2*(int)std::signbit(d2W_ip1h);
              // if (std::abs(sd2W_im3h + sd2W_im1h + sd2W_ip1h) < 3)
              // if (c > CRDPhysics::pressureIndex())
                {
                  Real minW_im3h = std::min(W_im2, W_im1);
                  Real maxW_im3h = std::max(W_im2, W_im1);
                  W_im3h = std::min(maxW_im3h, std::max(minW_im3h, W_im3h));
                  Real minW_im1h = std::min(W_im1, W_i);
                  Real maxW_im1h = std::max(W_im1, W_i);
                  W_im1h = std::min(maxW_im1h, std::max(minW_im1h, W_im1h));
                }

              // // Face ip1h uses d2W_ip3h instead of d2W_im3h
              // const int sd2W_ip3h = 1 - 2*(int)std::signbit(d2W_ip3h);
              // if (std::abs(sd2W_im1h + sd2W_ip1h + sd2W_ip3h) < 3)
              // if (c > CRDPhysics::pressureIndex())
                {
                  Real minW = std::min(W_i, W_ip1);
                  Real maxW = std::max(W_i, W_ip1);
                  W_ip1h = std::min(maxW, std::max(minW, W_ip1h));
                }

              // // Drop to 3rd order if not bounded
              // if ((W_im3h - W_im2)*(W_im3h - W_im1) > 0.)  // Not bounded
              //   {
              //     W_im3h = (9.*W_im1 - 3.*W_i + 2.*d2W_im3h)/6.;
              //   }
              // if ((W_im1h - W_im1)*(W_im1h - W_i)   > 0.)  // Not bounded
              //   {
              //     W_im1h = 0.5*(W_im1 + W_i) - d2W_im1h/6.;
              //   }
              // if ((W_ip1h - W_i)  *(W_ip1h - W_ip1) > 0.)  // Not bounded
              //   {
              //     W_ip1h = 0.5*(W_i + W_ip1) - d2W_ip1h/6.;
              //   }

              arrWFace[MD_OFFSETIX(i,-,o,c)] = W_im3h;
              arrWFace[MD_IX(i,c)]           = W_im1h;
              arrWFace[MD_OFFSETIX(i,+,o,c)] = W_ip1h;
            }
        }
      // Loop over center faces that are well away from boundaries
      MD_BOXLOOP(innerCenterFaces, i)  // Cell index i or face index i-1/2
        {
          // Cell-averaged values
          const Real W_im2 = arrW[MD_OFFSETIX(i,-,2*o,c)];
          const Real W_im1 = arrW[MD_OFFSETIX(i,-,1*o,c)];
          const Real W_i   = arrW[MD_IX(i,c)];
          const Real W_ip1 = arrW[MD_OFFSETIX(i,+,1*o,c)];
          const Real W_ip2 = arrW[MD_OFFSETIX(i,+,2*o,c)];
          const Real W_ip3 = arrW[MD_OFFSETIX(i,+,3*o,c)];

          CH_assert(W_im2 < hiTol && W_im2 > loTol);
          CH_assert(W_im1 < hiTol && W_im1 > loTol);
          CH_assert(W_i   < hiTol && W_i   > loTol);
          CH_assert(W_ip1 < hiTol && W_ip1 > loTol);
          CH_assert(W_ip2 < hiTol && W_ip2 > loTol);
          CH_assert(W_ip3 < hiTol && W_ip3 > loTol);

          // First get the 2nd derivatives on the faces at i-1/2, i+1/2, and
          // i+3/2)
          const Real d2W_im1h = 0.5*((W_ip1 - W_i)   - (W_im1 - W_im2));
          const Real d2W_ip1h = 0.5*((W_ip2 - W_ip1) - (W_i   - W_im1));
          const Real d2W_ip3h = 0.5*((W_ip3 - W_ip2) - (W_ip1 - W_i));

          // Compute the fifth-order face value
#if 1
          Real W_ip1h =
            0.5*(W_i + W_ip1) - (d2W_im1h + 18.*d2W_ip1h + d2W_ip3h)/120.;
#else
          // Compute the fourth-order face value
          Real W_ip1h = (-W_im1 + 7.*(W_i + W_ip1) - W_ip2)/12.;
#endif

          // // Not bad
          // const int sd2W_im1h = 1 - 2*(int)std::signbit(d2W_im1h);
          // const int sd2W_ip1h = 1 - 2*(int)std::signbit(d2W_ip1h);
          // const int sd2W_ip3h = 1 - 2*(int)std::signbit(d2W_ip3h);
          // if (std::abs(sd2W_im1h + sd2W_ip1h + sd2W_ip3h) < 3)
          //   {
          //     Real d2W_max = std::max(
          //       std::abs(d2W_im1h),
          //       std::max(std::abs(d2W_ip1h), std::abs(d2W_ip3h)));
          //     if (d2W_max > std::abs(W_i))
          // if (c > CRDPhysics::pressureIndex())
            {
              Real minW = std::min(W_i, W_ip1);
              Real maxW = std::max(W_i, W_ip1);
              W_ip1h = std::min(maxW, std::max(minW, W_ip1h));
            }
          // }

          // // Drop to 3rd order if not bounded
          // if ((W_ip1h - W_ip1)*(W_ip1h - W_i) > 0.)  // Not bounded
          //   {
          //     W_ip1h = 0.5*(W_i + W_ip1) - d2W_ip1h/6.;
          //   }

          CH_assert(a_WFace.box().contains(IntVect(D_DECL(i0+o0,
                                                          i1+o1,
                                                          i2+o2))));
          arrWFace[MD_OFFSETIX(i,+,o, c)] = W_ip1h;
        }
      if (hasHiFaces)
        {
          // Loop over faces near upper boundary
          MD_BOXLOOP(nearHiFaces, i)  // Cell index i or face index i-1/2
            {
              // Cell-averaged values
              const Real W_ip2 = arrW[MD_OFFSETIX(i,+,2*o,c)];
              const Real W_ip1 = arrW[MD_OFFSETIX(i,+,1*o,c)];
              const Real W_i   = arrW[MD_IX(i,c)];
              const Real W_im1 = arrW[MD_OFFSETIX(i,-,1*o,c)];
              const Real W_im2 = arrW[MD_OFFSETIX(i,-,2*o,c)];

              CH_assert(W_im2 < hiTol && W_im2 > loTol);
              CH_assert(W_im1 < hiTol && W_im1 > loTol);
              CH_assert(W_i   < hiTol && W_i   > loTol);
              CH_assert(W_ip1 < hiTol && W_ip1 > loTol);
              CH_assert(W_ip2 < hiTol && W_ip2 > loTol);

              // Linear extension to predict W_ip3
              const Real W_ip3 = 2.*W_ip2 - W_ip1;

              // 2nd derivatives on the faces at i+3/2 and i+5/2
              // If unconstrained (4 cells near high boundary)
              const Real d2W_ip5h =
                0.5*(5.*W_ip2 - 13.*W_ip1 + 11.*W_i - 3.*W_im1);
              const Real d2W_ip3h =
                0.5*(3.*W_ip2 -  7.*W_ip1 +  5.*W_i -    W_im1);
              // If constrained by d1W_ip5h = 0 (3 cells near high boundary)
              // const Real d2W_ip5h = ( -25.*W_ip2 + 32.*W_ip1 - 7.*W_i)/11.;
              // const Real d2W_ip3h = (     -W_ip2 -  4.*W_ip1 + 5.*W_i)/11.;
              // Centered d2W
              const Real d2W_ip1h = 0.5*((W_ip2 - W_ip1) - (W_i   - W_im1));
              const Real d2W_im1h = 0.5*((W_ip1 - W_i)   - (W_im1 - W_im2));

              // Compute fifth-order face values.  All use two cells and 2nd
              // derivatives
#if 1
              Real W_ip5h = 0.5*(3.*W_ip2 - W_ip1) -
                (-9.*d2W_ip5h - 32.*d2W_ip3h + d2W_ip1h )/120.;
              Real W_ip3h =
                0.5*(W_ip1 + W_ip2) - (d2W_ip1h + 18.*d2W_ip3h + d2W_ip5h)/120.;
              Real W_ip1h =
                0.5*(W_i   + W_ip1) - (d2W_im1h + 18.*d2W_ip1h + d2W_ip3h)/120.;
#else
              // Compute the second/third/fourth-order face value (full fourth-
              // order is commented out)
              Real W_ip5h = 0.5*(W_ip2 + W_ip1);
              // Real W_ip5h = (25.*W_ip2 - 23.*W_ip1 + 13.*W_i  - 3.*W_im1)/12.;
              Real W_ip3h = (2.*W_ip2 + 5.*W_ip1 - W_i)/6.;
              // Real W_ip3h = ( 3.*W_ip2 + 13.*W_ip1 -  5.*W_i  +    W_im1)/12.;
              Real W_ip1h = (   -W_ip2 + 7.*(W_ip1 +     W_i) -    W_im1)/12.;
#endif

              // // Not bad
              // const int sd2W_ip5h = 1 - 2*(int)std::signbit(d2W_ip5h);
              // const int sd2W_ip3h = 1 - 2*(int)std::signbit(d2W_ip3h);
              // const int sd2W_ip1h = 1 - 2*(int)std::signbit(d2W_ip1h);
              // if (std::abs(sd2W_ip5h + sd2W_ip3h + sd2W_ip1h) < 3)
              // if (c > CRDPhysics::pressureIndex())
                {
                  Real minW_ip5h = std::min(W_ip3, W_ip2);
                  Real maxW_ip5h = std::max(W_ip3, W_ip2);
                  W_ip5h = std::min(maxW_ip5h, std::max(minW_ip5h, W_ip5h));
                  Real minW_ip3h = std::min(W_ip2, W_ip1);
                  Real maxW_ip3h = std::max(W_ip2, W_ip1);
                  W_ip3h = std::min(maxW_ip3h, std::max(minW_ip3h, W_ip3h));
                }

              // // Face ip1h uses d2W_im1h instead of d2W_ip5h
              // const int sd2W_im1h = 1 - 2*(int)std::signbit(d2W_im1h);
              // if (std::abs(sd2W_im1h + sd2W_ip1h + sd2W_ip3h) < 3)
              // if (c > CRDPhysics::pressureIndex())
                {
                  Real minW = std::min(W_i, W_ip1);
                  Real maxW = std::max(W_i, W_ip1);
                  W_ip1h = std::min(maxW, std::max(minW, W_ip1h));
                }

              arrWFace[MD_OFFSETIX(i,+,o,c)]   = W_ip1h;
              arrWFace[MD_OFFSETIX(i,+,2*o,c)] = W_ip3h;
              arrWFace[MD_OFFSETIX(i,+,3*o,c)] = W_ip5h;
            }
        }
    }  // Loop over components
}

void
CRDutil::smoothTest(FArrayBox&           a_smoothTest,
                    const FArrayBox&     a_W,
                    const int&           a_numSlopes,
                    const int&           a_dir,
                    const Box&           a_cellBox,
                    const ProblemDomain& a_domain)
{
  Box box1cells = grow(a_cellBox, BASISV(a_dir));
  Box box1Dom = box1cells;
  box1Dom &= a_domain;
  CH_assert(a_W.box().contains(box1Dom));
  Box loBox, nextLoBox, hiBox, nextHiBox;
  Box centerBox, innerCenterBox, entireBox;
  int hasLo, hasHi;
  loHiCenter5(loBox, nextLoBox, hasLo,
              hiBox, nextHiBox, hasHi,
              centerBox, innerCenterBox, entireBox,
              box1cells, a_domain, a_dir);
  MD_ARRAY_RESTRICT(arrST, a_smoothTest);
  MD_ARRAY_RESTRICT(arrW, a_W);

  const Real hiTol = 1e100; // high tolerance to check if data is valid
  const Real loTol = -1e100; // low tolerance to check if data is valid

  const int MD_ID(o, a_dir);
  for (int c = 0; c != a_numSlopes; ++c)
    {
      if (hasLo)
        {
          MD_BOXLOOP(loBox, i)
            {
              const Real W_i   = arrW[MD_IX(i,c)];
              const Real W_ip1 = arrW[MD_OFFSETIX(i,+,1*o,c)];
              const Real W_ip2 = arrW[MD_OFFSETIX(i,+,2*o,c)];
              const Real W_ip3 = arrW[MD_OFFSETIX(i,+,3*o,c)];
              const Real W_ip4 = arrW[MD_OFFSETIX(i,+,4*o,c)];

              CH_assert(W_i   < hiTol && W_i   > loTol);
              CH_assert(W_ip1 < hiTol && W_ip1 > loTol);
              CH_assert(W_ip2 < hiTol && W_ip2 > loTol);
              CH_assert(W_ip3 < hiTol && W_ip3 > loTol);
              CH_assert(W_ip4 < hiTol && W_ip4 > loTol);

              const Real interpVal = (4.*(W_ip1 + W_ip3) - (6.*W_ip2 + W_ip4));
              Real E_i = interpVal/W_i - 1.;
              arrST[MD_IX(i,c)] = std::abs(E_i);
            }
          MD_BOXLOOP(nextLoBox, i)
            {
              const Real W_i   = arrW[MD_IX(i,c)];
              const Real W_im1 = arrW[MD_OFFSETIX(i,-,1*o,c)];
              const Real W_ip1 = arrW[MD_OFFSETIX(i,+,1*o,c)];
              const Real W_ip2 = arrW[MD_OFFSETIX(i,+,2*o,c)];
              const Real W_ip3 = arrW[MD_OFFSETIX(i,+,3*o,c)];

              CH_assert(W_im1 < hiTol && W_im1 > loTol);
              CH_assert(W_i   < hiTol && W_i   > loTol);
              CH_assert(W_ip1 < hiTol && W_ip1 > loTol);
              CH_assert(W_ip2 < hiTol && W_ip2 > loTol);
              CH_assert(W_ip3 < hiTol && W_ip3 > loTol);

              const Real interpVal = (W_im1 + 6.*W_ip1 - 4.*W_ip2 + W_ip3)/4.;
              Real E_i = interpVal/W_i - 1.;
              arrST[MD_IX(i,c)] = std::abs(E_i);
            }
        }
      MD_BOXLOOP(innerCenterBox, i)
        {
          const Real W_i   = arrW[MD_IX(i,c)];
          const Real W_im2 = arrW[MD_OFFSETIX(i,-,2*o,c)];
          const Real W_im1 = arrW[MD_OFFSETIX(i,-,1*o,c)];
          const Real W_ip1 = arrW[MD_OFFSETIX(i,+,1*o,c)];
          const Real W_ip2 = arrW[MD_OFFSETIX(i,+,2*o,c)];

          CH_assert(W_im2 < hiTol && W_im2 > loTol);
          CH_assert(W_im1 < hiTol && W_im1 > loTol);
          CH_assert(W_i   < hiTol && W_i   > loTol);
          CH_assert(W_ip1 < hiTol && W_ip1 > loTol);
          CH_assert(W_ip2 < hiTol && W_ip2 > loTol);

          const Real interpVal = (4.*(W_im1 + W_ip1) - (W_im2 + W_ip2))/6.;
          Real E_i = interpVal/W_i - 1.;
          arrST[MD_IX(i,c)] = std::abs(E_i);
        }
      if (hasHi)
        {
          MD_BOXLOOP(hiBox, i)
            {
              const Real W_i   = arrW[MD_IX(i,c)];
              const Real W_im1 = arrW[MD_OFFSETIX(i,-,1*o,c)];
              const Real W_im2 = arrW[MD_OFFSETIX(i,-,2*o,c)];
              const Real W_im3 = arrW[MD_OFFSETIX(i,-,3*o,c)];
              const Real W_im4 = arrW[MD_OFFSETIX(i,-,4*o,c)];

              CH_assert(W_im4 < hiTol && W_im4 > loTol);
              CH_assert(W_im3 < hiTol && W_im3 > loTol);
              CH_assert(W_im2 < hiTol && W_im2 > loTol);
              CH_assert(W_im1 < hiTol && W_im1 > loTol);
              CH_assert(W_i   < hiTol && W_i   > loTol);

              const Real interpVal = (4.*(W_im1 + W_im3) - (6.*W_im2 + W_im4));
              Real E_i = interpVal/W_i - 1.;
              arrST[MD_IX(i,c)] = std::abs(E_i);
            }
          MD_BOXLOOP(nextHiBox, i)
            {
              const Real W_i   = arrW[MD_IX(i,c)];
              const Real W_ip1 = arrW[MD_OFFSETIX(i,+,1*o,c)];
              const Real W_im1 = arrW[MD_OFFSETIX(i,-,1*o,c)];
              const Real W_im2 = arrW[MD_OFFSETIX(i,-,2*o,c)];
              const Real W_im3 = arrW[MD_OFFSETIX(i,-,3*o,c)];

              CH_assert(W_im3 < hiTol && W_im3 > loTol);
              CH_assert(W_im2 < hiTol && W_im2 > loTol);
              CH_assert(W_im1 < hiTol && W_im1 > loTol);
              CH_assert(W_i   < hiTol && W_i   > loTol);
              CH_assert(W_ip1 < hiTol && W_ip1 > loTol);

              const Real interpVal = (W_ip1 + 6.*W_im1 - 4.*W_im2 + W_im3)/4.;
              Real E_i = interpVal/W_i - 1.;
              arrST[MD_IX(i,c)] = std::abs(E_i);
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  C++ implementation of full PPM limiter
/**
 *//*-----------------------------------------------------------------*/

void
CRDutil::PPMLimiter(FArrayBox&           a_dWMinus,
                    FArrayBox&           a_dWPlus,
                    const FArrayBox&     a_W,
                    const int            a_numSlopes,
                    const int            a_dir,
                    const Box&           a_box,
                    const ProblemDomain& a_domain)
{
  CH_TIME("CRDutil::PPMLimiter");
  // We calculate d2Wfcf on a_box,
  // and we need a_dWMinus on a_box,
  // and we need a_dWPlus on a_box.
  // d2Wfcf[i] = 6 * (a_dWMinus[i] + a_dWPlus[i])
  //           = 6 * (thisFaceWDir[i-e/2] - a_cellW[i] +
  //                  thisFaceWDir[i+e/2] - a_cellW[i])
  FABSTACKTEMP(d2Wfcf, a_box, a_numSlopes);
  for (int c = 0; c != a_numSlopes; ++c)
    {
      MD_BOXLOOP(a_box, i)
        {
          d2Wfcf[MD_IX(i, c)] =
            6*(a_dWMinus[MD_IX(i, c)] + a_dWPlus[MD_IX(i, c)]);
        }
    }

  // In order to get a_dWMinus and a_dWPlus on a_box,
  // we need d2W on a_box grown by 3 in a_dir directions.
  Box box3 = a_box;
  box3.grow(a_dir, 3);
  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;
  loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox, box3,
             a_domain, a_dir);

  FABSTACKTEMP(d2W, entireBox, a_numSlopes);
  // On centerBox, use 3-point stencil for d2W;
  // on loBox and hiBox, copy result from neighbor.
  // In the case of centerBox == entireBox == box1,
  // where box1 is a_box grown by 1 in dimension a_dir:
  // We calculate d2W on a_box grown by 2 in dimension a_dir,
  // and we need a_W on a_box grown by 3 in dimension a_dir.
  CRDutil::getSecondDiff(d2W,
                         a_W,
                         a_numSlopes,
                         a_dir,
                         loBox,
                         hasLo,
                         hiBox,
                         hasHi,
                         centerBox);

  Box box1 = a_box;
  box1.grow(a_dir, 1);
  Box nextLoBox, nextHiBox, innerCenterBox;
  loHiCenter5(loBox, nextLoBox, hasLo,
              hiBox, nextHiBox, hasHi,
              centerBox, innerCenterBox, entireBox,
              box1, a_domain, a_dir);
  CH_assert(entireBox == a_box);

  const Real limitC = 1.25;
  const Real eps = 1.0e-12;
  const Real c3 = 0.1;
  const int useHOChecks = 1;        // Use 3rd-derivative checks in center
  const int useLoBoundHOChecks = 1; // Always use 3rd-derivative checks
  const int useHiBoundHOChecks = 1; // near hi and lo boundaries
  CRDutil::checkCubicLimiter(a_dWMinus,
                             a_dWPlus,
                             a_W,
                             d2W,
                             d2Wfcf,
                             a_numSlopes,
                             a_dir,
                             loBox,
                             nextLoBox,
                             hasLo,
                             hiBox,
                             nextHiBox,
                             hasHi,
                             innerCenterBox,
                             limitC,
                             c3,
                             eps,
                             useHOChecks,
                             useLoBoundHOChecks,
                             useHiBoundHOChecks);
}

/*--------------------------------------------------------------------*/
//  C++ implementation of getSecondDiff
/**
 *//*-----------------------------------------------------------------*/

void
CRDutil::getSecondDiff(FArrayBox&       a_d2W,
                       const FArrayBox& a_W,
                       const int        a_numSlopes,
                       const int        a_dir,
                       const Box&       a_loBox,
                       const int        a_hasLo,
                       const Box&       a_hiBox,
                       const int        a_hasHi,
                       const Box&       a_centerBox)
{
  CH_TIME("CRDutil::getSecondDiff");
  const int MD_ID(o, a_dir);
  const Real hiTol = 1e100; // high tolerance to check if data is valid
  const Real loTol = -1e100; // low tolerance to check if data is valid
  for (int c = 0; c != a_numSlopes; ++c)
    {
      MD_BOXLOOP(a_centerBox, i)
        {
          const Real W_l = a_W[MD_OFFSETIX(i,-,o, c)];
          const Real W_c = a_W[MD_IX(i, c)];
          const Real W_r = a_W[MD_OFFSETIX(i,+,o, c)];

          CH_assert(W_l < hiTol && W_l > loTol);
          CH_assert(W_c < hiTol && W_c > loTol);
          CH_assert(W_r < hiTol && W_r > loTol);

          const Real dW_l = W_c - W_l;
          const Real dW_r = W_r - W_c;
          a_d2W[MD_IX(i, c)] = dW_r - dW_l;
        }
      if (a_hasLo == 1)
        {
          // Left-sided coefficients:
          // 2, -5,  4, -1   which can be obtained from the sum of
          // 1, -2,  1       multiplied by 2, and
          //     1, -2,  1   multiplied by -1.
          MD_BOXLOOP(a_loBox, i)
            {
              a_d2W[MD_IX(i, c)] =
                2*a_d2W[MD_OFFSETIX(i,+,o, c)] - a_d2W[MD_OFFSETIX(i,+,2*o, c)];
              // This gives only FIRST-order approx to second derivative:
              // = a_d2W[MD_OFFSETIX(i,+,o, c)];

              CH_assert(a_d2W[MD_OFFSETIX(i,+,o, c)] < hiTol &&
                        a_d2W[MD_OFFSETIX(i,+,o, c)] > loTol);
              CH_assert(a_d2W[MD_OFFSETIX(i,+,2*o, c)] < hiTol &&
                        a_d2W[MD_OFFSETIX(i,+,2*o, c)] > loTol);
            }
        }
      if (a_hasHi == 1)
        {
          // Right-sided coefficients:
          // -1,  4, -5, 2   which can be obtained from the sum of
          //  1, -2,  1      multiplied by -1, and
          //      1, -2, 1   multiplied by 2.
          MD_BOXLOOP(a_hiBox, i)
            {
              a_d2W[MD_IX(i, c)] =
                2*a_d2W[MD_OFFSETIX(i,-,o, c)] - a_d2W[MD_OFFSETIX(i,-,2*o, c)];
              // This gives only FIRST-order approx to second derivative:
              // = a_d2W[MD_OFFSETIX(i,-,o, c)];

              CH_assert(a_d2W[MD_OFFSETIX(i,-,o, c)] < hiTol &&
                        a_d2W[MD_OFFSETIX(i,-,o, c)] > loTol);
              CH_assert(a_d2W[MD_OFFSETIX(i,-,2*o, c)] < hiTol &&
                        a_d2W[MD_OFFSETIX(i,-,2*o, c)] > loTol);
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  C++ implementation of checkCubicLimiter
/**
 *//*-----------------------------------------------------------------*/

void
CRDutil::checkCubicLimiter(FArrayBox&       a_dWMinus,
                           FArrayBox&       a_dWPlus,
                           const FArrayBox& a_W,
                           const FArrayBox& a_d2W,
                           const FArrayBox& a_dW2fcf,
                           const int        a_numSlopes,
                           const int        a_dir,
                           const Box&       a_loBox,
                           const Box&       a_nextLoBox,
                           const int        a_hasLo,
                           const Box&       a_hiBox,
                           const Box&       a_nextHiBox,
                           const int        a_hasHi,
                           const Box&       a_innerCenterBox,
                           const Real       a_limitC,
                           const Real       a_C3,
                           const Real       a_eps,
                           const int        a_useHOCheck,
                           const int        a_loBoundHOCheck,
                           const int        a_hiBoundHOCheck)
{
  CH_TIME("CRDutil::checkCubicLimiter");
  const int MD_ID(o, a_dir);
  const Real hiTol = 1e100; // high tolerance to check if data is valid
  const Real loTol = -1e100; // low tolerance to check if data is valid
  for (int c = 0; c != a_numSlopes; ++c)
    {
      MD_BOXLOOP(a_innerCenterBox, i)
        {
          // dWM = W[i - e/2] - W[i] == wfaceLo - wavg == -dwfm
          // dWP = W[i + e/2] - W[i] == wfaceHi - wavg == dwfp
          const Real dWM = a_dWMinus[MD_IX(i, c)];
          const Real dWP = a_dWPlus [MD_IX(i, c)];
          const bool bigM = std::abs(dWM) > 2*std::abs(dWP);
          const bool bigP = std::abs(dWP) > 2*std::abs(dWM);
          const Real WLL = a_W[MD_OFFSETIX(i,-,2*o, c)];
          const Real WC  = a_W[MD_IX(i, c)];
          const Real WRR = a_W[MD_OFFSETIX(i,+,2*o, c)];
          const Real dWavgM = WC - WLL;
          const Real dWavgP = WRR - WC;
          const Real prodE1 = dWM*dWP;
          const Real prodE2 = dWavgM*dWavgP;

          CH_assert(dWM < hiTol && dWM > loTol);
          CH_assert(dWP < hiTol && dWP > loTol);
          CH_assert(WLL < hiTol && WLL > loTol);
          CH_assert(WC  < hiTol && WC  > loTol);
          CH_assert(WRR < hiTol && WRR > loTol);

          // prodExtr1 == dwfm * dwfp == -dWM * dWP == -prodE1
          // prodExtr2 == dWavgM * dWavgP
          // condExtr == (prodExtr1 <= 0) | (prodExtr2 <= 0)
          //          == (-prodE1 <= 0) | (prodE2 <= 0)
          // Extremum check.
          if (prodE1 >= 0. || prodE2 <= 0.)
            {
              // This is an extremum.  Check if it is a potential discontinuity,
              // by checking relative sizes and signs of second derivatives.
              // d2WL == d2wm
              // d2WC == d2wm
              // d2WR == d2wp
              const Real d2WL  = a_d2W[MD_OFFSETIX(i,-,o, c)];
              const Real d2WC  = a_d2W[MD_IX(i, c)];
              const Real d2WR  = a_d2W[MD_OFFSETIX(i,+,o, c)];
              const Real atfcf = a_dW2fcf[MD_IX(i, c)];

              CH_assert(d2WL  < hiTol && d2WL  > loTol);
              CH_assert(d2WC  < hiTol && d2WC  > loTol);
              CH_assert(d2WR  < hiTol && d2WR  > loTol);
              CH_assert(atfcf < hiTol && atfcf > loTol);

              // Set d2Wlim = stuff, if signs are all same; 0, otherwise.
              // Set rho = d2Wlim/d2Wfcf, if |d2Wfcf| >= eps; 0, otherwise.
              // So rho is set to 0 unless signs are all same AND
              // |d2Wfcf| >= eps.
              Real d2Wlim = 0.;
              Real rho = 0.;
              if (std::abs(atfcf) >= a_eps)
                {
                  const int sd2WL  = (int)std::copysign((Real)1, d2WL);
                  const int sd2WC  = (int)std::copysign((Real)1, d2WC);
                  const int sd2WR  = (int)std::copysign((Real)1, d2WR);
                  const int sd2fcf = (int)std::copysign((Real)1, atfcf);
                  if (std::abs(sd2WL + sd2WC + sd2WR + sd2fcf) == 4)
                    // Signs are all the same, sd2WC
                    {
                      d2Wlim = sd2WC*std::min({
                          std::abs(atfcf),
                          a_limitC*std::abs(d2WL),
                          a_limitC*std::abs(d2WC),
                          a_limitC*std::abs(d2WR) });
                      rho = d2Wlim/atfcf;
                    }
                  // else, d2Wlim = 0, and rho = 0.
                }
              if (rho < (1. - a_eps))
                // Potential discontinuity:
                // Check if well-separated from a small perturbation of a cubic.
                {
                  const Real d2WLL = a_d2W[MD_OFFSETIX(i,-,2*o, c)];
                  const Real d2WRR = a_d2W[MD_OFFSETIX(i,+,2*o, c)];
                  // Third derivatives at faces.
                  const Real d3WLL = d2WL  - d2WLL;
                  const Real d3WL  = d2WC  - d2WL;
                  const Real d3WR  = d2WR  - d2WC;
                  const Real d3WRR = d2WRR - d2WR;
                  const Real d3Wmin = std::min({ d3WLL, d3WL, d3WR, d3WRR });
                  const Real d3Wmax = std::max({ d3WLL, d3WL, d3WR, d3WRR });
                  const Real prodD3 =
                    a_C3*std::max(std::abs(d3Wmax), std::abs(d3Wmin)) -
                    std::abs(d3Wmax - d3Wmin);
                  // Is the extremum is located close to where d2W is nearly
                  // linear?
                  if (prodD3 <= 0. || a_useHOCheck == 0)
                    // The answer is no.
                    {
                      // We are well-separated from a cubic, so we can apply
                      // the limiter.
                      // dwfm * dwfp == -dWM * dWP
                      if (prodE1 > 0.)
                        {
                          // extrapLo == wfaceLo == wavg - dwfm == wavg + dWM
                          // New extrapLo == wavg - dwfm*rho == wavg + dWM*rho
                          // Hence new dWM == extrapLo - wavg == dWM*rho
                          
                          // extrapHi == wfaceHi == wavg + dwfp == wavg + dWP
                          // New extrapHi == wavg + dwfp*rho == wavg + dWP*rho
                          // Hence new dWP == extrapHi - wavg == dWP*rho
                          a_dWMinus[MD_IX(i, c)] = dWM*rho;
                          a_dWPlus [MD_IX(i, c)] = dWP*rho;
                        }
                      else if (bigM)
                        {
                          // extrapLo == wfaceLo == wavg - dwfm == wavg + dWM
                          // New extrapLo
                          //   == wavg - (2*dwfp*(1 - rho) + dwfm*rho)
                          //   == wavg - (2*dWP*(1 - rho) - dWM*rho)
                          // Hence new dWM == extrapLo - wavg
                          //               == dWM*rho - 2*dWP*(1 - rho)
                          a_dWMinus[MD_IX(i, c)] = dWM*rho - 2*dWP*(1. - rho);
                        }
                      else if (bigP)
                        {
                          // extrapHi == wfaceHi == wavg + dwfp == wavg + dWP
                          // New extrapHi
                          //   == wavg + (2*dwfm*(1 - rho) + dwfp*rho)
                          //   == wavg + (-2*dWM*(1 - rho) + dWP*rho)
                          // Hence new dWP == extrapHi - wavg
                          //               == dWP*rho - 2*dWM*(1 - rho)
                          a_dWPlus [MD_IX(i, c)] = dWP*rho - 2*dWM*(1. - rho);
                        }
                    }
                }
            }  // If extremum
          else
            {
              // Not an extremum, so apply PPM limiter.
               if (bigM)
                 {
                   a_dWMinus[MD_IX(i, c)] = -2*dWP;
                 }
               if (bigP)
                 {
                   a_dWPlus [MD_IX(i, c)] = -2*dWM;
                 }
            }
        }  // Loop over innerCenterBox

//--Low side of box

      if (a_hasLo)
        {
          MD_BOXLOOP(a_loBox, i)
            {
              const Real dWM = a_dWMinus[MD_IX(i, c)];
              const Real dWP = a_dWPlus [MD_IX(i, c)];
              const bool bigM = std::abs(dWM) > 2*std::abs(dWP);
              const bool bigP = std::abs(dWP) > 2*std::abs(dWM);
              const Real prodE1 = dWM*dWP;
              if (prodE1 >= 0.)
                {
                  const Real d2WC  = a_d2W[MD_IX(i, c)];
                  const Real d2WR  = a_d2W[MD_OFFSETIX(i,+,o, c)];
                  const Real atfcf = a_dW2fcf[MD_IX(i, c)];
                  Real d2Wlim = 0.;
                  Real rho = 0.;
                  if (std::abs(atfcf) >= a_eps)
                    {
                      const int sd2WC  = (int)std::copysign((Real)1, d2WC);
                      const int sd2WR  = (int)std::copysign((Real)1, d2WR);
                      const int sd2fcf = (int)std::copysign((Real)1, atfcf);
                      if (std::abs(sd2WC + sd2WR + sd2fcf) == 3)
                        {
                          d2Wlim = sd2WC*std::min({
                              std::abs(atfcf),
                              a_limitC*std::abs(d2WC),
                              a_limitC*std::abs(d2WR) });
                          rho = d2Wlim/atfcf;
                        }
                    }
                  if (rho < (1. - a_eps))
                    {
                      const Real d2WRR = a_d2W[MD_OFFSETIX(i,+,2*o, c)];
                      const Real d3WR  = d2WR  - d2WC;
                      const Real d3WRR = d2WRR - d2WR;
                      const Real d3Wmin = std::min(d3WR, d3WRR);
                      const Real d3Wmax = std::max(d3WR, d3WRR);
                      const Real prodD3 =
                        a_C3*std::max(std::abs(d3Wmax), std::abs(d3Wmin)) -
                        std::abs(d3Wmax - d3Wmin);
                      if (prodD3 <= 0. || a_loBoundHOCheck == 0)
                        {
                          if (prodE1 > 0.)
                            {
                              a_dWMinus[MD_IX(i, c)] = dWM*rho;
                              a_dWPlus [MD_IX(i, c)] = dWP*rho;
                            }
                          else if (bigM)
                            {
                              a_dWMinus[MD_IX(i, c)] =
                                dWM*rho - 2*dWP*(1. - rho);
                            }
                          else if (bigP)
                            {
                              a_dWPlus [MD_IX(i, c)] =
                                dWP*rho - 2*dWM*(1. - rho);
                            }
                        }
                    }
                }
              else
                {
                  if (bigM)
                    {
                      a_dWMinus[MD_IX(i, c)] = -2*dWP;
                    }
                  if (bigP)
                    {
                      a_dWPlus [MD_IX(i, c)] = -2*dWM;
                    }
                }
            }  // Loop over low cells

//--Next to low side of box

          MD_BOXLOOP(a_nextLoBox, i)
            {
              const Real dWM = a_dWMinus[MD_IX(i, c)];
              const Real dWP = a_dWPlus [MD_IX(i, c)];
              const bool bigM = std::abs(dWM) > 2*std::abs(dWP);
              const bool bigP = std::abs(dWP) > 2*std::abs(dWM);
              const Real prodE1 = dWM*dWP;
              if (prodE1 >= 0.)
                {
                  const Real d2WL  = a_d2W[MD_OFFSETIX(i,-,o, c)];
                  const Real d2WC  = a_d2W[MD_IX(i, c)];
                  const Real d2WR  = a_d2W[MD_OFFSETIX(i,+,o, c)];
                  const Real atfcf = a_dW2fcf[MD_IX(i, c)];
                  Real d2Wlim = 0.;
                  Real rho = 0.;
                  if (std::abs(atfcf) >= a_eps)
                    {
                      const int sd2WL  = (int)std::copysign((Real)1, d2WL);
                      const int sd2WC  = (int)std::copysign((Real)1, d2WC);
                      const int sd2WR  = (int)std::copysign((Real)1, d2WR);
                      const int sd2fcf = (int)std::copysign((Real)1, atfcf);
                      if (std::abs(sd2WL + sd2WC + sd2WR + sd2fcf) == 4)
                        {
                          d2Wlim = sd2WC*std::min({
                              std::abs(atfcf),
                              a_limitC*std::abs(d2WL),
                              a_limitC*std::abs(d2WC),
                              a_limitC*std::abs(d2WR) });
                          rho = d2Wlim/atfcf;
                        }
                    }
                  if (rho < (1. - a_eps))
                    {
                      const Real d2WRR = a_d2W[MD_OFFSETIX(i,+,2*o, c)];
                      const Real d3WL  = d2WC  - d2WL;
                      const Real d3WR  = d2WR  - d2WC;
                      const Real d3WRR = d2WRR - d2WR;
                      const Real d3Wmin = std::min({ d3WL, d3WR, d3WRR });
                      const Real d3Wmax = std::max({ d3WL, d3WR, d3WRR });
                      const Real prodD3 =
                        a_C3*std::max(std::abs(d3Wmax), std::abs(d3Wmin)) -
                        std::abs(d3Wmax - d3Wmin);
                      if (prodD3 <= 0. || a_loBoundHOCheck == 0)
                        {
                          if (prodE1 > 0.)
                            {
                              a_dWMinus[MD_IX(i, c)] = dWM*rho;
                              a_dWPlus [MD_IX(i, c)] = dWP*rho;
                            }
                          else if (bigM)
                            {
                              a_dWMinus[MD_IX(i, c)] =
                                dWM*rho - 2*dWP*(1. - rho);
                            }
                          else if (bigP)
                            {
                              a_dWPlus [MD_IX(i, c)] =
                                dWP*rho - 2*dWM*(1. - rho);
                            }
                        }
                    }
                }
              else
                {
                  if (bigM)
                    {
                      a_dWMinus[MD_IX(i, c)] = -2*dWP;
                    }
                  if (bigP)
                    {
                      a_dWPlus [MD_IX(i, c)] = -2*dWM;
                    }
                }
            }  // Loop over next-to-low cells
        }  // Has low cells

//--High side of box

      if (a_hasHi)
        {
          MD_BOXLOOP(a_hiBox, i)
            {
              const Real dWM = a_dWMinus[MD_IX(i, c)];
              const Real dWP = a_dWPlus [MD_IX(i, c)];
              const bool bigM = std::abs(dWM) > 2*std::abs(dWP);
              const bool bigP = std::abs(dWP) > 2*std::abs(dWM);
              const Real prodE1 = dWM*dWP;
              if (prodE1 >= 0.)
                {
                  const Real d2WL  = a_d2W[MD_OFFSETIX(i,-,o, c)];
                  const Real d2WC  = a_d2W[MD_IX(i, c)];
                  const Real atfcf = a_dW2fcf[MD_IX(i, c)];
                  Real d2Wlim = 0.;
                  Real rho = 0.;
                  if (std::abs(atfcf) >= a_eps)
                    {
                      const int sd2WL  = (int)std::copysign((Real)1, d2WL);
                      const int sd2WC  = (int)std::copysign((Real)1, d2WC);
                      const int sd2fcf = (int)std::copysign((Real)1, atfcf);
                      if (std::abs(sd2WL + sd2WC + sd2fcf) == 3)
                        {
                          d2Wlim = sd2WC*std::min({
                              std::abs(atfcf),
                              a_limitC*std::abs(d2WL),
                              a_limitC*std::abs(d2WC) });
                          rho = d2Wlim/atfcf;
                        }
                    }
                  if (rho < (1. - a_eps))
                    {
                      const Real d2WLL = a_d2W[MD_OFFSETIX(i,-,2*o, c)];
                      const Real d3WLL = d2WL  - d2WLL;
                      const Real d3WL  = d2WC  - d2WL;
                      const Real d3Wmin = std::min(d3WLL, d3WL);
                      const Real d3Wmax = std::max(d3WLL, d3WL);
                      const Real prodD3 =
                        a_C3*std::max(std::abs(d3Wmax), std::abs(d3Wmin)) -
                        std::abs(d3Wmax - d3Wmin);
                      if (prodD3 <= 0. || a_hiBoundHOCheck == 0)
                        {
                          if (prodE1 > 0.)
                            {
                              a_dWMinus[MD_IX(i, c)] = dWM*rho;
                              a_dWPlus [MD_IX(i, c)] = dWP*rho;
                            }
                          else if (bigM)
                            {
                              a_dWMinus[MD_IX(i, c)] =
                                dWM*rho - 2*dWP*(1. - rho);
                            }
                          else if (bigP)
                            {
                              a_dWPlus [MD_IX(i, c)] =
                                dWP*rho - 2*dWM*(1. - rho);
                            }
                        }
                    }
                }
              else
                {
                  if (bigM)
                    {
                      a_dWMinus[MD_IX(i, c)] = -2*dWP;
                    }
                  if (bigP)
                    {
                      a_dWPlus [MD_IX(i, c)] = -2*dWM;
                    }
                }
            }  // Loop over high cells

//--Next to high side of box

          MD_BOXLOOP(a_nextHiBox, i)
            {
              const Real dWM = a_dWMinus[MD_IX(i, c)];
              const Real dWP = a_dWPlus [MD_IX(i, c)];
              const bool bigM = std::abs(dWM) > 2*std::abs(dWP);
              const bool bigP = std::abs(dWP) > 2*std::abs(dWM);
              const Real prodE1 = dWM*dWP;
              if (prodE1 >= 0.)
                {
                  const Real d2WL  = a_d2W[MD_OFFSETIX(i,-,o, c)];
                  const Real d2WC  = a_d2W[MD_IX(i, c)];
                  const Real d2WR  = a_d2W[MD_OFFSETIX(i,+,o, c)];
                  const Real atfcf = a_dW2fcf[MD_IX(i, c)];
                  Real d2Wlim = 0.;
                  Real rho = 0.;
                  if (std::abs(atfcf) >= a_eps)
                    {
                      const int sd2WL  = (int)std::copysign((Real)1, d2WL);
                      const int sd2WC  = (int)std::copysign((Real)1, d2WC);
                      const int sd2WR  = (int)std::copysign((Real)1, d2WR);
                      const int sd2fcf = (int)std::copysign((Real)1, atfcf);
                      if (std::abs(sd2WL + sd2WC + sd2WR + sd2fcf) == 4)
                        {
                          d2Wlim = sd2WC*std::min({
                              std::abs(atfcf),
                              a_limitC*std::abs(d2WL),
                              a_limitC*std::abs(d2WC),
                              a_limitC*std::abs(d2WR) });
                          rho = d2Wlim/atfcf;
                        }
                    }
                  if (rho < (1. - a_eps))
                    {
                      const Real d2WLL = a_d2W[MD_OFFSETIX(i,-,2*o, c)];
                      const Real d3WLL = d2WL  - d2WLL;
                      const Real d3WL  = d2WC  - d2WL;
                      const Real d3WR  = d2WR  - d2WC;
                      const Real d3Wmin = std::min({ d3WLL, d3WL, d3WR });
                      const Real d3Wmax = std::max({ d3WLL, d3WL, d3WR });
                      const Real prodD3 =
                        a_C3*std::max(std::abs(d3Wmax), std::abs(d3Wmin)) -
                        std::abs(d3Wmax - d3Wmin);
                      if (prodD3 <= 0. || a_hiBoundHOCheck == 0)
                        {
                          if (prodE1 > 0.)
                            {
                              a_dWMinus[MD_IX(i, c)] = dWM*rho;
                              a_dWPlus [MD_IX(i, c)] = dWP*rho;
                            }
                          else if (bigM)
                            {
                              a_dWMinus[MD_IX(i, c)] =
                                dWM*rho - 2*dWP*(1. - rho);
                            }
                          else if (bigP)
                            {
                              a_dWPlus [MD_IX(i, c)] =
                                dWP*rho - 2*dWM*(1. - rho);
                            }
                        }
                    }
                }
              else
                {
                  if (bigM)
                    {
                      a_dWMinus[MD_IX(i, c)] = -2*dWP;
                    }
                  if (bigP)
                    {
                      a_dWPlus [MD_IX(i, c)] = -2*dWM;
                    }
                }
            }  // Loop over next-to-high cells
        }  // Has high cells

    }  // Loop over components
}

/*--------------------------------------------------------------------*/
//  Apply extra limiting at some boundaries
/** \param[in]  a_WMinus
 *                      Difference from cell average to get the face
 *                      value on the left, e.g., a first-order scheme
 *                      would have a value of 0.  Stored on a cell-
 *                      centered data structure
 *  \param[out] a_WMinus
 *                      Extra limiting applied on faces orthogonal to
 *                      and adjacent to certain physical boundaries.
 *  \param[in]  a_WPlus Difference from cell average to get the face
 *                      value on the right, e.g., a first-order scheme
 *                      would have a value of 0.  Stored on a cell-
 *                      centered data structure
 *  \param[out] a_WPlus Extra limiting applied on faces orthogonal to
 *                      and adjacent to certain physical boundaries.
 *//*-----------------------------------------------------------------*/

void
CRDutil::extraBoundaryLimiting(FArrayBox&           a_WMinus,
                               FArrayBox&           a_WPlus,
                               const FArrayBox&     a_W,
                               const Box&           a_cellBox,
                               const int&           a_numSlopes,
                               const int&           a_dir,
                               const int&           a_level,
                               const Real&          a_time,
                               LevelGridMetrics&    a_gridMetrics,
                               const ProblemDomain& a_domain)
{
  CH_TIME("CRDutil::extraBoundaryLimiting");
  // Apply to first interior face adjacent to any inlet/outlet/farfield BC
  const int MD_ID(o, a_dir);
  CRDparam::DomainBCType limitBCType =
    static_cast<CRDparam::DomainBCType>(CRDparam::DomainBCTypeInOut);
  Vector<Box> boundBoxes;
  for (const auto side : EachSide)
    {
      boundBoxes.clear();
      CRDparam::g_CNSIBC->getBoundaryAdjCells(boundBoxes,
                                              a_cellBox,
                                              a_domain,
                                              a_dir,
                                              side,
                                              a_level,
                                              a_time,
                                              // A setting of 2 modifies 1 face
                                              2, // Return first two cells
                                              a_gridMetrics,
                                              limitBCType);

      for (int bcNum = 0; bcNum != boundBoxes.size(); ++bcNum)
        {
          Box curBox(boundBoxes[bcNum]);
          if (curBox.isEmpty()) continue;
          Box faceBox(curBox);
          faceBox.surroundingNodes(a_dir);
          // This grow shrinks on each side.  So you need at least two cells
          // to get a single face.
          faceBox.grow(-IntVect_basis(a_dir));
          a_WMinus.shiftHalf(a_dir, -1);
          a_WPlus .shiftHalf(a_dir, +1);
          for (int c = 0; c != a_numSlopes; ++c)
            {
              MD_BOXLOOP(faceBox, i)
                {
                  const Real W_im1  = a_W[MD_OFFSETIX(i,-,1*o,c)];
                  const Real W_i    = a_W[MD_IX(i,c)];
                  const Real W_im1h = 0.5*(W_im1 + W_i);

                  const Real Wm = a_WMinus[MD_IX(i,c)];
                  Real dWm_im1h = W_im1h - W_i;
                  if (dWm_im1h*Wm < 0.)
                    // The computed delta Wm is opposite in sign to the
                    // difference between cells
                    {
                      dWm_im1h = 0.;
                    }
                  if (std::abs(dWm_im1h) < std::abs(Wm))
                    // Take the smaller delta
                    {
                      a_WMinus[MD_IX(i,c)] = dWm_im1h;
                    }

                  const Real Wp = a_WPlus [MD_IX(i,c)];
                  Real dWp_im1h = W_im1h - W_im1;
                  if (dWp_im1h*Wp < 0.)
                    {
                    // The computed delta Wp is opposite in sign to the
                    // difference between cells
                      dWp_im1h = 0.;
                    }
                  if (std::abs(dWp_im1h) < std::abs(Wp))
                    // Take the smaller delta
                    {
                      a_WPlus[MD_IX(i,c)] = dWp_im1h;
                    }
                }
            }
          a_WMinus.shiftHalf(a_dir, +1);
          a_WPlus .shiftHalf(a_dir, -1);
          // // Severe
          // a_WMinus.setVal(0., curBox, 0, a_numSlopes);
          // a_WPlus.setVal(0., curBox, 0, a_numSlopes);
        }
    }
}

void
CRDutil::PPMUpwindFaceValues(FArrayBox&           a_WfaceMinus,
                             FArrayBox&           a_WfacePlus,
                             const FArrayBox&     a_WcellAvg,
                             const int&           a_numComp,
                             const int&           a_dir,
                             const Box&           a_box,
                             const ProblemDomain& a_domain)
{
  CH_TIME("CRDutil::PPMUpwindFaceValues");
  //**FIXME only valid for 64-bit
  const Real eps = 1.E11*std::numeric_limits<Real>::epsilon();
  CH_assert(a_WcellAvg.contains(a_box));
  CH_assert(a_WfaceMinus.contains(a_box));
  CH_assert(a_WfacePlus.contains(a_box));

  constexpr Real fact = 1./60.;
  const int stenSize = 2;

  IntVect growVect(IntVect::Zero);
  growVect[a_dir] = stenSize;
  const Box box2 = grow(a_box, growVect) & a_domain; // Touched cells a_WcellAvg
  Box box2Lo(a_box);
  box2Lo.growLo(a_dir, stenSize) & a_domain;
  Box box2Hi(a_box);
  box2Hi.growHi(a_dir, stenSize) & a_domain;
  growVect[a_dir] = -stenSize;
  const Box cenBox = grow(box2, growVect);     // Cells in result updated by
                                               // a centered stencil

  Box nearLoCells, farLoCells;
  Box nearHiCells, farHiCells;
  int hasLoCells = 0;
  int hasHiCells = 0;

  const int MD_ID(o, a_dir);

  CH_assert(!cenBox.isEmpty());
  if (cenBox.smallEnd(a_dir) != a_box.smallEnd(a_dir))
    {
      hasLoCells = 1;
      farLoCells = bdryLo(box2, a_dir, 1);
      farLoCells.shiftHalf(a_dir, 1);
      nearLoCells = adjCellHi(farLoCells, a_dir, 1);
    }
  if (cenBox.bigEnd(a_dir) != a_box.bigEnd(a_dir))
    {
      hasHiCells = 1;
      farHiCells = bdryHi(box2, a_dir, 1);
      farHiCells.shiftHalf(a_dir, -1);
      nearHiCells = adjCellLo(farHiCells, a_dir, 1);
    }

//--Box where centered stencil is applied

  for (int c = 0; c != a_numComp; ++c)
    {

//--Lower boundary adjacent cells

      if (hasLoCells)
        {
          // FIXME: Find proper limiting for boundary faces
          MD_BOXLOOP(farLoCells, i)
            {
              Real W_i   = a_WcellAvg[MD_IX(i, c)];
              Real W_ip1 = a_WcellAvg[MD_OFFSETIX(i,+,1*o, c)];
              Real W_ip2 = a_WcellAvg[MD_OFFSETIX(i,+,2*o, c)];
              Real W_ip3 = a_WcellAvg[MD_OFFSETIX(i,+,3*o, c)];
              Real W_ip4 = a_WcellAvg[MD_OFFSETIX(i,+,4*o, c)];
              Real vals[4] = {W_i, W_ip1, W_ip2, W_ip3};

              Real WfaceL = (12.*W_ip4 - 63.*W_ip3 - 163.*W_ip1
                             + 137.*(W_ip2 + W_i))*fact;
              Real WfaceL3 = (2.*W_ip2 - 7.*W_ip1 + 11.*W_i)/6.;

              Real WfaceR = (-3.*W_ip4 + 17.*W_ip3 - 43.*W_ip2
                             + 77.*W_ip1 + 12.*W_i)*fact;
              Real WfaceR3 = (-W_ip2 + 5.*W_ip1 + 2.*W_i)/6.;

              Real avg = 0.;
              for (int tc = 0; tc != 4; ++tc)
                {
                  avg += vals[tc];
                }
              avg /= 4.;
              Real varianceL = std::abs(avg - WfaceL);
              Real varianceL3 = std::abs(avg - WfaceL3);
              Real varianceR = std::abs(avg - WfaceR);
              Real varianceR3 = std::abs(avg- WfaceR3);

              // FIXME: This method selects the face construction
              // with the least amount of variance relative to adjacent cells
              if (varianceL < varianceL3)
                {
                  a_WfaceMinus[MD_IX(i, c)] = WfaceL;
                }
              else
                {
                  a_WfaceMinus[MD_IX(i, c)] = WfaceL3;
                }
              if (varianceR < varianceR3)
                {
                  a_WfacePlus[MD_IX(i, c)] = WfaceR;
                }
              else
                {
                  a_WfacePlus[MD_IX(i, c)] = WfaceR3;
                }
            }
          MD_BOXLOOP(nearLoCells, i)
            {
              Real W_im1 = a_WcellAvg[MD_OFFSETIX(i,-,1*o, c)];
              Real W_i   = a_WcellAvg[MD_IX(i, c)];
              Real W_ip1 = a_WcellAvg[MD_OFFSETIX(i,+,1*o, c)];
              Real W_ip2 = a_WcellAvg[MD_OFFSETIX(i,+,2*o, c)];
              Real W_ip3 = a_WcellAvg[MD_OFFSETIX(i,+,3*o, c)];

              const Real tol = eps*std::abs(W_i);

//--Left face (Wminus)

              // Real WfaceL = (-3.*W_ip3 + 17.*W_ip2 - 43.*W_ip1 + 77.*W_i
              //                + 12.*W_im1)*fact;
              Real WfaceL3 = (-W_ip1 + 5*W_i + 2*W_im1)/6.0;
              // if (std::abs(WfaceL - WfaceL3) > tol &&
              //     (WfaceL - W_im1)*(W_i - WfaceL) < 0.)
              //   {
              //     // From a quadratic fit to 3 cells
              //     const Real d2W_im1 = 2.*W_im1 - 5*W_i + 4.*W_ip1 - W_ip2;
              //     const Real d2W_i   = W_im1 - 2*W_i   + W_ip1;
              //     const Real d2W_ip1 = W_i - 2.*W_ip1 + W_ip2;
              //     //const Real d2W_ip1 = 2*W_ip1 - 5*W_i + 4*W_im1 - W_im2;
              //     if (d2W_im1*d2W_i < 0. || d2W_i*d2W_ip1 < 0.)
              //       {
              //         // From a cubic fit to four cells
              //         const Real d2W_im1h = 0.5*(3.*W_im1 - 7.*W_i
              //                                    + 5.*W_ip1 - W_ip2);
              //         WfaceL = d1Func(W_i, W_im1, d2W_im1h);
              //       }
              //   }
              // FIXME: Cannot seem to apply limiting or 5th-order and
              // remain stable
              a_WfaceMinus[MD_IX(i, c)] = WfaceL3;

//--Right face (Wplus)

              Real WfaceR = (2.*W_ip3 - 13.*W_ip2 + 47.*W_ip1
                             + 27.*W_i - 3.*W_im1)*fact;
              Real WfaceR3 = (-W_im1 + 5*W_i + 2*W_ip1)/6.0;
              if (std::abs(WfaceR - WfaceR3) > tol &&
                  (WfaceR - W_i)*(W_ip1 - WfaceR) < 0.)
                {
                  // From a quadratic fit to 3 cells
                  const Real d2W_im1 = 2*W_im1 - 15*W_i + 4*W_ip1 - W_ip2;
                  const Real d2W_i   = W_im1 - 2*W_i   + W_ip1;
                  const Real d2W_ip1 = W_i   - 2*W_ip1 + W_ip2;
                  const Real d2W_ip2 = 2*W_ip2 - 5*W_ip1 + 4*W_i - W_im1;
                  if (d2W_im1*d2W_i   < 0. ||
                      d2W_i  *d2W_ip1 < 0. ||
                      d2W_ip1*d2W_ip2 < 0.)
                    {
                      // From a cubic fit to four cells
                      const Real d2W_ip1h = 0.5*(W_im1 - W_i - W_ip1 + W_ip2);
                      WfaceR = d1Func(W_ip1, W_i, d2W_ip1h);
                    }
                }
              a_WfacePlus[MD_IX(i, c)] = WfaceR;
            }
        }

//--Interior cells

      MD_BOXLOOP(cenBox, i)
        {
          Real W_im2 = a_WcellAvg[MD_OFFSETIX(i,-,2*o, c)];
          Real W_im1 = a_WcellAvg[MD_OFFSETIX(i,-,1*o, c)];
          Real W_i   = a_WcellAvg[MD_IX(i, c)];
          Real W_ip1 = a_WcellAvg[MD_OFFSETIX(i,+,1*o, c)];
          Real W_ip2 = a_WcellAvg[MD_OFFSETIX(i,+,2*o, c)];

          const Real tol = std::abs(W_i)*eps;

//--Left face (Wminus)

          Real WfaceL = (2*W_ip2 - 13*W_ip1 + 47*W_i + 27*W_im1 - 3*W_im2)*fact;
          Real WfaceL3 = (-W_ip1 + 5*W_i + 2*W_im1)/6.0;
          if (std::abs(WfaceL - WfaceL3) > tol &&
              (WfaceL - W_im1)*(W_i - WfaceL) < 0.)
            {
              // From a quadratic fit to 3 cells
              const Real d2W_im2 = 2*W_im2 - 5*W_im1 + 4*W_i - W_ip1;
              const Real d2W_im1 = W_im2 - 2*W_im1 + W_i;
              const Real d2W_i   = W_im1 - 2*W_i   + W_ip1;
              const Real d2W_ip1 = 2*W_ip1 - 5*W_i + 4*W_im1 - W_im2;
              if (d2W_im2*d2W_im1 < 0. ||
                  d2W_im1*d2W_i   < 0. ||
                  d2W_i  *d2W_ip1 < 0.)
                {
                  // From a cubic fit to four cells
                  const Real d2W_im1h = 0.5*(W_im2 - W_im1 - W_i + W_ip1);
                  WfaceL = d1Func(W_i, W_im1, d2W_im1h);
                }
            }
          a_WfaceMinus[MD_IX(i, c)] = WfaceL;

//--Right face (Wplus)

          Real WfaceR = (2*W_im2 - 13*W_im1 + 47*W_i + 27*W_ip1 - 3*W_ip2)*fact;
          Real WfaceR3 = (-W_im1 + 5*W_i + 2*W_ip1)/6.0;
          if (std::abs(WfaceR - WfaceR3) > tol &&
              (WfaceR - W_i)*(W_ip1 - WfaceR) < 0.)
            {
              // From a quadratic fit to 3 cells
              const Real d2W_im1 = 2*W_im1 - 5*W_i + 4*W_ip1 - W_ip2;
              const Real d2W_i   = W_im1 - 2*W_i   + W_ip1;
              const Real d2W_ip1 = W_i   - 2*W_ip1 + W_ip2;
              const Real d2W_ip2 = 2*W_ip2 - 5*W_ip1 + 4*W_i - W_im1;
              if (d2W_im1*d2W_i   < 0. ||
                  d2W_i  *d2W_ip1 < 0. ||
                  d2W_ip1*d2W_ip2 < 0.)
                {
                  // From a cubic fit to four cells
                  const Real d2W_ip1h = 0.5*(W_im1 - W_i - W_ip1 + W_ip2);
                  WfaceR = d1Func(W_ip1, W_i, d2W_ip1h);
                }
            }
          a_WfacePlus[MD_IX(i, c)] = WfaceR;
        }

//--Upper boundary adjacent cells

      if (hasHiCells)
        {
          // FIXME: Find proper limiting for boundary faces
          MD_BOXLOOP(farHiCells, i)
            {
              Real W_i   = a_WcellAvg[MD_IX(i, c)];
              Real W_im1 = a_WcellAvg[MD_OFFSETIX(i,-,1*o, c)];
              Real W_im2 = a_WcellAvg[MD_OFFSETIX(i,-,2*o, c)];
              Real W_im3 = a_WcellAvg[MD_OFFSETIX(i,-,3*o, c)];
              Real W_im4 = a_WcellAvg[MD_OFFSETIX(i,-,4*o, c)];
              Real vals[4] = {W_i, W_im1, W_im2, W_im3};

              Real WfaceL = (-3.*W_im4 + 17.*W_im3 - 43.*W_im2
                             + 77.*W_im1 + 12.*W_i)*fact;
              Real WfaceL3 = (-W_im2 + 5.*W_im1 + 2.*W_i)/6.;

              Real WfaceR = (12.*W_im4 - 63.*W_im3 - 163.*W_im1
                             + 137.*(W_im2 + W_i))*fact;
              Real WfaceR3 = (2.*W_im2 - 7.*W_im1 + 11.*W_i)/6.;

              Real avg = 0.;
              for (int tc = 0; tc != 4; ++tc)
                {
                  avg += vals[tc];
                }
              avg /= 4.;

              Real varianceR = std::abs(avg - WfaceR);
              Real varianceR3 = std::abs(avg - WfaceR3);
              Real varianceL = std::abs(avg - WfaceL);
              Real varianceL3 = std::abs(avg - WfaceL3);

              if (varianceL < varianceL3)
                {
                  a_WfaceMinus[MD_IX(i, c)] = WfaceL;
                }
              else
                {
                  a_WfaceMinus[MD_IX(i, c)] = WfaceL3;
                }

              if (varianceR < varianceR3)
                {
                  a_WfacePlus[MD_IX(i, c)] = WfaceR;
                }
              else
                {
                  a_WfacePlus[MD_IX(i, c)] = WfaceR3;
                }
            }
          MD_BOXLOOP(nearHiCells, i)
            {
              Real W_ip1 = a_WcellAvg[MD_OFFSETIX(i,+,1*o, c)];
              Real W_i   = a_WcellAvg[MD_IX(i, c)];
              Real W_im1 = a_WcellAvg[MD_OFFSETIX(i,-,1*o, c)];
              Real W_im2 = a_WcellAvg[MD_OFFSETIX(i,-,2*o, c)];
              Real W_im3 = a_WcellAvg[MD_OFFSETIX(i,-,3*o, c)];

              const Real tol = eps*std::abs(W_i);

//--Left face (Wminus)

              Real WfaceL = (2.*W_im3 - 13.*W_im2 + 47.*W_im1
                             + 27.*W_i - 3.*W_ip1)*fact;
              Real WfaceL3 = (-W_ip1 + 5*W_i + 2*W_im1)/6.0;
              if (std::abs(WfaceL - WfaceL3) > tol &&
                  (WfaceL - W_im1)*(W_i - WfaceL) < 0.)
                {
                  // From a quadratic fit to 3 cells
                  const Real d2W_im2 = 2*W_im2 - 5*W_im1 + 4*W_i - W_ip1;
                  const Real d2W_im1 = W_im2 - 2*W_im1 + W_i;
                  const Real d2W_i   = W_im1 - 2*W_i   + W_ip1;
                  const Real d2W_ip1 = 2*W_ip1 - 5*W_i + 4*W_im1 - W_im2;
                  if (d2W_im2*d2W_im1 < 0. ||
                      d2W_im1*d2W_i   < 0. ||
                      d2W_i  *d2W_ip1 < 0.)
                    {
                      // From a cubic fit to four cells
                      const Real d2W_im1h = 0.5*(W_im2 - W_im1 - W_i + W_ip1);
                      WfaceL = d1Func(W_i, W_im1, d2W_im1h);
                    }
                }
              a_WfaceMinus[MD_IX(i, c)] = WfaceL;

//--Right face (Wplus)
              
              // Real WfaceR = (-3.*W_im3 + 17.*W_im2 - 43.*W_im1
              //                + 77.*W_i + 12.*W_ip1)*fact;
              Real WfaceR3 = (-W_im1 + 5*W_i + 2*W_ip1)/6.0;
              // if (std::abs(WfaceR - WfaceR3) > tol &&
              //     (WfaceR - W_i)*(W_ip1 - WfaceR) < 0.)
              //   {
              //     // From a quadratic fit to 3 cells
              //     const Real d2W_im1 = W_im2 - 2.*W_im1 + W_i;
              //     const Real d2W_i   = W_im1 - 2*W_i + W_ip1;
              //     const Real d2W_ip1 = 2.*W_ip1 - 5.*W_i + 4.*W_im1 - W_im2;
              //     //const Real d2W_ip2 = 2*W_ip2 - 5*W_ip1 + 4*W_i - W_im1;
              //     if (d2W_im1*d2W_i   < 0. || d2W_i*d2W_ip1 < 0.)
              //       {
              //         // From a cubic fit to four cells
              //         const Real d2W_ip1h = 0.5*(3.*W_ip1 - 7.*W_i
              //                                    + 5.*W_im1 - W_im2);
              //         WfaceR = d1Func(W_ip1, W_i, d2W_ip1h);
              //       }
              //   }
              // FIXME: Cannot seem to apply limiting or 5th-order and
              // remain stable
              a_WfacePlus[MD_IX(i, c)] = WfaceR3;
            }
        }
    }
}

void
CRDutil::getSecondDiff(FArrayBox&           a_d2Fab,
                       const FArrayBox&     a_vals,
                       const Box&           a_box,
                       const int            a_numComp,
                       const ProblemDomain& a_domain)
{
  CH_TIME("CRDutil::getSecondDiff");
  int hasLo, hasHi;
  Box loBox, hiBox, centerBox, entireBox;
  a_d2Fab.setVal(0.);
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      Box inBox(a_box);
      inBox.grow(dir,1);
      loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                 inBox, a_domain, dir);
      CH_assert(a_vals.contains(entireBox));
      FABSTACKTEMP(d2Dir, entireBox, a_numComp);
      FORT_GETSECONDDIFF(CHF_FRA(d2Dir),
                         CHF_CONST_FRA(a_vals),
                         CHF_CONST_INT(a_numComp),
                         CHF_CONST_INT(dir),
                         CHF_BOX(loBox),
                         CHF_CONST_INT(hasLo),
                         CHF_BOX(hiBox),
                         CHF_CONST_INT(hasHi),
                         CHF_BOX(centerBox));
      a_d2Fab.plus(d2Dir);
    }
}

// void
// CRDutil::convolve(FArrayBox& a_avgFab,
//                   const FArrayBox& a_pntFab,
//                   const Box& a_box,
//                   const int a_order,
//                   bool a_limit,
//                   int a_sign)
// {
//   MD_ARRAY_RESTRICT(arrAvgFab, a_avgFab);
//   MD_ARRAY_RESTRICT(arrPntFab, a_pntFab);
// }

inline static Real
limCen1stDiff2o(const Real W_im1,
                const Real W_i,
                const Real W_ip1) noexcept
{
  const Real d1L = W_i - W_im1;
  const Real d1R = W_ip1 - W_i;
  Real d1 = 0.5*(d1L + d1R);
  if (d1L*d1R < 0.)
    {
      d1 = 0.;
    }
  else if (d1L > 2*d1R)  // d1(3h/2) = 0
    {
      d1 = (3./2.)*d1R;
    }
  else if (d1R > 2*d1L)  // d1(-3h/2) = 0
    {
      d1 = (3./2.)*d1L;
    }
  return d1;
}

inline static Real
limFwd1stDiff2o(const Real W_im2,
                const Real W_im1,
                const Real W_i,
                const int  sgn) noexcept
{
  const Real d1L = W_im1 - W_im2;
  const Real d1R = W_i - W_im1;
  Real d1 = sgn*0.5*(3*W_i - 4*W_im1 + W_im2);
  if (d1L*d1R < 0.)
    {
      d1 = 0.;
    }
  else if (d1L > 2*d1R)  // d1(h/2) = 0
    {
      d1 = sgn*(0.5*d1R);
    }
  else if (d1R > 2*d1L)  // d1(-5*h/2) = 0
    {
      d1 = sgn*(1.25*d1R);
    }
  return d1;
}


/// Get the average kinetic energy in a cell using a triple product rule
/** This returns the average kinetic energy per unit volume
 */
void
CRDutil::avgKE(FArrayBox&           a_avgKEfab,
               const FArrayBox&     a_UavgFab,
               const Box&           a_box,
               const ProblemDomain& a_domain,
               bool                 a_limit)
{
  CH_assert(a_domain.contains(a_box));
  CH_assert(a_avgKEfab.contains(a_box));
  const Box box1 = grow(a_box, 1) & a_domain;  // Touched a_UavgFab cells
  CH_assert(a_UavgFab.contains(box1));

  constexpr int cRho   = CRDPhysics::densityIndex();
  constexpr int cMom0  = CRDPhysics::velocityInterval().begin();
  constexpr int cMomSD = CRDPhysics::velocityInterval().end();  // SpaceDim
  constexpr Real factor12 = 1./12.;

  a_avgKEfab.setVal(0.);

  struct PointDataKEOp
  {
    Real rho_i;
    Real mom_i;
    Real rhoInv;
  };
  PointDataKEOp pdataKEOp;

  for (int cMom = cMom0, cMom_end = cMomSD + 1; cMom != cMom_end; ++cMom)
    {
      // Pre-stencil point operation
      auto prePointKEOp =
        [&, &pdata = pdataKEOp]
        (MD_DECLIX(const int, i))
        {
          pdata.rho_i  = a_UavgFab[MD_IX(i, cRho)];
          pdata.mom_i  = a_UavgFab[MD_IX(i, cMom)];
          pdata.rhoInv = 1./pdata.rho_i;
        };

      // Product term for centered stencil
      auto centerStencilKEOp =
        [&, &pdata = pdataKEOp]
        (MD_DECLIX(const int, i),
         const int            dir,
         MD_DECLIX(const int, ii))
        {
          Real d1Rho = limCen1stDiff2o(a_UavgFab[MD_OFFSETIX(i,-,ii, cRho)],
                                       pdata.rho_i,
                                       a_UavgFab[MD_OFFSETIX(i,+,ii, cRho)]);
          Real d1Mom = limCen1stDiff2o(a_UavgFab[MD_OFFSETIX(i,-,ii, cMom)],
                                       pdata.mom_i,
                                       a_UavgFab[MD_OFFSETIX(i,+,ii, cMom)]);
          return d1Mom*(d1Mom - 2*pdata.rhoInv*pdata.mom_i*d1Rho);
        };

      // Product term for one-sided stencil
      auto offsetStencilKEOp =
        [&, &pdata = pdataKEOp]
        (const int            sgn,
         MD_DECLIX(const int, i),
         const int            dir,
         MD_DECLIX(const int, ii))
        {
          Real d1Rho = limFwd1stDiff2o(
            a_UavgFab[MD_OFFSETIX(i,-,2*sgn*ii, cRho)],
            a_UavgFab[MD_OFFSETIX(i,-,  sgn*ii, cRho)],
            pdata.rho_i, sgn);
          Real d1Mom = limFwd1stDiff2o(
            a_UavgFab[MD_OFFSETIX(i,-,2*sgn*ii, cMom)],
            a_UavgFab[MD_OFFSETIX(i,-,  sgn*ii, cMom)],
            pdata.mom_i, sgn);
          return d1Mom*(d1Mom - 2*pdata.rhoInv*pdata.mom_i*d1Rho);
        };

      // Post-stencil point operation for product rule
      auto postPointKEOp =
        [&, &pdata = pdataKEOp]
        (MD_DECLIX(const int, i),
         const Real           stencilTerm)
        {
          a_avgKEfab[MD_IX(i, 0)] += 0.5*(pdata.mom_i*pdata.mom_i*pdata.rhoInv +
                                          factor12*pdata.rhoInv*stencilTerm);
        };

      loHiCenter3ApplyOp<Real>(a_box,
                               a_domain,
                               prePointKEOp,
                               centerStencilKEOp,
                               offsetStencilKEOp,
                               postPointKEOp);
    }  // Loop over momentum components
}

/// Get the average internal energy from the kinetic energy
/** This returns the internal energy per unit mass
 *  \param[in]  a_avgE  Average kinetic energy in a cell
 *  \param[out] a_avgE  Average internal energy in a cell
 */
void
CRDutil::avgE(FArrayBox&           a_avgEfab,
              const FArrayBox&     a_barEfab,
              const FArrayBox&     a_UavgFab,
              const Box&           a_box,
              const ProblemDomain& a_domain,
              bool                 a_limit)
{
  CH_assert(a_domain.contains(a_box));
  CH_assert(a_avgEfab.contains(a_box));
  const Box box1 = grow(a_box, 1) & a_domain;  // Touched Uavg and barE cells
  CH_assert(a_barEfab.contains(box1));
  CH_assert(a_UavgFab.contains(box1));

  constexpr int cRho   = CRDPhysics::densityIndex();
  constexpr int cRhoET = CRDPhysics::energyFluxIndex();
  constexpr Real factor12 = 1./12.;

//--Determine <e> in each cell

  struct PointDataEOp
  {
    Real rho_i;
  };
  PointDataEOp pdataEOp;

  // Pre-stencil point operation for product rule
  auto prePointEOp =
    [&, &pdata = pdataEOp]
    (MD_DECLIX(const int, i))
  {
    pdata.rho_i = a_UavgFab[MD_IX(i, cRho)];
  };

  // Centered stencil for product term
  auto centerStencilEOp =
    [&, &pdata = pdataEOp]
    (MD_DECLIX(const int, i),
     const int            dir,
     MD_DECLIX(const int, ii))
  {
    const Real rho_im1 = a_UavgFab[MD_OFFSETIX(i,-,ii, cRho)];
    const Real rho_ip1 = a_UavgFab[MD_OFFSETIX(i,+,ii, cRho)];
          
    Real d1Rho  = limCen1stDiff2o(rho_im1, pdata.rho_i, rho_ip1);
    Real d1RhoE = limCen1stDiff2o(a_barEfab[MD_OFFSETIX(i,-,ii, 0)]*rho_im1,
                                  a_barEfab[MD_IX(i, 0)]           *pdata.rho_i,
                                  a_barEfab[MD_OFFSETIX(i,+,ii, 0)]*rho_ip1);
    return d1Rho*d1RhoE;
  };

  // One-sided stencil for product term
  auto offsetStencilEOp =
    [&, &pdata = pdataEOp]
    (const int            sgn,
     MD_DECLIX(const int, i),
     const int            dir,
     MD_DECLIX(const int, ii))
  {
    const Real rho_im2 = a_UavgFab[MD_OFFSETIX(i,-,2*sgn*ii, cRho)];
    const Real rho_im1 = a_UavgFab[MD_OFFSETIX(i,-,  sgn*ii, cRho)];
          
    Real d1Rho  = limFwd1stDiff2o(rho_im2, rho_im1, pdata.rho_i, sgn);
    Real d1RhoE = limFwd1stDiff2o(
      a_barEfab[MD_OFFSETIX(i,-,2*sgn*ii, 0)]*rho_im2,
      a_barEfab[MD_OFFSETIX(i,-,  sgn*ii, 0)]*rho_im1,
      a_barEfab[MD_IX(i, 0)]*pdata.rho_i, sgn);
    return d1Rho*d1RhoE;
  };

  // Post-stencil point operation for product rule
  auto postPointEOp =
    [&, &pdata = pdataEOp]
    (MD_DECLIX(const int, i),
     const Real           stencilTerm)
  {
    const Real rhoE  = a_UavgFab[MD_IX(i, cRhoET)] - a_avgEfab[MD_IX(i, 0)];
    const Real rhoInv = 1./pdata.rho_i;
    a_avgEfab[MD_IX(i, 0)] = rhoInv*(rhoE - factor12*rhoInv*stencilTerm);
  };

  loHiCenter3ApplyOp<Real>(a_box,
                           a_domain,
                           prePointEOp,
                           centerStencilEOp,
                           offsetStencilEOp,
                           postPointEOp);
}

void
CRDutil::limitTpnt(const FArrayBox&     a_UavgFab,
                   const FArrayBox&     a_WpntFab,
                   const Box&           a_box,
                   const ProblemDomain& a_domain,
                   bool                 a_limit)
{
  CH_assert(a_domain.contains(a_box));
  CH_assert(a_WpntFab.contains(a_box));
  constexpr int cRho   = CRDPhysics::densityIndex();
  constexpr int cRhoET = CRDPhysics::energyFluxIndex();
  constexpr int cMom0  = CRDPhysics::velocityInterval().begin();

  const Box box1 = grow(a_box, 1) & a_domain;
  const Box box2 = grow(a_box, 2) & a_domain;
  CH_assert(a_UavgFab.contains(box2));

//--Get the average kinetic energy in each cell

  FABSTACKTEMP(avgEfab, box1, 1);
  avgKE(avgEfab, a_UavgFab, box1, a_domain, a_limit);

//--Determine \bar{e} in each cell.  This is a second-order approximation of
//--<e> used for derivatives.  Note: multiply by <rho> to get <rho e> when
//--needed.

  FABSTACKTEMP(barEfab, box2, 1);
  MD_BOXLOOP(box2, i)
    {
      const Real rhoInv = 1./a_UavgFab[MD_IX(i, cRho)];
      const RealVect mom{D_DECL(a_UavgFab[MD_IX(i, cMom0 + 0)],
                                a_UavgFab[MD_IX(i, cMom0 + 1)],
                                a_UavgFab[MD_IX(i, cMom0 + 2)])};
      const Real rhoET = a_UavgFab[MD_IX(i, cRhoET)];
      barEfab[MD_IX(i, 0)] = (rhoET - 0.5*stc::sum(mom*mom)*rhoInv)*rhoInv;
    }

//--The the average intenal energy per unit mass

  avgE(avgEfab, barEfab, a_UavgFab, box1, a_domain, a_limit);

//--Determine \bar{T} in each cell

  FABSTACKTEMP(TbarFab, box1, 1);
  FABSTACKTEMP(dFdTbarFab, box1, 1);
  CRDparam::g_CRDPhysics->intermediateTandDiffF(TbarFab,
                                                dFdTbarFab,
                                                a_WpntFab,
                                                avgEfab,
                                                box1);

//--Determine T at cell centers

  // Get the Laplacian of \bar{e} ...
  FABSTACKTEMP(D2TbarFab, box1, 1);
  Laplacian(D2TbarFab,
            barEfab,
            box1,
            a_domain,
            Interval(0, 0),
            0);
  // ... and convert to a Laplacian of \bar{T}
  MD_BOXLOOP(a_box, i)
    {
      // The internal energy should always increase with T so dF/dT = de/dT
      // should always be greater than 0.
      CH_assert(dFdTbarFab[MD_IX(i, 0)] > 0.);
      D2TbarFab[MD_IX(i, 0)] /= dFdTbarFab[MD_IX(i, 0)];
    }

  // Get local min/max of \bar{T}
  FABSTACKTEMP(lclMinTbarFab, a_box, 1);
  FABSTACKTEMP(lclMaxTbarFab, a_box, 1);
  localBounds(lclMinTbarFab, lclMaxTbarFab, TbarFab, a_box, a_domain,
              Interval(0, 0), 0);

  constexpr Real factor24 = 1./24.;
  const int cT = CRDparam::g_CRDPhysics->temperatureIndex();
  MD_BOXLOOP(a_box, i)
    {
      Real D2T     = D2TbarFab[MD_IX(i, 0)];
      Real Tbar    = TbarFab[MD_IX(i, 0)];
      Real Tpnt    = Tbar - factor24*D2T;
      Real TbarMin = lclMinTbarFab[MD_IX(i, 0)];
      Real TbarMax = lclMaxTbarFab[MD_IX(i, 0)];
      if ((Tpnt - TbarMin)*(Tpnt - TbarMax) > 0.)
        // We are at an extrema of T.  Limit Laplacian instead.
        {
          Real lclD2Tmin;
          Real lclD2Tmax;
          localBounds(lclD2Tmin, lclD2Tmax, D2TbarFab, MD_GETIV(i), a_domain,
                      0);
          if (lclD2Tmin*lclD2Tmax <= 0.)
            {
              TbarMin = Tbar;
              TbarMax = Tbar;
            }
          // Otherwise d2T, lclD2Tmin, lclD2Tmax have the same sign
          else if (D2T < 0.)
            // The extremum is a maximum.  Allowed difference is closest
            // to zero which is maximum of negative D2T.
            {
              TbarMax  = Tbar - factor24*lclD2Tmax;
            }
          else if (D2T > 0.)
            // The extremum is a minimum.  Allowed difference is closest
            // to zero which is minimum of positive D2T.
            {
              TbarMin  = Tbar - factor24*lclD2Tmin;
            }
        }
      a_WpntFab[MD_IX(i, cT)] = 
        std::min(TbarMax, std::max(TbarMin, Tpnt));
    }
}

void
CRDutil::Laplacian(FArrayBox&           a_D2fab,
                   const FArrayBox&     a_D0fab,
                   const Box&           a_box,
                   const ProblemDomain& a_domain,
                   const Interval&      a_intv,
                   bool                 a_cSto)
{
  CH_TIME("CRDutil::Laplacian");
  CRD::msg << CRD::fv4 << "CRDutil::Laplacian" << CRD::end;
  CH_assert(a_domain.contains(a_box));
  CH_assert(a_D2fab.contains(a_box));
  const Box box1 = grow(a_box, 1) & a_domain;  // Touched D0fab cells
  CH_assert(a_D0fab.contains(box1));

  // Pre-stencil point operation
  auto prePointD2 =
    []
    (MD_DECLIX(const int, i))
  { };

  for (int iComp = a_intv.begin(), iComp_end = a_intv.end() + 1, cSto = a_cSto;
       iComp != iComp_end; ++iComp, ++cSto)
    {
      // Centered stencil for Laplacian
      auto centerStencilD2 =
        [&, iComp]
        (MD_DECLIX(const int, i),
         const int            dir,
         MD_DECLIX(const int, ii))
      {
        return     a_D0fab[MD_OFFSETIX(i,-,ii, iComp)]
               - 2*a_D0fab[MD_IX(i, iComp)]
               +   a_D0fab[MD_OFFSETIX(i,+,ii, iComp)];
      };

      // Offset stencil for Laplacian
      auto offsetStencilD2 =
        [&, iComp]
        (const int            sgn,
         MD_DECLIX(const int, i),
         const int            dir,
         MD_DECLIX(const int, ii))
      {
        return 2*a_D0fab[MD_IX(i, iComp)]
              -5*a_D0fab[MD_OFFSETIX(i,-,  sgn*ii, iComp)]
              +4*a_D0fab[MD_OFFSETIX(i,-,2*sgn*ii, iComp)]
              -1*a_D0fab[MD_OFFSETIX(i,-,3*sgn*ii, iComp)];
      };

      // Post-stencil point operation for Laplacian
      auto postPointD2 =
        [&, cSto]
        (MD_DECLIX(const int, i),
         const Real           stencilTerm)
      {
        a_D2fab[MD_IX(i, cSto)] = stencilTerm;
      };

      loHiCenter3ApplyOp<Real>(a_box,
                               a_domain,
                               prePointD2,
                               centerStencilD2,
                               offsetStencilD2,
                               postPointD2);
    }  // Loop over components
}

/*--------------------------------------------------------------------*/
//  Compute cell-centered values from cell-averaged values
/** \param[out] a_pntFab
 *                      Fab of cell-centered values
 *  \param[out] a_avgFab
 *                      Fab of cell-averaged values
 *  \param[in]  a_box   Box over which to compute cell-centered values
 *  \param[in]  a_domain
 *                      Problem domain for checking boundaries
 *  \param[in]  a_intv  Interval of components for computing
 *                      cell-centered values
 *  \param[in]  a_order Order of accuracy of deconvolution
 *  \param[in]  a_sign  Whether this is deconvolution (+) or 
 *                      convolution (-)
 *  \param[in]  a_limit Limit the deconvolution (if true)
 *  \param[in]  a_preCopy
 *                      a_pntFab only needs the laplacian subtracted
 *                      (if true)
 *//*-----------------------------------------------------------------*/

void
CRDutil::deconvolve(FArrayBox&           a_pntFab,
                    const FArrayBox&     a_avgFab,
                    const Box&           a_box,
                    const ProblemDomain& a_domain,
                    const Interval&      a_intv,
                    const int            a_order,
                    const int            a_sign,
                    bool                 a_limit,
                    bool                 a_preCopy)
{
  CH_TIME("CRDutil::deconvolve");
  CRD::msg << CRD::fv4 << "CRDutil::deconvolve" << CRD::end;
  CH_assert(a_pntFab.contains(a_box));
  CH_assert(a_pntFab.interval().contains(a_intv));
  CH_assert(a_avgFab.interval().contains(a_intv));
  CH_assert(a_order == 2 || a_order == 4);
  constexpr Real factor24 = 1./24.;

  if (a_order == 2)
    {
      CH_assert(a_avgFab.contains(a_box));
      a_pntFab.copy(a_avgFab, a_box, a_intv.begin(), a_box, a_intv.begin(),
                    a_intv.size());
      return;
    }

  Box avgBox = grow(a_box, 1 + (int)a_limit);
  avgBox &= a_domain;
  CH_assert(a_avgFab.contains(avgBox));

  Box D2box = grow(a_box, (int)a_limit);  // Box where we need the Laplacian
  D2box &= a_domain;

  FABSTACKTEMP(D2avgFab, D2box, 1);
  FABSTACKTEMP(lclMinAvgFab, a_box, 1);
  FABSTACKTEMP(lclMaxAvgFab, a_box, 1);
  for (int iComp = a_intv.begin(), iComp_end = a_intv.end() + 1;
       iComp != iComp_end; ++iComp)
    {
      // Get the Laplacian
      Laplacian(D2avgFab,
                a_avgFab,
                D2box,
                a_domain,
                Interval(iComp, iComp),
                0);
      if (a_limit && a_preCopy)
        {
          // Get local min/max of <U>
          localBounds(lclMinAvgFab, lclMaxAvgFab, a_avgFab, a_box, a_domain,
                      Interval(iComp, iComp), 0);
          MD_BOXLOOP(a_box, i)
            {
              Real D2U     = D2avgFab[MD_IX(i, 0)];
              Real Uavg    = a_avgFab[MD_IX(i, iComp)];
              Real Upnt    = a_pntFab[MD_IX(i, iComp)] - a_sign*factor24*D2U;
              Real UavgMin = lclMinAvgFab[MD_IX(i, 0)];
              Real UavgMax = lclMaxAvgFab[MD_IX(i, 0)];
              if ((Upnt - UavgMin)*(Upnt - UavgMax) > 0.)
                // We are at an extrema of U.  Limit Laplacian instead.
                {
                  Real lclD2Umin;
                  Real lclD2Umax;
                  localBounds(lclD2Umin, lclD2Umax, D2avgFab, MD_GETIV(i),
                              a_domain, 0);
                  if (lclD2Umin*lclD2Umax <= 0.)
                    {
                      UavgMin = Uavg;
                      UavgMax = Uavg;
                    }
                  // Otherwise d2U, lclD2Umin, lclD2Umax have the same sign
                  else if (D2U < 0.)
                    // The extremum is a maximum.  Allowed difference is closest
                    // to zero which is maximum of negative D2U.
                    {
                      UavgMax = Uavg - a_sign*factor24*lclD2Umax;
                    }
                  else if (D2U > 0.)
                    // The extremum is a minimum.  Allowed difference is closest
                    // to zero which is minimum of positive D2U.
                    {
                      UavgMin = Uavg - a_sign*factor24*lclD2Umin;
                    }
                }
              a_pntFab[MD_IX(i, iComp)] =
                std::min(UavgMax, std::max(UavgMin, Upnt));
            }
        }
      else if (a_limit)
        {
          // Get local min/max of <U>
          localBounds(lclMinAvgFab, lclMaxAvgFab, a_avgFab, a_box, a_domain,
                      Interval(iComp, iComp), 0);
          MD_BOXLOOP(a_box, i)
            {
              Real D2U     = D2avgFab[MD_IX(i, 0)];
              Real Uavg    = a_avgFab[MD_IX(i, iComp)];
              Real Upnt    = Uavg - a_sign*factor24*D2U;
              Real UavgMin = lclMinAvgFab[MD_IX(i, 0)];
              Real UavgMax = lclMaxAvgFab[MD_IX(i, 0)];
              if ((Upnt - UavgMin)*(Upnt - UavgMax) > 0.)
                // We are at an extrema of U.  Limit Laplacian instead.
                {
                  Real lclD2Umin;
                  Real lclD2Umax;
                  localBounds(lclD2Umin, lclD2Umax, D2avgFab, MD_GETIV(i),
                              a_domain, 0);
                  if (lclD2Umin*lclD2Umax <= 0.)
                    {
                      UavgMin = Uavg;
                      UavgMax = Uavg;
                    }
                  // Otherwise d2U, lclD2Umin, lclD2Umax have the same sign
                  else if (D2U < 0.)
                    // The extremum is a maximum.  Allowed difference is closest
                    // to zero which is maximum of negative D2U.
                    {
                      UavgMax = Uavg - a_sign*factor24*lclD2Umax;
                    }
                  else if (D2U > 0.)
                    // The extremum is a minimum.  Allowed difference is closest
                    // to zero which is minimum of positive D2U.
                    {
                      UavgMin = Uavg - a_sign*factor24*lclD2Umin;
                    }
                }
              a_pntFab[MD_IX(i, iComp)] =
                std::min(UavgMax, std::max(UavgMin, Upnt));
            }
        }
      else if (a_preCopy)  // No limiting
        {
          MD_BOXLOOP(a_box, i)
            {
              a_pntFab[MD_IX(i, iComp)] -=
                a_sign*factor24*D2avgFab[MD_IX(i, 0)];
            }
        }
      else // No limiting and no prefilling of a_pntFab
        {
          MD_BOXLOOP(a_box, i)
            {
              a_pntFab[MD_IX(i, iComp)] = a_avgFab[MD_IX(i, iComp)]
                - a_sign*factor24*D2avgFab[MD_IX(i, 0)];
            }
        }
    }  // Loop over components
}

/*--------------------------------------------------------------------*/
//  Compute face-centered values from face-averaged values
/** \param[out] a_pntFab
 *                      Fab of face-centered values
 *  \param[out] a_avgFab
 *                      Fab of face-averaged values
 *  \param[in]  a_box   Box over which to compute face-centered values
 *  \param[in]  a_domain
 *                      Problem domain for checking boundaries
 *  \param[in]  a_intv  Interval of components for computing
 *                      face-centered values
 *  \param[in]  a_dir   Face direction of a_pntFab
 *  \param[in]  a_order Order of accuracy of deconvolution
 *  \param[in]  a_limit Limit the deconvolution (if true)
 *  \param[in]  a_preCopy
 *                      a_pntFab only needs the laplacian subtracted
 *                      (if true)
 *//*-----------------------------------------------------------------*/

void
CRDutil::deconvolveFace(FArrayBox&           a_pntFab,
                        const FArrayBox&     a_avgFab,
                        const Box&           a_box,
                        const ProblemDomain& a_domain,
                        const Interval&      a_intv,
                        const int            a_dir,
                        const int            a_order,
                        bool                 a_limit,
                        bool                 a_preCopy)
{
  CH_TIME("CRDutil::deconvolveFace");
  CRD::msg << CRD::fv4 << "CRDutil::deconvolveFace" << CRD::end;
  //**FIXME: a lot of things need to be fixed
  //         1) add CH_asserts for box sizes and component counts
  //         2) need to change this to use a face-laplacian operator

  // We require that we're given a face-box in direction a_dir
  IntVect testBoxType = a_box.type();
  IntVect cellBoxType(IntVect::Zero);
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      if ((dir == a_dir) &&
          (testBoxType[dir] == cellBoxType[dir]))
        {
          CH_assert(false);
        }
      else if ((dir != a_dir) &&
               (testBoxType[dir] != cellBoxType[dir])) // not allowed
        {
          CH_assert(false);
        }
    }

  // Check that a_box is inside a_domain
  CH_assert(a_domain.contains(a_box));

  constexpr Real factor24 = 1./24.;

  if (a_order == 2)
    {
      if (a_preCopy)
        {
          return;
        }
      a_pntFab.copy(a_avgFab, a_box, a_intv.begin(), a_box, a_intv.begin(),
                    a_intv.size());
      return;
    }

  int limitStencil = (int)a_limit;
  IntVect growVect(limitStencil*IntVect::Unit);
  growVect[a_dir] = 0;
  Box D2box = grow(a_box, growVect);
  D2box &= a_domain;

  FABSTACKTEMP(D2avgFab, D2box, 1);
  FABSTACKTEMP(lclMinAvgFab, a_box, 1);
  FABSTACKTEMP(lclMaxAvgFab, a_box, 1);
  FABSTACKTEMP(lclD2UminFab, a_box, 1);
  FABSTACKTEMP(lclD2UmaxFab, a_box, 1);
  for (int iComp = a_intv.begin(), iComp_end = a_intv.end() + 1;
       iComp != iComp_end; ++iComp)
    {
      D2avgFab.setVal(0.);
      for (int tdir = 0; tdir != SpaceDim; ++tdir)
        {
          if (tdir != a_dir)
            {
              const int MD_ID(o, tdir);
              Box loBox, hiBox, centerBox, entireBox;
              int hasLo, hasHi;
              Box inputBox = D2box;
              inputBox.grow(tdir, 1);
              loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox,
                         entireBox, inputBox, a_domain, tdir);
              if (hasLo)
                {
                  loBox &= a_domain;
                }
              if (hasHi)
                {
                  hiBox &= a_domain;
                }
              MD_BOXLOOP(centerBox, i)
                {
                  D2avgFab[MD_IX(i, 0)] += (
                       a_avgFab[MD_OFFSETIX(i,-,o, iComp)]
                  - 2.*a_avgFab[MD_IX(i, iComp)]
                  +    a_avgFab[MD_OFFSETIX(i,+,o, iComp)]);
                }
              if (hasLo)
                {
                  MD_BOXLOOP(loBox, i)
                    {
                      D2avgFab[MD_IX(i, 0)] += (
                        2.*a_avgFab[MD_IX(i, iComp)]
                      - 5.*a_avgFab[MD_OFFSETIX(i,+,o, iComp)]
                      + 4.*a_avgFab[MD_OFFSETIX(i,+,2*o, iComp)]
                      -    a_avgFab[MD_OFFSETIX(i,+,3*o, iComp)]);
                    }
                }
              if (hasHi)
                {
                  MD_BOXLOOP(hiBox, i)
                    {
                      D2avgFab[MD_IX(i, 0)] += (
                        2.*a_avgFab[MD_IX(i, iComp)]
                      - 5.*a_avgFab[MD_OFFSETIX(i,-,o, iComp)]
                      + 4.*a_avgFab[MD_OFFSETIX(i,-,2*o, iComp)]
                      -    a_avgFab[MD_OFFSETIX(i,-,3*o, iComp)]);
                    }
                }
            }
        }
      if (a_limit && a_preCopy)
        {
          // Get local min/max of <U>
          localBoundsFace(lclMinAvgFab, lclMaxAvgFab, a_avgFab, a_box, a_domain,
                          Interval(iComp, iComp), 0, a_dir);
          // Get local min/max of D2<U>
          localBoundsFace(lclD2UminFab, lclD2UmaxFab, D2avgFab, a_box, a_domain,
                          Interval(0, 0), 0, a_dir);
          MD_BOXLOOP(a_box, i)
            {
              Real D2U       = D2avgFab[MD_IX(i, 0)];
              Real Uavg      = a_avgFab[MD_IX(i, iComp)];
              Real Upnt      = a_pntFab[MD_IX(i, iComp)] - factor24*D2U;
              Real UavgMin   = lclMinAvgFab[MD_IX(i, 0)];
              Real UavgMax   = lclMaxAvgFab[MD_IX(i, 0)];
              Real lclD2Umin = lclD2UminFab[MD_IX(i, 0)];
              Real lclD2Umax = lclD2UmaxFab[MD_IX(i, 0)];
              if ((Upnt - UavgMin)*(Upnt - UavgMax) > 0.)
                // We are at an extrema of U.  Limit Laplacian instead.
                {
                  if (lclD2Umin*lclD2Umax <= 0.)
                    {
                      UavgMin = Uavg;
                      UavgMax = Uavg;
                    }
                  // Otherwise d2U, lclD2Umin, lclD2Umax have the same sign
                  else if (D2U < 0.)
                    // The extremum is a maximum.  Allowed difference is closest
                    // to zero which is maximum of negative D2U.
                    {
                      UavgMax = a_pntFab[MD_IX(i, iComp)] - factor24*lclD2Umax;
                    }
                  else if (D2U > 0.)
                    // The extremum is a minimum.  Allowed difference is closest
                    // to zero which is minimum of positive D2U.
                    {
                      UavgMin = a_pntFab[MD_IX(i, iComp)] - factor24*lclD2Umin;
                    }
                }
              a_pntFab[MD_IX(i, iComp)] =
                std::min(UavgMax, std::max(UavgMin, Upnt));
            }
        }
      else if (a_limit)
        {
          // Get local min/max of <U>
          localBoundsFace(lclMinAvgFab, lclMaxAvgFab, a_avgFab, a_box, a_domain,
                          Interval(iComp, iComp), 0, a_dir);
          // Get local min/max of D2<U>
          localBoundsFace(lclD2UminFab, lclD2UmaxFab, D2avgFab, a_box, a_domain,
                          Interval(0, 0), 0, a_dir);
          MD_BOXLOOP(a_box, i)
            {
              Real D2U       = D2avgFab[MD_IX(i, 0)];
              Real Uavg      = a_avgFab[MD_IX(i, iComp)];
              Real Upnt      = a_avgFab[MD_IX(i, iComp)] - factor24*D2U;
              Real UavgMin   = lclMinAvgFab[MD_IX(i, 0)];
              Real UavgMax   = lclMaxAvgFab[MD_IX(i, 0)];
              Real lclD2Umin = lclD2UminFab[MD_IX(i, 0)];
              Real lclD2Umax = lclD2UmaxFab[MD_IX(i, 0)];
              if ((Upnt - UavgMin)*(Upnt - UavgMax) > 0.)
                // We are at an extrema of U.  Limit Laplacian instead.
                {
                  if (lclD2Umin*lclD2Umax <= 0.)
                    {
                      UavgMin = Uavg;
                      UavgMax = Uavg;
                    }
                  // Otherwise d2U, lclD2Umin, lclD2Umax have the same sign
                  else if (D2U < 0.)
                    // The extremum is a maximum.  Allowed difference is closest
                    // to zero which is maximum of negative D2U.
                    {
                      UavgMax = a_avgFab[MD_IX(i, iComp)] - factor24*lclD2Umax;
                    }
                  else if (D2U > 0.)
                    // The extremum is a minimum.  Allowed difference is closest
                    // to zero which is minimum of positive D2U.
                    {
                      UavgMin = a_avgFab[MD_IX(i, iComp)] - factor24*lclD2Umin;
                    }
                }
              a_pntFab[MD_IX(i, iComp)] =
                std::min(UavgMax, std::max(UavgMin, Upnt));
            }
        }
      else if (a_preCopy)  // No limiting
        {
          MD_BOXLOOP(a_box, i)
            {
              a_pntFab[MD_IX(i, iComp)] -=
                factor24*D2avgFab[MD_IX(i, 0)];
            }
        }
      else // No limiting and no prefilling of a_pntFab
        {
          MD_BOXLOOP(a_box, i)
            {
              a_pntFab[MD_IX(i, iComp)] = a_avgFab[MD_IX(i, iComp)]
                - factor24*D2avgFab[MD_IX(i, 0)];
            }
        }
    }  // Loop over components
}

/*--------------------------------------------------------------------*/
//  Compute face-centered values from face-averaged values
/** \param[out] a_pntFxb
 *                      Fab of face-centered values
 *  \param[out] a_avgFxb
 *                      Fab of face-averaged values
 *  \param[in]  a_box   Box over which to compute face-centered values
 *  \param[in]  a_domain
 *                      Problem domain for checking boundaries
 *  \param[in]  a_intv  Interval of components for computing
 *                      face-centered values
 *  \param[in]  a_order Order of accuracy of deconvolution
 *  \param[in]  a_limit Limit the deconvolution (if true)
 *  \param[in]  a_preCopy
 *                      a_pntFxb only needs the laplacian subtracted
 *                      (if true)
 *//*-----------------------------------------------------------------*/

void
CRDutil::deconvolveFace(FluxBox&             a_pntFxb,
                        const FluxBox&       a_avgFxb,
                        const Box&           a_box,
                        const ProblemDomain& a_domain,
                        const Interval&      a_intv,
                        const int            a_order,
                        bool                 a_limit,
                        bool                 a_preCopy)
{
  CH_TIME("CRDutil::deconvolveFace");
  CRD::msg << CRD::fv4 << "CRDutil::deconvolveFace" << CRD::end;
  //**FIXME: there's a lot to fix here (CH_asserts)
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      deconvolveFace(a_pntFxb[dir],
                     a_avgFxb[dir],
                     a_box,
                     a_domain,
                     a_intv,
                     dir,
                     a_order,
                     a_limit,
                     a_preCopy);
    }
}

/*--------------------------------------------------------------------*/
//  Compute face-averaged values from face-centered values
/** \param[out] a_avgFab
 *                      Fab of face-averaged values
 *  \param[out] a_pntFab
 *                      Fab of face-centered values
 *  \param[in]  a_box   Box over which to compute face-averaged values
 *  \param[in]  a_domain
 *                      Problem domain for checking boundaries
 *  \param[in]  a_intv  Interval of components for computing
 *                      face-averaged values
 *  \param[in]  a_dir   Face direction of a_avgFab
 *  \param[in]  a_order Order of accuracy of convolution
 *  \param[in]  a_interior
 *                      If false, one-sided convolutions are used at
 *                      domain boundaries
 *  \param[in]  a_limit Limit the convolution (if true)
 *  \param[in]  a_preCopy
 *                      a_avgFab only needs the laplacian added
 *                      (if true)
 *//*-----------------------------------------------------------------*/

void
CRDutil::convolveFace(FArrayBox&           a_avgFab,
                      const FArrayBox&     a_pntFab,
                      const Box&           a_box,
                      const ProblemDomain& a_domain,
                      const Interval&      a_intv,
                      const int            a_dir,
                      const int            a_order,
                      bool                 a_interior,
                      bool                 a_limit,
                      bool                 a_preCopy)
{
  CH_TIME("CRDutil::convolveFace");
  CRD::msg << CRD::fv4 << "CRDutil::convolveFace" << CRD::end;
  //**FIXME: a lot of things need to be fixed
  //         1) add CH_asserts for box sizes and component counts
  //            a) some odd ones that should be checked include a_box
  //               being inside a_domain's box
  //         2) need to add in the limiting (if necessary)
  //         3) need to change this to use a face-laplacian operator
  // We require that we're given a face-box in direction a_dir
  IntVect testBoxType = a_box.type();
  IntVect cellBoxType(IntVect::Zero);
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      if ((dir == a_dir) &&
          (testBoxType[dir] == cellBoxType[dir]))
        {
          CH_assert(false);
        }
      else if ((dir != a_dir) &&
               (testBoxType[dir] != cellBoxType[dir])) // not allowed
        {
          CH_assert(false);
        }
    }

  // Check that a_box is inside a_domain
  if (!a_interior)
    {
      CH_assert(a_domain.contains(a_box));
    }

  constexpr Real factor24 = 1./24.;

  if (a_order == 2)
    {
      if (a_preCopy)
        {
          return;
        }
      a_avgFab.copy(a_pntFab, a_box, a_intv.begin(), a_box, a_intv.begin(),
                    a_intv.size());
      return;
    }

  int limitStencil = (int)a_limit;
  IntVect growVect(limitStencil*IntVect::Unit);
  growVect[a_dir] = 0;
  Box D2box = grow(a_box, growVect);
  D2box &= a_domain;

  FABSTACKTEMP(D2pntFab, D2box, 1);
  for (int iComp = a_intv.begin(), iComp_end = a_intv.end() + 1;
       iComp != iComp_end; ++iComp)
    {
      D2pntFab.setVal(0.);
      for (int tdir = 0; tdir != SpaceDim; ++tdir)
        {
          if (tdir != a_dir)
            {
              const int MD_ID(o, tdir);
              if (!a_interior)
                {
                  Box loBox, hiBox, centerBox, entireBox;
                  int hasLo, hasHi;
                  Box inputBox = D2box;
                  inputBox.grow(tdir, 1);
                  loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox,
                             entireBox, inputBox, a_domain, tdir);
                  if (hasLo)
                    {
                      loBox &= a_domain;
                    }
                  if (hasHi)
                    {
                      hiBox &= a_domain;
                    }
                  MD_BOXLOOP(centerBox, i)
                    {
                      D2pntFab[MD_IX(i, 0)] += (
                        a_pntFab[MD_OFFSETIX(i,-,o, iComp)]
                        - 2.*a_pntFab[MD_IX(i, iComp)]
                        +    a_pntFab[MD_OFFSETIX(i,+,o, iComp)]);
                    }
                  if (hasLo)
                    {
                      MD_BOXLOOP(loBox, i)
                        {
                          D2pntFab[MD_IX(i, 0)] += (
                            2.*a_pntFab[MD_IX(i, iComp)]
                            - 5.*a_pntFab[MD_OFFSETIX(i,+,o, iComp)]
                            + 4.*a_pntFab[MD_OFFSETIX(i,+,2*o, iComp)]
                            -    a_pntFab[MD_OFFSETIX(i,+,3*o, iComp)]);
                        }
                    }
                  if (hasHi)
                    {
                      MD_BOXLOOP(hiBox, i)
                        {
                          D2pntFab[MD_IX(i, 0)] += (
                            2.*a_pntFab[MD_IX(i, iComp)]
                            - 5.*a_pntFab[MD_OFFSETIX(i,-,o, iComp)]
                            + 4.*a_pntFab[MD_OFFSETIX(i,-,2*o, iComp)]
                            -    a_pntFab[MD_OFFSETIX(i,-,3*o, iComp)]);
                        }
                    }
                }
              else // use centered stencil everywhere
                {
                  MD_BOXLOOP(D2box, i)
                    {
                      D2pntFab[MD_IX(i, 0)] += (
                        a_pntFab[MD_OFFSETIX(i,-,o, iComp)]
                        - 2.*a_pntFab[MD_IX(i, iComp)]
                        +    a_pntFab[MD_OFFSETIX(i,+,o, iComp)]);
                    }
                }
            }
        }
      if (a_preCopy)  // No limiting
        {
          MD_BOXLOOP(a_box, i)
            {
              a_avgFab[MD_IX(i, iComp)] += factor24*D2pntFab[MD_IX(i, 0)];
            }
        }
      else // No limiting and no prefilling of a_pntFab
        {
          MD_BOXLOOP(a_box, i)
            {
              a_avgFab[MD_IX(i, iComp)] = a_pntFab[MD_IX(i, iComp)]
                + factor24*D2pntFab[MD_IX(i, 0)];
            }
        }
    }  // Loop over components
}

/*--------------------------------------------------------------------*/
//  Compute face-averaged values from face-centered values
/** \param[out] a_avgFxb
 *                      FluxBox of face-averaged values
 *  \param[out] a_pntFxb
 *                      FluxBox of face-centered values
 *  \param[in]  a_box   Box over which to compute face-averaged values
 *  \param[in]  a_domain
 *                      Problem domain for checking boundaries
 *  \param[in]  a_intv  Interval of components for computing
 *                      face-averaged values
 *  \param[in]  a_order Order of accuracy of convolution
 *  \param[in]  a_interior
 *                      If false, one-sided convolutions are used at
 *                      domain boundaries
 *  \param[in]  a_limit Limit the convolution (if true)
 *  \param[in]  a_preCopy
 *                      a_avgFxb only needs the laplacian added
 *                      (if true)
 *//*-----------------------------------------------------------------*/

void
CRDutil::convolveFace(FluxBox&             a_avgFxb,
                      const FluxBox&       a_pntFxb,
                      const Box&           a_box,
                      const ProblemDomain& a_domain,
                      const Interval&      a_intv,
                      const int            a_order,
                      bool                 a_interior,
                      bool                 a_limit,
                      bool                 a_preCopy)
{
  CH_TIME("CRDutil::convolveFace");
  CRD::msg << CRD::fv4 << "CRDutil::convolveFace" << CRD::end;
  //**FIXME: there's a lot to fix here (CH_asserts)
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      const Box inputBox = surroundingNodes(a_box, dir);
      convolveFace(a_avgFxb[dir],
                   a_pntFxb[dir],
                   inputBox,
                   a_domain,
                   a_intv,
                   dir,
                   a_order,
                   a_interior,
                   a_limit,
                   a_preCopy);
    }
}

/*--------------------------------------------------------------------*/
//  Compute the min and max values in a local neighborhood of a cell
/** \param[out] a_localMin
 *                      Fab of min values over local regions
 *  \param[out] a_localMax
 *                      Fab of max values over local regions
 *  \param[in]  a_data  Fab of input data
 *  \param[in]  a_box   Box over which to check min and max values
 *  \param[in]  a_domain
 *                      Problem domain for checking boundaries
 *  \param[in]  a_intv  Interval of components for checking min and
 *                      max values in a_data
 *  \param[in]  a_localCompBeg
 *                      Beginning component of a_localMin/a_localMax
 *//*-----------------------------------------------------------------*/

void
CRDutil::localBounds(FArrayBox&           a_localMin,
                     FArrayBox&           a_localMax,
                     const FArrayBox&     a_data,
                     const Box&           a_box,
                     const ProblemDomain& a_domain,
                     const Interval&      a_intv,
                     const int            a_localCompBeg)
{
  CH_TIME("CRDutil::localBounds");
  CRD::msg << CRD::fv4 << "CRDutil::localBounds" << CRD::end;
  CH_assert(a_domain.contains(a_box));
  CH_assert(a_localMin.contains(a_box));
  CH_assert(a_localMax.contains(a_box));
  CH_assert(a_data.interval().contains(a_intv));
  const Interval localIntv(a_localCompBeg,
                           a_localCompBeg + a_intv.size() - 1);
  CH_assert(a_localMin.interval().contains(localIntv));
  CH_assert(a_localMax.interval().contains(localIntv));

  const Box box1 = grow(a_box, 1);
  Box scanBox = box1 & a_domain;  // Cells accessed for local min/max
  CH_assert(a_data.contains(scanBox));

  a_localMin.copy(a_data, a_box, a_intv.begin(), a_box, a_localCompBeg,
                  a_intv.size());
  a_localMax.copy(a_data, a_box, a_intv.begin(), a_box, a_localCompBeg,
                  a_intv.size());
  MD_ARRAY_RESTRICT(arrMin, a_localMin);
  MD_ARRAY_RESTRICT(arrMax, a_localMax);
  MD_ARRAY_RESTRICT(arrData, a_data);

  if (!scanBox.contains(box1))  // Have to deal with boundaries
    {
      for (int cData = a_intv.begin(), cData_end = a_intv.end() + 1;
           cData != cData_end; ++cData)  // c for data
        {
          // c for localMin/localMax
          const int cLocal = cData - a_intv.begin() + a_localCompBeg;
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              for (SideIterator side; side.ok(); ++side)
                {
                  if (scanBox.sideEnd(side())[dir] ==
                      a_box.sideEnd(side())[dir])
                    {
                      Box boBox = adjCellBox(a_box, dir, side(), -1);
                      MD_BOXLOOP(boBox, i)
                        {
                          Real minVal = arrMin[MD_IX(i, cLocal)];
                          Real maxVal = arrMax[MD_IX(i, cLocal)];
                          const IntVect iv{MD_EXPANDIX(i)};
                          Box boStencil(iv - IntVect_unit, iv + IntVect_unit);
                          boStencil &= scanBox;
                          MD_BOXLOOP(boStencil, ii)
                            {
                              const Real val = arrData[MD_IX(ii, cData)];
                              minVal = std::min(minVal, val);
                              maxVal = std::max(maxVal, val);
                            }
                          arrMin[MD_IX(i, cLocal)] = minVal;
                          arrMax[MD_IX(i, cLocal)] = maxVal;
                        }
                    }
                }
            }
        }
    }

  constexpr Box c_stencil(-IntVect_unit, IntVect_unit);

  // Shrink the scan box by 1 to get all cells where we can apply a completely
  // centered stencil
  scanBox.grow(-1);
  if (!scanBox.isEmpty())
    {
      for (int cData = a_intv.begin(), cData_end = a_intv.end() + 1;
           cData != cData_end; ++cData)  // c for data
        {
          // c for localMin/localMax
          const int cLocal = cData - a_intv.begin() + a_localCompBeg;
          MD_BOXLOOP(scanBox, i)
            {
              Real minVal = arrMin[MD_IX(i, cLocal)];
              Real maxVal = arrMax[MD_IX(i, cLocal)];
              MD_BOXLOOP(c_stencil, ii)
                {
                  const Real val = arrData[MD_OFFSETIX(i,+,ii, cData)];
                  minVal = std::min(minVal, val);
                  maxVal = std::max(maxVal, val);
                }
              arrMin[MD_IX(i, cLocal)] = minVal;
              arrMax[MD_IX(i, cLocal)] = maxVal;
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute the min and max values in a local neighborhood of a face
/** \param[out] a_localMin
 *                      Fab of min values over local regions
 *  \param[out] a_localMax
 *                      Fab of max values over local regions
 *  \param[in]  a_data  Fab of input data
 *  \param[in]  a_box   Box over which to check min and max values
 *  \param[in]  a_domain
 *                      Problem domain for checking boundaries
 *  \param[in]  a_intv  Interval of components for checking min and
 *                      max values in a_data
 *  \param[in]  a_localCompBeg
 *                      Beginning component of a_localMin/a_localMax
 *  \param[in]  a_dir   Direction of faces of a_localMin/a_localMax
 *//*-----------------------------------------------------------------*/

void
CRDutil::localBoundsFace(FArrayBox&           a_localMin,
                         FArrayBox&           a_localMax,
                         const FArrayBox&     a_data,
                         const Box&           a_box,
                         const ProblemDomain& a_domain,
                         const Interval&      a_intv,
                         const int            a_localCompBeg,
                         const int            a_dir)
{
  CH_TIME("CRDutil::localBoundsFace");
  CRD::msg << CRD::fv4 << "CRDutil::localBoundsFace" << CRD::end;
  //**FIXME: add in ch_asserts
  //**FIXME: move most of the framework for this to neighborApply or
  //         something similar
  //**FIXME: make a neighborApplyCell and neighborApplyFace
  //**FIXME: if done properly, this function should require 2 evaluations
  //         per evaluation-dimension
  for (int cData = a_intv.begin(), cData_end = a_intv.end() + 1;
       cData != cData_end; ++cData)  // c for data
    {
      // c for localMin/localMax
      const int cLocal = cData - a_intv.begin() + a_localCompBeg;

      // set a_localMin and a_localMax to the a_data value
      MD_BOXLOOP(a_box, i)
        {
          a_localMin[MD_IX(i, cLocal)] = a_data[MD_IX(i,cData)];
          a_localMax[MD_IX(i, cLocal)] = a_data[MD_IX(i,cData)];
        }

      const int rad = 1; // radius of apply box
      IntVect loSide(-rad*IntVect::Unit);
      loSide[a_dir] = 0;
      IntVect hiSide(rad*IntVect::Unit);
      hiSide[a_dir] = 0;
      Box offsetBox(loSide, hiSide);
      MD_BOXLOOP(offsetBox, i) // a box-loop to represent neighbor IntVects
        {
          const IntVect os(MD_GETIV(i)); // get the neighbor IntVect
          Box applyBox(a_box);
          applyBox.shift(os); // shift the box to the neighbor IntVect
          applyBox &= a_domain; // see if the neighbor IntVect exists
          applyBox.shift(-os); // shift the box back to original location
          if (os != IntVect::Zero)
            {
              MD_BOXLOOP(applyBox, j) // loop over the entire box
                {
                  Real minVal = a_localMin[MD_IX(j, cLocal)]; // local minimum
                  Real maxVal = a_localMax[MD_IX(j, cLocal)]; // local maximum
                  const Real val = a_data[MD_OFFSETIV(j,+,os,cData)];
                  a_localMin[MD_IX(j, cLocal)] = std::min(val, minVal);
                  a_localMax[MD_IX(j, cLocal)] = std::max(val, maxVal);
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute the min and max values in a local neighborhood of a cell
/** \param[out] a_localMin
 *                      Min value over local region
 *  \param[out] a_localMax
 *                      Max value over local region
 *  \param[in]  a_data  Fab of input data over local box
 *  \param[in]  a_iv    Intvect of center of local box over which
 *                      to check min and max values
 *  \param[in]  a_domain
 *                      Problem domain for checking boundaries
 *  \param[in]  a_comp  Component of a_data to evaluate
 *//*-----------------------------------------------------------------*/

void
CRDutil::localBounds(Real&                a_localMin,
                     Real&                a_localMax,
                     const FArrayBox&     a_data,
                     const IntVect&       a_iv,
                     const ProblemDomain& a_domain,
                     const int            a_comp)
{
  CH_TIME("CRDutil::localBounds");
  CRD::msg << CRD::fv4 << "CRDutil::localBounds" << CRD::end;
  CH_assert(a_domain.contains(a_iv));
  CH_assert(a_data.interval().contains(a_comp));

  const Box stencilBox = (grow(Box(a_iv, a_iv), 1)) & a_domain;
  CH_assert(a_data.contains(stencilBox));

  a_localMin = a_data[MD_IV(a_iv, a_comp)];
  a_localMax = a_localMin;
  MD_BOXLOOP(stencilBox, i)
    {
      const Real val = a_data[MD_IX(i, a_comp)];
      a_localMin = std::min(a_localMin, val);
      a_localMax = std::max(a_localMax, val);
    }
}

/*--------------------------------------------------------------------*/
//  Check positivity of cell values
/** \param[in]  a_inputFab
 *                      Fab of values to check positivity of
 *  \param[in]  a_dependentFab
 *                      Fab of values used to obtain a_inputFab
 *  \param[in]  a_box   Box on which to check a_inputFab
 *  \param[in]  a_domain
 *                      Problem domain for boundary checking
 *  \param[in]  a_localComp
 *                      Component of a_inputFab to check
 *  \param[in]  a_neighborRadius
 *                      Radius of faces of a_dependentFab on which
 *                      a_inputFab depends (print these to terminal)
 *  \param[in]  a_inputVarName
 *                      Name of input parameter to print to terminal
 *  \param[in]  a_dependentVarName
 *                      Name of dependent parameter to print
 *//*-----------------------------------------------------------------*/

void
CRDutil::checkCellPositivity(const FArrayBox&     a_inputFab,
                             const FArrayBox&     a_dependentFab,
                             const Box&           a_box,
                             const ProblemDomain& a_domain,
                             const int            a_localComp,
                             const int            a_neighborRadius,
                             const std::string&   a_inputVarName,
                             const std::string&   a_dependentVarName)
{
  CH_TIME("CRDutil::checkCellPositivity");
  CRD::msg << CRD::fv4 << "CRDutil::checkCellPositivity" << CRD::end;
  CH_assert(a_domain.contains(a_box));
  CH_assert(a_inputFab.contains(a_box));
  CH_assert(a_dependentFab.contains(a_box));
  CH_assert(a_inputFab.nComp() >= a_localComp);
  CH_assert(a_localComp >= 0);
  CH_assert(a_dependentFab.nComp() >= a_localComp);
  MD_BOXLOOP(a_box, i)
    {
      Real val    = a_inputFab[MD_IX(i, a_localComp)];
      Real depVal = a_dependentFab[MD_IX(i, a_localComp)];
      if ((val <= 0.) || std::isnan(val) || std::isinf(val))
        {
          CRD::msg << "WARNING: in cell " << MD_GETIV(i) << ": "
                   << a_inputVarName << " = " << val << ": "
                   << a_dependentVarName << " = " << depVal << CRD::warn;
          if (a_neighborRadius > 0)
            {
              Box neighborBox = grow(Box(MD_GETIV(i), MD_GETIV(i)),
                                     a_neighborRadius);
              neighborBox &= a_domain;
              MD_BOXLOOP(neighborBox, j)
                {
                  Real neighborVal = a_inputFab[MD_IX(j, a_localComp)];
                  Real neighborDependentVal =
                    a_dependentFab[MD_IX(j, a_localComp)];
                  CRD::msg << "  Neighbor cell " << MD_GETIV(j)
                           << ": " << a_inputVarName << " = " << neighborVal
                           <<  ": " << a_dependentVarName << " = "
                           << neighborDependentVal << CRD::warn;
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Check positivity of face values
/** \param[in]  a_inputFab
 *                      Fab of values to check positivity of
 *  \param[in]  a_dependentFab
 *                      Fab of values used to obtain a_inputFab
 *  \param[in]  a_box   Box on which to check a_inputFab
 *  \param[in]  a_domain
 *                      Problem domain for boundary checking
 *  \param[in]  a_localComp
 *                      Component of a_inputFab to check
 *  \param[in]  a_dir   Direction of faces of a_inputFab
 *  \param[in]  a_neighborRadius
 *                      Radius of faces of a_dependentFab on which
 *                      a_inputFab depends (print these to terminal)
 *  \param[in]  a_inputVarName
 *                      Name of input parameter to print to terminal
 *  \param[in]  a_dependentVarName
 *                      Name of dependent parameter to print
 *//*-----------------------------------------------------------------*/

void
CRDutil::checkFacePositivity(const FArrayBox&     a_inputFab,
                             const FArrayBox&     a_dependentFab,
                             const Box&           a_box,
                             const ProblemDomain& a_domain,
                             const int            a_localComp,
                             const int            a_dir,
                             const int            a_neighborRadius,
                             const std::string&   a_inputVarName,
                             const std::string&   a_dependentVarName)
{
  CH_TIME("CRDutil::checkFacePositivity");
  CRD::msg << CRD::fv4 << "CRDutil::checkFacePositivity" << CRD::end;
  CH_assert(a_domain.contains(a_box));
  CH_assert(a_inputFab.contains(a_box));
  CH_assert(a_dependentFab.contains(a_box));
  CH_assert(a_inputFab.nComp() >= a_localComp);
  CH_assert(a_localComp >= 0);
  CH_assert(a_dependentFab.nComp() >= a_localComp);
  MD_BOXLOOP(a_box, i)
    {
      Real val    = a_inputFab[MD_IX(i, a_localComp)];
      Real depVal = a_dependentFab[MD_IX(i, a_localComp)];
      if ((val <= 0.) || std::isnan(val) || std::isinf(val))
        {
          CRD::msg << "WARNING: on face " << MD_GETIV(i) << ": "
                   << a_inputVarName << " = " << val << ": "
                   << a_dependentVarName << " = " << depVal << CRD::warn;
          if (a_neighborRadius > 0)
            {
              Box neighborBox = grow(Box(MD_GETIV(i), MD_GETIV(i)),
                                     a_neighborRadius);
              neighborBox.grow(a_dir, -a_neighborRadius);
              neighborBox &= a_domain;
              MD_BOXLOOP(neighborBox, j)
                {
                  Real neighborVal = a_inputFab[MD_IX(j, a_localComp)];
                  Real neighborDependentVal =
                    a_dependentFab[MD_IX(j, a_localComp)];
                  CRD::msg << "  Neighbor face " << MD_GETIV(j)
                           << ": " << a_inputVarName << " = " << neighborVal
                           <<  ": " << a_dependentVarName << " = "
                           << neighborDependentVal << CRD::warn;
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Check positivity of face values
/** \param[in]  a_inputFxb
 *                      FluxBox of values to check positivity of
 *  \param[in]  a_dependentFxb
 *                      FluxBox of values used to obtain a_inputFxb
 *  \param[in]  a_box   Box on which to check a_inputFxb
 *  \param[in]  a_domain
 *                      Problem domain for boundary checking
 *  \param[in]  a_localComp
 *                      Component of a_inputFxb to check
 *  \param[in]  a_neighborRadius
 *                      Radius of faces of a_dependentFxb on which
 *                      a_inputFxb depends (print these to terminal)
 *  \param[in]  a_inputVarName
 *                      Name of input parameter to print to terminal
 *  \param[in]  a_dependentVarName
 *                      Name of dependent parameter to print
 *//*-----------------------------------------------------------------*/

void
CRDutil::checkFacePositivity(const FluxBox&       a_inputFab,
                             const FluxBox&       a_dependentFab,
                             const Box&           a_box,
                             const ProblemDomain& a_domain,
                             const int            a_localComp,
                             const int            a_neighborRadius,
                             const std::string&   a_inputVarName,
                             const std::string&   a_dependentVarName)
{
  CH_TIME("CRDutil::checkFacePositivity");
  CRD::msg << CRD::fv4 << "CRDutil::checkFacePositivity" << CRD::end;
  //**FIXME: there's a lot to fix here (CH_asserts)
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      checkFacePositivity(a_inputFab[dir],
                          a_dependentFab[dir],
                          a_box,
                          a_domain,
                          a_localComp,
                          dir,
                          a_neighborRadius,
                          a_inputVarName,
                          a_dependentVarName);
    }
}

/*--------------------------------------------------------------------*/
//  Modify velocity gradient to incorporate wall-shear-stress from model
/** \param[out] a_NGradUfacePntFxb
 *                      Physical-space velocity-gradient
 *  \param[in]  a_unitNormalsFxb
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive state
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_domain
 *                      Problem domain for the level
 *  \param[in]  a_box   Box over which wall velocity-gradients
 *                      are modified based on wall-model
 *  \param[in]  a_disjointBox
 *                      Disjoint box of this dataIndex
 *//*-----------------------------------------------------------------*/

void
CRDutil::enforceLawOfTheWallShearStress(
  FluxBox&                a_NGradUfacePntFxb,
  const FluxBox&          a_unitNormalsFxb,
  const FluxBox&          a_WfaceAvgFxb,
  const LevelGridMetrics& a_gridMetrics,
  const ProblemDomain&    a_domain,
  const Box&              a_box,
  const Box&              a_disjointBox)
{
  CH_TIME("CRDutil::enforceLawOfTheWallShearStress");

  // Loop over all face directions
  for (int dirSide = 0; dirSide != 2*SpaceDim; ++dirSide)
    {
      // The following few lines save us an indent by merging a "side loop"
      // with a "dir loop"
      // Translate the dirSide index into dir
      const int dir = dirSide/2;
      // Translate the dirSide index into side
      const int side = dirSide % 2;
      // Current side
      Side::LoHiSide whichSide = (side == 1) ? Side::Hi : Side::Lo;
      // Opposite side of current side (use for opposite shift direction)
      Side::LoHiSide oppSide = flip(whichSide);
      // Sign for shifting out of the domain: -1 if lo, +1 if hi
      const int normSign = sign(whichSide);
      // Sign for shifting into the domain: +1 if lo, -1 if hi
      const int offsetSign = sign(oppSide);
      // Constant indices
      const int cVel = CRDparam::g_CRDPhysics->velocityInterval().begin();
      const int cRho = CRDparam::g_CRDPhysics->densityIndex();
      // Check if there is a boundary box
      Box bndryFaceBox;
      CRDparam::g_CNSIBC->getBoundaryFaces(
        bndryFaceBox, a_box, a_domain, dir, whichSide);
      // Check if this is a no-slip wall
      if (!bndryFaceBox.isEmpty())
        {
          BoundaryIndex bcID;
          bcID.define(a_gridMetrics.getCoordSys().whichBlock(a_disjointBox),
                      dir, whichSide);
          BCInfo bc = CRDparam::g_CNSIBC->getDomainBC(bcID);
          if ((CRDparam::DomainBCTypeSWall & bc.m_type) &&
              !(CRDparam::DomainBCTypeMixed & bc.m_type))
            {
              const FArrayBox& unitNormalsFab = a_unitNormalsFxb[dir];
              const FArrayBox& WfaceAvgFab = a_WfaceAvgFxb[dir];
              // Get the block coordinate system
              const BlockCoordSys& blockCoordSys =
                *(a_gridMetrics.getCoordSys(a_disjointBox));
              // First, compute the wall shear-stress from Spalding's law
              FABSTACKTEMP(eta, bndryFaceBox, 1);
              eta.setVal(0.);
              // Get the physical-space coordinates of the wall face
              FABSTACKTEMP(faceXFab, bndryFaceBox, SpaceDim);
              FABSTACKTEMP(faceXiFab, bndryFaceBox, SpaceDim);
              CRDparam::g_CNSIBC->getFaceCoordinates(
                bndryFaceBox, faceXiFab, faceXFab, dir, blockCoordSys);
              // Get the physical-space coordinates of the 1st interior face
              Box firstIntFaceBox = bndryFaceBox;
              firstIntFaceBox.shift(dir, offsetSign);
              FABSTACKTEMP(face1XFab, firstIntFaceBox, SpaceDim);
              FABSTACKTEMP(face1XiFab, firstIntFaceBox, SpaceDim);
              CRDparam::g_CNSIBC->getFaceCoordinates(
                firstIntFaceBox, face1XiFab, face1XFab, dir, blockCoordSys);
              // Shift the second interior face coordinates to the wall
              face1XFab.shift(dir, normSign);
              // Get the physical-space coordinates of the 2nd interior face
              Box secondIntFaceBox = bndryFaceBox;
              secondIntFaceBox.shift(dir, 2*offsetSign);
              FABSTACKTEMP(face2XFab, secondIntFaceBox, SpaceDim);
              FABSTACKTEMP(face2XiFab, secondIntFaceBox, SpaceDim);
              CRDparam::g_CNSIBC->getFaceCoordinates(
                secondIntFaceBox, face2XiFab, face2XFab, dir, blockCoordSys);
              // Shift the second interior face coordinates to the wall
              face2XFab.shift(dir, 2*normSign);
              // Compute the distance from the wall to the 2nd interior face
              //**NOTE: This method isn't fully correct, but it's close.
              //        If necessary, use the surrounding faces to construct
              //        an estimate for the curvature of the face and the true
              //        distance to the wall.
              MD_BOXLOOP(bndryFaceBox, i)
                {
                  for (int comp = 0; comp != SpaceDim; ++comp)
                    {
                      face1XFab[MD_IX(i, comp)] =
                        offsetSign*(face1XFab[MD_IX(i, comp)]
                                  - faceXFab[MD_IX(i, comp)]);
                      faceXFab[MD_IX(i, comp)] =
                        offsetSign*(face2XFab[MD_IX(i, comp)]
                                  - faceXFab[MD_IX(i, comp)]);
                    }
                }
              // Transform both wall-distance vectors into normal-tangent space
              FORT_FORWARDTRANSFORMF(CHF_FRA(faceXFab),
                                     CHF_CONST_FRA(unitNormalsFab),
                                     CHF_BOX(bndryFaceBox));
              FORT_FORWARDTRANSFORMF(CHF_FRA(face1XFab),
                                     CHF_CONST_FRA(unitNormalsFab),
                                     CHF_BOX(bndryFaceBox));
              // Shift the 2nd interior-face velocity to the wall
              FABSTACKTEMP(interiorVel, secondIntFaceBox, SpaceDim);
              interiorVel.copy(WfaceAvgFab, secondIntFaceBox, cVel,
                               secondIntFaceBox, 0, SpaceDim);
              interiorVel.shift(dir, 2*normSign);
              // Transform the face-averaged velocity into normal-tangent space
              //**NOTE: This method isn't fully 4th-order. The velocity needs
              //        to be face-centered data in order to be 
              FORT_FORWARDTRANSFORMF(CHF_FRA(interiorVel),
                                     CHF_CONST_FRA(unitNormalsFab),
                                     CHF_BOX(bndryFaceBox));
              // Set the wall-normal velocity to zero
              interiorVel.setVal(0., dir);
              // Calculate the viscosity (and kappa unfortunately) at the wall
              FABSTACKTEMP(muFace, bndryFaceBox, 1);
              FABSTACKTEMP(kappaFace, bndryFaceBox, 1);
              CRDparam::g_CRDPhysics->calcCoeffKappaMu(
                bndryFaceBox, muFace, kappaFace, WfaceAvgFab);
              // Turn mu into nu
              MD_BOXLOOP(bndryFaceBox, i)
                {
                  muFace[MD_IX(i, 0)] /= WfaceAvgFab[MD_IX(i, cRho)];
                }
              // Compute eta_0 (wall-normal deriv of streamwise-velocity)
              FABSTACKTEMP(etaFab, bndryFaceBox, 1);
              FABSTACKTEMP(uTauFab, bndryFaceBox, 1);
              MD_BOXLOOP(bndryFaceBox, i)
                {
                  const Real nu = muFace[MD_IX(i, 0)];
                  const Real yVal = faceXFab[MD_IX(i, dir)];
                  RealVect streamwiseVel(D_DECL(interiorVel[MD_IX(i, 0)],
                                                interiorVel[MD_IX(i, 1)],
                                                interiorVel[MD_IX(i, 2)]));
                  Real u_tau = 1.e-10;
                  Real uVel = streamwiseVel.vectorLength();
                  if (uVel >= 0.) // If uVel = 0, u_tau is indeterminant
                    {
                      // Friction velocity (based on Spalding)
                      Real yPlusMaxGuess = 1000000.;
                      Real uTauMin = std::sqrt(uVel*nu/yVal);
                      Real uTauMax = yPlusMaxGuess*nu/yVal;
                      int iterBrent  = 0;
                      int errorBrent = 0;

                      const UTauSpaldingFunc& f =
                        UTauSpaldingFunc(yVal,nu,uVel);
                      u_tau = RootSolver::BrentER(
                        iterBrent,errorBrent,f,uTauMin,uTauMax);
                      if (errorBrent != 0 || u_tau != u_tau)
                        {
                          CRD::msg << "SSV: Bad uTau value: " << u_tau
                                   << CRD::error;
                        }
                    }
                  uTauFab[MD_IX(i, 0)] = u_tau;
                  etaFab[MD_IX(i, 0)] = u_tau*u_tau/nu;
                }
              // Normalize wall-tangent vector -- streamwise unit vector
              PatchMappedFunc::normalize(
                bndryFaceBox, interiorVel, Interval(0,SpaceDim-1));

              // Face-centered physical-space velocity-gradients
              FArrayBox& facePntVelGrad = a_NGradUfacePntFxb[dir];
              // Velocity gradients modified by LES wall-model
              FABSTACKTEMP(modelVelGrad, bndryFaceBox, SpaceDim*SpaceDim);
              // Everything except for normal derivs should be zero
              modelVelGrad.setVal(0.);
              // Get the face-normal vector at the boundary
              FABSTACKTEMP(normVect, bndryFaceBox, SpaceDim);
              normVect.setVal(0.);
              normVect.setVal(1., dir);
              // Transform normalVect into Cartesian physical space
              FORT_REVERSETRANSFORMF(CHF_FRA(normVect),
                                     CHF_CONST_FRA(unitNormalsFab),
                                     CHF_BOX(bndryFaceBox));
              // Compute wall-normal derivative of wall-normal velocity
              // Just extract it from a_NGradUfacePntFxb
              MD_BOXLOOP(bndryFaceBox, i)
                {
                  Real normDeriv = 0.;
                  for (int comp = 0; comp != SpaceDim; ++comp)
                    {
                      normDeriv += normVect[MD_IX(i, comp)]*(
                        D_TERM(
                          facePntVelGrad[MD_IX(i,comp*SpaceDim)]*
                          normVect[MD_IX(i,0)],
                          + facePntVelGrad[MD_IX(i,comp*SpaceDim+1)]*
                          normVect[MD_IX(i,1)],
                          + facePntVelGrad[MD_IX(i,comp*SpaceDim+2)]*
                          normVect[MD_IX(i,2)]));
                    }
                  modelVelGrad[MD_IX(i,dir*SpaceDim+dir)] = normDeriv;
                }

              // Fill normal-derivs of tangential velocities using model
              for (int velComp = 0; velComp != SpaceDim; ++velComp)
                {
                  if (velComp != dir)
                    {
                      const int cGrad = dir + SpaceDim*velComp;
                      MD_BOXLOOP(bndryFaceBox, i)
                        {
                          // Make sure to account for normal direction
                          modelVelGrad[MD_IX(i, cGrad)] =
                            offsetSign*interiorVel[MD_IX(i,velComp)]*
                            etaFab[MD_IX(i,0)];
                        }
                    }
                }
              // Transform back to Cartesian space (a lot here)
              FABSTACKTEMP(tempWallVelGrad, bndryFaceBox, SpaceDim*SpaceDim);
              // First, right multiply grad(u) by a_unitNormals
              for (int row = 0; row != SpaceDim; ++row)
                {
                  for (int col = 0; col != SpaceDim; ++col)
                    {
                      const int cGrad = col + SpaceDim*row;
                      MD_BOXLOOP(bndryFaceBox, i)
                        {
                          tempWallVelGrad[MD_IX(i, cGrad)] =
                            D_TERM(
                              modelVelGrad[MD_IX(i,row*SpaceDim)]*
                              unitNormalsFab[MD_IX(i,col)],
                              + modelVelGrad[MD_IX(i,row*SpaceDim+1)]*
                              unitNormalsFab[MD_IX(i,SpaceDim+col)],
                              + modelVelGrad[MD_IX(i,row*SpaceDim+2)]*
                              unitNormalsFab[MD_IX(i,2*SpaceDim+col)]);
                        }
                    }
                }
              // Next, left multiply by a_unitNormals^T
              FABSTACKTEMP(wallVelGrad, bndryFaceBox, SpaceDim*SpaceDim);
              for (int row = 0; row != SpaceDim; ++row)
                {
                  for (int col = 0; col != SpaceDim; ++col)
                    {
                      const int cGrad = col + SpaceDim*row;
                      MD_BOXLOOP(bndryFaceBox, i)
                        {
                          wallVelGrad[MD_IX(i, cGrad)] =
                            D_TERM(
                              unitNormalsFab[MD_IX(i,row)]*
                              tempWallVelGrad[MD_IX(i,col)],
                              + unitNormalsFab[MD_IX(i,SpaceDim+row)]*
                              tempWallVelGrad[MD_IX(i,SpaceDim+col)],
                              + unitNormalsFab[MD_IX(i,2*SpaceDim+row)]*
                              tempWallVelGrad[MD_IX(i,2*SpaceDim+col)]);
                        }
                    }
                }
              // Finally, copy the result into a_NGradUfacePntFxb
              for (int c = 0; c != SpaceDim*SpaceDim; ++c)
                {
                  MD_BOXLOOP(bndryFaceBox, i)
                    {
                      facePntVelGrad[MD_IX(i, c)] = wallVelGrad[MD_IX(i, c)];
                    }
                }

              // DEBUG ONLY
              // Try modifying the wall-normal derivatives at the first interior
              // face. Right now, the wall-normal derivs at the first interior
              // face are negative for coarse flat-plate simulations due to a
              // systemic oscillation (this should be investigated further)
              if (true)
                {
                  Box intFaceBox = bndryFaceBox;
                  intFaceBox.shift(dir, offsetSign);
                  // Velocity gradients modified by LES wall-model
                  FABSTACKTEMP(intVelGrad, intFaceBox, SpaceDim*SpaceDim);
                  intVelGrad.copy(facePntVelGrad);
                  // Transform intVelGrad into wall-normal space
                  intVelGrad.shift(dir,normSign); // Shift to wall temporarily
                  // Transform back to Cartesian space (a lot here)
                  FABSTACKTEMP(tempGrad1, bndryFaceBox, SpaceDim*SpaceDim);
                  // First, right multiply grad(u) by a_unitNormals^T
                  for (int row = 0; row != SpaceDim; ++row)
                    {
                      for (int col = 0; col != SpaceDim; ++col)
                        {
                          const int cGrad = col + SpaceDim*row;
                          MD_BOXLOOP(bndryFaceBox, i)
                            {
                              tempGrad1[MD_IX(i, cGrad)] =
                                D_TERM(
                                  intVelGrad[MD_IX(i,row*SpaceDim)]*
                                  unitNormalsFab[MD_IX(i,col*SpaceDim)],
                                  + intVelGrad[MD_IX(i,row*SpaceDim+1)]*
                                  unitNormalsFab[MD_IX(i,col*SpaceDim+1)],
                                  + intVelGrad[MD_IX(i,row*SpaceDim+2)]*
                                  unitNormalsFab[MD_IX(i,col*SpaceDim+2)]);
                            }
                        }
                    }
                  // Next, left multiply by a_unitNormals
                  FABSTACKTEMP(normGrad, bndryFaceBox, SpaceDim*SpaceDim);
                  for (int row = 0; row != SpaceDim; ++row)
                    {
                      for (int col = 0; col != SpaceDim; ++col)
                        {
                          const int cGrad = col + SpaceDim*row;
                          MD_BOXLOOP(bndryFaceBox, i)
                            {
                              normGrad[MD_IX(i, cGrad)] =
                                D_TERM(
                                  unitNormalsFab[MD_IX(i,row*SpaceDim)]*
                                  tempGrad1[MD_IX(i,col)],
                                  + unitNormalsFab[MD_IX(i,row*SpaceDim+1)]*
                                  tempGrad1[MD_IX(i,SpaceDim+col)],
                                  + unitNormalsFab[MD_IX(i,row*SpaceDim+2)]*
                                  tempGrad1[MD_IX(i,2*SpaceDim+col)]);
                            }
                        }
                    }
                  // Compute wall-normal derivatives at first interior face
                  FABSTACKTEMP(etaIntFab, bndryFaceBox, 1);
                  MD_BOXLOOP(bndryFaceBox, i)
                    {
                      const Real etaWall = etaFab[MD_IX(i,0)];
                      const Real y_val = face1XFab[MD_IX(i,dir)];
                      const Real nu = muFace[MD_IX(i, 0)];
                      const Real u_tau = uTauFab[MD_IX(i, 0)];
                      const Real y_p = y_val*u_tau/nu;
                      const Real k_1 = 0.41;
                      const Real c_1 = 0.001093;
                      const Real a_1 = k_1 + c_1*y_p*y_p;
                      const Real a_2 = a_1 + c_1*k_1*y_p*y_p*y_p;
                      etaIntFab[MD_IX(i, 0)] = etaWall*(a_1/a_2);
                    }
                  // Replace the two wall-normal components
                  for (int velComp = 0; velComp != SpaceDim; ++velComp)
                    {
                      if (velComp != dir)
                        {
                          const int cGrad = dir + SpaceDim*velComp;
                          MD_BOXLOOP(bndryFaceBox, i)
                            {
                              normGrad[MD_IX(i, cGrad)] =
                                offsetSign*interiorVel[MD_IX(i,velComp)]*
                                etaIntFab[MD_IX(i,0)];
                            }
                        }
                    }
                  // Transform back to Cartesian space
                  FABSTACKTEMP(tempGrad2, bndryFaceBox, SpaceDim*SpaceDim);
                  // First, right multiply grad(u) by a_unitNormals
                  for (int row = 0; row != SpaceDim; ++row)
                    {
                      for (int col = 0; col != SpaceDim; ++col)
                        {
                          const int cGrad = col + SpaceDim*row;
                          MD_BOXLOOP(bndryFaceBox, i)
                            {
                              tempGrad2[MD_IX(i, cGrad)] =
                                D_TERM(
                                  normGrad[MD_IX(i,row*SpaceDim)]*
                                  unitNormalsFab[MD_IX(i,col)],
                                  + normGrad[MD_IX(i,row*SpaceDim+1)]*
                                  unitNormalsFab[MD_IX(i,SpaceDim+col)],
                                  + normGrad[MD_IX(i,row*SpaceDim+2)]*
                                  unitNormalsFab[MD_IX(i,2*SpaceDim+col)]);
                            }
                        }
                    }
                  // Next, left multiply by a_unitNormals^T
                  for (int row = 0; row != SpaceDim; ++row)
                    {
                      for (int col = 0; col != SpaceDim; ++col)
                        {
                          const int cGrad = col + SpaceDim*row;
                          MD_BOXLOOP(bndryFaceBox, i)
                            {
                              intVelGrad[MD_IX(i, cGrad)] =
                                D_TERM(
                                  unitNormalsFab[MD_IX(i,row)]*
                                  tempGrad2[MD_IX(i,col)],
                                  + unitNormalsFab[MD_IX(i,SpaceDim+row)]*
                                  tempGrad2[MD_IX(i,SpaceDim+col)],
                                  + unitNormalsFab[MD_IX(i,2*SpaceDim+row)]*
                                  tempGrad2[MD_IX(i,2*SpaceDim+col)]);
                            }
                        }
                    }
                  // Shift back to interior face
                  intVelGrad.shift(dir, offsetSign);
                  // Finally, copy the result into a_NGradUfacePntFxb
                  for (int c = 0; c != SpaceDim*SpaceDim; ++c)
                    {
                      MD_BOXLOOP(intFaceBox, i)
                        {
                          facePntVelGrad[MD_IX(i, c)] = intVelGrad[MD_IX(i, c)];
                        }
                    }
                }
              // END DEBUG ONLY
            }
        }
    }
}
