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
 * \file CRDutil.H
 *
 * \brief Utilities for Chord
 *
 *//*+*************************************************************************/

#ifndef _CRDUTIL_H_
#define _CRDUTIL_H_

//----- Chombo Library -----//

#include "LevelGridMetrics.H"
#include "UsingNamespace.H"

//----- Internal -----//

//--Forward declaration from Chombo

#include "NamespaceHeader.H"
class Box;
class ProblemDomain;
class FArrayBox;
#include "NamespaceFooter.H"

/*******************************************************************************
 */
///  Utilities for Chord
/**
 *   This namespace is a random catch-all container for utilities that are
 *   used throughout Chord. Most of these utilities are specifically related
 *   to fourth-order components of the algorithm and are used for stability.
 *
 *//*+*************************************************************************/

namespace CRDutil
{

  /// PPM face values
  void PPMFaceValues(FArrayBox&           a_WFace,
                     const FArrayBox&     a_W,
                     const FArrayBox&     a_smoothTest,
                     const int&           a_numSlopes,
                     const bool&          a_useLimiting,
                     const int&           a_dir,
                     const Box&           a_faceBox,
                     const ProblemDomain& a_domain);

  /// Test the smoothness of an interpolant
  void smoothTest(FArrayBox&           a_smoothTest,
                  const FArrayBox&     a_W,
                  const int&           a_numSlopes,
                  const int&           a_dir,
                  const Box&           a_cellBox,
                  const ProblemDomain& a_domain);

  void PPMLimiter(FArrayBox&           a_dWMinus,
                  FArrayBox&           a_dWPlus,
                  const FArrayBox&     a_W,
                  const int            a_numSlopes,
                  const int            a_dir,
                  const Box&           a_box,
                  const ProblemDomain& a_domain);

  /// Compute the second-order second-derivative in a cell
  void getSecondDiff(FArrayBox&       a_d2W,
                     const FArrayBox& a_W,
                     const int        a_numSlopes,
                     const int        a_dir,
                     const Box&       a_loBox,
                     const int        a_hasLo,
                     const Box&       a_hiBox,
                     const int        a_hasHi,
                     const Box&       a_centerBox);

  void checkCubicLimiter(FArrayBox&       a_dWMinus,
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
                         const int        a_hiBoundHOCheck);

  /// Limiting near boundaries
  void extraBoundaryLimiting(FArrayBox&           a_WMinus,
                             FArrayBox&           a_WPlus,
                             const FArrayBox&     a_W,
                             const Box&           a_cellBox,
                             const int&           a_numSlopes,
                             const int&           a_dir,
                             const int&           a_level,
                             const Real&          a_time,
                             LevelGridMetrics&    a_gridMetrics,
                             const ProblemDomain& a_domain);

  /// PPM upwind face values
  void PPMUpwindFaceValues(FArrayBox&           a_WfaceMinus,
                           FArrayBox&           a_WfacePlus,
                           const FArrayBox&     a_WcellAvg,
                           const int&           a_numComp,
                           const int&           a_dir,
                           const Box&           a_box,
                           const ProblemDomain& a_domain);

  /// Compute the cell-averaged kinetic energy with a triple product rule
  void avgKE(FArrayBox&           a_avgKEfab,
             const FArrayBox&     a_UavgFab,
             const Box&           a_box,
             const ProblemDomain& a_domain,
             bool                 a_limit = true);

  /// Compute the cell-averaged internal energy from the kinetic energy
  void avgE(FArrayBox&           a_avgEfab,
            const FArrayBox&     a_barEfab,
            const FArrayBox&     a_UavgFab,
            const Box&           a_box,
            const ProblemDomain& a_domain,
            bool                 a_limit = true);

  /// Limit the deconvolution of temperature
  void limitTpnt(const FArrayBox&     a_UavgFab,
                 const FArrayBox&     a_WpntFab,
                 const Box&           a_box,
                 const ProblemDomain& a_domain,
                 bool                 a_limit = true);

  /// Compute the second-order Laplacian in a cell
  void Laplacian(FArrayBox&           a_D2fab,
                 const FArrayBox&     a_D0fab,
                 const Box&           a_box,
                 const ProblemDomain& a_domain,
                 const Interval&      a_intv,
                 bool                 a_cSto = 0);

  /// Compute the second-order second-derivative in a cell
  void getSecondDiff(FArrayBox&           a_d2Fab,
                     const FArrayBox&     a_vals,
                     const Box&           a_box,
                     const int            a_numComp,
                     const ProblemDomain& a_domain);

  /// Compute cell-averaged values from cell-centered values (NOT IMPLEMENTED)
  void convolve(FArrayBox& a_avgFab,
                const FArrayBox& a_pntFab,
                const Box& a_box,
                const int a_order,
                bool a_limit,
                int a_sign = 1);

  /// Compute cell-centered values from cell-averaged values
  void deconvolve(FArrayBox&           a_pntFab,
                  const FArrayBox&     a_avgFab,
                  const Box&           a_box,
                  const ProblemDomain& a_domain,
                  const Interval&      a_intv,
                  const int            a_order = 4,
                  const int            a_sign  = 1,
                  bool                 a_limit = true,
                  bool                 a_preCopy = false);

  /// Compute face-centered values from face-averaged values
  void deconvolveFace(FArrayBox&           a_pntFab,
                      const FArrayBox&     a_avgFab,
                      const Box&           a_box,
                      const ProblemDomain& a_domain,
                      const Interval&      a_intv,
                      const int            a_dir,
                      const int            a_order = 4,
                      bool                 a_limit = true,
                      bool                 a_preCopy = false);

  /// Compute face-centered values from face-averaged values
  void deconvolveFace(FluxBox&             a_pntFxb,
                      const FluxBox&       a_avgFxb,
                      const Box&           a_box,
                      const ProblemDomain& a_domain,
                      const Interval&      a_intv,
                      const int            a_order = 4,
                      bool                 a_limit = true,
                      bool                 a_preCopy = false);

  /// Compute face-averaged values from face-centered values
  void convolveFace(FArrayBox&           a_avgFab,
                    const FArrayBox&     a_pntFab,
                    const Box&           a_box,
                    const ProblemDomain& a_domain,
                    const Interval&      a_intv,
                    const int            a_dir,
                    const int            a_order = 4,
                    bool                 a_interior = false,
                    bool                 a_limit = true,
                    bool                 a_preCopy = false);

  /// Compute face-averaged values from face-centered values
  void convolveFace(FluxBox&             a_avgFxb,
                    const FluxBox&       a_pntFxb,
                    const Box&           a_box,
                    const ProblemDomain& a_domain,
                    const Interval&      a_intv,
                    const int            a_order = 4,
                    bool                 a_interior = false,
                    bool                 a_limit = true,
                    bool                 a_preCopy = false);

  /// Compute the min and max values in a local neighborhood of a cell
  void localBounds(FArrayBox&           a_localMin,
                   FArrayBox&           a_localMax,
                   const FArrayBox&     a_data,
                   const Box&           a_box,
                   const ProblemDomain& a_domain,
                   const Interval&      a_intv,
                   const int            a_localCompBeg = 0);

  /// Compute the min and max values in a local neighborhood of a face
  void localBoundsFace(FArrayBox&           a_localMin,
                       FArrayBox&           a_localMax,
                       const FArrayBox&     a_data,
                       const Box&           a_box,
                       const ProblemDomain& a_domain,
                       const Interval&      a_intv,
                       const int            a_localCompBeg,
                       const int            a_dir);

  /// Compute the min and max values in a local neighborhood of a cell
  void localBounds(Real&                a_localMin,
                   Real&                a_localMax,
                   const FArrayBox&     a_data,
                   const IntVect&       a_iv,
                   const ProblemDomain& a_domain,
                   const int            a_comp);

  /// Check positivity of cell values
  void checkCellPositivity(const FArrayBox&     a_inputFab,
                           const FArrayBox&     a_dependentFab,
                           const Box&           a_box,
                           const ProblemDomain& a_domain,
                           const int            a_localComp,
                           const int            a_neighborRadius,
                           const std::string&   a_inputVarName,
                           const std::string&   a_dependentVarName);

  /// Check positivity of face values
  void checkFacePositivity(const FArrayBox&     a_inputFab,
                           const FArrayBox&     a_dependentFab,
                           const Box&           a_box,
                           const ProblemDomain& a_domain,
                           const int            a_localComp,
                           const int            a_dir,
                           const int            a_neighborRadius,
                           const std::string&   a_inputVarName,
                           const std::string&   a_dependentVarName);

  /// Check positivity of face values
  void checkFacePositivity(const FluxBox&       a_inputFab,
                           const FluxBox&       a_dependentFab,
                           const Box&           a_box,
                           const ProblemDomain& a_domain,
                           const int            a_localComp,
                           const int            a_neighborRadius,
                           const std::string&   a_inputVarName,
                           const std::string&   a_dependentVarName);

  /// Modify the wall shear-stress to match logarithmic law-of-the-wall
  void enforceLawOfTheWallShearStress(
    FluxBox&                a_NGradUfacePntFxb,
    const FluxBox&          a_unitNormalsFxb,
    const FluxBox&          a_WfaceAvgFxb,
    const LevelGridMetrics& a_gridMetrics,
    const ProblemDomain&    a_domain,
    const Box&              a_box,
    const Box&              a_disjointBox);

/*--------------------------------------------------------------------*/
//  Non-linear function for estimating wall friction velocity
/**
 *  \param[in]  a_y     Wall-normal evaluation distance
 *  \param[in]  a_nu    Kinematic viscosity at wall
 *  \param[in]  a_u     van Driest velocity at evaluation point
 *//*-----------------------------------------------------------------*/

  struct UTauSpaldingFunc
  {
    UTauSpaldingFunc(const Real a_y,
                     const Real a_nu,
                     const Real a_u)
      :
      m_y(a_y),
      m_nu(a_nu),
      m_u(a_u)
      {}
    Real operator()(const Real& a_uTau) const
      {
        Real y_plus = m_y*a_uTau/m_nu;
        Real u_plus = std::abs(m_u)/a_uTau;
        Real b = 0.4*u_plus;
        Real b2 = b*b;
        Real b3 = b*b2;
        Real b4 = b*b3;
        Real c_0 = std::exp(b);
        Real c_1 = u_plus + 0.1108*(c_0 - 1. - b - 0.5*b2 - (b3/6.) - (b4/24.));
        return y_plus - c_1;
      }
    const Real m_y;
    const Real m_nu;
    const Real m_u;
  };

}  // namespace CRDutil

#endif /* ! defined _CRDUTIL_H_ */
