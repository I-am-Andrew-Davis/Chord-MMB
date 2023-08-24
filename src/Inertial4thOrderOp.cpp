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
 * \file Inertial4thOrderOp.cpp
 *
 * \brief Member functions for Inertial4thOrderOp
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "FluxBox.H"
#include "MOLUtilFunc.H"
#include "LevelGridMetrics.H"
#include "GodunovUtilitiesF_F.H"
#include "LoHiCenter.H"

//----- Internal -----//

#include "Inertial4thOrderOp.H"
#include "PatchCNSOp.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "CRDPhysics.H"
#include "CRDutil.H"
#include "CNSIBC.H"
#include "DataTemp.H"
#include "LGintegrator.H"  //**FIXME, use CRDPhysics instead
#include "DCFlattening.H"
#include "AMRIO.H"


/*******************************************************************************
 *
 * Class Inertial4thOrderOp: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_levelGridMetrics
 *                      Grid metrics for the level
 *//*-----------------------------------------------------------------*/

Inertial4thOrderOp::Inertial4thOrderOp(LevelGridMetrics& a_levelGridMetrics)
  :
  m_levelGridMetrics(a_levelGridMetrics)
{ }


/*==============================================================================
 * Public member functions
 *============================================================================*/

//**FIXME comments below
/*--------------------------------------------------------------------*/
//  Compute the inertial (hyperbolic) flux
/** \param[in]  a_box   Flux determined on faces of these cells.
 *                      'a_box' must be <= disjoint.
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[in]  a_fluxFaceAvgFace
 *                      Average flux on faces
 *  \param[out] a_fluxFaceAvgFace
 *                      Updated with contribution from inertial terms
 *  \param[out] a_fluxFromWavgFxb
 *                      Second-order estimate of flux
 *  \param[in]  a_WfaceAvgFab
 *                      The average primitive state on the faces.
 *                      This is defined on the disjoint box grown by 2
 *                      except across non-periodic boundaries where it
 *                      is only grown by 0.
 *  \param[out] a_WfacePntFxb
 *                      Point values of the full primtive state, W,
 *                      are set on all faces of 'a_boxWfacePnt'
 *  \param[in]  a_WcellAvgFab
 *                      The average primitive state in the cells.
 *                      This is defined on the disjoint box grown by 4
 *                      except across non-periodic boundaries where it
 *                      is only grown by 2.
 *  \param[out] a_WcellPntFab
 *                      The centered primitive state in the cells.
 *                      This is defined at the same places as a_WcellAvgFab
 *                      except is not needed across non-periodic boundaries.
 *  \param[in]  a_flattening
 *                      Flattening coefficients
 *  \param[out] a_flattening
 *                      Coefficients written if a_setFlattening is T
 *  \param[in]  a_bndryCellFab
 *                      Previous-time data for cells immediately adjacent
 *                      to boundaries
 *  \param[in]  a_bndryNtJ
 *                      NtJ on near-bndry faces necessary for LES wall-model
 *  \param[in]  a_UcellAvgFab
 *                      The average solution state in the cells in
 *                      physical space
 *  \param[in]  a_unitNormalFxb
 *                      Unit normals for transforming velocity space
 *                      to be normal to the faces
 *  \param[in]  a_dataIndx
 *                      Current data index
 *  \param[in]  a_time  Time (for applying BC)
 *  \param[in]  a_dt    Time-step size for current time-step and level
 *  \param[in]  a_prevDt
 *                      Time-step size for past time-step
 *  \param[in]  a_stageWeight
 *                      Weighting applied to contribution of flux
 *                      from this stage of RK4 to total flux
 *  \param[in]  a_setFlattening
 *                      T - Set flattening coefficients
 *  \param[in]  a_level Index of the AMR level
 *
 *  \note
 *  <ul>
 *    <li> Previous implementations gave the option to do
 *         convolutions.  Here, all convolutions are always performed
 *         (relates to pnt <-> avg values).
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
Inertial4thOrderOp::flux(const Box&                          a_disjointBox,
                         const Box&                          a_box,
                         const BlockDomain&                  a_domain,
                         FluxBox&                            a_fluxFaceAvgFxb,
                         FluxBox&                            a_fluxFromWavgFxb,
                         FluxBox&                            a_WfaceAvgFxb,
                         FluxBox&                            a_WfacePntFxb,
                         FArrayBox&                          a_WcellAvgFab,
                         const FArrayBox&                    a_WcellPntFab,
                         FluxBox&                            a_faceAvgPlotFxb,
                         FArrayBox&                          a_flattening,
                         stc::Vector<FArrayBox, 2*SpaceDim>& a_bndryCellFab,
                         stc::Vector<FArrayBox, 2*SpaceDim>& a_bndryNtJ,
                         const FArrayBox&                    a_UcellAvgFab,
                         const FluxBox&                      a_unitNormalFxb,
                         const DataIndex&                    a_dataIndx,
                         const Real                          a_time,
                         const Real                          a_dt,
                         const Real                          a_prevDt,
                         const Real                          a_stageWeight,
                         const bool                          a_setFlattening,
                         const int                           a_level) const
{
#ifndef NDEBUG
  // Face WpAvg - interpolate native primitive state and limit.  Together, these
  //   require 3 cells normal to the face
  // Note: CRDparam::NumGhostWfaceAvgNrm is used for number of all remaining
  //   stencil cells in the normal direction on faces for W state
  // Note: NumGhostWpFaceAvgTan == NumGhostWfaceAvgTan
  // Regarding the following ghosts, we generally need a full computation in
  // a cell normal to and outside the box.  I.e., the state must be interpolated
  // on both of of the faces.
  // For flattening, we need 1 plus half the flattening stencil.
  if (CRDparam::g_useFlattening)
    {
      CH_assert(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm) + 4 <=
      CRDparam::queryGhosts(CRDparam::NumGhostWfromUcellAvg));
    }
  // The number of ghosts needed in WcellAvg depend on the order of the face
  // interpolation.  In general, it is 1 plus half the stencil.  But see
  // exception below
  switch (CRDparam::g_faceInterpolationOrder)
    {
    case 1:
      CH_assert(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm) + 2 <=
                CRDparam::queryGhosts(CRDparam::NumGhostWcellAvg));
      break;
    case 2:
      CH_assert(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm) + 2 <=
                CRDparam::queryGhosts(CRDparam::NumGhostWcellAvg));
    case 4:
      if (CRDparam::g_limitFaceValues)
        {
          CH_assert(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm) + 4 <=
                    CRDparam::queryGhosts(CRDparam::NumGhostWcellAvg));
        }
      else
        {
          if (CRDparam::g_usePPMlimiter)
            {
              CH_assert(
                CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm) + 3 <=
                CRDparam::queryGhosts(CRDparam::NumGhostWcellAvg));
            }
          else
            {
              CH_assert(
                CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm) + 2 <=
                CRDparam::queryGhosts(CRDparam::NumGhostWcellAvg));
            }
        }
      break;
    case 5:
      CH_assert(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm) + 4 <=
                CRDparam::queryGhosts(CRDparam::NumGhostWcellAvg));
      break;
    }
#endif
  CH_assert(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgTan) <=
            CRDparam::queryGhosts(CRDparam::NumGhostWcellAvg));
  Box boxWpFaceAvgDom =
    grow(a_box,
         std::max(
           CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry),
           std::max(CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm),
                    CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgTan))));
  boxWpFaceAvgDom &= a_domain;

  // Face WpPnt - deconvolve WpAvg (requires 1 (normal) or 2 (limited) cells in
  //   tangential directions)
  // Note: NumGhostWpFacePntTan == NumGhostWfaceAvgTan
  CH_assert(CRDparam::queryGhosts(CRDparam::NumGhostWpFacePntTan) + 1 <=
            CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgTan));
  Box boxWpFacePntDom =
    grow(a_box,
         std::max(
           CRDparam::queryGhosts(CRDparam::NumGhostWfacePntTngBdry),
           std::max(CRDparam::queryGhosts(CRDparam::NumGhostWpFacePntNrm),
                    CRDparam::queryGhosts(CRDparam::NumGhostWpFacePntTan))));
  boxWpFacePntDom &= a_domain;

  // Face Wpnt - compute the full primitive state
  // Note: NumGhostWFacePntTan == NumGhostWfaceAvgTan
  const Box& boxWfacePntDom = boxWpFacePntDom;
  CH_assert(a_WfacePntFxb.box().contains(boxWfacePntDom));

  // Face Wavg (this is only used for viscous updates) - convolve Wpnt using
  // WpAvg, and extra primitives computed directly from that, for derivatives
  const Box& boxWfaceAvgDom = boxWpFacePntDom;  // Rename for clarity
  CH_assert(a_WfaceAvgFxb.box().contains(boxWfaceAvgDom));

  // Face Fpnt - compute the flux
  // Note: No faces are required in the normal direction for fluxes.  Suffix Tan
  //       is dropped.
  // Note: NumGhostInertialFfacePnt == NumGhostWfaceAvgTan
  Box boxFfacePntDom =
    grow(a_box, CRDparam::queryGhosts(CRDparam::NumGhostInertialFfacePnt));
  boxFfacePntDom &= a_domain;

  // Face Favg.  Convolve Fpnt but there are two options regarding the
  //   derivatives use both here and for the product rule Nt*Favg:
  //   a) Derivatives are always computed from FfromWpAvg.  This is used in the
  //      inertial-only solutions.
  //   b) Derivatives are computed from the respective state, e.g., Fpnt for
  //      the convolution and Favg for the product rule.  An extra ghost cell
  //      is required in each step.
  //   The number of tangential ghost cells required for Favg determines which
  //   approach is used: 0 for a) and 1 for b)

  // Face FfromWavg
  Box boxFfromWpAvgDom;  // Empty
  if (CRDparam::queryGhosts(CRDparam::NumGhostInertialFfaceAvg) == 0)
    {
      CH_assert(CRDparam::queryGhosts(CRDparam::NumGhostInertialFfaceAvg) + 1 <=
                CRDparam::queryGhosts(CRDparam::NumGhostFfromWpFaceAvg));
      // Recompute to neglect any extra cells in normal direction
      boxFfromWpAvgDom =
        grow(a_box, CRDparam::queryGhosts(CRDparam::NumGhostFfromWpFaceAvg));
      boxFfromWpAvgDom &= a_domain;
    }
  else
    {
      CH_assert(CRDparam::queryGhosts(CRDparam::NumGhostInertialFfaceAvg) + 1 <=
                CRDparam::queryGhosts(CRDparam::NumGhostInertialFfacePnt));
    }

  // Face Favg
  Box boxFfaceAvgDom =
    grow(a_box, CRDparam::queryGhosts(CRDparam::NumGhostInertialFfaceAvg));
  boxFfaceAvgDom &= a_domain;

  // Cells needed outside and in a tangential direction to a boundary
  // If ghost cells are not needed, set it equal to 
  Box boxWcellTngBdryDom =
    grow(a_box, std::max(
           CRDparam::queryGhosts(CRDparam::NumGhostWcellExtTngBdry),
           CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry)));
  boxWcellTngBdryDom &= grow(a_domain, 2);  // Restrict normal to Phys. Bound.
  CH_assert(a_WcellAvgFab.box().contains(boxWcellTngBdryDom));

//--Numbers of components  

  const int numUcomp  = CRDparam::g_CRDPhysics->numConservative();
  const int numWcomp  = CRDparam::g_CRDPhysics->numPrimitive();
  const int numFcomp  = CRDparam::g_CRDPhysics->numFluxes();
  CH_assert(a_UcellAvgFab.nComp() == numUcomp);
  CH_assert(a_WfaceAvgFxb.nComp() == numWcomp);
  CH_assert(a_WfacePntFxb.nComp() == numWcomp);
  CH_assert(a_WcellAvgFab.nComp() == numWcomp);
  CH_assert(a_fluxFaceAvgFxb.nComp() == SpaceDim*numFcomp);
  CH_assert(a_fluxFromWavgFxb.nComp() == SpaceDim*numFcomp);

//-- Start of work -------------------------------------------------------------

  // bool setFlattening = true;

  // Compute averages of W on faces
  // ---------------------------------------
  //  I/O           Var            Stencil  
  // =====  ====================  ==========
  // read   a_WcellAvgFab          8 per dir  //**FIXME this is either 6/8 or
  // read   WfromUavgFab          10 per dir  //**      8/10 depending on
  // read   a_flattening                   1  //**      PPMLimiter interface
  // write  a_flattening                   1
  // write  WpFaceAvgFxb                   1
  // ---------------------------------------
  // Fourth-order interpolation (if used):
  // WfaceAvg[i+e/2] = 7/12 * (a_WcellAvgFab[i] + a_WcellAvgFab[i+e])
  //                 - 1/12 * (a_WcellAvgFab[i-e] + a_WcellAvgFab[i+2e]).
  // Includes calls to limiter (if used) and flattening (if used),
  // and Riemann solver.
  // computeFaceAverage(a_box,
  //                    a_domain,
  //                    a_flattening,
  //                    a_WfaceAvgFxb,
  //                    a_WcellAvgFab,
  //                    a_WfromUavgFab,
  //                    a_unitNormalFxb,
  //                    a_time,
  //                    setFlattening,  //**FIXME or a_setFlattening?
  //                    a_level);

  // Compute point values of W, the full primitive state of W and <W>, and
  // apply BC
  computeFullFaceState(boxWfacePntDom,
                       boxWpFaceAvgDom,
                       boxWcellTngBdryDom,
                       a_box,
                       a_domain,
                       a_WcellAvgFab,
                       a_bndryCellFab,
                       a_bndryNtJ,
                       a_WfaceAvgFxb,
                       a_WfacePntFxb,
                       a_unitNormalFxb,
                       a_dataIndx,
                       a_time,
                       a_prevDt,
                       a_level);

  // Return if only BC are required
  // This conditional can be ignored for performance (e.g., loop-chain) work
  if (!(CRDparam::g_physicsModels & CRDparam::PhysicsInertial))
    {
      return;
    }

  // Compute the fluxes.  This include point fluxes on the faces computed
  // from WfacePntFxb and an approximation from WfaceAvgFxb (the latter is
  // only required to reduce the stencil if viscous terms are not included).
  // ---------------------------------------
  //  I/O           Var            Stencil  
  // =====  ====================  ==========
  // read   WfacePntFxb                    1
  // read   WfaceAvgFxb                    1
  // write  fluxFacePntFxb                 1
  // write  a_fluxFromWavgFxb              1
  // ---------------------------------------
  FLUXBOXSTACKTEMP(fluxFacePntFxb, boxFfacePntDom, SpaceDim*numFcomp);
  // Note: NumGhostInertialFfacePnt == NumGhostWfaceAvgTan
  const int numGhostFfacePnt =
    CRDparam::queryGhosts(CRDparam::NumGhostInertialFfacePnt);
  // Compute point values of fluxes
  computeAllFluxes(a_box,
                   a_domain,
                   fluxFacePntFxb,
                   a_WfacePntFxb,
                   numGhostFfacePnt);
  // Compute a flux from <W>
  if (!boxFfromWpAvgDom.isEmpty())// && (faceConvOrder == 4))
    {
      // Note: NumGhostFfromWpFaceAvg == NumGhostWpFaceAvgTan
      const int numGhostFfromWpFaceAvg =
        CRDparam::queryGhosts(CRDparam::NumGhostFfromWpFaceAvg);
      computeAllFluxes(a_box,
                       a_domain,
                       a_fluxFromWavgFxb,
                       a_WfaceAvgFxb,
                       numGhostFfromWpFaceAvg);
    }

  // Compute the average flux on the faces using Laplacian
  // ---------------------------------------------------
  //  I/O           Var            Stencil      Ghosts
  // =====  ====================  ==========  ==========
  // read   fluxFacePntFxb                 1           0
  // read   a_fluxFromWavgFxb     3 per tdir    1 in dom
  // write  a_fluxFaceAvgFxb               1           0
  // ---------------------------------------------------
  //**FIXME costly to copy and then perform operation.  Do the copy in the Op.
  a_fluxFaceAvgFxb.copy(fluxFacePntFxb, boxFfaceAvgDom);
  // Depending on number of ghosts available, may use Fpnt or FfromWAvg to
  // compute gradients of fluxes.
  FluxBox* fluxGrad =
    (CRDparam::queryGhosts(CRDparam::NumGhostInertialFfaceAvg) == 0) ?
    &a_fluxFromWavgFxb : &fluxFacePntFxb;
  if (CRDparam::g_faceConvolveFlatten < 2)
    {
      //**FIXME: are we deconvolving 1 layer of face-normal ghosts?
      MOLUtilFunc::deconvolveFace(a_fluxFaceAvgFxb,
                                  *fluxGrad,
                                  boxFfaceAvgDom,
                                  a_domain,
                                  1);
    }

  //**NOTE: We'll move this somewhere else soon
  // Plot squared variables (for RMS and stress calculations)
  int plotVariables = CRDparam::g_plotLoFaceAvgComps;
  plotVariables |= CRDparam::g_plotLoFaceAvgTimeAvgComps;
  if ((plotVariables & CRDparam::PlotTurbulentComps) &&
      (CNSIBC::s_firstRKStage))
    {
      int cStart = 0; // Yes, we'll use indexing within a plot class or
                      // another improved method --- as it is, this is annoying
      if (plotVariables & CRDparam::PlotPrimitive)
        {
          cStart = CRDparam::g_CRDPhysics->numPrimitive();
        }
      if ((plotVariables & CRDparam::PlotFluxes) &&
          (CRDparam::g_physicsModels & CRDparam::PhysicsInertial))
        {
          cStart += CRDparam::g_CRDPhysics->numFluxes();
        }
      if ((plotVariables & CRDparam::PlotFluxes) &&
          (CRDparam::g_physicsModels & CRDparam::PhysicsViscous))
        {
          cStart += CRDparam::g_CRDPhysics->numFluxes();
        }

      // Primitive state indices
      const int cRho  = CRDparam::g_CRDPhysics->densityIndex();
      const int cVel = CRDparam::g_CRDPhysics->velocityInterval().begin();
      const int cPress = CRDparam::g_CRDPhysics->pressureIndex();
      const int cTemp = CRDparam::g_CRDPhysics->temperatureIndex();
      // Compute face-avg state on a_box --- need pnt state on box1Dom
      Box box1Dom = grow(a_box, 1);
      box1Dom &= a_domain;
      int numTurbInviscidComps = 3 + SpaceDim*(SpaceDim - 1.)/2. + SpaceDim;
      FLUXBOXSTACKTEMP(turbFacePntFxb, box1Dom, numTurbInviscidComps);
      FLUXBOXSTACKTEMP(turbFaceAvgFxb, a_box, numTurbInviscidComps);
      for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
        {
          IntVect growVect = IntVect_unit;
          growVect[faceDir] = 0;
          Box faceBox = a_box;
          faceBox.grow(growVect); // Grow in tangential directions only
          faceBox &= a_domain;
          faceBox.surroundingNodes(faceDir);
          int cLoc = 0;
          // Compute <u_i*u_j>
          FArrayBox& turbFacePntFab = turbFacePntFxb[faceDir];
          FArrayBox& WfacePntFab = a_WfacePntFxb[faceDir];
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              for (int sComp = dir; sComp != SpaceDim; ++sComp)
                {
                  MD_BOXLOOP(faceBox, i)
                    {
                      Real uValOne = WfacePntFab[MD_IX(i, cVel+dir)];
                      Real uValTwo = WfacePntFab[MD_IX(i, cVel+sComp)];
                      turbFacePntFab[MD_IX(i,cLoc)] = uValOne*uValTwo;
                    }
                  ++cLoc;
                }
            }
          // Compute <p*p>
          MD_BOXLOOP(faceBox, i)
            {
              Real pressure = WfacePntFab[MD_IX(i, cPress)];
              turbFacePntFab[MD_IX(i,cLoc)] = pressure*pressure;
            }
          ++cLoc;
          // Compute <T*T>
          MD_BOXLOOP(faceBox, i)
            {
              Real temperature = WfacePntFab[MD_IX(i, cTemp)];
              turbFacePntFab[MD_IX(i,cLoc)] = temperature*temperature;
            }
          ++cLoc;
          // Compute <rho*rho>
          MD_BOXLOOP(faceBox, i)
            {
              Real density = WfacePntFab[MD_IX(i, cRho)];
              turbFacePntFab[MD_IX(i,cLoc)] = density*density;
            }
        }
      // Convolve the face-centered data
      turbFaceAvgFxb.copy(turbFacePntFxb, a_box);
      if (CRDparam::g_faceConvolveFlatten < 2)
        {
          MOLUtilFunc::deconvolveFace(
            turbFaceAvgFxb, turbFacePntFxb, a_box, a_domain, 1);
        }
      // Copy turbFaceAvgFxb into a_faceAvgPlotFxb
      a_faceAvgPlotFxb.copy(turbFaceAvgFxb, 0, cStart, numTurbInviscidComps);
    }

  //**This is only required for linear problems (fifth-order upwind is probably
  //**a better approach)
  // if (m_useArtificialDissipation)
  //   {
  //     // Add fifth-derivative correction to flux
  //     fluxCorrection(a_FfaceAvg, a_UcellAvgFab);
  //   }
}


/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Compute averages of W on faces
/** Sounds easy but this includes all stabilization methods including
 *  PPM, slope flattening, and solving the Riemann problem
 *  \param[out] a_boxWfaceAvgDom
 *                      Average values of \<W\> are computed on the
 *                      faces of this cell-centered box.  In each
 *                      direction, we take the surrounding nodes of
 *                      this box and then shrink by 1 in 'dir'.  I.e.,
 *                      there is one extra ghost in tangential
 *                      directions.
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[in]  a_flattening
 *                      Flattening coefficients
 *  \param[out] a_flattening
 *                      Coefficients written if a_setFlattening is T
 *  \param[out] a_WfaceAvgFxb
 *                      The average native (minimal) primitive state,
 *                      \<W\>, on faces.  'a_WfaceAvgFxb' should be
 *                      aliased to only have the native number of
 *                      primitive states.
 *  \param[in]  a_WcellAvgFab
 *                      The average native (minimal) primitive state,
 *                      \<W\>, in cells.  'a_WcellAvgFxb' should be
 *                      aliased to only have the native number of
 *                      primitive states.
 *  \param[in]  a_WfromUavgFab
 *                      The average native (minimal) primitive state,
 *                      \<W\>, in cells computed directly from \<U\>.
 *                      It is only second-order accurage.
 *                      'a_WfromUavgFab' should be aliased to only
 *                      have the native number of primitive states.
 *  \param[in]  a_unitNormalFxb
 *                      Unit normal basis for transforming velocity
 *                      space to be normal to the faces
 *  \param[in]  a_time  Time (for applying BC)
 *  \param[in]  a_dt    Time step size
 *  \param[in]  a_setFlattening
 *                      T - Set flattening coefficients
 *  \param[in]  a_level Index of the AMR level
 *
 *  \note
 *  <ul>
 *    <li> \<W\> is determined on 1 less face in normal directions
 *         This is so tangential gradients can be taken to compute
 *         point values of W on the faces.  In other words,
 *         a_boxWfacePntDom = grow(a_boxWfaceAvgDom, -1).
 *    <li> All arrays of W should only have the native (minimal)
 *         number of primitive data.
 *  </ul>
 *//*-----------------------------------------------------------------*/

// ---------------------------------------------------
//  I/O           Var            Stencil      Ghosts
// =====  ====================  ==========  ==========
// read   a_WcellAvgFab          4 per dir    6 in dom
// read   a_WfromUavgFab         7 per dir    6 in dom
// read   a_flattening                   1    3 in dom
// write  a_flattening                   1    3 in dom
// write  a_WfaceAvgFxb                  1    3 in dom
// ---------------------------------------------------
// Note: write to a_flattening only occurs if
// a_setFlattening is true.
// ---------------------------------------------------

void
Inertial4thOrderOp::computeFaceAverage(
  const Box&         a_box,
  const BlockDomain& a_domain,
  FArrayBox&         a_flattening,
  FluxBox&           a_WfaceAvgFxb,   // (a_faceW)
  const FArrayBox&   a_WcellAvgFab,   // (a_cellW)
  const FArrayBox&   a_WfromUavgFab,  // (a_WofUavg)
  const FluxBox&     a_unitNormalFxb,
  const Real         a_time,
  const bool         a_setFlattening,
  const int          a_level) const
{
  CH_TIME("Inertial4thOrderOp::computeFaceAverage");
  // Number of native (minimal) primitive states (same as F).  Ghost cells and
  // boxes are carefully named with Wp.  FArrayBoxes and FluxBoxes are generally
  // named with W but only Wp components are used.
  int numWcomp = CRDparam::g_CRDPhysics->numNativePrimitive();
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      numWcomp = CRDparam::g_CRDPhysics->numPrimitive();
    }

  // Number of normal faces on ghost cells to solve
  // Note: NumGhostWpFaceAvgNrm == NumGhostWfaceAvgNrm
  const int numGhostWpFaceAvgNrm =
    CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm);
  int numGhostWpFaceAvgTan =
    CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgTan);

  // Number of tangential faces on ghost cells to solve
  // If we are adjacent to any boundary, and there is at least one periodic
  // boundary, ensure numGhostWpFaceAvgTan supports NumGhostWcellTngBdry + 1
  if (!a_domain.contains(grow(a_box, 1)))// && a_domain.isPeriodic())
    {
      numGhostWpFaceAvgTan = std::max(
        numGhostWpFaceAvgTan,
        CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry));
      CH_assert(numGhostWpFaceAvgTan <=
                CRDparam::queryGhosts(CRDparam::NumGhostWcellAvg));
    }

  // Set flags for HO limiting
  // If no 3rd derivative checks are allowed
  // int loBoundDerivChecks = 1;
  // int hiBoundDerivChecks = 1;
  // if (CRDparam::g_noHOchecks)
  //   {
  //     HOchecks = 0;
  //     loBoundDerivChecks = 0;
  //     hiBoundDerivChecks = 0;
  //   }
  // // If no 3rd derivative checks near boundaries
  // // FIXME: These are split into upper and lower to allow for selective
  // // boundary adjustments
  // int loBoundDerivChecks = 1;
  // int hiBoundDerivChecks = 1;
  // if (CRDparam::g_extraBoundLim)
  //   {
  //     loBoundDerivChecks = 0;
  //     hiBoundDerivChecks = 0;
  //   }
  bool useFCOR = false;
  if (CRDparam::g_useFCOR)
    {
      useFCOR = true;
    }

  // Set flattening for all slopes if needed (done once per time-step).
  // Each CELL has a flattening coefficient, eta.
  // Flattening means replacing Wface on left and right of each cell
  // by eta*Wface + (1-eta)*Wcell.
  // = Wcell + eta * (Wface - Wcell).
  // That is:
  // WfaceL[+] = Wcell + eta * (WfaceL[+] - Wcell);
  // WfaceR[-] = Wcell + eta * (WfaceR[-] - Wcell).
  // We use a_WfromUAvgFab instead of a_WcellAvgFab because at second order,
  // it provides a more stable state to evaluate the physical, versus numerical,
  // conditions where flattening should be applied.  E.g., we are looking for a
  // steep shock, not an oscillation.
  if (CRDparam::g_useFlattening && a_setFlattening)
    {
      Interval velIntv   = CRDparam::g_CRDPhysics->velocityInterval();
      int pressureIdx    = CRDparam::g_CRDPhysics->pressureIndex();
      Real smallPressure = CRDparam::g_smallp;
      int bulkIdx        = CRDparam::g_CRDPhysics->bulkModulusIndex();
      // The flattening box is based off of NumGhostWpFaceAvgTan and
      // NumGhostWpFaceAvgNrm == NumGhostWfaceAvgNrm.  Tangential directions
      // are sufficient.  In the normal direction, grow by 1 for the low/high
      // face values "outside" the box.  Remember that this covers all
      // directions so we need a maximum.
      const int flattenGhost = std::max(numGhostWpFaceAvgNrm + 1,
                                        numGhostWpFaceAvgTan);
      Box flattenBox = grow(a_box, flattenGhost);
      flattenBox &= a_domain;
      CH_assert(a_flattening.box().contains(flattenBox));
      Box testFlattenBox = grow(flattenBox, 3);
      testFlattenBox &= a_domain;
      CH_assert(a_WfromUavgFab.box().contains(testFlattenBox));
      // ---------------------------------------
      //  I/O           Var            Stencil  
      // =====  ====================  ==========
      // read   a_WfromUAvgFab         7 per dir
      // write  a_flattening                   1
      // ---------------------------------------
      MOLUtilFunc::computeFlattening(a_flattening,
                                     a_WfromUavgFab,
                                     velIntv,
                                     pressureIdx,
                                     smallPressure,
                                     bulkIdx,
                                     flattenBox,
                                     a_domain);
    }

//--Primary loop over face directions

  for (const int dir : EachDir)
    {
      FArrayBox& WfaceAvgDirFab = a_WfaceAvgFxb[dir];

      // This is the box at which we ultimately need WpFaceAvg.
      IntVect growVect(numGhostWpFaceAvgTan*IntVect::Unit);
      growVect[dir] = numGhostWpFaceAvgNrm;
      Box cellBoxWpFaceAvgDom = grow(a_box, growVect);
      cellBoxWpFaceAvgDom &= a_domain;

      // Wminus and Wplus store +/- adjacent face values in cells.  Grow by 1
      // in normal direction to store outer face values and crop at domain.
      // This is similar to flattenBox except that it is directional.
      // This is also the box on which we will use the interior scheme to
      // interpolate face average values at the interior (see
      // faceBoxWpFaceAvgLRDom)
      Box cellBoxWpFaceAvgLRdom = cellBoxWpFaceAvgDom;
      if (CRDparam::g_usePPMlimiter ||
          (CRDparam::g_faceInterpolationOrder % 2) == 1)
        {
          cellBoxWpFaceAvgLRdom.grow(dir, 1);
        }
      cellBoxWpFaceAvgLRdom &= a_domain;
      Box eTolBox = cellBoxWpFaceAvgDom;
      eTolBox.grow(dir, 2);
      eTolBox &= a_domain;
      // PPM limiting needs L and R face values and we need to apply limiting
      // in one ghost cell in the normal direction to get the limited value
      // in the outside face.  The also required the opposing face so
      // faceBoxWpFaceAvgLRDom grows by 1 in the normal direction.
      Box faceBoxWpFaceAvgLRdom = surroundingNodes(cellBoxWpFaceAvgLRdom, dir);

      FABSTACKTEMP(Wminus, cellBoxWpFaceAvgLRdom, numWcomp);
      FABSTACKTEMP(Wplus, cellBoxWpFaceAvgLRdom, numWcomp);

//--Face interpolation (with possibly face limiting)

      switch (CRDparam::g_faceInterpolationOrder)
        {
        case 1:
          // Initially set Wminus and Wplus to the cell-average state
          // ---------------------------------------
          //  I/O           Var            Stencil  
          // =====  ====================  ==========
          // read   WfaceAvgDirFab                 2
          // write  Wminus                         1
          // write  Wplus                          1
          // ---------------------------------------
          Wminus.copy(a_WcellAvgFab,
                      cellBoxWpFaceAvgLRdom, 0,
                      cellBoxWpFaceAvgLRdom, 0,
                      numWcomp);
          Wplus .copy(a_WcellAvgFab,
                      cellBoxWpFaceAvgLRdom, 0,
                      cellBoxWpFaceAvgLRdom, 0,
                      numWcomp);
          break;
        case 2:
          // ---------------------------------------
          //  I/O           Var            Stencil  
          // =====  ====================  ==========
          // read   WfaceAvgDirFab                 2
          // write  Wminus                         1
          // write  Wplus                          1
          // ---------------------------------------
          solveSecondOrderVals(WfaceAvgDirFab,
                               a_WcellAvgFab,
                               numWcomp,
                               dir,
                               faceBoxWpFaceAvgLRdom,
                               a_domain);
          WfaceAvgDirFab.shiftHalf(dir, +1); // [i-e/2] becomes [i]
          Wminus.copy(WfaceAvgDirFab,
                      cellBoxWpFaceAvgLRdom, 0,
                      cellBoxWpFaceAvgLRdom, 0,
                      numWcomp);
          WfaceAvgDirFab.shift(dir, -1);     // [i+e/2] becomes [i]
          Wplus.copy( WfaceAvgDirFab,
                      cellBoxWpFaceAvgLRdom, 0,
                      cellBoxWpFaceAvgLRdom, 0,
                      numWcomp);
          WfaceAvgDirFab.shiftHalf(dir, +1); // final shift back
          break;
        case 4:
          if (CRDparam::g_limitFaceValues)
            {
              // FAB containing the smoothness tests
              FABSTACKTEMP(smoothTestFab, eTolBox, numWcomp);
              smoothTestFab.setVal(0.);
              // Perform the smoothness test
              // CRDutil::smoothTest(smoothTestFab,
              //                     a_WcellAvgFab,
              //                     numWcomp,
              //                     dir,
              //                     eTolBox,
              //                     a_domain);
              // ---------------------------------------
              //  I/O           Var            Stencil  
              // =====  ====================  ==========
              // read   a_WcellAvgFab                  6
              // write  WfaceAvgDirFab                 1
              // ---------------------------------------
              CRDutil::PPMFaceValues(WfaceAvgDirFab,
                                     a_WcellAvgFab,
                                     smoothTestFab,
                                     numWcomp,
                                     true,  // Not used
                                     dir,
                                     faceBoxWpFaceAvgLRdom,
                                     a_domain);
            }
          else
            {
              // Set WfaceAvgDirFab from a_WcellAvgFab by 4th-order
              // interpolation on all dir-faces of box1Dom with the formula
              // WfaceAvgDirFab[i+e/2] =
              //     7/12 * (a_WcellAvgFab[i  ] + a_WcellAvgFab[i+e ])
              //   - 1/12 * (a_WcellAvgFab[i-e] + a_WcellAvgFab[i+2e]).
              // ---------------------------------------
              //  I/O           Var            Stencil  
              // =====  ====================  ==========
              // read   a_WcellAvgFab                  4
              // write  WfaceAvgDirFab                 1
              // ---------------------------------------
              MOLUtilFunc::PPMFaceValues(WfaceAvgDirFab,
                                         a_WcellAvgFab,
                                         numWcomp,
                                         dir,
                                         faceBoxWpFaceAvgLRdom,
                                         a_domain);
            }
          // Initially set Wminus and Wplus to the face-average state
          // ---------------------------------------
          //  I/O           Var            Stencil  
          // =====  ====================  ==========
          // read   WfaceAvgDirFab                 2
          // write  Wminus                         1
          // write  Wplus                          1
          // ---------------------------------------
          WfaceAvgDirFab.shiftHalf(dir, +1); // [i-e/2] becomes [i]
          Wminus.copy(WfaceAvgDirFab,
                      cellBoxWpFaceAvgLRdom, 0,
                      cellBoxWpFaceAvgLRdom, 0,
                      numWcomp);
          WfaceAvgDirFab.shift(dir, -1);     // [i+e/2] becomes [i]
          Wplus.copy( WfaceAvgDirFab,
                      cellBoxWpFaceAvgLRdom, 0,
                      cellBoxWpFaceAvgLRdom, 0,
                      numWcomp);
          WfaceAvgDirFab.shiftHalf(dir, +1); // final shift back
          break;
        case 5:
          if (CRDparam::g_limitFaceValues)
            {
              // ---------------------------------------
              //  I/O           Var            Stencil  
              // =====  ====================  ==========
              // read   a_WcellAvgFab                  6
              // write  Wminus                         1
              // write  Wplus                          1
              // ---------------------------------------
              CRDutil::PPMUpwindFaceValues(Wminus,
                                           Wplus,
                                           a_WcellAvgFab,
                                           numWcomp,
                                           dir,
                                           cellBoxWpFaceAvgLRdom,
                                           a_domain);
            }
          else
            {
//          CH_assert(false);  //**FIXME
              // ---------------------------------------
              //  I/O           Var            Stencil  
              // =====  ====================  ==========
              // read   a_WcellAvgFab                  6
              // write  WfaceAvgDirFab                 1
              // ---------------------------------------
              solveFifthOrderVals(Wminus,
                                  Wplus,
                                  a_WcellAvgFab,
                                  numWcomp,
                                  dir,
                                  cellBoxWpFaceAvgLRdom,
                                  a_domain);
            }
          break;
        }

//--PPM cell limiting and flattening

      if (CRDparam::g_usePPMlimiter)
        {
          // Second-order plus and minus states
          FABSTACKTEMP(WminusLO, cellBoxWpFaceAvgLRdom, numWcomp);
          FABSTACKTEMP(WplusLO, cellBoxWpFaceAvgLRdom, numWcomp);
          if (useFCOR)
            {
              // Compute second-order face values and apply PPM limiter
              FABSTACKTEMP(W2ndFace, faceBoxWpFaceAvgLRdom, numWcomp);
              solveSecondOrderVals(W2ndFace, a_WfromUavgFab, numWcomp, dir,
                                   faceBoxWpFaceAvgLRdom, a_domain);
              W2ndFace.shiftHalf(dir, +1); // [i-e/2] becomes [i]
              WminusLO.copy(W2ndFace, cellBoxWpFaceAvgLRdom, 0,
                            cellBoxWpFaceAvgLRdom, 0, numWcomp);
              W2ndFace.shift(dir, -1);     // [i+e/2] becomes [i]
              WplusLO.copy(W2ndFace, cellBoxWpFaceAvgLRdom, 0,
                           cellBoxWpFaceAvgLRdom, 0, numWcomp);
              W2ndFace.shiftHalf(dir, +1); // final shift back
              WminusLO.minus(a_WfromUavgFab, cellBoxWpFaceAvgLRdom,
                             0, 0, numWcomp);
              WplusLO .minus(a_WfromUavgFab, cellBoxWpFaceAvgLRdom,
                             0, 0, numWcomp);
              CRDutil::PPMLimiter(WminusLO,
                                  WplusLO,
                                  a_WfromUavgFab,
                                  numWcomp,
                                  dir,
                                  cellBoxWpFaceAvgLRdom,
                                  a_domain);
              // MOLUtilFunc::PPMLimiter(WminusLO,
              //                         WplusLO,
              //                         a_WfromUavgFab,
              //                         numWcomp,
              //                         dir,
              //                         cellBoxWpFaceAvgLRdom,
              //                         a_domain,
              //                         HOchecks,
              //                         loBoundDerivChecks,
              //                         hiBoundDerivChecks);
            }
          // Wminus -= a_cellW;
          // Wplus  -= a_cellW;
          // Call minus() instead of the simpler -= because the domains do not
          // match.
          Wminus.minus(a_WcellAvgFab,
                       cellBoxWpFaceAvgLRdom,
                       0, 0, numWcomp);
          Wplus .minus(a_WcellAvgFab,
                       cellBoxWpFaceAvgLRdom,
                       0, 0, numWcomp);

          // Apply limiter to Wminus and Wplus
          // ---------------------------------------
          //  I/O           Var            Stencil  
          // =====  ====================  ==========
          // read   a_WcellAvgFab                  7
          // read   Wminus                         1
          // read   Wplus                          1
          // write  Wminus                         1
          // write  Wplus                          1
          // ---------------------------------------
          // Recall:
          // Wminus[i] = thisFaceWDir[i - e/2] - a_cellW.
          // Wplus[i]  = thisFaceWDir[i + e/2] - a_cellW.
          CRDutil::PPMLimiter(Wminus,
                              Wplus,
                              a_WcellAvgFab,
                              numWcomp,
                              dir,
                              // Use the box where Wminus and Wplus are
                              // defined
                              cellBoxWpFaceAvgLRdom,
                              a_domain);
          // MOLUtilFunc::PPMLimiter(Wminus,
          //                         Wplus,
          //                         a_WcellAvgFab,
          //                         numWcomp,
          //                         dir,
          //                         // Use the box where Wminus and Wplus are
          //                         // defined
          //                         cellBoxWpFaceAvgLRdom,
          //                         a_domain,
          //                         HOchecks,
          //                         loBoundDerivChecks,
          //                         hiBoundDerivChecks);
          if (useFCOR)
            {
              CRDparam::g_DCF->ppmFCOR(
                Wminus, Wplus, a_WcellAvgFab, a_WfromUavgFab, WminusLO, WplusLO,
                numWcomp, dir, cellBoxWpFaceAvgLRdom, a_domain);
            }
          if (CRDparam::g_useFlattening)
            {
              // Apply flattening to Wminus and Wplus
              // ---------------------------------------
              //  I/O           Var            Stencil  
              // =====  ====================  ==========
              // read   a_flattening                   1
              // read   Wminus                         1
              // read   Wplus                          1
              // write  Wminus                         1
              // write  Wplus                          1
              // ---------------------------------------
              CH_assert(a_flattening.box().contains(cellBoxWpFaceAvgLRdom));
              MOLUtilFunc::applyFlattening(Wminus,
                                           a_flattening,
                                           cellBoxWpFaceAvgLRdom);
              MOLUtilFunc::applyFlattening(Wplus,
                                           a_flattening,
                                           cellBoxWpFaceAvgLRdom);
            }
          if (CRDparam::g_extraBoundLim)
            {
              if (!a_domain.contains(
                    grow(cellBoxWpFaceAvgLRdom, IntVect_basis(dir))))
                {
                  CRDutil::extraBoundaryLimiting(Wminus,
                                                 Wplus,
                                                 a_WcellAvgFab,
                                                 cellBoxWpFaceAvgLRdom,
                                                 numWcomp,
                                                 dir,
                                                 a_level,
                                                 a_time,
                                                 m_levelGridMetrics,
                                                 a_domain);
                }
            }

          // Wminus[i] = WfaceAvgDirFab[i-e/2] on right after limiting
          //          Wminus += a_cellW;
          // Wplus[i] = WfaceAvgDirFab[i+e/2] on left after limiting
          //          Wplus  += a_cellW;
          Wminus.plus(a_WcellAvgFab, cellBoxWpFaceAvgLRdom, 0, 0, numWcomp);
          Wplus .plus(a_WcellAvgFab, cellBoxWpFaceAvgLRdom, 0, 0, numWcomp);
        }  // End if use PPM limiter

//--Solve the Riemann problem

      if (CRDparam::g_usePPMlimiter ||
          CRDparam::g_faceInterpolationOrder % 2 == 1)
        {

          // Modify WfaceAvgDirFab[i + e/2] based on
          // Wplus [i    ] = WleftExtrap [i + e/2] and
          // Wminus[i + e] = WrightExtrap[i + e/2].

          // The Riemann box is based on cellBoxWfaceAvgLRdom.  We use this
          // instead of cellBoxWpfaceAvgDom because it pulls one cell away from
          // boundaries.  This is the only difference from
          // cellBoxWpfaceAvgDom.surroundingNodes(dir)
          Box riemannBox(cellBoxWpFaceAvgLRdom);
          riemannBox.surroundingNodes(dir);
          riemannBox.grow(dir, -1);

          // A forward transform to get the velocities normal to the faces on
          // mapped grids
          // ---------------------------------------
          //  I/O           Var            Stencil  
          // =====  ====================  ==========
          // read   Wplus                          1
          // read   Wminus                         1
          // read   a_unitNormalFxb                1
          // write  Wplus                          1
          // write  Wminus                         1
          // ---------------------------------------
          // The values above consider Wplus and
          // Wminus shifted to be face-centered.
          // ---------------------------------------
          PatchCNSOp::preRiemann(Wplus,
                                 Wminus,
                                 a_unitNormalFxb,
                                 dir,
                                 riemannBox);

          // Solve Riemann problem
          // ---------------------------------------
          //  I/O           Var            Stencil  
          // =====  ====================  ==========
          // read   Wplus                          1
          // read   Wminus                         1
          // write  WfaceAvgDirFab                 1
          // ---------------------------------------
          // The values above consider Wplus and
          // Wminus shifted to be face-centered.
          // ---------------------------------------
          CRDparam::g_CRDPhysics->riemann(
            WfaceAvgDirFab,  // On dir-FACEs
            Wplus,           // On cells to left of idir-FACEs
            Wminus,          // On cells to right of idir-FACEs
            dir,
            riemannBox);

          // A reverse transform of the velocity components
          // ---------------------------------------
          //  I/O           Var            Stencil  
          // =====  ====================  ==========
          // read   WfaceAvgDirFab                 1
          // read   a_unitNormalFxb                1
          // write  WfaceAvgDirFab                 1
          // ---------------------------------------
          PatchCNSOp::postRiemann(WfaceAvgDirFab,
                                  a_unitNormalFxb,
                                  dir,
                                  riemannBox);
        }

//--Replace interpolated state on domain boundaries with limited state

      // Boundary face boxes for implementing BC.  These extend to 3 ghosts
      // tangential along the boundary unless there is a non-periodic domain BC
      // there
      Box boundaryFaceBoxLo;
      CRDparam::g_CNSIBC->getBoundaryFaces(
        boundaryFaceBoxLo, cellBoxWpFaceAvgDom, a_domain, dir, Side::Lo);
      Box boundaryFaceBoxHi;
      CRDparam::g_CNSIBC->getBoundaryFaces(
        boundaryFaceBoxHi, cellBoxWpFaceAvgDom, a_domain, dir, Side::Hi);

      if (!boundaryFaceBoxLo.isEmpty())
        {
          // Replace values of WfaceAvgDirFab on lowest face with Wminus
          Wminus.shiftHalf(dir, -1);  // [i] becomes [i-e/2]
          WfaceAvgDirFab.copy(Wminus, boundaryFaceBoxLo, 0, boundaryFaceBoxLo,
                              0, numWcomp);
          if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
            {
              CRDparam::g_CRDPhysics->density(WfaceAvgDirFab,
                                              boundaryFaceBoxLo);
            }
        }
      if (!boundaryFaceBoxHi.isEmpty())
        {
          // Replace values of WfaceAvgDirFab on highest face with Wplus
          Wplus.shiftHalf(dir, 1);  // [i] becomes [i+e/2]
          WfaceAvgDirFab.copy(Wplus, boundaryFaceBoxHi, 0, boundaryFaceBoxHi,
                              0, numWcomp);
          if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
            {
              CRDparam::g_CRDPhysics->density(WfaceAvgDirFab,
                                              boundaryFaceBoxHi);
            }
        }
    }  // Loop over directions
}

/*--------------------------------------------------------------------*/
//  Compute the full average and point primitive state on all the
//  faces
/** From the face-averaged native primitive state, as estimated by the
 *  interior scheme:
 *  <ul>
 *    <li> a) Compute point values of native state
 *    <li> b) Compute full primitive state at points and averages.
 *         E.g., compute temperature T, T from <W>, and <T>
 *    <li> c) Apply BC (modifies <W> on boundary faces)
 *    <li> d) Fix point values of W (including T) on boundary faces
 *  </ul>
 *  The above operations are applied independently to (1) interior
 *  faces, (2) lower boundary faces, and (3) upper boundary faces in
 *  each direction.  The groups of operations (1)--(3) are fully
 *  independent across all directions and can be separated into tasks.
 *
 *  \param[in]  a_cellBoxWfacePnt
 *                      Box on which face-centered values must be known
 *                      In general the average native primitive state, 
 *                      <Wp> is known on the faces of a box grown by 1 
 *                      in tangential directions.
 *  \param[in]  a_cellBoxWpFaceAvg
 *                      Box on which limited, native, face-averaged values 
 *                      are known and where the non-native face-averaged
 *                      state must be found on boundaries.
 *  \param[in]  a_cellBoxWcellBdry
 *                      Box describing extent to which boundary ghost
 *                      cells must be filled in tangential directions.
 *                      It defines the number of boundary faces,
 *                      outside the disjoint box, on which the full
 *                      primitive state <W> is required.  Generally,
 *                      this is greater than 'a_box'.
 *  \param[in]  a_disjointBox
 *                      A disjoint box for which the calculations
 *                      are being performed
 *  \param[in]  a_domain
 *                      The problem or block domain for the disjoint
 *                      box
 *  \param[in]  a_WcellAvgFab
 *                      The full averaged primitive state, <W> in the
 *                      cells
 *  \param[out] a_WcellAvgFab
 *                      The full averaged primitive state, <W>, is set
 *                      in two layers of ghost cells orthogonal to the
 *                      boundary (as required)
 *  \param[in]  a_bndryCellFab
 *                      Previous-time data for cells immediately adjacent
 *                      to boundaries
 *  \param[in]  a_bndryNtJ
 *                      NtJ on near-bndry faces necessary for LES wall-model
 *  \param[in]  a_WfaceAvgFxb
 *                      The native primitive state <W> is set on all
 *                      faces by the interior scheme (including
 *                      boundary faces).  In general this state
 *                      must be known in 'a_boxWfacePnt' grown by 1 in
 *                      tangential directions.
 *  \param[out] a_WfaceAvgFxb
 *                      Extra primitive states are set on all faces
 *                      of 'a_boxWfacePnt'.  Also, faces on boundaries
 *                      are corrected according to the BC
 *  \param[out] a_WfacePntFxb
 *                      Point values of the full primtive state, W,
 *                      are set on all faces of 'a_boxWfacePnt'
 *  \param[in]  a_unitNormalFxb
 *                      Unit normal basis for transforming velocity
 *                      space to be normal to the faces
 *  \param[in]  a_dataIndx
 *                      Current data index
 *  \param[in]  a_time  Time (for applying BC)
 *  \param[in]  a_prevDt
 *                      Time-step size for past time-step
 *  \param[in]  a_level Index of the AMR level
 *
 *  \note
 *  <ul>
 *    <li> Native primitive state is the minimal state, (e.g.,
 *         density, velocity, pressure).  This is denoted by Wp.  The
 *         full primitive state adds extra variables.  The extra
 *         variables are denoted by Wx.  The most common example of an
 *         extra primitive variable for CNS is temperature.  The full
 *         primitive state is denoted simply by W.
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
Inertial4thOrderOp::computeFullFaceState(
  const Box&                          a_cellBoxWfacePnt,
  const Box&                          a_cellBoxWpFaceAvg,
  const Box&                          a_cellBoxWcellBdry,
  const Box&                          a_disjointBox,
  const BlockDomain&                  a_domain,
  FArrayBox&                          a_WcellAvgFab,
  stc::Vector<FArrayBox, 2*SpaceDim>& a_bndryCellFab,
  stc::Vector<FArrayBox, 2*SpaceDim>& a_bndryNtJ,
  FluxBox&                            a_WfaceAvgFxb,
  FluxBox&                            a_WfacePntFxb,
  const FluxBox&                      a_unitNormalFxb,
  const DataIndex&                    a_dataIndx,
  const Real                          a_time,
  const Real                          a_prevDt,
  const int                           a_level) const
{
  CRD::msg << CRD::fv4 << "Inertial4thOrderOp::computeFullFaceState"
           << CRD::end;

  //**FIXME Why Vector<Vector<>> instead of Vector<> name[SpaceDim]
  Vector<Vector<BCInfo> > domainBCInfoLo(SpaceDim);
  Vector<Vector<Box> > boundaryBoxesLo(SpaceDim);
  Vector<Vector<BCInfo> > domainBCInfoHi(SpaceDim);
  Vector<Vector<Box> > boundaryBoxesHi(SpaceDim);

  //**FIXME
  /*
    There is sort of a mess in the code here where boundary face values are
    needed for:
    a) the inertial solution: taken from a_cellBoxWpFaceAvg
    b) extrapolations to boundary ghost cells needed for viscous operator:
       taken from a_cellBoxWcellBdry.
    The mess is the BC are imposed twice and we need to fix that.
    boundarySlipVelocity adds to the confusion because it would be used in
    both locations.  So it needs to be sized to the maximum of either.
  */

  /*
    Also, since boundary states are defined for Wp, the temperature needs to
    be solved or set there.
   */

  // We need to store the wall slip-velocities possibly given by a turb model
  // Since this is also used for the inertial solver, we need a box that is the
  // minbox of a_cellBoxWpFaceAvg and a_cellBoxWcellBdry.
  //**FIXME Do we want this to be the full box or just localized to the
  //        boundaries of the box?
  Box cellBoxSlipVel(a_cellBoxWcellBdry);
  cellBoxSlipVel.minBox(a_cellBoxWpFaceAvg);
  FLUXBOXSTACKTEMP(boundarySlipVelocity, cellBoxSlipVel, SpaceDim);
  boundarySlipVelocity.setVal(0.);

  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      // Customize tangential cell per direction based on BC
      int numGhostWpFaceAvgTan =
        CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgTan);
      int numGhostWpFacePntTan =
        CRDparam::queryGhosts(CRDparam::NumGhostWpFacePntTan);
      // If there is a boundary in this direction, and any other direction is
      // periodic, ensure number of tangential cells is sufficient.
      {
        Box bdryTestBox = a_disjointBox;
        bdryTestBox.grow(dir, 1);
        if (!a_domain.contains(bdryTestBox))  // Not periodic in dir
          {
            numGhostWpFaceAvgTan = std::max(
              numGhostWpFaceAvgTan,
              CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry));
            numGhostWpFacePntTan = std::max(
              numGhostWpFacePntTan,
              CRDparam::queryGhosts(CRDparam::NumGhostWfacePntTngBdry));
          }
      }

      FArrayBox& WfaceAvgDirFab = a_WfaceAvgFxb[dir];
      FArrayBox& WfacePntDirFab = a_WfacePntFxb[dir];
      // Alias for only native primitive components
      FArrayBox WpFaceAvgDirFab(CRDparam::g_CRDPhysics->nativePrimInterval(),
                                WfaceAvgDirFab);
      FArrayBox WpFacePntDirFab(CRDparam::g_CRDPhysics->nativePrimInterval(),
                                WfacePntDirFab);

      // The cell and face boxes where the interior scheme is used to compute
      // WpPnt (== Wpnt == Wavg)
      Box faceBoxWpFacePntDom;
      {
        IntVect growVect(numGhostWpFacePntTan*IntVect::Unit);
        growVect[dir] = CRDparam::queryGhosts(CRDparam::NumGhostWpFacePntNrm);
        Box cellBoxWpFacePntDom = grow(a_disjointBox, growVect);
        cellBoxWpFacePntDom &= a_domain;
        faceBoxWpFacePntDom = cellBoxWpFacePntDom;
      }
      // // The growing and shrinking crops faces on the boundary in direction
      // // 'dir'
      // faceBoxWpFacePntDom.grow(dir, 1);
      // faceBoxWpFacePntDom &= a_domain;
      faceBoxWpFacePntDom.surroundingNodes(dir);  // Switch to faces
      // faceBoxWpFacePntDom.grow(dir, -1);
      CH_assert(WfacePntDirFab.box().contains(faceBoxWpFacePntDom));

      // The cell and face boxes where we need to know WpAvg in order to compute
      // WpPnt (== Wpnt == Wavg)
      Box faceBoxWpFaceAvgDom;
      {
        IntVect growVect(numGhostWpFaceAvgTan*IntVect::Unit);
        growVect[dir] = CRDparam::queryGhosts(CRDparam::NumGhostWpFaceAvgNrm);
        Box cellBoxWpFaceAvgDom = grow(a_disjointBox, growVect);
        cellBoxWpFaceAvgDom &= a_domain;
        faceBoxWpFaceAvgDom = cellBoxWpFaceAvgDom;
      }
      // // The growing and shrinking crops faces on the boundary in direction
      // // 'dir'
      // faceBoxWpFaceAvgDom.grow(dir, 1);
      // faceBoxWpFaceAvgDom &= a_domain;
      faceBoxWpFaceAvgDom.surroundingNodes(dir);  // Switch to faces
      // faceBoxWpFaceAvgDom.grow(dir, -1);
      CH_assert(WfaceAvgDirFab.box().contains(faceBoxWpFaceAvgDom));

//-- (1a) Compute point values of native W on face centers using Laplacian.

      const int faceDeconvOrd = (CRDparam::g_faceDeconvolveFlatten < 2) ? 4 : 2;
      bool faceDeconvLimiting =
        (CRDparam::g_faceDeconvolveLimit == 1) ? true : false;
      WpFacePntDirFab.copy(WpFaceAvgDirFab, faceBoxWpFacePntDom);
      CRDutil::deconvolveFace(WpFacePntDirFab,
                              WpFaceAvgDirFab,
                              faceBoxWpFacePntDom,
                              a_domain,
                              Interval(0, WpFacePntDirFab.nComp() - 1),
                              dir,
                              faceDeconvOrd,
                              faceDeconvLimiting,
                              true);
      Box faceBoxTestDom;
      IntVect growTestVect(numGhostWpFacePntTan*IntVect::Unit);
      growTestVect[dir] = CRDparam::queryGhosts(CRDparam::NumGhostWpFacePntNrm);
      Box cellBoxTestDom = grow(a_disjointBox, growTestVect) & a_domain;
      faceBoxTestDom = surroundingNodes(cellBoxTestDom, dir);
      faceBoxTestDom.grow(dir, 1) & a_domain;
      faceBoxTestDom.grow(dir, -1);
      // const int rad = (faceDeconvLimiting) ? 2 : 1; // radius for face checking
      // is density > 0 ?
      // CRDutil::checkFacePositivity(WpFacePntDirFab, WpFaceAvgDirFab,
      //                              faceBoxTestDom, a_domain,
      //                              CRDparam::g_CRDPhysics->densityIndex(),
      //                              dir, rad, "pnt-density", "avg-density");
      // // is pressure > 0 ?
      // CRDutil::checkFacePositivity(WpFacePntDirFab, WpFaceAvgDirFab,
      //                              faceBoxTestDom, a_domain,
      //                              CRDparam::g_CRDPhysics->pressureIndex(),
      //                              dir, rad, "pnt-pressure", "avg-pressure");
//-- (1b) Get full primitive state at points and averages

      const int numWcomp = CRDparam::g_CRDPhysics->numPrimitive();
      const Interval extraPrimIntv =
        CRDparam::g_CRDPhysics->extraPrimInterval();
      const int numExtraPrim = extraPrimIntv.size();

      // Compute 2nd-order extra primitive state into WfaceAvgDirFab
      CRDparam::g_CRDPhysics->extraPrimitiveState(WfaceAvgDirFab,
                                                  faceBoxWpFaceAvgDom);

      // Compute face-centered extra primitive state
      //**FIXME: one of the extraPrimitiveState operations is unnecessary
      //         if face convolutions/deconvolutions are turned off
      //         Additionally, one of these may need to be expanded if
      //         face convolution/deconvolution limiting is on
      CRDparam::g_CRDPhysics->extraPrimitiveState(WfacePntDirFab,
                                                  faceBoxWpFacePntDom);
      // Now WfaceAvgDirFab contains face-averaged primitive values
      // and 2nd-order non-native primitive values. Compare the
      // face-centered values with these values in the NLVF
      if (CRDparam::g_faceDeconvolveFlatten == 1)
        {
          CRDparam::g_DCF->directDCFlatten(faceBoxWpFacePntDom,
                                           WfacePntDirFab,
                                           WfaceAvgDirFab,
                                           numWcomp);
        }
      // Make FAB containing the 2nd-order extra primitive variables
      FABSTACKTEMP(WxFromWpFaceAvgFab, faceBoxWpFaceAvgDom, numExtraPrim);
      // Copy 2nd-order extra primitive variables to FAB
      WxFromWpFaceAvgFab.copy(WfaceAvgDirFab,
                              faceBoxWpFaceAvgDom,
                              extraPrimIntv.begin(),
                              faceBoxWpFaceAvgDom,
                              0,
                              numExtraPrim);
      // Convolve to get averages of extra primitives, e.g. temperature <T>,
      // using aliases to only work on extra primitives
      WfaceAvgDirFab.copy(WfacePntDirFab,
                          faceBoxWpFacePntDom,
                          extraPrimIntv.begin(),
                          faceBoxWpFacePntDom,
                          extraPrimIntv.begin(),
                          numExtraPrim);
      // Use alias to access only the extra primitive variables
      FArrayBox WxFaceAvgDirFab(extraPrimIntv, WfaceAvgDirFab);
      if (CRDparam::g_faceConvolveFlatten < 2)
        {
          MOLUtilFunc::deconvolveCenterFace(WxFaceAvgDirFab,
                                            WxFromWpFaceAvgFab,
                                            faceBoxWpFacePntDom,
                                            a_domain,
                                            dir,
                                            1);
        }
      if (CRDparam::g_faceConvolveFlatten == 1)
        {
          CRDparam::g_DCF->directDCFlatten(faceBoxWpFacePntDom,
                                           WxFaceAvgDirFab,
                                           WxFromWpFaceAvgFab,
                                           numExtraPrim,
                                           false); // cn and vel checks
          // no longer apply
        }

      //-- (1c): Apply boundary conditions to faces

      Box boundaryFaceBoxLo;
      CRDparam::g_CNSIBC->getBoundaryFaces(
        boundaryFaceBoxLo, a_cellBoxWpFaceAvg, a_domain, dir, Side::Lo);
      Box boundaryFaceGhostBoxLo;
      CRDparam::g_CNSIBC->getBoundaryFaces(
        boundaryFaceGhostBoxLo, a_cellBoxWcellBdry, a_domain, dir, Side::Lo);
      if (!boundaryFaceBoxLo.isEmpty())
        {
          BCInfo domainBCInfoLo;
          Box boundaryBoxLo;
          // Create composite box from boundaryFaceBox and boundaryFaceGhostBox
          Box compositeBox = minBox(boundaryFaceBoxLo,boundaryFaceGhostBoxLo);
          // Slip-velocity on walls from turbulence wall-model
          FABSTACKTEMP(bndrySlipVelDirFab, compositeBox, SpaceDim);
          CRDparam::g_CNSIBC->setBoxWallVelocity(
            bndrySlipVelDirFab,
            a_unitNormalFxb[dir],
            compositeBox,
            a_disjointBox,
            m_levelGridMetrics,
            a_domain,
            a_time,
            a_level,
            dir,
            Side::Lo);
          // Near-bndry-face NtJ for LES wall-model
          const FArrayBox& bndryNtJFab =
            a_bndryNtJ[MultiBlockRegions::FaceTag::indexFace(dir,Side::Lo)];
          // Slip-velocity calculation for LES wall-model
          CRDparam::g_CRDPhysics->turbModelSlipVel(
            bndrySlipVelDirFab,      // slip-velocity fab on the boundaries
            a_WcellAvgFab,           // cell-avg W for setting face values
            a_WfaceAvgFxb,           // need this for tangential derivatives
            m_levelGridMetrics,
            a_domain,
            a_dataIndx,
            a_unitNormalFxb[dir],
            bndryNtJFab,
            boundaryFaceBoxLo,
            boundaryFaceGhostBoxLo,
            a_disjointBox,
            dir,
            Side::Lo);

          // Temporary FAB for wall primitive state
          FABSTACKTEMP(WfaceAvgBndryDirFab,
                       compositeBox,
                       CRDparam::g_CRDPhysics->numPrimitive());

          // One bndry tangent ghost-face is required for every exterior tangent
          // ghost cell. For viscous, this is 4, but only 3 tangent ghost faces
          // have been solved. So, work out extra one from interior cells.

          FArrayBox& bndryCellFab =
            a_bndryCellFab[MultiBlockRegions::FaceTag::indexFace(dir,Side::Lo)];

          // Set face-values using boundary condition information
          setBndryConditions(
            boundaryFaceBoxLo,       // box to set using WfaceAvgDirFab
            boundaryFaceGhostBoxLo,  // box to set using a_WcellAvgFab
            a_disjointBox,           // disjointBox
            a_domain,                // domain
            WfaceAvgBndryDirFab,     // temporary fab to store wall state in
            WfaceAvgDirFab,          // state to overwrite with WfaceAvgBndryFab
            WfacePntDirFab,          // pnt state from avg state
            a_WcellAvgFab,           // cell-avg W for setting face values
            bndrySlipVelDirFab,      // slip-velocity fab on the boundaries
            bndryCellFab,            // bndry-cell-avg from prev time-points
            domainBCInfoLo,          // domain bcs
            boundaryBoxLo,           // boundary boxes
            dir,                     // direction
            Side::Lo,                // side
            a_unitNormalFxb[dir],    // unit normals
            a_dataIndx,              // data index
            a_time,                  // time
            a_prevDt,                // previous time-step size
            a_level);                // level

          // Extrapolate ghost cells and apply other boundary conditions
          // that are not imposed, like Neumann conditions
          //**NOTE: IMPORTANT: only set ghost cells here. Boundary face values
          //        are set by setBoundaryConditions and never changed. If at
          //        some future point, ghost faces exterior to bndry in normal
          //        dir must be set, set these using another function named
          //        extrapolateDomainGhostFaces. Don't mix with
          //        extrapolateDomainGhostCells.

          //**NOTE: placing extrapolateDomainGhostCells in this loop because it
          //        should have no requirement for setBoundaryConditions to be
          //        finished looping over every direction. Plus, it reduces
          //        memory overhead. Also, it just keeps things less complex.
          extrapolateDomainGhostCells(boundaryFaceBoxLo,
                                      boundaryFaceGhostBoxLo,
                                      a_disjointBox,
                                      a_domain,
                                      a_WcellAvgFab,
                                      WfaceAvgBndryDirFab,
                                      WfaceAvgDirFab,
                                      domainBCInfoLo,
                                      boundaryBoxLo,
                                      dir,
                                      Side::Lo,
                                      a_unitNormalFxb[dir],
                                      a_dataIndx,
                                      a_time,
                                      a_level);

          // If necessary, place extrapolateDomainGhostFaces() here
        }

      Box boundaryFaceBoxHi;
      CRDparam::g_CNSIBC->getBoundaryFaces(
        boundaryFaceBoxHi, a_cellBoxWpFaceAvg, a_domain, dir, Side::Hi);
      Box boundaryFaceGhostBoxHi;
      CRDparam::g_CNSIBC->getBoundaryFaces(
        boundaryFaceGhostBoxHi, a_cellBoxWcellBdry, a_domain, dir, Side::Hi);
      if (!boundaryFaceBoxHi.isEmpty())
        {
          BCInfo domainBCInfoHi;
          Box boundaryBoxHi;
          // Create composite box from boundaryFaceBox and boundaryFaceGhostBox
          Box compositeBox = minBox(boundaryFaceBoxHi,boundaryFaceGhostBoxHi);
          // Slip-velocity on walls from turbulence wall-model
          FABSTACKTEMP(bndrySlipVelDirFab, compositeBox, SpaceDim);
          CRDparam::g_CNSIBC->setBoxWallVelocity(
            bndrySlipVelDirFab,
            a_unitNormalFxb[dir],
            compositeBox,
            a_disjointBox,
            m_levelGridMetrics,
            a_domain,
            a_time,
            a_level,
            dir,
            Side::Hi);
          // Near-bndry-face NtJ for LES wall-model
          const FArrayBox& bndryNtJFab =
            a_bndryNtJ[MultiBlockRegions::FaceTag::indexFace(dir,Side::Hi)];
          // Slip-velocity calculation for LES wall-model
          CRDparam::g_CRDPhysics->turbModelSlipVel(
            bndrySlipVelDirFab,      // slip-velocity fab on the boundaries
            a_WcellAvgFab,           // cell-avg W for setting face values
            a_WfaceAvgFxb,           // need this for tangential derivatives
            m_levelGridMetrics,
            a_domain,
            a_dataIndx,
            a_unitNormalFxb[dir],
            bndryNtJFab,
            boundaryFaceBoxHi,
            boundaryFaceGhostBoxHi,
            a_disjointBox,
            dir,
            Side::Hi);

          // Temporary FAB for wall primitive state
          FABSTACKTEMP(WfaceAvgBndryDirFab,
                       compositeBox,
                       CRDparam::g_CRDPhysics->numPrimitive());

          FArrayBox& bndryCellFab =
            a_bndryCellFab[MultiBlockRegions::FaceTag::indexFace(dir,Side::Hi)];

          // Set face-values using boundary condition information
          setBndryConditions(
            boundaryFaceBoxHi,       // box to set using WfaceAvgDirFab
            boundaryFaceGhostBoxHi,  // box to set using a_WcellAvgFab
            a_disjointBox,           // disjointBox
            a_domain,                // domain
            WfaceAvgBndryDirFab,     // temporary fab to store wall state in
            WfaceAvgDirFab,          // state to overwrite with WfaceAvgBndryFab
            WfacePntDirFab,          // pnt state from avg state
            a_WcellAvgFab,           // cell-avg W for setting face values
            bndrySlipVelDirFab,      // slip-velocity fab on the boundaries
            bndryCellFab,            // bndry-cell-avg from prev time-points
            domainBCInfoHi,          // domain bcs
            boundaryBoxHi,           // boundary boxes
            dir,                     // direction
            Side::Hi,                // side
            a_unitNormalFxb[dir],    // unit normals
            a_dataIndx,              // data index
            a_time,                  // time
            a_prevDt,                // previous time-step size
            a_level);

          // Extrapolate ghost cells and apply other boundary conditions
          // that are not imposed, like Neumann conditions
          extrapolateDomainGhostCells(boundaryFaceBoxHi,
                                      boundaryFaceGhostBoxHi,
                                      a_disjointBox,
                                      a_domain,
                                      a_WcellAvgFab,
                                      WfaceAvgBndryDirFab,
                                      WfaceAvgDirFab,
                                      domainBCInfoHi,
                                      boundaryBoxHi,
                                      dir,
                                      Side::Hi,
                                      a_unitNormalFxb[dir],
                                      a_dataIndx,
                                      a_time,
                                      a_level);
        }
    } // Loop over directions
}

/*--------------------------------------------------------------------*/
//  Set the boundary conditions on the domain boundary faces
/** NOTE: this function does not set ghost cells. That is taken care
 *   of in extrapolateDomainGhostCells
 *  \param[in]  a_faceBdryBox
 *                      Box on which the full state is required as
 *                      points and averages.
 *  \param[in]  a_cellBdryBox
 *                      Boundary faces to set adjacent ghost cells
 *  \param[in]  a_disjointBox
 *                      A disjoint box for which the calculations
 *                      are being performed
 *  \param[in]  a_domain
 *                      The problem or block domain for the disjoint
 *                      box
 *  \param[out] a_WfaceDirFab
 *                      Faces on boundaries are set based on BC
 *  \param[in]  a_WcellFab
 *                      Cell based values, averaged or 2nd order
 *  \param[in]  a_boundarySlipVelocity
 *                      Slip velocities on the boundary provided by
 *                      turbulence model wall-model
 *  \param[in]  a_bndryCellFab
 *                      Previous-time data for cells immediately adjacent
 *                      to boundaries
 *  \param[out] a_bcInfoLo
 *                      BC info for lower boundaries
 *  \param[out] a_bcInfoHi
 *                      BC info for upper boundaries
 *  \param[out] a_bcBoxesLo
 *                      BC boxes for lower boundaries
 *  \param[out] a_bcBoxesHi
 *                      BC boxes for upper boundaries
 *  \param[in]  a_dir   Direction normal to the boundary faces
 *  \param[in]  a_unitNormalFxb
 *                      Unit normal basis for transforming velocity
 *                      space to be normal to the faces
 *  \param[in]  a_time  Time (for applying BC)
 *  \param[in]  a_prevDt
 *                      Time-step size for past time-step
 *  \param[in]  a_level Index of the AMR level
 *//*-----------------------------------------------------------------*/

void
Inertial4thOrderOp::setBndryConditions(
  const Box&            a_boundaryFaceBox,
  const Box&            a_boundaryFaceGhostBox,
  const Box&            a_disjointBox,
  const BlockDomain&    a_domain,
  FArrayBox&            a_WfaceAvgBndryDirFab,
  FArrayBox&            a_WfaceAvgDirFab,
  FArrayBox&            a_WfacePntDirFab,
  FArrayBox&            a_WcellAvgFab,
  const FArrayBox&      a_bndryFaceAvgSlipVelDirFab,
  FArrayBox&            a_bndryCellFab,
  BCInfo&               a_bcInfo,
  Box&                  a_bcBox,
  const int             a_dir,
  const Side::LoHiSide& a_side,
  const FArrayBox&      a_unitNormalDirFab,
  const DataIndex&      a_dataIndx,
  const Real            a_time,
  const Real            a_prevDt,
  const int             a_level) const
{
  CH_TIME("Inertial4thOrderOp::SetBoundaryConditions");
  CRD::msg << CRD::fv4
           << "Inertial4thOrderOp::SetBoundaryConditions" << CRD::end;

  // Set the face-averaged primitive state on domain boundaries
  CRDparam::g_CNSIBC->setBCprimState(a_WfaceAvgBndryDirFab,
                                     a_WfaceAvgDirFab,
                                     a_WcellAvgFab,
                                     a_bndryFaceAvgSlipVelDirFab,
                                     a_bndryCellFab,
                                     a_boundaryFaceBox,
                                     a_boundaryFaceGhostBox,
                                     a_disjointBox,
                                     a_unitNormalDirFab,
                                     a_bcInfo,
                                     a_bcBox,
                                     a_dir,
                                     a_side,
                                     m_levelGridMetrics,
                                     a_time,
                                     a_prevDt,
                                     a_level);

  // Set the face-averaged turbulence variables on domain boundaries
  //**FIXME: SA model uses ghost cell velocity to compute wall grad-velocity
  //         However, ghost cells aren't set until extrapolateDomainGhostCells.
  //         Since ghost cell velocity isn't known just yet, we shouldn't have
  //         this here. Problem is, different turb models have different wall
  //         models. For some, only turb variables are changed, but in others,
  //         standard primitives are changed (e.g. velocity). We should probably
  //         move this to a later function called "setBndryTurbVariables" to
  //         differentiate btwn this and "turbModelSlipVel" (which probably
  //         should be renamed to setBndryTurbModelSlipVelocity).
  //**NOTE:  The other option is to assume (for now) that turb models always
  //         have no-slip walls and wall grad-velocity can be computed with
  //         interior data plus the wall-boundary condition (one-sided stencil).
  //         For now, this is what we're going to assume.
  //**NOTE:  Main laziness keeping us from moving this out right now is that we
  //         would need an extra deconvolution to get the point turbulent state.
  //         But, then again, this sort of laziness keeps us from implementing
  //         a new one-sided gradient operator as well . . .
  if (CRDparam::g_turbModelType)
    {
      CRDparam::g_CRDPhysics->applyWallModel(
        a_boundaryFaceBox,
        a_boundaryFaceGhostBox,
        a_disjointBox,
        a_WfaceAvgBndryDirFab,
        a_WfaceAvgDirFab,
        a_WcellAvgFab,
        a_unitNormalDirFab,
        a_bcInfo,
        a_bcBox,
        m_levelGridMetrics,
        a_dataIndx,
        a_dir,
        a_side,
        a_time,
        a_level);
    }

  // Set up some boxes just for deconvolutions
  //**FIXME: is this necessary? We already have a_boundaryFaceBox and
  //         a_boundaryFaceGhostBox, so why not use these?
  // Where boundaryFaceBox intersects with stored average state <W>
  Box bdryAvgBox(a_boundaryFaceBox);
  bdryAvgBox &= a_WfaceAvgDirFab.box();
  // Where boundaryFaceBoxLo intersects with stored state W at points
  Box bdryPntBox(a_boundaryFaceBox);
  bdryPntBox &= a_WfacePntDirFab.box();

  // Fill a_WfacePntDirFab with a_WfaceAvgDirFab in prep for deconvolution
  a_WfacePntDirFab.copy(a_WfaceAvgDirFab, bdryPntBox);

  // Set up the checks for deconvolutions of faces
  const int faceDeconvOrd = (CRDparam::g_faceDeconvolveFlatten < 2) ? 4 : 2;
  bool faceDeconvLimit = (CRDparam::g_faceDeconvolveLimit == 1) ? true : false;

  // Deconvolve face-avg W to get face-pnt W on domain boundaries
  // Only deconvolve face values if extra boundary limiting is off
  //**NOTE: previously, the implementation deconvolved native W if running
  //        only inertial physics and all W if running anything else. Here,
  //        we wrap this into the definition of the deconvolution interval
  if (!CRDparam::g_extraBoundLim)
    {
      Interval deconvInterval = CRDparam::g_CRDPhysics->nativePrimInterval();
      if ((CRDparam::g_physicsModels & CRDparam::PhysicsViscous) ||
          (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf))
        {
          deconvInterval = Interval(0,CRDparam::g_CRDPhysics->numPrimitive()-1);
        }
      CRDutil::deconvolveFace(a_WfacePntDirFab,
                              a_WfaceAvgDirFab,
                              bdryPntBox,
                              a_domain,
                              deconvInterval,
                              a_dir,
                              faceDeconvOrd,
                              faceDeconvLimit,
                              true);
    }

  //**NOTE: retaining this code as a reminder of directDCFlatten. We may want to
  //        see if this works. Otherwise, we should remove it from this section.
  //**FIXME: remember to uncomment this
  // if (CRDparam::g_faceDeconvolveFlatten == 1)
  //   {
  //     const int numWcomp = CRDparam::g_CRDPhysics->numPrimitive();
  //     CRDparam::g_DCF->directDCFlatten(bdryPntBox,
  //                                      a_WfacePntDirFab,
  //                                      a_WfaceAvgDirFab,
  //                                      numWcomp);
  //   }

  //**NOTE: retaining this code as a reminder of face positivity checks. We may
  //        want to see if this is helpful. Otherwise, we should remove it from
  //        this section.
  // const int rad = (faceDeconvLimiting) ? 2 : 1; // radius for face check
  // is density > 0 ?
  // CRDutil::checkFacePositivity(WpFacePntDirFab, WpFaceAvgDirFab,
  //                              bdryPntBox, a_domain,
  //                              CRDparam::g_CRDPhysics->densityIndex(),
  //                              a_dir, rad, "pnt-density", "avg-density");
  // // is pressure > 0 ?
  // CRDutil::checkFacePositivity(WpFacePntDirFab, WpFaceAvgDirFab,
  //                              bdryPntBox, a_domain,
  //                              CRDparam::g_CRDPhysics->pressureIndex(),
  //                              a_dir,rad,"pnt-pressure","avg-pressure");
}

/*--------------------------------------------------------------------*/
//  Apply the boundary conditions to domain ghost cells
/** \param[in]  a_boundaryFaceBox
 *                      Box on which the full state is required as
 *                      points and averages.
 *  \param[in]  a_boundaryFaceGhostBox
 *                      Boundary faces to set adjacent ghost cells
 *  \param[in]  a_disjointBox
 *                      A disjoint box for which the calculations
 *                      are being performed
 *  \param[in]  a_domain
 *                      The problem or block domain for the disjoint
 *                      box
 *  \param[in]  a_WcellAvgFab
 *                      The full averaged primitive state, <W> in the
 *                      cells
 *  \param[out] a_WcellAvgFab
 *                      The full averaged primitive state, <W>, is set
 *                      in two layers of ghost cells orthogonal to the
 *                      boundary (as required)
 *  \param[in]  a_WfaceAvgDirFab
 *                      The native primitive state <W> is set on all
 *                      faces by the interior scheme (including
 *                      boundary faces).  In general this state
 *                      must be known in 'a_boundaryFaceBox' grown by
 *                      1 in tangential directions.
 *  \param[out] a_WfaceAvgDirFab
 *                      Faces on boundaries are corrected according to
 *                      the BC
 *  \param[out] a_WfacePntDirFab
 *                      Point values of the full primtive state, W,
 *                      are set on all faces of 'a_boundaryFaceBox'
 *  \param[in]  a_boundarySlipVelocity
 *                      Slip velocities on the boundary provided by
 *                      turbulence model wall-model
 *  \param[in]  a_bcInfo
 *                      BC info
 *  \param[in]  a_bcBoxes
 *                      BC boxes
 *  \param[in]  a_dir   Direction normal to the boundary faces
 *  \param[in]  a_side  Side of the domain
 *  \param[in]  a_unitNormalFxb
 *                      Unit normal basis for transforming velocity
 *                      space to be normal to the faces
 *  \param[in]  a_time  Time (for applying BC)
 *  \param[in]  a_level Index of the AMR level
 *//*-----------------------------------------------------------------*/

void
Inertial4thOrderOp::extrapolateDomainGhostCells(
  const Box&            a_boundaryFaceBox,
  const Box&            a_boundaryFaceGhostBox,
  const Box&            a_disjointBox,
  const BlockDomain&    a_domain,
  FArrayBox&            a_WcellAvgFab,
  const FArrayBox&      a_WfaceAvgBndryDirFab,
  const FArrayBox&      a_WfaceAvgDirFab,
  const BCInfo&         a_bcInfo,
  const Box&            a_bcBox,
  const int             a_dir,
  const Side::LoHiSide& a_side,
  const FArrayBox&      a_unitNormalDirFab,
  const DataIndex&      a_dataIndx,
  const Real            a_time,
  const int             a_level) const
{
  CH_TIME("Inertial4thOrderOp::extrapolateDomainGhostCells");
  CRD::msg << CRD::fv4
           << "Inertial4thOrderOp::extrapolateDomainGhostCells" << CRD::end;

  //**NOTE: this function should only be called to fill ghost cells outside
  //        the domain. It doesn't need to be called for inertial-only cases
  //**NOTE: cleverly, this was previously setting faces in temp memory that were
  //        then used to fill ghost cells, but never actually copied to
  //        a_WfaceAvgDirFab except for characteristic BCs and imposed Neumann
  //        conditions (tricky to follow). Theoretically, this can lead to
  //        inconsistent treatment of domain ghost cells and boundary faces.
  //        This function now takes consistently set a_WfaceAvgBndryDirFab and
  //        uses it with interior cells to set domain ghosts
  //**NOTE: the previous implementation in applyBoundaryConditions was also not
  //        necessarily 4th-order as it used first-order face-values in combo
  //        with interior cell-values (the face and first-interior cell would
  //        interact in a way that could decrease the order-of-accuracy)
  if ((CRDparam::g_physicsModels & CRDparam::PhysicsViscous) ||
      (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf))
    {
      CRDparam::g_CNSIBC->extrapolateDomainGhostCells(
        a_WfaceAvgBndryDirFab,
        a_WfaceAvgDirFab,
        a_WcellAvgFab,
        a_boundaryFaceBox,
        a_boundaryFaceGhostBox,
        a_disjointBox,
        a_unitNormalDirFab,
        a_dir,
        a_side,
        a_bcInfo,
        a_bcBox,
        m_levelGridMetrics,
        a_time,
        a_level);
    }
}

/*--------------------------------------------------------------------*/
//  Compute fluxes on the faces from primitive state W on faces
/** This routine can be used to compute the flux and a point from
 *  point values of W or to compute a flux from the average primitive
 *  state \<W\>.  The latter is used to computed gradients for a
 *  product rule.
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[in]  a_domain
 *                      The problem or block domain for the disjoint
 *                      box
 *  \param[out] a_fluxFaceFxb
 *                      Flux on the faces
 *  \param[in]  a_WfaceFxb
 *                      Primitive state to use
 *  \param[in]  a_numTanGhost
 *                      Number of ghost cells in the tangent direction
 *                      to solve for the flux values
 *//*-----------------------------------------------------------------*/

// ---------------------------------------------------
//  I/O           Var            Stencil      Ghosts
// =====  ====================  ==========  ==========
// read   a_WfaceFxb                     1    see
// write  a_fluxFaceFxb                  1    caller
// ---------------------------------------------------

void
Inertial4thOrderOp::computeAllFluxes(const Box&         a_disjointBox,
                                     const BlockDomain& a_domain,
                                     FluxBox&           a_fluxFaceFxb,
                                     const FluxBox&     a_WfaceFxb,
                                     const int          a_numTanGhost) const
{
  CH_TIME("Inertial4thOrderOp::computeAllFluxes");

  CH_assert(a_fluxFaceFxb.box().contains(a_disjointBox));
  CH_assert(a_WfaceFxb.box().contains(a_disjointBox));

  const int numFcomp = CRDparam::g_CRDPhysics->numFluxes();

  // T - we need the full flux dyad on each face
  // F - we only need the flux in the direction orthogonal to the face
  const bool allFluxDir = (a_fluxFaceFxb.nComp() >= SpaceDim*numFcomp);

  for (int faceDir = 0; faceDir != SpaceDim; ++faceDir)
    {
      Box faceBox(a_disjointBox);
      IntVect growVect(a_numTanGhost*IntVect::Unit);
      // Only need flux values in tangent directions
      growVect[faceDir] = 0;
      faceBox.grow(growVect);
      faceBox &= a_domain;
      faceBox.surroundingNodes(faceDir);
                                                           // Old names
      const FArrayBox& WfaceDirFxb = a_WfaceFxb[faceDir];  //(WfaceCenDir)
      FArrayBox& fluxFaceDirFxb = a_fluxFaceFxb[faceDir];  //(FfaceCenDir)
      CH_assert(fluxFaceDirFxb.contains(faceBox));
      CH_assert(WfaceDirFxb.contains(faceBox));

      if (allFluxDir)
        {
          for (int fluxDir = 0; fluxDir < SpaceDim; ++fluxDir)
            {
              // Interval
              int fluxIntvLo = fluxDir*numFcomp;
              int fluxIntvHi = fluxIntvLo + numFcomp - 1;
              Interval fluxIntv(fluxIntvLo, fluxIntvHi);
              // Alias to pretend the fluxbox has only m_fluxes components
              FArrayBox fluxFaceDirFxb_fluxDir(fluxIntv, fluxFaceDirFxb);
              CRDparam::g_CRDPhysics->getFlux(fluxFaceDirFxb_fluxDir,
                                              WfaceDirFxb,
                                              fluxDir,
                                              faceBox);
            }
        }
      else
        {
          CH_assert(a_fluxFaceFxb.nComp() >= numFcomp);
          // Set a_fluxFacePntFxb = flux(a_WfacePntFxb) on all faceDir faces of
          // a_box
          CRDparam::g_CRDPhysics->getFlux(fluxFaceDirFxb,
                                          WfaceDirFxb,
                                          faceDir,
                                          faceBox);
        }
    }
}

/*--------------------------------------------------------------------*/
//  Solves for the fifth-order W- and W+ values
/** \param[out] a_Wminus
 *                      FAB containing the W- values, shifted to be
 *                      cell-centered
 *  \param[out] a_Wplus FAB containing the W+ values, shifted to be
 *                      cell-centered
 *  \param[in]  a_WcellAvgFab
 *                      FAB containing the cell-averaged primitive values
 *  \param[in]  a_numSlopes
 *                      Number of slopes to solve
 *  \param[in]  a_dir   Direction of faces
 *  \param[in]  a_box   Box to solve over
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *//*-----------------------------------------------------------------*/

void
Inertial4thOrderOp::solveFifthOrderVals(FArrayBox&         a_Wminus,
                                        FArrayBox&         a_Wplus,
                                        const FArrayBox&   a_WcellAvgFab,
                                        const int          a_numSlopes,
                                        const int          a_dir,
                                        const Box&         a_box,
                                        const BlockDomain& a_domain) const
{
  CH_TIME("Inertial4thOrderOp::solveFifthOrderVals");
  CRD::msg << CRD::fv4 << "Inertial4thOrderOp::solveFifthOrderVals" << CRD::end;
  const Real fact = 1./60.;
  const Real fact2 = 1./30.;
  Box loBox, nextLoBox, hiBox, nextHiBox, centerBox, innerCenterBox, entireBox;
  int hasLo, hasHi;
  Box inBox(a_box);
  inBox.grow(a_dir,1);
  loHiCenter5(loBox, nextLoBox, hasLo, hiBox, nextHiBox, hasHi, centerBox,
              innerCenterBox, entireBox, inBox, a_domain, a_dir);
  CH_assert(a_Wminus.nComp() == a_numSlopes);
  CH_assert(a_Wplus.nComp() == a_numSlopes);
  CH_assert(a_Wminus.box().contains(entireBox));
  CH_assert(a_Wplus.box().contains(entireBox));
  MD_ARRAY_RESTRICT(arrWminus, a_Wminus);
  MD_ARRAY_RESTRICT(arrWplus, a_Wplus);
  MD_ARRAY_RESTRICT(arrWCell, a_WcellAvgFab);
  const int MD_ID(o, a_dir);
  const Real alpha = CRDparam::g_fifthOrderBlendingCoef;
  for (int comp = 0; comp != a_numSlopes; ++comp)
    {
      MD_BOXLOOP(innerCenterBox, i)
        {
          Real W_im2 = arrWCell[MD_OFFSETIX(i,-,2*o,comp)];
          Real W_im1 = arrWCell[MD_OFFSETIX(i,-,o,comp)];
          Real W_i = arrWCell[MD_IX(i,comp)];
          Real W_ip1 = arrWCell[MD_OFFSETIX(i,+,o,comp)];
          Real W_ip2 = arrWCell[MD_OFFSETIX(i,+,2*o,comp)];
          // Compute a 1st-order 4th-derivative and use this for continuous
          // blending between 5th-order and 4th-order interpolation
          // Coefficient for 4th-deriv with centered stencils is 1
          Real d4 = fact2*(W_im2 - 4.*W_im1 + 6.*W_i - 4.*W_ip1 + W_ip2);
          arrWplus[MD_IX(i,comp)] = (2.*W_im2 - 13.*W_im1 + 47.*W_i + 27.*W_ip1
                                     - 3.*W_ip2)*fact - alpha*d4;
          arrWminus[MD_IX(i,comp)] = (2.*W_ip2 - 13.*W_ip1 + 47.*W_i + 27.*W_im1
                                      - 3.*W_im2)*fact - alpha*d4;
        }
      if (hasLo)
        {
          MD_BOXLOOP(loBox, i)
            {
              Real W_i = arrWCell[MD_IX(i,comp)];
              Real W_ip1 = arrWCell[MD_OFFSETIX(i,+,o,comp)];
              Real W_ip2 = arrWCell[MD_OFFSETIX(i,+,2*o,comp)];
              Real W_ip3 = arrWCell[MD_OFFSETIX(i,+,3*o,comp)];
              Real W_ip4 = arrWCell[MD_OFFSETIX(i,+,4*o,comp)];
              // Compute a 1st-order 4th-derivative and use this for continuous
              // blending between 5th-order and 4th-order interpolation
              // Coefficient for 4th-deriv with biased stencil is -3/2
              // Coefficient for 4th-deriv with one-sided stencil is 6
              Real d4 = fact2*(W_i - 4.*W_ip1 + 6.*W_ip2 - 4.*W_ip3 + W_ip4);
              arrWplus[MD_IX(i,comp)] =
                (-3.*W_ip4 + 17.*W_ip3 - 43.*W_ip2
                 + 77.*W_ip1 + 12.*W_i)*fact - (-3./2.)*alpha*d4;
              arrWminus[MD_IX(i,comp)] =
                (12.*W_ip4 - 63.*W_ip3 - 163.*W_ip1
                 + 137.*(W_ip2 + W_i))*fact - 6.*alpha*d4;
            }
          MD_BOXLOOP(nextLoBox, i)
            {
              Real W_im1 = arrWCell[MD_OFFSETIX(i,-,o,comp)];
              Real W_i = arrWCell[MD_IX(i,comp)];
              Real W_ip1 = arrWCell[MD_OFFSETIX(i,+,o,comp)];
              Real W_ip2 = arrWCell[MD_OFFSETIX(i,+,2*o,comp)];
              Real W_ip3 = arrWCell[MD_OFFSETIX(i,+,3*o,comp)];
              // Compute a 1st-order 4th-derivative and use this for continuous
              // blending between 5th-order and 4th-order interpolation
              // Coefficient for 4th-deriv with centered stencils is 1
              // Coefficient for 4th-deriv with biased stencil is -3/2
              Real d4 = fact2*(W_im1 - 4.*W_i + 6.*W_ip1 - 4.*W_ip2 + W_ip3);
              arrWplus[MD_IX(i,comp)] = (2.*W_ip3 - 13.*W_ip2 + 47.*W_ip1
                                         + 27.*W_i - 3.*W_im1)*fact - alpha*d4;
              arrWminus[MD_IX(i,comp)] =
                (-3.*W_ip3 + 17.*W_ip2 - 43.*W_ip1
                 + 77.*W_i + 12.*W_im1)*fact - (-3./2.)*alpha*d4;
            }
        }
      if (hasHi)
        {
          MD_BOXLOOP(hiBox, i)
            {
              Real W_i = arrWCell[MD_IX(i,comp)];
              Real W_im1 = arrWCell[MD_OFFSETIX(i,-,o,comp)];
              Real W_im2 = arrWCell[MD_OFFSETIX(i,-,2*o,comp)];
              Real W_im3 = arrWCell[MD_OFFSETIX(i,-,3*o,comp)];
              Real W_im4 = arrWCell[MD_OFFSETIX(i,-,4*o,comp)];
              // Compute a 1st-order 4th-derivative and use this for continuous
              // blending between 5th-order and 4th-order interpolation
              // Coefficient for 4th-deriv with one-sided stencil is 6
              // Coefficient for 4th-deriv with biased stencil is -3/2
              Real d4 = fact2*(W_i - 4.*W_im1 + 6.*W_im2 - 4.*W_im3 + W_im4);
              arrWplus[MD_IX(i,comp)] =
                (12.*W_im4 - 63.*W_im3 - 163.*W_im1
                 + 137.*(W_im2 + W_i))*fact - 6.*alpha*d4;
              arrWminus[MD_IX(i,comp)] =
                (-3.*W_im4 + 17.*W_im3 - 43.*W_im2
                 + 77.*W_im1 + 12.*W_i)*fact - (-3./2.)*alpha*d4;
            }
          MD_BOXLOOP(nextHiBox, i)
            {
              Real W_ip1 = arrWCell[MD_OFFSETIX(i,+,o,comp)];
              Real W_i = arrWCell[MD_IX(i,comp)];
              Real W_im1 = arrWCell[MD_OFFSETIX(i,-,o,comp)];
              Real W_im2 = arrWCell[MD_OFFSETIX(i,-,2*o,comp)];
              Real W_im3 = arrWCell[MD_OFFSETIX(i,-,3*o,comp)];
              // Compute a 1st-order 4th-derivative and use this for continuous
              // blending between 5th-order and 4th-order interpolation
              // Coefficient for 4th-deriv with biased stencil is -3/2
              // Coefficient for 4th-deriv with centered stencils is 1
              Real d4 = fact2*(W_ip1 - 4.*W_i + 6.*W_im1 - 4.*W_im2 + W_im2);
              arrWplus[MD_IX(i,comp)] =
                (-3.*W_im3 + 17.*W_im2 - 43.*W_im1
                 + 77.*W_i + 12.*W_ip1)*fact - (-3./2.)*alpha*d4;
              arrWminus[MD_IX(i,comp)] = (2.*W_im3 - 13.*W_im2 + 47.*W_im1
                                          + 27.*W_i - 3.*W_ip1)*fact - alpha*d4;
            }
        }
    } 
}

/*--------------------------------------------------------------------*/
//  Solves for the second-order W values
/** \param[out] a_WfaceAvgDirFab
 *                      FAB containing the W face-averaged values
 *  \param[out] a_WcellAvgFab 
 *                      FAB containing the W cell-averaged values
 *  \param[in]  a_numWcomp
 *                      Number of W face-averaged values to solve for
 *  \param[in]  a_dir   Direction of faces
 *  \param[in]  a_faceBox
 *                      Faces on which solution is needed
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *//*-----------------------------------------------------------------*/

void
Inertial4thOrderOp::solveSecondOrderVals(FArrayBox&         a_WfaceAvgDirFab,
                                         const FArrayBox&   a_WcellAvgFab,
                                         const int          a_numWcomp,
                                         const int          a_dir,
                                         const Box&         a_faceBox,
                                         const BlockDomain& a_domain) const
{
  CH_TIME("Inertial4thOrderOp::solveSecondOrderVals");
  CRD::msg << CRD::fv4 << "Inertial4thOrderOp::solveSecondOrderVals"
           << CRD::end;

  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;
  Box inBox = grow(a_faceBox, IntVect_basis(a_dir));
  inBox.enclosedCells();
  loHiCenterFace(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                 inBox, a_domain, a_dir);

  const int MD_ID(o, a_dir);
  const Real hiTol = 1e100; // high tolerance to check if data is valid
  const Real loTol = -1e100; // low tolerance to check if data is valid
  for (int comp = 0; comp != a_numWcomp; ++comp)
    {
      MD_BOXLOOP(centerBox, i)
        {
          Real W_im1 = a_WcellAvgFab[MD_OFFSETIX(i,-,o,comp)];
          Real W_i  = a_WcellAvgFab[MD_IX(i,comp)];
          a_WfaceAvgDirFab[MD_IX(i,comp)] = 0.5*(W_im1 + W_i);

          CH_assert(W_im1 < hiTol && W_im1 > loTol);
          CH_assert(W_i < hiTol && W_i > loTol);
        }
      if (hasLo)
        {
          MD_BOXLOOP(loBox, i)
            {
              Real W_ip1 = a_WcellAvgFab[MD_OFFSETIX(i,+,o,comp)];
              Real W_i   = a_WcellAvgFab[MD_IX(i,comp)];
              a_WfaceAvgDirFab[MD_IX(i,comp)] = 0.5*(3.*W_i - W_ip1);

              CH_assert(W_ip1 < hiTol && W_ip1 > loTol);
              CH_assert(W_i < hiTol && W_i > loTol);
            }
        }
      if (hasHi)
        {
          MD_BOXLOOP(hiBox, i)
            {
              Real W_im1 = a_WcellAvgFab[MD_OFFSETIX(i,-,o,comp)];
              Real W_im2 = a_WcellAvgFab[MD_OFFSETIX(i,-,2*o,comp)];
              a_WfaceAvgDirFab[MD_IX(i,comp)] = 0.5*(3.*W_im1 - W_im2);

              CH_assert(W_im2 < hiTol && W_im2 > loTol);
              CH_assert(W_im1 < hiTol && W_im1 > loTol);
            }
        }
    } 
}
