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
 * \file PatchCNSOp.cpp
 *
 * \brief Member functions for PatchCNSOp
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "LevelFluxRegister.H"
#include "DivergenceF_F.H"
#include "MOLUtilFunc.H"
#include "LevelGridMetrics.H"

//----- Internal -----//

#include "PatchCNSOp.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "CRDutil.H"
#include "PatchMappedFunc.H"
#include "CRDPhysics.H"
#include "DataTemp.H"
#include "CNSIBC.H"
#include "DCFlattening.H"


/*******************************************************************************
 *
 * Class PatchCNSOp: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_levelGridMetrics
 *                      Grid metrics for the level
 *--------------------------------------------------------------------*/

PatchCNSOp::PatchCNSOp(
  LevelGridMetrics& a_levelGridMetrics)
  :
  m_levelGridMetrics(a_levelGridMetrics),
  m_inertial4thOrderOp(a_levelGridMetrics),
  m_viscous4thOrderOp(a_levelGridMetrics)
{ }


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Compute the primitive cell point and cell average states from the
//  conservative cell average state
/** \param[out]  a_UcellPntFab
 *                      The conservative cell point state
 *  \param[out]  a_WcellPntFab
 *                      The primitive cell point state
 *  \param[out]  a_WcellAvgFab
 *                      The primitive cell average state
 *  \param[out]  a_dataIndx
 *                      Current data index
 *  \param[in]   a_box  Disjoint box region to compute the state
 *  \param[in]   a_domain
 *                      The problem or block domain for the box
 *  \param[in]   a_UcellAvgFab
 *                      The conservative cell average state
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::computeWcell(FArrayBox&         a_UcellPntFab,
                         FArrayBox&         a_WcellPntFab,
                         FArrayBox&         a_WcellAvgFab,
                         FArrayBox&         a_WfromUavgFab,
                         const DataIndex&   a_dataIndx,
                         const Box&         a_box,
                         const BlockDomain& a_domain,
                         const FArrayBox&   a_UcellAvgFab,
                         FArrayBox&         a_WOld) const
{
  CH_TIME("PatchCNSOp::computeWcell");
  CRD::msg << CRD::fv4 << "PatchCNSOp::computeWcell" << CRD::end;

/* When do we want a consistent primitive state?  Consistent can mean 3 things:
   a) p = \rho R T
   b) sum_j cn_j = 1
   c) 0 <= cn_j <= 1
   The answer may be never, just get the flux.  But running the Riemann solver
   without a consistent left or right state is confusing.
   Rules:
   1) Enforce a) at WcellPnt and WfacePnt.  There is no reason to expect it to
      hold for averages.  This means that the Riemann solution should be
      performed on points.  Historical practice recommends operating on
      averages though.  It would be nice to adjust R
   2) b) and c) can be more easily achieved if \rho = \sum_j (\rho*cn)_j
   3) We generally prefer obtaining \rho from p and T anyways
*/

//--Boxes (reminder: these are all cell-centered, even if they define data on
//--faces).  Face boxes are computed when looping over directions.

  // Cell Uavg and cell WfromUavg
  Box boxUcellAvgDom =
    grow(a_box, CRDparam::queryGhosts(CRDparam::NumGhostUcellAvg));
  CH_assert(a_UcellAvgFab.box().contains(boxUcellAvgDom));
  boxUcellAvgDom &= a_domain;
  const Box& boxWfromUavgDom = boxUcellAvgDom;  // Rename for clarity

  // Cell Upnt - deconvolve Uavg (requires 1 (normal) or 2 (limited) cells)
  // Note: NumGhostUcellPnt == NumGhostWcellAvg
  CH_assert(CRDparam::queryGhosts(CRDparam::NumGhostUcellPnt) + 1 <=
            CRDparam::queryGhosts(CRDparam::NumGhostUcellAvg));
  Box boxUcellPntDom =
    grow(a_box, CRDparam::queryGhosts(CRDparam::NumGhostUcellPnt));
  boxUcellPntDom &= a_domain;
  CH_assert(a_UcellPntFab.box().contains(boxUcellPntDom));

  // Cell Wpnt - compute from Upnt
  // Note: NumGhostWcellPnt == NumGhostWcellAvg
  const Box& boxWcellPntDom = boxUcellPntDom;  // Rename for clarity
  CH_assert(a_WcellPntFab.box().contains(boxWcellPntDom));

  // Cell Wavg - deconvolve Wpnt using WfromUavg for derivatives
  const Box& boxWcellAvgDom = boxUcellPntDom;  // Rename for clarity
  CH_assert(a_WcellAvgFab.box().contains(boxWcellAvgDom));

//--Numbers of components

  const int numUcomp  = CRDparam::g_CRDPhysics->numConservative();
  const int numWcomp  = CRDparam::g_CRDPhysics->numPrimitive();
  CH_assert(a_UcellAvgFab.nComp() == numUcomp);
  CH_assert(a_WcellAvgFab.nComp() == numWcomp);
  CH_assert(a_WcellPntFab.nComp() == numWcomp);
  CH_assert(a_WfromUavgFab.nComp() == numWcomp);

//-- Start of work -------------------------------------------------------------

  // NOTE: Stencil estimates do not include 1-sided operations

  // Compute a primitive state based on average of U (reduces stencil)
  // ---------------------------------------
  //  I/O           Var            Stencil  
  // =====  ====================  ==========
  // read   a_UcellAvgFab                  1
  // write  a_WfromUavgFab                 1
  // ---------------------------------------
//**FIXME The whole native versus full primitive state thing has ended up as
//**      being quite convoluted (in a non-mathematical sense).  It's probably
//**      worthwhile to dispense with the concept and just use the full
//**      primitive state everywhere.  That might be combined with a
//**      'make consistent perfect gas' function that ensures p = rho R T;
  CRDparam::g_CRDPhysics->consToPrim(a_WfromUavgFab,
                                     a_UcellAvgFab,
                                     boxWfromUavgDom,
                                     a_WOld);

  // Correct pressure using SGS KE estimate -- unnecessary for subsonic flows
  // ---------------------------------------
  //  I/O           Var            Stencil
  // =====  ====================  ==========
  // read   a_WfromUAvgFab                 0
  // write  a_WfromUavgFab                 0
  // ---------------------------------------
  CRDparam::g_CRDPhysics->consToPrimCorrectionLES(a_WfromUavgFab,
                                                  boxWfromUavgDom);
  
  // Add the extra primitive state (T or p) to WfromUavg
  // ---------------------------------------
  //  I/O           Var            Stencil  
  // =====  ====================  ==========
  // read   a_WfromUavgFab                 1
  // write  a_WfromUavgFab                 1
  // ---------------------------------------
  CRDparam::g_CRDPhysics->extraPrimitiveState(a_WfromUavgFab, boxWfromUavgDom);

  // Compute point values of U at cell centers using Laplacian
  // ---------------------------------------
  //  I/O           Var            Stencil  
  // =====  ====================  ==========
  // read   a_UcellAvgFab          3 per dir
  // write  UcellPntFab                    1
  // ---------------------------------------
  // Deconvolve to get Upnt from Uavg (with deconvolution limiting)
  a_UcellPntFab.copy(a_UcellAvgFab);
  const int cellDeconvOrd = (CRDparam::g_cellDeconvolveFlatten < 2) ? 4 : 2;
  bool cellDeconvLimiting =
    (CRDparam::g_cellDeconvolveLimit == 1) ? true : false;
  CRDutil::deconvolve(
    a_UcellPntFab,
    a_UcellAvgFab,
    boxUcellPntDom,
    a_domain,
    Interval(0, CRDparam::g_CRDPhysics->numConservative() - 1),
    cellDeconvOrd,
    1,
    cellDeconvLimiting,
    true);

  // Compute point values of W at cell centers from point values of U
  // ---------------------------------------
  //  I/O           Var            Stencil  
  // =====  ====================  ==========
  // read   UcellPntFab                    1
  // write  WcellPntFab                    1
  // ---------------------------------------
  CRDparam::g_CRDPhysics->consToPrim(a_WcellPntFab,
                                     a_UcellPntFab,
                                     boxWcellPntDom,
                                     a_WOld);

  // Correct pressure using SGS KE estimate -- unnecessary for subsonic flows
  // ---------------------------------------
  //  I/O           Var            Stencil
  // =====  ====================  ==========
  // read   a_WfromUAvgFab                 0
  // write  a_WfromUavgFab                 0
  // ---------------------------------------
  CRDparam::g_CRDPhysics->consToPrimCorrectionLES(a_WcellPntFab,
                                                  boxWcellPntDom);

  // Add the extra primitive state (T or p) to W
  // ---------------------------------------
  //  I/O           Var            Stencil  
  // =====  ====================  ==========
  // read   WcellPntFab                    1
  // write  WcellPntFab                    1
  // ---------------------------------------
  CRDparam::g_CRDPhysics->extraPrimitiveState(a_WcellPntFab, boxWcellPntDom);

  // const int rad = (cellDeconvLimiting) ? 2 : 1; // radius for cell checking
  // is density > 0 ? (from <U> --> U)
  // CRDutil::checkCellPositivity(UcellPntFab, a_UcellAvgFab, boxUcellPntDom,
  //                              a_domain, CRDparam::g_CRDPhysics->densityIndex(),
  //                              rad, "pnt-density", "avg-density");
  // // is pressure > 0 ? (from U --> W)
  // CRDutil::checkCellPositivity(a_WcellPntFab, a_WcellPntFab,
  //                              boxWcellPntDom, a_domain,
  //                              CRDparam::g_CRDPhysics->pressureIndex(),
  //                              0, "pnt-pressure", "pnt-pressure");

  // // Compute most of primitive state if thermally perfect, all otherwise
  // CRDparam::g_CRDPhysics->intermediateConsToPrim(a_WcellPntFab, 
  //                                                UcellPntFab,
  //                                                boxWcellPntDom);

  // if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
  //   {
  //     // Compute temperature
  //     CRDutil::limitTpnt(a_UcellAvgFab,
  //                        a_WcellPntFab,
  //                        boxWcellPntDom,
  //                        m_domain,
  //                        true);

  //     // Compute pressure
  //     CRDparam::g_CRDPhysics->pressure(a_WcellPntFab, boxWcellPntDom);
  //   }
  // else
  //   {
  //     CRDparam::g_CRDPhysics->extraPrimitiveState(a_WcellPntFab,
  //                                                 boxWcellPntDom);
  //   }

  // Compute cell average of W using Laplacian based on a_WfromUavgFab
  // ---------------------------------------
  //  I/O           Var            Stencil  
  // =====  ====================  ==========
  // read   WcellPntFab                    1
  // read   a_WfromUavgFab         3 per dir
  // write  a_WcellAvgFab                  1
  // ---------------------------------------
  a_WcellAvgFab.copy(a_WcellPntFab);
  if (CRDparam::g_cellConvolveFlatten < 2)
    {
      //**FIXME costly to copy and then perform operation. Do copy in Op.
      MOLUtilFunc::deconvolve(a_WcellAvgFab,  // Convolution despite the name
                              a_WfromUavgFab,
                              boxWcellAvgDom,
                              a_domain,
                              1);             // because of this 1
      // Mostly leave this state alone.  However, the species mass fractions
      // should sum to 1.
      CRDparam::g_CRDPhysics->normalizePrimSpecies(
        CRDPhysics::NormalizeTypeRedistribute,
        true,   // Bound
        true,   // Sort
        boxWcellAvgDom,
        a_WcellAvgFab);

      if (CRDparam::g_cellConvolveFlatten == 1)
        {
          CRDparam::g_DCF->directDCFlatten(boxWcellAvgDom,
                                           a_WcellAvgFab,
                                           a_WfromUavgFab,
                                           numWcomp);
        }
    }
  // Apply DC flattening and blend the high-order and low-order
  // cell-centered values
  if (CRDparam::g_cellDeconvolveFlatten == 1)
    {
      CRDparam::g_DCF->directDCFlatten(boxWcellAvgDom,
                                       a_WcellPntFab,
                                       a_WfromUavgFab,
                                       numWcomp);
    }
}

/*--------------------------------------------------------------------*/
//  Compute the primitive cell-centered state from the conservative
//  This is a more general function than computeWcell where
//  we assume that we're only computing the primitive state
//  (a) on the input box,
//  (b) at the order-of-accuracy specified by the input,
//  (c) as cell-centered or cell-averaged values based on an input
//  As opposed to computeWcell where we make a lot more decisions for
//  the user
/** \param[out]  a_WcellPntFab
 *                      The primitive cell point state
 *  \param[in]   a_box  Disjoint box region to compute the state
 *  \param[in]   a_domain
 *                      The problem or block domain for the box
 *  \param[in]   a_UcellAvgFab
 *                      The conservative cell average state
 *  \param[in]   a_cellAvgInput
 *                      If true, must deconvolve UcellAvg first
 *                      If false, nothing must be done
 *  \param[in]   a_fourthOrder
 *                      If true, see a_cellAvgInput
 *                      If false, disregard a_cellAvgInput and just
 *                      go ahead with the consToPrim conversion
 *  \param[in]   a_WOld An old primitive state used to help stabilize
 *                      constToPrim operations
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::computeWstate(FArrayBox&         a_WcellPntFab,
                          const FArrayBox&   a_UcellAvgFab,
                          FArrayBox&         a_WOld,
                          const Box&         a_box,
                          const BlockDomain& a_domain,
                          bool               a_cellAvgInput,
                          bool               a_fourthOrder) const
{
  CH_TIME("PatchCNSOp::computeWstate");
  CRD::msg << CRD::fv4 << "PatchCNSOp::computeWstate" << CRD::end;

  const int numUcomp  = CRDparam::g_CRDPhysics->numConservative();
  const int numWcomp  = CRDparam::g_CRDPhysics->numPrimitive();
  bool deconvCheck = ((a_fourthOrder) && (a_cellAvgInput) &&
                      (CRDparam::g_cellDeconvolveFlatten < 2)) ? true : false;
  bool cellDeconvLimiting =
    (CRDparam::g_cellDeconvolveLimit == 1) ? true : false;
  const int numGrowCells = deconvCheck ? (!cellDeconvLimiting ? 1 : 2) : 0;
  const Box consCellBox(grow(a_box, numGrowCells));

  CH_assert(a_WcellPntFab.nComp() == numWcomp);
  CH_assert(a_UcellAvgFab.nComp() == numUcomp);
  // CH_assert(a_UcellAvgFab.box().contains(consCellBox));
  CH_assert(a_WcellPntFab.box().contains(a_box));

  if (deconvCheck) // Deconvolve a_UcellAvgFab first
    {
      FABSTACKTEMP(UcellPntFab, a_box, numUcomp);
      CRDutil::deconvolve(UcellPntFab, a_UcellAvgFab, a_box, a_domain,
                          Interval(0, numUcomp-1), 4, 1, cellDeconvLimiting);

      CRDparam::g_CRDPhysics->consToPrim(a_WcellPntFab,
                                         UcellPntFab,
                                         a_box,
                                         a_WOld);
      CRDparam::g_CRDPhysics->consToPrimCorrectionLES(a_WcellPntFab, a_box);
      CRDparam::g_CRDPhysics->extraPrimitiveState(a_WcellPntFab, a_box);
    }
  else // No deconvolution necessary -- just go for it
    {
      CRDparam::g_CRDPhysics->consToPrim(a_WcellPntFab,
                                         a_UcellAvgFab,
                                         a_box,
                                         a_WOld);
      CRDparam::g_CRDPhysics->consToPrimCorrectionLES(a_WcellPntFab, a_box);
      CRDparam::g_CRDPhysics->extraPrimitiveState(a_WcellPntFab, a_box);
    }
}

/*--------------------------------------------------------------------*/
//  Compute the inertial (hyperbolic) flux
/** \param[in]  a_box   Flux determined on faces of these cells
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[in]  a_flux  Average flux on faces
 *  \param[out] a_flux  Updated with contribution from hyperbolic
 *                      terms
 *  \param[out] a_fluxFromWavg
 *                      Second-order estimate of flux based on average
 *                      of W on faces
 *  \param[in]  a_WfaceAvgFab
 *                      The average primitive state on the faces.
 *                      This is defined on the disjoint box grown by 2
 *                      except across non-periodic boundaries where it
 *                      is only grown by 0.
 *  \param[in]  a_WfaceAvgFab
 *  \param[out] a_WcellAvgFab
 *                      The average primitive state in the cells.
 *                      This is defined on the disjoint box grown by 4
 *                      except across non-periodic boundaries where it
 *                      is only grown by 2.
 *  \param[out] a_WcellAvgFab
 *  \param[out] a_WcellPntFab
 *                      The centered primitive state in the cells.
 *                      This is defined at the same places as a_WcellAvgFab
 *                      except is not needed across non-periodic boundaries.
 *  \param[in]  a_flattening
 *                      Flattening coefficients
 *  \param[out] a_flattening
 *                      Coefficients written if a_setFlattening is T
 *  \param[in]  a_bndryCellData
 *                      Previous-time data for cells immediately adjacent
 *                      to boundaries
 *  \param[in]  a_bndryNtJ
 *                      NtJ on near-bndry faces necessary for LES wall-model
 *  \param[in]  a_U     The average solution state in the cells in
 *                      physical space
 *  \param[in]  a_unitNormalFxb
 *                      Unit normal basis for transforming velocity
 *                      space to be normal to the faces
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
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::addInertialFlux(
  const Box&                          a_disjointBox,
  const Box&                          a_box,
  const BlockDomain&                  a_domain,
  FluxBox&                            a_flux,
  FluxBox&                            a_fluxFromWavg,
  FluxBox&                            a_WfaceAvgFxb,
  FluxBox&                            a_WfacePntFxb,
  FArrayBox&                          a_WcellAvgFab,
  const FArrayBox&                    a_WcellPntFab,
  FluxBox&                            a_faceAvgPlotFxb,
  FArrayBox&                          a_flattening,
  stc::Vector<FArrayBox, 2*SpaceDim>& a_bndryCellData,
  stc::Vector<FArrayBox, 2*SpaceDim>& a_bndryNtJ,
  const FArrayBox&                    a_U,
  const FluxBox&                      a_unitNormalFxb,
  const DataIndex&                    a_dataIndx,
  const Real                          a_time,
  const Real                          a_dt,
  const Real                          a_prevDt,
  const Real                          a_stageWeight,
  const bool                          a_setFlattening,
  const int                           a_level) const
{
  CH_TIME("PatchCNSOp::addInertialFlux");
  CRD::msg << CRD::fv4 << "PatchCNSOp::addInertialFlux (level: "
           << a_level << ")" << CRD::end;
  m_inertial4thOrderOp.flux(a_disjointBox,
                            a_box,
                            a_domain,
                            a_flux,
                            a_fluxFromWavg,
                            a_WfaceAvgFxb,
                            a_WfacePntFxb,
                            a_WcellAvgFab,
                            a_WcellPntFab,
                            a_faceAvgPlotFxb,
                            a_flattening,
                            a_bndryCellData,
                            a_bndryNtJ,
                            a_U,
                            a_unitNormalFxb,
                            a_dataIndx,
                            a_time,
                            a_dt,
                            a_prevDt,
                            a_stageWeight,
                            a_setFlattening,
                            a_level);
}

/*--------------------------------------------------------------------*/
//  Compute the face-average primitive state
/** This is the heart of the canonical inertial solve with limited
 *  reconstructionCompute the face-averaged
 *  \param[in]  a_box   Flux determined on faces of these cells
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[in]  a_WfaceAvgFab
 *                      The average primitive state on the faces.
 *                      This is defined on the disjoint box grown by 2
 *                      except across non-periodic boundaries where it
 *                      is only grown by 0.
 *  \param[in] a_WcellAvgFab
 *                      The average primitive state in the cells.
 *                      This is defined on the disjoint box grown by 4
 *                      except across non-periodic boundaries where it
 *                      is only grown by 2.
 *  \param[in]  a_WfromUavgFab
 *                      The average native (minimal) primitive state,
 *                      \<W\>, in cells computed directly from \<U\>.
 *                      It is only second-order accurage.
 *  \param[in]  a_flattening
 *                      Flattening coefficients
 *  \param[out] a_flattening
 *                      Coefficients written if a_setFlattening is T
 *  \param[in]  a_unitNormalFxb
 *                      Unit normal basis for transforming velocity
 *                      space to be normal to the faces
 *  \param[in]  a_dataIndx
 *                      Current data index
 *  \param[in]  a_time  Time (for applying BC)
 *  \param[in]  a_dt    Time step
 *  \param[in]  a_stageWeight
 *                      Weighting applied to contribution of flux
 *                      from this stage of RK4 to total flux
 *  \param[in]  a_setFlattening
 *                      T - Set flattening coefficients
 *  \param[in]  a_level Index of the AMR level
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::computeWfaceAvg(const Box&         a_disjointBox,
                            const Box&         a_box,
                            const BlockDomain& a_domain,
                            FluxBox&           a_WfaceAvgFxb,
                            const FArrayBox&   a_WcellAvgFab,
                            const FArrayBox&   a_WfromUavgFab,
                            FArrayBox&         a_flattening,
                            const FluxBox&     a_unitNormalFxb,
                            const DataIndex&   a_dataIndx,
                            const Real         a_time,
                            const Real         a_dt,
                            const Real         a_stageWeight,
                            const bool         a_setFlattening,
                            const int          a_level) const
{
  CH_TIME("PatchCNSOp::inertialFaceAverage");
  CRD::msg << CRD::fv4 << "PatchCNSOp::inertialFaceAverage (level: "
           << a_level << ")" << CRD::end;
  bool setFlattening = true;
  m_inertial4thOrderOp.computeFaceAverage(a_box,
                                          a_domain,
                                          a_flattening,
                                          a_WfaceAvgFxb,
                                          a_WcellAvgFab,
                                          a_WfromUavgFab,
                                          a_unitNormalFxb,
                                          a_time,
                                          setFlattening,
                                          a_level);
}

/*--------------------------------------------------------------------*/
//  Compute the viscous (elliptic) flux
/** \param[in]  a_box   Disjoint box
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[out] a_FluxfaceAvgFxb
 *                      Face-averaged viscous flux G
 *  \param[out] a_FluxfacePntFxb
 *                      Face-centered viscous flux G defined in domain +1 cells
 *  \param[out] a_turbSourceAvgFab
 *                      FAB containing mapped cell-averaged turbulent sources
 *  \param[out] a_invDtFab
 *                      Inverse time step values
 *  \param[in]  a_WcellAvgFab
 *                      Cell-averaged primitive values, is unchanged but
 *                      not const so we can use aliasing
 *  \param[in]  a_WcellPntFab
 *                      The cell-centered primitive variables, is unchanged but
 *                      not const so we can use aliasing
 *  \param[out] a_timeAvgFab
 *                      Time-averaged data for post-processing turbulent cases
 *  \param[in]  a_WfaceAvgFxb
 *                      Face-averaged primitive values, is unchanged but
 *                      not const so we can use aliasing
 *  \param[in]  a_WfacePntFxb
 *                      Face-centered primitive values, is unchanged but
 *                      not const so we can use aliasing
 *  \param[in]  a_facePntVelGradFxb
 *                      The face-centered velocity gradient
 *  \param[in]  a_facePntDeltaC
 *                      Face-centered cell-cutoff length for LES
 *  \param[in]  a_unitNormalsFxb
 *                      Face unit-normal vectors used for LES wall-model
 *  \param[in]  a_dataIndx
 *                      Data index for use with the metric terms
 *  \param[in]  a_dt    Time-step size on current level
 *  \param[in]  a_time  Time (for applying BC)
 *  \param[in]  a_timeAfterFilterInit
 *                      Total time after initializing time-averaging filter
 *  \param[in]  a_level Index of the AMR level
 *  \param[in]  a_minDiffDt
 *                      Previous minimum diffusive time step
 *  \param[out] a_minDiffDt
 *                      Updated minimum diffusive time step
 *  \param[out] a_minDiffDtCell
 *                      Updated cell with minimum diffusive time step
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::addViscousFlux(const Box&         a_box,
                           const BlockDomain& a_domain,
                           FluxBox&           a_FluxfaceAvgFxb,
                           FluxBox&           a_FluxfacePntFxb,
                           FArrayBox&         a_invDtFab,
                           FArrayBox&         a_turbSourceAvgFab,
                           FArrayBox&         a_WcellAvgFab,
                           FArrayBox&         a_WcellPntFab,
                           FArrayBox&         a_timeAvgFab,
                           FluxBox&           a_faceAvgPlotFxb,
                           FluxBox&           a_WfaceAvgFxb,
                           FluxBox&           a_WfacePntFxb,
                           FluxBox&           a_facePntVelGradFxb,
                           FluxBox&           a_stressFluxFxb,
                           const FluxBox&     a_facePntDeltaC,
                           const FluxBox&     a_faceCoord,
                           const FluxBox&     a_unitNormalsFxb,
                           const DataIndex&   a_dataIndx,
                           const Real         a_dt,
                           const Real         a_time,
                           const Real         a_timeAfterFilterInit,
                           const int          a_level,
                           Real&              a_minDiffDt,
                           IntVect&           a_minDiffDtCell) const
{
  CH_TIME("PatchCNSOp::addViscousFlux");
  CRD::msg << CRD::fv4 << "PatchCNSOp::addViscousFlux" << CRD::end;
  m_viscous4thOrderOp.flux(a_box,
                           a_domain,
                           a_FluxfaceAvgFxb,
                           a_FluxfacePntFxb,
                           a_invDtFab,
                           a_turbSourceAvgFab,
                           a_WcellAvgFab,
                           a_WcellPntFab,
                           a_timeAvgFab,
                           a_faceAvgPlotFxb,
                           a_WfaceAvgFxb,
                           a_WfacePntFxb,
                           a_facePntVelGradFxb,
                           a_stressFluxFxb,
                           a_facePntDeltaC,
                           a_faceCoord,
                           a_unitNormalsFxb,
                           a_dataIndx,
                           a_dt,
                           a_time,
                           a_timeAfterFilterInit,
                           a_level,
                           a_minDiffDt,
                           a_minDiffDtCell);
}

/*--------------------------------------------------------------------*/
//  Compute flux due to artificial viscosity and add to a_JUnew
/** \param[in]  a_box   Cell box of sites to be updated in a_JUnew.
 *                      Must be <= disjoint.
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[in]  a_JUnew
 *                      Solution in conserved state.
 *  \param[out] a_JUnew
 *                      Contribution of artificial viscosity added.
 *  \param[out] a_NtF   Flux due to artificial viscosity
 *  \param[in]  a_Uold  Solution at start of the time step
 *  \param[in]  a_N     Grid metrics on faces
 *  \param[in]  a_J     Metrics Jacobian in cells
 *  \param[in]  a_unitNormals
 *                      Unit normals on faces
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *  \param[in]  a_oldTime
 *                      Time at start of the time step
 *  \param[in]  a_dt    Weight for divergence operator
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::addArtificialViscosity(const Box&         a_box,
                                   const BlockDomain& a_domain,
                                   FArrayBox&         a_JUnew,
                                   FluxBox&           a_NtF,
                                   const FArrayBox&   a_Uold,
                                   const FArrayBox&   a_WOld,
                                   const FluxBox&     a_N,
                                   const FArrayBox&   a_J,
                                   const FluxBox&     a_unitNormals,
                                   const RealVect&    a_dx,
                                   const Real         a_dt,
                                   const Real         a_time,
                                   const int          a_level) const
{
  CH_TIME("PatchCNSOp::addArtificialViscosity");
  CRD::msg << CRD::fv4 << "PatchCNSOp::addArtificialViscosity" << CRD::end;
  CH_assert(CRDparam::g_useArtificialViscosity);
  Box box1Dom = grow(a_box, 1);
  box1Dom &= a_domain;

  // Compute the flux due to artificial viscosity
  CRDparam::g_CRDPhysics->artVisc(a_box,
                                  a_domain,
                                  a_NtF,
                                  a_Uold,
                                  a_WOld,
                                  a_N,
                                  a_J,
                                  a_unitNormals,
                                  a_dx,
                                  m_levelGridMetrics,
                                  a_time,
                                  a_level);

  // Get the divergence of the flux
  FABSTACKTEMP(divFlux, a_box, CRDparam::g_CRDPhysics->numFluxes());
  divFlux.setVal(0.);
  fluxDivergence(a_box, divFlux, a_NtF, a_dx);
  a_JUnew.plus(divFlux, a_dt);
}

/*--------------------------------------------------------------------*/
//  Compute flux due to artificial viscosity and add to a_JUnew
/** \param[in]  a_box   Cell box of sites to be updated in a_JUnew,
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[in]  a_JU    Solution in conserved state.
 *  \param[out] a_JU    Contribution of species correction
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::speciesCorrection(const Box&         a_box,
                              const BlockDomain& a_domain,
                              FArrayBox&         a_JU,
                              const RealVect&    a_dx) const
{
  CH_TIME("PatchCNSOp::speciesCorrection");
  CRD::msg << CRD::fv4 << "PatchCNSOp::speciesCorrection" << CRD::end;
  CRDparam::g_CRDPhysics->speciesCorrection(a_box,
                                            a_JU,
                                            a_domain,
                                            a_dx);
}

/*--------------------------------------------------------------------*/
//  Compute the divergence of a flux
/** \param[in]  a_box   Cell box of sites to be updated
 *  \param[out] a_rhs   Updated with divergence of the flux
 *  \param[in]  a_flux  Flux on faces
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::fluxDivergence(const Box&      a_box,
                           FArrayBox&      a_rhs,
                           const FluxBox&  a_flux,
                           const RealVect& a_dx) const
{
  CH_TIME("PatchCNSOp::fluxDivergence");
  CRD::msg << CRD::fv4 << "PatchCNSOp::fluxDivergence" << CRD::end;
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      const FArrayBox& fluxDir = a_flux[dir];
      // Increment a_rhs on all cells of a_box
      const Real dx = a_dx[0];  //**FIXME-dxVect
      FORT_DIVERGENCE(CHF_CONST_FRA(fluxDir),
                      CHF_FRA(a_rhs),
                      CHF_BOX(a_box),
                      CHF_CONST_REAL(dx),
                      CHF_INT(dir));
    }
  a_rhs.negate(a_box, 0, a_flux.nComp());
}

/*--------------------------------------------------------------------*/
//  Update flux registers
/** \param[out] a_crFluxRegister
 *                      Register with the next coarser level
 *  \param[out] a_fnFluxRegister
 *                      Register with the next finer level
 *  \param[in]  a_flux  Flux on the faces to add/subtract from
 *                      registers
 *  \param[in]  a_box   Disjoint box region to update the flux on
 *  \param[in]  a_dataIdx
 *                      Data iterator
 *  \param[in]  a_weight
 *                      Weight for divergence (essentially dt)
 *  \param[in]  a_hasCoarserGrid
 *                      T - A coarser level exists and is in use
 *  \param[in]  a_hasFinerGrid
 *                      T - A finer level exists and is in use
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::updateFluxRegisters(LevelFluxRegister& a_crFluxRegister,
                                LevelFluxRegister& a_fnFluxRegister,
                                const FluxBox&     a_flux,
                                const Box&         a_box,
                                const DataIndex&   a_dataIdx,
                                const Real         a_weight,
                                const bool         a_hasCoarserGrid,
                                const bool         a_hasFinerGrid) const
{
  CH_TIME("PatchCNSOp::updateFluxRegisters");
  CRD::msg << CRD::fv4 << "PatchCNSOp::updateFluxRegisters" << CRD::end;
  Interval fullFluxIntv(0, CRDparam::g_CRDPhysics->numFluxes() - 1);

  // Do flux register updates
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      // This gets shifted and shifted back, with no data change.
      FArrayBox& fluxDir = const_cast<FArrayBox&>(a_flux[dir]);

      // Increment fine flux registers between this level and the next
      // coarser level.
      if (a_hasCoarserGrid)
        {
          CH_assert(a_crFluxRegister.isDefined());
          a_crFluxRegister.incrementFine(fluxDir,
                                         a_weight,
                                         a_dataIdx,
                                         fullFluxIntv,
                                         fullFluxIntv,
                                         dir);
        }

      // Increment coarse flux register between this level and the next
      // finer level.
      if (a_hasFinerGrid)
        {
          CH_assert(a_fnFluxRegister.isDefined());
          a_fnFluxRegister.incrementCoarse(fluxDir,
                                           a_weight,
                                           a_dataIdx,
                                           fullFluxIntv,
                                           fullFluxIntv,
                                           dir,
                                           a_box);
        }
    }
}

/*--------------------------------------------------------------------*/
//  Add source terms
/** \param[out] a_sourceFab
 *                      Cell-centered source terms
 *  \param[out] a_invDtFab
 *                      Inverse of the time step
 *  \param[in]  a_Wcell Cell-centered primitive variables
 *  \param[in]  a_UcellAvg 
 *                      Cell-averaged conservative variables
 *  \param[in]  a_WfaceAvgFxb
 *                      Face-averaged primitive state
 *  \param[in]  a_time  Solution time
 *  \param[in]  a_stageWeight
 *                      Contribution of this stage to the final flux
 *                      over the time step
 *  \param[in]  a_level Current level
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[in]  a_domain
 *                      The problem or block domain for the box
 *  \param[in]  a_solveBox
 *                      Box to solve over
 *  \param[in]  a_globalKE
 *                      Global mean kinetic energy (turbulent forcing)
 *  \param[in]  a_globalHelicity
 *                      Global mean helicity (turbulent forcing)
 *//*-----------------------------------------------------------------*/

void 
PatchCNSOp::addSourceTerm(FArrayBox&         a_sourceFab,
                          FArrayBox&         a_invDtFab,
                          const FArrayBox&   a_Wcell,
                          const FArrayBox&   a_UcellAvg,
                          const FluxBox&     a_WfaceAvgFxb,
                          const Real         a_time,
                          const Real         a_stageWeight,
                          const int          a_level,
                          const Box&         a_disjointBox,
                          const BlockDomain& a_domain,
                          const Box&         a_solveBox,
                          const DataIndex&   a_dataIndx,
                          const Real         a_globalKE,
                          const Real         a_globalHelicity) const
{
  CH_TIME("PatchCNSOp::addSourceTerm");
  CRD::msg << CRD::fv4 << "PatchCNSOp::addSourceTerm" << CRD::end;
  CH_assert(a_invDtFab.contains(a_solveBox));
  CRDparam::g_CNSIBC->addSourceTerm(a_sourceFab,
                                    a_invDtFab,
                                    a_Wcell,
                                    a_UcellAvg,
                                    a_WfaceAvgFxb,
                                    a_domain,
                                    m_levelGridMetrics,
                                    a_time,
                                    a_stageWeight,
                                    a_level,
                                    a_disjointBox,
                                    a_solveBox,
                                    a_dataIndx,
                                    a_globalKE,
                                    a_globalHelicity);
}

/*--------------------------------------------------------------------*/
//  Solves for the mapped cell-averaged values from the cell-centered
//  values
/** \param[out] a_cellAvgJS
 *                      FAB containing the mapped cell-averaged values
 *  \param[in]  a_cellPntS
 *                      FAB containing the cell-centered values
 *  \param[out] a_cellPntS
 *                      FAB filled with mapped cell-centered values
 *  \param[in]  a_boxJS Box for solving <JS>
 *  \param[in]  a_boxS  Box for solving (JS), 1 bigger than a_boxJS
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[in]  a_dataIdx
 *                      Data iterator
 *  \param[in]  a_domain
 *                      The problem or block domain for the disjoint
 *                      box
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::solveJSfromS(FArrayBox&         a_cellAvgJS,
                         FArrayBox&         a_cellPntS,
                         const Box&         a_boxJS,
                         const Box&         a_boxS,
                         const DataIndex&   a_dataIndx,
                         const BlockDomain& a_domain,
                         const bool         a_fourthOrder) const
{
  CH_TIME("PatchCNSOp::solveJSfromS");
  CRD::msg << CRD::fv4 << "PatchCNSOp::solveJSfromS" << CRD::end;

  const int numComp = a_cellPntS.nComp();
  CH_assert(numComp == a_cellAvgJS.nComp());
  FABSTACKTEMP(cellAvgS, a_boxJS, numComp);

  // Convolve S to get <S>
  const int convOrder = (a_fourthOrder &&
                         (CRDparam::g_cellConvolveFlatten < 2)) ? 4 : 2;

//**FIXME THIS CONVOLUTION HAS TO BE LIMITED IF NOT USING CHEMICAL TIME STEP!!!
  CRDutil::deconvolve(cellAvgS, a_cellPntS, a_boxJS, a_domain,
                      Interval(0, numComp - 1), convOrder, -1, false);

  const FArrayBox& cellAvgJ = m_levelGridMetrics.m_J[a_dataIndx];
  FABSTACKTEMP(gradJ, a_boxJS, SpaceDim);
  FABSTACKTEMP(gradS, a_boxJS, numComp*SpaceDim);

  if (a_fourthOrder)
    {
      // Compute gradients of J
      PatchMappedFunc::undividedGrad2OCS(a_boxJS, gradJ, cellAvgJ, a_domain,
                                         Interval(0, 0));
      // Compute gradients of S
      PatchMappedFunc::undividedGrad2OCS(a_boxJS, gradS, a_cellPntS, a_domain,
                                         Interval(0, numComp - 1));
    }

  for (int comp = 0; comp != numComp; ++comp)
    {
      // Fill a_cellAvgJS with <J><S> + (1/12)*\sum_i((dJ/dx_i)*(dS/dx_i))
      MD_BOXLOOP(a_boxJS, i)
        {
          a_cellAvgJS[MD_IX(i, comp)] =
            cellAvgJ[MD_IX(i, 0)]*cellAvgS[MD_IX(i, comp)];
        }

      if (a_fourthOrder)
        {
          // Multiply the gradients of each component and add to a_cellAvgJS
          const Real factor = 1./12.;
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              const int sComp = comp*SpaceDim + dir;
              MD_BOXLOOP(a_boxJS, i)
                {
                  const Real dJ1 = gradJ[MD_IX(i, dir)];
                  const Real dS1 = gradS[MD_IX(i, sComp)];
                  a_cellAvgJS[MD_IX(i, comp)] += factor*dJ1*dS1;
                }
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Solves for the max norm of JU-\hat{JU}
/** \param[in]  a_JU    Solution vector
 *  \param[in]  a_JUhat Solution vector associated with the
 *                      embedded ARK scheme
 *  \param[in]  a_disjointBox
 *                      Disjoint box
 *  \param[in]  a_interval
 *                      Which components to compute the norm over
 *  \return The max norm of U-\hat{U}
 *//*-----------------------------------------------------------------*/

Real
PatchCNSOp::computeUminusUhatMaxNorm(const FArrayBox& a_JU,
                                     const FArrayBox& a_JUhat,
                                     const Box&       a_disjointBox,
                                     const Interval&  a_interval) const
{
  const int numComps = a_interval.size();
  std::vector<Real> speciesNorms(numComps);
  FABSTACKTEMP(UminusUHat, a_disjointBox, numComps);
  for (int comp = 0; comp != numComps; ++comp)
    {
      const int UcompIdx = a_interval.begin() + comp;
      MD_BOXLOOP(a_disjointBox, i)
        {
          UminusUHat[MD_IX(i, comp)] =
            a_JU[MD_IX(i, UcompIdx)]-a_JUhat[MD_IX(i, UcompIdx)];
        }

      speciesNorms[comp] = UminusUHat.norm(
          a_disjointBox, // box to compute norm on
          0, // Which norm to do (0 = max norm)
          comp, // Which component to compute norm of
          1); // number of components
    }

  // Now compute the largest magnitude of all of the components
  Real maxNorm = speciesNorms[0];
  for (int j = 0; j != numComps; ++j)
    {
      maxNorm = std::max(maxNorm, speciesNorms[j]);
    }

  return maxNorm;
}

/*--------------------------------------------------------------------*/
//  Things to do before calling the Riemann solver
/** Specifically, transform the velocity to be normal to the face.
 *  **FIXME Move to PatchMappedFunc and comments
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::preRiemann(FArrayBox&     a_WLeft,
                       FArrayBox&     a_WRight,
                       const FluxBox& a_unitNormalFxb,
                       const int      a_dir,
                       const Box&     a_box)
{
  CH_TIME("PatchCNSOp::preRiemann");
  CRD::msg << CRD::fv4 << "PatchCNSOp::preRiemann" << CRD::end;
  Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
  if (velIntv.begin() <= velIntv.end())
    {
      FArrayBox& shiftWLeft  = (FArrayBox&)a_WLeft;
      FArrayBox& shiftWRight = (FArrayBox&)a_WRight;
      shiftWLeft .shiftHalf(a_dir, 1);
      shiftWRight.shiftHalf(a_dir,-1);
      CH_assert(shiftWLeft .box().contains(a_box));
      CH_assert(shiftWRight.box().contains(a_box));

      FArrayBox velLeftFab(velIntv, shiftWLeft);    // Alias
      FArrayBox velRightFab(velIntv, shiftWRight);  // Alias
      PatchCNSOp::forwardBasisTransform(a_box,
                                        velLeftFab,
                                        a_unitNormalFxb,
                                        a_dir);
      PatchCNSOp::forwardBasisTransform(a_box,
                                        velRightFab,
                                        a_unitNormalFxb,
                                        a_dir);

      shiftWLeft .shiftHalf(a_dir,-1);
      shiftWRight.shiftHalf(a_dir, 1);
    }
}

/*--------------------------------------------------------------------*/
//  Things to do after calling the Riemann solver
/** Specifically, untransform the star velocity state.
 *  **FIXME Move to PatchMappedFunc and comments
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::postRiemann(FArrayBox&     a_Wface,
                                const FluxBox& a_unitNormalFxb,
                                const int      a_dir,
                                const Box&     a_box)
{
  CH_TIME("PatchCNSOp::postRiemann");
  CRD::msg << CRD::fv4 << "PatchCNSOp::postRiemann" << CRD::end;
  // Transform back velocity components in a_Wface.
  Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
  if (velIntv.begin() <= velIntv.end())
    {
      FArrayBox velFab(velIntv, a_Wface);
      PatchCNSOp::reverseBasisTransform(a_box, velFab, a_unitNormalFxb, a_dir);
    }
}

/*--------------------------------------------------------------------*/
//  Forward basis transform on mapped grids
//**FIXME Move to PatchMappedFunc
/**
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::forwardBasisTransform(const Box&     a_box,
                                  FArrayBox&     a_W,
                                  const FluxBox& a_unitNormalFxb,
                                  const int      a_dir)
{
  CH_TIME("PatchCNSOp::forwardBasisTransform");
  CRD::msg << CRD::fv4 << "PatchCNSOp::forwardBasisTransform" << CRD::end;
  const FArrayBox& unitNormalDirFab = a_unitNormalFxb[a_dir];
  CH_assert(unitNormalDirFab.box().contains(a_box));
  PatchMappedFunc::forwardTransform(a_W,
                                    unitNormalDirFab,
                                    a_box);
}

/*--------------------------------------------------------------------*/
//  Forward basis transform on mapped grids
//**FIXME Move to PatchMappedFunc
/**
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::reverseBasisTransform(const Box&     a_box,
                                  FArrayBox&     a_W,
                                  const FluxBox& a_unitNormalFxb,
                                  const int      a_dir)
{
  CH_TIME("PatchCNSOp::reverseBasisTransform");
  CRD::msg << CRD::fv4 << "PatchCNSOp::reverseBasisTransform" << CRD::end;
  const FArrayBox& unitNormalDirFab = a_unitNormalFxb[a_dir];
  CH_assert(unitNormalDirFab.box().contains(a_box));
  PatchMappedFunc::reverseTransform(a_W,
                                    unitNormalDirFab,
                                    a_box);
}

/*--------------------------------------------------------------------*/
//  Compute the primitive cell-centered state from the conservative
/** \param[out]  a_WcellPntFab
 *                      The primitive cell point state
 *  \param[in]   a_box  Disjoint box region to compute the state
 *  \param[in]   a_domain
 *                      The problem or block domain for the box
 *  \param[in]   a_UcellAvgFab
 *                      The conservative cell average state
 *  \param[in]   a_cellAvgInput
 *                      If true, must deconvolve UcellAvg first
 *                      If false, nothing must be done
 *  \param[in]   a_fourthOrder
 *                      If true, see a_cellAvgInput
 *                      If false, disregard a_cellAvgInput and just
 *                      go ahead with the consToPrim conversion
 *//*-----------------------------------------------------------------*/

void
PatchCNSOp::computeWpntCell(FArrayBox&         a_WcellPntFab,
                            const FArrayBox&   a_UcellAvgFab,
                            const FArrayBox&   a_WOld,
                            const Box&         a_box,
                            const BlockDomain& a_domain,
                            bool               a_cellAvgInput,
                            bool               a_fourthOrder)
{
  CH_TIME("PatchCNSOp::computeWpntCell");
  CRD::msg << CRD::fv4 << "PatchCNSOp::computeWpntCell" << CRD::end;

  const int numUcomp  = CRDparam::g_CRDPhysics->numConservative();
  const int numWcomp  = CRDparam::g_CRDPhysics->numPrimitive();
  bool deconvCheck = ((a_fourthOrder) && (a_cellAvgInput) &&
                      (CRDparam::g_cellDeconvolveFlatten < 2)) ? true : false;
  bool cellDeconvLimiting =
    (CRDparam::g_cellDeconvolveLimit == 1) ? true : false;
  const int numGrowCells = deconvCheck ? (cellDeconvLimiting ? 2 : 1) : 0;
  const Box consCellBox(grow(a_box, numGrowCells));

  CH_assert(a_WcellPntFab.nComp() == numWcomp);
  CH_assert(a_UcellAvgFab.nComp() == numUcomp);
  CH_assert(a_UcellAvgFab.box().contains(consCellBox));
  CH_assert(a_WcellPntFab.box().contains(a_box));

  if (deconvCheck) // Deconvolve a_UcellAvgFab first
    {
      FABSTACKTEMP(UcellPntFab, a_box, numUcomp);
      CRDutil::deconvolve(UcellPntFab, a_UcellAvgFab, a_box, a_domain,
                          Interval(0, numUcomp-1), 4, 1, cellDeconvLimiting);

      CRDparam::g_CRDPhysics->consToPrim(
        a_WcellPntFab, UcellPntFab, a_box, a_WOld);
      CRDparam::g_CRDPhysics->consToPrimCorrectionLES(a_WcellPntFab, a_box);
      CRDparam::g_CRDPhysics->extraPrimitiveState(a_WcellPntFab, a_box);
    }
  else // No deconvolution necessary -- just go for it
    {
      CRDparam::g_CRDPhysics->consToPrim(
        a_WcellPntFab, a_UcellAvgFab, a_box, a_WOld);
      CRDparam::g_CRDPhysics->consToPrimCorrectionLES(a_WcellPntFab, a_box);
      CRDparam::g_CRDPhysics->extraPrimitiveState(a_WcellPntFab, a_box);
    }
}
