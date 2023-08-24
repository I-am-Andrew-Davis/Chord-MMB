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
 * \file CNSIBCRecirculatingInletTFP.H
 *
 * \brief Initial and boundary conditions for recirculating turbulent flat plate
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCRECIRCULATINGINLETTFP_H_
#define _CNSIBCRECIRCULATINGINLETTFP_H_

//----- Internal -----//

#include "CNSIBCGeneralized.H"
#include "RootSolver.H"

/*******************************************************************************
 */
///  BC and initial conditions on the grid for turbulent flat plate
/**
 *//*+*************************************************************************/

class CNSIBCRecirculatingInletTFP : public CNSIBCGeneralized
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCRecirculatingInletTFP();

  /// Destructor
  virtual ~CNSIBCRecirculatingInletTFP();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCRecirculatingInletTFP(const CNSIBCRecirculatingInletTFP&);

//--Assignment not permitted

  CNSIBCRecirculatingInletTFP &operator=(const CNSIBCRecirculatingInletTFP&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Return a name describing the IBC
  virtual const char *const IBCName() const;

  /// Write any information related to the IBC to output
  virtual void writeIBCInfo() const;

  /// Set the tagging method if one can be associated with IBC
  virtual TagLevelFactory* setTagMethod(const int a_tagBufferSize);

  /// Initialize a level
  virtual void initialize(LevelData<FArrayBox>&      a_U,
                          LevelGridMetrics&          a_gridMetrics,
                          const LayoutData<FluxBox>& a_unitNormals,
                          const Real                 a_time,
                          const int                  a_level) const;

  /// Initialize the inlet region
  virtual void initializeInletDataStructures(
    const LevelData<FArrayBox>& a_U,
    const int                   a_level,
    const DisjointBoxLayout&    a_disjointBoxLayout,
    const LevelGridMetrics&     a_gridMetrics,
    const bool                  a_hasFinerGrid) const;

  /// Copy the interior plane to the inlet plane and rescale it
  virtual void copyInteriorToInletAndRescale(
    const LevelData<FArrayBox>& a_U,
    const int                   a_level,
    const Real                  a_t,
    const Real                  a_dt,
    const int                   a_stage,
    const DisjointBoxLayout&    a_disjointBoxLayout,
    const LevelGridMetrics&     a_gridMetrics) const;

  /// Add source term
  virtual void addSourceTerm(FArrayBox&           a_sourceFab,
                             FArrayBox&           a_invDtFab,
                             const FArrayBox&     a_Wcell,
                             const FArrayBox&     a_UcellAvg,
                             const FluxBox&       a_WfaceAvgFxb,
                             const ProblemDomain& a_domain,
                             LevelGridMetrics&    a_gridMetrics,
                             const Real           a_time,
                             const Real           a_stageWeight,
                             const int            a_level,
                             const Box&           a_disjointBox,
                             const Box&           a_solveBox,
                             const DataIndex&     a_dataIndx,
                             const Real           a_GlobalKE,
                             const Real           a_globalHelicity) const;

  /// Spatially average data if necessary
  virtual void spatiallyAverageData(
    LevelData<FArrayBox>& a_avgData,
    const int             a_stage,
    const Real            a_t) const;

//--Patch operators

  /// Does an exact solution exist?
  virtual bool haveExactSol() const;

/*--------------------------------------------------------------------*/
//  Non-linear function for wall friction-velocity solution (Urbin 2001 p. 1292)
/**
 *  \param[in]  a_a     Constant
 *  \param[in]  a_c     Constant
 *  \param[in]  a_k     Kolmogorov constant
 *  \param[in]  a_delta Local boundary layer thickness (delta_99)
 *  \param[in]  a_nu    Kinematic viscosity at wall
 *  \param[in]  a_uInf  Freestream velocity (untransformed)
 *  \param[in]  a_ColesParam
 *                      Coles' parameter for turbulent boundary layers
 *//*-----------------------------------------------------------------*/

  struct UTauMuskerFunc
  {
    UTauMuskerFunc(const Real a_delta,
                   const Real a_nu,
                   const Real a_uInf,
                   const Real a_ColesParam)
      :
      m_delta(a_delta),
      m_nu(a_nu),
      m_uInf(a_uInf),
      m_ColesParam(a_ColesParam)
      {}
    Real operator()(const Real& a_uTau) const
      {
        Real delta_plus = m_delta*a_uTau/m_nu;
        const Real b_0 = 5.424*std::atan2((2.*delta_plus - 8.15), 16.7);
        const Real b_1 = std::pow((delta_plus + 10.6), 9.6);
        const Real b_2 = (delta_plus*delta_plus - 8.15*delta_plus + 86);
        const Real b_3 = std::log10(b_1/(b_2*b_2));
        return a_uTau*(b_0 + b_3 - 3.52 + 4.88*m_ColesParam) - m_uInf;
      }
    const Real m_delta;
    const Real m_nu;
    const Real m_uInf;
    const Real m_ColesParam;
  };

  struct deltaFunc
  {
    deltaFunc(const Real a_delta,
              const Real a_nu,
              const Real a_uInf,
              const Real a_ColesParam,
              const Real a_yPlusMax)
      :
      m_delta99(a_delta),
      m_nu2(a_nu),
      m_uInf2(a_uInf),
      m_ColesParam2(a_ColesParam),
      m_yPlusMax2(a_yPlusMax)
      {}
    Real operator()(const Real& a_delta100) const
      {
        // First, solve for u_tau with the current estimate for a_delta100
        int iterBrent = 0;
        int errorBrent = 0;
        Real uTauMin = 0.;
        Real uTauMax = m_yPlusMax2*m_nu2/a_delta100;
        const UTauMuskerFunc& f =
          UTauMuskerFunc(a_delta100, m_nu2, m_uInf2, m_ColesParam2);
        Real u_tau =
          RootSolver::BrentER(iterBrent, errorBrent, f, uTauMin, uTauMax);
        if (errorBrent != 0 || u_tau != u_tau)
          {
            CRD::msg << "deltaFunc: Bad uTau value: " << u_tau << CRD::error;
          }
        // Next, evaluate the Musker profile at delta_99

        Real delta_plus = m_delta99*u_tau/m_nu2;
        const Real b_0 = 5.424*std::atan2((2.*delta_plus - 8.15), 16.7);
        const Real b_1 = std::pow((delta_plus + 10.6), 9.6);
        const Real b_2 = (delta_plus*delta_plus - 8.15*delta_plus + 86);
        const Real b_3 = std::log10(b_1/(b_2*b_2));
        // Add in the wake component
        const Real b_4 = m_delta99/a_delta100;
        const Real b_5 = b_4*b_4;
        const Real b_6 = b_5*b_4;
        const Real b_7 =
          2.44*(m_ColesParam2*(6.*b_5 - 4.*b_6) + b_5*(1. - b_4));
        // Evaluate the full profile at delta_99
        return u_tau*(b_0 + b_3 - 3.52 + b_7) - 0.99*m_uInf2;
      }
    const Real m_delta99;
    const Real m_nu2;
    const Real m_uInf2;
    const Real m_ColesParam2;
    const Real m_yPlusMax2;
  };

  struct yMatchingFunc
  {
    yMatchingFunc(const Real a_deltaInlet,
                  const Real a_deltaRecy,
                  const Real a_yInlet,
                  const Real a_uTauInlet,
                  const Real a_uTauRecy,
                  const Real a_nu,
                  const Real a_ColesParam)
      :
      m_delta100Inlet(a_deltaInlet),
      m_delta100Recy(a_deltaRecy),
      m_yInlet(a_yInlet),
      m_uTauInlet(a_uTauInlet),
      m_uTauRecy(a_uTauRecy),
      m_nu3(a_nu),
      m_ColesParam3(a_ColesParam)
      {}
    Real operator()(const Real& a_yRecy) const
      {
        Real y_p_inlet = m_yInlet*m_uTauInlet/m_nu3;
        Real y_p_recy = a_yRecy*m_uTauRecy/m_nu3;

        const Real a_0 = 5.424*std::atan2((2.*y_p_inlet - 8.15), 16.7);
        const Real a_1 = std::pow((y_p_inlet + 10.6), 9.6);
        const Real a_2 = (y_p_inlet*y_p_inlet - 8.15*y_p_inlet + 86);
        const Real a_3 = std::log10(a_1/(a_2*a_2));
        const Real a_4 = m_yInlet/m_delta100Inlet;
        const Real a_5 = a_4*a_4;
        const Real a_6 = a_5*a_4;
        const Real a_7 =
          2.44*(m_ColesParam3*(6.*a_5 - 4.*a_6) + a_5*(1. - a_4));

        const Real b_0 = 5.424*std::atan2((2.*y_p_recy - 8.15), 16.7);
        const Real b_1 = std::pow((y_p_recy + 10.6), 9.6);
        const Real b_2 = (y_p_recy*y_p_recy - 8.15*y_p_recy + 86);
        const Real b_3 = std::log10(b_1/(b_2*b_2));
        const Real b_4 = a_yRecy/m_delta100Recy;
        const Real b_5 = b_4*b_4;
        const Real b_6 = b_5*b_4;
        const Real b_7 =
          2.44*(m_ColesParam3*(6.*b_5 - 4.*b_6) + b_5*(1. - b_4));

        return m_uTauRecy*(b_0+b_3-3.52+b_7) - m_uTauInlet*(a_0+a_3-3.52+a_7);
      }
    const Real m_delta100Inlet;
    const Real m_delta100Recy;
    const Real m_yInlet;
    const Real m_uTauInlet;
    const Real m_uTauRecy;
    const Real m_nu3;
    const Real m_ColesParam3;
  };

/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:
  
  /// Set the imposed BC for primitive variables
  virtual void setImposedBCprimState(
    FArrayBox&                    a_Wface,
    const Box&                    a_boundaryFaceBox,
    const FArrayBox&              a_Wcell,
    const FArrayBox&              a_unitNormalBasisFab,
    const BoundaryIndex&          a_bcIdx,
    const Box&                    a_disjointBox,
    LevelGridMetrics&             a_gridMetrics,
    const Real                    a_time,
    const int                     a_level,
    const BCInfo&                 a_bcInfo) const;

  /// Set exterior face-state for relaxed characteristic boundary conditions
  virtual void setRelaxedCBCPrimState(
    FArrayBox&              a_WfaceAvgExterior,
    FArrayBox&              a_bndryCellFab,
    FArrayBox&              a_WcellAvgFab,
    const FArrayBox&        a_WfaceAvgDirFab,
    const FArrayBox&        a_unitNormalBasisFab,
    const BoundaryIndex&    a_bcIdx,
    const BCInfo&           a_bcInfo,
    const Box&              a_disjointBox,
    const Box&              a_totalFaceBox,
    const LevelGridMetrics& a_gridMetrics,
    const int               a_dir,
    const Side::LoHiSide    a_side,
    const Real              a_time,
    const Real              a_prevDt,
    const int               a_level) const;

/*==============================================================================
 * Private member functions
 *============================================================================*/

private:

  /// Read any information related to the IBC from input
  void readBCInfo();

  /// Space-time average the flow variables at the recycling plane
  void spaceTimeAverageVariables(const FArrayBox& a_WcellPntFab,
                                 const Box        a_samplePlaneBox,
                                 const Real       a_t,
                                 const Real       a_dt,
                                 const int        a_level) const;

  /// Compute the inlet recycled-state consisting of fluctuations
  void recycleFluctuations(FArrayBox&       a_inletPlaneW,
                           const FArrayBox& a_meanInletFab,
                           const FArrayBox& a_WcellPntFab,
                           const FArrayBox& a_WcellInletSampledPntFab,
                           const Box&       a_spBox,
                           const Box&       a_ipBox,
                           const int        a_level,
                           const Real       a_t) const;

  /// Compute an approximate initial condition based on Musker 1979
  void muskerIC(FArrayBox&       a_meanFab,
                FArrayBox&       a_etaFab,
                const FArrayBox& a_deltaFab,
                const FArrayBox& a_XFab,
                const FArrayBox& a_XNodeFab,
                const Real       a_virtualWallDy,
                const Box&       a_box) const;

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  std::vector<std::unique_ptr<LevelData<FArrayBox>>> m_singleBoxSampledPlane;
                                        ///< Data from sampled-plane residing in
                                        ///< the space of one box
  std::vector<
    std::unique_ptr<LevelData<FArrayBox>>> m_singleBoxSampledInletPlane;
                                        ///< Data from inlet-plane residing in
                                        ///< the space of one box
  std::vector<std::unique_ptr<LevelData<FArrayBox>>> m_singleBoxInletPlane;
                                        ///< Scaled data moved to inlet in the
                                        ///< space of one box
  std::vector<
    std::unique_ptr<LevelData<FArrayBox>>> m_singleBoxInletPlaneRelaxed;
                                        ///< Scaled data moved to inlet in the
                                        ///< space of one box
  std::vector<
    std::unique_ptr<LevelData<FArrayBox>>> m_singleBoxInletPlaneInst;
                                        ///< Scaled data moved to inlet in the
                                        ///< space of one box
  std::vector<std::unique_ptr<LevelData<FArrayBox>>> m_singleBoxMeanInletPlane;
                                        ///< Mean data imposed at inlet in the
                                        ///< space of one box -- Level Data
                                        ///< Note: this lives at the sampling
                                        ///< plane location for convenience
  std::vector<std::unique_ptr<LevelData<FArrayBox>>> m_multiBoxInletPlane;
                                        ///< Scaled data residing on normal
                                        ///< disjointBoxLayout of m_U data in
                                        ///< LevelCNSOp
  std::vector<
    std::unique_ptr<LevelData<FArrayBox>>> m_multiBoxInletPlaneRelaxed;
                                        ///< Scaled data residing on normal
                                        ///< disjointBoxLayout of m_U data in
                                        ///< LevelCNSOp
  std::vector<
    std::unique_ptr<LevelData<FArrayBox>>> m_multiBoxInletPlaneInst;
                                        ///< Scaled data residing on normal
                                        ///< disjointBoxLayout of m_U data in
                                        ///< LevelCNSOp
  std::vector<std::unique_ptr<Copier>> m_singleBoxSampledPlaneCopier;
                                        ///< Copier from global dataset to
                                        ///< m_singleBoxSampledPlane
  std::vector<std::unique_ptr<Copier>> m_singleBoxSampledInletPlaneCopier;
                                        ///< Copier from global dataset to
                                        ///< m_singleBoxSampledInletPlane
  std::vector<std::unique_ptr<Copier>> m_multiBoxInletPlaneCopier;
                                        ///< Copier from m_singleBoxInletPlane
                                        ///< to m_multiBoxInletPlane
  std::vector<
    std::unique_ptr<Copier>> m_multiBoxInletPlaneRelaxedCopier;
                                        ///< Copier from m_singleBoxInletPlane
                                        ///< to m_multiBoxInletPlane
  std::vector<
    std::unique_ptr<Copier>> m_multiBoxInletPlaneInstCopier;
                                        ///< Copier from m_singleBoxInletPlane
                                        ///< to m_multiBoxInletPlane
  std::vector<std::unique_ptr<Copier>> m_multiBoxInletExchangeCopier;
                                        ///< Exchange copier for inlet LevelData
  mutable std::vector<std::vector<Real>> m_uVelSTAvg;
                                        ///< X-velocity average at sample-plane
                                        ///< (space-time average)
  mutable std::vector<std::vector<Real>> m_uVelSAvg;
                                        ///< X-velocity average at sample-plane
                                        ///< (space average)
  mutable std::vector<std::vector<Real>> m_vVelSTAvg;
                                        ///< Y-velocity average at sample-plane
                                        ///< (space-time average)
  mutable std::vector<std::vector<Real>> m_vVelSAvg;
                                        ///< Y-velocity average at sample-plane
                                        ///< (space average)
  mutable std::vector<std::vector<Real>> m_wVelSAvg;
                                        ///< Z-velocity average at sample-plane
                                        ///< (space average)
  mutable std::vector<std::vector<Real>> m_rhoSTAvg;
                                        ///< Density average at sample-plane
                                        ///< (space-time average)
  mutable std::vector<std::vector<Real>> m_rhoSAvg;
                                        ///< Density average at sample-plane
                                        ///< (space average)
  mutable std::vector<std::vector<Real>> m_tempSTAvg;
                                        ///< Temperature average at sample-plane
                                        ///< (space-time average)
  mutable std::vector<std::vector<Real>> m_tempSAvg;
                                        ///< Temperature average at sample-plane
                                        ///< (space average)
  mutable std::vector<std::vector<Real>> m_yLoc;
                                        ///< Y-location in 1D vertical box
  mutable std::vector<std::vector<Real>> m_yLocMapped;
                                        ///< Y-location in 1D vertical box
  mutable std::vector<std::vector<Real>> m_zLoc;
                                        ///< Z-location in 1D horizontal box
  mutable std::vector<Real> m_deltaRecy;///< BL-thickness at recycling plane
                                        ///< on current level. Needed for
                                        ///< levels that don't extend through BL
  Real m_delta;                         ///< Boundary layer thickness (99%)
  Real m_alpha;                         ///< Boundary-layer blending coeff 1
  Real m_beta;                          ///< Boundary-layer blending coeff 2
  Real m_phi1;                          ///< Number of characteristic times to
                                        ///< use in start-up time-averaging
  Real m_phi2;                          ///< Number of characteristic times to
                                        ///< use in transition time-averaging
  Real m_psi1;                          ///< Number of characteristic times
                                        ///< before transitioning to m_phi2
  Real m_psi2;                          ///< Number of characteristic times
                                        ///< before transitioning to running
                                        ///< time-average (for rest of run)
  Real m_lambda;                        ///< Characteristic time-scale at inlet
  Real m_samplePlaneLoc;                ///< Location of sampling plane for
                                        ///< recycled inlet in terms of percent
                                        ///< of the domain x-direction length
  Real m_yPlusMaxGuess;                 ///< Guess of maximum likely y-plus
                                        ///< value in the first cell adjacent
                                        ///< to the wall
  Real m_velPerturb;                    ///< Magnitude of velocity perturbations
  RealVect m_perturbFreq;               ///< Frequency of velocity perturbations
  mutable std::vector<int> m_inletMeanInitialized;
                                        ///< If true, inlet mean is good to go!
  mutable std::vector<int> m_levelDefined;
                                        ///< If true, we can define curr level
  mutable std::vector<int> m_timeCounter;
                                        ///< Used for both time-dependent BCs
                                        ///< and time-averaging the solution
  mutable std::vector<Real> m_timeAfterInit;
                                        ///< Used for time-averaging solution
  mutable std::vector<int> m_stageCounter;
                                        ///< Used for both time-dependent BCs
                                        ///< and time-averaging the solution
  mutable std::vector<Real> m_currDt;
                                        ///< Specifically for computing the
                                        ///< backwards time-derivative for
                                        ///< improving the boundary conditions
  mutable std::vector<Real> m_prevDt;
                                        ///< Specifically for computing the
                                        ///< backwards time-derivative for
                                        ///< improving the boundary conditions
  int m_numLevels;                      ///< Maximum number of AMR levels.
                                        ///< Necessary for allocating recycling
                                        ///< plane data structures.
  int m_useAvgYVel;                     ///< Recycles mean wall-normal velocity
  int m_inflowMethod;                   ///< Inflow BC method
  int m_outflowMethod;                  ///< Outflow BC method
  int m_useWallModel;                   ///< 1 = using wall-model
  // Note: gcc ~6.3 won't tolerate emplacement of an FArrayBox so use a pointer
  // (using FArrayBox directly as follows does work on gcc ~8.3)
  // auto ins = m_hackRelax.insert(
  //   { a_boundaryFaceBox, FArrayBox(a_boundaryFaceBox, 1) });
  // This would make destruction very ugly for AMR...
  mutable std::unordered_map<Box, FArrayBox*, CH_Hash::google_CityHash<Box>>
  m_hackRelax;
};

#endif /* ! defined _CNSIBCRECIRCULATINGINLETTFP_H_ */
