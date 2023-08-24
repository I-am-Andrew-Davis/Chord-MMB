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
 * \file SA.H
 *
 * \brief Class for Large Eddy Simulation
 *
 *//*+*************************************************************************/

#ifndef _SA_H_
#define _SA_H_

//----- Chombo Library -----//

#include "LGintegrator.H"
#include "LevelGridMetrics.H"

//----- Internal -----//

#include "TurbModeling.H"

//--Forward declarations

/*******************************************************************************
 */
///  SA derived from turbulence modeling base class
/**
 *   This class is an implementation of the Spalart-Allmaras
 *
 *//*+*************************************************************************/

class SA : public TurbModeling
{

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  SA();

  /// Destructor
  ~SA();

/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  SA(const SA&);

//--Assignment not permitted

  SA &operator=(const SA&);

/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

//--Variables

  /// Return the number of turbulent variables
  virtual int numTurbVar() const;

//--Misc

  /// Return the name the variable names, index a_iComp only relates to
  /// turbulence variables, not any other variables
  virtual const char *const turbStateName(const int a_iComp) const;

  /// Output turbulence model info
  virtual void writeTurbModelInfo() const;

  /// Output turbulence model name
  virtual const char *const modelName() const;

//--Algorithms (on a patch)

  /// Calculate dynamic eddy viscosity and turbulent thermal conductivity
  virtual void calcCoeffKappatMut(const Box&              a_box,
                                  FArrayBox&              a_mutFab,
                                  FArrayBox&              a_kappatFab,
                                  const FArrayBox&        a_muFab,
                                  const FArrayBox&        a_kappaFab,
                                  const FArrayBox&        a_WfacePntFab,
                                  const FArrayBox&        a_strainRateFab,
                                  const LevelGridMetrics& a_gridMetrics,
                                  const int               a_dir) const;

  /// Calculate the turbulent convective fluxes
  virtual void addTurbConvectiveFlux(FArrayBox&       a_flux,
                                     const FArrayBox& a_WFace,
                                     const int&       a_dir,
                                     const Box&       a_box) const;

  /// Calculate the turbulent diffusive fluxes
  virtual void calcTurbDiffFlux(const Box&       a_box,
                                FArrayBox&       a_NGradNufacePntFab,
                                const FArrayBox& a_muFab,
                                const FArrayBox& a_kappaFab,
                                const FArrayBox& a_WfacePntFab,
                                const int        a_dir,
                                const int        a_wTurbComp,
                                const int        a_tComp) const;

  /// Calculate the turbulent source terms
  virtual void calcTurbSourceTerms(
    const Box&              a_box,
    const Box&              a_disjointBox,
    FArrayBox&              a_turbSourcePntFab,
    const FArrayBox&        a_GradWcellPntFab,
    const FArrayBox&        a_WcellPntFab,
    const DataIndex&        a_dataIndx,
    const LevelGridMetrics& a_gridMetrics) const;

  /// Initialize turbulent variables
  virtual void turbInitialize(FArrayBox&              a_U,
                              const FArrayBox&        a_W,
                              const LevelGridMetrics& a_gridMetrics,
                              const FluxBox&          a_unitNormals,
                              const DataIndex&        a_didx,
                              const Box&              a_disjointBox,
                              const Box&              a_box) const;

  /// Find primitive turbulent variables from conservative variables
  virtual void turbConsToPrim(FArrayBox&       a_W,
                              const FArrayBox& a_U,
                              const Box&       a_box) const;

  /// Set turbulence values at boundaries
  virtual void setTurbulentBC(
    const Box&                    a_boundaryFaceBox,
    FArrayBox&                    a_Wface,
    const FArrayBox&              a_Wcell,
    const FArrayBox&              a_unitNormalBasisFab,
    const int                     a_dir,
    const Side::LoHiSide&         a_side,
    const LevelGridMetrics&       a_gridMetrics,
    const Real                    a_time,
    const int                     a_level,
    const CRDparam::DomainBCType& a_domT) const;

  /// Apply wall model
  virtual void applyWallModel(
    const Box&                    a_boundaryFaceBox,
    const Box&                    a_boundaryFaceGhostBox,
    const Box&                    a_disjointBox,
    FArrayBox&                    a_WfaceBdryFab,
    FArrayBox&                    a_WfaceAvgDirFab,
    FArrayBox&                    a_WcellAvgFab,
    const FArrayBox&              a_unitNormalBasisFab,
    const CRDparam::DomainBCType& a_bcType,
    const Box&                    a_bcBox,
    const LevelGridMetrics&       a_gridMetrics,
    const DataIndex&              a_dataIndx,
    const int                     a_dir,
    const Side::LoHiSide&         a_side,
    const Real                    a_time,
    const int                     a_level) const;

/*==============================================================================
 * Data members
 *============================================================================*/

protected:
  Real m_Cb1;                         ///< Constant Cb1
  Real m_Cb2;                         ///< Constant Cb2
  Real m_Cw1;                         ///< Constant Cw1
  Real m_sigma;                       ///< Constant sigma
  Real m_k;                           ///< Constant k
  Real m_Cw2;                         ///< Constant Cw2
  Real m_Cw3;                         ///< Constant Cw3
  Real m_Cv1;                         ///< Constant Cv1
  Real m_Ct3;                         ///< Constant Ct3
  Real m_Ct4;                         ///< Constant Ct4
  Real m_Cprod;                       ///< From Fluent model
  Real m_Prt;                         ///< Turbulent Prandtl number
  Real m_initTurb;                    ///< Initial turbulent intensity
  Real m_inletTurb;                   ///< Inlet turbulent intensity
  Real m_refLen;                      ///< Turbulent length scale
  Real m_yPlusLam;                    ///< Laminar y+
  Real m_E;                           ///< Constant for wall model
};

/*******************************************************************************
 *
 * Class SA: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Number of turbulence variables
/** 
 *//*-----------------------------------------------------------------*/

inline int
SA::numTurbVar() const
{
  return 1;
}

/*--------------------------------------------------------------------*/
//  Return the name of the turbulence variables, default is nothing
/** \return             Name of the a_iComp variable
 *//*-----------------------------------------------------------------*/

inline const char *const
SA::turbStateName(const int a_iComp) const
{
  CH_assert(a_iComp == 0);
  return "mu_t";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the turbulence model
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

inline void
SA::writeTurbModelInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Cb1\n" << m_Cb1 << CRD::var;
  CRD::msg << "Cb2\n" << m_Cb2 << CRD::var;
  CRD::msg << "Cw1\n" << m_Cw1 << CRD::var;
  CRD::msg << "Sigma\n" << m_sigma << CRD::var;
  CRD::msg << "k\n" << m_k << CRD::var;
  CRD::msg << "Cw2\n" << m_Cw2 << CRD::var;
  CRD::msg << "Cw3\n" << m_Cw3 << CRD::var;
  CRD::msg << "Cprod\n" << m_Cprod << CRD::var;
  CRD::msg << "Turbulent Prandtl number\n" << m_Prt << CRD::var;
  CRD::msg << "Initial turbulent intensity\n" << m_initTurb << CRD::var;
  CRD::msg << "Inlet turbulent intensity\n" << m_inletTurb << CRD::var;
  CRD::msg << "Turbulent reference length\n" << m_refLen << CRD::var;
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Write the name of the turbulence model
/**
 *//*-----------------------------------------------------------------*/

inline const char *const
SA::modelName() const
{
  return "SA model";
}

#endif  /* ! defined _SA_H_ */
