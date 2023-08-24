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
 * \file CNSIBCRMI.H
 *
 * \brief Initial and boundary conditions for Richtmyer-Meshkov instability case
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCRMI_H_
#define _CNSIBCRMI_H_

//----- Internal -----//

#include "CNSIBCCombustionReference.H"


/*******************************************************************************
 */
///  Initialization and BC for Richtmyer-Meshkov instability case
/**
 *//*+*************************************************************************/

class CNSIBCRMI : public CNSIBCCombustionReference
{
public:

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCRMI();

  /// Destructor
  virtual ~CNSIBCRMI();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCRMI(const CNSIBCRMI&);

//--Assignment not permitted

  CNSIBCRMI &operator=(const CNSIBCRMI&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Return a name describing the IBC
  virtual const char *const IBCName() const;

  /// Write any information related to the IBC to output
  virtual void writeIBCInfo() const;

  /// Initialize a level
  virtual void initialize(LevelData<FArrayBox>&      a_U,
                          LevelGridMetrics&          a_gridMetrics,
                          const LayoutData<FluxBox>& a_unitNormals,
                          const Real                 a_time,
                          const int                  a_level) const;

//--Patch operators

  /// Does an exact solution exist?
  virtual bool haveExactSol() const;
  
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
    const int                     a_dir,
    const Side::LoHiSide&         a_side,
    const Box&                    a_disjointBox,
    LevelGridMetrics&             a_gridMetrics,
    const Real                    a_time,
    const int                     a_level,
    const CRDparam::DomainBCType& a_domT) const;

/*==============================================================================
 * Private member functions
 *============================================================================*/

private:

  /// Read any information related to the IBC from input
  void readBCInfo();


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  Real m_xs;                          ///< Initial location of the shock
  Real m_lambda;                      ///< Wavelength of the interface
  Real m_rhoo2;                       ///< Density of O2
  Real m_rhoh2;                       ///< Density of H2
  Real m_Patm;                        ///< Atmospheric pres in regions 2 and 3
  Real m_presRatio;                   ///< Pressure ratio across shock
  Real m_rhoRatio;                    ///< Density ratio across shock
  RealVect m_U0;                      ///< Interface velocity
  int m_O2comp;                       ///< O2 comp number
  int m_H2comp;                       ///< H2 comp number
  int m_loBC;                         ///< Low BC 0 - extrapolated
                                      ///< 1 - outflow, 2 - dirichlet
  int m_bcOrder;                      ///< BC order, either 1 or 4
  Real m_h0;                          ///< Initial amplitude of the interface
  Real m_Delta;                       ///< Interface thickness
};

#endif /* ! defined _CNSIBCRMI_H_ */
