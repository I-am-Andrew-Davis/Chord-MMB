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
 * \file CNSIBCVortex.H
 *
 * \brief Initial and boundary conditions for species mixing test case
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCVORTEX_H_
#define _CNSIBCVORTEX_H_

//----- Internal -----//

#include "CNSIBCReferenceCubeBC.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for species mixing test case.
/**
 *//*+*************************************************************************/

class CNSIBCVortex : public CNSIBCReferenceCubeBC
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCVortex();

  /// Destructor
  virtual ~CNSIBCVortex();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCVortex(const CNSIBCVortex&);

//--Assignment not permitted

  CNSIBCVortex &operator=(const CNSIBCVortex&);


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

  /// Set inlet profiles for velocity, temperature, species
  virtual void setCNSCBCProfiles(
    FArrayBox&                    a_BCProfile,
    const Box&                    a_boundaryBox,
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
  Real m_beta;                        ///< Vortex strength
  RealVect m_U;                       ///< Freestream velocity
  Real m_Pinf;                        ///< Stagnation pressure
  Real m_Tinf;                        ///< Stagnation temperature
  std::vector<Real> m_massFractions;  ///< Initial mass fractions in domain
};

#endif /* ! defined _CNSIBCVORTEX_H_ */
