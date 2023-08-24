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
 * \file CNSIBCFlame.H
 *
 * \brief Initial and boundary conditions for a general flame problem
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCFLAME_H_
#define _CNSIBCFLAME_H_

//----- Internal -----//

#include "CNSIBCGeneralizedSingleBlock.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for species mixing test case.
/**
 *                                 Periodic
 *             +===============================================
 *             |               !         !                     )
 *             |   rho1        !  rho2   !       rho1          )
 *             |   p1          !  p2     !       p1            )
 * Periodic    |   T1          !  T2     !       T1            )  Periodic
 *             |   u1          !  u2     !       u1            )
 *             |               !         !                     )
 *             +===============================================
 *                             x3        x2
 *
 *//*+*************************************************************************/

class CNSIBCFlame : public CNSIBCGeneralizedSingleBlock
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCFlame();

  /// Destructor
  virtual ~CNSIBCFlame();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCFlame(const CNSIBCFlame&);

//--Assignment not permitted

  CNSIBCFlame &operator=(const CNSIBCFlame&);


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
  
/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Set the imposed BC for primitive variables
  virtual void setImposedBCprimState(
    FArrayBox&           a_Wface,
    const Box&           a_boundaryFaceBox,
    const FArrayBox&     a_Wcell,
    const FArrayBox&     a_unitNormalBasisFab,
    const BoundaryIndex& a_bcIdx,
    const Box&           a_disjointBox,
    LevelGridMetrics&    a_gridMetrics,
    const Real           a_time,
    const int            a_level,
    const BCInfo&        a_domT) const;

  /// Set inlet profiles for velocity, temperature, species
  virtual void setCNSCBCProfiles(
    FArrayBox&           a_BCProfile,
    const Box&           a_boundaryBox,
    const FArrayBox&     a_Wcell,
    const FArrayBox&     a_unitNormalBasisFab,
    const BoundaryIndex& a_bcIdx,
    const Box&           a_disjointBox,
    LevelGridMetrics&    a_gridMetrics,
    const Real           a_time,
    const int            a_level,
    const BCInfo&        a_domT) const;

/*==============================================================================
 * Private member functions
 *============================================================================*/

private:

  /// Read any information related to the IBC from input
  void readBCInfo();

  /// Add inlet perturbations
  void addInletPerturbations(FArrayBox&            a_WFace,
                             const Box&            a_box,
                             const BoundaryIndex&  a_bcIdx,
                             const Box&            a_disjointBox,
                             LevelGridMetrics&     a_gridMetrics,
                             const Real            a_time,
                             const int             a_level) const;

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  int m_flameDir;                     ///< Direction of flame
  int m_varDir;                       ///< Direction perturbations occur
  Real m_pertMag;                     ///< Magnitude of perturbations
  Real m_pertPeriod;                  ///< Period of perturbations
  Real m_pertT;                       ///< Temperature in perturbation
  Real m_pertP;                       ///< Pressure in perturbation
  RealVect m_pertStartPoint;          ///< Start point for perturb initial
  RealVect m_pertEndPoint;            ///< End point for perturb initial
  Real m_pertRadius;                  ///< Radius of perturb
  Real m_pertCr;                      ///< Sharpness of pert interface
  Real m_pertStartTime;               ///< Start time for perturbing inlet
  Real m_pertEndTime;                 ///< End time for perturbing inlet
  Real m_Cr;                          ///< Determines sharpness of interfaces
  Real m_rho2;                        ///< Density in region 2
  Real m_p2;
  Real m_t2;
  RealVect m_loEnd2;                  ///< Lower end vector for region 2
  RealVect m_hiEnd2;                  ///< Upper end vector for region 2
  std::vector<Real> m_region2MassFraction;
                                      ///< Mass fractions in region 2
  int m_cirFlame;                     ///< 0 - rectangle, uses loEnd hiEnd
                                      ///< 1 - circular flame
};

#endif /* ! defined _CNSIBCFLAME_H_ */
