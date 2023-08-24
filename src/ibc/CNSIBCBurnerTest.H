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
 * \file CNSIBCBurnerTest.H
 *
 * \brief Initial and boundary conditions for species mixing test case
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCBURNERTEST_H_
#define _CNSIBCBURNERTEST_H_

//----- Internal -----//

#include "CNSIBCGeneralizedSingleBlock.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for burner test case
/**
 *                    Outlet
 *               -+-------------+-
 *                |             |
 *                |             |
 *         Slip   |             | Slip
 *         Wall   |             | Wall
 *                |             |
 *               _|             |____________________
 *                |             |   Top of the viewing window: m_boxTop
 *                |             |
 *     Air jets ->               <- Air jets
 *              ->               <-
 *                |             |
 *                |             |
 *               -|=============|-------------------- Interface
 *                |             |    Bottom of the viewing window: m_boxBottom
 *                +-------------+---------------------
 *                 ^^^^^^^^^^^^^
 *                     Inlet
 *//*+*************************************************************************/

class CNSIBCBurnerTest : public CNSIBCGeneralizedSingleBlock
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCBurnerTest();

  /// Destructor
  virtual ~CNSIBCBurnerTest();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCBurnerTest(const CNSIBCBurnerTest&);

//--Assignment not permitted

  CNSIBCBurnerTest &operator=(const CNSIBCBurnerTest&);


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

//--Patch operators

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
                             const Real           a_globalKE,
                             const Real           a_globalHelicity) const;
  
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

  /// Set boxes for inflow or outflow conditions in mixed boundaries
  virtual int setMixedBC(const Box&           a_boundaryFaceBox,
                         const BoundaryIndex& a_bcIdx,
                         const Box&           a_disjointBox,
                         LevelGridMetrics&    a_gridMetrics,
                         const Real           a_time,
                         const int            a_level,
                         Vector<Box>&         a_boxVect,
                         Vector<BCInfo>&      a_domainType) const override;

/*==============================================================================
 * Private member functions
 *============================================================================*/

private:

  /// Read any information related to the IBC from input
  void readBCInfo();

  /// Set jet boundary conditions
  void setJetBCVals(const Box&               a_box,
                    FArrayBox&               a_Wface,
                    const FArrayBox&         a_XFab,
                    const Real&              a_inletT,
                    const Real&              a_inletGamma,
                    const Real&              a_inletVel,
                    const Real&              a_Rgas,
                    const std::vector<Real>& a_inletCn,
                    const Real&              a_loLoc,
                    const Real&              a_hiLoc,
                    const int&               a_tanDir,
                    const int&               a_normDir) const;

  /// Need to find mapped coordinates
  /// FIXME: Not the best way to do this
  Real coordOut(const Real a_Xi);

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  int m_flowDir;                      ///< Direction of flow
  int m_jetDir;                       ///< Direction normal to the jets
  RealVect m_jetLoVel;                ///< Lower jet velocity
  RealVect m_jetHiVel;                ///< Upper jet velocity
  RealVect m_hotSpotCenter;           ///< Location of hot spot center
  Real m_hotSpotTemp;                 ///< Temperature of the hot spot
  Real m_hotSpotStartTime;            ///< Start time for hot spot to occur
  Real m_hotSpotRunTime;              ///< Time for the hot spot to exist
  bool m_energyIgnition;              ///< If true, energy will be added as
                                      ///< a source term
                                      ///< If false, ignition occurs through
                                      ///< lowering of the activation energy
  // FIXME: hotSpotIgn is a poor way to reduce time steps prior to applying
  // the source term
  mutable int m_hotSpotIgn;           ///< For setting dt before hot spot
  Real m_gravForce;                   ///< Gravity force
  Real m_jetLoc;                      ///< Jet location
  Real m_viewerTop;                   ///< Top of viewer window
  Real m_viewerBottom;                ///< Bottom of viewing window
  Real m_jetLength;                   ///< Length of jet
  std::vector<Real> m_jetLoLoc;       ///< Low end location of jet on each level
  std::vector<Real> m_jetHiLoc;       ///< Hi end location of jet on each level
  std::vector<Box> m_belowJets;       ///< Box of the cells below the jets
  std::vector<Box> m_atJets;          ///< Box of the cells at the jets
  std::vector<Box> m_aboveJets;       ///< Box of the cells above the jets
  std::vector<Box> m_slipWalls;       ///< Box of the top slip walls
  std::vector<Box> m_hotSpot;         ///< Box of the hot spot
  CRDparam::DomainBCType m_jetLoBC;   ///< BC type for lower jet
  CRDparam::DomainBCType m_jetHiBC;   ///< BC type for upper jet
  RealVect m_mappedShift;             ///< Shift operation for log stretch
                                      ///< mapping
};

#endif /* ! defined _CNSIBCBURNERTEST_H_ */
