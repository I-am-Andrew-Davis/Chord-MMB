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
 * \file CNSIBCShockBubble.H
 *
 * \brief Initial and boundary conditions for shock bubble case
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCSHOCKBUBBLE_H_
#define _CNSIBCSHOCKBUBBLE_H_

//----- Internal -----//

#include "CNSIBCCombustionReference.H"


/*******************************************************************************
 */
///  Initialization and BC for shock bubble case
/**
 *//*+*************************************************************************/

class CNSIBCShockBubble : public CNSIBCCombustionReference
{
public:

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCShockBubble();

  /// Destructor
  virtual ~CNSIBCShockBubble();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCShockBubble(const CNSIBCShockBubble&);

//--Assignment not permitted

  CNSIBCShockBubble &operator=(const CNSIBCShockBubble&);


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
  
  Real m_C;                           ///< Constant in initialization equation
  Real m_radius;                      ///< Radius of bubble
  Real m_shockLoc;                    ///< Location of shock
  Real m_t0;                          ///< Temperature upstream of shock
  Real m_t1;                          ///< Temperature downstream of shock
  Real m_p0;                          ///< Pressure upstream of shock
  Real m_p1;                          ///< Pressure downstream of shock
  RealVect m_U0;                      ///< Velocity upstream of shock
  RealVect m_U1;                      ///< Velocity downstream of shock
  RealVect m_Center;                  ///< Center location of bubble
  int m_O2comp;                       ///< O2 comp number
  int m_H2comp;                       ///< H2 comp number
  int m_N2comp;                       ///< N2 comp number
  int m_loBC;                         ///< Low and high BC's
  int m_hiBC;                         ///< 0 - extrapolated, 1 - outflow
  int m_bcOrder;                      ///< BC order, either 1 or 4
};

#endif /* ! defined _CNSIBCREACTIONADVECTION_H_ */
