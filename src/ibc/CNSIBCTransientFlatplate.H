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
 * \file CNSIBCTransientFlatplate.H
 *
 * \brief Initial and boundary conditions for transient boundary layer growth
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCTRANSIENTFLATPLATE_H_
#define _CNSIBCTRANSIENTFLATPLATE_H_

//----- Internal -----//

#include "CNSIBCReferenceCubeBC.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for transient boundary layer growth
///  !!! FIXME !!! this should derive from the generalized case
/**
 *//*+*************************************************************************/

class CNSIBCTransientFlatplate : public CNSIBCReferenceCubeBC
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
 CNSIBCTransientFlatplate();

  /// Destructor
  virtual ~CNSIBCTransientFlatplate();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCTransientFlatplate(const CNSIBCTransientFlatplate&);

//--Assignment not permitted

  CNSIBCTransientFlatplate &operator=(const CNSIBCTransientFlatplate&);


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

//--Level operators

  /// Initialize a level
  virtual void initialize(LevelData<FArrayBox>&      a_U,
                          LevelGridMetrics&          a_gridMetrics,
                          const LayoutData<FluxBox>& a_unitNormals,
                          const Real                 a_time,
                          const int                  a_level) const;
//--Patch operators

  /// Does an exact solution exist?
  virtual bool haveExactSol() const;

  /// Compute the exact solution state in the cells
  virtual int exactSol(FArrayBox&              a_Ux,
                       const Box&              a_box,
                       const Box&              a_disjointBox,
                       const LevelGridMetrics& a_gridMetrics,
                       const FluxBox&          a_unitNormals,
                       const DataIndex&        a_didx,
                       const Real              a_time,
                       const int               a_level) const;


/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Set the primitive state at wall BC
  virtual void setWallBCprimState(
    FArrayBox&                    a_Wface,
    const Box&                    a_boundaryFaceBox,
    const FArrayBox&              a_Wcell,
    const FArrayBox&              a_boundarySlipVelocity,
    const int                     a_dir,
    const Side::LoHiSide&         a_side,
    const Box&                    a_disjointBox,
    LevelGridMetrics&             a_gridMetrics,
    const Real                    a_time,
    const int                     a_level,
    const CRDparam::DomainBCType& a_domT) const;

  /// Set the boxes for the flat plate
  virtual int setMixedBC(const Box&            a_boundaryFaceBox,
                         const int             a_dir,
                         const Side::LoHiSide& a_side,
                         const Box&            a_disjointBox,
                         LevelGridMetrics&     a_gridMetrics,
                         const Real            a_time,
                         const int             a_level,
                         Vector<Box>&          a_boxVect,
                         Vector<BCInfo>&       a_domainType) const;
  

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

  RealVect m_farfieldVelocity;        ///< Velocity of the wall
  Real m_plateLength;                 ///< Length of the no-slip plate
  Real m_maxBLdisplacement;           ///< The maximum displacement of the
                                      ///< boundary layer at the end of the
                                      ///< plate
  int m_cellPlateLength;              ///< Plate length in terms of cells.
  int m_wallNormalDir;                ///< Direction normal to the wall
  int m_flowDir;                      ///< Direction of flow
  Box m_plateBox;                     ///< A box along the wall defining where
                                      ///< the flow is no-slip.  This is a cell
                                      ///< box extending 1 cell into the
                                      ///< interior in m_wallNormalDir
};

#endif /* ! defined _CNSIBCTRANSIENTFLATPLATE_H_ */
