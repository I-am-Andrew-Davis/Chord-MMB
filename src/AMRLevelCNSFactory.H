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
 * \file AMRLevelCNSFactory.H
 *
 * \brief Factory object for AMRLevelCNS and related strategies
 *
 *//*+*************************************************************************/

#ifndef _AMRLEVELCNSFACTORY_H_
#define _AMRLEVELCNSFACTORY_H_

//----- Chombo Library -----//

#include "AMRLevelFactory.H"
#include "UsingNamespace.H"
#include "NamespaceVar.H"

//----- Internal -----//

#include "CRDparam.H"

//-- Forward declarations from Chombo

#include "NamespaceHeader.H"
class MultiBlockCoordSys;
class MultiBlockCoordSysFactory;
#include "NamespaceFooter.H"

//-- Forward declarations

class TagLevel;
class TagLevelFactory;


/*******************************************************************************
 */
///  AMRLevelCNSFactory
/**
 *   This class derives from AMRLevelFactory.  It implements the single pure
 *   virtual function in AMRLevelFactory.  This class is also an aggregate
 *   of other factories used to construct objects required by the AMRLevelCNS
 *   class.
 *
 *//*+*************************************************************************/

//** COMMENT - Maybe a lot of this stuff should be block specific.  I.e.,
//**           tagging per block.

class AMRLevelCNSFactory : public CH_XDIR::AMRLevelFactory
{

//--Internal types

  enum FlagDefined : unsigned
  {
    FlagDefinedCls  = 1<<0,           ///< This class has been defined
    FlagDefinedPD   = 1<<1,           ///< Problem domain
    FlagDefinedMBCS = 1<<2,           ///< Coordinate system
    FlagDefinedTL   = 1<<3,           ///< Tag level factory
    FlagDefinedAll  = (1<<4) - 1
  };

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  AMRLevelCNSFactory();

  //**FIXME - what should this do? create preset MBCS and TL?
  // /// Construct with a known type
  // AMRLevelCNSFactory(CRDparam::ProblemType a_problemType);

  /// Weak construction
  void define();

  /// Destructor
  virtual ~AMRLevelCNSFactory();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  AMRLevelCNSFactory(const AMRLevelCNSFactory&);

//--Assignment not permitted

  AMRLevelCNSFactory &operator=(const AMRLevelCNSFactory&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Create a new AMRLevelCNS object
  virtual AMRLevel* new_amrlevel() const override;

  /// Returns objects required for defining an AMRLevelCNS object
  virtual void define_amrlevel(
    const int                          a_level,
    const int                          a_fnRefRatio,
    RealVect&                          a_dx,
    RefCountedPtr<MultiBlockCoordSys>& a_MBCoordSys,
    RefCountedPtr<TagLevel>&           a_tagLevel) const;

  /// Create a multiblock coordinate system object
  RefCountedPtr<MultiBlockCoordSys> new_MultiBlockCoordSys(
    const int       a_level) const;

  /// Create a new cell tagging strategy
  RefCountedPtr<TagLevel> new_TagLevel(
    const int a_level,
    const int a_fnRefRatio) const;

  /// define dx on a given level
  RealVect getLevelDx(const int a_level) const;
  
  /// Set or reset the problem domain
  void setBaseProblemDomain(const ProblemDomain& a_baseProblemDomain);

  /// Set or reset the coordinate system factory
  void setCoordinateSystemFactory(
    const MultiBlockCoordSysFactory *const a_multiBlockCoordSysFactory);

  /// Set or reset the tag level factory
  void setTagLevelFactory(
    const TagLevelFactory *const a_tagLevelFactory);


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  ProblemDomain m_baseProblemDomain;  ///< Problem domain on the base level.
                                      ///< For multiblock, this can be some
                                      ///< nominal set of cell numbers in each
                                      ///< direction that help define the grid
                                      ///< resolution.  How exactly they are
                                      ///< used depends on the coordinate
                                      ///< system.
  RefCountedPtr<const MultiBlockCoordSysFactory> m_multiBlockCoordSysFactory;
                                      ///< Factory for MultiBlockCoordSys
  mutable std::vector<RefCountedPtr<MultiBlockCoordSys>> m_MBCS;
                                      ///< A cache of MultiBlockCoordSys for
                                      ///< each level since these are requested
                                      ///< by many classes.
                                      //**FIXME There should be only one CS,
                                      //**      not one per level.
  RefCountedPtr<const TagLevelFactory> m_tagLevelFactory;
                                      ///< Factory for TagLevel
  unsigned m_flagDefined;             ///< Bit flag of factories that have
                                      ///< been defined.  Both
                                      ///< setCoordinateSystemFactory and
                                      ///< setTagLevelFactory must be called
                                      ///< before define_amrlevel
};

#endif  /* ! defined _AMRLEVELCNSFACTORY_H_ */
