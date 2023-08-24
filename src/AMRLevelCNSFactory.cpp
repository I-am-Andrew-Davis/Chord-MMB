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
 * \file AMRLevelCNSFactory.cpp
 *
 * \brief Member functions for AMRLevelCNSFactory
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "MultiBlockCoordSys.H"
#include "TimeInterpolatorRK4.H"
#include "MOLUtilities.H"

//----- Internal -----//

#include "AMRLevelCNSFactory.H"
#include "AMRLevelCNS.H"
#include "TagLevel.H"
#include "TagLevelFactory.H"


/*******************************************************************************
 *
 * Class AMRLevelCNSFactory: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default constructor
/**
 *//*-----------------------------------------------------------------*/

AMRLevelCNSFactory::AMRLevelCNSFactory()
:
  m_MBCS(),
  m_flagDefined(0)
{
}

/*--------------------------------------------------------------------*/
//  Weak construction
/** CRDparam::numAMRLevel() must be defined
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNSFactory::define()
{
  m_MBCS.assign(CRDparam::numAMRLevel(),
                RefCountedPtr<MultiBlockCoordSys>(nullptr));
  m_flagDefined |= FlagDefinedCls;
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

AMRLevelCNSFactory::~AMRLevelCNSFactory()
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Create a new AMRLevelCNS object
/** Deep in the library and difficult to change so weak construction
 *  is largely use to define the AMRLevel class within the application
 *  code
 *//*-----------------------------------------------------------------*/

AMRLevel*
AMRLevelCNSFactory::new_amrlevel() const
{
  return new AMRLevelCNS(this);
}

/*--------------------------------------------------------------------*/
//  Returns objects required for defining an AMRLevelCNS object
/** \param[in]  a_problemDomain
 *                      Problem domain defined on this level.  For
 *                      multiblock, it is the aggregate of all
 *                      domains.
 *  \param[in]  a_level Index of this level with coarsest defined
 *                      as 0.
 *  \param[in]  a_fnRefRatio
 *                      Refinement ratio between this level and the
 *                      next finer level.
 *  \param[in]  a_dx    Mesh spacing
 *  \param[out] a_MBCoordSys
 *                      A multiblock coordinate system created for the
 *                      level
 *  \param[out] a_tagLevel
 *                      A tagging strategy created for the level
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNSFactory::define_amrlevel(
  const int                          a_level,
  const int                          a_fnRefRatio,
  RealVect&                          a_dx,
  RefCountedPtr<MultiBlockCoordSys>& a_MBCoordSys,
  RefCountedPtr<TagLevel>&           a_tagLevel) const
{
  // Size of computational mesh spacing on this level.
  a_dx = getLevelDx(a_level);

  CH_assert(m_flagDefined == FlagDefinedAll);
  // !!!FIXME!!!  There should only be one coordinate system that exists
  a_MBCoordSys = new_MultiBlockCoordSys(a_level);
  a_tagLevel = new_TagLevel(a_level, a_fnRefRatio);
}

/*--------------------------------------------------------------------*/
//  Create a multiblock coordinate system object
/** The default implementation uses a MultiBlockCoordSysFactory object
 *  \param[in]  a_level Index of this level with coarsest defined
 *                      as 0.
 *  \return             RCP to multiblock coordinate system
 *//*-----------------------------------------------------------------*/

RefCountedPtr<MultiBlockCoordSys>
AMRLevelCNSFactory::new_MultiBlockCoordSys(const int a_level) const
{
  CH_assert(m_flagDefined & FlagDefinedCls);
  CH_assert(a_level < CRDparam::numAMRLevel());
  if (m_MBCS[a_level].getRefToThePointer() == nullptr)
    {
      // Size of computational mesh spacing on this level.
      const auto dx = getLevelDx(a_level);
  
      // Size of problem domain on this level
      const ProblemDomain problemDomain =
        refine(m_baseProblemDomain, CRDparam::g_refFromBase[a_level]);
      CH_assert(!problemDomain.isEmpty());

      m_MBCS[a_level] = RefCountedPtr<MultiBlockCoordSys>(
        m_multiBlockCoordSysFactory->getCoordSys(problemDomain, dx));
      m_MBCS[a_level]->setVerbosity(CRDparam::g_verbosity);
    }
  return m_MBCS[a_level];
}

/*--------------------------------------------------------------------*/
//  Create a new cell tagging strategy
/** \param[in]  a_level Index of this level with coarsest defined
 *                      as 0.
 *  \param[in]  a_fnRefRatio
 *                      Refinement ratio between this level and the
 *                      next finer level.
 *  \return             RCP to TagLevel object
 *//*-----------------------------------------------------------------*/

RefCountedPtr<TagLevel>
AMRLevelCNSFactory::new_TagLevel(const int a_level,
                                 const int a_fnRefRatio) const
{
  return m_tagLevelFactory->new_tagLevel(a_level, a_fnRefRatio);
}

/*--------------------------------------------------------------------*/
//  Solve for dx at a level
/** \param[in]  a_level Index of this level with coarsest defined
 *                      as 0.
 *  \return             dx corresponding to the given level
 *//*-----------------------------------------------------------------*/

RealVect
AMRLevelCNSFactory::getLevelDx(const int a_level) const
{
  if (CRDparam::g_coordSysType == CRDparam::CoordSysMultiBlockExternal ||
      CRDparam::g_coordSysType == CRDparam::CoordSysExtrudedMultiBlock)
    {
      // The grid has not yet been read by the CS class so we do not know the
      // domain size.  Dx is assumed to be 1 on the base grid.
      return RealVect(RealVect_unit/CRDparam::g_refFromBase[a_level]);
    }
  else
    {
      // Size of computational mesh spacing on this level.
      IntVect levelDomainSize =
        CRDparam::g_domainBaseSize*CRDparam::g_refFromBase[a_level];
      CH_assert(levelDomainSize > IntVect::Zero);
      CH_assert(CRDparam::g_domainLength > RealVect::Zero);
      return RealVect(CRDparam::g_domainLength/levelDomainSize);
    }
}

/*--------------------------------------------------------------------*/
//  Set or reset the problem domain
/** \param[in]  a_problemDomain
 *                      For multiblock, it is the aggregate of all
 *                      domains.
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNSFactory::setBaseProblemDomain(
  const ProblemDomain& a_baseProblemDomain)
{
  m_flagDefined |= FlagDefinedPD;
  m_baseProblemDomain = a_baseProblemDomain;
}

/*--------------------------------------------------------------------*/
//  Set or reset the coordinate system factory
/** \param[in]  a_multiBlockCoordSysFactory
 *                      Allocate with new and do not delete.
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNSFactory::setCoordinateSystemFactory(
  const MultiBlockCoordSysFactory *const a_multiBlockCoordSysFactory)
{
  m_flagDefined |= FlagDefinedMBCS;
  m_multiBlockCoordSysFactory =
    RefCountedPtr<const MultiBlockCoordSysFactory>(a_multiBlockCoordSysFactory);
}

/*--------------------------------------------------------------------*/
//  Set or reset the tag level factory
/** \param[in]  a_tagLevelFactory
 *                      Allocate with new and do not delete.
 *//*-----------------------------------------------------------------*/

void
AMRLevelCNSFactory::setTagLevelFactory(
  const TagLevelFactory *const a_tagLevelFactory)
{
  m_flagDefined |= FlagDefinedTL;
  m_tagLevelFactory = RefCountedPtr<const TagLevelFactory>(a_tagLevelFactory);
}
