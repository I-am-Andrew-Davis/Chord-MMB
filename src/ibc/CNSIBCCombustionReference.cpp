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
 * \file CNSIBCCombustionReference.cpp
 *
 * \brief Member functions for CNSIBCCombustionReference
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelGridMetrics.H"
#include "SetCentersF_F.H"
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"

//----- Internal -----//

#include "CNSIBCCombustionReference.H"
#include "CNSPhysics.H"
#include "ThermPhysics.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"
#include "TagMethodValue.H"
#include "TagMethodVMS.H"


/*******************************************************************************
 *
 * Class CNSIBCCombustionReference: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor set BC at domain extents
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCCombustionReference::CNSIBCCombustionReference()
  :
  CNSIBC(),
  m_tagVar(1,""),
  m_threshold(1,0.),
  m_loThreshold(1,0.),
  m_maxLevel(1,-1),
  m_tagType(1,0)
{
  ParmParse ppIBC("ibc");
  const int numTags = std::max(ppIBC.countval("tag_variable"),
                               ppIBC.countval("tag_type"));
  m_tagVar.resize(numTags);
  m_tagVar.assign(numTags, "");
  m_threshold.resize(numTags);
  m_threshold.assign(numTags,0.);
  m_maxLevel.resize(numTags);
  m_maxLevel.assign(numTags,-1);
  m_tagType.resize(numTags);
  m_tagType.assign(numTags,0);
  m_loThreshold.resize(numTags);
  m_loThreshold.assign(numTags, 0.);
  for (int tag = 0; tag != numTags; ++tag)
    {
      ppIBC.query("tag_max_level", m_maxLevel[tag], tag);
      if (m_maxLevel[tag] == 0)
        {
          CRD::msg << "'tag_max_level' must be either -1 (meaning all levels)"
                   << " or > 0!" << CRD::error;
        }
      std::string tagTypeStr("gradient");
      ppIBC.query("tag_type", tagTypeStr, tag);
      if (tagTypeStr == "vorticity")
        {
          m_tagType[tag] = 2;
          ppIBC.get("tag_threshold", m_threshold[tag], tag);
          m_loThreshold[tag] = -1.;
        }
      else if (tagTypeStr == "vms_vorticity")
        {
          m_tagType[tag] = 4;
          ppIBC.get("tag_threshold", m_threshold[tag], tag);
          m_loThreshold[tag] = -1.;
        }
      else
        { 
          if (tagTypeStr == "gradient" || tagTypeStr == "grad")
            {
              m_tagType[tag] = 0;
              m_loThreshold[tag] = 0.;
              ppIBC.query("lo_end_threshold", m_loThreshold[tag], tag);
            }
          else if (tagTypeStr == "value" || tagTypeStr == "val")
            {
              m_tagType[tag] = 1;
              ppIBC.get("tag_lo_threshold", m_loThreshold[tag], tag);
            }
          else if (tagTypeStr == "vms" || tagTypeStr == "VMS")
            {
              m_tagType[tag] = 3;
            }
          ppIBC.get("tag_variable", m_tagVar[tag], tag);
          ppIBC.get("tag_threshold", m_threshold[tag], tag);
        }
    }
  ppIBC.query("outflow_smoothing_beta", m_CBCbeta);
  ppIBC.query("outflow_smoothing_sigma", m_CBCsigma);
  ppIBC.query("inflow_smoothing", m_etaMax);
  ppIBC.query("inflow_smoothing_cn", m_etaCN);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCCombustionReference::~CNSIBCCombustionReference()
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Allocations of new physics states.  Customize for each derivative
//  class.
/**
 *//*-----------------------------------------------------------------*/

std::vector<CRDPhysics*>
CNSIBCCombustionReference::allocatePhysics()
{
  std::vector<CRDPhysics*> physics(1, nullptr);
  if (CRDparam::g_physicsModels & CRDparam::PhysicsThermPerf)
    {
      physics[0] = new ThermPhysics;
    }
  else
    {
      physics[0] = new CNSPhysics;
    }
  return physics;
}

/*--------------------------------------------------------------------*/
//  Return a name describing the IBC
/** \return             Name of IBC
 *//*-----------------------------------------------------------------*/

const char *const
CNSIBCCombustionReference::IBCName() const
{
  return "combustion reference";
}

/*--------------------------------------------------------------------*/
//  Set the tagging method if one can be associated with IBC
/**
 *  \param[in]  a_tagBufferSize
 *                      Requested tag buffer size (should be
 *                      respected).
 *
 *  \note
 *  <ul>
 *    <li> Allocate levels and methods with 'new' and do not delete
 *  </ul>
 *//*-----------------------------------------------------------------*/

TagLevelFactory*
CNSIBCCombustionReference::setTagMethod(const int a_tagBufferSize)
{
  TagLevel* tagLevel = new TagLevel;
  setTagMethodLevel(a_tagBufferSize, tagLevel);
  // Return the factory
  return new TagLevelFactory(tagLevel);
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set the tagging method if one can be associated with IBC
/**
 *  \param[in]  a_tagBufferSize
 *                      Requested tag buffer size (should be
 *                      respected).
 *
 *  \note
 *  <ul>
 *    <li> Allocate levels and methods with 'new' and do not delete
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
CNSIBCCombustionReference::setTagMethodLevel(const int        a_tagBufferSize,
                                             TagLevel*        a_tagLevel,
                                             std::vector<Box> a_restrictBox)
{
  const int numSpecies = CRDparam::g_numSpecies;
  const int numTags = m_tagVar.size();
  for (int tag = 0; tag != numTags; ++tag)
    {
      bool primTag = false;
      int tagComp = -1;
      if (m_tagType[tag] == 2)
        {
          tagComp = -2;
          const bool insideThreshold = false;
          a_tagLevel->appendTagMethod(
            new TagMethodValue(tagComp,
                               m_loThreshold[0],
                               m_threshold[0],
                               insideThreshold,
                               TagMethodValue::ValueTypeVorticity,
                               a_restrictBox));
          continue;
        }
      if (m_tagType[tag] == 4)
        {
          tagComp = -2;
          const bool insideThreshold = false;
          a_tagLevel->appendTagMethod(
            new TagMethodVMS(tagComp,
                             m_loThreshold[tag],
                             m_threshold[tag],
                             insideThreshold,
                             true));
          
          continue;
        }
      std::string tagVarSt = m_tagVar[tag];
      if (tagVarSt == "density")
        {
          tagComp = CRDparam::g_CRDPhysics->densityIndex();
        }
      else if (tagVarSt == "pressure")
        {
          tagComp = CRDparam::g_CRDPhysics->pressureIndex();
          primTag = true;
        }
      else if (tagVarSt == "x-velocity")
        {
          tagComp = CRDparam::g_CRDPhysics->velocityInterval().begin();
          primTag = true;
        }
      else if (tagVarSt == "y-velocity")
        {
          tagComp = CRDparam::g_CRDPhysics->velocityInterval().begin() + 1;
          primTag = true;
        }
      else if (tagVarSt == "z-velocity")
        {
          if (SpaceDim < 3)
            {
              CRD::msg << "Input (CombustionReference IBC): 'tag_variable'"
                       << " must exist!" << CRD::error;
            }
          tagComp = CRDparam::g_CRDPhysics->velocityInterval().begin() + 2;
          primTag = true;
        }
      else if (tagVarSt == "x-momentum")
        {
          tagComp = CRDparam::g_CRDPhysics->vectorFluxInterval().begin();
        }
      else if (tagVarSt == "y-momentum")
        {
          tagComp = CRDparam::g_CRDPhysics->vectorFluxInterval().begin() + 1;
        }
      else if (tagVarSt == "z-momentum")
        {
          if (SpaceDim < 3)
            {
              CRD::msg << "Input (CombustionReference IBC): 'tag_variable'"
                       << " must exist!" << CRD::error;
            }
          tagComp = CRDparam::g_CRDPhysics->vectorFluxInterval().begin() + 2;
        }
      else if (tagVarSt == "temperature")
        {
          tagComp = CRDparam::g_CRDPhysics->temperatureIndex();
          primTag = true;
        }
      else if (tagVarSt == "energy-density")
        {
          tagComp = CRDparam::g_CRDPhysics->energyFluxIndex();
        }
      else
        {
          // Check for primitive species mass fractions
          for (int i = 0; i != numSpecies; ++i)
            {
              if (tagVarSt == CRDparam::g_speciesNames[i])
                {
                  tagComp =
                    CRDparam::g_CRDPhysics->speciesPrimInterval().begin() + i;
                  primTag = true;
                }
            }
          // Check for conservative species mass fractions
          for (int i = 0; i != numSpecies; ++i)
            {
              std::string consvar;
              consvar = "rho" + CRDparam::g_speciesNames[i];
              if (tagVarSt == consvar)
                {
                  tagComp =
                    CRDparam::g_CRDPhysics->speciesConsInterval().begin() + i;
                  primTag = false;
                }
            }
        }
      if (tagComp == -1)
        {
          CRD::msg << "Tagging variable or type not found!" << CRD::error;
        }
      if (m_tagType[tag] == 0)
        {
          a_tagLevel->appendTagMethod(new TagMethodGradient(tagComp,
                                                            m_threshold[tag],
                                                            primTag,
                                                            m_loThreshold[tag],
                                                            a_restrictBox));
        }
      else if (m_tagType[tag] == 1)
        {
          const bool insideThreshold = true;
          TagMethodValue::ValueType vt = TagMethodValue::ValueTypeConservative;
          if (primTag)
            {
              vt = TagMethodValue::ValueTypePrimitive;
            }
          a_tagLevel->appendTagMethod(new TagMethodValue(tagComp,
                                                         m_loThreshold[tag],
                                                         m_threshold[tag],
                                                         insideThreshold,
                                                         vt,
                                                         a_restrictBox));
        }
      else if (m_tagType[tag] == 3)
        {
          const bool insideThreshold = false;
          const bool useVorticity = false;
          a_tagLevel->appendTagMethod(new TagMethodVMS(tagComp,
                                                       m_loThreshold[tag],
                                                       m_threshold[tag],
                                                       primTag,
                                                       insideThreshold,
                                                       useVorticity,
                                                       m_maxLevel[tag]));
        }
    }
  a_tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
  a_tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
}

/*--------------------------------------------------------------------*/
//  Display tagging information, this should be called in the writeIBCInfo()
//  function for any inhereted class. Also displays CNSCBC info when necessary
/** 
 *//*-----------------------------------------------------------------*/

void
CNSIBCCombustionReference::writeTagInfo() const
{
  // Check for CNSCBCOutflow
  bool CBCTest = false;
  BoundaryIndex bcIdx;
  bcIdx.m_block = 0;
  for (const auto& bcTypePair: m_blockBCInfo)
    {
      if (bcTypePair.second.m_type & CRDparam::DomainBCTypeCNSCBCOutflow)
        {
          CBCTest = true;
        }
    }
  if (CBCTest)
    {
      CRD::msg << "Outflow smoothing beta\n" << m_CBCbeta << CRD::var;
      CRD::msg << "Outflow smoothing sigma\n" << m_CBCsigma << CRD::var;
      CRD::msg << "Inflow smoothing\n" << m_etaMax << CRD::var;
      CRD::msg << "Inflow mass fractions smoothing\n" << m_etaCN << CRD::var;
    }
  if (CRDparam::numAMRLevel() == 1)
    {
      return;
    }

  CRD::msg.newline();
  const int numTags = m_tagVar.size();
  CRD::msg << "Tagging information:" << CRD::body;
  for (int tag = 0; tag != numTags; ++tag)
    {
      CRD::msg << "Tagging method " << tag+1 << "\n";
      if (m_tagType[tag] == 0)
        {
          CRD::msg << m_tagVar[tag] << " gradient" << CRD::var;
        }
      else if (m_tagType[tag] == 1)
        {
          CRD::msg << m_tagVar[tag] << " value" << CRD::var;
        }
      else if (m_tagType[tag] == 2)
        {
          CRD::msg << "vorticity" << CRD::var;
        }
      else if (m_tagType[tag] == 3)
        {
          CRD::msg << m_tagVar[tag] << " vms" << CRD::var;
        }
      else if (m_tagType[tag] == 4)
        {
          CRD::msg << "vms_vorticity" << CRD::var;
        }
      CRD::msg << "Tagging threshold " << tag+1 << "\n" << m_threshold[tag]
               << CRD::var;
    }
  CRD::msg.newline();
}
