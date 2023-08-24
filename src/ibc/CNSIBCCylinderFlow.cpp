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
 * \file CNSIBCCylinderFlow.cpp
 *
 * \brief Member functions for CNSIBCCylinderFlow
 *
 *//*+*************************************************************************/
#include <fstream>
//----- Chombo Library -----//
#include "CONSTANTS.H"
#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
 
//----- Internal -----//
#include "CNSIBCGeneralizedSingleBlock.H"
#include "CNSIBCCylinderFlow.H"
#include "ChordInput.H"
#include "CNSIBCF_F.H"
#include "CRDPhysics.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodValue.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"


/*******************************************************************************
 *
 * Class CNSIBCCylinderFlow: member definitions
 *
 *         (f)
 *     -----------                   ______          
 *     |         |             (f) /        \       
 *     |         |                /    _     \   
 * (p) |         | (p)           |    / \_____|
 *     |         |               |    \_/ (p) |
 *     |         |                \   (w)    /      
 *     |         |                 \ ______ /           
 *     -----------
 *         (w)
 *
 * (p) - periodic
 * (f) - farfield
 * (w) - wall
 * Used with the annulus mapping the idea is for the wall portion to be the
 * inner cylinder, the farfield is the outer, and the periodic connects it
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_farfieldVelocity
 *                      Velocity of the far field flow
 *//*-----------------------------------------------------------------*/

CNSIBCCylinderFlow::CNSIBCCylinderFlow()
  :
  CNSIBCGeneralized(),
  m_tagPerc(8)
{
  //--Read any BC info
  readBCInfo();
}


/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/
CNSIBCCylinderFlow::~CNSIBCCylinderFlow()
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Return a name describing the IBC
/** \return             Name of IBC
 *//*-----------------------------------------------------------------*/
const char *const
CNSIBCCylinderFlow::IBCName() const
{
  return "Navier-Stokes flow over cylinder";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/
void
CNSIBCCylinderFlow::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Tagging threshold\n(";
  for (int i = 0; i != CRDparam::numAMRLevel()-1; ++i)
    {
      if (i != 0)
        {
          CRD::msg << ',';
        }
      CRD::msg << m_threshold[i];
    }
  CRD::msg << ')' << CRD::var;
  CRD::msg << "Radial tagging ratio\n" << m_tagPerc << CRD::var;
  CRD::msg << "Radial in direction\n" << m_tagPerc << CRD::var;
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Set the tagging method if one can be associated with IBC
/** Fixed radial refinement, and tagging on vorticity
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
CNSIBCCylinderFlow::setTagMethod(const int a_tagBufferSize)
{ 
  const int maxAMRLevel = CRDparam::maxAMRLevel();
  TagLevel* tagLevel;
  std::vector<TagLevel*> tagLevelVec;
  std::vector<int> levelMapVec;
  tagLevelVec.resize(maxAMRLevel);
  levelMapVec.resize(maxAMRLevel);
  // still need to return an empty tag level even if no amr
  if (maxAMRLevel == 0)
    {
      tagLevelVec.resize(1);
      levelMapVec.resize(1);
      levelMapVec[0] = 0;
      tagLevelVec[0] = new TagLevel;
    }
  else{
    
    for (int l = 0; l != maxAMRLevel; ++l)
      {
        // create the tag level
        levelMapVec[l] = l;
        tagLevel = new TagLevel;
        // tag radial direction
        IntVect numCells = CRDparam::g_domainBaseSize;
        Box tagBox(IntVect::Zero, numCells);
        if (m_tagSide == 0)
          {
            tagBox.setBig(m_tagDir,
                          std::floor(numCells[m_tagDir]*m_tagPerc));
          }
        else
          {
            tagBox.setSmall(m_tagDir,
                            std::floor(numCells[m_tagDir]*(1.0-m_tagPerc)));
          }

        tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
      
        // Ensures no tagging occurs near the farfield boundary
        std::vector<Box> restrictBoxes(1);
        Box noTagBox(IntVect::Zero, numCells);
        if (m_tagSide == 0)
          {
            noTagBox.setSmall(m_tagDir,
                              std::floor(numCells[m_tagDir]*(1.0-m_tagPerc)));
          }
        else
          {
            noTagBox.setBig(m_tagDir,
                            std::floor(numCells[m_tagDir]*m_tagPerc));
          }
      
        // do vorticity taging above specified threshold
        tagLevel->appendTagMethod(new TagMethodValue(
                                    -1,
                                    m_threshold[l],
                                    1.0e+99,
                                    true,
                                    TagMethodValue::ValueTypeVorticity));
        // Add in tag buffer
        tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
        tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
        tagLevelVec[l] = tagLevel;
      }
  }
  // Return the factory
  return new TagLevelFactory(tagLevelVec, levelMapVec);
}

/*==============================================================================
 * Private member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Read any information related to the IBC from input
/** \note
 *  <ul>
 *    <li> No output should be printed aside from errors
 *  </ul>
 *//*-----------------------------------------------------------------*/
void
CNSIBCCylinderFlow::readBCInfo()
{
  //CNSIBCGeneralizedSingleBlock::readBCInfo();
  //CH_assert(!m_readInput);
  //CRD::msg.setTerminateOnError(false);  // Disable terminate on error
  ParmParse ppIBC("ibc");
  ParmParse ppCS("coordsys");

  Real diam;
  ppCS.get("inner_radius", diam);
  diam *= 2.0;

  //--Threshold per level of vorticity magnitude for refinement
  // size threshold vect
  m_threshold.resize(CRDparam::maxAMRLevel(), -1);
  if (ppIBC.contains("tag_threshold"))
    {
      ppIBC.getarr("tag_threshold", m_threshold, 0, CRDparam::maxAMRLevel());
      for (int i = 0; i != CRDparam::maxAMRLevel(); ++i)
        {
          if (m_threshold[i] < 0.)
            {
              CRD::msg << "Input (CylinderFlow IBC): 'tag_threshold' must be > 0!"
                       << CRD::error;
            }
        }
    }

  //--Percent to tag in radial direction
  ppIBC.query("radial_tagging", m_tagPerc);
  if (m_tagPerc < 0 || m_tagPerc > 1)
    {
      CRD::msg << "Input (CylinderFlow IBC): 'radial_tagging' must be > 0"
               << " and < 1!"<< CRD::error;
    }
  //--Define the radial direction
  ppIBC.query("tag_dir", m_tagDir);
  if (m_tagDir > SpaceDim || m_tagDir < 0)
    {
      CRD::msg << "Input (CylinderFlow IBC): 'tag_dir' must be in [0, SpaceDim]" << CRD::error;
    }
  //--Define the radial tagging side
  ppIBC.query("tag_side", m_tagSide);
  //--Enable terminate on error
  CRD::msg.setTerminateOnError(true);  
  if (CRD::msg.numErrors() > 0)
    {
      CRD::msg << "There " << CRD::wasNumError << " reading the ibc section of "
        "the input file.  Aborting now." << CRD::error;
    }
  m_readInput = true;
}
