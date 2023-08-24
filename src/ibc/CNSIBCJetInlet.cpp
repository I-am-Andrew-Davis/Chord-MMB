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
 * \file CNSIBCJetInlet.cpp
 *
 * \brief Member functions for CNSIBCJetInlet
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
#include "LoHiCenter.H"
 
//----- Internal -----//

#include "CNSIBCJetInlet.H"
#include "CNSIBCCombustionTestF_F.H"
#include "ViscousTensor4thOrderOpF_F.H"
#include "CNSIBCF_F.H"
#include "ChordInput.H"
#include "CRDPhysics.H"


/*******************************************************************************
 *
 * Class CNSIBCJetInlet: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCJetInlet::CNSIBCJetInlet()
  :
  CNSIBCGeneralizedSingleBlock(),
  m_flatplateSetup(false),
  m_slipWallFrac(-1.) // Default means there is no slip wall
{
  readBCInfo();
  int jetDir = 0;
  int jetWallDir = 1;
  // if (m_domainBC[1][0].type & CRDparam::DomainBCTypeCNSCBC ||
  //     m_domainBC[1][0].type & CRDparam::DomainBCTypeDirichlet ||
  //     m_domainBC[1][0].type & CRDparam::DomainBCTypeExtrapolated)
  //   {
  //     jetDir = 1;
  //     jetWallDir = 0;
  //   }
  setDomainBC(jetWallDir, 0, CRDparam::DomainBCTypeMixed);
  if (m_flatplateSetup)
    {
      setDomainBC(jetWallDir, 1, CRDparam::DomainBCTypeSlipWall);
    }
  else
    {
      setDomainBC(jetWallDir, 1, CRDparam::DomainBCTypeMixed);
    }
  if (SpaceDim == 3)
    {
      setDomainBC(2, 0, CRDparam::DomainBCTypeSlipWall);
      setDomainBC(2, 1, CRDparam::DomainBCTypeSlipWall);
    }
  IntVect loNoSlip(IntVect::Zero);
  m_slipWall.define(-3*IntVect::Unit, -IntVect::Unit);
  if (m_slipWallFrac > 0.)
    {
      // Define the slip wall
      IntVect lo(IntVect::Zero); 
      IntVect hi(CRDparam::g_domainBaseSize);
      hi[jetDir] *= m_slipWallFrac;
      m_slipWall.define(lo, hi);
      lo[jetDir] = m_slipWall.bigEnd()[jetDir] + 1;
      loNoSlip[jetDir] = lo[jetDir];
    }
  {
    // Define the high jet box
    IntVect hi(CRDparam::g_domainBaseSize);
    m_noSlipWall.define(loNoSlip, hi);
  }
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCJetInlet::~CNSIBCJetInlet()
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
CNSIBCJetInlet::IBCName() const
{
  return "Inlet jet case";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCJetInlet::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg << "Jet problem\n";
  if (m_flatplateSetup)
    {
      CRD::msg << "flatplate" << CRD::var;
    }
  else
    {
      CRD::msg << "channel" << CRD::var;
    }
  CRD::msg.setFloatDefault();
  CRD::msg << "Leading slip wall fraction\n";
  if (m_slipWallFrac > 1.)
    {
      CRD::msg << "entire wall" << CRD::var;
    }
  else if (m_slipWallFrac < 0.)
    {
      CRD::msg << "none" << CRD::var;
    }
  else
    {
      CRD::msg << m_slipWallFrac*100. << "%" << CRD::var;
    }
  CRD::msg.newline();
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set the boxes that are inflow or outflow if mixed BC is being used
/** \param[in]  a_boundaryFaceBox
 *                      Box of boundary being operated on
 *  \param[in]  a_dir   Direction of the boundary
 *  \param[in]  a_side  LoHi side of the boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[out] a_boxVect
 *                      Vector of boxes
 *//*-----------------------------------------------------------------*/

int
CNSIBCJetInlet::setMixedBC(
  const Box&           a_boundaryFaceBox,
  const BoundaryIndex& a_bcIdx,
  const Box&           a_disjointBox,
  LevelGridMetrics&    a_gridMetrics,
  const Real           a_time,
  const int            a_level,
  Vector<Box>&         a_boxVect,
  Vector<BCInfo>&      a_domainBC) const
{
  const int lohiSign = sign(a_bcIdx.m_side);
  a_boxVect.resize(2);
  a_domainBC.resize(2);
  Box slipWall(m_slipWall);
  Box noSlipWall(m_noSlipWall);
  slipWall.refine(CRDparam::g_refFromBase[a_level]);
  noSlipWall.refine(CRDparam::g_refFromBase[a_level]);
  slipWall.shiftHalf(a_bcIdx.m_dir, lohiSign);
  noSlipWall.shiftHalf(a_bcIdx.m_dir, lohiSign);
  a_boxVect[0] = slipWall & a_boundaryFaceBox;
  a_domainBC[0].m_type = CRDparam::DomainBCTypeSlipWall;
  a_domainBC[0].m_order = 1;
  a_boxVect[1] = noSlipWall & a_boundaryFaceBox;
  a_domainBC[1].m_type = CRDparam::DomainBCTypeAdiabaticWall;
  a_domainBC[1].m_order = 1;
  return 0;
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
CNSIBCJetInlet::readBCInfo()
{
  ParmParse ppIBC("ibc");
  ppIBC.query("percent_slip_wall", m_slipWallFrac);
  m_flatplateSetup = false;
  std::string setupName;
  ppIBC.query("jet_problem_type", setupName);
  if (setupName == "FlatPlate" || setupName == "flatplate" ||
     setupName == "Flatplate")
    {
      m_flatplateSetup = true;
    }
}
