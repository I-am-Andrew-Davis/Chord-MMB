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
 * \file CNSIBCShear.cpp
 *
 * \brief Member functions for CNSIBCShear
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
#include "CONSTANTS.H"

//----- Internal -----//

#include "CNSIBCShear.H"
#include "ChordInput.H"
#include "CRDPhysics.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodBaseBox.H"
#include "TagMethodGradient.H"


/*******************************************************************************
 *
 * Class CNSIBCShear: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCShear::CNSIBCShear()
  :
  CNSIBCGeneralizedSingleBlock(),
  m_velConst(100.*RealVect::Unit)
{
  //--Read any BC info

  readBCInfo();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCShear::~CNSIBCShear()
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
CNSIBCShear::IBCName() const
{
  return "Shear problem";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCShear::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg.setPrecFloatSN(4);
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Initialize \<U\>
/** Sets the initial state for a solution.  This routine must compute
 *  \<U\>.
 *  \param[out] a_U     State on the level to be initialized in this
 *                      routine
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_time  Initial time
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

void
CNSIBCShear::initialize(LevelData<FArrayBox>&      a_U,
                        LevelGridMetrics&          a_gridMetrics,
                        const LayoutData<FluxBox>& a_unitNormals,
                        const Real                 a_time,
                        const int                  a_level) const
{
  const int numSpecies = CRDparam::g_numSpecies;
  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int numGhosts = 2;
  const RealVect domLen = CRDparam::g_domainLength;
  const RealVect domStart = CRDparam::g_domainOrigin;
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get the domain for the block
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box2Dom = grow(box, numGhosts);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, numGhosts - 1);
      box1Dom &= blockDomain;
      
      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));
      // Get physical coordinates
      FABSTACKTEMP(XiFab, box2Dom, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, box2Dom, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(box2Dom, XiFab, XFab, blockCoordSys);
      // Average values of U (required on valid +1 except at physical
      // boundaries)
      FArrayBox& UFab = a_U[dit];
      // Pointwise values of U (required on valid +2 except at physical
      // boundaries)
      CH_assert(UFab.box().contains(box2Dom));
      FABSTACKTEMP(Wc, box2Dom, numWVar);
      Wc.setVal(0.);
      MD_ARRAY_RESTRICT(arrW, Wc);
      MD_ARRAY_RESTRICT(arrX, XFab);
      MD_BOXLOOP(box2Dom, i)
        {
          RealVect xloc(D_DECL(arrX[MD_IX(i,0)] - domStart[0],
                               arrX[MD_IX(i,1)] - domStart[1],
                               arrX[MD_IX(i,2)] - domStart[2]));
          Real rho = m_initRho;
          arrW[MD_IX(i,rhoIndx)] = rho;
          if (numSpecies > 0)
            {
              for (int spec = 0; spec != numSpecies; ++spec)
                {
                  int wComp = wCompStart + spec;
                  arrW[MD_IX(i,wComp)] = m_initMassFraction[spec];
                }
            }
          arrW[MD_IX(i,presIndx)] = m_initP;
          arrW[MD_IX(i,tempIndx)] = m_initT;
          D_TERM(arrW[MD_IX(i, velIndx)] = m_initVel[0] +
                 m_velConst[0]*std::cos(2.*Pi*xloc[1]/domLen[1]);,
                 arrW[MD_IX(i, velIndx+1)] = m_initVel[1] +
                 m_velConst[1]*std::cos(2.*Pi*xloc[0]/domLen[0]);,
                 arrW[MD_IX(i, velIndx+2)] = 0.;);
        }
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         box2Dom);
      fourthOrderAverageCell(a_U[dit], blockDomain, box1Dom);
    }
}

/*--------------------------------------------------------------------*/
//  Set the tagging method if one can be associated with IBC
/**
 *  \param[in]  a_tagBufferSize
 *                      Requested tag buffer size (should be
 *                      respected).
 *  \note
 *  <ul>
 *    <li> Allocate levels and methods with 'new' and do not delete
 *  </ul>
 *//*-----------------------------------------------------------------*/

TagLevelFactory*
CNSIBCShear::setTagMethod(const int a_tagBufferSize)
{
  TagLevel* tagLevel = new TagLevel;
  IntVect tagLo(IntVect::Zero);
  IntVect tagHi(CRDparam::g_domainBaseSize);
  tagLo = tagHi/4;
  tagHi = tagHi*3/4 - IntVect::Unit;
  // Set tag box
  Box tagBox(tagLo, tagHi);
  tagLevel->appendTagMethod(new TagMethodBaseBox(tagBox));
  return new TagLevelFactory(tagLevel);
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

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
CNSIBCShear::readBCInfo()
{
  ParmParse ppIBC("ibc");
  std::vector<Real> inputVel(SpaceDim, 0.);
  ppIBC.queryarr("vel_factor", inputVel, 0, SpaceDim);
  SpaceDimArray<Real, Real>::loadFromArray(m_velConst.dataPtr(),
                                           &inputVel.front());
  m_readInput = true;
}
