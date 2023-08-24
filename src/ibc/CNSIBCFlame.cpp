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
 * \file CNSIBCFlame.cpp
 *
 * \brief Member functions for CNSIBCFlame
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "ParmParse.H"
#include "LevelGridMetrics.H"
#include "DataTemp.H"
#include "FourthOrderUtil.H"
#include "CONSTANTS.H"
 
//----- Internal -----//

#include "CNSIBCFlame.H"
#include "CNSIBCGeneralizedSingleBlock.H"
#include "CNSIBCFlameF_F.H"
#include "CRDPhysics.H"
#include "ChordInput.H"


/*******************************************************************************
 *
 * Class CNSIBCFlame: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** 
 *//*-----------------------------------------------------------------*/

CNSIBCFlame::CNSIBCFlame()
  :
  CNSIBCGeneralizedSingleBlock(),
  m_flameDir(0),
  m_varDir(1),
  m_pertMag(0.),
  m_pertPeriod(1.),
  m_pertT(-1.),
  m_pertP(-1.),
  m_pertStartPoint(-9.*RealVect::Unit),
  m_pertEndPoint(-10.*RealVect::Unit),
  m_pertRadius(0.),
  m_pertCr(1.E8),
  m_pertStartTime(0.),
  m_pertEndTime(-1.),
  m_Cr(1.E8),
  m_rho2(-1.),
  m_p2(-1.),
  m_t2(-1.),
  m_loEnd2(-1.*RealVect::Unit),
  m_hiEnd2(-1.*RealVect::Unit),
  m_cirFlame(0)
{
  readBCInfo();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCFlame::~CNSIBCFlame()
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
CNSIBCFlame::IBCName() const
{
  return "General flame problem";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCFlame::writeIBCInfo() const
{
  CNSIBCCombustionReference::writeTagInfo();
  CRD::msg << "Region 1 temperature\n" << m_initT << CRD::var;
  CRD::msg << "Region 2 temperature\n" << m_t2 << CRD::var;
  CRD::msg << "Region 1 initial mass fractions\n(" << m_initMassFraction[0];
  for (int species = 1; species != CRDparam::g_numSpecies; ++species)
    {
      CRD::msg << ", " << m_initMassFraction[species];
    }
  CRD::msg << ")" << CRD::var;
  CRD::msg << "Region 2 initial mass fractions\n(" << m_region2MassFraction[0];
  for (int species = 1; species != CRDparam::g_numSpecies; ++species)
    {
      CRD::msg << ", " << m_region2MassFraction[species];
    }
  CRD::msg << ")" << CRD::var;
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Region 1 pressure\n" << m_initP << CRD::var;
  CRD::msg << "Region 2 pressure\n" << m_p2 << CRD::var;
  CRD::msg << "Interface sharpness\n" << m_Cr << CRD::var;
  CRD::msg.setFloatDefault();
  CRD::msg << "Region 2 end\n(" << m_loEnd2 << ", " << m_hiEnd2
           << ")" << CRD::var;
  if (m_pertStartTime < m_pertEndTime)
    {
      CRD::msg << "Inlet perturbation info" << CRD::h2;
      CRD::msg << "Start time\n" << m_pertStartTime << CRD::var;
      CRD::msg << "End time\n" << m_pertEndTime << CRD::var;
      CRD::msg << "Period\n" << m_pertPeriod << CRD::var;
      CRD::msg << "Magnitude\n" << m_pertMag << CRD::var;
    }
  else if (m_pertStartPoint[m_flameDir] < m_pertEndPoint[m_flameDir])
    {
      CRD::msg << "Initial perturbation info" << CRD::h2;
      CRD::msg << "Start location\n" << m_pertStartPoint << CRD::var;
      CRD::msg << "Radius\n" << m_pertRadius << CRD::var;
      CRD::msg << "End location\n" << m_pertEndPoint << CRD::var;
      CRD::msg << "Temperature\n" << m_pertT << CRD::var;
      CRD::msg << "Pressure\n" << m_pertP << CRD::var;
    }
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
CNSIBCFlame::initialize(LevelData<FArrayBox>&      a_U,
                        LevelGridMetrics&          a_gridMetrics,
                        const LayoutData<FluxBox>& a_unitNormals,
                        const Real                 a_time,
                        const int                  a_level) const
{
  // Set the initial values
  const int numSpecies = CRDparam::g_numSpecies;
  const int numWVar = CRDparam::g_CRDPhysics->numPrimitive();
  const int WvelIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
  const int wCompStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  int numGhosts = CRDparam::queryGhosts(CRDparam::NumGhostUcellAvg);
  Real rhoVal = m_initRho;
  Real sameCnTest = 10.;
  bool sameCnBool = true;
  for (int spec = 0; spec != numSpecies; ++spec)
    {
      sameCnTest = std::min(
        m_initMassFraction[spec] - m_region2MassFraction[spec], sameCnTest);
    }
  if (sameCnTest == 0.)
    {
      sameCnBool = false;
    }
  // Iterate over all boxes on this level
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box box = a_U.disjointBoxLayout()[dit];

      // Get coordinate system and domain for the block
      const BlockCoordSys& blockCoordSys = *(a_gridMetrics.getCoordSys(box));
      const BlockDomain& blockDomain =
        a_gridMetrics.getCoordSys().problemDomain(box);

      // Working set boxes
      Box box2Dom = grow(box, numGhosts);
      box2Dom &= blockDomain;
      Box box1Dom = grow(box, numGhosts - 1);
      box1Dom &= blockDomain;
      Box initBox = box2Dom;

      // Get physical coordinates
      FABSTACKTEMP(XiFab, initBox, SpaceDim);  // Cartesian coordinates
      FABSTACKTEMP(XFab, initBox, SpaceDim);   // Physical coordinates
      this->CNSIBC::getCellCoordinates(initBox, XiFab, XFab, blockCoordSys);

      // Average values of U (required on valid +1 except at physical
      // boundaries)
      FArrayBox& UFab = a_U[dit];
      UFab.setVal(0.);
      // Create a FAB of the initial primitive variables
      FABSTACKTEMP(Wc, initBox, numWVar);
      Wc.setVal(0.);
      // Set velocity to zero
      Wc.setVal(0., initBox, WvelIndx, SpaceDim);
      if (m_cirFlame == 1)
        {
          const RealVect center = CRDparam::g_domainOrigin +
            CRDparam::g_domainLength/2.;
          const Real Len = CRDparam::g_domainLength[m_flameDir];
          MD_ARRAY_RESTRICT(arrW, Wc);
          MD_ARRAY_RESTRICT(arrX, XFab);
          MD_BOXLOOP(initBox, i)
            {
              RealVect curLoc(D_DECL(arrX[MD_IX(i, 0)],
                                     arrX[MD_IX(i, 1)],
                                     arrX[MD_IX(i, 2)]));
              RealVect cirPoint = curLoc - center;
              Real curRadius = std::sqrt(
                D_TERM(cirPoint[0]*cirPoint[0],
                       + cirPoint[1]*cirPoint[1],));
              Real perc = 0.5*
                (1. + std::tanh((m_pertRadius - curRadius)*m_Cr/Len));
              Real summf = 0.;
              Real sumside = 0.;
              for (int spec = 0; spec != numSpecies; ++spec)
                {
                  const int wComp = wCompStart + spec;
                  const Real cn1 = m_initMassFraction[spec];
                  const Real cn2 = m_region2MassFraction[spec];
                  Real cnval = cn1 + perc*(cn2 - cn1);
                  arrW[MD_IX(i, wComp)] = cnval;
                  summf += (cnval - cn2)/(cn1 - cn2);
                  sumside += 1.;
                }
              summf /= sumside;
              Real tempVal = m_initT;
              Real presVal = m_initP;
              Real rhoVal = m_initRho;
              if (sameCnBool)
                {
                  tempVal = 1./(summf/m_initT + (1. - summf)/m_t2);
                  presVal = 1./(summf/m_initP + (1. - summf)/m_p2);
                  rhoVal = 1./(summf/m_initRho + (1. - summf)/m_rho2);
                }
              else
                {
                  tempVal = m_initT + perc*(m_t2 - m_initT);
                  presVal = m_initP + perc*(m_p2 - m_initP);
                  rhoVal = m_initRho + perc*(m_rho2 - m_initRho);
                }
              arrW[MD_IX(i, tempIndx)] = tempVal;
              arrW[MD_IX(i, presIndx)] = presVal;
              arrW[MD_IX(i, rhoIndx)] = rhoVal;
              D_TERM(arrW[MD_IX(i, WvelIndx)] = m_initVel[0];,
                     arrW[MD_IX(i, WvelIndx+1)] = m_initVel[1];,
                     arrW[MD_IX(i, WvelIndx+2)] = m_initVel[2];);
            }
        }
      else
        {
          const Real Len = CRDparam::g_domainLength[m_flameDir];
          MD_ARRAY_RESTRICT(arrW, Wc);
          MD_ARRAY_RESTRICT(arrX, XFab);
          MD_BOXLOOP(initBox, i)
            {
              RealVect curLoc(D_DECL(arrX[MD_IX(i, 0)],
                                     arrX[MD_IX(i, 1)],
                                     arrX[MD_IX(i, 2)]));
              Real locdir = curLoc[m_flameDir];
              Real yloc = curLoc[m_varDir];
              Real diff1 = std::abs(m_loEnd2[m_flameDir] - locdir);
              Real diff2 = std::abs(m_hiEnd2[m_flameDir] - locdir);
              Real tLen = 0.;
              if (diff1 < diff2)
                {
                  tLen = locdir - m_loEnd2[m_flameDir];
                }
              else
                {
                  tLen = m_hiEnd2[m_flameDir] - locdir;
                }
              Real perc = 0.5*(1. + std::tanh(m_Cr*tLen/Len));
              Real summf = 0.;
              Real sumside = 0.;
              for (int spec = 0; spec != numSpecies; ++spec)
                {
                  const int wComp = wCompStart + spec;
                  const Real cn1 = m_initMassFraction[spec];
                  const Real cn2 = m_region2MassFraction[spec];
                  Real cnval = cn1 + perc*(cn2 - cn1);
                  arrW[MD_IX(i, wComp)] = cnval;
                  if (cn2 != 0 && cn1 != cn2)
                    {
                      summf += (cnval - cn2)/(cn1 - cn2);
                      sumside += 1.;
                    }
                }
              summf /= std::max(1., sumside);
              Real meanTemp = 1./(summf/m_initT + (1. - summf)/m_t2);
              arrW[MD_IX(i, tempIndx)] = meanTemp;
              Real meanPres = 1./(summf/m_initP + (1. - summf)/m_p2);
              arrW[MD_IX(i, presIndx)] = meanPres;
              RealVect cirPoint = curLoc - m_pertStartPoint;
              Real curRadius = std::sqrt(
                D_TERM(cirPoint[0]*cirPoint[0],
                       + cirPoint[1]*cirPoint[1],));
              Real pertRadius = 0.5*
                (1. + std::tanh((m_pertRadius - curRadius)*m_pertCr));
              Real otherP = m_initP;
              Real otherT = m_initT;
              if (m_pertStartPoint < m_hiEnd2 && m_pertStartPoint > m_loEnd2)
                {
                  otherP = m_p2;
                  otherT = m_t2;
                }
              if (pertRadius > 0. && pertRadius <= 1.)
                {
                  arrW[MD_IX(i, presIndx)] = otherP +
                    pertRadius*(m_pertP - otherP);
                  arrW[MD_IX(i, tempIndx)] = otherT +
                    pertRadius*(m_pertT - otherT);
                }
              if (locdir > m_pertStartPoint[m_flameDir] &&
                 locdir < m_pertEndPoint[m_flameDir] &&
                 yloc > m_pertStartPoint[m_varDir] &&
                 yloc < m_pertEndPoint[m_varDir])
                {
                  arrW[MD_IX(i, presIndx)] = m_pertP;
                  arrW[MD_IX(i, tempIndx)] = m_pertT;
                }
              if (rhoVal < 0.)
                {
                  arrW[MD_IX(i, rhoIndx)] = -1.;
                }
              else
                {
                  arrW[MD_IX(i, rhoIndx)] = rhoVal + perc*(m_rho2 - rhoVal);
                }
              D_TERM(arrW[MD_IX(i, WvelIndx)] = m_initVel[0];,
                     arrW[MD_IX(i, WvelIndx+1)] = m_initVel[1];,
                     arrW[MD_IX(i, WvelIndx+2)] = m_initVel[2];);
            }
        }
      CRDparam::g_CRDPhysics->initialize(UFab,
                                         Wc,
                                         a_gridMetrics,
                                         a_unitNormals[dit],
                                         dit(),
                                         box,
                                         initBox);
      fourthOrderAverageCell(UFab, blockDomain, box1Dom);
    }
}

/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Set the imposed (exterior or farfield) primitive state at flow BCC
//  Also adds time and space dependent perturbations
/** \param[in]  a_Wface Primitive state on wall from interior scheme
 *                      (face-centered)
 *  \param[out] a_Wface Primitive state at exterior to boundary
 *  \param[in]  a_boundaryFaceBox
 *                      Box of faces on boundary to adjust
 *  \param[in]  a_Wcell Primitive state in cells at interior of
 *                      domain.  The cells have been shifted by half
 *                      so that the first interior layer of cells
 *                      overlaps a_Wface.
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces.  (This is the same basis used
 *                      for solving the Riemann problem on mapped
 *                      grids
 *  \param[in]  a_dir   Direction of the boundary
 *  \param[in]  a_side  LoHi side of the boundary
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current boundary condition being set (for mixed BCs)
 *//*-----------------------------------------------------------------*/

void
CNSIBCFlame::setImposedBCprimState(
  FArrayBox&           a_Wface,
  const Box&           a_boundaryFaceBox,
  const FArrayBox&     a_Wcell,
  const FArrayBox&     a_unitNormalBasisFab,
  const BoundaryIndex& a_bcIdx,
  const Box&           a_disjointBox,
  LevelGridMetrics&    a_gridMetrics,
  const Real           a_time,
  const int            a_level,
  const BCInfo&        a_domT) const
{
  CNSIBCGeneralizedSingleBlock::setImposedBCprimState(a_Wface,
                                           a_boundaryFaceBox,
                                           a_Wcell,
                                           a_unitNormalBasisFab,
                                           a_bcIdx,
                                           a_disjointBox,
                                           a_gridMetrics,
                                           a_time,
                                           a_level,
                                           a_domT);
  // Add inlet perturbations
  if (a_time > m_pertStartTime && a_time < m_pertEndTime &&
     (a_domT.m_type & CRDparam::DomainBCTypeInflow))
    {
      addInletPerturbations(a_Wface,
                            a_boundaryFaceBox,
                            a_bcIdx,
                            a_disjointBox,
                            a_gridMetrics,
                            a_time,
                            a_level);
    }
}

/*--------------------------------------------------------------------*/
//  Fill in the profiles for the CNSCBC boundary condition
//  Also adds time and space dependent perturbations
/** \param[out] a_BCProfile
 *                      Profile of boundary values
 *  \param[in]  a_boundaryBox
 *                      Box of locations at the boundary
 *  \param[in]  a_Wcell Cell-averaged primitive variables
 *  \param[in]  a_unitNormalBasis
 *                      A basis with the 'a_dir'th row normal to the
 *                      'a_dir' faces
 *  \param[in]  a_dir   Normal direction to the boundary
 *  \param[in]  a_side  LoHi side of boundary
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *  \param[in]  a_domT  Current CNSCBC boundary condition
 *//*-----------------------------------------------------------------*/

void
CNSIBCFlame::setCNSCBCProfiles(
  FArrayBox&           a_BCProfile,
  const Box&           a_boundaryBox,
  const FArrayBox&     a_Wcell,
  const FArrayBox&     a_unitNormalBasisFab,
  const BoundaryIndex& a_bcIdx,
  const Box&           a_disjointBox,
  LevelGridMetrics&    a_gridMetrics,
  const Real           a_time,
  const int            a_level,
  const BCInfo&        a_domT) const
{
  CNSIBCGeneralizedSingleBlock::setCNSCBCProfiles(a_BCProfile,
                                       a_boundaryBox,
                                       a_Wcell,
                                       a_unitNormalBasisFab,
                                       a_bcIdx,
                                       a_disjointBox,
                                       a_gridMetrics,
                                       a_time,
                                       a_level,
                                       a_domT);
  // Add inlet perturbations
  if (a_time > m_pertStartTime && a_time < m_pertEndTime &&
     (a_domT.m_type & CRDparam::DomainBCTypeCNSCBCInflow))
    {
      addInletPerturbations(a_BCProfile,
                            a_boundaryBox,
                            a_bcIdx,
                            a_disjointBox,
                            a_gridMetrics,
                            a_time,
                            a_level);
    }
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
CNSIBCFlame::readBCInfo()
{
  // FIXME: This IBC file is a mess, should be cleaned up
  int numSpecies = CRDparam::g_numSpecies;
  ParmParse ppIBC("ibc");
  // Interface between regions 1 and 2
  ppIBC.query("flame_dir", m_flameDir);
  ppIBC.query("pert_dir", m_varDir);
  ppIBC.query("interface_sharpness", m_Cr);
  // Mass fractions in region 2
  m_region2MassFraction.resize(numSpecies);
  m_region2MassFraction.assign(numSpecies,0.);
  // Possible to use alpha_h2, alpha_co2, and phi instead of specifying
  // each mass fraction
  if (ppIBC.contains("alpha_h2"))
    {
      m_initMassFraction.assign(numSpecies, 0.);
      m_region2MassFraction.assign(numSpecies, 0.);
      Real alphaH2;
      Real alphaCO2;
      Real phi;
      ppIBC.get("alpha_h2",alphaH2);
      ppIBC.get("alpha_co2",alphaCO2);
      ppIBC.get("equiv_ratio",phi);
      // Retrieve specific components
      int h2comp = -1;
      int cocomp = -1;
      int co2comp = -1;
      int o2comp = -1;
      int n2comp = -1;
      const Real mwO2 = 31.9988;
      const Real mwN2 = 28.01348;
      const Real mwCO = 28.0104;
      const Real mwCO2 = 44.0098;
      const Real mwH2 = 2.01588;
      for (int i = 0; i != numSpecies; ++i)
        {
          if (CRDparam::g_speciesNames[i] == "H2")
            {
              h2comp = i;
            }
          else if (CRDparam::g_speciesNames[i] == "CO")
            {
              cocomp = i;
            }
          else if (CRDparam::g_speciesNames[i] == "CO2")
            {
              co2comp = i;
            }
          else if (CRDparam::g_speciesNames[i] == "O2")
            {
              o2comp = i;
            }
          else if (CRDparam::g_speciesNames[i] == "N2")
            {
              n2comp = i;
            }
        }
      if (n2comp == -1 || h2comp == -1)
        {
          CRD::msg << "Input (Flame IBC): must be "
                   << "have N2, O2, and H2!" << CRD::error;
        }
      else if (alphaH2 < 1. || alphaCO2 > 0.)
        {
          if (co2comp == -1 || cocomp == -1)
            {
              CRD::msg << "Input (Flame IBC): must be "
                       << "have CO and CO2!" << CRD::error;
            }
        }
      const Real nCO2 = alphaCO2;
      const Real nH2 = (1. - alphaCO2)*alphaH2;
      const Real nCO = (1. - alphaCO2)*(1. - alphaH2);
      const Real nO2 = (nH2 + nCO)/(2.*phi);
      const Real nN2 = nO2*3.76;
      const Real mwTot = mwO2*nO2 + mwH2*nH2 + mwCO*nCO + mwCO2*nCO2 + mwN2*nN2;
      m_initMassFraction[h2comp] = nH2*mwH2/mwTot;
      m_region2MassFraction[h2comp] = nH2*mwH2/mwTot;
      
      m_initMassFraction[o2comp] = nO2*mwO2/mwTot;
      m_region2MassFraction[o2comp] = nO2*mwO2/mwTot;
      if (alphaH2 < 1. || alphaCO2 > 0.)
        {
          m_initMassFraction[co2comp] = nCO2*mwCO2/mwTot;
          m_region2MassFraction[co2comp] = nCO2*mwCO2/mwTot;
      
          m_initMassFraction[cocomp] = nCO*mwCO/mwTot;
          m_region2MassFraction[cocomp] = nCO*mwCO/mwTot;
        }

      m_initMassFraction[n2comp] = nN2*mwN2/mwTot;
      m_region2MassFraction[n2comp] = nN2*mwN2/mwTot;
      // Set the inlet mass fractions to be the initial mass fractions
      BoundaryIndex bcIdx;
      bcIdx.m_block = 0;
      for (const auto side : EachSide)
        {
          bcIdx.m_side = side;
          for (const auto dir : EachDir)
            {
              bcIdx.m_dir = dir;
              if (getDomainBC(bcIdx).m_type & CRDparam::DomainBCTypeInlet)
                {
                  // auto bc = getDomainBC(bcIdx);
                  // FIXME update to CRDstate
                  // BCstate inletState = getDomainBCstate(bc);
                  // inletState.m_MF = m_initMassFraction;
                  // FIXME update to CRDstate
                  //setDomainBCstate(bc.m_stateIdx, inletState);
                }
            }
        }
    }
  if (ppIBC.contains("region_2_mfs"))
    {
      int massTest = assignMassFractions(m_region2MassFraction,
                                         "region_2_specs",
                                         "region_2_mfs");
      if (massTest == 1)
        {
          CRD::msg << "Input (Flame IBC): 'region_2_mass_fractions'"
                   << " must be equal to 1!" << CRD::error;
        }
    }
  m_p2 = m_initP;
  m_t2 = m_initT;
  m_rho2 = m_initRho;
  ppIBC.query("circular_flame", m_cirFlame);
  if (ppIBC.contains("lo_end") || m_cirFlame)
    {
      ppIBC.query("p2", m_p2);
      ppIBC.query("rho2", m_rho2);
      ppIBC.query("t2", m_t2);
      if ((m_t2*m_p2*m_rho2 > 0.) || (m_t2+m_p2+m_rho2 < 0))
        {
          CRD::msg << "Input (Flame IBC): 1 of the variables (rho, p, T) "
                   << "in region 2 must be < 0" << CRD::error;
        }
      m_pertT = m_t2;
      m_pertP = m_p2;
      if (m_cirFlame)
        {
          ppIBC.get("circle_radius", m_pertRadius);
        }
    }
  m_loEnd2 = -1.*RealVect::Unit;
  m_hiEnd2 = RealVect::Zero;
  if (ppIBC.contains("lo_end"))
    { 
      std::vector<Real> IBCloEnd(SpaceDim);
      ppIBC.queryarr("lo_end", IBCloEnd, 0, SpaceDim);
      SpaceDimArray<Real, Real>::loadFromArray(m_loEnd2.dataPtr(),
                                               &IBCloEnd.front());
      if (ppIBC.contains("hi_end"))
        { 
          std::vector<Real> IBChiEnd(SpaceDim);
          ppIBC.queryarr("hi_end", IBChiEnd, 0, SpaceDim);
          SpaceDimArray<Real, Real>::loadFromArray(m_hiEnd2.dataPtr(),
                                                   &IBChiEnd.front());
        }
    }
  if (ppIBC.contains("pert_start_time"))
    {
      ppIBC.get("pert_start_time", m_pertStartTime);
      m_pertEndTime = 1.E6;
      ppIBC.query("pert_end_time", m_pertEndTime);
      ppIBC.query("inlet_pert_mag", m_pertMag);
      ppIBC.query("inlet_pert_period", m_pertPeriod);
      if (m_pertEndTime < m_pertStartTime)
        {
          CRD::msg << "Input (Flame IBC): 'pert_start_time' cannot be after"
                   << " 'pert_end_time'" << CRD::error;
        }
    }
  else if (ppIBC.contains("pert_start_point") ||
          ppIBC.contains("pert_radius"))
    {
      std::vector<Real> IBCpertStart(SpaceDim);
      ppIBC.queryarr("pert_start_point", IBCpertStart, 0, SpaceDim);
      SpaceDimArray<Real, Real>::loadFromArray(m_pertStartPoint.dataPtr(),
                                               &IBCpertStart.front());
      std::vector<Real> IBCpertEnd(SpaceDim);
      ppIBC.queryarr("pert_end_point", IBCpertEnd, 0, SpaceDim);
      SpaceDimArray<Real, Real>::loadFromArray(m_pertEndPoint.dataPtr(),
                                               &IBCpertEnd.front());
      m_pertCr = m_Cr;
      ppIBC.query("pert_radius", m_pertRadius);
      ppIBC.query("pert_cr", m_pertCr);
      ppIBC.query("pert_temperature", m_pertT);
      ppIBC.query("pert_pressure", m_pertP);
    }
  m_readInput = true;
}

/*--------------------------------------------------------------------*/
//  Add inlet perturbations
/** \param[out] a_WFace Boundary values with perturbations added
 *  \param[in]  a_box   Box of locations at the boundary
 *  \param[in]  a_dir   Normal direction to the boundary
 *  \param[in]  a_side  LoHi side of boundary
 *  \param[in]  a_time  Current time
 *  \param[in]  a_level Grid level
 *//*-----------------------------------------------------------------*/

void
CNSIBCFlame::addInletPerturbations(FArrayBox&            a_WFace,
                                   const Box&            a_box,
                                   const BoundaryIndex&  a_bcIdx,
                                   const Box&            a_disjointBox,
                                   LevelGridMetrics&     a_gridMetrics,
                                   const Real            a_time,
                                   const int             a_level) const
{
  const int varDir = m_varDir;
  const Real domLength = CRDparam::g_domainLength[varDir];
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();
  Box intBox = adjCellBox(a_box, a_bcIdx.m_dir, Side::flip(a_bcIdx.m_side), 1);
  // FIXME - this should be single-block so any box should work.  Ideally
  //         the parent block domain should be an annotation in the box
  // FIXME - this should work as long as min box size is >= 8
  const int numGhostW =
    CRDparam::queryGhosts(CRDparam::NumGhostWfaceAvgTngBdry);
  intBox.grow(a_bcIdx.m_dir, numGhostW);
  intBox.grow(-numGhostW);
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  FABSTACKTEMP(XiFab, a_box, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, a_box, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(a_box, XiFab, XFab, blockCoordSys);
  MD_ARRAY_RESTRICT(arrWface, a_WFace);
  MD_ARRAY_RESTRICT(arrX, XFab);
  MD_BOXLOOP(a_box, i)
    {
      Real loc = arrX[MD_IX(i, varDir)];
      Real lRat = loc/domLength;
      Real pert = m_pertMag*std::sin((a_time/m_pertPeriod + lRat)*2.*Pi);
      arrWface[MD_IX(i, velIndx + varDir)] += pert;
    }
}
