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
 * \file TagMethodValue.cpp
 *
 * \brief Member functions for TagMethodValue
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "IntVectSet.H"
#include "Box.H"
#include "DataIndex.H"
#include "LoHiCenter.H"
#include "MappedLevelData.H"
#include "LevelGridMetrics.H"
#include "MultiBlockCoordSys.H"

//----- Internal -----//

#include "TagMethodValue.H"
#include "CRDparam.H"
#include "PolytropicPhysicsF_F.H"
#include "DataTemp.H"
#include "CRDPhysics.H"
#include "PatchMappedFunc.H"


/*******************************************************************************
 *
 * Class TagMethodValue: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_comp  Component to use for value
 *  \param[in]  a_loThreshold
 *                      Lower threshold for refinement.  Only cells with
 *                      values above this value will be tagged.
 *  \param[in]  a_upThreshold
 *                      Upper threshold for refinement.  Only cells with
 *                      values below this value will be tagged.
 *  \param[in]  a_usePrimitive
 *                      T - Values of primitive variables (default)
 *                      F - Values of conservative variables
 *  \param[in]  a_insideThreshold
 *                      T - Tag cells with values inside the thresholds
 *                      F - Tag cells with values outside the thresholds
 *  \param[in]  a_valueType
 *                      Type of value to tag:
 *                        ValueTypeConservative - conservative state variable
 *                        ValueTypePrimitive    - primitive state variables
 *                        ValueTypeVorticity    - vorticity magnitude
 *                        ValueTypeHydrocarbonFlame
 *                          - hydrocarbon flame surface detected via product of
 *                            mass fractions: OH*CH2O
 *  \param[in]  a_restrictBox
 *                      Vector of boxes for which refinement does not occur
 *                      These boxes should only be for the base grid
 *//*-----------------------------------------------------------------*/

TagMethodValue::TagMethodValue(const int        a_comp,
                               const Real       a_loThreshold,
                               const Real       a_upThreshold,
                               const bool       a_insideThreshold,
                               const ValueType  a_valueType,
                               std::vector<Box> a_restrictBox)
  :
  m_comp(a_comp),
  m_loThreshold(a_loThreshold),
  m_upThreshold(a_upThreshold),
  m_insideThreshold(a_insideThreshold),
  m_valueType(a_valueType),
  m_restrictBox(a_restrictBox),
  m_cOH(-1),
  m_cCH2O(-1)
{
  if (m_valueType == ValueTypeHydrocarbonFlame)
    {
      m_cOH   = CRDparam::g_CRDPhysics->speciesConsIndex("OH");
      if (m_cOH == -1)
        {
          CRD::msg << "Cannot find OH species required for tagging a "
            "hydrocarbon flame!" << CRD::error;
        }
      m_cCH2O = CRDparam::g_CRDPhysics->speciesConsIndex("CH2O");
      if (m_cOH == -1)
        {
          CRD::msg << "Cannot find CH2O species required for tagging a "
            "hydrocarbon flame!" << CRD::error;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

TagMethodValue::~TagMethodValue()
{
}

/*--------------------------------------------------------------------*/
//  Tag all cells with values either inside or outside of the thresholds
/** \param[out] a_tags  IntVectSet of cells to refine local to
 *                      'a_box'.
 *  \param[in]  a_box   Disjoint box in which cells should be marked
 *                      for refinement
 *  \param[in]  a_didx  Data index of the current box
 *  \param[in]  a_data  Data on the level, both mapped and physical.
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[out] a_gridMetrics
 *                      Only caches should be modified
 *  \param[in]  a_MBCoordSys
 *                      Coordinate system describing the mapping
 *  \param[in]  a_time  Current solution time
 *  \param[in]  a_level Index of the AMR level
 *//*-----------------------------------------------------------------*/

void
TagMethodValue::operator()(IntVectSet&               a_tags,
                           const Box&                a_box,
                           const DataIndex&          a_didx,
                           const MappedLevelData&    a_data,
                           LevelGridMetrics&         a_gridMetrics,
                           const MultiBlockCoordSys& a_MBCoordSys,
                           const Real                a_time,
                           const int                 a_level)
{
  const Interval velIntv = CRDparam::g_CRDPhysics->velocityInterval();
  const int compRho = CRDparam::g_CRDPhysics->densityIndex();
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_box);
  const Box box1 = grow(a_box, 1);
  Box box1Dom(box1);
  box1Dom &= blockDomain;
  int useComp = m_comp;

//--Get the data

  const FArrayBox& Ufab = a_data.getU(1,1)[a_didx];
  const FArrayBox* phiPtr = &Ufab;
  // Primitive state
  FABSTACKTEMP(Wfab, box1Dom, CRDparam::g_CRDPhysics->numPrimitive());
  // Scalar value (aliases primitive state, use one or the other)
  FArrayBox phiFab(Interval(0, 0), Wfab);

  switch (m_valueType)
    {
    case ValueTypeConservative:
      break;

//--Use a primitive state value

    case ValueTypePrimitive:
      CRDparam::g_CRDPhysics->consToPrim(Wfab, Ufab, box1Dom, Wfab);
      if (CRDparam::g_CRDPhysics->extraPrimInterval().contains(m_comp))
        {
          CRDparam::g_CRDPhysics->extraPrimitiveState(Wfab, box1Dom);
        }
      phiPtr = &Wfab;
      break;

//--Use vorticity

    case ValueTypeVorticity:
    {
      useComp = 0;
      FABSTACKTEMP(vort, a_box, PatchMappedFunc::m_numCurlComps);
      FABSTACKTEMP(vel, box1Dom, velIntv.size());
      // compute velocity
      PatchMappedFunc::divideVec(box1Dom,
                                 vel,
                                 Ufab,
                                 velIntv,
                                 compRho,
                                 Interval(0, velIntv.size()-1));
      // compute mapped vorticity
      PatchMappedFunc::curl2OPS(a_box,
                                vort,
                                vel,
                                a_gridMetrics.m_N[a_didx],
                                a_gridMetrics.m_J[a_didx],
                                blockDomain,
                                Interval(0, velIntv.size()-1),
                                a_gridMetrics.dxVect());
          
      // compute magnitude of vorticity
      PatchMappedFunc::magnitude(a_box,
                                 phiFab,
                                 0,
                                 vort,
                                 vort.interval());
      // set the tag pointer
      phiPtr = &phiFab;
      break;
    }

//--Hydrocarbon flame edge

    case ValueTypeHydrocarbonFlame:
    {
      useComp = 0;
      CH_assert(m_cOH >= 0 && m_cCH2O >= 0);
      MD_BOXLOOP(a_box, i)
        {
          const Real rho = Ufab[MD_IX(i, compRho)];
          const Real mfOH   = Ufab[MD_IX(i, m_cOH  )]/rho;
          const Real mfCH2O = Ufab[MD_IX(i, m_cCH2O)]/rho;
          phiFab[MD_IX(i, 0)] = mfOH*mfCH2O;
        }
      phiPtr = &phiFab;
      break;
    }
    }

  // Tag values
  const FArrayBox& phi = *phiPtr;
  // For when we don't have boxes to restrict refinement
  if (m_restrictBox.size() == 0)
    {
      // contained inside threshold
      if (m_insideThreshold)
        {
          MD_BOXLOOP(a_box, i)
            {
              const Real phiVal = phi[MD_IX(i, useComp)];
              if (phiVal >= m_loThreshold && phiVal <= m_upThreshold)
                {
                  a_tags |= MD_GETIV(i);
                }
            }
        }
      // contained outside threshold
      else
        {
          MD_BOXLOOP(a_box, i)
            {
              const Real phiVal = phi[MD_IX(i, useComp)];
              if (phiVal <= m_loThreshold || phiVal >= m_upThreshold)
                {
                  a_tags |= MD_GETIV(i);
                }
              // Some ongoing debugging for BluffBodyC3H8Air.inputs
              // if (a_box.contains(IntVect_zero))
              //   {
              //     if (i1 == 7 || i1 == 15)
              //       {
              //         pout() << "Tagging for cell " << MD_GETIV(i) << ": "
              //                << phiVal << ' ' << m_upThreshold << std::endl;
              //       }
              //   }
            }
        }
    }
  else
    {
      // contained inside threshold
      if (m_insideThreshold)
        {
          MD_BOXLOOP(a_box, i)
            {
              const Real phiVal = phi[MD_IX(i, useComp)];
              if (phiVal >= m_loThreshold && phiVal <= m_upThreshold)
                {
                  bool nonRestricted = true;
                  for (int iv = 0; iv != m_restrictBox.size(); ++iv)
                    {
                      Box checkBox = m_restrictBox[iv];
                      checkBox.refine(CRDparam::g_refFromBase[a_level]);
                      if (checkBox.contains(MD_GETIV(i)))
                        {
                          nonRestricted = false;
                          continue;
                        }
                    }
                  // If not inside restricted boxes
                  if (nonRestricted)
                    {
                      a_tags |= MD_GETIV(i);
                    }
                }
            }
        }
      // contained outside threshold
      else
        {
          MD_BOXLOOP(a_box, i)
            {
              const Real phiVal = phi[MD_IX(i, useComp)];
              if (phiVal <= m_loThreshold || phiVal >= m_upThreshold)
                {
                  bool nonRestricted = true;
                  for (int iv = 0; iv != m_restrictBox.size(); ++iv)
                    {
                      Box checkBox = m_restrictBox[iv];
                      checkBox.refine(CRDparam::g_refFromBase[a_level]);
                      if (checkBox.contains(MD_GETIV(i)))
                        {
                          nonRestricted = false;
                          continue;
                        }
                    }
                  // If not inside restricted boxes
                  if (nonRestricted)
                    {
                      a_tags |= MD_GETIV(i);
                    }
                }
            }
        }
    }
}
