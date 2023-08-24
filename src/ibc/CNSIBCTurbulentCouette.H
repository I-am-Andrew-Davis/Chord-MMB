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
 * \file CNSIBCTurbulentCouette.H
 *
 * \brief Initial and boundary conditions for the turbulent Couette flow case
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCTURBULENTCOUETTE_H_
#define _CNSIBCTURBULENTCOUETTE_H_

//----- Internal -----//

#include "CNSIBCGeneralized.H"


/*******************************************************************************
 */
///  Initialization and BC for the turbulent Couette flow case
/**
 *//*+*************************************************************************/

class CNSIBCTurbulentCouette : public CNSIBCGeneralized
{
public:

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  CNSIBCTurbulentCouette();

  /// Destructor
  virtual ~CNSIBCTurbulentCouette();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCTurbulentCouette(const CNSIBCTurbulentCouette&);

//--Assignment not permitted

  CNSIBCTurbulentCouette &operator=(const CNSIBCTurbulentCouette&);


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Return a name describing the IBC
  virtual const char *const IBCName() const;

  /// Set the tagging method if one can be associated with IBC
  virtual TagLevelFactory* setTagMethod(const int a_tagBufferSize);

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

  // Real solTermC4(Real yp, int n, Real f)
  //   {
  //     Real vp = std::sqrt(2/f);
  //     return ((2.5/vp)*std::log(2.0*yp)+0.5)*std::sin(n*Pi*yp);
  //   }
  struct anlySolTermC4 : public std::unary_function<Real, Real>
  {
    anlySolTermC4(const Real a_f, const int a_n)
      :
      m_n(a_n),
      m_vp(std::sqrt(2./a_f))
      { }
    Real operator()(const Real& a_yp) const
      {
        return ((2.5/m_vp)*std::log(2.0*a_yp)+0.5)*std::sin(m_n*Pi*a_yp);
      }
  private:
    int  m_n;
    const Real m_vp;
  };
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

  Box m_loTagBox;                     ///< Box for tagging near lower wall
  Box m_hiTagBox;                     ///< Box for tagging near upper wall
  Real m_tagPct;                      ///< Percentage to tag away from boundary
};

#endif /* ! defined _CNSIBCTURBULENTCOUETTE_H_ */
