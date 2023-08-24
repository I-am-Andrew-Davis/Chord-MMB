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
 * \file CNSIBCTransientCouette.H
 *
 * \brief Initial and boundary conditions for transient Couette flow
 *
 *//*+*************************************************************************/

#ifndef _CNSIBCTRANSIENTCOUETTE_H_
#define _CNSIBCTRANSIENTCOUETTE_H_

//----- Internal -----//

#include "CNSIBCGeneralized.H"


/*******************************************************************************
 */
///  BC and initial conditions on the grid for transient Couette flow.
/**
 *//*+*************************************************************************/

class CNSIBCTransientCouette : public CNSIBCGeneralized
{
public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  /// Constructor initializes
  /*  Note that Re*nu/L is the velocity (V) according to the Reynolds number.
   *  For SpaceDim > 2, the wall moves in all tangential directions, set by
   *  Vi = sqrt(V^2/(SpaceDim-1))
   */
  CNSIBCTransientCouette();

  /// Destructor
  virtual ~CNSIBCTransientCouette();


/*==============================================================================
 * Private constructors and destructors
 *============================================================================*/

private:

//--Copy not permitted

  CNSIBCTransientCouette(const CNSIBCTransientCouette&);

//--Assignment not permitted

  CNSIBCTransientCouette &operator=(const CNSIBCTransientCouette&);


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

  virtual TagLevelFactory* setTagMethodSteadyState(const int a_tagBufferSize);

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

  virtual int steadyState(
                       FArrayBox&              a_Ux,
                       const Box&              a_box,
                       const Box&              a_disjointBox,
                       const LevelGridMetrics& a_gridMetrics,
                       const Real              a_time,
                       const int               a_level) const;

/*==============================================================================
 * Private member functions
 *============================================================================*/

private:

  /// Read any information related to the IBC from input
  void readBCInfo();

  /// Equation for turbulent couette profile approximation
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
 * Data members
 *============================================================================*/

protected:

  int m_wallNormalDir;                ///< Direction normal to the moving wall
  int m_movWall;                      ///< 0 - low wall moves, 1 - high wall
  std::string m_outputFileName;       ///< Name of file with exact solution
  Box m_loTagBox;                     ///< Box for tagging near lower wall
  Box m_hiTagBox;                     ///< Box for tagging near upper wall
  Real m_tagPct;                      ///< Percentage to tag away from boundary
  Real m_threshold;                   ///< Tag refinement threshold
  bool m_steadyState;                 ///< True: do steady state, False: do transient case
  bool m_turbulent;                   ///< Use turbulent initial conditions
};

#endif /* ! defined _CNSIBCTRANSIENTCOUETTE_H_ */
