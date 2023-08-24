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
 * \file CNSIBC.cpp
 *
 * \brief Member functions for CNSIBC
 *
 *//*+*************************************************************************/

//----- Chombo Library -----//

#include "CONSTANTS.H"
#include "ParmParse.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelGridMetrics.H"
#include "NewFourthOrderCoordSys.H"
#include "FourthOrderUtil.H"

//----- Internal -----//

#include "CNSIBCEulerAdvectionCube.H"
#include "CRDparam.H"
#include "CRDPhysics.H"
#include "DataTemp.H"
#include "ChordInput.H"
#include "TagLevelFactory.H"
#include "TagLevel.H"
#include "TagMethodBuffer.H"
#include "TagMethodGradient.H"
#include "TagMethodAnalyticAdvection.H"


/*******************************************************************************
 *
 * Class CNSIBCEulerAdvectionCube: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCEulerAdvectionCube::CNSIBCEulerAdvectionCube()
  :
  CNSIBC(),
  // All invalid values which must be corrected
  m_velocity(IntVect::Zero),
  m_center(-1.*RealVect::Unit),
  m_deltaRho(0.),
  m_sizeGaussian(-1.),
  m_radMax(-1.)
{

//--Read any BC info

  readBCInfo();

//--Set BC Type

  // BC are all periodic by default
  //setAllDomainBC(CRDparam::DomainBCTypePeriodic);
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

CNSIBCEulerAdvectionCube::~CNSIBCEulerAdvectionCube()
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
CNSIBCEulerAdvectionCube::IBCName() const
{
  return "Euler advection cube";
}

/*--------------------------------------------------------------------*/
//  Write any information related to the IBC to output
/** Terminate with empty line
 *//*-----------------------------------------------------------------*/

void
CNSIBCEulerAdvectionCube::writeIBCInfo() const
{
  CRD::msg.setPrecFloatSN(4);
  CRD::msg << "Delta density\n" << m_deltaRho << CRD::var;
  CRD::msg << "Velocity\n(";
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      if (dir != 0)
        {
          CRD::msg << ',';
        }
      CRD::msg << m_velocity[dir];
    }
  CRD::msg << ')' << CRD::var;
  if (!m_useSharpProfiles)
    {
      CRD::msg << "Size of Gaussian (domain Dim.)\n" << m_sizeGaussian << CRD::var;
      CRD::msg << "Center of Gaussian (domain  Dim.)\n(";
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          if (dir != 0)
            {
              CRD::msg << ',';
            }
          CRD::msg << m_center[dir];
        }
      CRD::msg << ')' << CRD::var;
    }
  else
    {
      CRD::msg << "Ball's and Jack's profiles\n";
    }
  CRD::msg << "Maximum radius (domain Dim.)\n" << m_radMax << CRD::var;
  CRD::msg << "Tagging threshold\n" << m_threshold << CRD::var;
  CRD::msg.setFloatDefault();
  CRD::msg.newline();
}

/*--------------------------------------------------------------------*/
//  Set the tagging method if one can be associated with IBC
/** For the Gaussian profile, tags are set based on gradients of
 *  density
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
CNSIBCEulerAdvectionCube::setTagMethod(const int a_tagBufferSize)
{
  TagLevel* tagLevel = new TagLevel;
  tagLevel->FIXMEsetTagBufferSize(a_tagBufferSize);
  // Tag density gradient
  tagLevel->appendTagMethod(
    new TagMethodGradient(CRDparam::g_CRDPhysics->densityIndex(),
                          m_threshold,
                          true));
  // Add in tag buffer
  tagLevel->appendTagMethod(new TagMethodBuffer(a_tagBufferSize));
  // Analytic tagging
  //tagLevel->appendTagMethod(new TagMethodAnalyticAdvection());
  //tagLevel->appendTagMethod(new TagMethodAnalyticAdvection(0.2694206));
  // Return the factory
  return new TagLevelFactory(tagLevel);
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
CNSIBCEulerAdvectionCube::initialize(LevelData<FArrayBox>&      a_U,
                                     LevelGridMetrics&          a_gridMetrics,
                                     const LayoutData<FluxBox>& a_unitNormals,
                                     const Real                 a_time,
                                     const int                  a_level) const
{
  // Iterate over all boxes on this level and set to the exact solution
  for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
    {
      // The disjoint box
      const Box dbox = a_U.disjointBoxLayout()[dit];
 
      // Working set boxes
      Box box2Dom = grow(dbox, 2);
      //box2Dom &= blockDomain; // periodic so unnecessary 

      // set the exact solution
      exactSol(a_U[dit],
               box2Dom,
               dbox,
               a_gridMetrics,
               a_unitNormals[dit],
               dit(),
               a_time,
               a_level);
    }
}

/*--------------------------------------------------------------------*/
//  Does an exact solution exist?
/** \return             T - Exact solution is available through
 *                          'exactSol'
 *                      F - No exact solution is available
 *//*-----------------------------------------------------------------*/

bool
CNSIBCEulerAdvectionCube::haveExactSol() const
{
  return true;
}

/*--------------------------------------------------------------------*/
//  Compute the exact solution state \<U\> in the cells
/** The exact solution is obtained through the initialization function
 *  \param[out] a_Ux    Exact solution
 *  \param[in]  a_box   Cell locations to update with \<U\>.  One
 *                      assumes that this is larger than the disjoint
 *                      box so that \<JU\> can be computed.
 *  \param[in]  a_disjointBox
 *                      The disjoint box
 *  \param[in]  a_gridMetrics
 *                      Grid metrics for the level
 *  \param[in]  a_unitNormals
 *                      Face unit-normal vectors used for turbulence
 *                      variable initialization
 *  \param[in]  a_didx  Current DataIndex on the current disjoint box
 *  \param[in]  a_time  Current solution time
 *  \param[in]  a_level Grid level
 *  \return             0 - Successfully computed exact solution
 *                      1 - Exact solution is not known
 *//*-----------------------------------------------------------------*/

int 
CNSIBCEulerAdvectionCube::exactSol(FArrayBox&              a_Ux,
                                   const Box&              a_box,
                                   const Box&              a_disjointBox,
                                   const LevelGridMetrics& a_gridMetrics,
                                   const FluxBox&          a_unitNormals,
                                   const DataIndex&        a_didx,
                                   const Real              a_time,
                                   const int               a_level) const
{
  const Real gamma      = CRDparam::g_gamma;
  const Real p0         = CRDparam::g_rho*CRDparam::g_R*CRDparam::g_T;
  const Real rho0       = CRDparam::g_rho;

  // Get coordinate system and domain for the block
  const BlockCoordSys& blockCoordSys =
    *(a_gridMetrics.getCoordSys(a_disjointBox));
  const BlockDomain& blockDomain =
    a_gridMetrics.getCoordSys().problemDomain(a_disjointBox);

  // Working set boxes
  Box box1Dom = grow(a_box, 1);
  box1Dom &= blockDomain;

  // Get physical coordinates
  FABSTACKTEMP(XiFab, box1Dom, SpaceDim);  // Cartesian coordinates
  FABSTACKTEMP(XFab, box1Dom, SpaceDim);   // Physical coordinates
  this->CNSIBC::getCellCoordinates(box1Dom, XiFab, XFab, blockCoordSys);

  // Pointwise values of U (required on a_box +1 except at physical boundaries)
  CH_assert(a_Ux.box().contains(box1Dom));
  RealVect timeCenter = m_center + m_velocity*a_time;
  D_TERM(timeCenter[0] -= std::floor(timeCenter[0]);,
         timeCenter[1] -= std::floor(timeCenter[1]);,
         timeCenter[2] -= std::floor(timeCenter[2]);)
    const Real rad0sqr = m_sizeGaussian*m_sizeGaussian;
  const Real ke = 0.5*stc::dot(m_velocity, m_velocity);
  MD_BOXLOOP(box1Dom, i)
    {
      // relative point location
      const RealVect pt{D_DECL(XFab[MD_IX(i, 0)] - timeCenter[0],
                               XFab[MD_IX(i, 1)] - timeCenter[1],
                               XFab[MD_IX(i, 2)] - timeCenter[2])};
      // gaussian profile
      Real rho = rho0;
      if (!m_useSharpProfiles)
        {
          const Real radsqr = stc::dot(pt, pt);
          const Real radnrm = std::sqrt(radsqr)/m_radMax;
          Real smoo = 0.;
          if (std::fabs(radnrm) < 1.0)
            {
              smoo = std::pow(std::cos(0.5*Pi*radnrm), 6);
            }
          rho = rho0 + m_deltaRho*std::exp(-(radsqr/rad0sqr))*smoo;
        }
      // Sharp "ball and jack's" profile
      else
        {
          // small shell
          {
            const RealVect ptC = pt - RealVect{D_DECL(0.1, -0.1, 0)};
            const Real rad = std::sqrt(stc::dot(ptC, ptC));
            Real outerRad = 0.07;
            Real width = 0.03;
            if ((rad <= outerRad) && (rad >= outerRad - width))
              {
                rho = rho0 + m_deltaRho;
              }
          }
          // large shell
          {
            const RealVect ptC = pt - RealVect{D_DECL(0.1, 0.1, 0)};
            const Real rad = std::sqrt(stc::dot(ptC, ptC));
            Real outerRad = 0.10;
            Real width = 0.03;
            if ((rad <= outerRad) && (rad >= outerRad - width))
              {
                rho = rho0 + m_deltaRho;
              }
          }
          // Jack
          {
            const RealVect ptC = pt - RealVect{D_DECL(-0.11, -0.11, 0)};
            Real width = 0.03;
            Real length = 0.21;
            Real offset = 0.07;
            Real rotate = Pi/4.;
            for (const auto dir : EachDir)
              {
                RealVect line; // a unit vector defining the direction of each spike
                if (dir == 0)
                  {
                    line = RealVect{D_DECL(std::cos(rotate), std::sin(rotate), 0)};
                  }
                else if (dir == 1)
                  {
                    line = RealVect{D_DECL(-std::sin(rotate), std::cos(rotate), 0)};
                  }
                else
                  {
                    line = RealVect_basis(dir);
                  }
                // break into normal and tangent lengths
                const Real sProj = stc::dot(line, ptC);
                const Real norm = stc::mag(ptC - sProj*line);
                if ((sProj <= length-offset) && (sProj >= -offset) && (norm < width/2))
                  {
                    rho = rho0 + m_deltaRho;
                  }
              }
          }
        }
      // Store state
      a_Ux[MD_IX(i, CRDPhysics::densityIndex())] = rho;
      stc::forEachElement<SpaceDim>(
        [&, rho=rho]
        (const stc::array_size_type a_idx)
        { 
          a_Ux[MD_IX(i, CRDPhysics::velocityInterval().begin() + a_idx)] =
            m_velocity[a_idx]*rho;
        });
      a_Ux[MD_IX(i, CRDPhysics::bulkModulusIndex())] =
        p0/(gamma - 1.) + rho*ke;
    }

  // Average values of U (required on a_box except at physical boundaries)
  if (!m_useSharpProfiles)
    {
      fourthOrderAverageCell(a_Ux, blockDomain, a_box);
    }
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
CNSIBCEulerAdvectionCube::readBCInfo()
{
  CH_assert(!m_readInput);
  CRD::msg.setTerminateOnError(false);  // Disable terminate on error
  ParmParse ppIBC("ibc");

  // Many parameters are input as normalized
  const int minWidthDir = CRDparam::g_domainLength.minDir(true);
  const Real minWidth = CRDparam::g_domainLength[minWidthDir];

//--Use "ball and jacks" instead of guassian if true (default false)

  m_useSharpProfiles = false;
  ppIBC.query("sharp_profiles", m_useSharpProfiles);

//--Change in density from center of Gaussian to ambient (default 0.14)

  m_deltaRho = 0.14;
  ppIBC.query("delta_density", m_deltaRho);

//--Constant velocity field in the domain (default (1.0, 0.5, 1.0))

  m_velocity = RealVect(D_DECL(1.0, 0.5, 1.0));
  if (ppIBC.contains("velocity"))
    {
      std::vector<Real> IBCvelocity(SpaceDim);
      ppIBC.getarr("velocity", IBCvelocity, 0, SpaceDim);
      SpaceDimArray<Real, Real>::loadFromArray(m_velocity.dataPtr(),
                                               &IBCvelocity.front());
    }

//--Size of the Gaussian given as one standard deviation, on a domain from 0
//--to 1 (default 0.1)

  m_sizeGaussian = 0.1;
  ppIBC.query("size_gaussian", m_sizeGaussian);
  if (m_sizeGaussian <= 0. || m_sizeGaussian > 1.0)
    {
      CRD::msg << "Input (AdvectionCube IBC): 'size_gaussian' must be > 0.0 "
        "and < 1.0!" << CRD::error;
    }
  m_sizeGaussian *= minWidth;

//--Center of the Gaussian at time t=0, on a domain from 0 to 1
//--(default (0.5,0.5,0.5))

  m_center = 0.5*RealVect::Unit;
  if (ppIBC.contains("center"))
    {
      std::vector<Real> IBCcenter(SpaceDim);
      ppIBC.getarr("center", IBCcenter, 0, SpaceDim);
      SpaceDimArray<Real, Real>::loadFromArray(m_center.dataPtr(),
                                               &IBCcenter.front());
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          if (m_center[dir] < 0. || m_center[dir] > 1.)
            {
              CRD::msg << "Input (AdvectionCube IBC): 'center' must be >= 0.0 "
                "and < 1.0!" << CRD::error;
            }
        }
      m_center = m_center*CRDparam::g_domainLength;
    }

/*  Since a Gaussian never goes to zero, it is forced to zero at boundaries by
 *  applying smoothing.  'outer_flat_buffer' allows one to specify the size of
 *  a flat buffer near boundaries.
 */

//--Maximum radius before the profile should be made flat (through smoothing)
//--This is given as the region to keep flat near the domain extents, on a
//--domain from 0 to 1 (default 0, meaning flat right at the boundary)

  m_radMax = 0.;
  ppIBC.query("outer_flat_buffer", m_radMax);
  if (m_radMax < 0. || m_radMax > 1.)
    {
      CRD::msg << "Input (AdvectionCube IBC): 'outer_flat_buffer' must be >= "
        "0.0 and <= 1.0!" << CRD::error;
    }
  m_radMax = 0.5*(1.0 - m_radMax)*minWidth;

//--Threshold of relative density gradient for refinement

  m_threshold = 0.02;
  ppIBC.query("tag_threshold", m_threshold);
  if (m_threshold < 0.)
    {
      CRD::msg << "Input (AdvectionCube IBC): 'tag_threshold' must be > 0.0!"
               << CRD::error;
    }

  CRD::msg.setTerminateOnError(true);  // Enable terminate on error
  if (CRD::msg.numErrors() > 0)
    {
      CRD::msg << "There " << CRD::wasNumError << " reading the ibc section of "
        "the input file.  Aborting now." << CRD::error;
    }
  m_readInput = true;
}
