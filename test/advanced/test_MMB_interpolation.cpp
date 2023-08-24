#ifdef CH_LANG_CC
/*
 *      _______               __
 *     / ___/ /  ___  __  ___/ /
 *    / /__/ _ \/ _ \/ _\/ _  /
 *    \___/_//_/\___/_/  \_._/
 *    Please refer to Copyright.txt, in Chord's root directory.
 */
#endif

//----- Standard Library -----//

//----- System -----//

//----- Chombo Library -----//

#include "CH_Timer.H"
#include "REAL.H"
#include "parstream.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelGridMetrics.H"
#include "AMR.H"
#include "FourthOrderUtil.H"

//----- Internal -----//

#include "chord.H"
#include "CRDPhysics.H"
#include "CNSIBCFactory.H"
#include "CNSIBCEulerAdvectionCube.H"
#include "DataTemp.H"
#include "parseTestOptions.H"

#include "UsingNamespace.H"

/// Prototypes:
int
test_MMB_interpolation();

/// Global variables for handling output:
static const char *pgmname = "test_MMB_interpolation";
static const char *indent = "   ", *indent2 = "      ";
static bool verbose = true;

int
main(int argc, const char* argv[])
{
  parseTestOptions(argc, argv, verbose);

  if (verbose)
    {
      pout() << indent2 << "Beginning " << pgmname << " ..." << std::endl;
    }

//--Run the tests

  int ret = test_MMB_interpolation();
  if (ret == 0)
    {
      pout() << indent << pgmname << " passed all tests" << std::endl;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)"
             << std::endl;
    }

  return ret;
}


/*******************************************************************************
 */
/// Custom IBC
/*
 ******************************************************************************/

class CNSIBC_test_MMB_interpolation : public CNSIBCEulerAdvectionCube
{
public:
  /// Constructor initializes
  CNSIBC_test_MMB_interpolation()
    :
    CNSIBCEulerAdvectionCube()
    {
      // m_velocity = IntVect_zero;
    }

  /// Destructor
  virtual ~CNSIBC_test_MMB_interpolation()
    { }

//--Copy and assignment not permitted

  CNSIBC_test_MMB_interpolation(
    const CNSIBC_test_MMB_interpolation&) = delete;
  CNSIBC_test_MMB_interpolation &operator=(
    const CNSIBC_test_MMB_interpolation&) = delete;

  /// Exact solution is used for initialization
  virtual int exactSol(FArrayBox&        a_Ux,
                       const Box&        a_box,
                       const Box&        a_disjointBox,
                       LevelGridMetrics& a_gridMetrics,
                       const Real        a_time,
                       const int         a_level) const
    {
      constexpr int init = 1;
      if (init == 0)
        {
          constexpr Real rho = 1.;
          constexpr Real eT = 1.;
          constexpr Real vel = 0.;

          a_Ux.setVal(rho, a_box, CRDPhysics::densityIndex());
          a_Ux.setVal(vel, a_box,
                      CRDPhysics::velocityInterval().begin(),
                      CRDPhysics::velocityInterval().size());
          a_Ux.setVal(eT,  a_box, CRDPhysics::bulkModulusIndex());
        }
      else if (init == 1)
        {
          const Real gamma = CRDparam::g_gamma;
          const Real p0    = CRDparam::g_rho*CRDparam::g_R*CRDparam::g_T;
          const Real rho0  = CRDparam::g_rho;

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

          // Pointwise values of U (required on a_box +1 except at physical
          // boundaries)
          CH_assert(a_Ux.box().contains(box1Dom));
          RealVect timeCenter = m_center + m_velocity*a_time;
          D_TERM(timeCenter[0] -= std::floor(timeCenter[0]);,
                 timeCenter[1] -= std::floor(timeCenter[1]);,
                 timeCenter[2] -= std::floor(timeCenter[2]);)
          const Real ke = 0.5*stc::dot(m_velocity, m_velocity);

          // Vary the solution orthogonal to the velocity.
          if (SpaceDim == 3)
            {
              CH_assert(m_velocity[2] == 0.);
            }

          // Choose a vector along which the solution varies
          RealVect uvec = RealVect_basis(0);
          // if (stc::mag(m_velocity) > 1.E-6)
          //   {
          //     uvec = stc::unit(RealVect{m_velocity[1], -m_velocity[0], 0.});
          //   }

          MD_BOXLOOP(box1Dom, i)
            {
              const RealVect pt{D_DECL(XFab[MD_IX(i, 0)] - timeCenter[0],
                                       XFab[MD_IX(i, 1)] - timeCenter[1],
                                       XFab[MD_IX(i, 2)] - timeCenter[2])};
              Real s = stc::dot(pt, uvec);
              if (s > 0.5)
                {
                  s = 1 - s;
                }
              else if (s < -0.5)
                {
                  s = -1 - s;
                }
              constexpr Real a = 0.;
              constexpr Real b = 0.25;
              constexpr Real c = 0.;
              const Real rho = rho0 + s*(a + s*(b + c*s));

              // Store state
              a_Ux[MD_IX(i, CRDPhysics::densityIndex())] = rho;
              stc::forEachElement<SpaceDim>(
                [&, rho=rho]
                (const stc::array_size_type a_idx)
                { 
                  a_Ux[MD_IX(i, CRDPhysics::velocityInterval().begin()
                             + a_idx)] = m_velocity[a_idx]*rho;
                });
              a_Ux[MD_IX(i, CRDPhysics::bulkModulusIndex())] =
                p0/(gamma - 1.) + rho*ke;
            }

          // Average values of U (required on a_box except at physical
          // boundaries)
          fourthOrderAverageCell(a_Ux, blockDomain, a_box);
        }
      return 0;
    }
};

/// Factory for IBC
class CNSIBCFactory_test_MMB_interpolation : public CNSIBCFactory
{
public:
  virtual ~CNSIBCFactory_test_MMB_interpolation() = default;

  /// Make a CNSIBC object
  virtual CNSIBC* new_CNSIBC() const override
    {
      return new CNSIBC_test_MMB_interpolation;
    }
};

/*******************************************************************************
 */
/// Construct the fixed hierarch of boxes for the level
/*
 ******************************************************************************/

Vector<Vector<Box>> makeAMRBoxes(const AMR& a_AMR)
{
  const Vector<AMRLevel*> AMRLevels = const_cast<AMR*>(&a_AMR)->getAMRLevels();
  const int numLevel = AMRLevels.size();
  Vector<Vector<Box>> boxes(numLevel);
  for (int idxLvl = 0; idxLvl != numLevel; ++idxLvl)
    {
      const IntVect deltaToCenter =
        CRDparam::g_domainBaseSize/(2*CRDparam::g_refFromBase[idxLvl]);
      IntVect lo = CRDparam::g_domainBaseSize/2 - deltaToCenter;
      IntVect hi = CRDparam::g_domainBaseSize + deltaToCenter;
      lo *= CRDparam::g_refFromBase[idxLvl];
      hi *= CRDparam::g_refFromBase[idxLvl];
      const Box valid(lo*IntVect_unit, (hi - 1)*IntVect_unit);
      std::pair<int, const Box*> domains =
        AMRLevels[idxLvl]->blockDomainVector();
      for (int idxDom = 0, idxDom_end = domains.first; idxDom != idxDom_end;
           ++idxDom)
        {
          const Box validDomain = (domains.second)[idxDom] & valid;
          a_AMR.makeBoxesInDomain(boxes[idxLvl], validDomain);
        }
    }
  return boxes;
}


/*******************************************************************************
 */
/// Routine test_CRDutil_Laplacian
/*
 ******************************************************************************/

int
test_MMB_interpolation()
{
  CH_TIME("test_MMB_interpolation");
  int status = 0;

  // Set up command-line arguments
  const int argc = 2;
  const char* argv[] = { "chord", "test_MMB_interpolation.inputs" };

  // Custom IBC Factory
  CNSIBCFactory_test_MMB_interpolation ibcFactory1;
  CNSIBCFactory* ibcFactory = &ibcFactory1;

  // Custom fixed grid
  std::function<Vector<Vector<Box>>(const AMR&)> makeAMRBoxesFunc =
    makeAMRBoxes;

  chord(argc, argv, ibcFactory, nullptr, &makeAMRBoxesFunc);

  return status;
}
