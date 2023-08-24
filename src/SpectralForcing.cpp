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
 * \file SpectralForcing.cpp
 *
 * \brief Member functions for SpectralForcing
 *
 *//*+*************************************************************************/

#if defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3

//----- Chombo -----//

#include "CONSTANTS.H"

//----- Internal -----//

#include "SpectralForcing.H"
#include "CRDmsg.H"
#include "CRDparam.H"
#include "CRDPhysics.H"

//----- Standard Library -----//

#include <chrono>

/*******************************************************************************
 *
 * Class SpectralForcing: member definitions
 *
 ******************************************************************************/


/*-------------------------------------------------------------------------*//**
 *  \brief Define the class (weak construction)
 *  \param[in]  a_problemDomain
 *                      Problem domain must be single-level periodic
 *                      rectangular domain
 *  \param[in]  a_boxes Disjoint box layout of the state to filter
 *//*-------------------------------------------------------------------------*/
void
SpectralForcing::define(const ProblemDomain&     a_problemDomain,
                        const DisjointBoxLayout& a_boxes)
{
  // SpectralForcing always operates on a velocity vector
  const int ncomps = SpaceDim;

  // Spectral filtering might work with 32-bit floats but has not been tested
  assert(sizeof(Real) == 8);

  // Requirements:
  // 1) Single level with rectangular domain (we do not check that a_boxes
  //    cover the domain although this is a requirement).
  // 2) Periodic domain in all directions.
  for (int dir = 0; dir != SpaceDim; ++dir)
      {
        CH_assert(a_problemDomain.isPeriodic(dir));
      }

  // fftw use C (row-ordered) storage.  We simply reverse the indices of our
  // array.

  // dims gives dimensions of the complex array.  Note modification to
  // dims[SpaceDim-1] to account for complex vars.  For the real array cast in
  // the same memory, r_n0 = 2*dims[SpaceDim-1] and others remain the same.

  // Prefix c means applied to complex variables only and r for real only.
  // If c and r are both absent, it means the value it is the same for both.
  m_c_numI = a_problemDomain.size(0)/2 + 1;
  m_numJ   = a_problemDomain.size(1);
  m_numK   = a_problemDomain.size(2);

  ptrdiff_t c_dims[SpaceDim];
  c_dims[SpaceDim - 1] = m_c_numI;
  for (int dir = 1; dir != SpaceDim; ++dir)
    {
      c_dims[SpaceDim - dir - 1] = a_problemDomain.size(dir);
    }

  // Except for dims, in the following, we refer to our normal ordering of
  // dimensions.  E.g., 0 = x-direction.  The last direction, SpaceDim - 1,
  // is denoted by K and the second last by J.
  // Remember that l_alloc is the number of complex variables that must be
  // locally allocated
  ptrdiff_t lc_alloc = fftw_mpi_local_size_many(
    SpaceDim, c_dims, ncomps,
    FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD, &m_l_numK, &m_l_idxKbeg);

  // Swap 0 and 1 indices to get transpose shape
  {
    ptrdiff_t tmp = c_dims[0];
    c_dims[0] = c_dims[1];
    c_dims[1] = tmp;
    ptrdiff_t dummy   = fftw_mpi_local_size_many(
      SpaceDim, c_dims, ncomps,
      FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD, &m_lt_numJ, &m_lt_idxJbeg);
    (void)dummy;
  }

  // Construct fftw boxes
  std::vector<int> all_numK(numProc());
  MPI_Allgather(&m_l_numK,       1, MPI_INT,
                all_numK.data(), 1, MPI_INT, Chombo_MPI::comm);

  // Count the number of non-empty boxes
  int fftwNumBox = 0;
  for (int iProc = 0, iProc_end = numProc(); iProc != iProc_end; ++iProc)
    {
      if (all_numK[iProc] > 0) ++fftwNumBox;
    }

  Vector<int> procIDs(fftwNumBox);
  Vector<Box> fftwBoxes(fftwNumBox);

  const int r_loI = a_problemDomain.domainBox().smallEnd(0);
  const int   loK = a_problemDomain.domainBox().smallEnd(SpaceDim - 1);
  int r_maxI = r_loI;
  int idxK = loK;
  for (int iProc = 0; iProc != fftwNumBox; ++iProc)
    {
      CH_assert(all_numK[iProc] > 0);
      procIDs[iProc] = iProc;
      Box box = a_problemDomain.domainBox();
      box.setBig(0, r_loI + 2*m_c_numI - 1);
      box.setSmall(SpaceDim - 1, idxK);
      idxK += all_numK[iProc];
      box.setBig  (SpaceDim - 1, idxK - 1);
      fftwBoxes[iProc] = box;
      if (procID() == iProc)
        {
          CH_assert(box.numPts() < 2*lc_alloc);
        }
      r_maxI = std::max(r_maxI, box.bigEnd(0));
    }

  CH_assert(idxK == m_numK);
  CRD::msg << CRD::fvar << CRD::verb(c_verbosity) << "Local FFTW box\n";
  if (m_l_numK > 0)
    {
      CRD::msg << CRD::verb(c_verbosity) << fftwBoxes[procID()] << CRD::end;
    }
  else
    {
      CRD::msg << CRD::verb(c_verbosity) << "empty" << CRD::end;
    }

  // We have to make the problem domain match the fftwBoxes.  No need for
  // periodic
  Box fftwDomainBox = a_problemDomain.domainBox();
  fftwDomainBox.setBig(0, r_maxI);
  ProblemDomain fftwProblemDomain(fftwDomainBox);
  // std::cout << fftwProblemDomain << std::endl;
  DisjointBoxLayout fftw_box_layout(fftwBoxes, procIDs, fftwProblemDomain);
  fftw_box_layout.close();

  // l_alloc is used for the fabs in the LevelData and memory is allocated
  // with fftw_malloc (although our own alignment support would have been
  // sufficient)
  m_fftwDataFactory.define(lc_alloc);
  m_fftwData.define(fftw_box_layout,
                    ncomps,
                    IntVect::Zero,
                    m_fftwDataFactory);

  Real* l_data;
  if (m_l_numK > 0)
    {
      // Grab the data from the fab
      l_data = m_fftwData[m_fftwData.dataIterator()].dataPtr();
    }
  else
    {
      // Directly get the data from the factory (there is no box on this
      // processor but fftw may want us to allocate some memory anyways).
      l_data = m_fftwDataFactory.dataPtr();
    }
  fftw_complex *lc_data = reinterpret_cast<fftw_complex*>(l_data);

  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      // This is the size of the transform, both for real and complex.  Just
      // reusing c_dims.
      c_dims[SpaceDim - dir - 1] = a_problemDomain.size(dir);
    }

  using Clock = std::chrono::steady_clock;
  using TimePoint = typename Clock::time_point;
  TimePoint startTime = Clock::now();
  CRD::msg << CRD::verb(c_verbosity) << "Developing FFTW plans..." << CRD::end;

  m_fftwForwardPlan  = fftw_mpi_plan_many_dft_r2c(
    SpaceDim, c_dims, ncomps,
    FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, l_data, lc_data,
    MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);

  m_fftwBackwardPlan = fftw_mpi_plan_many_dft_c2r(
    SpaceDim, c_dims, ncomps,
    FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, lc_data, l_data,
    MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);

  CRD::msg << CRD::fvar << CRD::verb(c_verbosity) << "Time to develop plans\n"
           << std::chrono::duration_cast<std::chrono::duration<double>>(
             Clock::now() - startTime).count() << " s" << CRD::end;
  CRD::msg.newline();

  // In copier
  ProblemDomain notPeriodicDomain(a_problemDomain.domainBox());
  m_fftwInCopier.define(a_boxes, fftw_box_layout, notPeriodicDomain);

  // Out copier
  m_fftwOutCopier = m_fftwInCopier;
  m_fftwOutCopier.reverse();
}


/*-------------------------------------------------------------------------*//**
 *  \brief calculate new spectral forcing field from current solution field
 *  \param[in]  a_boxes     DisjointBoxLayout of solution data on level 0
 *  \param[in]  a_U         Conservative solution state on level 0
 *  \param[out] a_source    Unscaled spectral forcing field
 *//*-------------------------------------------------------------------------*/
void
SpectralForcing::calcSpectralForce(const DisjointBoxLayout& a_boxes,
                                   LevelData<FArrayBox>&    a_U,
                                   LevelData<FArrayBox>&    a_source)
{
  CH_TIME("SpectralForcing::calcSpectralForce");
  CRD::msg << CRD::fv1 << "SpectralForcing::calcSpectralForce" << CRD::end;

  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int velIndx = CRDparam::g_CRDPhysics->velocityInterval().begin();

  //compute raw a_source
  for (DataIterator dit = a_boxes.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = a_boxes[dit];
      FArrayBox& UFab = a_U[dit];
      FArrayBox& SFab = a_source[dit];
      MD_ARRAY_RESTRICT(arrU, UFab);
      MD_ARRAY_RESTRICT(arrS, SFab);
      MD_BOXLOOP(box, i)
        {
          const Real sqrtrho = std::sqrt(arrU[MD_IX(i, rhoIndx)]);
          D_TERM(arrS[MD_IX(i, 0)] = arrU[MD_IX(i, velIndx)]/sqrtrho;,
                 arrS[MD_IX(i, 1)] = arrU[MD_IX(i, velIndx+1)]/sqrtrho;,
                 arrS[MD_IX(i, 2)] = arrU[MD_IX(i, velIndx+2)]/sqrtrho;);
        }
    }

  // filter and project a_source
  rfft(a_source, m_intv);
  filter_and_project(CRDparam::g_spectralForcingInterval);
  irfft(a_source, m_intv);
}


/*-------------------------------------------------------------------------*//**
 *  \brief Apply forward real-to-complex FFT to data using fftw-mpi
 *  \param[in]  a_U     Real-space state vector to copy and transform
 *  \param[in]  a_intv  Component interval in state vector to copy
 *//*-------------------------------------------------------------------------*/
void
SpectralForcing::rfft(LevelData<FArrayBox>& a_U,
                      const Interval&       a_intv)
{
  CH_TIME("SpectralForcing::rfft");
  CRD::msg << CRD::fv3 << "SpectralForcing::rfft" << CRD::end;
  CH_assert(a_intv.begin() >= 0 && a_intv.end() < a_U.nComp());
  CH_assert(a_intv.size() == SpaceDim);

  // Copy SoA layout into AoS layout for fftw
  a_U.copyTo(a_intv,
             m_fftwData,
             Interval(0, SpaceDim-1),
             m_fftwInCopier,
             LDOperatorAoSIn<FArrayBox>{});

  // Forward
  fftw_execute(m_fftwForwardPlan);
}


/*-------------------------------------------------------------------------*//**
 *  \brief Apply inverse complex-to-real FFT to data using fftw-mpi
 *  \param[out]  a_U     Real-space state vector to copy into
 *  \param[in]   a_intv  Component interval in state vector to copy back into
 *//*-------------------------------------------------------------------------*/
void
SpectralForcing::irfft(LevelData<FArrayBox>& a_U,
                       const Interval&       a_intv)
{
  CH_TIME("SpectralForcing::irfft");
  CRD::msg << CRD::fv3 << "SpectralForcing::irfft" << CRD::end;
  CH_assert(a_intv.begin() >= 0 && a_intv.end() < a_U.nComp());
  CH_assert(a_intv.size() == SpaceDim);

  const Box domainBox = a_U.disjointBoxLayout().physDomain().domainBox();
  const int numTotal = domainBox.numPts();

  // Reverse
  fftw_execute(m_fftwBackwardPlan);

  // Copy AoS layout into standard SoA layout
  m_fftwData.copyTo(Interval(0, SpaceDim-1),
                    a_U,
                    a_intv,
                    m_fftwOutCopier,
                    LDOperatorAoSOut<FArrayBox>{});
}


/*-------------------------------------------------------------------------*//**
 *  \brief Filter and project the state stored in m_fftwData
 *  \param[in]  a_intv  bandpass filter interval, [kmin, kmax)
 *//*-------------------------------------------------------------------------*/
void
SpectralForcing::filter_and_project(const Interval& a_intv)
{
  CH_TIME("SpectralForcing::filter");
  CRD::msg << CRD::fv3 << "SpectralForcing::filter" << CRD::end;

  // CH_assert that kmin and kmax are acceptable
  CH_assert(a_intv.begin() >= 1 && a_intv.end() > a_intv.begin());

  // using kSq avoids the cost of sqrt inside the loops
  const int kSqMin = a_intv.begin() * a_intv.begin();
  const int kSqMax = a_intv.end() * a_intv.end();

  // These should be stored as an array and created in the
  // SpectralForcing::define method or during problem setup or something.
  // But, for now, as long as we're doing basic HIT, they'll always be unity.
  const double Scale0 = 1.0;
  const double Scale1 = 1.0;
  const double Scale2 = 1.0;

  // Only processes with data participate in the filtering
  if (m_l_numK > 0)
  {
    // You cannot use the BaseFab directly since the components have unit
    // stride. This is a transposed array.  In 3-D, the c-storage-order
    // complex-variable dimensions are [ltc_numJ][numK][c_numI][nC]
    const int c_numI = m_c_numI;  // Must be on stack for VLA-reshape
    const int numK   = m_numK;
    const int numJ   = m_numJ;
    const int nC     = SpaceDim;

    // Reshape data into transposed form
    FArrayBox& fftwFab = m_fftwData[m_fftwData.dataIterator()];
    auto ltc_data3d = (Real (*)[numK][c_numI][nC][2])fftwFab.dataPtr();
    // Assert this memory is used
    assert(&(fftwFab(fftwFab.smallEnd(), 0)) == &(ltc_data3d[0][0][0][0][0]));

    for (int i1 = 0; i1 != m_lt_numJ; ++i1) {
      for (int i2 = 0; i2 != numK; ++i2) {
        for (int i0 = 0; i0 != c_numI; ++i0)
        {
          const int k1i = i1 + m_lt_idxJbeg;
          const double k1 = Scale2 * (k1i - numJ*((2*k1i)/numJ));
          const double k2 = Scale1 * (i2 - numK*((2*i2)/numK));
          const double k0 = Scale0 * i0;
          const double kmagSq = k0*k0 + k1*k1 + k2*k2;
          const int pass = (kmagSq >= kSqMin) && (kmagSq < kSqMax);

          for (int c = 0; c != nC; ++c)
          {
            ltc_data3d[i1][i2][i0][c][0] *= pass;
            ltc_data3d[i1][i2][i0][c][1] *= pass;
          }

          double dot0 = 0.0;
          double dot1 = 0.0;

          dot0 += ltc_data3d[i1][i2][i0][0][0]*k0;
          dot0 += ltc_data3d[i1][i2][i0][1][0]*k1;
          dot0 += ltc_data3d[i1][i2][i0][2][0]*k2;

          dot1 += ltc_data3d[i1][i2][i0][0][1]*k0;
          dot1 += ltc_data3d[i1][i2][i0][1][1]*k1;
          dot1 += ltc_data3d[i1][i2][i0][2][1]*k2;

          double inv_ksq = 0.0;
          if (kmagSq > 0.0)
            inv_ksq = 1.0/kmagSq;

          ltc_data3d[i1][i2][i0][0][0] -= dot0*k0*inv_ksq;
          ltc_data3d[i1][i2][i0][1][0] -= dot0*k1*inv_ksq;
          ltc_data3d[i1][i2][i0][2][0] -= dot0*k2*inv_ksq;

          ltc_data3d[i1][i2][i0][0][1] -= dot1*k0*inv_ksq;
          ltc_data3d[i1][i2][i0][1][1] -= dot1*k1*inv_ksq;
          ltc_data3d[i1][i2][i0][2][1] -= dot1*k2*inv_ksq;
        } // i0
      } // i2
    } // i1

  } // if (m_l_numK > 0)

}


/*-------------------------------------------------------------------------*//**
 *  \brief calculate new random initial condition with specified spectrum
 *  \param[in]  a_U         Conservative solution state on level 0
 *//*-------------------------------------------------------------------------*/
void
SpectralForcing::calcRandomIC(LevelData<FArrayBox>& a_U,
                              const long int        a_seed)
{
  CH_TIME("SpectralForcing::calcRandomIC");
  CRD::msg << CRD::fv1 << "SpectralForcing::calcRandomIC" << CRD::end;

  generateRandomSpectrum(a_seed);
  irfft(a_U, CRDparam::g_CRDPhysics->velocityInterval());
}


void
SpectralForcing::generateRandomSpectrum(const long int a_seed)
{
  CH_TIME("SpectralForcing::generateRandomSpectrum");
  CRD::msg << CRD::fv3 << "SpectralForcing::generateRandomSpectrum" << CRD::end;

  // These should be stored as an array and created in the
  // SpectralForcing::define method or during problem setup or something.
  // But, for now, as long as we're doing basic HIT, they'll always be unity.
  const double Scale0 = 1.0;
  const double Scale1 = 1.0;
  const double Scale2 = 1.0;

  std::srand(a_seed);

  // Only MPI tasks with data participate
  if (m_l_numK > 0)
  {
    // You cannot use the BaseFab directly since the components have unit
    // stride. This is a transposed array.  In 3-D, the c-storage-order
    // complex-variable dimensions are [ltc_numJ][numK][c_numI][nC]
    const int c_numI = m_c_numI;  // Must be on stack for VLA-reshape
    const int numK   = m_numK;
    const int numJ   = m_numJ;
    const int nC     = SpaceDim;

    const double numI = static_cast<double>(c_numI-1);

    // For now just hard-code some logic regarding a Gamie-Ostriker spectrum:
    //        u_hat(k) ~ k^(expo-1) * exp(-k/kpeak),
    const double expo = -1./3.; // gives E(k) ~ k^(-2/3)
    const double kmax_iso = std::fmin(std::fmin(
                              Scale0*numI, Scale1*numJ/2.0), Scale2*numK/2.0);
    const double kpeak = kmax_iso/8.;

    // Reshape data into transposed form
    FArrayBox& fftwFab = m_fftwData[m_fftwData.dataIterator()];
    auto ltc_data3d = (Real (*)[numK][c_numI][nC][2]) fftwFab.dataPtr();

    // Assert this memory is used
    assert(&(fftwFab(fftwFab.smallEnd(), 0)) == &(ltc_data3d[0][0][0][0][0]));

    for (int i1 = 0; i1 != m_lt_numJ; ++i1) {
      for (int i2 = 0; i2 != numK; ++i2) {
        for (int i0 = 0; i0 != c_numI; ++i0)
        {
          const int k1i = i1 + m_lt_idxJbeg;
          const double k1 = Scale2 * (k1i - numJ*((2*k1i)/numJ));
          const double k2 = Scale1 * (i2 - numK*((2*i2)/numK));
          const double k0 = Scale0 * i0;
          const double kmag = std::sqrt(k0*k0 + k1*k1 + k2*k2);

          // scale-factor for hard-coded Gamie-Ostriker spectrum
          double kscale = 0.0;
          if (kmag > 0.0 && kmag < kmax_iso)
            kscale = std::pow(kmag, expo-1.0) * std::exp(-kmag/kpeak);

          // generate scaled random IC
          for (int c = 0; c != nC; ++c)
          {
            const double q0 = static_cast<double>(std::rand())/RAND_MAX;
            const double q1 = static_cast<double>(std::rand())/RAND_MAX;
            const double q2 = static_cast<double>(std::rand())/RAND_MAX;
            const double q3 = std::sqrt(-2*std::log(q0+1e-20)) * std::cos(2*PI*q1);
            ltc_data3d[i1][i2][i0][c][0] = kscale*q3*std::cos(2.0*PI*q2);
            ltc_data3d[i1][i2][i0][c][1] = kscale*q3*std::sin(2.0*PI*q2);
          }

          // solenoidally-project random IC
          double dot0 = 0.0;
          double dot1 = 0.0;

          dot0 += ltc_data3d[i1][i2][i0][0][0]*k0;
          dot0 += ltc_data3d[i1][i2][i0][1][0]*k1;
          dot0 += ltc_data3d[i1][i2][i0][2][0]*k2;

          dot1 += ltc_data3d[i1][i2][i0][0][1]*k0;
          dot1 += ltc_data3d[i1][i2][i0][1][1]*k1;
          dot1 += ltc_data3d[i1][i2][i0][2][1]*k2;

          double inv_ksq = 0.0;
          if (kmag > 0.0)
            inv_ksq = 1.0/(kmag*kmag);

          ltc_data3d[i1][i2][i0][0][0] -= dot0*k0*inv_ksq;
          ltc_data3d[i1][i2][i0][1][0] -= dot0*k1*inv_ksq;
          ltc_data3d[i1][i2][i0][2][0] -= dot0*k2*inv_ksq;

          ltc_data3d[i1][i2][i0][0][1] -= dot1*k0*inv_ksq;
          ltc_data3d[i1][i2][i0][1][1] -= dot1*k1*inv_ksq;
          ltc_data3d[i1][i2][i0][2][1] -= dot1*k2*inv_ksq;
        } // i0
      } // i2
    } // i1

  } // if (m_l_numK > 0)
}

#endif  /* ! (defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3) */
