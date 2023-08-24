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
 * \brief Spectral operations for use by Chombo applications
 *
 * Makes use of fftw for FFTs.
 *
 * REMINDER ABOUT FFTW OUTPUT WAVENUMBERS (using integer division):
 * where i is the array index with size ni, and k is the wavenumber index
 * with 2-sided maximum nk = ni/2
 * i = [0, 1, ..., ni/2-1, ni/2, ..., ni-1]
 * k = [0, 1, ...,   nk-1,  -nk, ...,   -1] if ni is even
 * k = [0, 1, ...,     nk,  -nk, ...,   -1] if ni is odd (ni/2 truncates)
 *
 * the indicial transform is:
 * k = i - ni*((2*i)/ni), for ni even OR odd!
 * (again noting this is all integer math)
 *
 *//***************************************************************************/

#ifndef _SPECTRALFORCING_H_
#define _SPECTRALFORCING_H_

//----- Chombo -----//

#include "FArrayBox.H"
#include "LevelData.H"
#include "BaseFabMacros.H"
#include "SpectralUtil.H"

//----- External -----//

#if defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3
#include "fftw3-mpi.h"
#endif


/******************************************************************************/
/**
 * \brief Class that implements spectral forcing using fftw
 *
 *//***************************************************************************/

class SpectralForcing
{

/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

public:

  /// Construct with undefined state
  SpectralForcing()
    { }

#if defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3

  /// Full construction
  SpectralForcing(const ProblemDomain&     a_problemDomain,
                  const DisjointBoxLayout& a_boxes)
    { define(a_problemDomain, a_boxes); }

  /// Destructor
  ~SpectralForcing() = default;

  /// Copy constructor and copy assignment are deleted
  SpectralForcing(const SpectralForcing&) = delete;
  SpectralForcing& operator=(const SpectralForcing&) = delete;

  /// Move constructor and move assignment are default
  //  Beware that underlying objects may not yet have proper move semantics.
  //  Test extensively if you use these.
  SpectralForcing(SpectralForcing&&) = default;
  SpectralForcing& operator=(SpectralForcing&&) = default;

/*==============================================================================
 * Member functions
 *============================================================================*/

  /// Weak construction
  void define(const ProblemDomain&     a_problemDomain,
              const DisjointBoxLayout& a_boxes);

  /// Forward (R2C) 3D FFT
  void rfft(LevelData<FArrayBox>& a_U,
            const Interval&       a_intv);

  /// Inverse (C2R) 3D FFT
  void irfft(LevelData<FArrayBox>& a_U,
             const Interval&       a_intv);

  /// Apply an isotropic bandpass filter and Leray-Hopf Projection,
  /// which is equivalent to a Helmholtz decomposition in physical space
  void filter_and_project(const Interval& a_intv);

  /// Generate a new isotropic, solenoidal, and large-scale perturbation field
  void calcSpectralForce(const DisjointBoxLayout& a_boxes,
                         LevelData<FArrayBox>&    a_U,
                         LevelData<FArrayBox>&    a_perturb);

  /// Gnerate a random, isotropic, solenoidal velocity field and output to a_U
  void calcRandomIC(LevelData<FArrayBox>& a_U,
                    const long int        a_seed);

  /// Generate a random, isotropic, solenoidal velocity spectrum in m_fftwData
  //    uses a fixed Gammie-Ostriker spectrum for now
  void generateRandomSpectrum(const long int a_seed);

#else  /* Null operation if FFTW not available */

  void rfft(LevelData<FArrayBox>& a_U,
            const Interval&       a_intv)
    { }

  void irfft(LevelData<FArrayBox>& a_U,
             const Interval&       a_intv)
    { }

  void filter_and_project(const Interval& a_intv)
    { }

  void calcSpectralForce(const DisjointBoxLayout& a_boxes,
                         LevelData<FArrayBox>&    a_U,
                         LevelData<FArrayBox>&    a_perturb)
    { }

  void calcRandomIC(LevelData<FArrayBox>& a_U,
                    const long int        a_seed)
    { }

  void generateRandomSpectrum(const long int a_seed)
    { }

#endif

/*==============================================================================
 * Data members
 *============================================================================*/

public:

  static constexpr int c_verbosity = 1;

#if defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3

protected:

  const Interval m_intv = Interval(0, SpaceDim-1);

  /* Prefix c means applied to complex variables only and r for real only.
   * If c and r are both absent, it means the value it is the same for both.
   * Prefix l means local to this processor
   * Prefix t means this value applies to transposed form
   */
  ptrdiff_t m_c_numI;                 ///< Number of complex variables in unit
                                      ///< stride direction

  ptrdiff_t m_numJ;                   ///< Number of real/complex variables in
                                      ///< next to least-contiguous direction

  ptrdiff_t m_numK;                   ///< Number of real/complex variables in
                                      ///< least-contiguous direction

  ptrdiff_t m_l_numK;                 ///< Number of local real/complex
                                      ///< variables in least-contiguous
                                      ///< direction on this processor

  ptrdiff_t m_l_idxKbeg;              ///< Begin 'K' index on this processor for
                                      ///< local partition of array

  ptrdiff_t m_lt_numJ;                ///< Number of local real/complex
                                      ///< variables in next to least-contiguous
                                      ///< direction on this processor when in
                                      ///< transposed form

  ptrdiff_t m_lt_idxJbeg;             ///< Begin 'J' index on this processor
                                      ///< when in transposed form for local
                                      ///< partition of array

  FFTWDataFactory<FArrayBox> m_fftwDataFactory;
                                      ///< Factory using fftw_malloc

  LevelData<FArrayBox> m_fftwData;    ///< Layouts that match fftw's
                                      ///< requirements

  fftw_plan m_fftwForwardPlan;        ///< Plans do to the FFTs (forward)
  fftw_plan m_fftwBackwardPlan;       ///< Plans do to the FFTs (reverse)

  Copier m_fftwInCopier;              ///< Copy from m_boxes to fftw layout
  Copier m_fftwOutCopier;             ///< Copy from fftw layout to m_boxes

#endif  /* ! (defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3) */
};

#endif  /* ! defined _SPECTRALFORCING_H_ */
