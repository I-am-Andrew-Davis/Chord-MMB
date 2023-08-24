#ifdef CH_LANG_CC
/*
 *      _______               __
 *     / ___/ /  ___  __  ___/ /
 *    / /__/ _ \/ _ \/ _\/ _  /
 *    \___/_//_/\___/_/  \_._/
 *    Please refer to Copyright.txt, in Chord's root directory.
 */
#endif

#include <iostream>
#include <cstring>
#include <limits>

#include <cstdlib>
#include <ctime>

// #ifdef CH_MPI
// #include <mpi.h>
// #endif

#include <REAL.H>
#include <Misc.H>
#include <parstream.H>
// #include <CH_Timer.H>
#include <Vector.H>

#include "CRDparam.H"
#include "CRDmsg.H"

#include "NVector_LDFAB.H"

#include "UsingNamespace.H"

/// Prototypes:
// class PointerPool;

/*--------------------------------------------------------------------*
 * Definition of class global parameters declared in CRDparam.H and
 * CRDmsg.H
 *--------------------------------------------------------------------*/

CRDparam::CRDparamVar CRDparam::CRDP;
CRD::Msg& CRD::msg = CRD::Msg::getInstance();

int test_NVector_LDFAB_PointerPool();

void util_fill_vector(Vector<Real>& a_vec);

// Compares that the vector and array are equal
bool util_compare_vector(const Vector<Real>& a_vec,
                         const Real* const a_data,
                         const int a_prec);

void util_copy_vector(const Vector<Real>& a_vec,
                      Real* a_data);

Real util_rand_float();

/// Global variables for handling output:
static const char *pgmname = "test_NVector_LDFAB";
static const char *indent = "   ", *indent2 = "      ";
static bool verbose = true;

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // parseTestOptions(argc, argv, verbose);

  if (verbose)
    {
      pout() << indent2 << "Beginning " << pgmname << " ..." << std::endl;
    }

//--Run the tests

  int ret = test_NVector_LDFAB_PointerPool();
  if (ret == 0)
    {
      pout() << indent << pgmname << " passed all tests" << std::endl;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)"
             << std::endl;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}


void util_fill_vector(Vector<Real>& a_vec)
{
  for(int i=0; i != a_vec.size(); ++i)
    {
      a_vec[i] = util_rand_float();
    }
}

bool util_compare_vector(const Vector<Real>& a_vec,
                         const Real* const a_data,
                         const int a_prec)
{
  for(int i=0; i != a_vec.size(); ++i)
    {
      if(!Misc::compare1(a_vec[i], a_data[i], a_prec))
        {
          return false;
        }
    }
  return true;
}

Real util_rand_float()
{
  return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}

void util_copy_vector(const Vector<Real>& a_vec,
                      Real* a_data)
{
  if(a_data == nullptr)
    {
      pout() << "Error: copying to uninitialized pointer in util_copy_vector" << std::endl;
      return;
    }
  for(int i=0; i != a_vec.size(); ++i)
    {
      a_data[i] = a_vec[i];
    }
}

int test_NVector_LDFAB_PointerPool()
{
  // CH_TIMERS("test_CRDutil_Laplacian");
  const int prec = std::numeric_limits<Real>::digits10 - 1;
  int status = 0;

  std::srand(std::time(nullptr));

  Vector<Real> vec1(5);
  util_fill_vector(vec1);
  Vector<Real> vec2(7);
  util_fill_vector(vec2);
  Vector<Real> vec3(6);
  util_fill_vector(vec3);

/*--------------------------------------------------------------------*
 * Test storing
 *--------------------------------------------------------------------*/

  // CH_TIMER("PointerPool:storePointers", t1);
  // CH_START(t1);
  Real* data1 = new Real[vec1.size()];
  util_copy_vector(vec1, data1);
  Real* data2 = new Real[vec2.size()];
  util_copy_vector(vec2, data2);
  Real* data3 = new Real[vec3.size()];
  util_copy_vector(vec3, data3);

  PointerPool::addPointer(data1);
  PointerPool::addPointer(data2);
  PointerPool::addPointer(data3);
  //   CH_STOP(t1);
  {
    int lstat = 0;
    if(PointerPool::m_ptrs[0] != data1)
      ++lstat;
    if(PointerPool::m_ptrs[1] != data2)
      ++lstat;
    if(PointerPool::m_ptrs[2] != data3)
      ++lstat;

    pout() << "data1: " << data1 << " data2: " << data2 << " data3: " << data3 << std::endl;
    pout() << "ptr0: " << PointerPool::m_ptrs[0] << " ptr1: " << PointerPool::m_ptrs[1] << " ptr2: " << PointerPool::m_ptrs[2] << std::endl;

    bool equal = util_compare_vector(vec1, PointerPool::m_ptrs[0], prec);
    if(!equal)
      ++lstat;
    equal = util_compare_vector(vec2, PointerPool::m_ptrs[1], prec);
    if(!equal)
      ++lstat;
    equal = util_compare_vector(vec3, PointerPool::m_ptrs[2], prec);
    if(!equal)
      ++lstat;

    PointerPool::freePointers();

    if(PointerPool::m_ptrs[0] != nullptr)
      ++lstat;
    if(PointerPool::m_ptrs[1] != nullptr)
      ++lstat;
    if(PointerPool::m_ptrs[2] != nullptr)
      ++lstat;

  }

  return status;
}
