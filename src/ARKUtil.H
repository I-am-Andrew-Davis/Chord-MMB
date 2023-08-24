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
 * \file ARKUtil.H
 *
 * \brief Utility functions for Additive Runge Kutta
 *
 *//*+*************************************************************************/

#ifndef _ARKUTIL_H_
#define _ARKUTIL_H_

//----- STL -----//
#include <deque>
#include <queue>

//----- Chombo -----//
#include "LevelData.H"

/// Sets the constraints on components for an ARK solve
void setConstraints(LevelData<FArrayBox>& a_constraints);

/// Sets the scale on components for an ARK solve
void setScale(LevelData<FArrayBox>& a_scale,
              const LevelData<FArrayBox>& a_JU);

/// Sets test data for the Jacobian
void setJacTestData(FArrayBox& a_data);

// Writes out the jacobian to files so we can do processing with python
void writeJacobianData(const FArrayBox& a_jacobian, const FArrayBox& a_data, const Box& a_disjointBox);

// Writes out the jacobian to files so we can do processing with python
void writeJacobianData(const LevelData<FArrayBox>& a_jacobian, const LevelData<FArrayBox>& a_WcellAvg);

/*--------------------------------------------------------------------*/
///  A fixed length queue stores only the last MaxLen elements pushed onto it
/**
 *   This class implements a fixed length queue. The size of the queue is zero
 *   upon initialization, and then increases by one for every element that is
 *   "pushed" into it. If the queue is at max length, then the the next call to
 *   push will remove the oldest element to maintain the maximum length.
 *
 *   Index 0 is the oldest element in the queue.
 *
 *//*-----------------------------------------------------------------*/
template <typename T, int MaxLen, typename Container=std::deque<T>>
class FixedQueue : public std::queue<T, Container>
{
public:

  /*--------------------------------------------------------------------*/
  //  Pushes a new element onto the queue
  //  If the queue is at max length, the oldes value is removed before
  //  the new value is pushed.
  /** \param[in]  a_value The value to push onto the queue
   *//*-----------------------------------------------------------------*/
  void push(const T& a_value)
  {
    if (this->size() == MaxLen)
      {
        this->c.pop_front();
      }
    std::queue<T, Container>::push(a_value);
  }

  /*--------------------------------------------------------------------*/
  //  Returns the element at a_index
  //  The oldest element is at index 0. The newest element is at index
  //  this->size()-1
  /** \param[in]  a_index Index of element to retrieve
   *  \returns The element at the provided index
   *//*-----------------------------------------------------------------*/
  const T& operator[](size_t a_index)
  {
    return this->c[a_index];
  }
};

#endif // _ARKUTIL_H_
