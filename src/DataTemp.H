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
 * \file DataTemp.H
 *
 * \brief Abstractions for temporary data
 *
 *//*+*************************************************************************/

#ifndef _DATATEMP_H_
#define _DATATEMP_H_

//----- Standard Library -----//

#include <memory>
#include <vector>

//----- System -----//

#include <unistd.h>

//----- Chombo Library -----//

#include "FArrayBox.H"
#include "CHArray.H"
#include "UsingNamespace.H"

//----- Internal -----//

#include "CRDmsg.H"

// Defined by testing for verbose output
// #define DEBUGSTACKALLOC
#ifdef DEBUGSTACKALLOC
  #define STACKALLOCDBG(x) x
#else
  #define STACKALLOCDBG(x) (void)0
#endif

//--Forward declarations

class StackMem;


/*******************************************************************************
 */
/// Object left on stack to release attached memory on closing brace
/** \tparam T           Type of memory
 *  The default sentinel uses tls_stack
 *
 *//*+*************************************************************************/

template <typename T>
class DefaultStackSentinel
{
public:

  /// Constructor (a_data obtained using tls_stack.get())
  DefaultStackSentinel(std::size_t a_num, void *const a_data)
    :
    m_data(nullptr)
    {
      STACKALLOCDBG(CRD::msg << CRD::fv1 << "DefaultStackSentinel::constructor "
                    << this << " for addr " << m_data << CRD::end);
      if (a_data == nullptr)
        {
          CRD::msg << "DefaultStackSentinel::constructor insufficient stack "
            "memory!" << CRD::error;
        }
      // We _must_ use placement new to begin the lifetime of new objects, even
      // if they are trivial or POD.  You do not need to destroy the objects
      // unless the destructor has side effects.
      switch (a_num)
        {
        case 0:
          m_data = static_cast<T*>(a_data);
          break;
        case 1:
          m_data = new(a_data) T;
          break;
        default:
          m_data = new(a_data) T[a_num];
          // Some compilers might pad placement new array.  A workaround in
          // that case is to separately call placement new on each element.
          // That would really be miserable since most often T is POD.
          CH_assert(m_data == a_data);
        }
    }

  // All other constructors prohibited
  DefaultStackSentinel(const DefaultStackSentinel&) = delete;
  DefaultStackSentinel(DefaultStackSentinel&&) = delete;
  DefaultStackSentinel& operator=(const DefaultStackSentinel&) = delete;
  DefaultStackSentinel& operator=(DefaultStackSentinel&&) = delete;

  /// Destructor releases memory to tls_stack
  ~DefaultStackSentinel();

  /// Return the memory location
  T* data() const
    {
      return m_data;
    }

private:

  T* m_data;                          ///< Memory retrieved from tls_stack
};


/*******************************************************************************
 */
/// Maintains a stack of memory
/**
 *//*+*************************************************************************/

class StackMem
{

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  StackMem()
    :
    m_stackBeg(nullptr),
    m_stack(nullptr),
    m_space(0),
    m_maxUsed(0)
    { }

  /// Construct with stack size a_numBytes (allocation in units of pages)
  StackMem(const std::size_t a_numBytes)
    :
    m_stackBeg(nullptr),
    m_stack(nullptr),
    m_space(0),
    m_maxUsed(0)
    {
      m_stackBeg = alloc(a_numBytes, m_space);
      m_stack = m_stackBeg;
      CRD::msg << "Stack allocation on thread from " << m_stackBeg << " to "
               << static_cast<void*>(static_cast<unsigned char*>(m_stackBeg) +
                                     m_space)
               << ": ";
      CRD::msg.putByteUnits(m_space);
      CRD::msg << CRD::h2;
    }

  // All other constructors prohibited
  StackMem(const StackMem&) = delete;
  StackMem(StackMem&&) = delete;
  StackMem& operator=(const StackMem&) = delete;
  StackMem& operator=(StackMem&&) = delete;

  /// Destructor
  ~StackMem()
    {
      CRD::msg << "Stack deallocation on thread at " << m_stackBeg
               << " with maximum usage of ";
      CRD::msg.putByteUnits(m_maxUsed);
      CRD::msg << CRD::h2;
      free();
    }


/*==============================================================================
 * Public member functions
 *============================================================================*/

  /// Print state of the stack
  void print() const
    {
      CRD::msg << "StackMem valid: (" << m_stackBeg << ':' <<
        static_cast<void*>(static_cast<unsigned char*>(m_stack) + m_space)
               << ") current: " << m_stack << ", remaining bytes: " << m_space
               << CRD::end;
    }

  /// Get the page size
  static std::size_t pageSize()
    {
      return s_pageSize;
    }

  /// Weak construction
  void define(const std::size_t a_numBytes)
    {
      free();
      m_stackBeg = alloc(a_numBytes, m_space);
      m_stack = m_stackBeg;
    }

  /// Get some memory from the stack
  /** WARNING: This returns raw uninitialized memory.  You *must* use placement
   *  new before using the memory to begin the lifetime of objects, even if the
   *  objects are trivial or POD
   */
  template <typename T = unsigned char>
  void* get(const std::size_t a_num, std::size_t a_align = 1)
    {
      a_align = std::max(a_align, alignof(T));
      // Require that a_align % alignof(T) == 0
      a_align = ((a_align + alignof(T) - 1)/alignof(T))*alignof(T);
      // Note that sizeof(T) % alignof(T) == 0 always!
      std::size_t numBytes = a_num*sizeof(T);
      if (a_num == 0)  // Supported, return 1 byte with no alignment
        {
          numBytes = 1;
          a_align = 1;
        }
      // Note: m_stack and m_space are adjusted as needed for alignment
      void* addr = std::align(a_align, numBytes, m_stack, m_space);
      if (addr == nullptr) return nullptr;
      m_stack =
        static_cast<void*>(static_cast<unsigned char*>(m_stack) + numBytes);
      m_space -= numBytes;
      m_maxUsed =
        std::max(m_maxUsed,
                 static_cast<std::size_t>(
                   static_cast<unsigned char*>(m_stack) -
                   static_cast<unsigned char*>(m_stackBeg)));
      return addr;
    }

  /// Return the memory to the stack
  void release(void *const a_addr)
    {
      // We do not destroy objects here.  The using class should do that.  Also,
      // lifetime will end when storage is reused.
      if (a_addr < m_stackBeg || a_addr >= m_stack)
        {
          CRD::msg << "Error releasing stack memory.  Address " << a_addr
                   << " is not in range of stack memory [" << m_stackBeg
                   << ':' << m_stack << ")!" << CRD::error;
        }
      m_space += (static_cast<unsigned char*>(m_stack) -
                  static_cast<unsigned char*>(a_addr));
      m_stack = a_addr;
    }

/*==============================================================================
 * Protected member functions
 *============================================================================*/

  /// Allocate memory for the stack
  static void* alloc(const std::size_t a_numBytes, std::size_t& a_allocBytes)
    {
      CH_assert(a_numBytes > 0);
      long errPageSize = sysconf(_SC_PAGESIZE);
      if (errPageSize <= 0)
        {
          CRD::msg << "Error obtaining page size!" << CRD::error;
        }
      std::size_t pageSize = static_cast<std::size_t>(errPageSize);
      CH_assert(pageSize == s_pageSize);
      std::size_t numPage = (a_numBytes + pageSize - 1)/pageSize;
      a_allocBytes = numPage*pageSize;
      STACKALLOCDBG(CRD::msg << CRD::fv1 << "StackMem::alloc pageSize: "
                    << pageSize << ", numBytes " << a_allocBytes << CRD::end);
      // Not until C++17 (or can use C11 but play it safe and use posix for now)
      // return aligned_alloc(pageSize, a_allocBytes);
      void* addr = nullptr;
      int err = CH_System::memalign(&addr, pageSize, a_allocBytes);
      STACKALLOCDBG(CRD::msg << CRD::fv1 << "StackMem::alloc err: "
                    << err << ", addr " << addr << CRD::end);
      switch (err)
        {
        case 0: break;
        case EINVAL:
          CRD::msg << "StackMem::alloc invalid alignment of " << pageSize
                   << " for stack allocation!" << CRD::error;
          break;
        case ENOMEM:
          CRD::msg << "StackMem::alloc insufficient memory for stack "
            "allocation of " << a_allocBytes << " bytes!" << CRD::error;
          break;
        default:
          CRD::msg << "StackMem::alloc unknown error (this shouldn't happpen)!"
                   << CRD::error;
        }
#ifndef CHDEF_SYSTEM_HAVE_POSIXMEMALIGN
      if (!Misc::isAligned(addr, pageSize))
        {
          CRD::msg << "Stack memory not aligned to page as intended due to "
            "missing posix_memalign (run configure in top-level directory of "
            "Chombo and re-compile)" << CRD::warn;
        }
#endif
      return addr;
    }

  /// Free memory used by the stack
  void free()
    {
      if (m_stackBeg != nullptr)
        {
#ifdef CHDEF_SYSTEM_HAVE_POSIXMEMALIGN
          if (m_stack != m_stackBeg)
            {
              CRD::msg << "Error freeing stack memory.  Memory still in use!"
                       << CRD::error;
            }
          // Otherwise alignment may cause the first allocation to be offset
          // from m_stackBeg.  A way to fix this is to adjust m_stackBeg if
          // addr != m_stackBeg when m_space is full.  But I'd rather save the
          // cycles and find bugs with a configured system instead.
#endif
          std::free(m_stackBeg);
          m_stackBeg = nullptr;
          m_stack = nullptr;
          m_space = 0;
          m_maxUsed = 0;
        }
    }

/*==============================================================================
 * Data members
 *============================================================================*/

private:

  void* m_stackBeg;                   ///< Beginning of stack memory
  void* m_stack;                      ///< Current stack pointer (one past last
                                      ///< memory chunk handed out)
  std::size_t m_space;                ///< Available space in bytes
  std::size_t m_maxUsed;              ///< Maximum usage of the stack in bytes
  static std::size_t s_pageSize;      ///< Page size on system in bytes
};


/*******************************************************************************
 */
/// Allocator allows using StackMem with some (not all) STL containers
/**
 *  Expectations are that a container is created, used, and then destroyed.
 *  Creating copies or trying to move or assign them is not well tested.
 *  The memory used by the containers must be predicted in advance.  Discarding
 *  and allocating elements out of order will fragment and use up excess memory.
 *
 *  The associated macro should be used for creating containers
 *
 *  Reverting to using the heap if the allocated stack runs out of space is
 *  not implemented
 *
 *//*+*************************************************************************/

template <typename T>
class StackAllocImpl : public std::allocator<T>
{


/*==============================================================================
 * Types
 *============================================================================*/

public:

  // Propagate types from std::allocator
  using size_type       = typename std::allocator<T>::size_type;
  using difference_type = typename std::allocator<T>::difference_type;
  using pointer         = typename std::allocator<T>::pointer;
  using const_pointer   = typename std::allocator<T>::const_pointer;
  using reference       = typename std::allocator<T>::reference;
  using const_reference = typename std::allocator<T>::const_reference;
  using value_type      = typename std::allocator<T>::value_type;

  // Void pointer types
  using void_pointer       = void*;
  using const_void_pointer = const void*;

  /// Bind this allocator to another type (ArrayClassIndex propagates)
  template<typename U>
  struct rebind
  {
    using other = StackAllocImpl<U>;
  };

  // No propagation (not well tested)
  using propagate_on_container_copy_assignment = std::false_type;
  using propagate_on_container_move_assignment = std::false_type;
  using propagate_on_container_swap = std::false_type;
  using is_always_equal = std::false_type;


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

  /// Constructor
  StackAllocImpl(StackMem&       a_stackMem,
                 const size_type a_size,
                 const bool      a_useHeap = false)
    :
    m_stackMem(a_stackMem),
    m_data(nullptr),  // No data yet
    m_size(a_size),
    m_index(0),
    m_useHeap(a_useHeap)
    {
      STACKALLOCDBG(CRD::msg << CRD::fv1 << "StackAlloc::constructor (as "
                    << this << ')' << CRD::end);
    }

  /// Copy constructor
  StackAllocImpl(const StackAllocImpl& a_other)
    :
    m_stackMem(a_other.m_stackMem),
    // Reserve on copy
    m_data(static_cast<pointer>(m_stackMem.get<T>(a_other.m_size))),
    m_size(a_other.m_size),
    m_index(0),
    m_useHeap(a_other.m_useHeap)
    {
      STACKALLOCDBG(CRD::msg << CRD::fv1 << "StackAlloc::copy constructor (by "
                    << this << "): " << sizeof(T) << ' ' << m_data << CRD::end);
      if (m_data == nullptr)
        {
          CRD::msg << "StackAlloc insufficient stack memory!" << CRD::error;
        }
    }

  /// Construct from a stack allocator for another type
  template <typename U>
  StackAllocImpl(const StackAllocImpl<U>& a_other)
    :
    m_stackMem(a_other.m_stackMem),
    // You can either do this here or in allocate.  Allocate would be preferable
    // but the user may defer the call that ultimately reserves memory on the
    // stack, separating the stack access from the local variable.
    // Reserve on copy
    m_data(static_cast<pointer>(m_stackMem.get<T>(a_other.m_size))),
    m_size(a_other.m_size),
    m_index(0),
    m_useHeap(a_other.m_useHeap)
    {
      STACKALLOCDBG(CRD::msg << CRD::fv1 << "StackAlloc::copy U constructor "
                    "(by " << this << "): " << sizeof(T) << ' ' << m_data
                    << CRD::end);
      if (m_data == nullptr)
        {
          CRD::msg << "StackAlloc insufficient stack memory!" << CRD::error;
        }
    }

  /// Move constructor
  StackAllocImpl(StackAllocImpl&& a_other)
    :
    m_stackMem(a_other.m_stackMem),
    m_data(a_other.m_data),
    m_size(a_other.m_size),
    m_index(a_other.m_index),
    m_useHeap(a_other.m_useHeap)
    {
      STACKALLOCDBG(CRD::msg << CRD::fv1 << "StackAlloc::move constructor "
                    "(by " << this << "): from: " << &a_other << CRD::end);
      if (m_index > 0)
        {
          // Moves can only occur during construction in containers
          CRD::msg << "StackAlloc moving after creating objects on the stack "
            "is not allowed!" << CRD::error;
        }
      a_other.m_data = nullptr;
      a_other.m_index = 0;
    }

  // Assignment is prohibited
  StackAllocImpl& operator=(const StackAllocImpl&) = delete;
  StackAllocImpl& operator=(StackAllocImpl&&) = delete;

  /// Destructor
  ~StackAllocImpl()
    {
      STACKALLOCDBG(CRD::msg << CRD::fv1 << "StackAlloc::destructor (by "
                    << this << "): " << m_data << CRD::end);
      if (m_data != nullptr)
        {
          m_stackMem.release(m_data);
        }
    }

  /// Friend with other StackAlloc
  template <typename U>
  friend class StackAllocImpl;


/*==============================================================================
 * Member functions
 *============================================================================*/

  /// Allocate uninitialized memory (no construction)
  pointer allocate(size_type a_size, const_void_pointer a_hint = nullptr)
    {
      if (a_size <= m_size - m_index)
        {
          pointer p = m_data + m_index;
          m_index += a_size;
          STACKALLOCDBG(CRD::msg << CRD::fv1 << "StackAlloc::allocate on stack "
                        " (by " << this << ") size: " << a_size << " at: "
                        << p << CRD::end);
          return p;
        }
      else if (m_useHeap)
        {
          pointer p = std::allocator<T>::allocate(a_size, a_hint);
          STACKALLOCDBG(CRD::msg << CRD::fv1 << "StackAlloc::allocate on heap "
                        " (by " << this << ") size: " << a_size << " at: "
                        << p << CRD::end);
          return p;
        }
      else
        {
          // Not a bad alloc since we are not necessarily out of memory in
          // any way, but m_size is insufficient.
          CRD::msg << "Exceeded reservation for stack allocator" << CRD::error;
          return nullptr;
        }
    }

  /// Deallocate storage
  void deallocate(pointer a_addr, size_type a_size)
    {
      if (a_addr < m_data + m_size)
        {
          STACKALLOCDBG(CRD::msg << CRD::fv1 << "StackAlloc::deallocate on "
                        "stack (by " << this << ") size: " << a_size << " at: "
                        << a_addr << CRD::end);
          // Only if this was the last allocation can we backup the stack
          // reservation
          if (a_addr + a_size == m_data + m_index)
            {
              m_index -= a_size;
            }
          // Otherwise we cannot deallocate (not an error, that memory is
          // unusable until the stack unwinds)
        }
      else if (m_useHeap)
        {
          STACKALLOCDBG(CRD::msg << CRD::fv1 << "StackAlloc::deallocate on "
                        "heap (by " << this << ") size: " << a_size << " at: "
                        << a_addr << CRD::end);
          std::allocator<T>::deallocate(a_addr, a_size);
        }
      else
        {
          throw std::out_of_range("Invalid deallocate to stack allocator");
        }
    }


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  StackMem& m_stackMem;               ///< Stack memory object
  pointer m_data;                     ///< Our reservation on the stack
  const size_type m_size;             ///< Number of objects reserved
  size_type m_index;                  ///< Index into objects reserved on stack
  bool m_useHeap;                     ///< T : Use the heap when the stack
                                      ///<     reservation is full
                                      ///< F : Fail if stack reservation is
                                      ///<     full and need more memory
};

/// Allocators are not equal becuase they cannot manage each others resources
template <typename T1, typename T2>
inline constexpr bool
operator==(const StackAllocImpl<T1>&,
           const StackAllocImpl<T2>&) noexcept
{
  return false;
}

/// Allocators are inequal becuase they cannot manage each others resources
template <typename T1, typename T2>
inline constexpr bool
operator!=(const StackAllocImpl<T1>&,
           const StackAllocImpl<T2>&) noexcept
{
  return true;
}


/*******************************************************************************
 *
 * The stack for each thread
 *
 ******************************************************************************/

extern thread_local StackMem tls_stack;


/*******************************************************************************
 *
 * Class DefaultStackSentinel: inline member definitions
 *
 ******************************************************************************/

template <typename T>
inline
DefaultStackSentinel<T>::~DefaultStackSentinel()
{
  STACKALLOCDBG(CRD::msg << CRD::fv1 << "DefaultStackSentinel::destructor "
                << this << " for addr " << m_data << CRD::end);
  tls_stack.release(m_data);
}

/*******************************************************************************
 *
 * Macros for creating temporary objects
 *
 ******************************************************************************/

#define USE_STACK
#ifdef USE_STACK

template <typename T>
using StackAlloc = StackAllocImpl<T>;

#define ALIGNFABTOPAGE

// Get some memory of a certain type on the stack.  Do not manual release it!
// The sentinel releases the memory when it goes out of scope
#define STACKGET(name, T, num)                                          \
  DefaultStackSentinel<T> name ## _ss ((num), tls_stack.get<T>((num))); \
  T* name = name ## _ss.data()

// Defines an STL vector on the stack with a reservation
#define VECSTACKTEMP(name, T, size)                                     \
  std::vector<T, StackAlloc<T>> name(StackAlloc<T>(tls_stack, (size))); \
  name.reserve(size)

// Defines an STL vector on the stack with a size
#define VECSTACKTEMPSIZE(name, T, size)                                 \
  std::vector<T, StackAlloc<T>> name(                                   \
    size, StackAlloc<T>(tls_stack, (size)));

// Defines an STL vector on the stack an assigned value
#define VECSTACKTEMPASSIGN(name, T, size, val)                          \
  std::vector<T, StackAlloc<T>> name(                                   \
    size, (T)(val), StackAlloc<T>(tls_stack, (size)));

//--Temporary FAB

// Defines FAB on the stack
#ifdef ALIGNFABTOPAGE
#define FABSTACKTEMP(fabname, box, ncomp)                               \
  DefaultStackSentinel<Real> fabname ## _ss (                           \
    box.size().product()*(ncomp),                                       \
    tls_stack.get<Real>(box.size().product()*(ncomp), StackMem::pageSize())); \
  FArrayBox fabname(box, (ncomp), fabname ## _ss.data())
#else
#define FABSTACKTEMP(fabname, box, (ncomp))                             \
  DefaultStackSentinel<Real> fabname ## _ss (                           \
    box.size().product()*(ncomp),                                       \
    tls_stack.get<Real>(box.size().product()*(ncomp)));                 \
  FArrayBox fabname(box, (ncomp), fabname ## _ss.data())
#endif

//--Fluxbox

#ifdef ALIGNFABTOPAGE
#define FLUXBOXSTACKTEMP(flxname, box, ncomp)                           \
  D_TERM6(                                                              \
  DefaultStackSentinel<Real> flxname ## _0_ss (                         \
    surroundingNodes(box, 0).size().product()*(ncomp),                  \
    tls_stack.get<Real>(                                                \
      surroundingNodes(box, 0).size().product()*(ncomp),                \
      StackMem::pageSize()));,                                          \
  DefaultStackSentinel<Real> flxname ## _1_ss (                         \
    surroundingNodes(box, 1).size().product()*(ncomp),                  \
    tls_stack.get<Real>(surroundingNodes(box, 1).size().product()*(ncomp)));, \
  DefaultStackSentinel<Real> flxname ## _2_ss (                         \
    surroundingNodes(box, 2).size().product()*(ncomp),                  \
    tls_stack.get<Real>(surroundingNodes(box, 2).size().product()*(ncomp)));, \
  DefaultStackSentinel<Real> flxname ## _3_ss (                         \
    surroundingNodes(box, 3).size().product()*(ncomp),                  \
    tls_stack.get<Real>(surroundingNodes(box, 3).size().product()*(ncomp)));, \
  DefaultStackSentinel<Real> flxname ## _4_ss (                         \
    surroundingNodes(box, 4).size().product()*(ncomp),                  \
    tls_stack.get<Real>(surroundingNodes(box, 4).size().product()*(ncomp)));, \
  DefaultStackSentinel<Real> flxname ## _5_ss (                         \
    surroundingNodes(box, 5).size().product()*(ncomp),                  \
    tls_stack.get<Real>(surroundingNodes(box, 5).size().product()*(ncomp)));) \
  FluxBox flxname(box, (ncomp), D_DECL6(flxname ## _0_ss.data(),        \
                                        flxname ## _1_ss.data(),        \
                                        flxname ## _2_ss.data(),        \
                                        flxname ## _3_ss.data(),        \
                                        flxname ## _4_ss.data(),        \
                                        flxname ## _5_ss.data()))
#else
#define FLUXBOXSTACKTEMP(flxname, box, ncomp)                           \
  D_TERM6(                                                              \
  DefaultStackSentinel<Real> flxname ## _0_ss (                         \
    surroundingNodes(box, 0).size().product()*(ncomp),                  \
    tls_stack.get<Real>(surroundingNodes(box, 0).size().product()*(ncomp)));, \
  DefaultStackSentinel<Real> flxname ## _1_ss (                         \
    surroundingNodes(box, 1).size().product()*(ncomp),                  \
    tls_stack.get<Real>(surroundingNodes(box, 1).size().product()*(ncomp)));, \
  DefaultStackSentinel<Real> flxname ## _2_ss (                         \
    surroundingNodes(box, 2).size().product()*(ncomp),                  \
    tls_stack.get<Real>(surroundingNodes(box, 2).size().product()*(ncomp)));, \
  DefaultStackSentinel<Real> flxname ## _3_ss (                         \
    surroundingNodes(box, 3).size().product()*(ncomp),                  \
    tls_stack.get<Real>(surroundingNodes(box, 3).size().product()*(ncomp)));, \
  DefaultStackSentinel<Real> flxname ## _4_ss (                         \
    surroundingNodes(box, 4).size().product()*(ncomp),                  \
    tls_stack.get<Real>(surroundingNodes(box, 4).size().product()*(ncomp)));, \
  DefaultStackSentinel<Real> flxname ## _5_ss (                         \
    surroundingNodes(box, 5).size().product()*(ncomp),                  \
    tls_stack.get<Real>(surroundingNodes(box, 5).size().product()*(ncomp)));) \
  FluxBox flxname(box, (ncomp), D_DECL6(flxname ## _0_ss.data(),        \
                                        flxname ## _1_ss.data(),        \
                                        flxname ## _2_ss.data(),        \
                                        flxname ## _3_ss.data(),        \
                                        flxname ## _4_ss.data(),        \
                                        flxname ## _5_ss.data()))
#endif

//--CHArray

// Defines CHArray (BNN) on the stack
#define CHARRAYBNNSTACKTEMP(arrname, T, rank, config, box, dimc1, dimc0) \
  DefaultStackSentinel<T> arrname ## _ss (                              \
    box.size().product()*(dimc1)*(dimc0),                               \
    tls_stack.get<T>(box.size().product()*(dimc1)*(dimc0)));            \
  CHArray<T, rank, config> arrname;                                     \
  arrname.define(arrname ## _ss.data(), box, (dimc1), (dimc0))

// Defines CHArray (BN) on the stack
#define CHARRAYBNSTACKTEMP(arrname, T, rank, config, box, dimc0)        \
  DefaultStackSentinel<T> arrname ## _ss (                              \
    box.size().product()*(dimc0),                                       \
    tls_stack.get<T>(box.size().product()*(dimc0)));                    \
  CHArray<T, rank, config> arrname;                                     \
  arrname.define(arrname ## _ss.data(), box, (dimc0))

// Defines CHArray (NNB) on the stack
#define CHARRAYNNBSTACKTEMP(arrname, T, rank, config, dimc1, dimc0, box) \
  DefaultStackSentinel<T> arrname ## _ss (                              \
    (dimc1)*(dimc0)*box.size().product(),                               \
    tls_stack.get<T>((dimc1)*(dimc0)*box.size().product()));            \
  CHArray<T, rank, config> arrname;                                     \
  arrname.define(arrname ## _ss.data(), (dimc1), (dimc0), box)

// Defines CHArray (NB) on the stack
#define CHARRAYNBSTACKTEMP(arrname, T, rank, config, dimc0, box)        \
  DefaultStackSentinel<T> arrname ## _ss (                              \
    (dimc0)*box.size().product(),                                       \
    tls_stack.get<T>((dimc0)*box.size().product()));                    \
  CHArray<T, rank, config> arrname;                                     \
  arrname.define(arrname ## _ss.data(), (dimc0), box)

#else

template <typename T>
using StackAlloc = std::allocator<T>;  // May be different in C++17

// Defines an STL vector with a reservation
#define VECSTACKTEMP(name, T, size)                                     \
  std::vector<T> name;                                                  \
  name.reserve(size);

// Defines an STL vector with a size
#define VECSTACKTEMPSIZE(name, T, size)                                 \
  std::vector<T> name(size);

// Defines an STL vector on the stack an assigned value
#define VECSTACKTEMPASSIGN(name, T, size, val)                          \
  std::vector<T> name(size, (T)(val))

// Allocate FAB on the heap
#define FABSTACKTEMP(fabname, box, ncomp)                               \
  FArrayBox fabname(box, ncomp)

//--Fluxbox

#define FLUXBOXSTACKTEMP(flxname, box, ncomp)                           \
  FluxBox flxname(box, ncomp)

//--CHArray

// Defines CHArray (BNN) on the stack
#define CHARRAYBNNSTACKTEMP(arrname, T, rank, config, box, dimc1, dimc0) \
  CHArray<T, rank, config> arrname(box, (dimc1), (dimc0))

// Defines CHArray (BN) on the stack
#define CHARRAYBNSTACKTEMP(arrname, T, rank, config, box, dimc0)        \
  CHArray<T, rank, config> arrname(box, (dimc0))

// Defines CHArray (NNB) on the stack
#define CHARRAYNNBSTACKTEMP(arrname, T, rank, config, dimc1, dimc0, box) \
  CHArray<T, rank, config> arrname((dimc1), (dimc0), box)

// Defines CHArray (NB) on the heap
#define CHARRAYNBSTACKTEMP(arrname, T, rank, config, dimc0, box)        \
  CHArray<T, rank, config> arrname((dimc0), box)

#endif

#endif /* ! defined _DATATEMP_H_ */
