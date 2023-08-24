/*******************************************************************************
 mpicxx -o filter_2d_mpi filter_2d_mpi.cpp -lfftw3_mpi -lfftw3 -lm
 mpirun -np 2 filter_2d_mpi

 Test spectral filtering of a scalar in 2-D with MPI.  Gnuplot is used to view
 results but the code will still compile and run if Gnuplot is not installed.

 Note: this code includes the optimization of transformed output from r2c and
 transformed input to c2r.  In other words, index 0 and index 1 are swapped
 (which in multidim are index {rank-1} and index {rank-2}.  Look for t_
 prefixes for transposed indices.

 Explore modifying the number of grid points "n", and the number of grids points
 you want the solution to be resolved on "nC".  Set n == nC for identity
 transform.  Set nC == 1 to retain only "DC".

*******************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <limits>

#include <fftw3-mpi.h>

/*--------------------------------------------------------------------*/
/// Comparison with limit tol as x and y -> 0.
/** \param[in]  x       First operand
 *  \param[in]  y       Second operand
 *  \param[in]  prec    Number of base 10 significant digits to compare
 *  \return             T - Not equal
 *//*-----------------------------------------------------------------*/

template <typename T>
inline bool compare1(const T &x, const T &y, int prec)
{
  const T tol = std::pow(10., -std::abs(prec));
  return std::fabs(x - y) >
  std::min(std::fabs(x), std::fabs(y))*tol + tol;
}

/*--------------------------------------------------------------------*/
/// Gather all of state onto processor 0
/**
 *//*-----------------------------------------------------------------*/

void gather_u(const double* l_data, const int count, const int r_offset,
              MPI_Win& win)
{
  MPI_Win_fence(0, win);
  MPI_Put(l_data, count, MPI_DOUBLE, 0, r_offset, count, MPI_DOUBLE, win);
  MPI_Win_fence(MPI_MODE_NOSTORE | MPI_MODE_NOPUT, win);
}

/*============================================================================*/
/// Example of spectral filtering in 2-D with MPI
/**
 *//*=========================================================================*/

int main(int argc, char **argv)
{
  constexpr double pi = 4*std::atan(1.);
  constexpr int n = 8;   // Number of grid points
  constexpr int nC = 4;  // Number of grid points used to resolve (1 <= nC <= n)
  constexpr int fw = n - nC + 1;  // a.k.a filter width
  assert(fw <= n);

//--MPI setup

  MPI_Init(&argc, &argv);
  int numProc;
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  int procID;
  MPI_Comm_rank(MPI_COMM_WORLD, &procID);

  // Prepare location and window to gather data for printing
  double* g_u = nullptr;
  MPI_Aint datasize = 0;
  if (procID == 0)
    {
      datasize = n*2*(n/2+1)*sizeof(double);
    }
  MPI_Win win;
  MPI_Win_allocate(datasize, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD,
                   &g_u, &win);
  fftw_complex *g_c_u = reinterpret_cast<fftw_complex*>(g_u);
  double (*g_u2d)[2*(n/2+1)] = (double (*)[2*(n/2+1)])g_u;
  // Note that i_{rank-1} (of size n) and i_{rank-2} (of size n/2+1 in this 2D
  // case) are transposed in the following
  fftw_complex (*g_t_c_u2d)[n] = (fftw_complex (*)[n])g_c_u;

//--FFTW setup

  // Allocations
  ptrdiff_t l_n1, l_i1Beg;
  ptrdiff_t l_t_n0, l_t_i0Beg;
  ptrdiff_t alloc_local = fftw_mpi_local_size_2d_transposed(
    n, n/2+1, MPI_COMM_WORLD, &l_n1, &l_i1Beg, &l_t_n0, &l_t_i0Beg);
  for (int iProc = 0; iProc != numProc; ++iProc)
    {
      if (procID == iProc)
        {
          // Use this for output so I/O from processes do not get mixed.
          std::ostringstream ost;
          ost << "Proc " << procID
              << "\n         alloc: " << alloc_local
              << "\n            i1: " << l_i1Beg << ' ' << l_n1
              << "\n  transpose i0: " << l_t_i0Beg << ' ' << l_t_n0
              << std::endl;
          std::cout << ost.str();
          ost.str("");
        }
      MPI_Barrier(MPI_COMM_WORLD);
    }

  // Data allocation and reshape to 2-D
  double *l_u = static_cast<double*>(
    fftw_malloc(sizeof(fftw_complex)*alloc_local));
  fftw_complex *l_c_u = reinterpret_cast<fftw_complex*>(l_u);

  double (*l_u2d)[2*(n/2+1)] = (double (*)[2*(n/2+1)])l_u;
  fftw_complex (*l_c_u2d)[n/2+1] = (fftw_complex (*)[n/2+1])l_c_u;
  fftw_complex (*l_t_c_u2d)[n] = (fftw_complex (*)[n])l_c_u;

  // FFTW plans
  fftw_plan forwardplan_u  = fftw_mpi_plan_dft_r2c_2d(
    n, n, l_u, l_c_u, MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT);
  fftw_plan backwardplan_u = fftw_mpi_plan_dft_c2r_2d(
    n, n, l_c_u, l_u, MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);

  // Initialize
  constexpr double L = 2*pi;
  constexpr double dl = L/n;
  for (int i1 = 0; i1 != l_n1; ++i1)
    {
      // const double y = (i1 + l_i1Beg + 0.5)*dl;
      const double y = (i1 + l_i1Beg)*dl;  // Allows for "exact" transform
      for (int i0 = 0; i0 != n; ++i0)
        {
          // const double x = (i0 + 0.5)*dl;
          const double x = (i0)*dl;        // Allows for "exact" transform
          double u = 1.0;
          l_u2d[i1][i0] = u;
          for (int k = 1; k < (n/2+1); ++k)
            {
              // c decreases from ~1 to 1/(n/2)
              const double c = u*((double)(n - k - n/2 + 1))/(n/2);
              l_u2d[i1][i0] += c*(std::cos(k*x) +
                                  std::sin(k*y) +
                                  std::cos(k*x)*std::sin(k*y));
            }
        }
    }

//--If master proc, print original and compute desired filtered state

  double umin = std::numeric_limits<double>::max();
  double umax = std::numeric_limits<double>::min();
  double *uref;
  double (*uref2d)[n];
  if (procID == 0)
    {
      FILE *file = fopen("filter_2d_mpi_uorig.dat", "w");
      printf("original u, 2d, mpi:\n", fw);
      fprintf(file, "#original u, 2d, mpi:\n", fw);
      uref = static_cast<double*>(fftw_malloc(sizeof(double)*n*n));
      uref2d = (double (*)[n])uref;
      for (int i1 = 0; i1 != n; ++i1)
        {
          // const double y = (i1 + 0.5)*dl;
          const double y = (i1)*dl;
          for (int i0 = 0; i0 != n; ++i0)
            {
              // const double x = (i0 + 0.5)*dl;
              const double x = (i0)*dl;
              double uDC = 1.0;
              double upnt = uDC;
              for (int k = 1; k < (n/2+1); ++k)
                {
                  // c decrease from ~1 to 1/(n/2)
                  const double c = uDC*((double)(n - k - n/2 + 1))/(n/2);
                  if (i1 == 0 && i0 == 0)
                    {
                      std::cout << "k: " << k << ", c: " << c << std::endl;
                    }
                  upnt += c*(std::cos(k*x) +
                             std::sin(k*y) +
                             std::cos(k*x)*std::sin(k*y));
                }
              umin = std::min(umin, upnt);
              umax = std::max(umax, upnt);
              if (i0 != 0)
                {
                  putc(' ', stdout);
                  putc(' ', file);
                }
              printf("%6.2f", upnt);
              fprintf(file, "%15.8e", upnt);
              // The filtered state should reproduce this
              uref2d[i1][i0] = uDC;
              for (int k = 1; k < (n/2+1) - (fw/2); ++k)
                {
                  const double c = uDC*((double)(n - k - n/2 + 1))/(n/2);
                  uref2d[i1][i0] += c*(std::cos(k*x) +
                                       std::sin(k*y) +
                                       std::cos(k*x)*std::sin(k*y));
                }
            }
          putc('\n', stdout);
          putc('\n', file);
        }
      fclose(file);
    }

//--Do the transforms

  // Forward
  fftw_execute(forwardplan_u);

  // Print
  gather_u(l_u, l_t_n0*2*(n), l_t_i0Beg*2*(n), win);
  if (procID == 0)
    {
      printf("F(u):\n");
      // We print out the transposed array as is
      for (int i0 = 0; i0 != (n/2+1); ++i0)
        {
          printf("(%6.2f,%6.2f)",
                 g_t_c_u2d[i0][0][0]/(n*n),
                 g_t_c_u2d[i0][0][1]/(n*n));
          for (int i1 = 1; i1 != n; ++i1)
            {
              printf(" (%6.2f,%6.2f)",
                     g_t_c_u2d[i0][i1][0]/(n*n),
                     g_t_c_u2d[i0][i1][1]/(n*n));
            }
          printf("\n");
        }
    }

  // Filter
  // if (fw > 1) ... logic still works for 0 <= fw <= 1.
  {
    constexpr int isEven = 1 - n%2;
    for (int i0 = 0; i0 != l_t_n0; ++i0)
      for (int i1 = (n/2) - (fw/2-1); i1 <= (n/2) + (fw/2-isEven); ++i1)
        {
          l_t_c_u2d[i0][i1][0] = 0.;
          l_t_c_u2d[i0][i1][1] = 0.;
        }
    if (l_t_i0Beg + l_t_n0 - 1 >= (n/2+1) - (fw/2))
      {
        for (int i0 = std::max((ptrdiff_t)0, (n/2+1) - (fw/2) - l_t_i0Beg);
             i0 != l_t_n0; ++i0)
          for (int i1 = 0; i1 < n; ++i1)
            {
              l_t_c_u2d[i0][i1][0] = 0.;
              l_t_c_u2d[i0][i1][1] = 0.;
            }
      }
  }

  // Print
  gather_u(l_u, l_t_n0*2*(n), l_t_i0Beg*2*(n), win);
  if (procID == 0)
    {
      printf("filtered F(u):\n");
      // We print out the transposed array as is
      for (int i0 = 0; i0 != (n/2+1); ++i0)
        {
          printf("(%6.2f,%6.2f)",
                 g_t_c_u2d[i0][0][0]/(n*n),
                 g_t_c_u2d[i0][0][1]/(n*n));
          for (int i1 = 1; i1 != n; ++i1)
            {
              printf(" (%6.2f,%6.2f)",
                     g_t_c_u2d[i0][i1][0]/(n*n),
                     g_t_c_u2d[i0][i1][1]/(n*n));
            }
          printf("\n");
        }
    }

  // Reverse
  fftw_execute(backwardplan_u);

  // Normalize
  for (int i1 = 0; i1 != l_n1; ++i1)
    {
      for (int i0 = 0; i0 != n; ++i0)
        {
          l_u2d[i1][i0] /= (n*n);
        }
    }

  // Print
  gather_u(l_u, l_n1*2*(n/2+1), l_i1Beg*2*(n/2+1), win);
  if (procID == 0)
    {
      FILE *file = fopen("filter_2d_mpi_ufilt.dat", "w");
      printf("filtered u, 2d, mpi, width=%d:\n", fw);
      fprintf(file, "#filtered u, 2d, mpi, width=%d:\n", fw);
      for (int i1 = 0; i1 != n; ++i1)
        {
          printf("%6.2f", g_u2d[i1][0]);
          fprintf(file, "%15.8e", g_u2d[i1][0]);
          for (int i0 = 1; i0 != n; ++i0)
            {
              printf(" %6.2f", g_u2d[i1][i0]);
              fprintf(file, " %15.8e", g_u2d[i1][i0]);
            }
          printf("\n");
          fprintf(file, "\n");
        }
      fclose(file);
    }

//--Test pass or fail

  int status = 0;
  if (procID == 0)
    {
      for (int i1 = 0; i1 != n; ++i1)
        {
          for (int i0 = 0; i0 != n; ++i0)
            {
              if (compare1(g_u2d[i1][i0], uref2d[i1][i0], 10))
                {
                  printf("diff (%d,%d): %17.10e %17.10e\n", i1, i0,
                         g_u2d[i1][i0], uref2d[i1][i0]);
                  ++status;
                }
            }
        }
      static const char* const statLbl[] = { "failed", "passed" };
      std::cout << std::left << std::setw(40) << "filter_2d_mpi"
                << statLbl[(status == 0)] << std::endl;
    }

//--Clean up

  if (procID == 0)
    {
      fftw_free(uref);
    }
  fftw_destroy_plan(forwardplan_u);
  fftw_destroy_plan(backwardplan_u);
  fftw_free(l_u);
  MPI_Win_free(&win);
  MPI_Finalize();

//--Write the gnuplot commands and execute

  if (procID == 0)
    {
      std::ofstream plot("filter_2d_mpi_plot.plt");
      plot << "#!/usr/bin/gnuplot -persist\n";
      plot << "set terminal wxt 0 size 2000,1000 enhanced title "
        "\"filter_2d_mpi\" persist\n";
      plot << "unset key\n";
      plot << "set multiplot layout 1,2\n";
      plot << "set hidden3d\n";
      plot << "set palette model RGB rgb 30,31,32\n";
      plot << "set cbrange [" << umin << ':' << umax << "]\n";
      plot << "set zrange [" << umin-0.5 << ':' << umax+0.5 << "]\n";
      plot << "splot \"filter_2d_mpi_uorig.dat\" matrix with pm3d\n";
      plot << "splot \"filter_2d_mpi_ufilt.dat\" matrix with pm3d\n";
      plot << "unset multiplot\n";
      plot.close();
      system("gnuplot filter_2d_mpi_plot.plt");
    }

  return status;
}
