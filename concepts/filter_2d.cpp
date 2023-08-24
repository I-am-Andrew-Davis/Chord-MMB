/*******************************************************************************
 g++ -o filter_2d filter_2d.cpp -lfftw3 -lm

 Test spectral filtering of a scalar in 2-D.  Gnuplot is used to view results
 but the code will still compile and run if Gnuplot is not installed.

 Explore modifying the number of grid points "n", and the number of grids points
 you want the solution to be resolved on "nC".  Set nC == n for identity
 transform.  Set nC == 1 to retain only "DC".

*******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <limits>

#include <fftw3.h>

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

/*============================================================================*/
/// Example of spectral filtering in 2-D
/**
 *//*=========================================================================*/


int main(int argc, char **argv)
{
  constexpr double pi = 4*std::atan(1.);
  constexpr int n = 8;   // Number of grid points
  constexpr int nC = 4;  // Number of grid points used to resolve (1 <= nC <= n)
  constexpr int fw = n - nC + 1;  // a.k.a filter width
  assert(fw <= n);

//--FFTW setup

  // Data allocation and reshape to 2-D
  constexpr int compstride = sizeof(fftw_complex) * n*(n/2+1);
  double *l_u = static_cast<double*>(fftw_malloc(compstride));
  fftw_complex *l_c_u = reinterpret_cast<fftw_complex*>(l_u);
  double (*l_u2d)[2*(n/2+1)] = (double (*)[2*(n/2+1)])l_u;
  fftw_complex (*l_c_u2d)[n/2+1] = (fftw_complex (*)[n/2+1])l_c_u;

  // FFTW plans
  fftw_plan forwardplan_u  = fftw_plan_dft_r2c_2d(n, n, l_u, l_c_u,
                                                  FFTW_ESTIMATE);
  fftw_plan backwardplan_u = fftw_plan_dft_c2r_2d(n, n, l_c_u, l_u,
                                                  FFTW_ESTIMATE);

  // Initialize, print original, and compute desired filtered state
  FILE *file = fopen("filter_2d_uorig.dat", "w");
  printf("original u, 2d:\n", fw);
  fprintf(file, "#original u, 2d:\n", fw);
  double *uref = static_cast<double*>(fftw_malloc(sizeof(double)*n*n));
  double (*uref2d)[n] = (double (*)[n])uref;
  double umin = std::numeric_limits<double>::max();
  double umax = std::numeric_limits<double>::min();
  constexpr double L = 2*pi;
  constexpr double dl = L/n;
  for (int i1 = 0; i1 != n; ++i1)
    {
      // const double y = (i1 + 0.5)*dl;
      const double y = (i1)*dl;  // Allows for "exact" transform
      for (int i0 = 0; i0 != n; ++i0)
        {
          // const double x = (i0 + 0.5)*dl;
          const double x = (i0)*dl;  // Allows for "exact" transform
          double uDC = 1.0;
          l_u2d[i1][i0] = uDC;
          for (int k = 1; k < (n/2+1); ++k)
            {
              // c decrease from ~1 to 1/(n/2)
              const double c = uDC*((double)(n - k - n/2 + 1))/(n/2);
              if (i1 == 0 && i0 == 0)
                {
                  std::cout << "k: " << k << ", c: " << c << std::endl;
                }
              l_u2d[i1][i0] += c*(std::cos(k*x) +
                                  std::sin(k*y) +
                                  std::cos(k*x)*std::sin(k*y));
            }
          umin = std::min(umin, l_u2d[i1][i0]);
          umax = std::max(umax, l_u2d[i1][i0]);
          if (i0 != 0)
            {
              putc(' ', stdout);
              putc(' ', file);
            }
          printf("%6.2f", l_u2d[i1][i0]);
          fprintf(file, "%15.8e", l_u2d[i1][i0]);
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

//--Do the transforms

  fftw_execute(forwardplan_u);

  printf("F(u):\n");
  for (int i1 = 0; i1 != n; ++i1)
    {
      for (int i0 = 0; i0 != (n/2+1); ++i0)
        {
          printf(" (%6.2f,%6.2f)",
                 l_c_u2d[i1][i0][0]/(n*n), l_c_u2d[i1][i0][1]/(n*n));
        }
      printf("\n");
    }

  // Filter
  // if (fw > 1) ... logic still works for 0 <= fw <= 1.
    {
      constexpr int isEven = 1 - n%2;
      for (int i1 = (n/2) - (fw/2-1); i1 <= (n/2) + (fw/2-isEven); ++i1)
        for (int i0 = 0; i0 != (n/2+1); ++i0)
          {
            l_c_u2d[i1][i0][0] = 0.;
            l_c_u2d[i1][i0][1] = 0.;
          }
      for (int i1 = 0; i1 < n; ++i1)
        for (int i0 = (n/2+1) - (fw/2); i0 != (n/2+1); ++i0)
          {
            l_c_u2d[i1][i0][0] = 0.;
            l_c_u2d[i1][i0][1] = 0.;
          }
    }

  printf("filtered F(u):\n");
  for (int i1 = 0; i1 != n; ++i1)
    {
      for (int i0 = 0; i0 != (n/2+1); ++i0)
        {
          printf(" (%6.2f,%6.2f)",
                 l_c_u2d[i1][i0][0]/(n*n), l_c_u2d[i1][i0][1]/(n*n));
        }
      printf("\n");
    }

  // Reverse
  fftw_execute(backwardplan_u);

  // Normalize
  for (int i1 = 0; i1 != n; ++i1)
    {
      for (int i0 = 0; i0 != n; ++i0)
        {
          l_u2d[i1][i0] /= (n*n);
        }
    }

  // Print
  file = fopen("filter_2d_ufilt.dat", "w");
  printf("filtered u, 2d, width=%d:\n", fw);
  fprintf(file, "#filtered u, 2d, width=%d:\n", fw);
  for (int i1 = 0; i1 != n; ++i1)
    {
      printf("%6.2f", l_u2d[i1][0]);
      fprintf(file, "%15.8e", l_u2d[i1][0]);
      for (int i0 = 1; i0 != n; ++i0)
        {
          printf(" %6.2f", l_u2d[i1][i0]);
          fprintf(file, " %15.8e", l_u2d[i1][i0]);
        }
      printf("\n");
      fprintf(file, "\n");
    }
  fclose(file);

//--Test pass or fail

  int status = 0;
  for (int i1 = 0; i1 != n; ++i1)
    {
      for (int i0 = 0; i0 != n; ++i0)
        {
          if (compare1(l_u2d[i1][i0], uref2d[i1][i0], 10))
            {
              printf("diff (%d,%d): %17.10e %17.10e\n", i1, i0,
                     l_u2d[i1][i0], uref2d[i1][i0]);
              ++status;
            }
        }
    }
  static const char* const statLbl[] = { "failed", "passed" };
  std::cout << std::left << std::setw(40) << "filter_2d"
            << statLbl[(status == 0)] << std::endl;

//--Clean up

  fftw_destroy_plan(forwardplan_u);
  fftw_destroy_plan(backwardplan_u);
  fftw_free(l_u);

//--Write the gnuplot commands and execute

  std::ofstream plot("filter_2d_plot.plt");
  plot << "#!/usr/bin/gnuplot -persist\n";
  plot << "set terminal wxt 0 size 2000,1000 enhanced title \"filter_2d\" "
  "persist\n";
  plot << "unset key\n";
  plot << "set multiplot layout 1,2\n";
  plot << "set hidden3d\n";
  plot << "set palette model RGB rgb 30,31,32\n";
  plot << "set cbrange [" << umin << ':' << umax << "]\n";
  plot << "set zrange [" << umin-0.5 << ':' << umax+0.5 << "]\n";
  plot << "splot \"filter_2d_uorig.dat\" matrix with pm3d\n";
  plot << "splot \"filter_2d_ufilt.dat\" matrix with pm3d\n";
  plot << "unset multiplot\n";
  plot.close();
  system("gnuplot filter_2d_plot.plt");
}


