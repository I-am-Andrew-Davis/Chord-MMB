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
 * \file ARKUtil.cpp
 *
 * \brief Utility functions for Additive Runge Kutta
 *
 *//*+*************************************************************************/

//----- Standard -----//
#include <random>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

//----- Chombo -----//
#include "Interval.H"
#include "FArrayBox.H"

#include "CRDmsg.H"
#include "CRDparam.H"
#include "CRDPhysics.H"

// Development of constraints
//   (1) We start with Ax = b
//   (2) Next, we extend our x vector to 3 times its original length
//       by adding x_low and x_high as variables to constrain --- these
//       vectors are each the same size as the original x and are equal to
//       x_low  = x - min_x_val >= 0
//       x_high = x - max_x_val <= 0
//   (3) Now, we place the constraints on x_low and x_high
//       x_low  = 1 or 2
//       x_high = -1 or -2
//   (4) We have to appropriately define min_x_val and max_x_val
//       (a) we know that rho > 0 is always true
//       (b) we know that pressure > 0 is always true
//       (c) we know that 200 < temperature < 6000 is always true
//       (d) we think that rho is unchanged by this update (?)
//       (e) we think that rho*u is unchanged by this update (?)
//       (f) we think that rho*e is unchanged by this update (?)
//       (g) we think that, in principle, rho*c_n should be the only
//           conservative variable changed in this update (?)
//       (h) we think that temperature and pressure can be changed
//           in this update
//       (i) we think that velocity remains unchanged in this update
//       (j) we assume that we must specify a small epsilon*rho as the
//           margin of variation for rho (and similarly for other
//           conservative variables not being changed)
//       (k) we assume that the temperature bound is a strict one
//       (l) we know that pressure comes from temperature and species
//       (m) we know that rho*e comes from pressure and temperature
//   (5) We are pretty sure at this point that the constraints should be
//       (a) rho_0 - epsilon     <= rho_star     <= rho_0 + epsilon
//       (b) rho*u_0 - epsilon   <= rho*u_star   <= rho*u_0 + epsilon
//       (c) rho*e_0 - epsilon   <= rho*e_star   <= rho*e_0 + epsilon
//       (d) rho*c_n_0 - epsilon <= rho*c_n_star <= rho*c_n_0 + epsilon
//       This should take care of all constraints naturally and consistently
//
// Sponge value layout:
//    (spongeIdxBegin+DensityIdx) Density Hi
//    (..) X-velocity Hi
//    (..) Y-velocity Hi
//    ....
//    (spongeIdxBegin+DensityIdx+numComponents) Density Lo
//    (..) X-velocity Lo
//    (..) Y-velocity Lo
//
// Essentially all components Hi first, then all components Lo
void setConstraints(LevelData<FArrayBox>& a_constraints)
{
  // LevelData<FArrayBox> constraints(dbl, spongeNumComps, ghostvect);
  // int spongeHiIdxBegin = CRDparam::g_CRDPhysics->numConservative();
  // int spongeHiIdxEnd = spongeHiIdxBegin+CRDparam::g_CRDPhysics->numConservative()-1;
  // int spongeLoIdxBegin = spongeHiIdxEnd+1;
  // int spongeLoIdxEnd = spongeLoIdxBegin+CRDparam::g_CRDPhysics->numConservative()-1;
  // int consIdxSize = CRDparam::g_CRDPhysics->numConservative();
  // std::cout << "SpongeHiIdxBegin: " << spongeHiIdxBegin
  //           << " SpongeHiIdxEnd: " << spongeHiIdxEnd
  //           << " spongeLoIdxBegin: " << spongeLoIdxBegin
  //           << " spongeLoIdxEnd: " << spongeLoIdxEnd
  //           << " consIdxSize: " << consIdxSize
  //           << " spongeNumComps: " << spongeNumComps
  //           << std::endl;
  // Set each component to 0, +-1, or +-2
  // Real loConstraint = 1.;
  // Real hiConstraint = -1.;

  // Not set up for turbulence yet... sorry guys!
  // std::cout << "Setting unconstrained: start: " << densityIdx
  //           << " num: " << consIdxSize << std::endl;
  // std::cout << "Setting hiConstraint start: " << spongeHiIdxBegin
  //           << " numComps: " << consIdxSize << std::endl;
  // std::cout << "Setting loConstraint start: " << spongeLoIdxBegin
  //           << " numComps: " << consIdxSize << std::endl;
  // for(DataIterator dit=constraints.dataIterator(); dit.ok(); ++dit)
  //   {
  //     FArrayBox& fab = constraints[dit];
  //     fab.setVal(0., fab.box(), densityIdx, consIdxSize); // Conserved values
  //     fab.setVal(hiConstraint, fab.box(), spongeHiIdxBegin, consIdxSize);
  //     fab.setVal(loConstraint, fab.box(), spongeLoIdxBegin, consIdxSize);
  //   }
  // N_Vector NVEC_constraints;
  // NVEC_constraints = N_VMake_LDFAB(&constraints);
  // flag = KINSetConstraints(m_kmem, NVEC_constraints);

  // Initial sponge values are zero for non-species
  // std::cout << "Setting spongeIntervalHi Begin: " << spongeHiIdxBegin << " end: " << spongeHiIdxEnd << std::endl;
  // std::cout << "Setting spongeIntervalLo Begin: " << spongeLoIdxBegin << " end: " << spongeLoIdxEnd << std::endl;
  // int spongeHiSpecStart = spongeHiIdxBegin+specIdxBegin;
  // int spongeLoSpecStart = spongeLoIdxBegin+specIdxBegin;
  // std::cout << "spongeHiSpecStart: " << spongeHiSpecStart
  //           << " spongeLoSpecStart: " << spongeLoSpecStart << std::endl;
  // Interval spongeIntervalHiNonSpec(spongeHiIdxBegin, spongeHiSpecStart-1);
  // Interval spongeIntervalHiSpec(spongeHiSpecStart, spongeHiIdxEnd);
  // Interval spongeIntervalLoNonSpec(spongeLoIdxBegin, spongeLoSpecStart-1);
  // Interval spongeIntervalLoSpec(spongeLoSpecStart, spongeLoIdxEnd);
  // std::cout << "hinonspec.start: " << spongeIntervalHiNonSpec.begin()
  //           << " hinonspec.end: " << spongeIntervalHiNonSpec.end()
  //           << std::endl;
  // std::cout << "hispec.start: " << spongeIntervalHiSpec.begin()
  //           << " hispec.end: " << spongeIntervalHiSpec.end()
  //           << std::endl;
  // std::cout << "lononspec.start: " << spongeIntervalLoNonSpec.begin()
  //           << " lononspec.end: " << spongeIntervalLoNonSpec.end()
  //           << std::endl;
  // std::cout << "hispec.start: " << spongeIntervalLoSpec.begin()
  //           << " hispec.end: " << spongeIntervalLoSpec.end()
  //           << std::endl;
  // for(DataIterator dit=data_constrained.dataIterator(); dit.ok(); ++dit)
  //   {
  //     FArrayBox& fab = data_constrained[dit];
  //     // Set non-species constraint variables to zero
  //     fab.setVal(0., fab.box(), spongeIntervalHiNonSpec.begin(), spongeIntervalHiNonSpec.size());
  //     fab.setVal(0., fab.box(), spongeIntervalLoNonSpec.begin(), spongeIntervalLoNonSpec.size());

  //     // Now deal with the nasty species
  //     // MD_ARRAY_RESTRICT(fabmd, fab);
  //     for(int comp = spongeIntervalHiSpec.begin(); comp != spongeIntervalHiSpec.end()+1; ++comp)
  //       {
  //         int spec_comp = comp - spongeHiIdxBegin;
  //         MD_BOXLOOP(dbl[dit], cell)
  //           {
  //             fab[MD_IX(cell, comp)] = fab[MD_IX(cell, spec_comp)] - fab[MD_IX(cell, densityIdx)] - epsilon;
  //           }
  //       }
  //     for(int comp = spongeIntervalLoSpec.begin(); comp != spongeIntervalLoSpec.end()+1; ++comp)
  //       {
  //         int spec_comp = comp - spongeLoIdxBegin;
  //         MD_BOXLOOP(dbl[dit], cell)
  //           {
  //             fab[MD_IX(cell,comp)] = fab[MD_IX(cell, densityIdx)] - fab[MD_IX(cell, spec_comp)] + epsilon;
  //           }
  //       }
  //   }

}

void setScale(LevelData<FArrayBox>& a_scale, const LevelData<FArrayBox>& a_JU)
{
  CH_assert(a_scale.nComp() == a_JU.nComp());
  const DisjointBoxLayout& dbl = a_scale.getBoxes();

  // Scale components to be approximately equal in magnitude
  for(DataIterator dit=a_scale.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& scaleFab = a_scale[dit];
      const FArrayBox& JUFab = a_JU[dit];
      Box disjointBox = dbl[dit];
      // fab.setVal(1.0);
      for(int comp = 0; comp != a_scale.nComp(); ++comp)
      {
        Real min = Abs(JUFab.min(disjointBox, comp));
        Real max = Abs(JUFab.max(disjointBox, comp));
        Real absMax = Max(min,max);
        scaleFab.setVal(1./absMax, comp);
      }
    }
}


void setJacTestData(FArrayBox& a_data)
{

  D_TERM(int num_rows = a_data.box().size(0);,
         int num_cols = a_data.box().size(1);,
         (void)0;);
  
  const int primSpecStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
  const int numSpecies = CRDparam::g_numSpecies;
  const int presIndx = CRDparam::g_CRDPhysics->pressureIndex();
  const int rhoIndx = CRDparam::g_CRDPhysics->densityIndex();
  const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();


  // X is constant species mass fractions (i.e. each cell in a row has the same species mass fractions)
  // Y is constant temperature (i.e. each cell in a column has the same temperature)
  a_data.setVal(1.225,rhoIndx);
  
  a_data.setVal(1.,presIndx);

  for(int col = 0; col != num_cols; ++col)
  {
    Box subBox(IntVect(D_DECL(col,0,0)),IntVect(D_DECL(col,num_rows-1,0)));
    Real temp = 300+50*col;
    a_data.setVal(temp,subBox,tempIndx);
  }

  // Our random number generator
  std::random_device rd;
  std::mt19937 gen(rd()); // The classic Mersenne Twister!
  std::uniform_real_distribution<Real> dis(0,std::nextafter(1., std::numeric_limits<Real>::max()));

  // Set first row to all zeros (except one species)
  Box subBox(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(num_cols-1,0,0)));
  a_data.setVal(1, subBox, primSpecStart);
  for(unsigned int sp = 1; sp != numSpecies; ++sp)
  {
    int specIdx = sp+primSpecStart;
    a_data.setVal(0, subBox, specIdx);
  }

  for(int row = 1; row != num_rows; ++row)
  {
    // First generate species mass fracs that add up to one
    Vector<Real> spVals(static_cast<unsigned int>(numSpecies));
    Real sum = 0;
    for(unsigned int sp = 0; sp != numSpecies; ++sp)
    {
      Real rndVal = dis(gen);
      spVals[sp] = rndVal;
      sum += rndVal;
    }

    Box subBox(IntVect(D_DECL(0,row,0)),IntVect(D_DECL(num_cols-1,row,0)));
    // Now normalize and set
    for(unsigned int sp = 0; sp != numSpecies; ++sp)
    {
      int specIdx = sp+primSpecStart;
      spVals[sp] /= sum;
      a_data.setVal(spVals[sp], subBox, specIdx);
    }
  }
}

void writeJacobianData(const FArrayBox& a_jacobian,
                       const FArrayBox& a_data,
                       const Box&       a_disjointBox)
{
 
  std::string folder("matrixdata/");

  const int numNativePrim = CRDparam::g_CRDPhysics->numNativePrimitive();
  const int jacComps = a_jacobian.nComp();
  
  CH_assert(numNativePrim*numNativePrim == jacComps);

  Box box = a_disjointBox;
  MD_ARRAY_RESTRICT(arrData, a_data);
  MD_ARRAY_RESTRICT(arrJac, a_jacobian);
  MD_BOXLOOP(box, i)
  {
    std::stringstream matrixfilename;
    matrixfilename << folder << "matrix." << D_TERM(std::setfill('0') << std::setw(6) << i0,
             << "." << std::setfill('0') << std::setw(6) << i1,
             << "." << std::setfill('0') << std::setw(6) << i2)
             << ".dat";
    std::stringstream datafilename;
    datafilename << folder << "data." << D_TERM(std::setfill('0') << std::setw(6) << i0,
             << "." << std::setfill('0') << std::setw(6) << i1,
             << "." << std::setfill('0') << std::setw(6) << i2)
             << ".dat";

    std::ofstream matrixfile;
    matrixfile.open(matrixfilename.str());
    std::ofstream datafile;
    datafile.open(datafilename.str());

    const int numSpecies = CRDparam::g_numSpecies;
    const int tempIndx = CRDparam::g_CRDPhysics->temperatureIndex();
    const int primSpecStart = CRDparam::g_CRDPhysics->speciesPrimInterval().begin();
    datafile D_TERM(<< "XLoc\t" << i0, << "\nYLoc\t" << i1, << "\nZLoc\t" << i2) << "\nTemp\t" << arrData[MD_IX(i,tempIndx)];
    for (int sp = 0; sp != numSpecies; ++sp)
    {
      int specIndx = sp + primSpecStart;
      datafile << "\nspecies" << sp << "\t" << arrData[MD_IX(i,specIndx)];
    }

    for (int row=0; row != numNativePrim; ++row)
    {
      for(int col = 0; col != numNativePrim; ++col)
      {
        int idx = row*numNativePrim+col;
        matrixfile << arrJac[MD_IX(i,idx)];
        if(col != numNativePrim-1)
        {
          matrixfile << " ";
        }
      }
      if(row != numNativePrim)
        matrixfile << "\n";
    }
  }
}

void writeJacobianData(const LevelData<FArrayBox>& a_jacobian, const LevelData<FArrayBox>& a_WcellAvg)
{
  for(DataIterator dit=a_jacobian.dataIterator(); dit.ok(); ++dit)
  {
    const FArrayBox& jfab = a_jacobian[dit];
    const FArrayBox& wfab = a_WcellAvg[dit];
    const Box& disjointBox = a_WcellAvg.getBoxes()[dit];
    writeJacobianData(jfab, wfab, disjointBox);
  }
}
