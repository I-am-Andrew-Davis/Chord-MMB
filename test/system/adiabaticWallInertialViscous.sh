#!/usr/bin/env bash
set -euo pipefail
# set -o xtrace

### Simple test to start with

### NOTE: At some point, I would like to automatically test building the exe.
###       Then, I would like to add the capability to detect what the name is.
###       After this, it would be nice to be able to have just run this without
###       any user interaction and see some nice results (tables, plots, etc.).

### Adiabatic wall, inertial, Cartesian, 4th-order central, no limiting

### NOTE: At some point, it would be nice to be able to run a simple test case,
###       check for the smallest time-step, and then use this to specify the
###       time-step size of each of the other tests. This could help make
###       these tests more resilient against changes in the code and stability
###       regressions. We could also hard-code a known time-step size and add
###       a warning if the time-step size deviates too much from this.

### NOTE: Rather than trust that an input file is going to be unchanging (hah),
###       all of the inputs should be specified here. This will just increase
###       the long-term resiliency of this test

### make DIM=2 DEBUG=FALSE OPT=HIGH MPI=TRUE USE_CGNS=TRUE
exe1=chord2d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex

### Copy the executable to this directory and then test here
cp ../../exec/$exe1 ./

### Load the necessary modules -- assuming the modules match
module load cgns/4.1.1_openmpi_gcc-8.3.0

echo ""

### Input file
inputFile=MMS_test2.inputs

### Copy the input file to this directory
cp ../../exec/$inputFile ./

### Order-of-accuracy testing
OoA_dir=OoA
mkdir -p $OoA_dir

### Number of processors
numProc=4

### 16^2 grid

echo -e "****  16^2 mesh  ****\n\n\n"

plot_dir=plot/2D_16_I_C_4_NoLim_Wall
check_dir=checkpoint/2D_16_I_C_4_NoLim_Wall
mkdir -p $plot_dir
mkdir -p $check_dir

mpirun -np $numProc ./$exe1 \
       $inputFile \
       file.plot_interval = 100 \
       file.plot_prefix = $plot_dir/MMS.plot. \
       file.checkpoint_prefix = $check_dir/MMS.check. \
       file.plot_error = 1 \
       physics.fluid_models = source inertial viscous \
       ibc.periodicity = 1 0 1 \
       grid.num_cells = 16 16 16 \
       sol.max_step = 50 \
       sol.fixed_dt = 1.6E-6 \
       sol.cfl = 1. \
       limit.method = FourthOrderCentral \
       amr.max_grid_size = 8 \

cp pout.0 $OoA_dir/pout_16
mv pout* $plot_dir
mv time* $plot_dir

### 32^2 grid

echo -e "****  32^2 mesh  ****\n\n\n"

plot_dir=plot/2D_32_I_C_4_NoLim_Wall
check_dir=checkpoint/2D_32_I_C_4_NoLim_Wall
mkdir -p $plot_dir
mkdir -p $check_dir

mpirun -np $numProc ./$exe1 \
       $inputFile \
       file.plot_interval = 100 \
       file.plot_prefix = $plot_dir/MMS.plot. \
       file.checkpoint_prefix = $check_dir/MMS.check. \
       file.plot_error = 1 \
       physics.fluid_models = source inertial viscous  \
       ibc.periodicity = 1 0 1 \
       grid.num_cells = 32 32 32 \
       sol.max_step = 100 \
       sol.fixed_dt = 8.0E-7 \
       sol.cfl = 1. \
       limit.method = FourthOrderCentral \
       amr.max_grid_size = 16 \

cp pout.0 $OoA_dir/pout_32
mv pout* $plot_dir
mv time* $plot_dir

### 64^2 grid

echo -e "****  64^2 mesh  ****\n\n\n"

plot_dir=plot/2D_64_I_C_4_NoLim_Wall
check_dir=checkpoint/2D_64_I_C_4_NoLim_Wall
mkdir -p $plot_dir
mkdir -p $check_dir

mpirun -np $numProc ./$exe1 \
       $inputFile \
       file.plot_interval = 100 \
       file.plot_prefix = $plot_dir/MMS.plot. \
       file.checkpoint_prefix = $check_dir/MMS.check. \
       file.plot_error = 1 \
       physics.fluid_models = source inertial viscous  \
       ibc.periodicity = 1 0 1 \
       grid.num_cells = 64 64 64 \
       sol.max_step = 200 \
       sol.fixed_dt = 4.0E-7 \
       sol.cfl = 1. \
       limit.method = FourthOrderCentral \
       amr.max_grid_size = 32 \

cp pout.0 $OoA_dir/pout_64
mv pout* $plot_dir
mv time* $plot_dir

### 128^2 grid

echo -e "****  128^2 mesh  ****\n\n\n"

plot_dir=plot/2D_128_I_C_4_NoLim_Wall
check_dir=checkpoint/2D_128_I_C_4_NoLim_Wall
mkdir -p $plot_dir
mkdir -p $check_dir

mpirun -np $numProc ./$exe1 \
       $inputFile \
       file.plot_interval = 100 \
       file.plot_prefix = $plot_dir/MMS.plot. \
       file.checkpoint_prefix = $check_dir/MMS.check. \
       file.plot_error = 1 \
       physics.fluid_models = source inertial viscous  \
       ibc.periodicity = 1 0 1 \
       grid.num_cells = 128 128 128 \
       sol.max_step = 400 \
       sol.fixed_dt = 2.0E-7 \
       sol.cfl = 1. \
       limit.method = FourthOrderCentral \
       amr.max_grid_size = 64 \

cp pout.0 $OoA_dir/pout_128
mv pout* $plot_dir
mv time* $plot_dir

### NOTE: At some point, I would like to add the capability for comparing
###       HDF5 files from runs using different levels of optimization, serial
###       versus parallel, and different box sizes.

### NOTE: We also need to add automatic order-of-accuracy testing here
python3 ChordExactConvergance.py

### Clean up after everything is done

rm $inputFile # remove any input files

rm $exe1 # remove any executables
