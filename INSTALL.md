# Overview

Chord is used to solve the compressible reacting Navier-Stokes equations using solution-adaptive mesh refinement.  Chord is based on the Chombo parallel AMR library, which itself has prerequisites of HDF5 and CGNS as input/output libraries.  The build system is in Chombo and both Chord and Chombo libraries are built at the same time.

We generally recommend compiling and installing serial and parallel versions of HDF5 and CGNS from source.  These are the only prerequisites.

# Obtaining the source code

Make a top-level directory that will contain both Chord and Chombo

    >$ mkdir $HOME/chombo
    >$ cd chombo

The source code is hosted at https://cfd-repo.engr.colostate.edu/ in mercurial repositories.  To clone the source (substitute username for your own):

    >$ hg clone https://username@cfd-repo.engr.colostate.edu/hgroot/chombo/Chombo
    >$ hg clone https://username@cfd-repo.engr.colostate.edu/hgroot/chombo/Chord
    >$ cd Chombo
    >$ hg update MMB
    >$ cd ../Chord
    >$ hg update MMB
    >$ cd ..

To obtain the source codes and install scripts for HDF5 and CGNS, either visit the site https://cfd-repo.engr.colostate.edu/share/iolib_local.tar.bz2 in a web browser or from the command line (note, this file is 11MB in size):

    >$ wget https://cfd-repo.engr.colostate.edu/share/iolib_local.tar.bz2

These instructions are in $HOME/chombo/Chord/INSTALL.md

# System prerequisites

Standard build tools are required for the system, especially a compiler and MPI library.  On a Debian system, we install:

     >$ apt install build-essential liblapack-dev doxygen csh autoconf mercurial gcc g++ cpp gfortran gdb valgrind cmake cmake-curses-gui libopenmpi-dev openmpi-bin libtool automake hwloc libhwloc-dev numactl libnuma-dev zlib1g-dev 

Similar packages should be available for any distribution.  Chord/Chombo should work with any MPI implementation.  Also consider installing clang++ as the C++ compiler if that is your preference (we normally use g++).

# Installing HDF5 and CGNS

HDF5 and CGNS can be installed either to /usr/local (if you have sufficient write privileges) or to $HOME/local (otherwise).  We recommend installing both libraries from source except on supercomputers where HDF5 is probably already installed and accessible via environment-modules or a similar technology.  We do not recommend using source packages from a distribution.

## Installing to a local file system

    >$ tar -xvjf iolib_local.tar.bz2
    >$ cd iolib_local/scripts

Use an editor to modify the file install_local.sh

    >$ emacs install_local.sh

The file is currently set up for Debian10 but should work easily on any Linux OS.  Adjust 'gcc_version', 'mpiccompiler', 'mpitag', 'cctag', as desired.  Set 'installdir' to your preferred installation directory.  It is currently set to $HOME/local, i.e., a directory name 'local' in your home directory.  To install for the system, change this to /usr/local, but you will need write permissions or must run the script as root.  Make sure the install directory exists.  E.g.,

    >$ mkdir $HOME/local

The srcdir variable must point to the source libraries.  If you have been following this guide closely, they are in '$HOME/chombo/iolib_local/src'

By default, the libraries are built directly in system memory for speed.  To alter the build location, change the variable 'scratchdir' to point to e.g., /tmp.

By default, CGNS tools are built as they can be very useful for inspecting and modifying CGNS files.  However, additional dependencies are required to use graphics:

    >$ apt install libxi-dev libxmu-dev freeglut3-dev libf2c2-dev libxt-dev tk-dev

If you do not have these packages, and do not have permission to install them, you can disable the building of CGNS tools which makes the compilation much more straightforward.  Do this on line 117 of the script by setting -DCGNS_BUILD_CGNSTOOLS:BOOL= to OFF.  Note, do not build CGNS tools in parallel.

Run the script to build and install HDF5 and CGNS

    >$ ./install_local.sh

It will take some time, up to 10 min.  The script cleans up after itself, but only if successful.  If not, be sure to manually remove the build directory (given by the variable 'builddir'.  Normally this is

    >$ rm -rf /dev/shm/install_tmp

If you have issues, first try disabling the building of CGNS tools.  Check if you have szip, serial hdf5, serial cgns, parallel hdf5, and parallel hdf5.  What failed might guide you on missing requirements.  Feel free to reach out if you have difficulties.

To check the libraries, go to their location and use 'ldd' to check that all dependencies on libraries are available.  You should see something like the following:

    >$ ldd $HOME/local/cgns/4.2.0_gcc-8.3.0/lib/libcgns.so
        linux-vdso.so.1 (0x00007ffef49fd000)
        libhdf5.so.200 => /usr/local/hdf5/1.12.0_gcc-8.3.0/lib/libhdf5.so.200 (0x00007f803a71d000)
        libsz.so.2 => /usr/local/szip/2.1.1/lib/libsz.so.2 (0x00007f803a708000)
        libz.so.1 => /lib/x86_64-linux-gnu/libz.so.1 (0x00007f803a4ba000)
        libdl.so.2 => /lib/x86_64-linux-gnu/libdl.so.2 (0x00007f803a4b5000)
        libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007f803a332000)
        libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007f803a16f000)
        /lib64/ld-linux-x86-64.so.2 (0x00007f803ac18000)

    >$ ldd $HOME/local/cgns/4.2.0_openmpi_gcc-8.3.0/lib/libcgns.so
        linux-vdso.so.1 (0x00007ffd363ec000)
        libhdf5.so.200 => /usr/local/hdf5/1.12.0_gcc-8.3.0/lib/libhdf5.so.200 (0x00007f06f1257000)
        libsz.so.2 => /usr/local/szip/2.1.1/lib/libsz.so.2 (0x00007f06f1242000)
        libz.so.1 => /lib/x86_64-linux-gnu/libz.so.1 (0x00007f06f0ff4000)
        libdl.so.2 => /lib/x86_64-linux-gnu/libdl.so.2 (0x00007f06f0fef000)
        libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007f06f0e6c000)
        libmpi.so.40 => /lib/x86_64-linux-gnu/libmpi.so.40 (0x00007f06f0d61000)
        libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007f06f0ba0000)
        /lib64/ld-linux-x86-64.so.2 (0x00007f06f175b000)
        libopen-rte.so.40 => /lib/x86_64-linux-gnu/libopen-rte.so.40 (0x00007f06f0ae8000)
        libopen-pal.so.40 => /lib/x86_64-linux-gnu/libopen-pal.so.40 (0x00007f06f0a3b000)
        librt.so.1 => /lib/x86_64-linux-gnu/librt.so.1 (0x00007f06f0a31000)
        libutil.so.1 => /lib/x86_64-linux-gnu/libutil.so.1 (0x00007f06f0a2c000)
        libhwloc.so.5 => /lib/x86_64-linux-gnu/libhwloc.so.5 (0x00007f06f09e8000)
        libevent-2.1.so.6 => /lib/x86_64-linux-gnu/libevent-2.1.so.6 (0x00007f06f0792000)
        libevent_pthreads-2.1.so.6 => /lib/x86_64-linux-gnu/libevent_pthreads-2.1.so.6 (0x00007f06f058f000)
        libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007f06f056e000)
        libnuma.so.1 => /lib/x86_64-linux-gnu/libnuma.so.1 (0x00007f06f0560000)
        libltdl.so.7 => /lib/x86_64-linux-gnu/libltdl.so.7 (0x00007f06f0553000)

Note above that the hdf5 libraries are the ones we installed and not from a system location.

# Installing VisIt

VisIt is required to view the results.  VisIt is a free and open-source software available at https://visit-dav.github.io/visit-website/.  While it can be built from source, we recommend installing the prebuilt binaries found on the DOWNLOADS->RELEASES tab.  VisIt has customizations for Chombo output.  To use these, launch VisIt and then select 'Options->Host profiles and configuration setup...'  In the pop-up window, under 'Select default configuration', choose 'Chombo Users'.  This will simply add some useful macros and not significantly alter the way VisIt behaves.  Select the 'Install' button and these macros should be available on the left pane the next time you start VisIt.

# Configuring Chombo and Chord

Chombo/lib/mk/Make.defs.local is the configuration that governs all options for building Chombo and Chord.  The build is done in place and combinations of optimized/debug, serial/parallel, 2-D/3-D can co-exist without conflict.

## Hardware

However, there is an optional step to detect hardware which permits some optimizations.  Run configure on the machine you will use for computation:

    >$ cd $HOME/chombo/Chombo
    >$ ./configure

All results from hardware detection are placed in the file Chombo/lib/src/BaseTools/CH_config.H.  Note that configure does not do any software configuration and you only need to run it once per system.

## Software

Start by using one of our templates for machine configuration.

    >$ cd $HOME/chombo/Chombo/lib/mk
    >$ cp local/Make.defs.CSUlocal ./Make.defs.local

Make.defs.local is where you can select the compiler and many other options.  It should mostly be set up correctly.  However, let's double check that the locations of the HDF5 and CGNS libraries are defined correctly.  About halfway down the file, you will see lines like:
```
_szip_root=${HOME}/local/szip/2.1.1
ifeq ($(MPI),TRUE)
  _hdf_root=${HOME}/local/hdf5/1.12.0_openmpi_gcc-8.3.0
  HDFMPIINCFLAGS = -I$(_hdf_root)/include
  HDFMPILIBFLAGS = -L$(_hdf_root)/lib -lhdf5 -L$(_szip_root)/lib -lsz -lz
else
  _hdf_root=${HOME}/local/hdf5/1.12.0_gcc-8.3.0
endif
```
and
```
ifeq ($(USE_CGNS),TRUE)
  ifeq ($(MPI),TRUE)
    _cgns_root=${HOME}/local/cgns/4.2.0_openmpi_gcc-8.3.0
  else
    _cgns_root=${HOME}/local/cgns/4.2.0_gcc-8.3.0
  endif
  CGNSINCFLAGS = -I$(_cgns_root)/include
  CGNSLIBFLAGS = -L$(_cgns_root)/lib -lcgns
  XTRALDFLAGS += -Wl,-rpath,$(_cgns_root)/lib
endif
```
Make sure _hdf_root points to your parallel (if MPI=TRUE) or serial hdf5 directories that you just installed, and similar for _cgns_root.  Also check _szip_root.  The location you set should have the include and lib directories.  For example, if I type 

    >$ cd ${HOME}/local/hdf5/1.12.0_openmpi_gcc-8.3.0
    >$ ls
    bin  hdf5_config.log  include  lib  share 

# Building Chombo and Chord

All building is done from within Chord; Chombo libraries are automatically built.

## Serial build

    >$ cd ${HOME}/chombo/Chord/exec
    >$ make -j16 example

The '-j16' performs a parallel build.  We suggest using 2-to-3 times the number of physical processing cores available on your computer.  E.g., if the machine has 6 CPU cores, use -j12 or -j18.

Check if you have an executable.

    >$ ls *.ex
    chord2d.Linux.64.g++.gfortran.DEBUG.OPT.ex

Note that configuration options (compiler.DEBUG.OPT) are included in the file name.  If you do not have an executable, try running the same make command again (sometimes one has to run it twice).  Next, check that the executable is using your HDF5 and CGNS libraries:

    >$ ldd chord2d.Linux.64.g++.gfortran.DEBUG.OPT.ex
        linux-vdso.so.1 (0x00007ffe558ee000)
        libcgns.so.4.2 => /home/sguzik/local/cgns/4.2.0_gcc-8.3.0/lib/libcgns.so.4.2 (0x00007ffa50ea9000)
        libhdf5.so.200 => /home/sguzik/local/hdf5/1.12.0_gcc-8.3.0/lib/libhdf5.so.200 (0x00007ffa50a80000)
        libsz.so.2 => /home/sguzik/local/szip/2.1.1/lib/libsz.so.2 (0x00007ffa50a6b000)
        libz.so.1 => /lib/x86_64-linux-gnu/libz.so.1 (0x00007ffa5081d000)
        ...

Here, you can see that libcgns.so.* and libhdf5.so.* are pointing to the local install I did to my home directory.

Once you have an executable, try a test problem.

    >$ ./chord2d.Linux.64.g++.gfortran.DEBUG.OPT.ex AMRAdvectionCube.inputs

That should run fairly fast (about 2 min), and place plot files in the plot directory.  Open the plot files in VisIt to check the results.  Launch VisIt and then select 'File->Open file...'.  From the pop-up window, navigate to Chord/exec/plot and from the right pane, select the complete set of 'gaussianAdvec.plot.*.2d.hdf5' files.  Under 'Plots', select the icon that has overlapping circles.  Cleck the box next to 'Mesh' to enable all levels and then the Apply button followed by the Dismiss button.  Click the 'Boxes' button under Macros (lower left pane).  Also click the button labeled 'Mapping on' to displace the grid to match physical space.  Next click the 'Draw' button under 'Plots' to make those changes.

Finally, click the play button under 'Time' to see how the solution evolves.  You can also show the mesh by highlighting the 'Mesh' plot and clicking Hide/Show.

## Parallel build

Now let's attempt to run Chord in parallel.  The standard configuration options are:

    DIM=[2|3]
    OPT=[TRUE|FALSE|HIGH]
    DEBUG=[TRUE|FALSE]
    MPI=[TRUE|FALSE]

This build will be in parallel and have full optimization

    >$ make -j16 MPI=TRUE OPT=HIGH DEBUG=FALSE example

As before, check for the executable and that it is using the proper parallel HDF5 and CGNS libraries.

    >$ ls *.ex
    chord2d.Linux.64.g++.gfortran.DEBUG.OPT.ex
    chord2d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex

    >$ ldd chord2d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex
        linux-vdso.so.1 (0x00007ffeb39ac000)
        libcgns.so.4.2 => /home/sguzik/local/cgns/4.2.0_openmpi_gcc-8.3.0/lib/libcgns.so.4.2 (0x00007efcc7a56000)
        libhdf5.so.200 => /home/sguzik/local/hdf5/1.12.0_openmpi_gcc-8.3.0/lib/libhdf5.so.200 (0x00007efcc7602000)
        libsz.so.2 => /home/sguzik/local/szip/2.1.1/lib/libsz.so.2 (0x00007efcc75ed000)
        libz.so.1 => /lib/x86_64-linux-gnu/libz.so.1 (0x00007efcc739f000)
        ...

This time, we will run the double-Mach reflection problem.

    >$ mpirun -np 6 chord2d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex AMRDoubleMachReflection.inputs

You will not see any output while the parallel job is running.  Instead, each processor will be writing independent output to a pout.* file where * is the processor number.  You can open a new terminal, navigate to the 'exec' directory, and view any pout.* file.

    >$ cd $HOME/chombo/Chord/exec
    >$ tail -f pout.0

The run should take only 30 sec.  View the plot files in VisIt using the same techniques.  If curious, open the file 'AMRDoubleMachReflection.inputs' and change amr.max_level from 2 to 3.  This will produce much more detailed physics but also take much longer.