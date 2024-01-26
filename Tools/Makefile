# Set this to the shell you use.  I need this because I have gfortran aliased in my .tcshrc
#SHELL = /bin/tcsh

# In the next two lines, specify how to invoke the fortran compilers on your machine.

ifeq ($(OSTYPE),darwin)
   fortran = gfortran
else ifeq ($(HOST),linux166.nas.nasa.gov)
   fortran = gfortran
else
   fortran = ifort
endif

$(info The variable fortran is $(fortran))

# If you make spectra then you will need the FFTW library.  Below enter the necessary paths.
fftw_include_path = $(HOME)/homebrew/opt/fftw/include
fftw-library-path = $(HOME)/homebrew/opt/fftw/lib

# The targets for this makefile are:
#    clean : Delete all objects, executables, and compiled modules.

#    partition_tool : A serial code which should be used to see whether your grid
#                sizes can be split into a given number of processors.  It will
#                ask you to enter the grid dimensions and number of processors.

#    time_average_profiles: Performs a time average of a bunch of profile files.
#                           Look at the source code time_average_profiles.f90
#                           for instructions on how to use it.
#    stretched: Tool for experimenting with parameters for creating a stretched
#               mesh.
#    spectrum: Spectrum w.r.t. of velocity fluctuations.
#    enstrophy_spectrum as a function of kr using a horizontal plane.
#    plot_Melander: Outputs data to plot the Melander C_infty function.
#    time_phi_Favre_averages.g90 : Takes phi Reynolds averages for many time instants
#                            and computes phi-time Favre averages including stresses.
#    post_processes_stresses.f90 : Take tecplot files phi_time_Favre_means.tec and phi_time_Favre_stresses
#                                  and computes the radial velocity at the midplane from the (phi,r) Favre
#                                  stress

#
# You don't have to read any further if you don't want to.

# The extension .mpio is used for object files pertaining to the mpi code.
# This is to allow both the serial and mpi codes to be built in the same
# directory.
.SUFFIXES : .f .f90 .o

# Run time check:
ifeq ($(fortran),gfortran)
   run_time_check = -fcheck=all
else
   run_time_check = -check bounds
endif

# List of possible compiler options I like to use for gfortran:
init_nan = -finit-real=nan
warn_un_init = -Wuninitialized
# Warn about array temporaries created by the compiler
warn_array_temps = -Warray-temporaries
# This has the effect of putting implicit none at the start of every subroutine
# in case I forget.
implicit_none = -fimplicit-none

optimization_level = 3

# For sanity checks, the code checks for uninitialized variables by comparing against NaN, so do
# not remove the init_nan or init=snan option below:

# Fortran options I always use:
ifeq ($(fortran), gfortran)
   fortran_options = $(init_nan) $(warn_un_init) $(implicit_none) $(compiled_modules_dir) $(enable_preprocessing) $(warn_array_temps)
else ifeq ($(fortran), ifort)
   fortran_options = -init=snan $(enable_preprocessing) $(compiled_modules_dir)
endif

# Add some options for debug mode:
ifeq ($(mode), debug)
   fortran_options += $(run_time_check)
endif

! For spectrum we need the FFTW library:
ifeq ($(MAKECMDGOALS),vel_spectrum)
   fortran_options += -I$(fftw_include_path)
   library_links += -L$(fftw-library-path) -lfftw3 -lm
endif

! For spectrum we need the FFTW library:
ifeq ($(MAKECMDGOALS),spectrum)
   fortran_options += -I$(fftw_include_path)
   library_links += -L$(fftw-library-path) -lfftw3 -lm
endif

! For spectrum we need the FFTW library:
ifeq ($(MAKECMDGOALS),enstrophy_spectrum)
   fortran_options += -I$(fftw_include_path)
   library_links += -L$(fftw-library-path) -lfftw3 -lm
   $(info Goal is enstrophy_spectrum.  fortran_options are $(fortran_options))
endif


exec1 = partition_tool
exec2 = time_average_profiles
exec3 = stretched
exec4 = spectrum
exec5 = lines
exec6 = plot_Melander
exec7 = post_process_conservation_diagnostics
exec8 = phi_time_averages
exec9 = post_process_stresses

# Computes phi-z-time averages from a list of files that contain phi-z averages.
# This is useful for runs in which the z direction is periodic as for
# Taylor-Couette flow.
exec10 = pzt_averages

exec11 = enstrophy_spectrum
exec12 = vel_spectrum
exec13 = make_log_log_line

objects1 = partition_tool.o
objects2 = time_average_profiles.o
objects3 = stretched.o
objects4 = spectrum.o
objects5 = lines.o
objects6 = plot_Melander.o
objects7 = post_process_conservation_diagnostics.o
objects8 = phi_time_averages.o
objects9 = post_process_stresses.o
objects10 = pzt_averages.o
objects11 = enstrophy_spectrum.o
objects12 = vel_spectrum.o
objects13 = make_log_log_line.o

$(exec1):	$(objects1)
		$(fortran) -o $(exec1) $(objects1)

$(exec2):	$(objects2)
		$(fortran) -o $(exec2) $(objects2)

$(exec3):	$(objects3)
		$(fortran) -o $(exec3) $(objects3)

$(exec4):	$(objects4)
		$(fortran) -o $(exec4) $(objects4) $(library_links)

$(exec5):	$(objects5)
		$(fortran) -o $(exec5) $(objects5)

$(exec6):	$(objects6)
		$(fortran) -o $(exec6) $(objects6)

$(exec7):	$(objects7)
		$(fortran) -o $(exec7) $(objects7)

$(exec8):	$(objects8)
		$(fortran) -o $(exec8) $(objects8)

$(exec9):	$(objects9)
		$(fortran) -o $(exec9) $(objects9)

$(exec10):	$(objects10)
		$(fortran) -o $(exec10) $(objects10)

$(exec11):	$(objects11)
		$(fortran) -o $(exec11) $(objects11) $(library_links)

$(exec12):	$(objects12)
		$(fortran) -o $(exec12) $(objects12) $(library_links)

$(exec13):	$(objects13)
		$(fortran) -o $(exec13) $(objects13)

# The dashes below make the command continue even after an error, e.g., the file does not exist.
clean:
	-\rm *.mod
	-\rm *.o
	-\rm $(exec1) $(exec2) $(exec3) $(exec4) $(exec5) $(exec6) $(exec7) $(exec8) $(exec9) $(exec10) $(exec11)

.f90.o :
	$(info fortran options are $(fortran_options))
	$(fortran) -O$(optimization_level) $(fortran_options) -c $<

.f.o :
	$(fortran) -O$(optimization_level) $(fortran_options) -c $<








