##############################################################################80

# To make the executable for the pade1.1 type:
# make pade1.1

# Command line options: print=yes/no checks=yes/no

# This section must be set-up by the user.

# Set this to yes if you want to perform the test of uniform advection in z and phi:
advection_test = no

# If parallel = yes, an mpi code will be made.  Else a serial code will be made.
parallel = no

# If hdf5 = yes, then the code will be linked with the hdf5 library permitting hdf5 file output.
hdf5 = no

# This is an option for certain nodes at NASA Ames, specifically, Electra Skylake
# and Aitken Cascade Lake.  It turns on the Advanced Vector Extensions optimizations.
AVX512 = no

AMD_Rome = no

# Turn on gdb debugger information for gfortran.  This will only work if you compile with
# checks = yes, i.e., you type "make pade checks = yes"
GDB_INFO = no

# Traceback (with -g which provides symbolic information) for the ifort compiler.  This makes is costly.
ifort_traceback = no

# In the next two lines, specify how to invoke the serial and parallel fortran
# compilers on your machine.  If you only want to make the serial code
# then don't worry about the mpi_fortran variable.

# Note: Currently, this makefile sets compiler options (via the variable $(fortran_options))
# for only the gfortran and ifort compilers.  

ifeq ($(OSTYPE),darwin)
#  For my mac:
   fortran_compiler = gfortran
   mpi_fortran = mpif90
else ifeq ($(HOST),linux261.nas.nasa.gov)
#  For my linux box:
   fortran_compiler = gfortran
   mpi_fortran = mpif90
else ifeq ($(MYHOST),pfe)
#  For the Pleiades supercomputer at NASA Ames:
# Note: In my .tcshrc I have setenv MYHOST pfe
   fortran_compiler = ifort
   mpi_fortran = mpif90
endif

# In the next line, select the fft routine to be used by the corrected FARGO method.
# The two choices are fftw or rogallo.  For fftw you need to specify the include paths
# below for the fftw library.  The rogallo fft comes with the pade code(see rogallo.f90).
#fft = rogallo
fft = fftw

# Set this for your situation (needed only if you chose fftw above).
ifeq ($(fft),fftw)
   ifeq ($(OSTYPE),darwin)
      fftw_include_path = $(HOME)/homebrew/opt/fftw/include
      fftw-library-path = $(HOME)/homebrew/opt/fftw/lib
   else ifeq ($(HOST),linux261.nas.nasa.gov)
      fftw_include_path = $(HOME)/WORK/fftw-3.3.8/include
      fftw-library-path = $(HOME)/WORK/fftw-3.3.8/lib
   else
#     Pleiades:
      fftw_include_path = $(HOME)/Applications/fftw-3.3.7/include
      fftw-library-path = $(HOME)/Applications/fftw-3.3.7/lib
   endif
endif

# Set this for your situation (needed only if you chose hdf5 = yes above).
ifeq ($(hdf5),yes)
   ifeq ($(OSTYPE),darwin)
      hdf5_include_path = $(HOME)/homebrew/include
      hdf5_library_path = $(HOME)/homebrew/lib
   else
#     Pleiades:
      hdf5_include_path = $(HOME)/Applications/fftw-3.3.7/include
      hdf5_include_path = $(HOME)/Applications/fftw-3.3.7/lib
   endif
   library_links += -L$(hdf5_library_path) -lhdf5
endif

# End of section to be set-up by the user.

############# make commands you may type #################################

# You can type the following commands to the terminal while in SRC 
# directory or its copy.
#
# (1) make pade ! To make the executable called pade
#
# (2) make clean ! Deletes the executable (pade), all object
# (.o) files and all compiled modules (.mod) files.
#
# (3) Options for debugging: print=yes; checks=yes
#     e.g. make pade print=yes checks=yes
#     (a) print=yes gives detailed output of what the code is doing.
#     (b) checks=yes enables various run-time checks, traceback, traps, etc.
#     Both of these will consume time and so should not be used for production runs.

# (4) Timing option.
#     make pade transpose_timing=yes
#     This will perform cpu timing for transposes.  The result will be output to
#     stdout at the end of every step and at the end of the run.

# You don't have to read any further if you don't want to.

#########################################################################

# I have currently commented out the call to a Bessel function for Seligman's
# vortex so that the user need not worry about having the math77 library.

# If you uncomment the call to the bessel function then you need to set math77 = yes
# and provide the path for the math77 library.
# If you set math77 to yes then you need to supply the location of the library.
# Note: Pleiades has the math77 library on the system so you don't need to do anything.
math77 = no
ifeq ($(math77),yes)
   ifeq ($(OSTYPE),darwin)
      math77-library-path = /Users/kshariff/MY_APPLICATIONS/math77-master/src/MATH77
   else ifeq ($(HOST),linux261.nas.nasa.gov)
      math77-library-path = $(HOME)/WORK/MATH77
   endif
   library_links += -L$(math77-library-path) -lmath77
endif

# Set this to yes if you will invoke subroutines to gather phi averaged statistics and/or
# history diagnostics (velocity rms in the domain and/or averaged fluctuation kinetic energy
# in the domain)
phi_averaged_statistics_and_history_diagnostics = yes

.SUFFIXES : .f .f90 .o

ifeq ($(parallel), yes)
   fortran = $(mpi_fortran)
else
   fortran = $(fortran_compiler)
endif

# Preprocessor options:
enable_preprocessing = -cpp
ifeq ($(advection_test), yes)
   fortran_options += -Dadvection_test 
endif

ifeq ($(print), yes)
   fortran_options += -Ddebug_print 
endif

ifeq ($(parallel), yes)
   fortran_options += -Dmpi_code
endif

ifeq ($(fft), rogallo)
   fortran_options += -Drogallo_fft
endif

# Run time check for debug mode.
ifeq ($(fortran_compiler),gfortran)
#   run_time_checks = -fbounds-check -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace
   run_time_checks = -fbounds-check -fcheck=all -fbacktrace
   ifeq ($(GDB_INFO),yes)
      run_time_checks += -g
   endif
else ifeq ($(fortran_compiler),ifort)
   run_time_checks = -check all -fpe0 -g -traceback
endif

# List of possible compiler options I like to use for gfortran:
init_nan = -finit-real=nan
warn_un_init = -Wuninitialized
# Warn about array temporaries created by the compiler
warn_array_temps = -Warray-temporaries
# This has the effect of putting implicit none at the start of every subroutine
# in case I forget.
implicit_none = -fimplicit-none
# In case needed for linking with c objects.
no_underscoreing = -fno-underscoring -fno-second-underscoring
max_errors = -fmax-errors=10
no_max_line_len = -ffree-line-length-none

#  This is currently not used.
# Directory in which .mod files will be kept.  I simply do not like to see all those
# .mod files mixed in with the .f90 files.
#ifeq ($(fortran_compiler),gfortran)
#   compiled_modules_dir = -JCOMPILED_MODULES
#else ifeq ($(fortran),ifort)
#   compiled_modules_dir = -module COMPILED_MODULES
#endif

optimization_level = 3

# For sanity checks, the code checks for uninitialized variables by comparing against NaN, so
# DO NOT remove the init_nan or init=snan option below:

# Fortran options I always use:
ifeq ($(fortran_compiler),gfortran)
   fortran_options += $(init_nan) $(warn_un_init) $(implicit_none) $(compiled_modules_dir) $(enable_preprocessing) $(warn_array_temps) $(max_errors) \
                     $(no_max_line_len)
else ifeq ($(fortran_compiler),ifort)
   fortran_options += -init=snan $(enable_preprocessing) $(compiled_modules_dir)
endif

ifeq ($(checks), yes)
   fortran_options += $(run_time_checks) 
endif

ifeq ($(fortran_compiler), ifort)
ifeq ($(ifort_traceback), yes)
   fortran_options += -g -traceback 
endif
endif


# Add pre-processor variable if needed:
ifeq ($(transpose_timing), yes)
   fortran_options += -Dtranspose_timing
endif

# Add options for certain nodes at NASA Ames:
ifeq ($(MYHOST),pfe)
ifeq ($(AVX512),yes)
   fortran_options += -xCORE-AVX512 -qopt-zmm-usage=high
endif
endif

ifeq ($(MYHOST),pfe)
ifeq ($(AMD_Rome),yes)
   fortran_options += -march=core=avx2
endif
endif

# This is for gfortran:
####fortran_options += -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace

# For FFTW to implement fargo shifts:
ifeq ($(fft), fftw)
   fortran_options += -I$(fftw_include_path)
   library_links += -L$(fftw-library-path) -lfftw3 -lm
#  This adds an ifdef preprocessor option called fftw for the code to use:
   fortran_options += -Dfftw
   ifeq ($(parallel), yes)
      library_links += -L$(fftw-library-path) -lfftw3_mpi
   endif
endif

ifeq ($(hdf5), yes)
   fortran_options += -I$(hdf5_include_path)
   library_links += -L$(hdf5_library_path) -lhdf5_fortran
#  This adds an ifdef preprocessor option cfor the code to use:
   fortran_options += -Dhdf5
   #ifeq ($(parallel), yes)
   #   library_links += -L$(hdf5_library_path) -lhdf5_mpi
   #endif
endif

#fortran_options = -fbacktrace
# If we are using the ifort compiler then the code needs to issue a
# use IFPORT
# the ifort portability library in order than the gfortran random number
# generator can be used.  This pre-processor variable allows us to do that.
ifeq ($(fortran_compiler), ifort)
   fortran_options += -Difort
endif

exec1 = pade1.1

# These were used for micro testing purposes while the code was being developed.
exec3 = pade_diff_test
exec4 = Shu_Osher_test
exec5 = phi_advection_test
exec6 = hyperslab_by_col

# For micro testing purposes as the code was being developed.
objects3 = pade_diff_test.o pade_diff.o tridiagonal_solvers.o
objects4 = Shu_Osher_test.o
objects5 = phi_advection_test.o flux_diff.o tridiagonal_solvers.o
objects6 = hyperslab_by_col.o

# The order of the .o files is important to ensure that the modules required by succeeding
# objects have been compiled.

# For some reason I have to declare these dependences for files that contain a module containing subroutines.
# Otherwise, if I modify these files, they do not get re-compiled.
activate_routines.o:	activate_routines.f90
			$(fortran) -O$(optimization_level) $(fortran_options) -c activate_routines.f90

set_up_routines.o:	set_up_routines.f90
			$(fortran) -O$(optimization_level) $(fortran_options) -c set_up_routines.f90

grid.o:			grid.f90
			$(fortran) -O$(optimization_level) $(fortran_options) -c grid.f90


# For some reason I had to change the name of sponge.f90 to sponge1.f90
# otherwise "make" does not recognize changes to the file.
# Objects for pade code.  The main program is in pade.f90
objects_for_pade =	modules.o \
			pade.o \
			basic_state.o \
			smooth.o \
			rogallo_fft.o \
			fargo_and_plotting_shift.o \
			gravity.o \
			boundary_conditions.o \
			viscous_terms.o \
			artificial_pressure.o \
			set_up_routines.o \
			sponge1.o \
			pade_diff.o \
			grid.o \
			domain_mesh_and_partition.o \
			add_derivative_routines.o \
			rotating_frame.o \
			rhs.o \
			control_and_check_routines.o \
			time_stepping_routines.o \
			tridiagonal_solvers.o \
			transpose_routines.o \
			plotting_output.o \
			restart_and_save_file_operations.o \
			vorticity_and_dilatation.o \
			baroclinic_term.o \
			Cassen_Moosman_infall.o \
			zeroin.o \
			conservation_diagnostics.o \
			vorticity_equation_terms.o \
			activate_routines.o

ifeq ($(phi_averaged_statistics_and_history_diagnostics), yes)
   objects_for_pade += phi_avg_statistics_and_history_diagnostics.o
endif

ifeq ($(advection_test), yes)
   objects_for_pade += app_advection_test.o
endif

# The appropriate subroutine is invoked in the main program in
# source file pade.f90 according to the i_run_type variable
# specified in the input file.
application_objects =	app_user.o \
			app_euler1d_tests.o \
			app_hydrostatic_test.o \
			app_homentropic_solid_body_rotation_test.o \
			app_vsi_3D.o \
			app_single_vortex_fargo_test.o \
			app_vortex_pair.o \
			app_taylor_couette.o \
			app_vsi_with_dt_output.o

objects_for_pade += $(application_objects)

$(exec1):	$(objects_for_pade)
		$(fortran) $(compiled_modules_dir) -o $(exec1) $(objects_for_pade) $(library_links)

$(exec3):	$(objects3)
		$(fortran) $(compiled_modules_dir) -o $(exec3) $(objects3)

$(exec4):	$(objects4)
		$(fortran) $(compiled_modules_dir) -o $(exec4) $(objects4)

$(exec5):	$(objects5)
		$(fortran) $(compiled_modules_dir) -o $(exec5) $(objects5)

$(exec6):	$(objects6)
		$(fortran) $(compiled_modules_dir) -o $(exec6) $(objects6) $(library_links)

$(info parallel is $(parallel))
$(info fortran is $(fortran))
$(info mpi_fortran is $(mpi_fortran))

# FORD is an automatic html documentation generator for modern (1990 onward) Fortran code.
# To create updated documentation, simply type
# make ford_doc
# You can navigate through the documentation by going to your browser and starting with
# /PADE_CODE/FORD_DOCUMENTATION/index.html

ford_doc:
		ford -I$(fftw_include_path) ../project-file.md

# The dashes below make the command continue even after an error, e.g., if the file does not exist.
clean:
	-\rm COMPILED_MODULES/*.mod
	-\rm *.mod
	-\rm *.o
	-\rm $(exec1) $(exec2) $(exec3) $(exec4) $(exec5) $(exec6)

.f90.o :
	$(fortran) -O$(optimization_level) $(fortran_options) -c $<

.f.o :
	$(fortran) -O$(optimization_level) $(fortran_options) -c $<








