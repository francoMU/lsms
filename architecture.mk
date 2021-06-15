export USE_OPENMP=1
export USE_LIBXC=1
export HAS_BACKTRACE=1


export LIBS +=
export ADD_LIBS += -L/usr/lib/x86_64-linux-gnu/ \
                    -lhdf5_serial \
				-llapack \
					-lblas \
                    -lgfortran

export INC_PATH += -I/usr/include/hdf5/serial/



ifdef USE_OPENMP
  export CXX=mpicxx.openmpi -g -std=c++14 -fopenmp -O2 $(OPT_DEFINES) -mtune=skylake -march=skylake -Wall
  export F77=gfortran -g -fopenmp -O2 -fbacktrace -mtune=skylake -march=skylake -finit-local-zero
else
  export CXX=mpicxx.openmpi -g -O2 -std=c++14  $(OPT_DEFINES) -mtune=skylake -march=skylake -Wall
  export F77=gfortran -g -O2 -fbacktrace -mtune=skylake -march=skylake -finit-local-zero
endif

# JSON++ requires flex and bison (version>=2.4)
export FLEX=flex
export BISON=bison

export LUACXX = $(CXX)
