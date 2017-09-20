#PLATFORM is either mac, machome or linux

SOURCE_DIR=./

# Binary type
ifeq (${BUILD64},0)
	BINTYPE=32bit
else
	BINTYPE=64bit
endif

ifeq (${PLATFORM},Darwin)
    MACPORTS_DIR=/opt/local
    GCC_PRE = g++
    GCC_COMP = clang++
    PYTHON = python
    ARPACK_LIB_DIR=${MACPORTS_DIR}/lib
    ARPACK_LIB=-larpack
    LAPACK_BLAS_LIBS=-llapack -lblas -latlas
    LAPACK_BLAS_LIB_DIR=/usr/lib
    GFORT_SYSTEM_LIB_DIR=${MACPORTS_DIR}/lib/gcc48
    GFORT_SYSTEM_LIBS=-lgfortran -lgomp
    PLPLOT_INCLUDE_DIR=${MACPORTS_DIR}/include
    HAVE_PLPLOT_LIBS=0
    PLPLOT_LIB_DIR=${MACPORTS_DIR}/lib
    PLPLOT_LIBS=-lplplotcxxd -lplplotd
    METIS_INCLUDE_DIR=${MACPORTS_DIR}/include/metis
    METIS_LIB_DIR=${MACPORTS_DIR}lib
    METIS_LIBS=-lmetis
    TBB_LIB_DIR=${MACPORTS_DIR}/tbb41_20121003oss/lib
    TBB_INCLUDE_DIR=${MACPORTS_DIR}/tbb41_20121003oss/include
    TBB_LIBS=-ltbb -ltbbmalloc
    BOOST_DIR=${MACPORTS_DIR}/include
    LIBMESH_DIR=../../../libmesh
    TRILINOS_DIR=/usr/local/trilinos
    VTK_INCLUDE_DIR=${MACPORTS_DIR}/include/vtk-5.10
    MPI_INCLUDE_DIR=${MACPORTS_DIR}/include/openmpi
    MPI_LIB_DIR=${MACPORTS_DIR}/lib/openmpi
    GCC_FLAGS_DEBUG = -g3 -O0 -Wall 
    GCC_FLAGS_OPT = -DNDEBUG -O2 -felide-constructors -Qunused-arguments 
    GCC_FLAGS_PRO = -DNDEBUG -O2 -g -felide-constructors -Qunused-arguments 
endif

ifeq (${PLATFORM},Linux)
    GCC_PRE = g++
    GCC_COMP = g++
    PYTHON = python
    ARPACK_LIB_DIR=./
    ARPACK_LIB=./
    LAPACK_BLAS_LIBS=./
    LAPACK_BLAS_LIB_DIR=./
    GFORT_SYSTEM_LIB_DIR=./
    GFORT_SYSTEM_LIBS=./
    PLPLOT_INCLUDE_DIR=./
    HAVE_PLPLOT_LIBS=0
    PLPLOT_LIB_DIR=./
    PLPLOT_LIBS=./
    METIS_INCLUDE_DIR=/home/bhatiam/usr/include
    METIS_LIB_DIR=/home/bhatiam/usr/lib
    METIS_LIBS=./
    HDF_LIB_DIR=/apps/global/oss/hdf5/1.8.0/lib
    TBB_LIB_DIR=/usr/local/intel/2013-u1/tbb/lib/intel64
    TBB_INCLUDE_DIR=/usr/local/intel/2013-u1/tbb/include
    TBB_LIBS=-ltbb -ltbbmalloc
    BOOST_DIR=./
    LIBMESH_DIR=../../../libmesh
    TRILINOS_DIR=./
    VTK_LIB_DIR=/home/bhatiam/usr/vtk/lib/vtk-5.10
    VTK_LIBS=-lvtkIO -lvtkCommon -lvtkFiltering
    VTK_INCLUDE_DIR=/home/bhatiam/usr/vtk/include/vtk-5.10
    MPI_INCLUDE_DIR=/usr/local/openmpi/openmpi-1.6.3/x86_64/ib/gcc/include
    MPI_LIB_DIR=/usr/local/openmpi/openmpi-1.6.3/x86_64/ib/gcc/lib
    GCC_FLAGS_DEBUG = -g3 -O0 -Wall -std=c++0x
    GCC_FLAGS_OPT = -DNDEBUG -O2 -felide-constructors -std=c++0x
    GCC_FLAGS_PRO = -DNDEBUG -O2 -g -felide-constructors -std=c++0x
endif

ifeq (${PLATFORM},script)
    ARPACK_LIB_DIR=${BASEDIR}/install/${KERNEL}-${BINTYPE}/lib
    LAPACK_BLAS_LIB_DIR=${BASEDIR}/install/${KERNEL}-${BINTYPE}/lib
    GFORT_SYSTEM_LIB_DIR=${GCCHOME}/lib
    GFORT_SYSTEM_LIBS=-lgfortran -lgomp
    HAVE_PLPLOT_LIBS=0
    PLATFORM=linux
endif


#GCC_FLAGS_OPT = -g0 -Os 

D      = cd
ECHO    = echo
LN      = ln
LNFLAGS = -s
MAKE    = make
RM      = rm
RMFLAGS = -f
SHELL   = /bin/sh
ARFLAGS = rv
RANLIB   = ranlib
AR = ar


LIBMESH_DIR=../contrib/libmesh
DEALII_DIR=../contrib/deal.II
SRC_SUBDIR=../src
APP_SUBDIR=../applications
BUILD_DIR=${PWD}

EXTERNAL_INCLUDE_DIRS+=-I${}
EXTERNAL_INCLUDE_DIRS+=-I${SRC_SUBDIR}
EXTERNAL_INCLUDE_DIRS+=-I${SRC_SUBDIR}/../examples

EXTERNAL_LIB_DIRS=-L${SOURCE_DIR}

EXTERNAL_LIBS+=${ARPACK_LIB}
ifeq (${HAVE_PLPLOT_LIBS},1)
    EXTERNAL_LIBS+=${PLPLOT_LIBS}
endif
LIBMESH_LIBS=-lmesh_${METHOD}
MPI_LIBS=-lmpi_cxx -lmpi


ifeq (${HAVE_PLPLOT_LIBS},1)
    COMPILATION_DEFINES=-DHAVE_PLPLOT
else
    COMPILATION_DEFINES=
endif


FESYSTEM_LIB_SOURCES = $(wildcard ${SRC_SUBDIR}/*.cpp ${SRC_SUBDIR}/*.C ${SRC_SUBDIR}/*/*.cpp ${SRC_SUBDIR}/*/*.C ${SRC_SUBDIR}/*/*/*.cpp ${SRC_SUBDIR}/*/*/*.C)
DEP_FILES= $(patsubst %.cpp,%-cpp.dep,${FESYSTEM_LIB_SOURCES} ${FESYSTEM_APP_SOURCES} ${FESYSTEM_EXAMPLE_SOURCES})

all: fesystem-${METHOD} euler-${METHOD}-${SCALAR}

fesystem-dbg: $(patsubst %.cpp,%-${METHOD}-cpp.o,${FESYSTEM_EXAMPLE_SOURCES})  libfesystem-dbg.a 
	${GCC_COMP} ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $^ ${EXTERNAL_LIB_DIRS} ${EXTERNAL_LIBS} ${MPI_LIBS} -o $@

fesystem-opt: $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_EXAMPLE_SOURCES})  libfesystem-opt.a 
	${GCC_COMP} ${GCC_FLAGS_OPT} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $^ ${EXTERNAL_LIB_DIRS} ${EXTERNAL_LIBS} ${MPI_LIBS} -o $@

%-cpp.dep: %.cpp 
	${GCC_PRE} -M  ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp > $@
	${PYTHON} ${BUILD_DIR}/modify_dep_file.py $@ ${*D}

%-dbg-cpp.o: %.cpp %.C
	${GCC_COMP} -c ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-opt-cpp.o: %.cpp %.C
	${GCC_COMP} -c ${GCC_FLAGS_OPT} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-dbg-real-cpp.o: %.cpp %.C
	${GCC_COMP} -c ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-opt-real-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_OPT} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-oprof-real-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_PRO} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-dbg-complex-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-opt-complex-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_OPT} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-oprof-complex-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_PRO} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@


libfesystem-${METHOD}.a: $(patsubst %.cpp,%-${METHOD}-cpp.o,${FESYSTEM_LIB_SOURCES}) 
	ar -rv $@ $^ 

clean:
	rm -rf $(patsubst %.cpp,%-${METHOD}-cpp.o,${FESYSTEM_EXAMPLE_SOURCES}) $(patsubst %.cpp,%-${METHOD}-cpp.o,${FESYSTEM_LIB_SOURCES})  libfesystem-${METHOD}.a  fesystem-${METHOD}
	rm -rf $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_APP_SOURCES})  euler-${METHOD}-${SCALAR}
	rm -f ${DEP_FILES}

-include ${DEP_FILES}

