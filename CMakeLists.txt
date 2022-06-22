#
#  This is a sample CMakeLists.txt suitable for direct use with a new PETSc application or
#  you can add fragments of material below to an existing applications CMakeLists.txt
#
#  Usage:  To build ex1 from ex1.c
#    rm -fr build
#    mkdir build
#    cd build
#    cmake ..
#    make
#    ./ex1
#
#  By default it gets both the compiler information and the library information from PETSc
#
cmake_minimum_required(VERSION 3.1.0)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# set root of location to find PETSc's pkg-config
set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig:$ENV{SLEPC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig)
#set(MKL_DIR $ENV{MKLROOT}/lib/cmake/mkl)
#set(SLEPC $ENV{SLEPC_DIR}/$ENV{PETSC_ARCH})
#set(ENV{PKG_CONFIG_PATH} ${SLEPC}/lib/pkgconfig)

# Remove the lines below if you do not wish to have PETSc determine the compilers
execute_process ( COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
#execute_process ( COMMAND pkg-config SLEPc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
SET(CMAKE_C_COMPILER ${C_COMPILER})
execute_process ( COMMAND pkg-config PETSc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
#execute_process ( COMMAND pkg-config SLEPc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
if (CXX_COMPILER)
  SET(CMAKE_CXX_COMPILER ${CXX_COMPILER})
endif (CXX_COMPILER)
execute_process ( COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
#execute_process ( COMMAND pkg-config SLEPc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
if (FORTRAN_COMPILER)
  SET(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
  enable_language(Fortran)
endif (FORTRAN_COMPILER)
add_compile_options(-cpp)
#  tells CMake to build the application ex1 from the source file ex1.c
#  this must appear AFTER the compilers are set
project(START_LNS)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
SET(CMAKE_INSTALL_PREFIX $ENV{HOME})
find_package(PkgConfig REQUIRED)
file(GLOB  SRC_FILES
    "${PROJECT_SOURCE_DIR}/src/*.f90"
    "${PROJECT_SOURCE_DIR}/src/*.F90"
    "${PROJECT_SOURCE_DIR}/src/third_party/cfgio/*.f90"
    "${PROJECT_SOURCE_DIR}/src/third_party/PENF/*.F90"
    "${PROJECT_SOURCE_DIR}/src/third_party/BeFoR64/*.F90"
    "${PROJECT_SOURCE_DIR}/src/third_party/StringiFor/*.F90"
    )
add_executable(START_LNS ${SRC_FILES})

#set(SLEPC $ENV{SLEPC_DIR}/$ENV{PETSC_ARCH})
#set(ENV{PKG_CONFIG_PATH} ${SLEPC}/lib/pkgconfig)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSC)
pkg_search_module(SLEPC REQUIRED IMPORTED_TARGET SLEPC)
#find_package(MKL)
target_link_libraries(START_LNS PkgConfig::PETSC PkgConfig::SLEPC lapack)# MKL::MKL)
install(TARGETS START_LNS DESTINATION bin)
