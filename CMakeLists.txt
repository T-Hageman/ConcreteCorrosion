cmake_minimum_required(VERSION 3.12...3.27)
project(CorrosionFEM   VERSION 1.0
                    DESCRIPTION "Remind me to type this"
                    LANGUAGES CXX)

set(CMAKE_BUILD_TYPE "Release")   #Debug or Release

if(DEFINED ENV{CurrentEnv})
  set(CurrentEnv2 "hpc")
else ()
  set(CurrentEnv2 "local")
endif()

option(MI_XMALLOC           "Enable abort() call on memory allocation failure by default" ON)

## Get all files in folder and subfolders, excluding build folder
set (EXCLUDE_DIR "/CMakeFiles/")
set (EXCLUDE_DIR2 "/build/")
file (GLOB_RECURSE source_list "*.cpp" "*.h")
foreach (TMP_PATH ${source_list})
    string (FIND ${TMP_PATH} ${EXCLUDE_DIR} EXCLUDE_DIR_FOUND)
    string (FIND ${TMP_PATH} ${EXCLUDE_DIR2} EXCLUDE_DIR2_FOUND)
    if (NOT ${EXCLUDE_DIR_FOUND} EQUAL -1)
        list (REMOVE_ITEM source_list ${TMP_PATH})
    endif ()
    if (NOT ${EXCLUDE_DIR2_FOUND} EQUAL -1)
        list (REMOVE_ITEM source_list ${TMP_PATH})
    endif ()
endforeach(TMP_PATH)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(CurrentEnv2 MATCHES "local")
  add_definitions(-D LOCALENVIRONMENT=1)

  find_package(PkgConfig REQUIRED)

  # set root of location to find PETSc's pkg-config
  set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
  set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)
  pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)

  # Remove the lines below if you do not wish to have PETSc determine the compilers
  execute_process ( COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
  SET(CMAKE_C_COMPILER ${C_COMPILER})
  execute_process ( COMMAND pkg-config PETSc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
  if (CXX_COMPILER)
    SET(CMAKE_CXX_COMPILER ${CXX_COMPILER})
  endif (CXX_COMPILER)
  execute_process ( COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
  if (FORTRAN_COMPILER)
    SET(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
    enable_language(Fortran)
  endif (FORTRAN_COMPILER)

  find_package(VTK COMPONENTS 
    ChartsCore
    CommonColor
    CommonCore
    CommonDataModel
    FiltersSources
    InteractionStyle
    RenderingAnnotation
    RenderingContext2D
    RenderingContextOpenGL2
    RenderingCore
    RenderingFreeType
    RenderingGL2PSOpenGL2
    RenderingOpenGL2
    TestingRendering
    ViewsContext2D
    CommonTransforms
    FiltersCore
    FiltersGeneral
    FiltersPoints
    ImagingHybrid
    InteractionStyle
  )

  include_directories("$ENV{LIBRARY_DIR}/eigen")
  include_directories("$ENV{LIBRARY_DIR}/rapidjson/include")
  include_directories("$ENV{LIBRARY_DIR}/HighFive/include")
  link_directories("$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/externalpackages/hdf5-1.14.3-p1/src/.libs")

else ()
  add_compile_options(-diag-disable=10441)
  find_package(PkgConfig REQUIRED)

  #eigen
  include_directories("$ENV{LIBRARY_DIR}/eigen")
  include_directories("$ENV{LIBRARY_DIR}/rapidjson/include")
  include_directories("$ENV{LIBRARY_DIR}/HighFive/include")
  

  #hdf5
  link_directories("$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib")

  #petsc
  set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
  set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)
  pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)

  # Remove the lines below if you do not wish to have PETSc determine the compilers
  execute_process ( COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
  SET(CMAKE_C_COMPILER ${C_COMPILER})
  execute_process ( COMMAND pkg-config PETSc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
  if (CXX_COMPILER)
    SET(CMAKE_CXX_COMPILER ${CXX_COMPILER})
  endif (CXX_COMPILER)
  execute_process ( COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
  if (FORTRAN_COMPILER)
    SET(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
    enable_language(Fortran)
  endif (FORTRAN_COMPILER)

endif()

add_executable(CorrosionFEM ${source_list})

target_link_libraries(CorrosionFEM PkgConfig::PETSC)
target_link_libraries(CorrosionFEM -lhdf5)
if (CurrentEnv2 MATCHES "local")
  target_link_libraries(CorrosionFEM ${VTK_LIBRARIES})
  vtk_module_autoinit(
    TARGETS CorrosionFEM
    MODULES ${VTK_LIBRARIES}
  )
else()
  add_compile_options(-diag-disable=10441)
endif()

