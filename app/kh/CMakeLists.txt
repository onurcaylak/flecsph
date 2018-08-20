#----------------------------------------------------------------------------#
# Copyright (c) 2017 Los Alamos National Security, LLC
# All rights reserved.
#----------------------------------------------------------------------------#

 #------------------------------------------------------------------------------#
# Debug and release flags
#------------------------------------------------------------------------------#

set(CMAKE_CXX_FLAGS_DEBUG
  "-g -O2 -Wall -Wno-return-type -Wno-unused -Wno-comment -Wno-parentheses")
set(CMAKE_CXX_FLAGS_RELEASE
  "-O3 -Wall -Wno-return-type -Wno-unused -Wno-comment -Wno-parentheses")

# includes

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/physics)
include_directories(${CMAKE_SOURCE_DIR}/app/kh/include)
include_directories(${CMAKE_SOURCE_DIR}/mpisph)
include_directories(${CMAKE_SOURCE_DIR}/third-party-libraries/ScalingFramework/IOTests)

#------------------------------------------------------------------------------#
# Executables
#------------------------------------------------------------------------------#

# MPILEGION version
add_executable(kh
  main.cc
  main_driver.cc
  ${FleCSI_RUNTIME}/runtime_driver.cc
)

# MPILEGION version
add_executable(kh_generator
  generator/main.cc
)

target_link_libraries(kh_generator ${MPI_LIBRARIES}
  ${H5hut_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(kh ${FleCSPH_LIBRARIES})

add_definitions(-DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT)
add_definitions(-DFLECSI_OVERRIDE_DEFAULT_SPECIALIZATION_DRIVER)

#------------------------------------------------------------------------------#
# kh test, call the default parameter file
#------------------------------------------------------------------------------#

cinch_add_unit(kh_test
  SOURCES
    test/kh.cc
    ${CMAKE_SOURCE_DIR}/app/kh/main_driver.cc
    ${FleCSI_RUNTIME}/runtime_driver.cc
  INPUTS
    ${PROJECT_SOURCE_DIR}/data/kh_sqn100.par
  LIBRARIES
    ${FleCSPH_LIBRARIES}
)

if(ENABLE_UNIT_TESTS)
# run the generator with default file
add_custom_target(
    kh_generator_test
    COMMAND kh_generator ${PROJECT_SOURCE_DIR}/data/kh_sqn100.par
    COMMAND mv kh_sqn100.h5part ../../test/kh
)
add_dependencies(kh_test kh_generator_test)
endif()
