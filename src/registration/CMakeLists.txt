## Version: $Id$
##
######################################################################
##
### Commentary:
##
######################################################################
##
### Change Log:
##
######################################################################
##
### Code:

project(mrtrix3Registration)

## #################################################################
## Sources
## #################################################################

set(${PROJECT_NAME}_HEADERS
  linear.h
  multi_resolution_lmax.h
  nonlinear.h
  metric/cross_correlation.h
  metric/demons4D.h
  metric/demons.h
  metric/difference_robust.h
  metric/evaluate.h
  metric/mean_squared.h
  metric/normalised_cross_correlation.h
  metric/params.h
  metric/robust_estimators.h
  metric/thread_kernel.h
  transform/affine.h
  transform/base.h
  transform/convergence_check.h
  transform/initialiser.h
  transform/initialiser_helpers.h
  transform/reorient.h
  transform/rigid.h
  transform/search.h
  warp/compose.h
  warp/convert.h
  warp/helpers.h
  warp/invert.h
)

set(${PROJECT_NAME}_SOURCES
  linear.cpp
  nonlinear.cpp
  transform/affine.cpp
  transform/convergence_check.cpp
  transform/initialiser.cpp
  transform/initialiser_helpers.cpp
  transform/rigid.cpp
)

## #################################################################
## Build rules
## #################################################################

add_library(${PROJECT_NAME} SHARED
  ${${PROJECT_NAME}_SOURCES}
  ${${PROJECT_NAME}_HEADERS})

## ###################################################################
## Link rules
## ###################################################################

## #################################################################
## Target properties
## #################################################################

set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_RPATH 0)
set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION   ${mrtrix3_VERSION}
                                                 SOVERSION ${mrtrix3_VERSION_MAJOR})

## #################################################################
## Export header file
## #################################################################

generate_export_header(${PROJECT_NAME} EXPORT_FILE_NAME "${CMAKE_BINARY_DIR}/${PROJECT_NAME}Export.h")
generate_export_header(${PROJECT_NAME} EXPORT_FILE_NAME "${CMAKE_BINARY_DIR}/${PROJECT_NAME}Export")

set(${PROJECT_NAME}_HEADERS
  ${${PROJECT_NAME}_HEADERS}
 "${CMAKE_BINARY_DIR}/${PROJECT_NAME}Export"
 "${CMAKE_BINARY_DIR}/${PROJECT_NAME}Export.h")

## ###################################################################
## Install rules - files
## ###################################################################

install(FILES ${${PROJECT_NAME}_HEADERS}
  DESTINATION include/${PROJECT_NAME}
  COMPONENT registration)

## ###################################################################
## Install rules - targets
## ###################################################################

install(TARGETS ${PROJECT_NAME}
  EXPORT mrtrix3Depends
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          COMPONENT   registration)

######################################################################
### CMakeLists.txt ends here
