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

project(mrtrix3Core)

## #################################################################
## Sources
## #################################################################

set(${PROJECT_NAME}_HEADERS
  ../exec_version.h
  app.h
  apply.h
  axes.h
  bitset.h
  command.h
  datatype.h
  debug.h
  exception.h
  header.h
  image.h
  memory.h
  mrtrix.h
  phase_encoding.h
  progressbar.h
  raw.h
  signals.h
  signal_handler.h
  stats.h
  stride.h
  thread.h
  timer.h
  transform.h
  types.h
  adapter/base.h
  adapter/extract.h
  adapter/gaussian1D.h
  adapter/gradient1D.h
  adapter/gradient3D.h
  adapter/jacobian.h
  adapter/median.h
  adapter/neighbourhood3D.h
  adapter/normalise3D.h
  adapter/permute_axes.h
  adapter/replicate.h
  adapter/reslice.h
  adapter/subset.h
  adapter/warp.h
  algo/copy.h
  algo/histogram.h
  algo/iterator.h
  algo/loop.h
  algo/min_max.h
  algo/neighbourhooditerator.h
  algo/random_loop.h
  algo/random_threaded_loop.h
  algo/stochastic_threaded_loop.h
  algo/threaded_copy.h
  algo/threaded_loop.h
  connectome/connectome.h
  connectome/enhance.h
  connectome/lut.h
  connectome/mat2vec.h
  file/config.h
  file/copy.h
  file/entry.h
  file/gz.h
  file/json.h
  file/json_utils.h
  file/key_value.h
  file/mgh.h
  file/mmap.h
  file/name_parser.h
  file/nifti1.h
  file/nifti1_utils.h
  file/nifti2.h
  file/nifti2_utils.h
  file/nifti_utils.h
  file/ofstream.h
  file/path.h
  file/tiff.h
  file/utils.h
  file/dicom/csa_entry.h
  file/dicom/definitions.h
  file/dicom/element.h
  file/dicom/image.h
  file/dicom/mapper.h
  file/dicom/patient.h
  file/dicom/quick_scan.h
  file/dicom/series.h
  file/dicom/study.h
  file/dicom/tree.h
  fixel/helpers.h
  fixel/keys.h
  fixel/loop.h
  fixel/legacy/fixel_image.h
  fixel/legacy/fixel_keys.h
  fixel/legacy/fixel_metric.h
  formats/list.h
  formats/mrtrix_utils.h
  image_io/base.h
  image_io/default.h
  image_io/fetch_store.h
  image_io/gz.h
  image_io/mosaic.h
  image_io/pipe.h
  image_io/ram.h
  image_io/scratch.h
  image_io/tiff.h
  image_io/sparse.h
  math/average_space.h
  math/bessel.cpp
  math/bessel.h
  math/cauchy.h
  math/chebyshev.h
  math/check_gradient.h
  math/constrained_least_squares.h
  math/cubic_spline.h
  math/gaussian.h
  math/golden_section_search.h
  math/gradient_descent_bb.h
  math/gradient_descent.h
  math/hermite.h
  math/least_squares.h
  math/legendre.h
  math/math.h
  math/median.h
  math/quadratic_line_search.h
  math/rician.h
  math/rng.h
  math/sech.h
  math/SH.h
  math/sinc.h
  math/Sn_scale_estimator.h
  math/sphere.h
  math/versor.h
  math/ZSH.h
  math/stats/glm.h
  math/stats/permutation.h
  math/stats/typedefs.h
  surface/freesurfer.h
  surface/mesh.h
  surface/mesh_multi.h
  surface/polygon.h
  surface/scalar.h
  surface/types.h
  surface/utils.h
  surface/algo/image2mesh.h
  surface/algo/mesh2image.h
  surface/filter/base.h
  surface/filter/smooth.h
  surface/filter/vertex_transform.h)

set(${PROJECT_NAME}_SOURCES
  ../exec_version.cpp
  app.cpp
  axes.cpp
  bitset.cpp
  datatype.cpp
  exception.cpp
  header.cpp
  mrtrix.cpp
  phase_encoding.cpp
  progressbar.cpp
  signal_handler.cpp
  stats.cpp
  stride.cpp
  thread.cpp
  version.cpp
  adapter/reslice.cpp
  algo/histogram.cpp
  connectome/connectome.cpp
  connectome/enhance.cpp
  connectome/lut.cpp
  file/config.cpp
  file/json_utils.cpp
  file/key_value.cpp
  file/mgh.cpp
  file/mmap.cpp
  file/name_parser.cpp
  file/nifti1_utils.cpp
  file/nifti2_utils.cpp
  file/nifti_utils.cpp
  file/ofstream.cpp
  file/tiff.cpp
  file/dicom/dict.cpp
  file/dicom/element.cpp
  file/dicom/image.cpp
  file/dicom/mapper.cpp
  file/dicom/patient.cpp
  file/dicom/quick_scan.cpp
  file/dicom/select_cmdline.cpp
  file/dicom/series.cpp
  file/dicom/study.cpp
  file/dicom/tree.cpp
  formats/list.cpp
  formats/mrtrix_utils.cpp
  formats/analyse.cpp
  formats/dicom.cpp
  formats/list.cpp
  formats/list.h
  formats/mgh.cpp
  formats/mgz.cpp
  formats/mri.cpp
  formats/mrtrix.cpp
  formats/mrtrix_gz.cpp
  formats/mrtrix_sparse_legacy.cpp
  formats/mrtrix_utils.cpp
  formats/nifti1.cpp
  formats/nifti1_gz.cpp
  formats/nifti2.cpp
  formats/nifti2_gz.cpp
  formats/pipe.cpp
  formats/ram.cpp
  formats/tiff.cpp
  formats/xds.cpp
  image_io/base.cpp
  image_io/default.cpp
  image_io/fetch_store.cpp
  image_io/gz.cpp
  image_io/mosaic.cpp
  image_io/pipe.cpp
  image_io/ram.cpp
  image_io/scratch.cpp
  image_io/sparse.cpp
  image_io/tiff.cpp
  math/average_space.cpp
  math/SH.cpp
  math/stats/glm.cpp
  math/stats/permutation.cpp
  surface/algo/mesh2image.cpp
  surface/filter/base.cpp
  surface/filter/smooth.cpp
  surface/filter/vertex_transform.cpp
  surface/freesurfer.cpp
  surface/mesh.cpp
  surface/mesh_multi.cpp
  surface/polygon.cpp
  surface/scalar.cpp
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

target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads)
target_link_libraries(${PROJECT_NAME} PRIVATE ${ZLIB_LIBRARIES})

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
  COMPONENT core)

## ###################################################################
## Install rules - targets
## ###################################################################

install(TARGETS ${PROJECT_NAME}
  EXPORT mrtrix3Depends
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          COMPONENT   core)

######################################################################
### CMakeLists.txt ends here
