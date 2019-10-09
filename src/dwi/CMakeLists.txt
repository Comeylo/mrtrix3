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

project(mrtrix3DWI)

## #################################################################
## Sources
## #################################################################

set(${PROJECT_NAME}_HEADERS
  bootstrap.h
  gradient.h
  fmls.h
  shells.h
  tensor.h
  directions/file.h
  directions/mask.h
  directions/predefined.h
  directions/set.h
  sdeconv/csd.h
  sdeconv/msmt_csd.h
  tractography/file_base.h
  tractography/file.h
  tractography/properties.h
  tractography/rng.h
  tractography/roi.h
  tractography/scalar_file.h
  tractography/streamline.h
  tractography/weights.h
  tractography/ACT/act.h
  tractography/ACT/gmwmi.h
  tractography/ACT/method.h
  tractography/ACT/shared.h
  tractography/ACT/tissues.h
  tractography/algorithms/calibrator.h
  tractography/algorithms/fact.h
  tractography/algorithms/iFOD1.h
  tractography/algorithms/iFOD2.h
  tractography/algorithms/nulldist.h
  tractography/algorithms/sd_stream.h
  tractography/algorithms/seedtest.h
  tractography/algorithms/tensor_det.h
  tractography/algorithms/tensor_prob.h
  tractography/connectome/connectome.h
  tractography/connectome/exemplar.h
  tractography/connectome/extract.h
  tractography/connectome/mapped_track.h
  tractography/connectome/mapper.h
  tractography/connectome/matrix.h
  tractography/connectome/metric.h
  tractography/connectome/streamline.h
  tractography/connectome/tck2nodes.h
  tractography/editing/editing.h
  tractography/editing/loader.h
  tractography/editing/receiver.h
  tractography/editing/worker.h
  tractography/GT/energy.h
  tractography/GT/externalenergy.h
  tractography/GT/gt.h
  tractography/GT/internalenergy.h
  tractography/GT/mhsampler.h
  tractography/GT/particlegrid.h
  tractography/GT/particle.h
  tractography/GT/particlepool.h
  tractography/GT/spatiallock.h
  tractography/mapping/buffer_scratch_dump.h
  tractography/mapping/fixel_td_map.h
  tractography/mapping/loader.h
  tractography/mapping/mapper.h
  tractography/mapping/mapper_plugins.h
  tractography/mapping/mapping.h
  tractography/mapping/twi_stats.h
  tractography/mapping/voxel.h
  tractography/mapping/writer.h
  tractography/mapping/gaussian/mapper.h
  tractography/mapping/gaussian/voxel.h
  tractography/resampling/arc.h
  tractography/resampling/downsampler.h
  tractography/resampling/endpoints.h
  tractography/resampling/fixed_num_points.h
  tractography/resampling/fixed_step_size.h
  tractography/resampling/resampling.h
  tractography/resampling/upsampler.h
  tractography/seeding/base.h
  tractography/seeding/basic.h
  tractography/seeding/dynamic.h
  tractography/seeding/gmwmi.h
  tractography/seeding/list.h
  tractography/seeding/seeding.h
  tractography/SIFT/fixel.h
  tractography/SIFT/gradient_sort.h
  tractography/SIFT/model_base.h
  tractography/SIFT/model.h
  tractography/SIFT/output.h
  tractography/SIFT/proc_mask.h
  tractography/SIFT/sifter.h
  tractography/SIFT/sift.h
  tractography/SIFT/track_contribution.h
  tractography/SIFT/track_index_range.h
  tractography/SIFT/types.h
  tractography/SIFT2/coeff_optimiser.h
  tractography/SIFT2/fixel.h
  tractography/SIFT2/fixel_updater.h
  tractography/SIFT2/line_search.h
  tractography/SIFT2/reg_calculator.h
  tractography/SIFT2/regularisation.h
  tractography/SIFT2/streamline_stats.h
  tractography/SIFT2/tckfactor.h
  tractography/tracking/early_exit.h
  tractography/tracking/exec.h
  tractography/tracking/generated_track.h
  tractography/tracking/method.h
  tractography/tracking/shared.h
  tractography/tracking/tractography.h
  tractography/tracking/types.h
  tractography/tracking/write_kernel.h
)

set(${PROJECT_NAME}_SOURCES
  gradient.cpp
  fmls.cpp
  shells.cpp
  directions/file.cpp
  directions/mask.cpp
  directions/predefined.cpp
  directions/set.cpp
  sdeconv/csd.cpp
  tractography/file_base.cpp
  tractography/properties.cpp
  tractography/rng.cpp
  tractography/roi.cpp
  tractography/weights.cpp
  tractography/algorithms/iFOD2.cpp
  tractography/ACT/gmwmi.cpp
  tractography/ACT/act.cpp
  tractography/connectome/connectome.cpp
  tractography/connectome/exemplar.cpp
  tractography/connectome/extract.cpp
  tractography/connectome/matrix.cpp
  tractography/connectome/tck2nodes.cpp
  tractography/editing/editing.cpp
  tractography/editing/receiver.cpp
  tractography/editing/worker.cpp
  tractography/GT/externalenergy.cpp
  tractography/GT/gt.cpp
  tractography/GT/internalenergy.cpp
  tractography/GT/mhsampler.cpp
  tractography/GT/particle.cpp
  tractography/GT/particlegrid.cpp
  tractography/mapping/fixel_td_map.cpp
  tractography/mapping/mapper.cpp
  tractography/mapping/mapper_plugins.cpp
  tractography/mapping/mapping.cpp
  tractography/mapping/twi_stats.cpp
  tractography/mapping/voxel.cpp
  tractography/mapping/writer.cpp
  tractography/mapping/gaussian/mapper.cpp
  tractography/resampling/arc.cpp
  tractography/resampling/downsampler.cpp
  tractography/resampling/endpoints.cpp
  tractography/resampling/fixed_num_points.cpp
  tractography/resampling/fixed_step_size.cpp
  tractography/resampling/resampling.cpp
  tractography/resampling/upsampler.cpp
  tractography/seeding/basic.cpp
  tractography/seeding/dynamic.cpp
  tractography/seeding/gmwmi.cpp
  tractography/seeding/list.cpp
  tractography/seeding/seeding.cpp
  tractography/SIFT/gradient_sort.cpp
  tractography/SIFT/proc_mask.cpp
  tractography/SIFT/sift.cpp
  tractography/SIFT/sifter.cpp
  tractography/SIFT/track_contribution.cpp
  tractography/SIFT/track_index_range.cpp
  tractography/SIFT2/coeff_optimiser.cpp
  tractography/SIFT2/fixel_updater.cpp
  tractography/SIFT2/line_search.cpp
  tractography/SIFT2/reg_calculator.cpp
  tractography/SIFT2/streamline_stats.cpp
  tractography/SIFT2/tckfactor.cpp
  tractography/tracking/early_exit.cpp
  tractography/tracking/method.cpp
  tractography/tracking/shared.cpp
  tractography/tracking/tractography.cpp
  tractography/tracking/write_kernel.cpp
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

target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)

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
  COMPONENT dwi)

## ###################################################################
## Install rules - targets
## ###################################################################

install(TARGETS ${PROJECT_NAME}
  EXPORT mrtrix3Depends
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          COMPONENT   dwi)

######################################################################
### CMakeLists.txt ends here
