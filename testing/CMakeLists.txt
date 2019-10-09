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

## #################################################################
## Build rules
## #################################################################

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(testing_diff_dir cmd/testing_diff_dir.cpp)
target_link_libraries(testing_diff_dir mrtrix3Core)

add_executable(testing_diff_fixel cmd/testing_diff_fixel.cpp)
target_link_libraries(testing_diff_fixel mrtrix3Core)

add_executable(testing_diff_fixel_old cmd/testing_diff_fixel_old.cpp)
target_link_libraries(testing_diff_fixel_old mrtrix3Core)

add_executable(testing_diff_image cmd/testing_diff_image.cpp)
target_link_libraries(testing_diff_image mrtrix3Core)

add_executable(testing_diff_matrix cmd/testing_diff_matrix.cpp)
target_link_libraries(testing_diff_matrix mrtrix3Core)

add_executable(testing_diff_mesh cmd/testing_diff_mesh.cpp)
target_link_libraries(testing_diff_mesh mrtrix3Core)

add_executable(testing_diff_peaks cmd/testing_diff_peaks.cpp)
target_link_libraries(testing_diff_peaks mrtrix3Core)

add_executable(testing_diff_tck cmd/testing_diff_tck.cpp)
target_link_libraries(testing_diff_tck mrtrix3Core)
target_link_libraries(testing_diff_tck mrtrix3DWI)

add_executable(testing_diff_tsf cmd/testing_diff_tsf.cpp)
target_link_libraries(testing_diff_tsf mrtrix3Core)
target_link_libraries(testing_diff_tsf mrtrix3DWI)

add_executable(testing_gen_data cmd/testing_gen_data.cpp)
target_link_libraries(testing_gen_data mrtrix3Core)

## ###################################################################
## Link rules
## ###################################################################

######################################################################
### CMakeLists.txt ends here
