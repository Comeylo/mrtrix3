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

project(mrtrix3Gui)

## #################################################################
## Sources
## #################################################################

set(CMAKE_AUTORCC ON)

set(${PROJECT_NAME}_HEADERS
  color_button.h
  crosshair.h
  cursor.h
  gui.h
  lighting_dock.h
  projection.h
  dialog/dialog.h
  dialog/dicom.h
  dialog/file.h
  dialog/image_properties.h
  dialog/list.h
  dialog/opengl.h
  dialog/progress.h
  dialog/report_exception.h
  dwi/renderer.h
  dwi/render_frame.h
  icons/mrview.qrc
  opengl/font.h
  opengl/gl_core_3_3.h
  opengl/gl.h
  opengl/lighting.h
  opengl/shader.h
  opengl/transformation.h
  mrview/adjust_button.h
  mrview/colourmap_button.h
  mrview/colourmap.h
  mrview/combo_box_error.h
  mrview/displayable.h
  mrview/gui_image.h
  mrview/icons.h
  mrview/spin_box.h
  mrview/volume.h
  mrview/window.h
  mrview/mode/base.h
  mrview/mode/lightbox_gui.h
  mrview/mode/lightbox.h
  mrview/mode/list.h
  mrview/mode/ortho.h
  mrview/mode/slice.h
  mrview/mode/volume.h
  mrview/tool/base.h
  mrview/tool/list.h
  mrview/tool/list_model_base.h
  mrview/tool/overlay.h
  mrview/tool/screen_capture.h
  mrview/tool/view.h
  mrview/tool/connectome/colourmap_observers.h
  mrview/tool/connectome/connectome.h
  mrview/tool/connectome/edge.h
  mrview/tool/connectome/file_data_vector.h
  mrview/tool/connectome/matrix_list.h
  mrview/tool/connectome/node.h
  mrview/tool/connectome/node_list.h
  mrview/tool/connectome/node_overlay.h
  mrview/tool/connectome/selection.h
  mrview/tool/connectome/shaders.h
  mrview/tool/connectome/types.h
  mrview/tool/fixel/base_fixel.h
  mrview/tool/fixel/directory.h
  mrview/tool/fixel/fixel.h
  mrview/tool/fixel/image4D.h
  mrview/tool/fixel/vector_structs.h
  mrview/tool/odf/item.h
  mrview/tool/odf/model.h
  mrview/tool/odf/odf.h
  mrview/tool/odf/preview.h
  mrview/tool/odf/type.h
  mrview/tool/roi_editor/item.h
  mrview/tool/roi_editor/model.h
  mrview/tool/roi_editor/roi.h
  mrview/tool/roi_editor/undoentry.h
  mrview/tool/tractography/track_scalar_file.h
  mrview/tool/tractography/tractogram_enums.h
  mrview/tool/tractography/tractogram.h
  mrview/tool/tractography/tractography.h
  shapes/cube.h
  shapes/cylinder.h
  shapes/halfsphere.h
  shapes/sphere.h
  shview/icons.h
  shview/render_window.h
  )

set(${PROJECT_NAME}_SOURCES
  color_button.cpp
  crosshair.cpp
  cursor.cpp
  gui.cpp
  lighting_dock.cpp
  projection.cpp
  dialog/dialog.cpp
  dialog/dicom.cpp
  dialog/file.cpp
  dialog/image_properties.cpp
  dialog/list.cpp
  dialog/opengl.cpp
  dialog/progress.cpp
  dialog/report_exception.cpp
  dwi/renderer.cpp
  dwi/render_frame.cpp
  opengl/font.cpp
  opengl/gl_core_3_3.cpp
  opengl/gl.cpp
  opengl/lighting.cpp
  opengl/shader.cpp
  mrview/adjust_button.cpp
  mrview/colourmap_button.cpp
  mrview/colourmap.cpp
  mrview/combo_box_error.cpp
  mrview/displayable.cpp
  mrview/gui_image.cpp
  mrview/volume.cpp
  mrview/window.cpp
  mrview/mode/base.cpp
  mrview/mode/lightbox.cpp
  mrview/mode/ortho.cpp
  mrview/mode/slice.cpp
  mrview/mode/volume.cpp
  mrview/tool/base.cpp
  mrview/tool/overlay.cpp
  mrview/tool/screen_capture.cpp
  mrview/tool/view.cpp
  mrview/tool/connectome/colourmap_observers.cpp
  mrview/tool/connectome/connectome.cpp
  mrview/tool/connectome/edge.cpp
  mrview/tool/connectome/file_data_vector.cpp
  mrview/tool/connectome/matrix_list.cpp
  mrview/tool/connectome/node.cpp
  mrview/tool/connectome/node_list.cpp
  mrview/tool/connectome/node_overlay.cpp
  mrview/tool/connectome/selection.cpp
  mrview/tool/connectome/shaders.cpp
  mrview/tool/fixel/base_fixel.cpp
  mrview/tool/fixel/directory.cpp
  mrview/tool/fixel/fixel.cpp
  mrview/tool/fixel/image4D.cpp
  mrview/tool/fixel/legacy.cpp
  mrview/tool/odf/model.cpp
  mrview/tool/odf/item.cpp
  mrview/tool/odf/odf.cpp
  mrview/tool/odf/preview.cpp
  mrview/tool/roi_editor/item.cpp
  mrview/tool/roi_editor/model.cpp
  mrview/tool/roi_editor/roi.cpp
  mrview/tool/roi_editor/undoentry.cpp
  mrview/tool/tractography/track_scalar_file.cpp
  mrview/tool/tractography/tractogram.cpp
  mrview/tool/tractography/tractography.cpp
  shapes/cube.cpp
  shapes/cylinder.cpp
  shapes/halfsphere.cpp
  shapes/sphere.cpp
  shview/render_window.cpp
  )

## #################################################################
## Build rules
## #################################################################p

add_library(${PROJECT_NAME} SHARED
  ${${PROJECT_NAME}_SOURCES}
  ${${PROJECT_NAME}_HEADERS})

## ###################################################################
## Link rules
## ###################################################################
target_link_libraries(${PROJECT_NAME} mrtrix3Core)
target_link_libraries(${PROJECT_NAME} mrtrix3DWI)

target_link_libraries(${PROJECT_NAME} Qt5::Core)
target_link_libraries(${PROJECT_NAME} Qt5::OpenGL)
target_link_libraries(${PROJECT_NAME} Qt5::Widgets)

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
  COMPONENT gui)

## ###################################################################
## Install rules - targets
## ###################################################################

install(TARGETS ${PROJECT_NAME}
  EXPORT mrtrix3Depends
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          COMPONENT   gui)

######################################################################
### CMakeLists.txt ends here
