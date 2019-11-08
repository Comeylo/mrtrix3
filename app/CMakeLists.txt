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

add_executable(5tt2gmwmi cmd/5tt2gmwmi.cpp)
target_link_libraries(5tt2gmwmi mrtrix3Core)
target_link_libraries(5tt2gmwmi mrtrix3DWI)
target_link_libraries(5tt2gmwmi Threads::Threads)

add_executable(5tt2vis cmd/5tt2vis.cpp)
target_link_libraries(5tt2vis mrtrix3Core)
target_link_libraries(5tt2vis mrtrix3DWI)
target_link_libraries(5tt2vis Threads::Threads)

add_executable(5ttcheck cmd/5ttcheck.cpp)
target_link_libraries(5ttcheck mrtrix3Core)
target_link_libraries(5ttcheck mrtrix3DWI)
target_link_libraries(5ttcheck Threads::Threads)

add_executable(5ttedit cmd/5ttedit.cpp)
target_link_libraries(5ttedit mrtrix3Core)
target_link_libraries(5ttedit mrtrix3DWI)
target_link_libraries(5ttedit Threads::Threads)

add_executable(afdconnectivity cmd/afdconnectivity.cpp)
target_link_libraries(afdconnectivity mrtrix3Core)
target_link_libraries(afdconnectivity mrtrix3DWI)
target_link_libraries(afdconnectivity Threads::Threads)

add_executable(amp2response cmd/amp2response.cpp)
target_link_libraries(amp2response mrtrix3Core)
target_link_libraries(amp2response mrtrix3DWI)
target_link_libraries(amp2response Threads::Threads)

add_executable(amp2sh cmd/amp2sh.cpp)
target_link_libraries(amp2sh mrtrix3Core)
target_link_libraries(amp2sh mrtrix3DWI)
target_link_libraries(amp2sh Threads::Threads)

add_executable(connectome2tck cmd/connectome2tck.cpp)
target_link_libraries(connectome2tck mrtrix3Core)
target_link_libraries(connectome2tck mrtrix3DWI)
target_link_libraries(connectome2tck Threads::Threads)

add_executable(connectomestats cmd/connectomestats.cpp)
target_link_libraries(connectomestats mrtrix3Core)
target_link_libraries(connectomestats mrtrix3Stats)
target_link_libraries(connectomestats Threads::Threads)

add_executable(dcmedit cmd/dcmedit.cpp)
target_link_libraries(dcmedit mrtrix3Core)

add_executable(dcminfo cmd/dcminfo.cpp)
target_link_libraries(dcminfo mrtrix3Core)

add_executable(dirflip cmd/dirflip.cpp)
target_link_libraries(dirflip mrtrix3Core)
target_link_libraries(dirflip mrtrix3DWI)
target_link_libraries(dirflip Threads::Threads)

add_executable(dirgen cmd/dirgen.cpp)
target_link_libraries(dirgen mrtrix3Core)
target_link_libraries(dirgen mrtrix3DWI)
target_link_libraries(dirgen Threads::Threads)

add_executable(dirmerge cmd/dirmerge.cpp)
target_link_libraries(dirmerge mrtrix3Core)
target_link_libraries(dirmerge mrtrix3DWI)

add_executable(dirorder cmd/dirorder.cpp)
target_link_libraries(dirorder mrtrix3Core)
target_link_libraries(dirorder mrtrix3DWI)

add_executable(dirsplit cmd/dirsplit.cpp)
target_link_libraries(dirsplit mrtrix3Core)
target_link_libraries(dirsplit mrtrix3DWI)
target_link_libraries(dirsplit Threads::Threads)

add_executable(dirstat cmd/dirstat.cpp)
target_link_libraries(dirstat mrtrix3Core)
target_link_libraries(dirstat mrtrix3DWI)

add_executable(dwi2adc cmd/dwi2adc.cpp)
target_link_libraries(dwi2adc mrtrix3Core)
target_link_libraries(dwi2adc mrtrix3DWI)
target_link_libraries(dwi2adc Threads::Threads)

add_executable(dwi2fod cmd/dwi2fod.cpp)
target_link_libraries(dwi2fod mrtrix3Core)
target_link_libraries(dwi2fod mrtrix3DWI)
target_link_libraries(dwi2fod Threads::Threads)

add_executable(dwi2mask cmd/dwi2mask.cpp)
target_link_libraries(dwi2mask mrtrix3Core)
target_link_libraries(dwi2mask mrtrix3DWI)
target_link_libraries(dwi2mask Threads::Threads)

add_executable(dwi2tensor cmd/dwi2tensor.cpp)
target_link_libraries(dwi2tensor mrtrix3Core)
target_link_libraries(dwi2tensor mrtrix3DWI)
target_link_libraries(dwi2tensor Threads::Threads)

add_executable(dwidenoise cmd/dwidenoise.cpp)
target_link_libraries(dwidenoise mrtrix3Core)
target_link_libraries(dwidenoise mrtrix3DWI)
target_link_libraries(dwidenoise Threads::Threads)

add_executable(dwiextract cmd/dwiextract.cpp)
target_link_libraries(dwiextract mrtrix3Core)
target_link_libraries(dwiextract mrtrix3DWI)
target_link_libraries(dwiextract Threads::Threads)

add_executable(dwinormalise cmd/dwinormalise.cpp)
target_link_libraries(dwinormalise mrtrix3Core)
target_link_libraries(dwinormalise mrtrix3DWI)
target_link_libraries(dwinormalise Threads::Threads)

add_executable(fixel2sh cmd/fixel2sh.cpp)
target_link_libraries(fixel2sh mrtrix3Core)
target_link_libraries(fixel2sh Threads::Threads)

add_executable(fixel2tsf cmd/fixel2tsf.cpp)
target_link_libraries(fixel2tsf mrtrix3Core)
target_link_libraries(fixel2tsf mrtrix3DWI)
target_link_libraries(fixel2tsf Threads::Threads)

add_executable(fixel2voxel cmd/fixel2voxel.cpp)
target_link_libraries(fixel2voxel mrtrix3Core)
target_link_libraries(fixel2voxel Threads::Threads)

add_executable(fixelcfestats cmd/fixelcfestats.cpp)
target_link_libraries(fixelcfestats mrtrix3Core)
target_link_libraries(fixelcfestats mrtrix3DWI)
target_link_libraries(fixelcfestats mrtrix3Stats)
target_link_libraries(fixelcfestats Threads::Threads)

add_executable(fixelconvert cmd/fixelconvert.cpp)
target_link_libraries(fixelconvert mrtrix3Core)
target_link_libraries(fixelconvert Threads::Threads)

add_executable(fixelcorrespondence cmd/fixelcorrespondence.cpp)
target_link_libraries(fixelcorrespondence mrtrix3Core)
target_link_libraries(fixelcorrespondence Threads::Threads)

add_executable(fixelcrop cmd/fixelcrop.cpp)
target_link_libraries(fixelcrop mrtrix3Core)
target_link_libraries(fixelcrop Threads::Threads)

add_executable(fixelreorient cmd/fixelreorient.cpp)
target_link_libraries(fixelreorient mrtrix3Core)
target_link_libraries(fixelreorient Threads::Threads)

add_executable(fod2dec cmd/fod2dec.cpp)
target_link_libraries(fod2dec mrtrix3Core)
target_link_libraries(fod2dec mrtrix3DWI)
target_link_libraries(fod2dec Threads::Threads)

add_executable(fod2fixel cmd/fod2fixel.cpp)
target_link_libraries(fod2fixel mrtrix3Core)
target_link_libraries(fod2fixel mrtrix3DWI)
target_link_libraries(fod2fixel Threads::Threads)

add_executable(label2colour cmd/label2colour.cpp)
target_link_libraries(label2colour mrtrix3Core)
target_link_libraries(label2colour Threads::Threads)

add_executable(label2mesh cmd/label2mesh.cpp)
target_link_libraries(label2mesh mrtrix3Core)
target_link_libraries(label2mesh Threads::Threads)

add_executable(labelconvert cmd/labelconvert.cpp)
target_link_libraries(labelconvert mrtrix3Core)
target_link_libraries(labelconvert Threads::Threads)

add_executable(maskdump cmd/maskdump.cpp)
target_link_libraries(maskdump mrtrix3Core)
target_link_libraries(maskdump Threads::Threads)

add_executable(maskfilter cmd/maskfilter.cpp)
target_link_libraries(maskfilter mrtrix3Core)
target_link_libraries(maskfilter Threads::Threads)

add_executable(mesh2voxel cmd/mesh2voxel.cpp)
target_link_libraries(mesh2voxel mrtrix3Core)
target_link_libraries(mesh2voxel Threads::Threads)

add_executable(meshconvert cmd/meshconvert.cpp)
target_link_libraries(meshconvert mrtrix3Core)

add_executable(meshfilter cmd/meshfilter.cpp)
target_link_libraries(meshfilter mrtrix3Core)

add_executable(mraverageheader cmd/mraverageheader.cpp)
target_link_libraries(mraverageheader mrtrix3Core)
target_link_libraries(mraverageheader Threads::Threads)

add_executable(mrcalc cmd/mrcalc.cpp)
target_link_libraries(mrcalc mrtrix3Core)
target_link_libraries(mrcalc mrtrix3DWI)
target_link_libraries(mrcalc Threads::Threads)

add_executable(mrcat cmd/mrcat.cpp)
target_link_libraries(mrcat mrtrix3Core)
target_link_libraries(mrcat mrtrix3DWI)
target_link_libraries(mrcat Threads::Threads)

add_executable(mrcheckerboardmask cmd/mrcheckerboardmask.cpp)
target_link_libraries(mrcheckerboardmask mrtrix3Core)
target_link_libraries(mrcheckerboardmask Threads::Threads)

add_executable(mrclusterstats cmd/mrclusterstats.cpp)
target_link_libraries(mrclusterstats mrtrix3Core)
target_link_libraries(mrclusterstats mrtrix3Stats)
target_link_libraries(mrclusterstats Threads::Threads)

add_executable(mrconvert cmd/mrconvert.cpp)
target_link_libraries(mrconvert mrtrix3Core)
target_link_libraries(mrconvert mrtrix3DWI)
target_link_libraries(mrconvert Threads::Threads)

add_executable(mrcrop cmd/mrcrop.cpp)
target_link_libraries(mrcrop mrtrix3Core)
target_link_libraries(mrcrop Threads::Threads)

add_executable(mrdegibbs cmd/mrdegibbs.cpp)
target_link_libraries(mrdegibbs mrtrix3Core)
target_link_libraries(mrdegibbs Threads::Threads)

add_executable(mrdump cmd/mrdump.cpp)
target_link_libraries(mrdump mrtrix3Core)
target_link_libraries(mrdump Threads::Threads)

add_executable(mredit cmd/mredit.cpp)
target_link_libraries(mredit mrtrix3Core)
target_link_libraries(mredit Threads::Threads)

add_executable(mrfilter cmd/mrfilter.cpp)
target_link_libraries(mrfilter mrtrix3Core)
target_link_libraries(mrfilter Threads::Threads)

add_executable(mrhistmatch cmd/mrhistmatch.cpp)
target_link_libraries(mrhistmatch mrtrix3Core)
target_link_libraries(mrhistmatch Threads::Threads)

add_executable(mrhistogram cmd/mrhistogram.cpp)
target_link_libraries(mrhistogram mrtrix3Core)
target_link_libraries(mrhistogram Threads::Threads)

add_executable(mrinfo cmd/mrinfo.cpp)
target_link_libraries(mrinfo mrtrix3Core)
target_link_libraries(mrinfo mrtrix3DWI)

add_executable(mrmath cmd/mrmath.cpp)
target_link_libraries(mrmath mrtrix3Core)
target_link_libraries(mrmath mrtrix3DWI)
target_link_libraries(mrmath Threads::Threads)

add_executable(mrmetric cmd/mrmetric.cpp)
target_link_libraries(mrmetric mrtrix3Core)
target_link_libraries(mrmetric Threads::Threads)

add_executable(mrpad cmd/mrpad.cpp)
target_link_libraries(mrpad mrtrix3Core)
target_link_libraries(mrpad Threads::Threads)

add_executable(mrregister cmd/mrregister.cpp)
target_link_libraries(mrregister mrtrix3Core)
target_link_libraries(mrregister mrtrix3DWI)
target_link_libraries(mrregister mrtrix3Registration)
target_link_libraries(mrregister Threads::Threads)

add_executable(mrresize cmd/mrresize.cpp)
target_link_libraries(mrresize mrtrix3Core)
target_link_libraries(mrresize Threads::Threads)

add_executable(mrstats cmd/mrstats.cpp)
target_link_libraries(mrstats mrtrix3Core)
target_link_libraries(mrstats mrtrix3Stats)
target_link_libraries(mrstats Threads::Threads)

add_executable(mrthreshold cmd/mrthreshold.cpp)
target_link_libraries(mrthreshold mrtrix3Core)
target_link_libraries(mrthreshold Threads::Threads)

add_executable(mrtransform cmd/mrtransform.cpp)
target_link_libraries(mrtransform mrtrix3Core)
target_link_libraries(mrtransform mrtrix3DWI)
target_link_libraries(mrtransform Threads::Threads)

add_executable(mrview cmd/mrview.cpp)
target_link_libraries(mrview mrtrix3Core)
target_link_libraries(mrview mrtrix3Gui)

add_executable(mtnormalise cmd/mtnormalise.cpp)
target_link_libraries(mtnormalise mrtrix3Core)
target_link_libraries(mtnormalise Threads::Threads)

add_executable(peaks2amp cmd/peaks2amp.cpp)
target_link_libraries(peaks2amp mrtrix3Core)
target_link_libraries(peaks2amp Threads::Threads)

add_executable(sh2amp cmd/sh2amp.cpp)
target_link_libraries(sh2amp mrtrix3Core)
target_link_libraries(sh2amp mrtrix3DWI)
target_link_libraries(sh2amp Threads::Threads)

add_executable(sh2peaks cmd/sh2peaks.cpp)
target_link_libraries(sh2peaks mrtrix3Core)
target_link_libraries(sh2peaks Threads::Threads)

add_executable(sh2power cmd/sh2power.cpp)
target_link_libraries(sh2power mrtrix3Core)
target_link_libraries(sh2power Threads::Threads)

add_executable(sh2response cmd/sh2response.cpp)
target_link_libraries(sh2response mrtrix3Core)
target_link_libraries(sh2response Threads::Threads)

add_executable(shbasis cmd/shbasis.cpp)
target_link_libraries(shbasis mrtrix3Core)
target_link_libraries(shbasis Threads::Threads)

add_executable(shconv cmd/shconv.cpp)
target_link_libraries(shconv mrtrix3Core)
target_link_libraries(shconv Threads::Threads)

add_executable(shview cmd/shview.cpp)
target_link_libraries(shview mrtrix3Core)
target_link_libraries(shview mrtrix3Gui)

add_executable(tck2connectome cmd/tck2connectome.cpp)
target_link_libraries(tck2connectome mrtrix3Core)
target_link_libraries(tck2connectome mrtrix3DWI)
target_link_libraries(tck2connectome Threads::Threads)

add_executable(tck2fixel cmd/tck2fixel.cpp)
target_link_libraries(tck2fixel mrtrix3Core)
target_link_libraries(tck2fixel mrtrix3DWI)
target_link_libraries(tck2fixel Threads::Threads)

add_executable(tckconvert cmd/tckconvert.cpp)
target_link_libraries(tckconvert mrtrix3Core)
target_link_libraries(tckconvert mrtrix3DWI)

add_executable(tckdfc cmd/tckdfc.cpp)
target_link_libraries(tckdfc mrtrix3Core)
target_link_libraries(tckdfc mrtrix3DWI)
target_link_libraries(tckdfc Threads::Threads)

add_executable(tckedit cmd/tckedit.cpp)
target_link_libraries(tckedit mrtrix3Core)
target_link_libraries(tckedit mrtrix3DWI)
target_link_libraries(tckedit Threads::Threads)

add_executable(tckgen cmd/tckgen.cpp)
target_link_libraries(tckgen mrtrix3Core)
target_link_libraries(tckgen mrtrix3DWI)
target_link_libraries(tckgen Threads::Threads)

add_executable(tckglobal cmd/tckglobal.cpp)
target_link_libraries(tckglobal mrtrix3Core)
target_link_libraries(tckglobal mrtrix3DWI)
target_link_libraries(tckglobal Threads::Threads)

add_executable(tckinfo cmd/tckinfo.cpp)
target_link_libraries(tckinfo mrtrix3Core)
target_link_libraries(tckinfo mrtrix3DWI)

add_executable(tckmap cmd/tckmap.cpp)
target_link_libraries(tckmap mrtrix3Core)
target_link_libraries(tckmap mrtrix3DWI)
target_link_libraries(tckmap Threads::Threads)

add_executable(tckresample cmd/tckresample.cpp)
target_link_libraries(tckresample mrtrix3Core)
target_link_libraries(tckresample mrtrix3DWI)
target_link_libraries(tckresample Threads::Threads)

add_executable(tcksample cmd/tcksample.cpp)
target_link_libraries(tcksample mrtrix3Core)
target_link_libraries(tcksample mrtrix3DWI)
target_link_libraries(tcksample Threads::Threads)

add_executable(tcksift cmd/tcksift.cpp)
target_link_libraries(tcksift mrtrix3Core)
target_link_libraries(tcksift mrtrix3DWI)
target_link_libraries(tcksift Threads::Threads)

add_executable(tcksift2 cmd/tcksift2.cpp)
target_link_libraries(tcksift2 mrtrix3Core)
target_link_libraries(tcksift2 mrtrix3DWI)
target_link_libraries(tcksift2 Threads::Threads)

add_executable(tckstats cmd/tckstats.cpp)
target_link_libraries(tckstats mrtrix3Core)
target_link_libraries(tckstats mrtrix3DWI)
target_link_libraries(tckstats mrtrix3Stats)

add_executable(tcktransform cmd/tcktransform.cpp)
target_link_libraries(tcktransform mrtrix3Core)
target_link_libraries(tcktransform mrtrix3DWI)
target_link_libraries(tcktransform Threads::Threads)

add_executable(tensor2metric cmd/tensor2metric.cpp)
target_link_libraries(tensor2metric mrtrix3Core)
target_link_libraries(tensor2metric Threads::Threads)

add_executable(transformcalc cmd/transformcalc.cpp)
target_link_libraries(transformcalc mrtrix3Core)

add_executable(transformcompose cmd/transformcompose.cpp)
target_link_libraries(transformcompose mrtrix3Core)
target_link_libraries(transformcompose Threads::Threads)

add_executable(transformconvert cmd/transformconvert.cpp)
target_link_libraries(transformconvert mrtrix3Core)

add_executable(tsfdivide cmd/tsfdivide.cpp)
target_link_libraries(tsfdivide mrtrix3Core)
target_link_libraries(tsfdivide mrtrix3DWI)

add_executable(tsfinfo cmd/tsfinfo.cpp)
target_link_libraries(tsfinfo mrtrix3Core)
target_link_libraries(tsfinfo mrtrix3DWI)

add_executable(tsfmult cmd/tsfmult.cpp)
target_link_libraries(tsfmult mrtrix3Core)
target_link_libraries(tsfmult mrtrix3DWI)

add_executable(tsfsmooth cmd/tsfsmooth.cpp)
target_link_libraries(tsfsmooth mrtrix3Core)
target_link_libraries(tsfsmooth mrtrix3DWI)

add_executable(tsfthreshold cmd/tsfthreshold.cpp)
target_link_libraries(tsfthreshold mrtrix3Core)
target_link_libraries(tsfthreshold mrtrix3DWI)

add_executable(tsfvalidate cmd/tsfvalidate.cpp)
target_link_libraries(tsfvalidate mrtrix3Core)
target_link_libraries(tsfvalidate mrtrix3DWI)

add_executable(vectorstats cmd/vectorstats.cpp)
target_link_libraries(vectorstats mrtrix3Core)
target_link_libraries(vectorstats mrtrix3Stats)
target_link_libraries(vectorstats Threads::Threads)

add_executable(voxel2fixel cmd/voxel2fixel.cpp)
target_link_libraries(voxel2fixel mrtrix3Core)
target_link_libraries(voxel2fixel Threads::Threads)

add_executable(voxel2mesh cmd/voxel2mesh.cpp)
target_link_libraries(voxel2mesh mrtrix3Core)
target_link_libraries(voxel2mesh Threads::Threads)

add_executable(warp2metric cmd/warp2metric.cpp)
target_link_libraries(warp2metric mrtrix3Core)
target_link_libraries(warp2metric Threads::Threads)

add_executable(warpconvert cmd/warpconvert.cpp)
target_link_libraries(warpconvert mrtrix3Core)
target_link_libraries(warpconvert Threads::Threads)

add_executable(warpcorrect cmd/warpcorrect.cpp)
target_link_libraries(warpcorrect mrtrix3Core)
target_link_libraries(warpcorrect Threads::Threads)

add_executable(warpinit cmd/warpinit.cpp)
target_link_libraries(warpinit mrtrix3Core)
target_link_libraries(warpinit Threads::Threads)

add_executable(warpinvert cmd/warpinvert.cpp)
target_link_libraries(warpinvert mrtrix3Core)
target_link_libraries(warpinvert Threads::Threads)

## #################################################################
## Install rules
## #################################################################

install(TARGETS mrview
   BUNDLE DESTINATION bin
   RUNTIME DESTINATION bin)

######################################################################
### CMakeLists.txt ends here
