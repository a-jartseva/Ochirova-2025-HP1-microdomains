image_basename = "90X_HILO_488_OD4_640narrow_OD0p8_3ms_1_MMStack_Pos0_DNA";

dir_raw = "/home/sasha/PhD/Laue_lab/Microscopy/HP1b_live_blinkingDyes/" +
"20221018_Halo-JF639b_Hoechst-SPY505/TopHat_segmentation/PNGs_ManThreshold_raw/";
dir_edited = "/home/sasha/PhD/Laue_lab/Microscopy/HP1b_live_blinkingDyes/" +
"20221018_Halo-JF639b_Hoechst-SPY505/TopHat_segmentation/PNGs_ManThreshold_edited/";
dir_pseudo = "/home/sasha/PhD/Laue_lab/Microscopy/HP1b_live_blinkingDyes/" +
"20221018_Halo-JF639b_Hoechst-SPY505/TopHat_segmentation/PNGs_pseudo/";

////Scale up
run("Scale...", "x=10 y=10 z=1.0 interpolation=Bilinear average process create");

////Run segmentation
image = image_basename + "-1.tif";
setBatchMode(true);
n = nSlices();
for (i=1; i<n+1; i++) {
	selectWindow(image);
	setSlice(i);
	run("Duplicate...", " ");
	rename("Dup");
	run("Morphological Filters", "operation=[White Top Hat] element=Square radius=30");
	rename("Tophat" + d2s(i, 0));
	close("Dup");
}
run("Images to Stack", "name=Tophat_stack title=Tophat use");
setBatchMode(false);
selectWindow("Tophat_stack");
run("Threshold...");
waitForUser;
//run("Auto Threshold", "method=Yen white use_stack_histogram");
run("Analyze Particles...", "size=70-Infinity show=Masks exclude clear stack");
run("Invert LUT");
rename("Mask");
close("Tophat_stack");

////Save as png
run("Image Sequence... ", "dir=" + dir_raw + " format=PNG" +
" name=" + image_basename + "_foci digits=2");

//for (i=1; i<15; i++) {
//	setSlice(i);
//	setForegroundColor(255,255,255);
//	run("Fill", "slice");
//}

////Edit manually and save
waitForUser;
run("Image Sequence... ", "dir=" + dir_edited + " format=PNG" +
" name=" + image_basename + "_foci digits=2");

////Colour foci
//Duplicate FinalMask -> RGB
//run("Analyze Particles...", "exclude clear add stack"); //in the 8-bit mask
//image = "Mask-1";
//nROIs=roiManager("count");
//ROIs = newArray();
//for(r=0; r<nROIs; r++){
//	selectWindow(image);
//  roiManager("select", r);
//	setForegroundColor(0, 255, 255);
//	run("Fill", "slice");
//	ROIs = Array.concat(ROIs, r);
//}
//roiManager("Select", ROIs);
//roiManager("Delete");

//run("Image Sequence... ", "dir=/home/sasha/PhD/Laue_lab/Microscopy/HP1b_live_blinkingDyes/" +
//"20221018_Halo-JF639b_Hoechst-SPY505/TopHat_segmentation/PNGs_ManThreshold_coloured/ format=PNG" +
//" name=Pos5_foci digits=2");

//Pseudo mask
selectImage("Mask");
run("Analyze Particles...", "exclude clear add stack"); //in the 8-bit mask
//Manually translate the ROIs
waitForUser;
newImage("PseudoMask", "8-bit black", 1040, 1180, 8);
image = "PseudoMask";
nROIs=roiManager("count");
ROIs = newArray();
for(r=0; r<nROIs; r++){
	selectWindow(image);
    roiManager("select", r);
	setForegroundColor(255, 255, 255);
	run("Fill", "slice");
	ROIs = Array.concat(ROIs, r);
}
roiManager("Select", ROIs);
roiManager("Delete");

//Pseudo mask
run("Image Sequence... ", "dir=" + dir_pseudo + " format=PNG" +
" name=" + image_basename + "_foci digits=2");


run("Close All");
