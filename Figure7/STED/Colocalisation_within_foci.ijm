folder = "/home/sasha/Science/PhD/Laue_lab/Microscopy/STED/HP1b_DNA_differentiation_Devina/";
imgdir = "STED/";
filename = "dna-JF646_HP1-af594_48h_cell";
savename = "DNA_HP1_48h_cell";
for (f = 1; f < 8; f++) {
	open(folder + imgdir + filename + d2s(f, 0) + "_STED_Gauss1.tif");
	run("Split Channels");
	selectImage("C1-" + filename + d2s(f, 0) + "_STED_Gauss1.tif");
	resetMinAndMax();
	run("Apply LUT");
	rename("Channel1");
	selectImage("C2-" + filename + d2s(f, 0) + "_STED_Gauss1.tif");
	resetMinAndMax();
	run("Apply LUT");
	rename("Channel2");
	imageCalculator("Add create 32-bit", "Channel1","Channel2");
	rename("Sum");
	setAutoThreshold("MaxEntropy dark no-reset");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	rename("Mask");
	run("Analyze Particles...", "size=100-Infinity pixel show=Masks clear");
	rename("Mask2");
	run("Invert");
	waitForUser;
	selectImage("Mask2");
	run("Dilate");
	run("Dilate");
	run("Dilate");
	run("Dilate");
	run("Dilate");
	run("Fill Holes");
	saveAs("tiff", folder + "Colocalisation_analysis/Masks/" + savename + d2s(f, 0) + ".tif");
	run("Analyze Particles...", "size=100-Infinity pixel show=Nothing clear add");
	roiManager("Save", folder + "Colocalisation_analysis/Masks/" + savename + d2s(f, 0) + "_roi.zip");
	nROIs = roiManager("count");
	R = "";
	for (i=0; i<nROIs; i++){
		selectWindow("Channel1");
		roiManager("Select", i);
		run("Coloc 2", "channel_1=Channel1" + " channel_2=Channel2" +
		" roi_or_mask=[ROI(s) in channel 1] threshold_regression=Bisection psf=3 costes_randomisations=10");
		string = getInfo("log");
		R = R + d2s(i, 0) +" " + substring(string, indexOf(string, "Pearson's R value (no threshold), ")+lengthOf("Pearson's R value (no threshold), "),
		indexOf(string, "Pearson's R value (below threshold)"));
		print("\\Clear");
	}
	File.saveString(R, folder + "Colocalisation_analysis/Coloc/" + savename + d2s(f, 0) + "_coloc.txt");
	
	selectWindow("ROI Manager");
	run("Close");
	run("Close All");
}

