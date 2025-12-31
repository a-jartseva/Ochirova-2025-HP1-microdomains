folder = "/home/sasha/Science/PhD/Laue_lab/Microscopy/Confocal/Histone_marks/20220201_HP1b-JF646_Hoechst_H3K9me3-Cy3_CAICAiry/";
basename = "20220201_Pos";
for (f = 1; f < 8; f++) {
	//f=1;
	open(folder + basename + f + "_Airy.czi");
	run("Split Channels");
	selectImage("C1-" + basename + f + "_Airy.czi");
	rename("HP1");
	selectImage("C2-" + basename + f + "_Airy.czi");
	rename("H3K9me3");
	selectImage("C3-" + basename + f + "_Airy.czi");
	rename("DNA");
	
	//Correct translation
	selectWindow("DNA");
	getDimensions(width, height,a,a,a);
	makeRectangle(0, 6, width, height-6);
	run("Crop");
	selectWindow("HP1");
	getDimensions(width, height,a,a,a);
	makeRectangle(0, 0, width, height-6);
	run("Crop");
	selectWindow("H3K9me3");
	getDimensions(width, height,a,a,a);
	makeRectangle(0, 0, width, height-6);
	run("Crop");
	
	imageCalculator("Add create 32-bit", "DNA", "HP1");
	rename("Sum");
	run("Duplicate...", " ");
	setAutoThreshold("Li dark no-reset");
	run("Convert to Mask");
	rename("Mask");
	run("Erode");
	run("Erode");
	run("Erode");
	run("Dilate");
	run("Dilate");
	run("Set Measurements...", "area mean min max area_fraction redirect=None decimal=3");
	run("Analyze Particles...", "size=200-Infinity pixel clear add");
	waitForUser;
	//folder = "/home/sasha/Science/PhD/Laue_lab/Microscopy/Confocal/Histone_marks/20220201_HP1b-JF646_Hoechst_H3K9me3-Cy3_CAICAiry/";
	//f=4;
	roiManager("Save", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_roi.zip");
	nROIs = roiManager("count");
	for (j=0; j<nROIs; j++){
		selectImage("DNA");
		run("Select None");
		run("Duplicate...", " ");
		roiManager("Select", j);
		setBackgroundColor(0, 0, 0);
		run("Clear Outside");
		setAutoThreshold("MaxEntropy dark"); //different from run("Auto Threshold") - different binning for 16-bit images
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Analyze Particles...", "size=40-Infinity pixel show=Masks clear");
		rename("DNA_mask" + j);
		run("Invert");
		run("Dilate");
		run("Erode");
//		run("Erode");
//		run("Dilate");
		waitForUser;
		if (j==0) {
			selectImage("DNA_mask"+j);
			rename("DNA_mask_old");
		} else {
			imageCalculator("Add create 32-bit", "DNA_mask"+j,"DNA_mask_old");
			rename("DNA_mask_new");
			close("DNA_mask_old");
			selectImage("DNA_mask_new");
			rename("DNA_mask_old");
		}
//			selectImage("HP1");
//			run("Duplicate...", " ");
//			roiManager("Select", j);
//			setBackgroundColor(0, 0, 0);
//			run("Clear Outside");
//			run("Auto Threshold", "method=MaxEntropy ignore_black white");
//			run("Analyze Particles...", "size=20-Infinity pixel show=Masks clear");
//			rename("HP1_mask" + j);
//			run("Invert");
//			waitForUser;
	}
	selectImage("DNA_mask_old");
	run("Convert to Mask");
	saveAs("tiff", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_DNA_mask.tif");
	rename("DNA_mask");
	
	selectImage("HP1");
	run("Duplicate...", " ");
	roiManager("Combine");
	setBackgroundColor(0, 0, 0);
	run("Clear Outside");
	setAutoThreshold("Triangle dark"); //different from run("Auto Threshold") - different binning for 16-bit images
	setOption("BlackBackground", true);
	run("Convert to Mask");
	//run("Auto Threshold", "method=Triangle ignore_black white");
	run("Analyze Particles...", "size=10-Infinity pixel show=Masks clear");
	rename("HP1_mask");
	run("Invert");
	waitForUser;
	selectImage("HP1_mask");
	run("Convert to Mask");
	saveAs("tiff", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_HP1_mask.tif");
	rename("HP1_mask");
	
	//Measurements per DNA focus
	run("Set Measurements...", "area mean min max redirect=H3K9me3 decimal=3");
	nROIs = roiManager("count");
	for (j=0; j<nROIs; j++){
		selectWindow("DNA_mask");
		roiManager("Select", j);
		run("Analyze Particles...", "size=20-Infinity pixel display summarize");
	}
	selectWindow("Summary");
	saveAs("Results", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_H3K9me3_DNAfocusSummary.csv");
	NumberofRows=Table.size("Pos" + d2s(f, 0) + "_H3K9me3_DNAfocusSummary.csv");
	Table.deleteRows(0, NumberofRows-1, "Pos" + d2s(f, 0) + "_H3K9me3_DNAfocusSummary.csv");
	selectWindow("Results");
	saveAs("Results", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_H3K9me3_perDNAFocus.csv");
	run("Clear Results");
	//Per DNA focus, adj by DNA intensity
	imageCalculator("Divide create 32-bit", "H3K9me3","DNA");
	rename("H3K9me3_norm");
	run("Set Measurements...", "area mean min max redirect=H3K9me3_norm decimal=3");
	nROIs = roiManager("count");
	for (j=0; j<nROIs; j++){
		selectWindow("DNA_mask");
		roiManager("Select", j);
		run("Analyze Particles...", "size=20-Infinity pixel display summarize");
	}
	selectWindow("Summary");
	saveAs("Results", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_H3K9me3norm_DNAfocusSummary.csv");
	NumberofRows=Table.size("Pos" + d2s(f, 0) + "_H3K9me3norm_DNAfocusSummary.csv");
	Table.deleteRows(0, NumberofRows-1, "Pos" + d2s(f, 0) + "_H3K9me3norm_DNAfocusSummary.csv");
	selectWindow("Results");
	saveAs("Results", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_H3K9me3norm_perDNAFocus.csv");
	run("Clear Results");
	
	//Measurements per HP1 focus
	run("Set Measurements...", "area mean min max redirect=H3K9me3 decimal=3");
	nROIs = roiManager("count");
	for (j=0; j<nROIs; j++){
		selectWindow("HP1_mask");
		roiManager("Select", j);
		run("Analyze Particles...", "size=10-Infinity pixel display summarize");
	}
	selectWindow("Summary");
	saveAs("Results", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_H3K9me3_HP1focusSummary.csv");
	NumberofRows=Table.size("Pos" + d2s(f, 0) + "_H3K9me3_HP1focusSummary.csv");
	Table.deleteRows(0, NumberofRows-1, "Pos" + d2s(f, 0) + "_H3K9me3_HP1focusSummary.csv");
	selectWindow("Results");
	saveAs("Results", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_H3K9me3_perHP1Focus.csv");
	run("Clear Results");
	//Per HP1 focus, adj by DNA intensity
	imageCalculator("Divide create 32-bit", "H3K9me3","DNA");
	rename("H3K9me3_norm");
	run("Set Measurements...", "area mean min max redirect=H3K9me3_norm decimal=3");
	nROIs = roiManager("count");
	for (j=0; j<nROIs; j++){
		selectWindow("HP1_mask");
		roiManager("Select", j);
		run("Analyze Particles...", "size=10-Infinity pixel display summarize");
	}
	selectWindow("Summary");
	saveAs("Results", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_H3K9me3norm_HP1focusSummary.csv");
	NumberofRows=Table.size("Pos" + d2s(f, 0) + "_H3K9me3norm_HP1focusSummary.csv");
	Table.deleteRows(0, NumberofRows-1, "Pos" + d2s(f, 0) + "_H3K9me3norm_HP1focusSummary.csv");
	selectWindow("Results");
	saveAs("Results", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_H3K9me3norm_perHP1Focus.csv");
	run("Clear Results");
	
	waitForUser;
	
	//Measurements per nucleus outside of foci
	imageCalculator("Add create", "DNA_mask", "HP1_mask");
	run("Convert to Mask");
	rename("Sum_mask");
	run("Set Measurements...", "area mean min max redirect=H3K9me3 decimal=3");
	nROIs = roiManager("count");
	for (j=0; j<nROIs; j++){
		selectWindow("Sum_mask");
		run("Select None");
		run("Duplicate...", " ");
		roiManager("Select", j);
		run("Invert");
		run("Clear Outside");
		run("Create Selection");
		run("Measure");
		close();
	}
	selectWindow("Results");
	saveAs("Results", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_H3K9me3_outsideFoci.csv");
	run("Clear Results");
	//Adj by DNA intensity
	imageCalculator("Divide create 32-bit", "H3K9me3","DNA");
	rename("H3K9me3_norm");
	run("Set Measurements...", "area mean min max redirect=H3K9me3_norm decimal=3");
	nROIs = roiManager("count");
	for (j=0; j<nROIs; j++){
		selectWindow("Sum_mask");
		run("Select None");
		run("Duplicate...", " ");
		roiManager("Select", j);
		run("Invert");
		run("Clear Outside");
		run("Create Selection");
		run("Measure");
		close();
	}
	selectWindow("Results");
	saveAs("Results", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_H3K9me3norm_outsideFoci.csv");
	run("Clear Results");
	
	//Correlations per cell
	nROIs = roiManager("count");
	R_DNA_HP1 = "";
	R_DNA_H3K9me3 = "";
	R_HP1_H3K9me3 = "";
	for (j=0; j<nROIs; j++){
		selectWindow("DNA");
		roiManager("Select", j);
		run("Coloc 2", "channel_1=DNA" + " channel_2=HP1" +
		" roi_or_mask=[ROI(s) in channel 1] threshold_regression=Bisection psf=3 costes_randomisations=10");
		string = getInfo("log");
		R_DNA_HP1 = R_DNA_HP1 + d2s(j, 0) + " " + substring(string, indexOf(string, "Pearson's R value (no threshold), ")+lengthOf("Pearson's R value (no threshold), "),
		indexOf(string, "Pearson's R value (below threshold)"));
		print("\\Clear");
		
		selectWindow("DNA");
		roiManager("Select", j);
		run("Coloc 2", "channel_1=DNA" + " channel_2=H3K9me3" +
		" roi_or_mask=[ROI(s) in channel 1] threshold_regression=Bisection psf=3 costes_randomisations=10");
		string = getInfo("log");
		R_DNA_H3K9me3 = R_DNA_H3K9me3 + d2s(j, 0) + " " + substring(string, indexOf(string, "Pearson's R value (no threshold), ")+lengthOf("Pearson's R value (no threshold), "),
		indexOf(string, "Pearson's R value (below threshold)"));
		print("\\Clear");
		
		selectWindow("HP1");
		roiManager("Select", j);
		run("Coloc 2", "channel_1=HP1" + " channel_2=H3K9me3" +
		" roi_or_mask=[ROI(s) in channel 1] threshold_regression=Bisection psf=3 costes_randomisations=10");
		string = getInfo("log");
		R_HP1_H3K9me3 = R_HP1_H3K9me3 + d2s(j, 0) + " " + substring(string, indexOf(string, "Pearson's R value (no threshold), ")+lengthOf("Pearson's R value (no threshold), "),
		indexOf(string, "Pearson's R value (below threshold)"));
		print("\\Clear");
	}
	File.saveString(R_DNA_HP1, folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_DNA_HP1_coloc.txt");
	File.saveString(R_DNA_H3K9me3, folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_DNA_H3K9me3_coloc.txt");
	File.saveString(R_HP1_H3K9me3, folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_HP1_H3K9me3_coloc.txt");
	
	run("Close All");
	run("Clear Results");
	roiManager("Delete");
}



////////////////////////////////////////
//Correction: due to outliers in Rep2, better to calc avg DNA and then divide avg H3K9me3 by avg DNA,
//rather than divide pixel by pixel.
//Thus, need to remeasure avg DNA intensities of foci and outfoci

folder = "/home/sasha/Science/PhD/Laue_lab/Microscopy/Confocal/Histone_marks/20220201_HP1b-JF646_Hoechst_H3K9me3-Cy3_CAICAiry/";
basename = "20220201_Pos";
for (f = 1; f < 4; f++) {
	roiManager("Open", folder + "Automatic_quantification/Pos" + f + "_roi.zip");
	
	//f=1;
	open(folder + basename + f + "_Airy.czi");
	run("Split Channels");
	selectImage("C1-" + basename + f + "_Airy.czi");
	rename("HP1");
	selectImage("C2-" + basename + f + "_Airy.czi");
	rename("H3K9me3");
	selectImage("C3-" + basename + f + "_Airy.czi");
	rename("DNA");
	
	//Correct translation
	selectWindow("DNA");
	getDimensions(width, height,a,a,a);
	makeRectangle(0, 6, width, height-6);
	run("Crop");
	selectWindow("HP1");
	getDimensions(width, height,a,a,a);
	makeRectangle(0, 0, width, height-6);
	run("Crop");
	selectWindow("H3K9me3");
	getDimensions(width, height,a,a,a);
	makeRectangle(0, 0, width, height-6);
	run("Crop");
	
	open(folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_DNA_mask.tif");
	rename("DNA_mask");
	open(folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_HP1_mask.tif");
	rename("HP1_mask");
	
	//Measurements per DNA focus
	run("Set Measurements...", "area mean min max redirect=DNA decimal=3");
	nROIs = roiManager("count");
	for (j=0; j<nROIs; j++){
		selectWindow("DNA_mask");
		roiManager("Select", j);
		run("Analyze Particles...", "size=20-Infinity pixel display");
	}
	selectWindow("Results");
	saveAs("Results", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_DNA_perDNAFocus.csv");
	run("Clear Results");
	
	//Measurements per HP1 focus
	run("Set Measurements...", "area mean min max redirect=DNA decimal=3");
	nROIs = roiManager("count");
	for (j=0; j<nROIs; j++){
		selectWindow("HP1_mask");
		roiManager("Select", j);
		run("Analyze Particles...", "size=10-Infinity pixel display");
	}
	selectWindow("Results");
	saveAs("Results", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_DNA_perHP1Focus.csv");
	run("Clear Results");
	
	//Measurements per nucleus outside of foci
	imageCalculator("Add create", "DNA_mask", "HP1_mask");
	run("Convert to Mask");
	rename("Sum_mask");
	run("Set Measurements...", "area mean min max redirect=DNA decimal=3");
	nROIs = roiManager("count");
	for (j=0; j<nROIs; j++){
		selectWindow("Sum_mask");
		run("Select None");
		run("Duplicate...", " ");
		roiManager("Select", j);
		run("Invert");
		run("Clear Outside");
		run("Create Selection");
		run("Measure");
		close();
	}
	selectWindow("Results");
	saveAs("Results", folder + "Automatic_quantification/Pos" + d2s(f, 0) + "_DNA_outsideFoci.csv");
	run("Clear Results");
	
	run("Close All");
	run("Clear Results");
	roiManager("Deselect");
	roiManager("Delete");
}
