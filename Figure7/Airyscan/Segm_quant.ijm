folder = "/home/sasha/Science/PhD/Laue_lab/Microscopy/Confocal/20241101_differentiation_HP1b_DNA_Airyscan_Devina/";
hl = newArray("2iL", "24h", "48h");
for (i = 0; i < lengthOf(hl); i++) {
	for (f = 5; f < 6; f++) {
		open(folder + "MAX_" + hl[i] + "_" + d2s(f, 0) + ".tif");
		run("Split Channels");
		selectImage("C2-MAX_" + hl[i] + "_" + d2s(f, 0) + ".tif");
		rename("DNA");
		selectImage("C1-MAX_" + hl[i] + "_" + d2s(f, 0) + ".tif");
		rename("HP1");
		imageCalculator("Add create 32-bit", "DNA","HP1");
		rename("Sum");
		run("Duplicate...", " ");
		setAutoThreshold("Li dark no-reset");
		run("Convert to Mask");
		rename("Mask");
		run("Set Measurements...", "area mean min max area_fraction redirect=None decimal=3");
		run("Analyze Particles...", "clear add");
		waitForUser;
		roiManager("Save", folder + "Automatic_quantification/" + hl[i] + "_Pos" + d2s(f, 0) + "_roi.zip");
		nROIs = roiManager("count");
		for (j=0; j<nROIs; j++){
			selectImage("DNA");
			run("Select None");
			run("Duplicate...", " ");
			roiManager("Select", j);
			setBackgroundColor(0, 0, 0);
			run("Clear Outside");
			run("Auto Threshold", "method=MaxEntropy ignore_black white");
			run("Analyze Particles...", "size=20-Infinity pixel show=Masks clear");
			rename("DNA_mask" + j);
			run("Invert");
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
		saveAs("tiff", folder + "Automatic_quantification/" + hl[i] + "_Pos" + d2s(f, 0) + "_DNA_mask.tif");
		rename("DNA_mask");
		
		selectImage("HP1");
		run("Duplicate...", " ");
		roiManager("Combine");
		setBackgroundColor(0, 0, 0);
		run("Clear Outside");
		run("Auto Threshold", "method=Triangle ignore_black white");
		run("Analyze Particles...", "size=10-Infinity pixel show=Masks clear");
		rename("HP1_mask");
		run("Invert");
		waitForUser;
		selectImage("HP1_mask");
		run("Convert to Mask");
		saveAs("tiff", folder + "Automatic_quantification/" + hl[i] + "_Pos" + d2s(f, 0) + "_HP1_mask.tif");
		rename("HP1_mask");
		
		//Measurements per cell
		nROIs = roiManager("count");
		run("Set Measurements...", "area mean min max area_fraction redirect=None decimal=3");
		for (j=0; j<nROIs; j++){
			selectWindow("DNA_mask");
			roiManager("Select", j);
			run("Measure");
			selectWindow("HP1_mask");
			roiManager("Select", j);
			run("Measure");
			selectWindow("DNA");
			roiManager("Select", j);
			run("Measure");
			selectWindow("HP1");
			roiManager("Select", j);
			run("Measure");
		}
		saveAs("Results", folder + "Automatic_quantification/" + hl[i] + "_Pos" + d2s(f, 0) + "_perCell.csv");
		run("Clear Results");
		
		//Measurements per DNA focus
		run("Set Measurements...", "area mean min max redirect=DNA decimal=3");
		nROIs = roiManager("count");
		for (j=0; j<nROIs; j++){
			selectWindow("DNA_mask");
			roiManager("Select", j);
			run("Analyze Particles...", "size=20-Infinity pixel display summarize");
		}
		selectWindow("Summary");
		saveAs("Results", folder + "Automatic_quantification/" + hl[i] + "_Pos" + d2s(f, 0) + "_focusSummary_DNA.csv");
		NumberofRows=Table.size(hl[i] + "_Pos" + d2s(f, 0) + "_focusSummary_DNA.csv");
		Table.deleteRows(0, NumberofRows-1, hl[i] + "_Pos" + d2s(f, 0) + "_focusSummary_DNA.csv");
		selectWindow("Results");
		saveAs("Results", folder + "Automatic_quantification/" + hl[i] + "_Pos" + d2s(f, 0) + "_perFocus_DNA.csv");
		run("Clear Results");
		
		//Measurements per HP1 focus
		run("Set Measurements...", "area mean min max redirect=HP1 decimal=3");
		nROIs = roiManager("count");
		for (j=0; j<nROIs; j++){
			selectWindow("HP1_mask");
			roiManager("Select", j);
			run("Analyze Particles...", "size=10-Infinity pixel display summarize");
		}
		selectWindow("Summary");
		saveAs("Results", folder + "Automatic_quantification/" + hl[i] + "_Pos" + d2s(f, 0) + "_focusSummary_HP1.csv");
		NumberofRows=Table.size(hl[i] + "_Pos" + d2s(f, 0) + "_focusSummary_HP1.csv");
		Table.deleteRows(0, NumberofRows-1, hl[i] + "_Pos" + d2s(f, 0) + "_focusSummary_HP1.csv");
		selectWindow("Results");
		saveAs("Results", folder + "Automatic_quantification/" + hl[i] + "_Pos" + d2s(f, 0) + "_perFocus_HP1.csv");
		run("Clear Results");
		
		waitForUser; //why???
		
		//Measurements per nucleus outside of DNA foci
		run("Set Measurements...", "area mean min max redirect=DNA decimal=3");
		nROIs = roiManager("count");
		for (j=0; j<nROIs; j++){
			selectWindow("DNA_mask");
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
		saveAs("Results", folder + "Automatic_quantification/" + hl[i] + "_Pos" + d2s(f, 0) + "_outsideFoci_DNA.csv");
		run("Clear Results");
		
		//Measurements per nucleus outside of HP1 foci
		run("Set Measurements...", "area mean min max redirect=HP1 decimal=3");
		nROIs = roiManager("count");
		for (j=0; j<nROIs; j++){
			selectWindow("HP1_mask");
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
		saveAs("Results", folder + "Automatic_quantification/" + hl[i] + "_Pos" + d2s(f, 0) + "_outsideFoci_HP1.csv");
		run("Clear Results");
		
		run("Close All");
		run("Clear Results");
		roiManager("Delete");
	}
}