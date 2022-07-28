// Macro to calculate contractility parameters from image stack
//
// Author: Francis Burton (Francis.Burton (at) glasgow.ac.uk)
// Created 14/5/15 FB  -  Original version
// Modified 28/5/15 FB -  Made showing of individual transients a separate option.
//                        For generality, now refer to "transients" rather than "contractions".
//                        Added initial image scaling step (optional).
//                        Calling macro with argument "noask" skips interactive menu and uses
//                        previously entered values (stored in the preferences file, and can be
//                        set using 'call("ij.Prefs.set",...').
//                        To make above work, all settings (including checkboxes) are now stored
//                        between runs.
// Modified 29/5/15 FB -  Added option to save final profile data to text file automatically.
//                        Modified manual mode to allow user to mark transient start points on
//                        profile plot rather than having to type them in.
//                        Added version date to initial dialog.
// Modified 31/5/15 FB -  Changed default algorithm to find start of transient: "Start offset"
//                        now specifies number of frames before shoulder/bend point to search for
//                        a lower point, rather than simply offsetting the start by the specified
//                        number of frames. However, the latter (old) behaviour can be achieved
//                        by specifying a negative start offset.
//                        Added mean, SD, min & max statistics to results.
//                        Now report all transient times (TU10 to TD90) relative to T0.
//                        Added optional reporting of times in ms rather than s.
//						  Added option to not smooth transient before calculating times.
//                        Add option to save results table to csv file.
//                        (Profile and results data are saved to same folder as original images.)
//                        Rearranged dialog items more logically.
// Modified 1/6/15 FB -   Added option to calculate average trace and corresponding statistics,
//                        and save average trace to file.
// Modified 2/6/15 FB -   Added various checks of transient quality.
// Modified 31/7/15 FB -  Changed settings file parameter names to include prepended "c1_" to
//                        distinguish them from parameters in other analyses and prevent crosstalk.
// Modified 2/11/15 FB -  Added option to search forward from peak for downstroke levels (invert
//                        downstroke search).
//                        Added option to filter/smooth average trace before analysis and export.
// Modified 14/4/22 FB -  Added averaging after aligning with blip markers in recording
// Modified 15/6/22 FB -  Added blip detection level adjustment
// Modified 16/6/22 FB -  Added blip frame normalization option & retired manual mode
//

// Installation: Copy this file into the folder called "macros" that is inside the folder where ImageJ
// or Fiji is located, possibly replacing a previous version.
//
// To add a button labelled C1 to the ImageJ toolbar that invokes this macro, append the following 3 
// lines starting with '//macro' to the file called StartupMacros (StartupMacros.fiji.ijm for Fiji, or
// StartupMacros.txt for plain ImageJ), deleting the initial '//'.
//macro "Contractility Action Tool - C800C080C880C008D11D12D13D14D15D16D17D18D19D1aD1bD1cD1dD1eD21D22D23D2cD2dD2eD31D32D34D35D36D37D38D39D3aD3bD3dD3eD41D43D44D45D46D47D48D49D4aD4bD4cD4eD51D53D54D55D56D57D58D59D5aD5bD5cD5eD61D62D64D65D66D67D68D69D6aD6bD6dD6eD71D72D73D74D75D76D77D78D79D7aD7bD7cD7dD7eD81D82D83D84D85D86D87D88D89D8aD8bD8cD8dD8eD91D92D94D95D96D97D98D99D9aD9bD9cD9eDa1Da3Da4Da5Da6Da7Da8Da9DaaDabDacDaeDb1DbeDc1Dc2Dc3Dc4Dc5Dc6Dc7Dc8Dc9DcaDcbDccDceDd1Dd2Dd3Dd4Dd5Dd6Dd7Dd8Dd9DdaDdbDdcDdeDe1De2De3De4De5De6De7De8De9DeaDebDecDedDeeC808D00D01D02D03D04D05D06D07D08D09D0aD0bD0cD0dD0eD0fD10D1fD20D2fD30D3fD40D4fD50D5fD60D6fD70D7fD80D8fD90D9fDa0DafDb0DbfDc0DcfDd0DdfDe0DefDf0Df1Df2Df3Df4Df5Df6Df7Df8Df9DfaDfbDfcDfdDfeDffC088CcccC888Cf00C0f0Cff0C00fCf0fC0ffD24D25D26D27D28D29D2aD2bD33D3cD42D4dD52D5dD63D6cD93D9dDa2DadDb2Db3Db4Db5Db6Db7Db8Db9DbaDbbDbcDbdDcdDdd"{
//     runMacro("Contractility");
//}

// Global variables
var	VERSION_DATE = "16/6/22";
var SHOW_EXTENDED_PARAMETERS = false;
var BATCH_MODE = false;
var SAVE_PLOTS = false;
var FILE_ROOT = "";
var BLIP_MODE = 5;
var MIN_BLIP_AMPLITUDE = 100;
var MIN_BLIP_INTERVAL = 480;
var MAX_BLIP_INTERVAL = 3020;
var BLIP_TOLERANCE_PC = 5;
var BLIP_NORMALIZATION_PC = 99;

macro "Contractility" {
	
	requires("1.49t");

	// Read saved parameter values from settings file (IJ_Prefs.txt), using sensible
	// default values if absent from settings file.
	default_imagescaling = readSettingValue("c1.imagescaling", 100);
	default_frameinterval = readSettingValue("c1.frameinterval", 0.01);
	default_windowlength = readSettingValue("c1.windowlength", default_frameinterval*10);
	default_nwin = readSettingValue("c1.nwin", 0);
	default_startoffset = readSettingValue("c1.startoffset", 2);
	default_qualitythreshold = readSettingValue("c1.qualitythreshold", 50);
	default_cvthreshold = readSettingValue("c1.cvthreshold", 5);
	default_blipthreshold = readSettingValue("c1.blipthreshold", 50);
	default_averagemode = readSettingValue("c1.averagemode", 0);
	default_normalizeblipframes = readSettingValue("c1.normalizeblipframes", 0);
	default_filterdata = readSettingBoolean("c1.filterdata", false);
	default_showintermediate = readSettingBoolean("c1.showintermediate", false);
	default_showtransients = readSettingBoolean("c1.showtransients", false);
	default_extendedparameters = readSettingBoolean("c1.extendedparameters", false);
	default_smoothtransient = readSettingBoolean("c1.smoothtransient", true);
	default_smoothaverage = readSettingBoolean("c1.smoothaverage", false);
	default_invertdownstrokesearch = readSettingBoolean("c1.invertdownstrokesearch", false); 
	default_reportmilliseconds = readSettingBoolean("c1.reportmilliseconds", false);
	default_saveprofile = readSettingBoolean("c1.saveprofile", false);
	default_saveresults = readSettingBoolean("c1.saveresults", false);
	default_saveplots = readSettingBoolean("c1.saveplots", false);

    print("nImages  " + nImages);
	// Import image sequence if nothing is open already
	if (nImages==0){
		
		run("Image Sequence...");
	}

	// Get stack dimensions
	getDimensions(width, height, channels, slices, frames);
	
	// Check for suitable stack
	if (slices<2)
		exit("This macro requires a stack");
	if (slices<20)
		exit("This macro requires a stack of at least 20 slices");
	if (channels!=1)
		exit("This macro doesn't work with multiple-channel images");
		
	// Remember original stack name for creating profile file later
	dirname = getDirectory("image");
	dataname = getTitle();
	FILE_ROOT = dirname + dataname.replace(".tif", "");

	 print("FILE_ROOT  " + FILE_ROOT);

		
	if (getArgument()=="noask") {	// For batch mode operation
		BATCH_MODE = true;
		image_scaling = default_imagescaling;
		frame_interval = default_frameinterval;
		wind_length = default_windowlength;
		nwin = default_nwin;
		start_offset = default_startoffset;
		quality_threshold = default_qualitythreshold;
		cv_threshold = default_cvthreshold;
		blip_threshold = default_blipthreshold;
		average_mode = default_averagemode;
		normalize_blip_frames = default_normalizeblipframes;
		detection_filter = default_filterdata;
		show_intermediate_results = default_showintermediate;
		show_transients = default_showtransients;
		SHOW_EXTENDED_PARAMETERS = default_extendedparameters;
		smooth_transient = default_smoothtransient;
		smooth_average = default_smoothaverage;
		invert_downstroke_search = default_invertdownstrokesearch;
		report_milliseconds = default_reportmilliseconds;
		save_profile = default_saveprofile;
		save_results = default_saveresults;
		SAVE_PLOTS = default_saveplots;
	} else {		
		// Get analysis parameters via dialog box, using settings file (or default) values
		showStatus("Get analysis parameters");
		Dialog.create("Analysis parameters ("+VERSION_DATE+")");
		Dialog.addNumber("Image scaling", default_imagescaling, 0, 4, "%");
		Dialog.addNumber("Frame interval", default_frameinterval, 2, 6, "s");
		Dialog.addNumber("Window duration", default_windowlength, 2, 6, "s");
		Dialog.addNumber("(or window length", default_nwin, 0, 4, "frames)");
		Dialog.addNumber("Start offset", default_startoffset, 0, 4, "frames");
		Dialog.addNumber("Quality threshold", default_qualitythreshold, 0, 4, "%");
		Dialog.addNumber("CV threshold", default_cvthreshold, 0, 4, "%");
		Dialog.addNumber("Blip threshold", default_blipthreshold, 0, 4, "%");
		average_items = newArray("None", "Align starts", "Align 10% ups", "Align 50% ups", "Align peaks", "Align blips");
		Assert((average_items[BLIP_MODE] == "Align blips"), "BLIP_MODE has incorrect value");
		Dialog.addChoice("Calculate average", average_items, average_items[default_averagemode]); 
		Dialog.addCheckbox("Normalize blip frames", default_normalizeblipframes);
		Dialog.addCheckbox("Detection filter (optional)", default_filterdata);
		Dialog.addCheckbox("Show intermediate results", default_showintermediate);
		Dialog.addCheckbox("Show individual transients", default_showtransients);
		Dialog.addCheckbox("Show extended parameters", default_extendedparameters);
		Dialog.addCheckbox("Analyze smoothed transients", default_smoothtransient);
		Dialog.addCheckbox("Smooth average", default_smoothaverage);
		Dialog.addCheckbox("Invert downstroke search", default_invertdownstrokesearch);
		Dialog.addCheckbox("Report parameters as ms", default_reportmilliseconds);
		Dialog.addCheckbox("Save profile data (x y) as txt", default_saveprofile);
		Dialog.addCheckbox("Save results table as csv", default_saveresults);
		Dialog.addCheckbox("Save plots as png", default_saveplots);
		Dialog.show();
		
		image_scaling = Dialog.getNumber();
		frame_interval = Dialog.getNumber();
		wind_length = Dialog.getNumber();
		nwin = Dialog.getNumber();
		start_offset = Dialog.getNumber();
		quality_threshold = Dialog.getNumber();
		cv_threshold = Dialog.getNumber();
		blip_threshold = Dialog.getNumber();
		average_mode_string = Dialog.getChoice();
		for (i=0; i<average_items.length; i++)
			if (average_items[i]==average_mode_string)
				average_mode = i;
		normalize_blip_frames = Dialog.getCheckbox();
		detection_filter = Dialog.getCheckbox();
		show_intermediate_results = Dialog.getCheckbox();
		show_transients = Dialog.getCheckbox();
		SHOW_EXTENDED_PARAMETERS = Dialog.getCheckbox();
		smooth_transient = Dialog.getCheckbox();
		smooth_average = Dialog.getCheckbox();
		invert_downstroke_search = Dialog.getCheckbox();
		report_milliseconds = Dialog.getCheckbox();
		save_profile = Dialog.getCheckbox();
		save_results = Dialog.getCheckbox();
		SAVE_PLOTS = Dialog.getCheckbox();

		// Save new parameter values to settings file for next time
		writeSettingValue("c1.imagescaling", image_scaling, 1);
		writeSettingValue("c1.frameinterval", frame_interval, 2);
		writeSettingValue("c1.windowlength", wind_length, 2);
		writeSettingValue("c1.nwin", nwin, 0);
		writeSettingValue("c1.startoffset", start_offset, 0);
		writeSettingValue("c1.qualitythreshold", quality_threshold, 1);
		writeSettingValue("c1.cvthreshold", cv_threshold, 1);
		writeSettingValue("c1.blipthreshold", blip_threshold, 1);
		writeSettingValue("c1.averagemode", average_mode, 0);
		writeSettingBoolean("c1.normalizeblipframes", normalize_blip_frames);
		writeSettingBoolean("c1.filterdata", detection_filter);
		writeSettingBoolean("c1.showintermediate", show_intermediate_results);
		writeSettingBoolean("c1.showtransients", show_transients);
		writeSettingBoolean("c1.extendedparameters", SHOW_EXTENDED_PARAMETERS);
		writeSettingBoolean("c1.smoothtransient", smooth_transient );
		writeSettingBoolean("c1.smoothaverage", smooth_average);
		writeSettingBoolean("c1.invertdownstrokesearch", invert_downstroke_search);
		writeSettingBoolean("c1.reportmilliseconds", report_milliseconds);
		writeSettingBoolean("c1.saveprofile", save_profile);
		writeSettingBoolean("c1.saveresults", save_results);
		writeSettingBoolean("c1.saveplots", SAVE_PLOTS);
	}
	
	// Check parameters have sensible values
	if (frame_interval<=0)
		exit("Invalid frame interval: must be positive");
	if (nwin == 0) {
		if (wind_length<=0)
			exit("Invalid window length: must be positive");
		nwin = round(wind_length/frame_interval);
		}
	if (nwin < 1) {
		nwin = 1;
		print("Detection window length set to 1 frame");
		}
	nframes = nSlices;
	if (nwin > nframes/2) {
		nwin = round(nframes/2);
		print("Detection window length set to ", nwin);
		}

	// Optionally resize stack images
	if (image_scaling!=100) {
		factor = image_scaling/100;
		sf = d2s(image_scaling/100, 4);
		s = "x=" + sf + " y=" + sf + " z=1.0";
		s = s + " width=" + i2s(round(width*factor));
		s = s + " height=" + i2s(round(height*factor));
		s = s + " depth=" + i2s(slices);
		s = s + " interpolation=Bilinear average process create title=Scaled";
		run("Scale...", s);
		closeWindow(dataname);
		}
		
	// Rename stack
	rename("Seq");
	
	setBatchMode(true);

	// The "blipper" is a pulsed light source that allows the time of an electrical
	// stimulus to be stored in the recording as single lighter or darker frames.
	use_blipper = (average_mode == BLIP_MODE);
	
	if (use_blipper) {
		// Look for blips (short excursions/pulses in mean image amplitude)
		blip_start = getBlips(blip_threshold);
		nblips = blip_start.length;
		if (nblips==0)
			exit("No blips detected");
			
		checkBlipGoodness(blip_start, frame_interval);
		
		// Experimental feature to remove blip artefact from the lighter or darker frame
		if (normalize_blip_frames) {
			// Remove blips from image data
			if (show_intermediate_results) {
				fs = "frame";
				if (nblips>1)
					fs = fs+"s";
				print("Normalizing " + i2s(nblips) + " blip " + fs);
				}

			selectWindow("Seq");
			run("Select All");

			for (i=0; i<nblips; i++) {
				kblip = blip_start[i]+1;
				kref = kblip-1;
				if (kblip==0) {
					kref = 2;
					kblip = 1;
					}
				setSlice(kref);
				href = getPercentileLevel(BLIP_NORMALIZATION_PC);
				setSlice(kblip);
				hblip = getPercentileLevel(BLIP_NORMALIZATION_PC);
				scale_level = href/hblip;
				sl = d2s(scale_level,4);
				
				if (show_intermediate_results) {
					print("Blip " + i2s(i+1) + "/" + i2s(nblips));
					print("Blip frame: " + i2s(kblip) + " (Ref frame: " + i2s(kref) + ")");
					print("Blip level: " + d2s(hblip,2) + ", Ref level: " + d2s(href,2)	+ ", Scale: " + sl);
					}
					
				s = "value=" + sl + " slice";
				run("Multiply...", s);
				}
							
				if (show_intermediate_results)
					print(" ");
			}
		}
	
	// Choose a suitable diastolic still frame (ref_frame) for global difference calculation:
	// 1) calculate running window differences
	showStatus("Calculating running window differences");
	getStackDifferenceProfile("Seq", "A", "B", 1, nframes-nwin, 1+nwin, nframes);
	
	// 2) find lowest point in running window differences
	Plot.getValues(x, y);
	yf = smoothData(y);
	minima = Array.findMinima(yf, 1);
	ref_frame = minima[0] + floor(nwin/2);

	if (show_intermediate_results && !use_blipper) {
		xmark = newArray(1);
		xmark[0] = ref_frame;
		plotTrace("Running diff", x, yf, xmark, "red");
		savePlot("Running diff");
		setBatchMode("show");
		
		selectWindow("Result of A");
		rename("Running diff");
		}
	else
		closeWindow("Result of A");

	closeWindow("Result of A-0-0");
	closeWindow("A");
	closeWindow("B");
	setBatchMode("exit and display");
	setBatchMode(true);
	
	// Calculate global difference trace
	showStatus("Calculating global difference trace");
	getImageDifferenceProfile("Seq", "Ref", ref_frame);
	
	Plot.getValues(x, y);

	n = y.length;
	// Adjust zero point at reference frame (where Ref-Ref=0)
	for (i=0; i<n; i++) {
		if (y[i]==0) {
			if (i==0)	// first point: use point after
				y[i] = y[i+1];
			else
				y[i] = y[i-1];
			}
		}
		
	if (use_blipper) {
		if (!normalize_blip_frames) {
			// Remove blips (large amplitude excursions) from trace data
			npoints = y.length;
			ysmooth = smoothData(y);
		
			for (i=0; i<nblips; i++) {
				k = blip_start[i];
				if (k>1)
					y[k-1] = y[k-2];
				else
					y[k-1] = ysmooth[k-1];
				
				if (k>0)
					y[k] = y[k-1];
				else
					y[k] = ysmooth[k];
				
				if (k<npoints-2)
					y[k+1] = (y[k] + y[k+2])/2;
				else
					y[k+1] = ysmooth[k+1];
				}
			}
		
		// Note that a negative start_offset will shift start later
		transient_start = Array.copy(blip_start);
		for (i=0; i<nblips; i++)
			transient_start[i] = blip_start[i] - start_offset;
		
		}
	else {	
		if (detection_filter)
			y = smoothData(y);
		
		Array.getStatistics(y, ymin, ymax, ymean, ysd);

		// Calculate 1st derivative
		dy = newArray(n);
		dy[0] = 0;
		for (i=0; i<n-1; i++)
			dy[i+1] = y[i+1]-y[i];
		dyf = smoothData(dy);   

		// Find 1st derivative peaks
		Array.getStatistics(dyf, dymin, dymax, dymean, dysd);
		eps = (dymax-dymin)/2;
		maxima_dyf = Array.findMaxima(dyf, eps);
		Array.sort(maxima_dyf);
		nevents = maxima_dyf.length;

		if (show_intermediate_results) {
			plotTrace("1st derivative smoothed", x, dyf, maxima_dyf, "#00aa00");
			savePlot("1stDS");
			selectWindow("Result of Seq");
			rename("Global diff");
			}
		else
			closeWindow("Result of Seq");
		closeWindow("Result of Seq-0-0");
		closeWindow("Ref");
		setBatchMode("exit and display");
		setBatchMode(true);

		if (nevents<2)
			exit("Need at least 2 peaks in 1st derivative");

		// Find trace minima between peaks of 1st derivative
		minima_y = newArray(nevents);
		for (i=0; i<nevents; i++) {
			if (i==0)
				j1 = 0;
			else
				j1 = maxima_dyf[i-1];
			j2 = maxima_dyf[i];
			ymin = y[j1];
			jmin = 0;
			for (j=j1; j<=j2; j++) {
				if (y[j]<=ymin) {
					ymin = y[j];
					jmin = j;
					}
				}
			if (jmin==0)
				jmin = -1;
			minima_y[i] = jmin;
			}

		// Find shoulder (bend) at start of upstroke as maximum perpendicular
		// distance from trace points to line joining where the trace is lowest
		// to where upstroke is steepest
		bend_y = newArray(nevents);
		for (i=0; i<nevents; i++) {
			j1 = minima_y[i];
			j2 = maxima_dyf[i];
			if (j1<=0 || j1==j2)
				bend_y[i] = -1;
			else {
				y1 = y[j1];
				y2 = y[j2];
				dmax = 0;
				jmax = -1;
				for (j=j1; j<=j2; j++) {
					d = perpDistance(j1,y1, j2,y2, j,y[j]);
					if (d>dmax) {
						dmax = d;
						jmax = j;
						}
					}
				bend_y[i] = jmax;
				}
			}

		transient_start = newArray(bend_y.length);
		
		// Shift transient start by specified offset, either back by a fixed amount
		// or back to the minimum of the trace within the specified number of samples
		for (i=0; i<bend_y.length; i++) {
			if (start_offset <= 0)
				transient_start[i] = bend_y[i] + start_offset;
			else {
				j = bend_y[i];
				j0 = j - start_offset;
				if (j0 >= 0) {
					ymin = y[j];
					jmin = j;
					while (j >= j0) {
						if (y[j] < ymin) {
							ymin = y[j];
							jmin = j;
							}
						j--;
						}
					transient_start[i] = jmin;
					}
				else
					transient_start[i] = -1;
				}
			}

		}
		
	// Cull invalid negative start points and points near ends of trace
	temp = newArray(0);
	for (i=0; i<transient_start.length; i++) {
		if (transient_start[i]>=0)
			temp = Array.concat(temp, transient_start[i]);
		}
	transient_start = Array.copy(temp);

	nevents = transient_start.length;
	if (nevents<2)
		exit("Need at least two transients");
		
	if (show_intermediate_results) {
		if (use_blipper)
			plotTrace("Signal", x, y, transient_start, "blue"); 
		else
			plotTraceWithCrosses("Signal", x, y, transient_start, "#aa0000", 
				bend_y, "red", minima_y, "blue", maxima_dyf, "#00aa00");
		savePlot("Signal");
		setBatchMode("show");
		}
			
	showStatus("Baseline correcting individual transients");
	
	// Remove baseline by subtracting a piecewise linear line between successive
	// transient starts
	k = 0;
	yc = newArray(0);
	start_y = newArray(0);
	for (i=0; i<nevents-1; i++) {
		j1 = transient_start[i];
		j2 = transient_start[i+1];
		if (use_blipper) {
			j1 = j1+1;
			j2 = j2+1;
			}
		yt = Array.slice(y, j1, j2);
		y2 = yt[yt.length-1];
		y1 = yt[0];
		for (j=0; j<yt.length; j++)
			yt[j] = yt[j] - ((((y2-y1)/(j2-j1))*j)+y1);
		yc = Array.concat(yc, yt);
		start_y = Array.concat(start_y, k);
		k = k + (j2-j1);
		}

	// Append index of end of signal profile
	k = yc.length-1;
	start_y = Array.concat(start_y, k);

	xc = Array.copy(yc);
	t = 0;
	for (i=0; i<xc.length; i++) {
		xc[i] = t;
		t = t+frame_interval;
		}
	
	// Save corrected signal profile to text file if requested
	if (save_profile) {
		showStatus("Writing profile data to file");
		f = File.open(FILE_ROOT+".txt");
		for (i=0; i<xc.length; i++)
			print(f, xc[i]+"\t"+yc[i]);
		File.close(f);
		}

	plotTrace("Corrected signal", xc, yc, start_y, "#0000cc");
	savePlot("Corrected");
	setBatchMode("show");
	
	run("Clear Results");
	showStatus("Analysing transients");

	start_avg = newArray(0);
	tvals = newArray(0);
	
	// Analyze each transient in turn
	for (i=0; i<start_y.length-1; i++) {
		j1 = start_y[i];
		j2 = start_y[i+1];
		xx = Array.slice(xc, j1, j2);
		yy = Array.slice(yc, j1, j2);
		checkTransientGoodness(i+1, yy, quality_threshold);
		atimes = analyzeTransient(i+1, xx, yy, smooth_transient, report_milliseconds, show_transients);
		if (use_blipper)
			start_avg = Array.concat(start_avg, start_y[i]);
		else if (average_mode>0)
			start_avg = Array.concat(start_avg, round(atimes[average_mode-1]/frame_interval));
		tvals = Array.concat(tvals, xx[xx.length-1]-xx[0]);
		}
	
	updateResults();
	run("Summarize");		// Add mean, SD, min & max statistics
	
	// Check interval variation to highlight gross irregularities in signal
	Array.getStatistics(tvals, tmin, tmax, tmean, tsd);
	tcv = (tsd/tmean)*100;
	if (tcv>cv_threshold)
		print("Coefficient of variation of transient intervals", tcv, " exceeds threshold", cv_threshold, "%");
	
	// Optionally calculate average transient and statistics
	ntransients = start_avg.length; 
	if (average_mode>0 && ntransients>1) {
		showStatus("Calculating average");
		// start_y has one extra element
		Assert((start_y.length==start_avg.length+1), "Start array length mismatch");
		// Find smallest pre- and post-alignment times, the sum of which will
		// determine the duration of the averaged signal.
		// (Both start_y and start_avg contain indices not times.)
		minpre = xc.length;
		for (i=0; i<ntransients; i++) {
			ipre = start_avg[i] - start_y[i];
			Assert((ipre>=0), "Negative pre-time");
			if (ipre<minpre)
				minpre = ipre;
			}
		minpost = xc.length;
		for (i=0; i<ntransients; i++) {
			ipost = start_y[i+1] - start_avg[i];
			Assert((ipost>0), "Negative or zero post-time");
			if (ipost<minpost)
				minpost = ipost;
			}
		navg = minpre + minpost;
		if (navg<5)
			exit("Average has fewer than 5 points");
		xavg = newArray(navg);
		yavg = newArray(navg);

		// Cumulate sum of individual transients, aligned by start_avg
		for (i=0; i<navg; i++) {
			xavg[i] = (i-minpre)*frame_interval;
			yavg[i] = 0;
			}

		for (i=0; i<ntransients; i++) {
			j1 = start_avg[i] - minpre;
			j2 = start_avg[i] + minpost;
			k = 0;
			for (j=j1; j<j2; j++) {
				yavg[k] = yavg[k] + yc[j];
				k++;
				}
			}

		for (i=0; i<navg; i++)
			yavg[i] = yavg[i] / ntransients;
			
		if (smooth_average)
			yavg = smoothData(yavg);

		// Analyze and plot average transient
		// (The "false" argument inhibits smoothing; the "true" arguments forces plotting.)
		analyzeTransient(0, xavg, yavg, false, report_milliseconds, true);

		// Adjust labels for clarity
		lastrow = nResults-1;
		setResult("Label", lastrow, "Average");
		setResult("T0", lastrow, "");
		setResult("Interval", lastrow, "");
		
		// Save average profile to text file
		f = File.open(FILE_ROOT+"_avg"+".txt");
		for (i=0; i<navg; i++)
			print(f, xavg[i]+"\t"+yavg[i]);
		File.close(f);
		}
	
	if (save_results) {
		updateResults();
		saveAs("results", FILE_ROOT+".csv");
		}
		
	setBatchMode(false);
		
	showStatus("Done");
	
	if (isOpen("Log"))
		showMessage("See Log window for messages");

} // End macro "Contractility"


function getStackDifferenceProfile(sname, aname, bname, a1, a2, b1, b2) {
	selectWindow(sname);
	run("Make Substack...", "  slices=" + i2s(a1) + "-" + i2s(a2));
	rename(aname);
	selectWindow(sname);
	run("Make Substack...", "  slices=" + i2s(b1) + "-" + i2s(b2));
	rename(bname);
	imageCalculator("Difference create 32-bit stack", aname, bname);
	run("Select All");
	run("Plot Z-axis Profile");
	}
	
function getImageDifferenceProfile(sname, rname, iframe) {
	selectWindow(sname);
	run("Make Substack...", "  slices=" + i2s(iframe));
	rename(rname);
	imageCalculator("Difference create 32-bit stack", sname, rname);
	run("Select All");
	run("Plot Z-axis Profile");
	}
	
function i2s(i) {
	return d2s(i, 0);
	}
	
function plotTrace(title, x, y, xmarkers, markcolor) {
	Plot.create(title, "Time (s)", "y", x, y);
	Plot.setFormatFlags("1100000011");
	Plot.setFontSize(9);
	if (xmarkers.length>0)
		plotVertLines(x, y, xmarkers, markcolor);
	Plot.setFontSize(11);
	Plot.setColor("black");
	Plot.setFrameSize(800, 300);
	Plot.show();
	}
	
function plotTraceWithCrosses(title, x, y, xm0, c0, xm1, c1, xm2, c2, xm3, c3) {
	Plot.create(title, "Time (s)", "y", x, y);
	Plot.setFormatFlags("1100000011");
	Plot.setFontSize(9);
	plotVertLines(x, y, xm0, c0);
	Plot.setLineWidth(2);
	plotCrosses(x, y, xm1, c1);
	plotCrosses(x, y, xm2, c2);
	plotCrosses(x, y, xm3, c3);
	Plot.setFontSize(11);
	Plot.setColor("black");
	Plot.setLineWidth(1);
	Plot.setFrameSize(800, 300);
	Plot.show();
	}
	
function plotVertLines(x, y, xmarkers, markcolor) {
	Plot.setColor(markcolor);
	Array.getStatistics(y, ymin, ymax, ymean, ysd);
	ymargin = (ymax-ymin)/20;
	for (i=0; i<xmarkers.length; i++) {
		xm = xmarkers[i];
		if (xm>=0 && xm<x.length) {
			mx = x[xm];
			Plot.drawLine(mx, ymin-ymargin, mx, ymax+ymargin);
			}
		}
	}
	
function plotHorzLine(x1, x2, y, linecolor) {
	Plot.setColor(linecolor);
	Plot.drawLine(x1, y, x2, y);
	}
	
function plotCrosses(x, y, xm, c) {
	mx = newArray(0);
	my = newArray(0);
	for (i=0; i<xm.length; i++) {
		k = xm[i];
		if (k>=0 && k<x.length) {
			mx = Array.concat(mx, x[k]);
			my = Array.concat(my, y[k]);
			}
		}		
	Plot.setColor(c);
	Plot.add("crosses", mx, my);
	}
	
function checkTransientGoodness(nth, y, qthresh) {
	ms = "Transient " + i2s(nth);
	// Must have at least 3 points
	ny = y.length;
	if (ny<3)
		print(ms, "has fewer than 3 points");
	else {
		// Find peak
		ymax = y[0];
		imax = 0;
		for (i=0; i<ny; i++) {
			if (y[i]>ymax) {
				ymax = y[i];
				imax = i;
				}
			}
		if (imax==0)
			print(ms, "maximum is at start");
		else if (imax==ny-1)
			print(ms, "maximum is at end");
		else {
			y0 = y[0];
			yn = y[ny-1];
			dup = abs(y0-ymax);
			ddn = abs(yn-ymax);
			dmx = maxOf(dup, ddn);
			if (dup/dmx < 0.5)
				print(ms, "has abnormally high start");
			else if (ddn/dmx < 0.5)
				print(ms, "has abnormally high end");
			else {
				// Scan from start to peak
				sum_good = 0;
				sum_bad = 0;
				ydiff = 0;
				penalty = 1;
				for (i=0; i<=imax; i++) {
					if (i>0)
						ydiff = y[i] - y[i-1];
					if (ydiff >= 0) {
						sum_good = sum_good + ydiff;
						penalty = 1;
						}
					else {
						sum_bad = sum_bad - ydiff*penalty;
						penalty++;
						}
					}
				bgu = sum_good/(sum_good+sum_bad) * 100;		// 100%=perfect, 0%=rubbish
				if (bgu<qthresh)
					print(ms, "has poor quality upstroke, index", d2s(bgu,1));

				// Scan from peak to end
				sum_good = 0;
				sum_bad = 0;
				ydiff = 0;
				penalty = 1;
				for (i=imax; i<ny; i++) {
					if (i>imax)
						ydiff = y[i] - y[i-1];
					if (ydiff <= 0) {
						sum_good = sum_good - ydiff;
						penalty = 1;
						}
					else {
						sum_bad = sum_bad + ydiff*penalty;
						penalty++;
						}
					}
				bgd = sum_good/(sum_good+sum_bad) * 100;
				if (bgd<qthresh)
					print(ms, "has poor quality downstroke, index", d2s(bgd,1));
				}
			}
		}
	}

function analyzeTransient(nth, x, y, smot, rpms, plotflag) {
	// Takes (x,y) arrays for one transient profile,
	// calculates transient parameters
	// and adds parameters to a new line in Results table
	x0 = x[0];
	xlast = x[x.length-1];
	if (smot)
		yf = smoothData(y);
	else
		yf = Array.copy(y);
	Array.getStatistics(yf, ymin, y100, ymean, ystd);
	// Start of transient is zero by definition
	m = Array.findMaxima(yf, y100*0.9);
	if (m.length>1)
		print("Transient ", nth, " has more than one peak");
		
	// Following times are absolute; 
	// for time relative to start of transient, subtract x0.
	
	// Peak time
	i100 = m[0];
	x100 = x[i100];
	
	// Find upstroke times
	i0 = 0;
	u10 = findPercentLevelCrossing(x, y, 10, 1, i0);
	u20 = findPercentLevelCrossing(x, y, 20, 1, i0);
	u25 = findPercentLevelCrossing(x, y, 25, 1, i0);
	u30 = findPercentLevelCrossing(x, y, 30, 1, i0);
	u40 = findPercentLevelCrossing(x, y, 40, 1, i0);
	u50 = findPercentLevelCrossing(x, y, 50, 1, i0);
	u60 = findPercentLevelCrossing(x, y, 60, 1, i0);
	u70 = findPercentLevelCrossing(x, y, 70, 1, i0);
	u75 = findPercentLevelCrossing(x, y, 75, 1, i0);
	u80 = findPercentLevelCrossing(x, y, 80, 1, i0);
	u90 = findPercentLevelCrossing(x, y, 90, 1, i0);
	
	// Find downstroke times
	if (invert_downstroke_search) {
		i0 = i100;
		kdir = 1;
		}
	else {
		i0 = y.length-2;
		kdir = -1;
		}
	d90 = findPercentLevelCrossing(x, y, 10, kdir, i0);
	d80 = findPercentLevelCrossing(x, y, 20, kdir, i0);
	d75 = findPercentLevelCrossing(x, y, 25, kdir, i0);
	d70 = findPercentLevelCrossing(x, y, 30, kdir, i0);
	d60 = findPercentLevelCrossing(x, y, 40, kdir, i0);
	d50 = findPercentLevelCrossing(x, y, 50, kdir, i0);
	d40 = findPercentLevelCrossing(x, y, 60, kdir, i0);
	d30 = findPercentLevelCrossing(x, y, 70, kdir, i0);
	d25 = findPercentLevelCrossing(x, y, 75, kdir, i0);
	d20 = findPercentLevelCrossing(x, y, 80, kdir, i0);
	d10 = findPercentLevelCrossing(x, y, 90, kdir, i0);
	
	if (plotflag) {
		if (nth>0)
			plot_title = "Transient "+i2s(nth);
		else {
			if (smot)
				plot_title = "Average transient (smoothed)";
			else
				plot_title = "Average transient";
			}
		Plot.create(plot_title, "Time (s)", "y", x, yf);
		Plot.setColor("red");
		Plot.drawLine(x100, 0, x100, y100);
		plotHorzLine(x0, xlast, 0, "blue");
		plotHorzLine(u10, d90, 0.1*y100, "#4400bb"); 
		plotHorzLine(u20, d80, 0.2*y100, "#5500bb"); 
		plotHorzLine(u25, d75, 0.25*y100, "#007744"); 
		plotHorzLine(u30, d70, 0.3*y100, "#7700bb"); 
		plotHorzLine(u40, d60, 0.4*y100, "#aa00bb"); 
		plotHorzLine(u50, d50, 0.5*y100, "#bb00aa"); 
		plotHorzLine(u60, d40, 0.6*y100, "#bb0077"); 
		plotHorzLine(u70, d30, 0.7*y100, "#bb0066"); 
		plotHorzLine(u75, d25, 0.75*y100, "#447700"); 
		plotHorzLine(u80, d20, 0.8*y100, "#bb0055"); 
		plotHorzLine(u90, d10, 0.9*y100, "#bb0044"); 
		Plot.setColor("black");
		Plot.setFrameSize(250, 300);
		Plot.show();
		savePlot(plot_title);
		setBatchMode("show");
		}

	// Add calculated parameters to Results table window
	if (rpms)
		m = 1000;
	else
		m = 1;
	k = nResults;
	setResult("T0", k, m*x0);
	setResult("Interval", k, m*(xlast-x0));
	setResult("Amp", k, y100);
	setResult("TTP", k, m*(x100-x0));
	setResult("Up90", k, m*(u90-u10));
	setResult("Dn90", k, m*(d90-d10));
	setResult("CD10", k, m*(d90-u10));
	setResult("CD25", k, m*(d75-u25));
	setResult("CD50", k, m*(d50-u50));
	setResult("CD75", k, m*(d25-u75));
	setResult("CD90", k, m*(d10-u90));
	if (SHOW_EXTENDED_PARAMETERS) {
		setResult("TU10", k, m*(u10-x0));
		setResult("TU20", k, m*(u20-x0));
		setResult("TU25", k, m*(u25-x0));
		setResult("TU30", k, m*(u30-x0));
		setResult("TU40", k, m*(u40-x0));
		setResult("TU50", k, m*(u50-x0));
		setResult("TU60", k, m*(u60-x0));
		setResult("TU70", k, m*(u70-x0));
		setResult("TU75", k, m*(u75-x0));
		setResult("TU80", k, m*(u80-x0));
		setResult("TU90", k, m*(u90-x0));
		setResult("TPk", k, m*(x100-x0));
		setResult("TD10", k, m*(d10-x0));
		setResult("TD20", k, m*(d20-x0));
		setResult("TD25", k, m*(d25-x0));
		setResult("TD30", k, m*(d30-x0));
		setResult("TD40", k, m*(d40-x0));
		setResult("TD50", k, m*(d50-x0));
		setResult("TD60", k, m*(d60-x0));
		setResult("TD70", k, m*(d70-x0));
		setResult("TD75", k, m*(d75-x0));
		setResult("TD80", k, m*(d80-x0));
		setResult("TD90", k, m*(d90-x0));
		}
		
	atimes = newArray(x0, u10, u50, x100);
	return atimes;
	}

function findPercentLevelCrossing(x, y, pc, kdir, i0) {
	if((pc<1)||(pc>99))
		exit("Percent level must be in range 1-99");
	// Calculate level from data and pc (%)
	Array.getStatistics(y, ymin, y100, ymean, ystd);
	// ymin will be zero for baseline-corrected trace
	level = y100*(pc/100.0);
	sign0 = signDiff(y[i0], level);
	found = 0;
	if (kdir==1) {		// search forwards from starting index
		i = i0;
		while (i<y.length && found==0) {
			if (signDiff(y[i], level) != sign0) {
				if (y[i-1]==y[i])
					xcross = 0;
				else
					xcross = x[i] - (y[i]-level) * (x[i]-x[i-1]) / (y[i]-y[i-1]);
				found = 1;
			}
			i++;
		}
	} else {			// search backwards from starting index
		i = i0;
		while (i>=0 && found==0) {
			if (signDiff(y[i], level) != sign0) {
				if (y[i+1]==y[i])
					xcross = 0;
				else
					xcross = x[i] - (level-y[i]) * (x[i+1]-x[i]) / (y[i]-y[i+1]);
				found = 1;
			}
			i--;
		}
	}
	if (found==0)
		exit("findPercentLevelCrossing failed");

	return xcross;
	}
	
function signDiff(a, b) {
	if (a < b)
		return -1;
	else if (a > b)
		return 1;
	else
		return 0;
	}

function smoothData(y) {
    yf = filterMedian5(y);
    yf = filterMedian3(yf);
    yf = filterHanning(yf);
    return yf;
    }
    
function filterMedian5(y) {
    n = y.length;
    yfilt = newArray(n);
    
    yfilt[0] = y[0];
    yfilt[1] = y[1];
    yfilt[n-2] = y[n-2];
    yfilt[n-1] = y[n-1];
    
    for (i=2; i<=n-3; i++)
     	yfilt[i] = getMedian5(y[i-2], y[i-1], y[i], y[i+1], y[i+2]);
    
    return yfilt;
    }
    
function getMedian5(y0, y1, y2, y3, y4) {
    if (y0 > y1) { t = y0; y0 = y1; y1 = t; }
    if (y3 > y4) { t = y3; y3 = y4; y4 = t; }
    if (y0 > y3) { t = y0; y0 = y3; y3 = t; }
    if (y1 > y4) { t = y1; y1 = y4; y4 = t; }
    if (y1 > y2) { t = y1; y1 = y2; y2 = t; }
    if (y2 > y3) { t = y2; y2 = y3; y3 = t; }
    if (y1 > y2) { t = y1; y1 = y2; y2 = t; }
    return y2;
    }
   
function filterMedian3(y) {
    n = y.length;
    yfilt = newArray(n);
    
    yfilt[0] = y[0];
    yfilt[n-1] = y[n-1];
    
    for (i=1; i<=n-2; i++)
     	yfilt[i] = getMedian3(y[i-1], y[i], y[i+1]);
    
    return yfilt;
    }
    
function getMedian3(y0, y1, y2) {
    if (y0 < y1) {
     	if (y0 > y2)
    		return y0;
    	else {
      	  	if (y1 > y2)
          		return y2;
        	else
          	return y1;
        	}
      	}
    else {
    	if (y1 > y2)
    		return y1;
    	else {
        	if (y0 > y2)
          		return y2;
        	else
    			return y0;
    		}
    	}
    }
   
function filterHanning(y) {
    n = y.length;
    yfilt = newArray(n);
    
    yfilt[0] = y[0];
    yfilt[n-1] = y[n-1];
    
    for (i=1; i<=n-2; i++)
     	yfilt[i] = y[i]*0.5 + (y[i-1]+y[i+1])*0.25;
    
    return yfilt;
    }
 
function getBlips(thresholdpc) {
	selectWindow("Seq");
	run("Select All");
	run("Plot Z-axis Profile");
	Plot.getValues(x, z);

	n = lengthOf(z);
	Array.getStatistics(z, min, max, mean, sd);

	if (show_intermediate_results) {
		print("Whole frame amplitude: " + i2s(min) + "-" + i2s(max) + " = " + i2s(max-min));
		xm = newArray(0);
		plotTrace("Whole frame amplitude", x, z, xm, "");
		savePlot("Amplitude");
		setBatchMode("show");
		}
		
	if (max-min < MIN_BLIP_AMPLITUDE)
		exit("No blips detected (too small amplitude)!");
		
	// Zero baseline (assuming baseline is level)
	for (i=0; i<n; i++)
		z[i] = z[i] - min;
	max -= min;
	mean -= min;
	min -= min;
	
	// Invert if blips are downward going
	if (abs(mean-min) > abs(mean-max)) {
		for (i=0; i<n; i++)
			z[i] = max - z[i];
		}
	blips = Array.findMaxima(z, max*thresholdpc/100);
	Array.sort(blips);	// because findMaxima returns in peak size not time order
	return blips;
	}
    
function checkBlipGoodness(blips, frame_interval) {
	n = blips.length;
	if (n == 0)
		exit("No blips detected (peak detection failed)!");
	else if (n == 1)
		exit("Only one blip detected!");
	else {
		blip_intervals = newArray(n-1);
		for (i=0; i<n-1; i++)
			blip_intervals[i] = (blips[i+1] - blips[i]) * frame_interval * 1000;
		Array.getStatistics(blip_intervals, bmin, bmax, bmean, bsd);
		bmed = arrayMedian(blip_intervals);
		minpcdev = (bmed-bmin)/bmed*100;
		maxpcdev = (bmax-bmed)/bmed*100;
		if (show_intermediate_results) {
			print("Blips detected: " + i2s(n));
			print("Min interval: " + i2s(bmin) + "ms");
			print("Max interval: " + i2s(bmax) + "ms");
			print("Median interval: " + i2s(bmed) + "ms");
			print("Short deviation: " + d2s(minpcdev, 1) + "%");
			print("Long deviation: " + d2s(maxpcdev, 1) + "%");
			}
			
		if (bmin < MIN_BLIP_INTERVAL)
			print("Blip interval " + i2s(bmin) + "ms less than minimum " + i2s(MIN_BLIP_INTERVAL));
		if (bmax < MIN_BLIP_INTERVAL)
			print("Blip interval " + i2s(bmax) + "ms greater than maximum " + i2s(MAX_BLIP_INTERVAL));
		if (minpcdev > BLIP_TOLERANCE_PC)
			print("One or more short blip intervals detected");
		if (maxpcdev > BLIP_TOLERANCE_PC)
			print("One or more long blip intervals detected");
		}
	}
    
function arrayMedian(y) {
	n = y.length;
	if (n == 0)
		return 0;
	else if (n == 1)
		return y[0];
	else {
		a = Array.copy(y);	
		Array.sort(a);
		k = floor(n/2);
		if (n%2 == 0)
			return (a[k-1] + a[k]) / 2;
		else
			return a[k];
		}
	}
	
function getPercentileLevel(pc) {
	resetMinAndMax(); 
	getHistogram(values, counts, 256); 
	
	// find cumulative sum 
	npixels = 0; 
	for (i = 0; i<counts.length; i++) 
		npixels += counts[i];
	  
	nbelowthreshold = npixels * pc/100; 
	sum = 0; 
	for (i = 0; i<counts.length; i++) { 
		sum = sum + counts[i]; 
		if (sum >= nbelowthreshold) {
			return(values[i]);
	  		} 
	  	}
	}
   
function closeWindow(name) {
	selectWindow(name);
	close();
	}
	
function closeSpecialWindow(name) {
	if (isOpen(name)) { 
		selectWindow(name); 
		run("Close"); 
		}
	} 

function readSettingValue(name, defval) {
	s = call("ij.Prefs.get", "contractility."+name, "?");
	if (s=="?")
		return defval;
	else
		return parseFloat(s);
	}

function writeSettingValue(name, val, ndecpl) {
	call("ij.Prefs.set", "contractility."+name, d2s(val,ndecpl));
	}

function readSettingBoolean(name, defval) {
	s = call("ij.Prefs.get", "contractility."+name, "?");
	if (s=="?")
		return defval;
	else if (s=="0")
		return false;
	else
		return true;
	}

function writeSettingBoolean(name, val) {
	s = "contractility."+name;
	if (val==false)
		call("ij.Prefs.set", s, "0");
	else
		call("ij.Prefs.set", s, "1");
	}

// perpendicular distance from point (x0,y0) to line joining (x1,y1) and (x2,y2)
function perpDistance(x1, y1, x2, y2, x0, y0) {
	Assert(((x1!=x2)||(y1!=y2)), "Points are equal in function perpDistance");
	return abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1) / sqrt((y2-y1)*(y2-y1) + (x2-x1)*(x2-x1));
	}
	
function inputValues(prompt, endval, vmin, vmax) {
	vals = newArray(0);
	done = false;
	while (!done) {
		v = getNumber(prompt, endval);
		if (v==endval)
			done = true;
		else if (v>=vmin && v<=vmax)
			vals = Array.concat(vals, v);
		else
			waitForUser("Value out of range "+i2s(vmin)+".."+i2s(vmax));
		}
	return vals;
	}
	
function Assert(cond, msg) {
	if (cond==0)
		exit("Assertion failed: " + msg);
	}

function getXPoints() {
	Plot.getLimits(xmin, xmax, ymin, ymax);
	Plot.setFrameSize(800, 300);
	// We assume a fixed size border round plot axes for set frame size
	LEFT = 60;
	RIGHT = 860;
	BOTTOM = 313;
	TOP = 13;
	leftButton = 16;
	rightButton = 4;
	xfac = (RIGHT-LEFT) / (xmax-xmin);
	xoff = LEFT - xmin*xfac;
	yfac = (TOP-BOTTOM) / (ymax-ymin);
	yoff = BOTTOM - ymin*yfac;
	run("Select None");
	showMessageWithCancel("Mark transient start times", "Left click to mark; right click to delete; click outside axes to finish");
	setOption("DisablePopupMenu", true);
	setOption("ExpandableArrays", true);
	xpoints = newArray; 
	ypoints = newArray; 
	npoints = 0;
	done = false;
	getCursorLoc(x2, y2, z, flags);
	moved = true;
	while (!done) {
		getCursorLoc(x, y, z, flags);
		ns = d2s(npoints,0) + " point";
		if (npoints!=1)
			ns = ns + "s";
		showStatus("("+x+","+y+") "+ns+" marked");
		if (flags&rightButton!=0 && moved) {	// delete last selected point
			if (npoints>0) {
				npoints--;
				xpoints = Array.slice(xpoints, 0, xpoints.length-1);
				ypoints = Array.slice(ypoints, 0, ypoints.length-1);
				makeSelection("point", xpoints, ypoints);
				}
			else {
				npoints = 0;
				xpoints = newArray; 
				ypoints = newArray; 
				run("Select None");
				}
			x2 = x;
			y2 = y;
			moved = false;
			wait(50);
			}
		else if (flags&leftButton!=0 && moved) {	// add new point
			if (x<LEFT-2 || x>RIGHT+2 || y<TOP-2 || y>BOTTOM+2)
				done = true;
			else {
				xpoints[npoints] = x;
				ypoints[npoints] = y;
				npoints++;
				x2 = x;
				y2 = y;
				moved = false;
				makeSelection("point", xpoints, ypoints);
				wait(50);
				}
			}
		else if (flags&leftButton==0) {
			// prevent multiple points or deletes if mouse moved while button down
			moved = ((x!=x2) || (y!=y2));
			}
		}

	run("Select None");
	setOption("DisablePopupMenu", false);
	setOption("ExpandableArrays", false);

	// convert from pixels to plot coordinates
	for (i=0; i<npoints; i++) {
		xpoints[i] = (xpoints[i]-xoff) / xfac;
		ypoints[i] = (ypoints[i]-yoff) / yfac;
	}

	return xpoints;
	}

function savePlot(name) {
	if (SAVE_PLOTS) {
		fn = FILE_ROOT + "_" + name + ".png";
		saveAs("png", fn);
		}
	}