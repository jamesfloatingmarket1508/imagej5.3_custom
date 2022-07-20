// ImageJ macro to convert multiple TIFF files that comprise an image sequence
// into a single TIFF stack file. Do this for all folders in a tree selected by
// the user.
//
// This macro is intended to allow image sequences collected by HCIImageLive
// to be used in ContractilityTool's batch processing procedure.
//
// FB 30/11/17

requires("1.50c");

inDir = getDirectory("Choose root folder to look for TIFF sequences");
outDir = getDirectory("Choose root folder for saving TIFF stacks");
if (inDir==outDir)
	exit("Input and output folders must be different!");
	
setBatchMode(true);
print("Running macro to convert image sequences to TIFF stacks");
minSeq = 10;
tiffExt = ".TIF";
nseq = 0;
processFiles(inDir, outDir, "");

if (nseq==1)
	plural = "";
else
	plural = "s";
print("Converted " + nseq + " sequence" + plural + " to stack format");
print("Finished processing");

function processFiles(inBase, outBase, subFolder) {
	if (subFolder!="")
		print("Scanning folder: " + subFolder);
		
	list = getFileList(inBase + subFolder);
	list = Array.sort(list);
	File.makeDirectory(outBase + subFolder);
	
	// Count TIFF files and get index of first one
	ntiff = 0;
	tiff1 = -1;
	for (i=0; i<list.length; i++) {
		if (endsWith(toUpperCase(list[i]), tiffExt)) {
			ntiff++;
			if (tiff1<0)
				tiff1 = i;
			}
		}
		
	// Assume folder contains image sequence if many TIFF files present
	if (ntiff>minSeq)
		processSequence(inBase, outBase, subFolder, list[tiff1], ntiff);
		
	// Continue to look for subfolders
	for (i=0; i<list.length; i++) {
		path = subFolder + list[i];
		if (endsWith(path, "/")) {
			// recurse into subdirectories
			processFiles(inBase, outBase, path);
	    	}
		}
	}

function processSequence(in, out, sub, firstfile, nfiles) {
	inspec = in + sub + firstfile;
	print("Reading sequence of " + nfiles + " TIFF files");
	cmd = "open=[" + inspec + "] number=" + nfiles;
	cmd = cmd + " starting=1 increment=1 sort"; 
	run("Image Sequence...", cmd);
	if (bitDepth()!=16) {
		print("Converting to 16 bit");
		run("32-bit");
		run("16-bit");
		}
	outspec = out + sub + makeCtFileName(firstfile);
	saveAs("Tiff", outspec);
	close();
	nseq++;
	}
	
function makeCtFileName(fn) {
	fn = replace(fn, " ", "");
	parts = split(fn, ".");
	np = lengthOf(parts);
	s = noRaw(parts[0]);
	// replace any dots in filename with underscores
	if (np>2) {
		for (i=1; i<np-1; i++)
			s = s + "_" + noRaw(parts[i]);
		}
	s = s + "_Raw.tif";
	return s;
	}
	
function noRaw(s) {
	return replace(s, "_Raw", "");
	}
