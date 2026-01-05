// -----------------------------------------------------------------------------
// ImageCytox_Quantification.ijm
// Batch quantification of endpoint cytotoxicity assay images (Fiji/ImageJ macro).
//
// Input:
//   - Multichannel TIFF hyperstacks (one file per field of view).
//   - One TOP folder containing one subfolder per experimental condition.
//   - A TOP/Results folder is created automatically for outputs (and is ignored on re-runs).
//
// Output written to TOP/Results/:
//   - ImageCytox_Results.csv
//   - Segmentation_QC_stack.tif   (RGB stack with raw channels + classification outlines)
//
// QC stack contents:
//   - Background (raw channels):
//       nuclei = cyan, PI = magenta, optional T marker = yellow
//   - Classification outlines (ROI Manager stroke color):
//       live tumor = white (line width 2)
//       dead tumor (PI+) = magenta (line width 1)
//       T cells = orange (line width 1) [only if T marker channel is provided]
//       fragments/unclassified = gray (line width 1)
//
// Requirements:
//   - Fiji/ImageJ with StarDist 2D + CSBDeep installed.
// -----------------------------------------------------------------------------

// -------------------- global variables ---------------------------------------
nucleiChan    = 0;     // Hoechst channel index (1-based)
piChan        = 0;     // PI channel index (1-based)
tChan         = 0;     // optional T marker channel index (1-based; 0 = none)

tCellFactor   = 1.2;   // T cell marker positivity factor (mean / background); only used if T marker channel present
minTumorArea  = 0;     // minimum tumor nucleus area (pixels)
maxNucleusArea = 0;     // optional maximum nucleus area (pixels; 0 = off)
doSubtractBG   = false; // optional: Subtract Background on nuclei channel before StarDist
bgRolling      = 50;    // rolling ball radius for background subtraction (pixels)
doCondensed    = false; // optional: exclude condensed Hoechst-bright nuclei
piFactor      = 0;     // dead if meanPI >= backgroundMean * piFactor
hoechstBrightFactor = 0; // 0 = off; Hoechst-bright, smallish nuclei are treated as "Fragments / other"

csvPath       = "";
topResultsDir = "";

resultsStackTitle = "Segmentation_QC_stack";

// -----------------------------------------------------------------------------
// MAIN MACRO
// -----------------------------------------------------------------------------
macro "ImageCytox_Quantification" {

    topDir = getDirectory("Select TOP folder with condition subfolders");
    if (topDir == "") exit("No folder selected.");

    topResultsDir = topDir + "Results" + File.separator;
    File.makeDirectory(topResultsDir);

    Dialog.create("ImageCytox quantification");

    // -------------------- required inputs --------------------
    Dialog.addMessage("Required channels (1-based indices as shown in Fiji). Use 0 if a channel is not present.");
    Dialog.addString("Nuclei channel (Hoechst) index", "1");
    Dialog.addString("PI channel index (dead nuclei; 0 = none)", "2");

    Dialog.addMessage("Optional marker channels");
    Dialog.addString("T cell marker channel index (0 = none)", "0");
    Dialog.addNumber("T cell positivity factor (mean / background)", 1.2);
    Dialog.addNumber("Tcell max nucleus area (px) (0 = auto)", 0);

    // -------------------- nuclei segmentation / gating --------------------
    Dialog.addMessage("Nuclei size gating (pixels)");
    Dialog.addNumber("Minimum nucleus area (px)", 600);
    Dialog.addNumber("Maximum nucleus area (px) (0 = off)", 0);

    // -------------------- PI classification --------------------
    Dialog.addMessage("Dead nuclei classification (PI)");
    Dialog.addNumber("PI dead threshold factor (mean / background)", 1.2);

    // -------------------- optional preprocessing --------------------
    Dialog.addMessage("Optional preprocessing (recommended for dim nuclei or high background)");
    Dialog.addCheckbox("Subtract background (rolling ball) before segmentation", true);
    Dialog.addNumber("Rolling ball radius (px)", 100);

    // -------------------- optional exclusion rules --------------------
    Dialog.addMessage("Optional exclusion: condensed Hoechst-bright nuclei");
    Dialog.addCheckbox("Exclude condensed Hoechst-bright nuclei (small but very bright; apoptotic / T cells)", true);
    Dialog.addNumber("Condensed brightness factor (× reference mean of largest nuclei)", 1.7);

    Dialog.show();

    nucleiChan   = parseInt(Dialog.getString());
    piChan       = parseInt(Dialog.getString());
    tChan        = parseInt(Dialog.getString());
    tCellFactor  = Dialog.getNumber();
    tCellMaxArea = Dialog.getNumber();
    minTumorArea = Dialog.getNumber();
    maxNucleusArea = Dialog.getNumber();
    if (tCellMaxArea <= 0) tCellMaxArea = minTumorArea / 2;
    piFactor     = Dialog.getNumber();
    doSubtractBG = Dialog.getCheckbox();
    bgRolling    = Dialog.getNumber();
    doCondensed  = Dialog.getCheckbox();
    hoechstBrightFactor = Dialog.getNumber();
    if (!doCondensed) hoechstBrightFactor = 0;


    if (nucleiChan < 1) exit("Nuclei channel index must be >= 1.");
    if (piChan < 1)     exit("PI channel index must be >= 1.");
    if (tChan < 0)      exit("T cell marker channel index must be >= 0.");

    // Create a unique CSV path (avoid overwriting and avoid Excel file locks)
    baseCsv = topResultsDir + "ImageCytox_Results";
    csvPath = baseCsv + ".csv";
    idx = 1;
    while (File.exists(csvPath)) {
        csvPath = baseCsv + "_" + idx + ".csv";
        idx++;
    }
    if (File.exists(csvPath)) File.delete(csvPath);

    // Column names chosen for broad, biologist-friendly compatibility
    File.append("Condition,Position,Survivors,,TotalROIs,Excluded:,DeadNuclei,T-Cells,HoechstBright,Fragments/Other,,SanityCheck:,NucleiMeanIntensity,NucleiIntThreshold,PI_BG_Intensity,MedianNucleiSize\n", csvPath);
setBatchMode(true);
    debugPrinted = false;

    condList = getFileList(topDir);
    for (c = 0; c < condList.length; c++) {

        condName = condList[c];
        condPath = topDir + condName;

        if (!File.isDirectory(condPath)) continue;

        // Robustly exclude TOP/Results (handles trailing separators and odd naming)
        condLabel = normalizeFolderName(getFolderName(condPath));
        condNameClean = normalizeFolderName(condName);
        condLower = toLowerCase(condLabel);
        nameLower = toLowerCase(condNameClean);
        if (condLower == "results" || startsWith(condLower, "results")) continue;
        if (nameLower == "results" || startsWith(nameLower, "results")) continue;

        if (!endsWith(condPath, File.separator))
            condPath = condPath + File.separator;

        print("=== Condition: " + condLabel + " ===");
        processConditionFolder(condPath, condLabel);
    }

    setBatchMode(false);

    if (isOpen(resultsStackTitle)) {
        selectWindow(resultsStackTitle);
        saveAs("Tiff", topResultsDir + "Segmentation_QC_stack.tif");
        print("Saved QC stack to: " + topResultsDir + "Segmentation_QC_stack.tif");
    } else {
        print("Warning: QC stack was not created (no images processed?).");
    }

    closeAllButLogAndQCStack();
    print("ImageCytox_Quantification finished.");
}


// -----------------------------------------------------------------------------
// processConditionFolder
// -----------------------------------------------------------------------------
function processConditionFolder(folder, condLabel) {

    // Safety: never process a Results folder
    if (endsWith(toLowerCase(folder), File.separator + "results" + File.separator) ||
        endsWith(toLowerCase(folder), "results" + File.separator))
        return;

    fileList = getFileList(folder);
    images = getFilesByExtension(fileList, ".tif");
    images = appendArray(images, getFilesByExtension(fileList, ".tiff"));

    if (images.length == 0) {
        print("  No .tif or .tiff images in folder: " + folder);
        return;
    }

    for (i = 0; i < images.length; i++) {

        filename = images[i];
        fullPath = folder + filename;

        // Skip directories (getFileList can return subfolders with trailing separator)
        if (File.isDirectory(fullPath)) continue;

        // Safety: never analyze files inside any Results subfolder
        lp = toLowerCase(fullPath);
        if (indexOf(lp, "\\results\\") >= 0 || indexOf(lp, "/results/") >= 0 ||
            indexOf(lp, File.separator + "results" + File.separator) >= 0)
            continue;

        print("  Processing: " + fullPath);

        open(fullPath);
        origTitle = getTitle();
        baseName  = stripExtension(filename);

        nucTitle = baseName + "_Nuc";
        piTitle  = baseName + "_PI";
        tTitle   = "";

        // Split channels (window names can vary across Fiji/OS; Windows may show "C1-\well-...tif")
        selectWindow(origTitle);
        run("Split Channels");
        wait(150);

        nucWin = findSplitWindowFlexible(nucleiChan, baseName);
        piWin  = findSplitWindowFlexible(piChan,     baseName);

        if (nucWin == "") exit("Could not find nuclei split window for baseName '" + baseName + "'.");
        if (piWin  == "") exit("Could not find PI split window for baseName '" + baseName + "'.");

        selectWindow(nucWin); rename(nucTitle);
        selectWindow(piWin);  rename(piTitle);

        if (tChan > 0) {
            tWin = findSplitWindowFlexible(tChan, baseName);
            if (tWin != "") {
                tTitle = baseName + "_T";
                selectWindow(tWin); rename(tTitle);
            } else {
                print("    Warning: T marker channel " + tChan + " not found for " + filename + " (continuing without marker).");
                tTitle = "";
            }
        }

        // Close any remaining split windows related to this file (channels we did not keep)
        closeSplitWindowsForBase(baseName, nucTitle, piTitle, tTitle);

        // Close original image if still open
        if (isOpen(origTitle)) { selectWindow(origTitle); run("Close"); }

        // --- StarDist segmentation (nuclei ROIs) -----------------------------
        // Optional preprocessing: subtract background on nuclei channel before StarDist
        if (doSubtractBG) {
            selectWindow(nucTitle);
            run("Subtract Background...", "rolling=" + bgRolling);
            if (tTitle != "") { selectWindow(tTitle); run("Subtract Background...", "rolling=" + bgRolling); }
        }
        runStarDistOnNuclei(nucTitle);

        // Background in T-marker channel (mean intensity outside ALL nuclei ROIs)
        tBG = 0;
        if (tTitle != "") { tBG = estimateBackgroundOutsideROIs(tTitle); }

        // Optional max nucleus area: do NOT delete large ROIs.
        // Large ROIs are kept and displayed as thin gray outlines in the QC overlay, and counted as excluded.
        totalROI = roiManager("count");
        // Median (50th percentile) of tumor-eligible nucleus areas in this image (used only if Hoechst-bright rule is enabled)
        // Median nucleus ROI area across ALL detected nuclei (includes small ones & fragments)
        medianNucleiSize = computeAllNucleiAreaMedian(nucTitle, maxNucleusArea);
        // Median (50th percentile) tumor-eligible nucleus area (used to define the "smaller half" for condensed-nuclei rule)
        medianTumorArea = computeTumorAreaMedian(nucTitle, tTitle, minTumorArea, maxNucleusArea, tBG, tCellFactor);
        // Reference intensity: mean Hoechst intensity of the TOP 25% largest tumor-eligible nuclei (by area)
        nucRefMeanTop25 = computeRefMeanNucleiIntensityTop25(nucTitle, tTitle, minTumorArea, maxNucleusArea, tBG, tCellFactor);

        // Sanity-check outputs
        nucleiMeanIntensity = nucRefMeanTop25;
        nucleiIntThreshold = "";
        if (hoechstBrightFactor > 0 && nucRefMeanTop25 > 0) {
            nucleiIntThreshold = "" + (nucRefMeanTop25 * hoechstBrightFactor);
        }
        if (!debugPrinted && hoechstBrightFactor > 0) {
            print("Condensed Hoechst-bright exclusion enabled: factor=" + hoechstBrightFactor + " | medianArea=" + medianTumorArea + " | refMeanTop25=" + nucRefMeanTop25);
            debugPrinted = true;
        }

        selectWindow(nucTitle);
        run("Select All");
        getStatistics(_aNimg, nucRefMedian, _miNimg, _maNimg, _sNimg);
        run("Select None");

        if (totalROI == 0) {
            print("    No nuclei detected (StarDist produced 0 ROIs). Skipping image: " + baseName);
            closeAllButLogAndQCStack();
            continue;
        }

        // --- Robust PI background estimate: pixels outside nuclei ROIs -------
        piBG = estimatePIBackgroundOutsideROIs(piTitle);

        // --- QC background creation (merge DUPLICATES; never touch originals) -
        nucDisp = baseName + "_Nuc_disp";
        piDisp  = baseName + "_PI_disp";
        tDisp   = "";

        selectWindow(nucTitle);
        run("Select None");
        run("Duplicate...", "title=" + nucDisp);

        selectWindow(piTitle);
        run("Select None");
        run("Duplicate...", "title=" + piDisp);

        if (tTitle != "") {
            tDisp = baseName + "_T_disp";
            selectWindow(tTitle);
            run("Select None");
            run("Duplicate...", "title=" + tDisp);
        }

        qcTitle = baseName + "_QC_RGB";

        // Merge Channels requires identical canvas sizes. If sizes differ, skip multi-channel merge for this image.
        selectWindow(nucDisp); wN = getWidth(); hN = getHeight();
        selectWindow(piDisp);  wP = getWidth(); hP = getHeight();

        if (wN != wP || hN != hP) {
            print("    Warning: channel canvas size mismatch for " + baseName + " (Nuc " + wN + "x" + hN + ", PI " + wP + "x" + hP + ").");
            print("    Skipping multi-channel QC merge for this image (quantification will continue).");

            // Fallback QC: nuclei only
            selectWindow(nucDisp);
            run("8-bit");
            run("Enhance Contrast", "saturated=0.35");
            run("Cyan");
            run("RGB Color");
            rename(qcTitle);

        } else {
            if (tDisp != "") {
                selectWindow(tDisp); wT = getWidth(); hT = getHeight();
                if (wT != wN || hT != hN) {
                    print("    Warning: T channel canvas size mismatch for " + baseName + " (T " + wT + "x" + hT + ", expected " + wN + "x" + hN + ").");
                    print("    Skipping T channel in QC merge for this image.");
                    run("Merge Channels...", "c1=[" + nucDisp + "] c2=[" + piDisp + "] create");
                } else {
                    run("Merge Channels...", "c1=[" + nucDisp + "] c2=[" + piDisp + "] c3=[" + tDisp + "] create");
                }
            } else {
                run("Merge Channels...", "c1=[" + nucDisp + "] c2=[" + piDisp + "] create");
            }
            rename(qcTitle);
        }

        // Close display duplicates we no longer need (qcTitle is created)
        if (isOpen(nucDisp)) { selectWindow(nucDisp); run("Close"); }
        if (isOpen(piDisp))  { selectWindow(piDisp);  run("Close"); }
        if (tDisp != "" && isOpen(tDisp)) { selectWindow(tDisp); run("Close"); }

        // Standardize QC display and convert to true RGB pixels for robust drawing
        selectWindow(qcTitle);
        run("Make Composite");

        // Contrast adjustment is for visualization only.
        // Important: do NOT auto-stretch PI (or T marker) when the signal is absent; it inflates background noise.
        // We only enhance the nuclei channel for a consistent QC background.
        Stack.setChannel(1); run("Enhance Contrast", "saturated=0.35");
Stack.setChannel(1); run("Cyan");
        Stack.setChannel(2); run("Magenta");
        if (tTitle != "") { Stack.setChannel(3); run("Yellow"); }

        run("RGB Color");
        setOption("BlackBackground", false);

        survivors = 0;
        deadNuclei = 0;
        tCells    = 0;
        fragments = 0;
        otherCount = 0;
        hoechstBrightCount = 0;


        // Per-ROI styling for QC overlay (applied after classification)
        roiColors = newArray(totalROI);
        roiWidths = newArray(totalROI);
        // Default outline thickness = 1; only live tumor ROIs will be drawn with width 3
        defaultLW = 1;

        // --- classification + ROI drawing -----------------------------------
        for (n = 0; n < totalROI; n++) {

            // Nucleus measurements (size gate)
            selectWindow(nucTitle);
            roiManager("select", n);
            getStatistics(areaNuc, meanNuc, minNuc, maxNuc, stdNuc);

            // Exclude overly large ROIs (but keep them for QC overlay as thin gray outlines)
            if (maxNucleusArea > 0 && areaNuc > maxNucleusArea) {
                otherCount++;
                roiColors[n] = "gray";
                roiWidths[n] = 1;
                continue;
            }

            // PI measurements for same ROI
            selectWindow(piTitle);
            roiManager("select", n);
            getStatistics(areaPI, meanPI, minPI, maxPI, stdPI);

            // Optional T marker intensity
            tMean = 0;
            if (tTitle != "") {
                selectWindow(tTitle);
                roiManager("select", n);
                getStatistics(areaT, meanT, minT, maxT, stdT);
                tMean = meanT;
            }

            isTcell = false;
            isDeadTumor = false;
            isLiveTumor = false;
            isHoechstBright = false;
            isOther = false;

            // T cell exclusion:
            // If NO T marker channel is provided (tChan = 0), do NOT classify T cells by size alone.
            if (tTitle != "" && tCellFactor > 0 && tBG > 0) {
                if (tMean >= (tBG * tCellFactor)) {
                    if (areaNuc <= tCellMaxArea) isTcell = true;
                }
            }

            // Tumor classification (only if not T cell):
            // - only consider as tumor if nucleus is large enough
            // - dead tumor: PI mean above threshold
            // - live tumor: PI mean below threshold
            if (!isTcell) {
                if (areaNuc >= minTumorArea) {

                    // Primary dead call: PI-positive nuclei
                    if (meanPI >= piBG * piFactor) {
                        isDeadTumor = true;
                    } else {

                        // Optional exclusion of condensed nuclei (OFF by default):
                        // Smallish nuclei (lower 50% by area among tumor-eligible nuclei in this image) that are
                        // unusually Hoechst-bright are treated as "Fragments / other" (gray) rather than surviving tumor.
                        // Threshold: meanHoechst >= (referenceMeanTop25Largest * factor).
                        if (hoechstBrightFactor > 0 && medianTumorArea > 0 &&
                            areaNuc <= medianTumorArea &&
                            nucRefMeanTop25 > 0 && meanNuc >= (nucRefMeanTop25 * hoechstBrightFactor)) {

                            isHoechstBright = true; // separate excluded category
                        } else {
                            isLiveTumor = true;
                        }
                    }
                }
            }
            // Store ROI styling for QC overlay (we burn overlays via Show All + Flatten later)
            if (isTcell) {
                tCells++;
                roiColors[n] = "orange";
                roiWidths[n] = 1;
            } else if (isDeadTumor) {
                deadNuclei++;
                roiColors[n] = "magenta";
                roiWidths[n] = 1;
            } else if (isHoechstBright) {
                hoechstBrightCount++;
                roiColors[n] = "cyan";
                roiWidths[n] = 1;
            } else if (isLiveTumor) {
                survivors++;
                roiColors[n] = "white";
                roiWidths[n] = 2;
            } else {
                // Excluded: too small, too large (if max enabled), or otherwise unclassifiable
                // For reporting, fragments+otherCount are combined as "Fragments/Other" in the CSV.
                if (areaNuc < minTumorArea) fragments++;
                else otherCount++;
                roiColors[n] = "gray";
                roiWidths[n] = 1;
            }
}

        // --- outputs ----------------------------------------------------------
        
        // Apply per-ROI stroke settings for QC overlay (matches manual workflow: set properties -> Show All -> Flatten)
        for (n = 0; n < totalROI; n++) {
            roiManager("select", n);
            roiManager("Set Color", roiColors[n]);
            roiManager("Set Line Width", roiWidths[n]);
            roiManager("Update");
        }

        // Activate overlay display on the QC image (burned in appendToQCStack via Flatten)
        selectWindow(qcTitle);
        roiManager("Show All");

        csvLine = condLabel + "," + baseName + "," + survivors + ",," +
            totalROI + "," + "" + "," + deadNuclei + "," + tCells + "," + hoechstBrightCount + "," + (fragments + otherCount) + ",," +
            "SanityCheck:," + nucleiMeanIntensity + "," + nucleiIntThreshold + "," + piBG + "," + medianNucleiSize + "\n";
        File.append(csvLine, csvPath);
// Append QC slice to the consolidated QC stack
        labelText = condLabel + " | " + baseName;
        appendToQCStack(qcTitle, labelText);

        closeAllButLogAndQCStack();
    }
}


// -----------------------------------------------------------------------------
// Append to QC stack
// -----------------------------------------------------------------------------
function appendToQCStack(qcTitle, labelText) {

    selectWindow(qcTitle);
    // Ensure ROI overlay is visible before flattening
    roiManager("Show All");
    // Ensure a true RGB canvas (prevents LUT/state effects on text color)
    run("RGB Color");

    // Draw label text in white
    setFont("SansSerif", 24, "bold");
    setColor(255, 255, 255);
    setForegroundColor(255, 255, 255);
    drawString(labelText, 10, 32);
run("Flatten");
    run("RGB Color");
    rename("QC_slice");

    if (isOpen(resultsStackTitle) == 0) {
        rename(resultsStackTitle);
        return;
    }

    run("Concatenate...",
        "title=[" + resultsStackTitle + "] " +
        "image1=[" + resultsStackTitle + "] " +
        "image2=[QC_slice] " +
        "image3=[-- None --]");

    closeAllButLogAndQCStack();
}


// -----------------------------------------------------------------------------
// StarDist call
// -----------------------------------------------------------------------------
function runStarDistOnNuclei(nucTitle) {

    if (isOpen("ROI Manager")) {
        selectWindow("ROI Manager");
        roiManager("Reset");
    } else {
        run("ROI Manager...");
        roiManager("Reset");
    }

    if (!isOpen(nucTitle)) exit("Nuclei image window not found: " + nucTitle);

    selectWindow(nucTitle);

    run("Command From Macro",
        "command=[de.csbdresden.stardist.StarDist2D], " +
        "args=['input':'" + nucTitle + "', " +
        "'modelChoice':'Versatile (fluorescent nuclei)', " +
        "'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', " +
        "'probThresh':'0.5', 'nmsThresh':'0.4', 'outputType':'Both', 'nTiles':'1', " +
        "'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', " +
        "'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
}


// -----------------------------------------------------------------------------
// Robust PI background estimate: mean PI in pixels outside nuclei ROIs
// -----------------------------------------------------------------------------
function estimateBackgroundOutsideROIs(imgTitle) {

    // Paint all current ROI Manager ROIs into a mask image, then measure mean intensity
    // in the inverse (background) selection on the ORIGINAL image.

    selectWindow(imgTitle);
    run("Duplicate...", "title=__bgmask_tmp");

    selectWindow("__bgmask_tmp");
    run("8-bit");
    run("Select All");
    setForegroundColor(0, 0, 0);
    run("Fill");
    run("Select None");

    nR = roiManager("count");
    if (nR > 0) {
        setForegroundColor(255, 255, 255);
        for (i = 0; i < nR; i++) {
            roiManager("select", i);
            run("Fill");
        }
    }

    run("Create Selection");
    run("Make Inverse");

    // Add background ROI temporarily
    roiManager("Add");
    bgIndex = roiManager("count") - 1;

    selectWindow(imgTitle);
    roiManager("select", bgIndex);
    getStatistics(bgArea, bgMean, bgMin, bgMax, bgStd);

    // Remove background ROI, keep nuclei ROIs
    roiManager("select", bgIndex);
    roiManager("Delete");

    // Clear lingering selection
    selectWindow(imgTitle);
    run("Select None");

    if (isOpen("__bgmask_tmp")) { selectWindow("__bgmask_tmp"); run("Close"); }

    // Fallback: if selection failed, measure whole-image mean
    if (bgArea <= 0) {
        selectWindow(imgTitle);
        run("Select All");
        getStatistics(a2, m2, mi2, ma2, s2);
        run("Select None");
        return m2;
    }

    return bgMean;
}

function estimatePIBackgroundOutsideROIs(piTitle) {
    return estimateBackgroundOutsideROIs(piTitle);
}




// -----------------------------------------------------------------------------
// Compute 75th percentile of tumor-eligible nucleus areas in the current image.
// Tumor-eligible = not a T cell (if T marker provided) and area >= minTumorArea.
// Returns 0 if not enough ROIs.
// -----------------------------------------------------------------------------
function computeTumorAreaMedian(nucTitle, tTitle, minTumorArea, maxNucleusArea, tBG, tCellFactor) {

    nAll = roiManager("count");
    if (nAll <= 0) return 0;

    areas = newArray();

    for (i = 0; i < nAll; i++) {

        selectWindow(nucTitle);
        roiManager("select", i);
        getStatistics(aN, mN, miN, maN, sN);

        // Apply max nucleus area only for statistics (do not delete ROIs)
        if (maxNucleusArea > 0 && aN > maxNucleusArea) continue;

        // Exclude T cells based on marker intensity (if provided)
        isT = false;
        if (tTitle != "" && tCellFactor > 0 && tBG > 0) {
            selectWindow(tTitle);
            roiManager("select", i);
            getStatistics(aT, mT, miT, maT, sT);
            if (mT >= (tBG * tCellFactor)) isT = true;
        }

        if (!isT && aN >= minTumorArea) {
            areas = appendToArray(aN, areas);
        }
    }

    if (areas.length < 4) return 0;

    areas = sortArrayAscending(areas);

    idx = floor(0.50 * (areas.length - 1));
    if (idx < 0) idx = 0;
    if (idx >= areas.length) idx = areas.length - 1;

    return areas[idx];
}


// -----------------------------------------------------------------------------
// Compute 75th percentile of per-nucleus mean Hoechst intensities among tumor-eligible ROIs.
// Tumor-eligible = not a T cell (if T marker provided) and area >= minTumorArea.
// Returns 0 if not enough ROIs.
// -----------------------------------------------------------------------------
function computeTumorMeanNucMedian(nucTitle, tTitle, minTumorArea, 0) {

    nAll = roiManager("count");
    if (nAll <= 0) return 0;

    means = newArray();

    for (i = 0; i < nAll; i++) {

        selectWindow(nucTitle);
        roiManager("select", i);
        getStatistics(aN, mN, miN, maN, sN);

        isT = false;
        if (tTitle != "") {
            selectWindow(tTitle);
            roiManager("select", i);
            getStatistics(aT, mT, miT, maT, sT);
            if (aN <= 0 && mT > 0) isT = true;
        }

        if (!isT && aN >= minTumorArea) {
            means = appendToArray(mN, means);
        }
    }

    if (means.length < 4) return 0;

    means = sortArrayAscending(means);

    idx = floor(0.50 * (means.length - 1));
    if (idx < 0) idx = 0;
    if (idx >= means.length) idx = means.length - 1;

    return means[idx];
}


// -----------------------------------------------------------------------------
// Sort numeric array ascending (insertion sort; avoids non-short-circuit && issues)
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Remove overly large ROIs from ROI Manager (optional max nucleus area filter)
// -----------------------------------------------------------------------------
function filterNucleiByMaxArea(maxArea) {
    total = roiManager("count");
    // Delete from the end to keep indices valid
    for (i = total - 1; i >= 0; i--) {
        roiManager("select", i);
        getStatistics(area, mean, min, max, std);
        if (area > maxArea) {
            roiManager("Delete");
        }
    }
}

// -----------------------------------------------------------------------------
// Median ROI area across ALL nuclei ROIs (includes fragments/small objects)
// -----------------------------------------------------------------------------
function computeAllNucleiAreaMedian(nucTitle, maxNucleusArea) {

    nAll = roiManager("count");
    if (nAll <= 0) return 0;

    areas = newArray();

    for (i = 0; i < nAll; i++) {
        selectWindow(nucTitle);
        roiManager("select", i);
        getStatistics(aN, mN, miN, maN, sN);
        if (maxNucleusArea > 0 && aN > maxNucleusArea) continue;
        areas = appendToArray(aN, areas);
    }

    if (areas.length < 1) return 0;

    areas = sortArrayAscending(areas);
    idx = floor(0.50 * (areas.length - 1));
    if (idx < 0) idx = 0;
    if (idx >= areas.length) idx = areas.length - 1;
    return areas[idx];
}


// -----------------------------------------------------------------------------
// Reference mean Hoechst intensity: mean of TOP 25% largest tumor-eligible nuclei by area
// Tumor-eligible = area >= minTumorArea and not classified as T cell.
// -----------------------------------------------------------------------------
function computeRefMeanNucleiIntensityTop25(nucTitle, tTitle, minTumorArea, maxNucleusArea, tBG, tCellFactor) {

    // Mean Hoechst intensity of the largest 25% of tumor-eligible nuclei (by area).
    // Tumor-eligible: area >= minTumorArea, (optional) area <= maxNucleusArea, and NOT T-marker positive.

    nAll = roiManager("count");
    if (nAll <= 0) return 0;

    // collect (area, intensity) pairs for eligible nuclei
    areas = newArray();
    means = newArray();

    for (i = 0; i < nAll; i++) {

        selectWindow(nucTitle);
        roiManager("select", i);
        getStatistics(aN, mN, miN, maN, sN);

        if (aN < minTumorArea) continue;
        if (maxNucleusArea > 0 && aN > maxNucleusArea) continue;

        // Exclude T cells based on marker intensity (if provided)
        isT = false;
        if (tTitle != "" && tCellFactor > 0 && tBG > 0) {
            selectWindow(tTitle);
            roiManager("select", i);
            getStatistics(aT, mT, miT, maT, sT);
            if (mT >= (tBG * tCellFactor)) isT = true;
        }
        if (isT) continue;

        areas = appendToArray(aN, areas);
        means = appendToArray(mN, means);
    }

    n = areas.length;
    if (n < 4) return 0;

    // sort indices by area
    idx = newArray(n);
    for (k = 0; k < n; k++) idx[k] = k;

    // insertion sort by areas[idx[*]] ascending (macro-safe)
    for (k = 1; k < n; k++) {
        key = idx[k];
        j = k - 1;
        while (j >= 0) {
            if (areas[idx[j]] > areas[key]) {
                idx[j+1] = idx[j];
                j--;
            } else {
                break;
            }
        }
        idx[j+1] = key;
    }

    // top 25% count (ceil)
    kTop = floor(n * 0.25);
    if (kTop * 4 < n) kTop = kTop + 1; // round up if needed
    if (kTop < 1) kTop = 1;
    if (kTop > n) kTop = n;

    sum = 0;
    for (p = n - kTop; p < n; p++) {
        ii = idx[p];
        sum += means[ii];
    }

    return sum / kTop;
}


function sortArrayAscending(arr) {
    out = arr;
    for (i = 1; i < out.length; i++) {
        key = out[i];
        j = i - 1;
        while (j >= 0) {
            if (out[j] > key) {
                out[j + 1] = out[j];
                j--;
            } else {
                break;
            }
        }
        out[j + 1] = key;
    }
    return out;
}


// -----------------------------------------------------------------------------
// Close everything except Log and QC stack
// -----------------------------------------------------------------------------
function closeAllButLogAndQCStack() {
    list = getList("image.titles");
    for (i = 0; i < list.length; i++) {
        title = list[i];
        if (title != "Log" && title != resultsStackTitle) {
            if (isOpen(title)) {
                selectWindow(title);
                run("Close");
            }
        }
    }
}


// -----------------------------------------------------------------------------
// Find split-channel window flexibly (handles naming quirks like "\well-...")
// -----------------------------------------------------------------------------
function findSplitWindowFlexible(chanIdx, baseName) {
    list = getList("image.titles");

    prefix1 = "C" + chanIdx + "-";
    for (i = 0; i < list.length; i++) {
        t = list[i];
        if (startsWith(t, prefix1) && indexOf(t, baseName) >= 0) return t;
    }

    // Alternative prefix variant (rare)
    prefix2 = "C" + chanIdx + " ";
    for (i = 0; i < list.length; i++) {
        t = list[i];
        if (startsWith(t, prefix2) && indexOf(t, baseName) >= 0) return t;
    }

    return "";
}


// -----------------------------------------------------------------------------
// Close remaining split-channel windows related to baseName
// -----------------------------------------------------------------------------
function closeSplitWindowsForBase(baseName, keep1, keep2, keep3) {
    list = getList("image.titles");
    for (i = 0; i < list.length; i++) {
        t = list[i];
        if (t == keep1 || t == keep2 || t == keep3 || t == resultsStackTitle || t == "Log") continue;

        if (startsWith(t, "C") && indexOf(t, baseName) >= 0) {
            if (isOpen(t)) { selectWindow(t); run("Close"); }
        }
    }
}


// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
function appendArray(a1, a2) {
    n1 = a1.length;
    n2 = a2.length;
    out = newArray(n1 + n2);
    for (i = 0; i < n1; i++) out[i] = a1[i];
    for (i = 0; i < n2; i++) out[n1 + i] = a2[i];
    return out;
}

function getFilesByExtension(fileList, extension) {
    extension = toLowerCase(extension);
    files = newArray();
    for (i = 0; i < fileList.length; i++) {
        name = fileList[i];
        nameLower = toLowerCase(name);
        if (endsWith(nameLower, extension)) {
            files = appendToArray(name, files);
        }
    }
    return files;
}

function appendToArray(value, array) {
    tempArray = newArray(array.length + 1);
    for (r = 0; r < array.length; r++) tempArray[r] = array[r];
    tempArray[array.length] = value;
    return tempArray;
}

function stripExtension(name) {
    dot = lastIndexOf(name, ".");
    if (dot < 0) return name;
    return substring(name, 0, dot);
}

function getFolderName(path) {
    p = path;
    if (endsWith(p, File.separator))
        p = substring(p, 0, lengthOf(p) - 1);
    idx = lastIndexOf(p, File.separator);
    if (idx < 0) return p;
    return substring(p, idx + 1);
}

function normalizeFolderName(name) {
    n = name;
    while (endsWith(n, "/") || endsWith(n, "\\") || endsWith(n, File.separator)) {
        n = substring(n, 0, lengthOf(n) - 1);
    }
    return n;
}
