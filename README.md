# ImageCytox

ImageCytox is an ImageJ/Fiji macro for nuclei-based quantification of image-based 3D cytotoxicity assays.

The macro segments nuclei, classifies them based on size and marker intensities, generates quality-control (QC) overlay images and exports per-image summary statistics as a CSV file. It is designed for robust, reproducible quantification of cytotoxicity in 3D imaging assays.

## Features
- Automated nuclei segmentation using StarDist
- Classification of nuclei into survivors, dead nuclei, T cells and excluded objects based on size (fragments and merged ROIs)
- Optional background subtraction to improve segmentation and intensity measurements
- QC overlay images visualizing classification results
- CSV output with summary metrics and sanity-check parameters

## Requirements
- Fiji / ImageJ
- StarDist plugin (2D)
- Image data with a nuclear stain (e.g. Hoechst); optional additional channels for T cell markers

## Usage
1. Open Fiji.
2. Run the ImageCytox macro (`.ijm` file).
3. Select the input directory containing the images.
4. Adjust analysis parameters in the dialog if needed.
5. The macro processes all images, writes results to a CSV file and saves a RGB stack with overlays color-coded based on their classification for quality control of segmentation and classification.

A detailed user guide and example data are provided in this repository.

## Citation
If you use ImageCytox in your work, please cite: 
Weigelin, B. (Year). ImageCytox: An ImageJ macro for nuclei-based quantification of image-based 3D cytotoxicity assays. Zenodo. https://doi.org/XXXX

## License
See the LICENSE file for details.
