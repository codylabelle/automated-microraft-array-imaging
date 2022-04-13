# automated-microraft-array-imaging

**Description and Intellectual Property Notice**
This repository contains the MATLAB code (2019a) to run an automated microraft array platform for imaging, image analysis, and single cell isolation. The code was written in the lab of Nancy Allbritton by Cody LaBelle as part of research at the University of North Carolina at Chapel Hill. All code within this repository is the intellectual property of Cody LaBelle, Nancy Allbritton, and the University of North Carolina at Chapel Hill.

**GUIs**:
The GUIs created for this platform were developed using MATLAB's App Designer (mlapp files) which cannot be viewed on GITHUB, however are downloadable along with all other files located in this repository.

The primary GUI used to run the automated microraft array platform is named InvertedMVXGUI. All other GUIs are iniated within the main GUI and begin with GUI in their name.

**Image Processing**:
The primary image processing function is named ANALYZE_IMG_t0, which is used to analyze individual images in both brightfield and fluorescence channels. Image data is extracted from each image using this function.
