# "AgeingHRF" from https://github.com/RikHenson/AgeingHRF/

Derived data and matlab scripts for paper on effects of Age on HRF

The raw fMRI images are available on request from this website: 

    https://camcan-archive.mrc-cbu.cam.ac.uk/dataaccess/

    Specifically, select the "preproc_fMRI_DARTEL" images that have 
    been realigned, slice-time corrected and DARTEL normalised 
    (see above website for more details). 

    Note these are large files, and include images from two other
    runs (Rest and Movie) that you can delete.

Some "derived" data available in folders of this repository:

    "Onsets" (BIDS events files)
    "Confounds" (realignment parameters and CSF/WM timeseries) 
    "ROI_data" (some mask images, example SPM.mat file, external 
                measurements like MEG, updated canonical HRF)

(some files in these folders need unzipping once downloaded - see 
README.txt within each folder).

You need Matlab installed (www.mathworks.com; the present scripts used 
MATLAB Version 9.9.0.1495850 (R2020b) Update 1), and then download 
SPM12 freeware (present scripts used version r7771) from: 

     https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ 

You'll also need the repository folder "Matlab_Utils" higher in the matlab
path (using "addpath" after starting SPM in Matlab) because it contains 
some minor but important changes to SPM12 functions (as well as some other
functions used in Matlab scripts below).

Finally, for HDM fitting, you'll need the "HDM-toolbox" from 

    https://github.com/pzeidman/HDM-toolbox

The file "participants.csv" in root folder has ID & Age of each participant.

The Matlab scripts from raw data to figures in paper are:

    s1_fit_SPM_FIR.m - fit individual participant SPM FIR models to raw 
        fMRI images, and estimate a group GLM

    s2_get_VOI.m - extract timeseries for each ROI 

    s3_fit_hdm.m - fit the HDM and apply PEB+BMR

    s4_refit_ROIs.m - create CSV files with parameters and fits from each 
        model in the paper

    s5_produce_figures.m - create TIF files for majority of figures in
        paper (excluding schematics)

For the mediation analyses in s5_produce_figures.m, you'll also need to 
download the Matlab mediation toolbox from 
https://github.com/canlab/MediationToolbox, which uses some utility 
functions that you can download from https://github.com/canlab/CanlabCore

The folder "ManualGraphics" contains schematic figures created by hand.

For users without Matlab, some additional derived gzipped CSV files with 
fMRI timeseries and model parameters for each of 4 ROIs described in the 
paper are in the repository folder "ROI_data".

rik.henson@mrc-cbu.cam.ac.uk                Revised Jan 2024

with thanks to Wiktor Olszowy for careful code checking and improvements,
and Peter Zeidman, Kamen Tsvetanov and Pranay Yadav (and Wiktor again) 
for conceptual help with the accompanying paper.

