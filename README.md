# "AgeingHRF" from https://github.com/RikHenson/AgeingHRF/

Derived data and matlab scripts for paper on effects of Age on HRF

The raw fMRI images are available on request from this website: 

    https://camcan-archive.mrc-cbu.cam.ac.uk/dataaccess/

    Specifically, select the "preproc_fMRI_Smooth" images that have 
    been realigned, slice-time corrected, DARTEL normalised and 
    smoothed (see above website for more details). 

    Note these are large files, and include images from two other
    runs (Rest and Movie) that you can delete.

Some "derived" data available in folders of this repository (which need
unzipping when uploaded):

    "Onsets" (BIDS events files)
    "Confounds" (realignment parameters and CSF/WM timeseries) 

You need Matlab and then to download the SPM12 freeware from: 

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

The folder "Graphics" contains schematic figures created by hand.


For users without Matlab, some additional derived gzipped CSV files with 
fMRI timeseries and model parameters for each of 4 ROIs described in the 
paper are in the repository folder "ROI_data".

rik.henson@mrc-cbu.cam.ac.uk                Jan 2023
