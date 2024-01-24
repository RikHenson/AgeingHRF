The gzipped tarball "all_rois.tar.gz" contains gzipped CSVs for each of the 4 ROIs (lAC, bVC, lMC, rMC):

    <roi>_BOLD.csv.gz = BOLD timeseries (first temporal singular vector over voxels within ROI), in units of percent change relative to mean over all scans and all voxels in brain

    and also for each of the 4 models (FIR, CAN, NLF, HDM)

        <roi>_<model>_beta.csv.gz = estimates for each participant and parameter of model
        <roi>_<model>_fit.csv.gz  = fitted response for each participant and timepoint
        <roi>_<model>_RMSE.csv.gz = root-mean-sum of squared residuals for each participant

For the matlab scripts to work, you need to extract these files, eg with "tar xfz all_rois.tar.gz" in linux.

Note for the FIR model, the <roi>_FIR_fit.csv.gz file does not exist because it is identical to <roi>_FIR_beta.csv.gz.

There is a Matlab file "SPM_CC110037_Stim.mat" containing the 1st-level fMRI analysis in SPM12 from an example participant (containing FIR design matrix in SPM.xX.X, as well as other meta data, eg TR etc)

There is a binary grey-matter mask image dartel_GM_mask.nii.gz (obtained by thresholding >0.5 the CamCAN DARTEL GM template).

There are also four MNI nii mask images for the four ROIs (obtained by thresholding the Group SPM analyses, as described in paper)

There are also CSV files for latent vascular factors (LVFs) and MEG ERF energies for each participant for mediation analyses in paper

Finally there is a CSV for a revised canonical HRF and its partial derivatives in revised_canonical_3bf_across3ROIs.csv
