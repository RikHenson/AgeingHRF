The gzipped CSVs here contain, for each of the 4 ROIs (lAC, bVC, lMC, rMC):

    <roi>_BOLD.csv.gz = BOLD timeseries (first temporal singular vector over voxels within ROI), in units of percent change relative to mean over all scans and all voxels in brain

    and also for each of the 4 models (FIR, CAN, NLF, HDM)

        <roi>_<model>_beta.csv.gz = estimates for each participant and parameter of model
        <roi>_<model>_fit.csv.gz  = fitted response for each participant and timepoint
        <roi>_<model>_SSE.csv.gz  = sum of squared residuals for each participant

Note for the FIR model, the <roi>_FIR_fit.csv.gz file does not exist because it is identical to <roi>_FIR_beta.csv.gz.

There is also a Matlab file "SPM_CC110037_Stim.mat" containing the 1st-level fMRI analysis in SPM12 from an example participant (containing FIR design matrix in SPM.xX.X, as well as other meta data, eg TR etc)

There are also 4 NII binary mask images containing voxels that define lAC, bVC, lMC and rMC.
