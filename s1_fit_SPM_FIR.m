% Matlab script for running SPM 1st-level FIR models on CamCAN Phase II
% Sensorimotor task (SMT) fMRI data
%
%  Download/clone to "bas_dir" below:
%
%       1. "AgeingHRF" from https://github.com/RikHenson/AgeingHRF/
%
%       2. Request preprocessed data from this website: 
%               https://camcan-archive.mrc-cbu.cam.ac.uk/dataaccess/
%          Specifically select the "aamod_norm_write_dartel_00001" images that have 
%          been realigned, slice-time corrected, DARTEL normalised 
%          (see above website for more details). 
%          Note these are large files, and include images from two other
%          runs (Rest and Movie) that you can delete.
%
%   Link to local SPM12 directory, or install from 
%       https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ 
%
% rik.henson@mrc-cbu.cam.ac.uk, Jan 2023
% Revised Nov 2023

clear

bas_dir = '/imaging/rh01/CamCAN/700/HRF_SMT/Revision' % Change to your working directory

raw_dir = '/imaging/camcan/cc700/mri/pipeline/release004/data_fMRI/aamod_norm_write_dartel_00001/'; % Change to where you downloaded processed SMT images from CamCAN website above

git_dir = fullfile(bas_dir,'AgeingHRF-main') % https://github.com/RikHenson/AgeingHRF

spm_dir = '/imaging/local/software/spm_cbu_svn/releases/spm12_latest/' % Your local installation of SPM12
addpath(spm_dir);

addpath(fullfile(git_dir,'Matlab_Utils')) % Some SPM12 updates needed for below (eg to properly turn-off orthogonalisation of basis functions), and for batch_spm_anova.m function below
spm('Defaults','fMRI')

bwd = fullfile(bas_dir,'outputs'); % Where results will go
try mkdir(bwd); end

ana_dir = fullfile(bwd,'SPM_FIR'); 
try mkdir(ana_dir); end

participants = spm_load(fullfile(git_dir,'participants.csv'));
nparticipants = length(participants.CCID)

glms = {'Stim-locked','Resp-locked'};

for sr = 1:length(glms)
    owd = fullfile(ana_dir,glms{sr})
    try mkdir(owd); end; cd(owd)
    
    parfor s = 1:nparticipants
        ccid = sprintf('CC%d',participants.CCID(s)); disp(ccid)            
        cd(owd); try mkdir(ccid); end; cd(ccid)
        
        SPM = [];
        SPM.nscan = 261;
        SPM.xY.RT = 1.97;
        SPM.xBF.T = 32;
        SPM.xBF.T0 = 17; % assumes data slice-time corrected to middle slice in space
        SPM.xBF.name = 'Finite Impulse Response';
        SPM.xBF.UNITS = 'secs';
        SPM.xBF.length = 32; SPM.xBF.order = 32;
        SPM.xBF.Volterra = 1;
        
        SPM.xY.P = spm_select('ExtFPList',fullfile(raw_dir,ccid,'SMT'),'^swauf.*nii',Inf); % No smoothing ("s" here is just 1-vox smoothing for DARTEL artifacts)
        
        ons = spm_load(fullfile(git_dir,'Onsets',sprintf('sub-%s_epi_smt_onsets.csv',ccid)));
        
        if strcmp(glms{sr},'Stim-locked')
            f = find(strncmp(ons.trial_type,'AudVid',6) | strncmp(ons.trial_type,'VidOnly',7) | strncmp(ons.trial_type,'AudOnly',7));
            if length(f)~=128; error('!'); end
        elseif strcmp(glms{sr},'Resp-locked')   
            f = find(strncmp(ons.trial_type,'button',6)); % number varies
        else
            error('!')
        end        
        
        SPM.Sess(1).U(1).name{1} = 'AllEvents';
        SPM.Sess(1).U(1).ons = ons.onset(f);
        SPM.Sess(1).U(1).dur = 0;
        SPM.Sess(1).U(1).P.name = 'none';
        SPM.Sess(1).U(1).orth = 0;   
         
        cof = spm_load(fullfile(git_dir,'Confounds',sprintf('sub-%s_epi_smt_confounds.csv',ccid)));
       
        SPM.Sess(1).C.name = fields(cof);
        SPM.Sess(1).C.C = struct2array(cof);
        
        SPM.xX.K(1).HParam = 128;
        SPM.xVi.form = 'AR(0.2)';
        
        SPM.xGX.iGXcalc = 'none';
        
        SPM = spm_fmri_spm_ui(SPM);
        
        SPM = spm_spm(SPM);
        
        try SPM = rmfield(SPM,'xCon'); end;
        c = [eye(SPM.xBF.order) zeros(SPM.xBF.order,size(SPM.Sess(1).C.C,2)+1)]';
        SPM.xCon(1)   = spm_FcUtil('Set','Effects of Interest','F','c',c,SPM.xX.xKXs);        
        spm_contrasts(SPM);
    end
end


%% Group Results for FIR

nBF = 32;

mc_age = detrend(participants.Age,0);

for sr = 1:length(glms)
    owd = fullfile(ana_dir,glms{sr})
    cd(owd)
    
    G = [];
    G.outdir = 'Group';
    G.sub_effects = 0; % not with age effects below
    G.nsph_flag = 0; % too long to estimate
    
    G.mask = fullfile(git_dir,'ROI_data','dartel_GM_mask.nii');
    
    G.imgfiles{1} = {};
    for s = 1:nparticipants
        ccid = sprintf('CC%d',participants.CCID(s)); 
        
        P = spm_select('ExtFPList',fullfile(owd,ccid),'^beta.*nii');
        
        G.imgfiles{1}{end+1} = P(1:nBF,:);
     end

     G.user_regs{1} = [kron(eye(nBF),mc_age) kron(eye(nBF),detrend(mc_age.^2,1))];
      
     G.contrasts{1}.name = 'HRF';
     G.contrasts{1}.c = [eye(nBF) zeros(nBF,2*nBF)];
     G.contrasts{1}.type = 'F';
     
     G.contrasts{2}.name = 'Age (Lin)';
     G.contrasts{2}.c = [zeros(nBF) eye(nBF) zeros(nBF)];
     G.contrasts{2}.type = 'F';
     
     G.contrasts{3}.name = 'Age (Quad)';
     G.contrasts{3}.c = [zeros(nBF,2*nBF) eye(nBF)];
     G.contrasts{3}.type = 'F';
     
     G.contrasts{4}.name = 'Age (Lin+Quad)';
     G.contrasts{4}.c = [zeros(2*nBF,nBF) eye(2*nBF)];
     G.contrasts{4}.type = 'F';
     
     batch_spm_anova(G);    
end

return
