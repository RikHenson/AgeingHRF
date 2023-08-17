% Matlab script for fitting HDM to ROI timeseries
%
%   Download/clone to "bas_dir" below:
%       1. "AgeingHRF" from https://github.com/RikHenson/AgeingHRF/
%       2. "HDM-toolbox" from https://github.com/pzeidman/HDM-toolbox
%   
%   Link to local SPM12 directory, or install from 
%       https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ 
%
%   You need SPM.mat and VOI.mat files from s1_fit_SPM_FIR.m and
%   s2_get_VOI.m
%
% rik.henson@mrc-cbu.cam.ac.uk, Jan 2023

clear

bas_dir = '/imaging/rh01/CamCAN/700/HRF_SMT' % Change to your working directory

git_dir = fullfile(bas_dir,'AgeingHRF') % https://github.com/RikHenson/AgeingHRF

spm_dir = '/imaging/local/software/spm_cbu_svn/releases/spm12_latest/' % Your local installation of SPM12
addpath(spm_dir);

addpath(fullfile(git_dir,'Matlab_Utils')) % Some SPM12 updates needed for below (eg to properly turn-off orthogonalisation of basis functions)
spm('Defaults','fMRI')

out_dir = fullfile(bas_dir,'outputs'); % Where results will go

ana_dir = fullfile(out_dir,'SPM_FIR'); 

hdm_dir = fullfile(out_dir,'HDM_fits');
try mkdir(hdm_dir); end

participants = spm_load(fullfile(git_dir,'participants.csv'));
nparticipants = length(participants.CCID)

roi_names = {'lAC','bVC','lMC','rMC'};
nrois = length(roi_names);

glms = {'Stim-locked','Resp-locked'};
glm_type = [1 1 2 2]; % 1 = stim-locked, 2 = resp-locked

addpath(fullfile(bas_dir,'HDM-toolbox-master','toolbox')) % https://github.com/pzeidman/HDM-toolbox

% Set binary vector with on / off for each SPM regressor (we'll use the
% first FIR bin as the driving input)
u_idx = 1; 
  
cd(hdm_dir)

params = {};
params{1} = {'efficacy','decay','transit'};
params{2} = {'efficacy','decay','transit','alpha'};
params{3} = {'efficacy','decay','transit','alpha','feedback'};
params{4} = {'efficacy','decay','transit','alpha','feedback','E0'}; % E0 doesn't change enough with age?

TR = 1.97;

for p = 1:length(params)
    param  = params{p};
    nparam = length(param);
    
    % HDM options
    options = struct();
    options.TE = 0.03;
    options.delays = TR/2; % TR/2 because middle slice (and not dealing with downsampled BFs here) (Empty means defaults to TR/2)
    %options.delays = TR/2 + TR/2; % additional TR/2 because priors based on this?
    eval(sprintf('options.priorfun = @spm_hdm_priors_hdm%d',nparam));
    
    for r = 1:nrois
        bwd = fullfile(ana_dir,glms{glm_type(r)});
        
        % Specify models
        GCM = cell(nparticipants,1);
        parfor s = 1:nparticipants
            % Subject's name
            name = ['CC' num2str(participants.CCID(s))];
            
            % Choose SPM file for timing
            SPM = fullfile(bwd, name, 'SPM.mat');
            
            % Choose the VOI file for data
            xY = fullfile(bwd, name, sprintf('VOI_%s_1.mat', roi_names{r}));
            
            % Specify HDM
            HDM = spm_hdm_specify(SPM,xY,u_idx,options);
            HDM.subject = name;
            
            % Estimate
            GCM{s} = spm_hdm_estimate(HDM);
        end
        
        % Save
        fn = fullfile(hdm_dir,sprintf('GCM_HDM%d_%s.mat',nparam,roi_names{r}));
        save(fn,'GCM');
    end
end


%% Run PEB for each ROI
M = struct();
M.Q = 'fields';
M.X = [ones(nparticipants,1) zscore(participants.Age)];

for p = 1:length(params)
    param  = params{p};
    nparam = length(param);
    
    PEBs = cell(1,nrois); BMAs = PEBs; GCMs_PEB = PEBs;
    for r = 1:nrois
        load(fullfile(hdm_dir,sprintf('GCM_HDM%d_%s.mat',nparam,roi_names{r})));
        
        [PEBs{r}, GCMs_PEB{r}] = spm_dcm_peb(GCM,M,param);
        
        BMAs{r} = spm_dcm_peb_bmc(PEBs{r});
    end
    save(fullfile(hdm_dir,sprintf('BMAs_%d.mat',nparam)),'BMAs','PEBs','roi_names','param','M');
    save(fullfile(hdm_dir,sprintf('GCMs_PEB_%d.mat',nparam)),'GCMs_PEB','roi_names','param','M');
end

% Compare F across models (cannot comapre across ROIs)
allF = [];
for p = 1:length(params)
    param  = params{p};
    nparam = length(param);
    load(fullfile(hdm_dir,sprintf('BMAs_%d.mat',nparam)),'BMAs','PEBs','roi_names','param','M');    
    allF(p,:) = [PEBs{1}.F PEBs{2}.F PEBs{3}.F PEBs{4}.F];
end
allF

