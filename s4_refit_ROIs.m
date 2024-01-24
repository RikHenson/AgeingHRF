% Matlab script for creating ROI timeseries
%
%   Download/clone to "bas_dir" below:
%       1. "AgeingHRF" from https://github.com/RikHenson/AgeingHRF/
%
%   Link to local SPM12 directory, or install from 
%       https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ 
%
% rik.henson@mrc-cbu.cam.ac.uk, Jan 2023


%% Finish RMSE of timeseries, including for NLF, and add RMSE for FIR fit?

clear

bas_dir = '/imaging/rh01/CamCAN/700/HRF_SMT/Revision' % Change to your working directory

git_dir = fullfile(bas_dir,'AgeingHRF-main') % https://github.com/RikHenson/AgeingHRF

spm_dir = '/imaging/local/software/spm_cbu_svn/releases/spm12_latest/' % Your local installation of SPM12
addpath(spm_dir);

addpath(fullfile(git_dir,'Matlab_Utils')) % Some SPM12 updates needed for below (eg to properly turn-off orthogonalisation of basis functions)
spm('Defaults','fMRI')

out_dir = fullfile(bas_dir, 'outputs');

ana_dir = fullfile(out_dir,'SPM_FIR'); % From s1_fit_SPM_FIR.m

roi_dir = fullfile(git_dir,'ROI_data') % CSVs (or could use VOI*mat files)

participants = spm_load(fullfile(git_dir,'participants_include.csv'));
nparticipants = length(participants.CCID)

roi_names = {'lAC','bVC','lMC','rMC'};
nrois = length(roi_names);

glms = {'Stim-locked','Resp-locked'};
glm_type = [1 1 2 2]; % 1 = stim-locked, 2 = resp-locked

% Load an example SPM 1st level model (from first participant)
load(fullfile(roi_dir,sprintf('SPM_CC%d_Stim.mat',participants.CCID(1))),'SPM');

sf      = max(SPM.xX.X(:,1)); % scaling of regressor
nscans  = SPM.nscan;
nparam  = SPM.xBF.order;
pst     = round(1000*[0.5:1:nparam] * SPM.xBF.length / SPM.xBF.order); % in ms to avoid writing decimals to CSV header!
npst    = length(pst); % = nparam for FIR
FIRpst  = pst; % (need for NLF below)
TR      = SPM.xY.RT;

RMSE = nan(4,nparticipants,nrois);


%% Get fitted timeseries from the SPM FIR model (bins NOT orthogonalised)
for r = 1:nrois
    fn_fit = fullfile(roi_dir,sprintf('%s_FIR%d_fit.csv',roi_names{r},nparam));
    fp_fit = fopen(fn_fit,'w');
    for b = 1:(npst-1); fprintf(fp_fit,'t%d,',pst(b)); end; fprintf(fp_fit,sprintf('t%d\n',pst(npst)));

    fn_RMSE = fullfile(roi_dir,sprintf('%s_FIR%d_RMSE.csv',roi_names{r},nparam));
    fp_RMSE = fopen(fn_RMSE,'w');
    fprintf(fp_RMSE,'RMSE\n');
    
    % Load the prewhitened adjusted timeseries for each participant
    xY  = struct2array(spm_load(fullfile(roi_dir,sprintf('%s_BOLD.csv.gz', roi_names{r}))));
    
    glm_dir = glms{glm_type(r)};
    
    for s = 1:nparticipants 
        name = sprintf('CC%d',participants.CCID(s));
        swd = fullfile(ana_dir,glm_dir,name);
        
        % Load the SPM design matrix X
        SPM = load(fullfile(swd,'SPM.mat')); SPM = SPM.SPM;
                     
        % Estimate FIR Betas
        Betas = SPM.xX.pKX*xY(s,:)';
        Betas = Betas * sf; % Scale to percent signal change
        for b = 1:(nparam-1); fprintf(fp_fit,'%6.5f,',Betas(b)); end; fprintf(fp_fit,'%6.5f\n',Betas(nparam));
        
        % Calculate residual sum of squares (RMSE) (scaling data to match rescaling of Betas)
        e = xY(s,:)' - (SPM.xX.xKXs.X*Betas/sf);
        RMSE(1,s,r) = sqrt(e'*e/nscans);

        fprintf(fp_RMSE,'%6.5f\n',RMSE(1,s,r));
        fprintf('.')
    end
    fprintf('\n')
    
    fclose(fp_fit); fclose(fp_RMSE);
    
    gzip(fn_fit); gzip(fn_RMSE)
    delete(fn_fit); delete(fn_RMSE); 
end


%% Get parameter estimates for SPM Canonical Basis Set

orthflag = 0; % do not orthogonalise

% Load example SPM.mat to get PST (same for all participants)
load(fullfile(roi_dir,sprintf('SPM_CC%d_Stim',participants.CCID(1))))
SPM.xBF.order = 3;
SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
SPM.Sess(1).U(1).orth = orthflag; 
SPM = spm_fMRI_design(SPM, 0); % 0 = don't save
nparam = SPM.xBF.order;
pst = round(1000*[0:SPM.xBF.dt:(SPM.xBF.length-SPM.xBF.dt)]);
npst = length(pst);
sf   = full(max(SPM.Sess(1).U.u));

params = {};
params{1} = {'hrf'};
params{2} = {'hrf (with time derivative)'};
params{3} = {'hrf (with time and dispersion derivatives)'};
nparams = [1 2 3];

for r = 1:nrois
    glm_dir = glms{glm_type(r)};
    
    % Load the prewhitened adjusted timeseries for each participant
    xY  = struct2array(spm_load(fullfile(roi_dir,sprintf('%s_BOLD.csv.gz', roi_names{r}))));
    
    for p = 1:length(params)
        param  = params{p};
        nparam = nparams(p);
        
        fn_beta = fullfile(roi_dir,sprintf('%s_CAN%d_beta.csv',roi_names{r},nparam));
        fp_beta = fopen(fn_beta,'w');
        for b = 1:(nparam-1); fprintf(fp_beta,'CAN_beta%d,',b); end; fprintf(fp_beta,sprintf('CAN_beta%d\n',nparam));
        
        fn_fit = fullfile(roi_dir,sprintf('%s_CAN%d_fit.csv',roi_names{r},nparam));
        fp_fit = fopen(fn_fit,'w');
        for b = 1:(npst-1); fprintf(fp_fit,'t%d,',pst(b)); end; fprintf(fp_fit,sprintf('t%d\n',pst(npst)));
        
        fn_RMSE = fullfile(roi_dir,sprintf('%s_CAN%d_RMSE.csv',roi_names{r},nparam));
        fp_RMSE = fopen(fn_RMSE,'w');
        fprintf(fp_RMSE,'RMSE\n');
        
        for s = 1:nparticipants
            name = sprintf('CC%d',participants.CCID(s));
            swd = fullfile(ana_dir,glm_dir,name);
            
            % Load the SPM design matrix X
            SPM = load(fullfile(swd,'SPM.mat')); SPM = SPM.SPM;
            
            % redefine the basis set
            SPM.xBF.order = nparam;
            SPM.xBF.name = param{1};
            SPM.Sess(1).U(1).orth = orthflag;
            SPM = spm_fMRI_design(SPM, 0); % 0 = don't save
            
            % Estimate CAN Betas (after prewhitening model to match data)
            SPM.xX.xKXs.X = spm_filter(SPM.xX.K,SPM.xX.W*SPM.xX.X);
            SPM.xX.pKX = pinv(SPM.xX.xKXs.X);
            Betas = SPM.xX.pKX*xY(s,:)';
            for b = 1:(nparam-1); fprintf(fp_beta,'%6.5f,',Betas(b)); end; fprintf(fp_beta,'%6.5f\n',Betas(nparam));
            
            % Calculate fit
            Y = SPM.xBF.bf * Betas(1:nparam) * sf;
            for b = 1:(npst-1); fprintf(fp_fit,'%6.5f,',Y(b)); end; fprintf(fp_fit,'%6.5f\n',Y(npst));
            
            % Calculate residual sum of squares (RMSE) (scaling data to match rescaling of Betas)
            e = xY(s,:)' - (SPM.xX.xKXs.X*Betas);
            RMSE(2,s,r) = sqrt(e'*e/nscans);

            fprintf(fp_RMSE,'%6.5f\n',RMSE(2,s,r));
            fprintf('.');
        end
        fprintf('\n')
        
        fclose(fp_beta); fclose(fp_fit); fclose(fp_RMSE);
        
        gzip(fn_beta); gzip(fn_fit); gzip(fn_RMSE)
        delete(fn_beta); delete(fn_fit); delete(fn_RMSE);
    end
end


%% Get parameter estimates for NLF

param = {'lat_off','lat_scl','amp_off','amp_scl','R2'};
nparam = 4; % R2 not a parameter but fit metric
 
pst = FIRpst;
npst = length(pst);

% invariant arguments for NLF_HRF below
P = [];
P.pst = pst/1000;
P.t0 = 0;
P.trange = [0 32];
%P.doplot = 2;
P.doplot = 0;
P.stepsize = [1 0.1]; 
P.difflim = 1e-6;
P.AbsCor = 1; % necessary for rMC deactivations
P.GridShift = linspace(-5,5,20);
P.GridStrch = linspace(.3,2,20);

for r = 1:nrois
 
    % Load up FIR
    fs = fullfile(roi_dir,sprintf('%s_FIR32_fit.csv.gz',roi_names{r}));  
    Y = struct2array(spm_load(fs));
    
    % Create template (SVD works better than mean for rMC)
    %P.yt  = squeeze(mean(Y(:,:,r),1));
    [~,~,V] = spm_svd(squeeze(Y));
    P.yt  = full(V(:,1));
    [~,i] = max(abs(P.yt)); if sign(P.yt(i))<0; P.yt  = -P.yt; end % ensure template always positive (eg for rMC)

    fn_beta = fullfile(roi_dir,sprintf('%s_NLF%d_beta.csv',roi_names{r},nparam));
    fp_beta = fopen(fn_beta,'w');
    for b = 1:(nparam-1); fprintf(fp_beta,'%s,',param{b}); end; fprintf(fp_beta,sprintf('%s\n',param{nparam}));

    fn_fit = fullfile(roi_dir,sprintf('%s_NLF%d_fit.csv',roi_names{r},nparam));
    fp_fit = fopen(fn_fit,'w');
    for b = 1:(npst-1); fprintf(fp_fit,'t%d,',pst(b)); end; fprintf(fp_fit,sprintf('t%d\n',pst(npst)));

%     fn_RMSE = fullfile(roi_dir,sprintf('%s_NLF%d_RMSE.csv',roi_names{r},nparam)); % Note this over FIR bins
%     fp_RMSE = fopen(fn_RMSE,'w');
%     fprintf(fp_RMSE,'RMSE\n');
    
    for s = 1:nparticipants
        P.y = Y(s,:)';
        NLF = NLF_HRF(P);
        
        % NLF parameters
        for b = 1:(nparam-1); fprintf(fp_beta,'%6.5f,',NLF.(param{b})); end; fprintf(fp_beta,'%6.5f\n',NLF.(param{nparam})); 
        
        % NLF fit
        for b = 1:(npst-1); fprintf(fp_fit,'%6.5f,',NLF.fit(b)); end; fprintf(fp_fit,'%6.5f\n',NLF.fit(npst));       
        
%         % Calculate residual sum of squares (RMSE) (scaling data to match rescaling of Betas)
%         fprintf(fp_RMSE,'%6.5f\n',sqrt(NLF.SSE/npst));
        fprintf('.')
    end
    fprintf('\n')
    
    fclose(fp_beta); fclose(fp_fit); %fclose(fp_RMSE);
    
    gzip(fn_beta); gzip(fn_fit); %gzip(fn_RMSE)
    delete(fn_beta); delete(fn_fit); %delete(fn_RMSE); 
end


%% re-insert NLF4 fit into GLM to get timeseries residuals 

for r = 1:nrois
    % Get participant-specific HRFs (NLF fits)
    SS_HRF = struct2array(spm_load(fullfile(roi_dir,sprintf('%s_NLF4_fit.csv.gz',roi_names{r}))));
    glm_dir = glms{glm_type(r)};
    
    parfor s=1:nparticipants
        name = sprintf('CC%d',participants.CCID(s));
        swd = fullfile(ana_dir,glm_dir,name);
        cd(swd)
        % Load the SPM design matrix X
        SPM = rikload('SPM.mat'); SPM = SPM.SPM;
        
        xY = rikload(fullfile(swd,sprintf('VOI_%s_1.mat', roi_names{r})));
        y  = xY.Y;      
          
        SPM.xBF.order = 1;
        SPM.xBF.name = 'User';
        SPM.xBF.bf = resample(SS_HRF(s,:)',round(1/SPM.xBF.dt),1);
        SPM = spm_fMRI_design(SPM, 0); % 0 = don't save

        SPM.xX.xKXs.X = spm_filter(SPM.xX.K,SPM.xX.W*SPM.xX.X);
        SPM.xX.pKX = pinv(SPM.xX.xKXs.X);
        
        b = SPM.xX.pKX*y; 
        
        % Calculate residual
        e = y - (SPM.xX.xKXs.X*b); % but have to re-estimate upsampled BF
        RMSE(3,s,r) = sqrt(e'*e/nscans);
    end
end

% Write out NLF RMSE
for r = 1:nrois
    fn_RMSE = fullfile(roi_dir,sprintf('%s_NLF%d_RMSE.csv',roi_names{r},nparam)); % Note this over FIR bins
    fp_RMSE = fopen(fn_RMSE,'w');
    fprintf(fp_RMSE,'RMSE\n');
    
    for s = 1:nparticipants
        fprintf(fp_RMSE,'%6.5f\n',RMSE(3,s,r));
    end
    fclose(fp_RMSE);    
    gzip(fn_RMSE);
    delete(fn_RMSE); 
end


%% Get parameter estimates for HDM

% Need to have created GCM.mat files from s3_fit_hdm.m
hdm_dir = fullfile(out_dir,'HDM_fits')

params = {};
params{1} = {'efficacy','decay','transit'};
params{2} = {'efficacy','decay','transit','alpha'};
params{3} = {'efficacy','decay','transit','alpha','feedback'};
params{4} = {'efficacy','decay','transit','alpha','feedback','E0'}; % E0 doesn't change enough with age?

hdm_sf = []; % Scaling of 1st-order kernel lost, so scale to match max-min of FIR fit
sfK2K1 = []; % Relative scaling of 2nd-order kernel vs 1st-order (ie degree of nonlinearity)
for p = 1:length(params)
    param  = params{p};
    nparam = length(param);
    is_logscale = ones(1,nparam); is_logscale(1) = 0; % assumes "efficiency" first

    for r = 1:nrois
%        GCM = GCMs_PEB{r};

        load(fullfile(hdm_dir,sprintf('GCM_HDM%d_%s.mat',nparam,roi_names{r}))); % if want pre-PEB params

        % Scaling of HDM fit lost, so scale to max-min FIR
        Y = spm_load(fullfile(roi_dir,sprintf('%s_FIR32_fit.csv.gz',roi_names{r})));
        Y = struct2array(Y);
       
        % Used for re-scaling to match (FIR) data below
        hdm_sf(r,1) = mean(mean(Y));
        hdm_sf(r,2) = max(mean(Y,1)) - min(mean(Y,1));
        
        % Relative scaling of 2nd- to 1st-order kernels, as estimate of degree of nonlinearity
        sf = [];
        for s = 1:length(GCM)
            sf(s,1) = max(GCM{s}.K1) - min(GCM{s}.K1);
            sf(s,2) = max(GCM{s}.K2(:)) - min(GCM{s}.K2(:));
        end
        sfK2K1(p,r,:) = sf(:,2)./sf(:,1);
        
        pst = round(1000*[1:GCM{1}.M.N]*GCM{1}.M.dt);
        npst = length(pst);
        
        % Identify free parameters (same for all participants)
        d=[];
        for v = 1:nparam
            d(v,1) = GCM{1}.M.De.(param{v});
        end
        
        fn_beta = fullfile(roi_dir,sprintf('%s_HDM%d_beta.csv',roi_names{r},nparam));
        fp_beta = fopen(fn_beta,'w');
        for b = 1:(nparam-1); fprintf(fp_beta,'%s,',param{b}); end; fprintf(fp_beta,sprintf('%s\n',param{nparam}));
        
        fn_fit = fullfile(roi_dir,sprintf('%s_HDM%d_fit.csv',roi_names{r},nparam));
        fp_fit = fopen(fn_fit,'w');
        for b = 1:(npst-1); fprintf(fp_fit,'t%d,',pst(b)); end; fprintf(fp_fit,sprintf('t%d\n',pst(npst)));
        
        fn_RMSE = fullfile(roi_dir,sprintf('%s_HDM%d_RMSE.csv',roi_names{r},nparam));
        fp_RMSE = fopen(fn_RMSE,'w');
        fprintf(fp_RMSE,'RMSE\n');
        
        % Fits
        Y = cell2mat(cellfun(@(HDM)HDM.K1', GCM, 'UniformOutput', false));
        
        % Now scale to FIR
        %Y = Y - mean(mean(Y)) + hdm_sf(r,1); 
        Yr = max(mean(Y,1)) - min(mean(Y,1));
        Y = Y * hdm_sf(r,2) / Yr;
         
        for s = 1:nparticipants
            Ep = [];
            for v = 1:nparam
                if is_logscale(v)
                    Ep(v) = d(v).*exp(GCM{s}.Ep.(param{v}));
                else
                    Ep(v) = d(v).*GCM{s}.Ep.(param{v});
                end
            end
            
            % HDM parameters
            for b = 1:(nparam-1); fprintf(fp_beta,'%6.5f,',Ep(b)); end; fprintf(fp_beta,'%6.5f\n',Ep(nparam));
            
            % HDM Fit
            for b = 1:(npst-1); fprintf(fp_fit,'%6.5f,',Y(s,b)); end; fprintf(fp_fit,'%6.5f\n',Y(s,npst));
            
            R  = GCM{s}.R;
            RMSE(4,s,r) = sqrt(R'*R/nscans);
            
            fprintf(fp_RMSE,'%6.5f\n', RMSE(4,s,r));
            fprintf('.')
        end
        fprintf('\n')
        
        fclose(fp_beta); fclose(fp_fit); fclose(fp_RMSE);
        
        gzip(fn_beta); gzip(fn_fit); gzip(fn_RMSE)
        delete(fn_beta); delete(fn_fit); delete(fn_RMSE);
    end
end

drange(sfK2K1(1,:,:));
for r = 1:nrois
    drange(sfK2K1(1,r,:));
end

return
