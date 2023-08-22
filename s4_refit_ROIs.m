% Matlab script for creating ROI timeseries
%
%   Download/clone to "bas_dir" below:
%       1. "AgeingHRF" from https://github.com/RikHenson/AgeingHRF/
%
%   Link to local SPM12 directory, or install from 
%       https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ 
%
% rik.henson@mrc-cbu.cam.ac.uk, Jan 2023

clear

bas_dir = '/imaging/rh01/CamCAN/700/HRF_SMT' % Change to your working directory

git_dir = fullfile(bas_dir,'AgeingHRF') % https://github.com/RikHenson/AgeingHRF

spm_dir = '/imaging/local/software/spm_cbu_svn/releases/spm12_latest/' % Your local installation of SPM12
addpath(spm_dir);

addpath(fullfile(git_dir,'Matlab_Utils')) % Some SPM12 updates needed for below (eg to properly turn-off orthogonalisation of basis functions)
spm('Defaults','fMRI')

out_dir = fullfile(bas_dir, 'outputs');

ana_dir = fullfile(out_dir,'SPM_FIR'); % From s1_fit_SPM_FIR.m

roi_dir = fullfile(git_dir,'ROI_data') % CSVs (or could use VOI*mat files)

participants = spm_load(fullfile(git_dir,'participants.csv'));
nparticipants = length(participants.CCID)

roi_names = {'lAC','bVC','lMC','rMC'};
nrois = length(roi_names);

glms = {'Stim-locked','Resp-locked'};
glm_type = [1 1 2 2]; % 1 = stim-locked, 2 = resp-locked

%% Load an example SPM 1st level model (from first participant)

load(fullfile(roi_dir,sprintf('SPM_CC%d_Stim.mat',participants.CCID(1))),'SPM');

sf      = max(SPM.xX.X(:,1)); % scaling of regressor
nscans  = SPM.nscan;
nparam  = SPM.xBF.order;
pst     = round(1000*[0.5:1:nparam] * SPM.xBF.length / SPM.xBF.order); % in ms to avoid writing decimals to CSV header!
npst    = length(pst); % = nparam for FIR
FIRpst  = pst; % (need for NLF below)
TR      = SPM.xY.RT;


%% Get fitted timeseries from the SPM FIR model (bins NOT orthogonalised)

for r = 1:nrois
    fn_fit = fullfile(roi_dir,sprintf('%s_FIR%d_fit.csv',roi_names{r},nparam));
    fp_beta = fopen(fn_fit,'w');
    for b = 1:(npst-1); fprintf(fp_beta,'t%d,',pst(b)); end; fprintf(fp_beta,sprintf('t%d\n',pst(npst)));

    fn_SSE = fullfile(roi_dir,sprintf('%s_FIR%d_SSE.csv',roi_names{r},nparam));
    fp_SSE = fopen(fn_SSE,'w');
    fprintf(fp_SSE,'SSE\n');
    
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
        for b = 1:(nparam-1); fprintf(fp_beta,'%6.5f,',Betas(b)); end; fprintf(fp_beta,'%6.5f\n',Betas(nparam));
        
        % Calculate residual sum of squares (SSE) (scaling data to match rescaling of Betas)
        e = xY(s,:)' - (SPM.xX.xKXs.X*Betas/sf);
        fprintf(fp_SSE,'%6.5f\n',e'*e);
        fprintf('.')
    end
    fprintf('\n')
    
    fclose(fp_beta); fclose(fp_SSE);
    
    gzip(fn_fit); gzip(fn_SSE)
    delete(fn_fit); delete(fn_SSE); 
end


%% Get parameter estimates for SPM Canonical Basis Set

orthflag = 0; % do not orthogonalise

% Load example SPM.mat to get PST (same for all participants)
load(fullfile(ana_dir,'Stim-locked',sprintf('CC%d',participants.CCID(1)),'SPM.mat')); % Just get one SPM.mat
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
        
        fn_SSE = fullfile(roi_dir,sprintf('%s_CAN%d_SSE.csv',roi_names{r},nparam));
        fp_SSE = fopen(fn_SSE,'w');
        fprintf(fp_SSE,'SSE\n');
        
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
            
            % Calculate residual sum of squares (SSE) (scaling data to match rescaling of Betas)
            e = xY(s,:)' - (SPM.xX.xKXs.X*Betas);
            fprintf(fp_SSE,'%6.5f\n',e'*e);
            fprintf('.')
        end
        fprintf('\n')
        
        fclose(fp_beta); fclose(fp_fit); fclose(fp_SSE);
        
        gzip(fn_beta); gzip(fn_fit); gzip(fn_SSE)
        delete(fn_beta); delete(fn_fit); delete(fn_SSE);
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

    fn_SSE = fullfile(roi_dir,sprintf('%s_NLF%d_SSE.csv',roi_names{r},nparam)); % Note this is not SSE over volumes but over FIR bins
    fp_SSE = fopen(fn_SSE,'w');
    fprintf(fp_SSE,'SSE\n');
    
    for s = 1:nparticipants
        P.y = Y(s,:)';
        NLF = NLF_HRF(P);
        
        % NLF parameters
        for b = 1:(nparam-1); fprintf(fp_beta,'%6.5f,',NLF.(param{b})); end; fprintf(fp_beta,'%6.5f\n',NLF.(param{nparam})); 
        
        % NLF fit
        for b = 1:(npst-1); fprintf(fp_fit,'%6.5f,',NLF.fit(b)); end; fprintf(fp_fit,'%6.5f\n',NLF.fit(npst));       
        
        % Calculate residual sum of squares (SSE) (scaling data to match rescaling of Betas)
        fprintf(fp_SSE,'%6.5f\n',NLF.SSE);
        fprintf('.')
    end
    fprintf('\n')
    
    fclose(fp_beta); fclose(fp_fit); fclose(fp_SSE);
    
    gzip(fn_beta); gzip(fn_fit); gzip(fn_SSE)
    delete(fn_beta); delete(fn_fit); delete(fn_SSE); 
end


%% Get parameter estimates for HDM

% Need to have created GCM.mat files from s3_fit_hdm.m
hdm_dir = fullfile(out_dir,'HDM_fits')

params = {};
params{1} = {'efficacy','decay','transit'};
params{2} = {'efficacy','decay','transit','alpha'};
params{3} = {'efficacy','decay','transit','alpha','feedback'};
params{4} = {'efficacy','decay','transit','alpha','feedback','E0'}; % E0 doesn't change enough with age?

for p = 1:length(params)
    param  = params{p};
    nparam = length(param);
    is_logscale = ones(1,nparam); is_logscale(1) = 0; % assumes "efficiency" first
    
%    load(fullfile(hdm_dir,sprintf('GCMs_PEB_%d.mat',nparam))); % if want post-PEB params
    for r = 1:nrois
%        GCM = GCMs_PEB{r};
        load(fullfile(hdm_dir,sprintf('GCM_HDM%d_%s.mat',nparam,roi_names{r}))); % if want pre-PEB params
        
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
        
        fn_SSE = fullfile(roi_dir,sprintf('%s_HDM%d_SSE.csv',roi_names{r},nparam));
        fp_SSE = fopen(fn_SSE,'w');
        fprintf(fp_SSE,'SSE\n');
        
        % Fits
        Y = cell2mat(cellfun(@(HDM)HDM.K1', GCM, 'UniformOutput', false));
        
        % Load up FIR just for scaling
        FIR = struct2array(spm_load(fullfile(roi_dir,sprintf('%s_FIR32_fit.csv.gz',roi_names{r}))));
        
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
            fprintf(fp_SSE,'%6.5f\n', R'*R);
            fprintf('.')
        end
        fprintf('\n')
        
        fclose(fp_beta); fclose(fp_fit); fclose(fp_SSE);
        
        gzip(fn_beta); gzip(fn_fit); gzip(fn_SSE)
        delete(fn_beta); delete(fn_fit); delete(fn_SSE);
    end
end

return


%% old code to re-insert NLF4 fit into GLM to get residuals and subject-specific HRF fit
% 
% Y = cell(nparticipants,1); 
% R = Y; 
% nparams = 1;
% parfor s=1:nparticipants
%     Y{s} = nan(nparams,nrois);
%     R{s} = nan(1,nrois);
% 
%     for r = 1:nrois
%         % Choose the correct GLM (time-locked to motor or sensory inputs)
%         if is_motor(r)
%             glm_dir = glm_motor_dir;
%         else
%             glm_dir = glm_sensory_dir;
%         end
%         
%         name = sprintf('CC%d',participants.CCID(s));
%         roi_dir = fullfile(fwd,glm_dir);
%         swd = fullfile(roi_dir,name);
%         cd(swd)
% 
%         % Load the SPM design matrix X
%         SPM = rikload('SPM.mat'); SPM = SPM.SPM;
%         
%         xY = rikload(fullfile(roi_dir,sprintf('VOI_%s_Sm1_%s_SS1_1.mat', roi_names{r}, name)));
%         y  = xY.Y;% 
%         
%         SPM.xBF.order = nparams;
%         SPM.xBF.name = 'User';
%         SPM.xBF.bf = resample(squeeze(Fits(s,:,r))',round(1/SPM.xBF.dt),1);
%         SPM = spm_fMRI_design(SPM, 0); % 0 = don't save
% 
%         SPM.xX.xKXs.X = spm_filter(SPM.xX.K,SPM.xX.W*SPM.xX.X);
%         SPM.xX.pKX = pinv(SPM.xX.xKXs.X);
%         
%         b = SPM.xX.pKX*y; 
%         
%         % Store betas
%         Y{s}(:,r) = b(1:nparams)';
%         
%         % Calculate residual
% %         e = y - (X*b); % ignores prewhitening of design matrix
%         e = y - (SPM.xX.xKXs.X*b); % but have to re-estimate upsampled BF
%         R{s}(r) = e'*e;
%     end
% end
% Y = permute(reshape(cat(1,Y{:}),[nparams nparticipants nrois]),[2 1 3]);
% R = cat(1,R{:});
% 
% mean(R)
% 
% % Calculate peri-stimulus time labels for FIR bins
% tY = (0.5:1:size(SPM.xBF.bf,1))*dt;
% 
% % Calculate times for observations
% dt = SPM.xY.RT(1);
% tR = (0:(SPM.nscan-1)) * dt + (SPM.xBF.T0/SPM.xBF.T);
% 
% save(fullfile(out_dir,'SSF_betas.mat'),'tY','Y','tR','R','roi_names');

% figure, hold on
% for r = 1:nrois
%     subplot(2,2,r),hold on
%     for s = 1:nparticipants
%         plot(tY,resample(squeeze(Fits(s,:,r))',round(1/SPM.xBF.dt),1) * squeeze(mean(Y(:,:,r)))');
%     end
% end
