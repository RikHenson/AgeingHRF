% Matlab script for producing figures in paper
%
% To use:
%   Download/clone to "bas_dir" below:
%       1. "AgeingHRF" from https://github.com/RikHenson/AgeingHRF/
%
%   Link to local SPM12 directory, or install from 
%       https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ 
%
% rik.henson@mrc-cbu.cam.ac.uk, Jan 2023

clear

bas_dir = '/imaging/rh01/CamCAN/700/HRF_SMT/Revision' % Change to wherever you downloaded/cloned "AgeingHRF" and "HDM-toolbox"

git_dir = fullfile(bas_dir,'AgeingHRF-main')

spm_dir = '/imaging/local/software/spm_cbu_svn/releases/spm12_latest/' % Your local installation of SPM12
addpath(spm_dir);

% Some SPM12 updates needed for below (eg to properly turn-off orthogonalisation of basis functions, and for variance explained by K2 etc)
addpath(fullfile(git_dir,'Matlab_Utils')) 
spm('Defaults','fMRI')

out_dir = fullfile(bas_dir,'outputs'); % Where results will go

hdm_dir = fullfile(out_dir,'HDM_fits');

try mkdir(fullfile(out_dir,'Graphics')); end % for all figures

participants = spm_load(fullfile(git_dir,'participants_include.csv'));
nparticipants = length(participants.CCID)
age_range = [min(participants.Age) max(participants.Age)]

roi_names = {'lAC','bVC','lMC','rMC'};
nrois = length(roi_names);
roi_dir = fullfile(git_dir,'ROI_data');

models = {'FIR32','CAN3','NLF4','HDM3'};
nmods  = length(models);

addpath(fullfile(bas_dir,'HDM-toolbox-master','toolbox')) % https://github.com/pzeidman/HDM-toolbox

% Get FIR parameters from example SPM fMRI design matrix
load(fullfile(roi_dir,sprintf('SPM_CC%d_Stim',participants.CCID(1))))
nFIR = SPM.xBF.order;
tFIR = (0.5:1:nFIR) * SPM.xBF.length / SPM.xBF.order;
sFIR = max(SPM.xX.X(:,1)); % Scaling factor for percent signal change (height of FIR bin ~ 16)
nscan = SPM.nscan;
max_pst = 16; % for plotting below


%% Figure 1 - show basis sets for FIR, Can and NLF (in higher time resolution)
f1 = figure; 
subplot(3,1,1)
%strips(SPM.xBF.bf)

% FIR
dt  = 1/round(1/SPM.xBF.dt); % because spm_get_bf rounds to nearest bin
pst = [dt:dt:SPM.xBF.length];
plot([-dt pst pst(end)+dt]',[zeros(1,nFIR); SPM.xBF.bf; zeros(1,nFIR)]) % bin edges are not vertical because of plotting, but helps visualisation
axis([-1 33 -0.1 1.2])
grid on
set(gca,'YTick',[0],'XTick',[0:8:32],'XTickLabel',[0:8:32])

% Can
SPM.xBF.order = 3;
SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
SPM.Sess(1).U(1).orth = 0;
SPM = spm_fMRI_design(SPM, 0); % 0 = don't save
subplot(3,1,2)
bf = SPM.xBF.bf;
for b=1:size(bf,2); bf(:,b) = bf(:,b)/(max(bf(:,b))-min(bf(:,b))); end
%strips(bf)
dt = SPM.xBF.dt;
t  = [0:dt:(SPM.xBF.length-dt)]';
plot(t,bf)
axis([-1 33 -0.7 1])
grid on
set(gca,'YTick',[0],'XTick',[0:8:32],'XTickLabel',[0:8:32])

% FIR SVD (used in NLF)
fs = fullfile(roi_dir,sprintf('%s_FIR32_fit.csv.gz',roi_names{1}));
Y = struct2array(spm_load(fs));
[~,~,V] = spm_svd(squeeze(Y));
bf  = full(V(:,1));
subplot(3,1,3)
plot(tFIR,bf)
axis([-1 33 -0.2 0.7])
grid on
set(gca,'YTick',[0],'XTick',[0:8:32],'XTickLabel',[0:8:32])

eval(sprintf('print -dtiff -f%d %s',f1.Number,fullfile(out_dir,'Graphics','Fig1_BFs.tif')))

% Last graphic of HDM from Friston paper, organised in powerpoint and saved in "ManualGraphics" folder


%% Figure 2 - MIPs and FIRs from SPM Group Analyses

% Prepared manually by exporting parts from SPM windows, organising them within powerpoint, and saving in "ManualGraphics" folder


%% Figure 3 - HRFs heatmap for each model (rows) and ROI (columns)
s3 = figure('OuterPosition',[100 100 1100 1100]);  % Supplementary Figure 4
tiledlayout(nmods,nrois,'TileSpacing','Compact'); 
ysample = [1:100:nparticipants];
smooth_flag = 1; pst_tick = [0.5 5 10 15];
FIR_RMSE = nan(nparticipants,nmods,nrois);
minY = 0; maxY = 0;
for r = 1:nrois
    for m = 1:nmods
        nexttile((m-1)*nmods+r)

        Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{m})));
        pst = strvcat(fields(Y)); pst = str2num(pst(:,2:end))/1000;
        Y = struct2array(Y);
        ind = find(pst <= max_pst);
        Y = Y(:,ind);
       
        if m==1
            FIR_Y = Y;
            FIR_T = pst(ind);
            PST_ind = FIR_T;
        else            
            PST_ind = [];
            for t = 1:length(FIR_T) % Note that last bin will be same bin for HDM since max(pst)=24
                [~,PST_ind(end+1)] = min(abs(pst-FIR_T(t)));
            end
            
            FIR_RMSE(:,m,r) = sqrt(mean((Y(:,PST_ind) - FIR_Y).^2,2));
        end
 
        if smooth_flag % Smooth the plot across participants 
            for t = 1:size(Y,2)
                Y(:,t) = smooth(Y(:,t),5);
            end
        end
        
        imagesc(Y);
                
        minY = min([minY min(Y(:))]);
        maxY = max([maxY max(Y(:))]);

        % resample for x-axis ticks
        PST_ind = [];
        for t = 1:length(pst_tick) % Note that last bin will be same bin for HDM where max(pst)=24
            [~,PST_ind(end+1)] = min(abs(pst-pst_tick(t)));
        end

        if m == 1
            title(roi_names{r})
        end
        if m  == nmods
            set(gca,'XTick',PST_ind,'XTickLabel',pst_tick,'FontSize',10);
            xlabel('PST (s)','FontSize',12);
        else
            set(gca,'XTick',[]);
        end
        if r == 1
            set(gca,'YTick',ysample,'YTickLabel',participants.Age(ysample),'FontSize',10)
            ylabel(sprintf('%s: Age (years)',models{m}),'FontSize',12)
        else
            set(gca,'YTick',[]);
        end
        
        colormap('parula');
        colorbar        
        axis square;
        drawnow;
    end
end

for r = 1:nrois
    for m = 1:nmods
        nexttile((m-1)*nmods+r)
        caxis manual; caxis([minY maxY]); % ensure same scale for all ROIs
    end
end

eval(sprintf('print -dtiff -f%d %s',s3.Number,fullfile(out_dir,'Graphics','HRF_heatmaps.tif')))


%% Figure 4 - HRFs by age tertile for each model (rows) and ROI (columns)
p = prctile(participants.Age,[33 67]);
pp{1} = find(participants.Age <= p(1));
pp{2} = find(participants.Age > p(1) & participants.Age <= p(2));
pp{3} = find(participants.Age > p(2));
%cp = {'r-','g-','b-'}; % fits

%c = colormap('parula'); for p = 1:3; cp{p} = c((p-1)*127+1,:); end
for p = 1:3; cp{p} = ones(1,3)*(3-p)/3; end % grayscale

f4 = figure('OuterPosition',[100 100 1100 1100]); 
minY = zeros(1,nrois); maxY = zeros(1,nrois);
for m = 1:nmods
    for r = 1:nrois
        subplot(nmods,nrois,(m-1)*nmods+r) ; hold on
        yline(0)
        grid on, axis square
       
        % Load FIR ("pst data")
        Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{1})));
        pst = strvcat(fields(Y)); pst = str2num(pst(:,2:end))/1000;
        Y = struct2array(Y);
        ind = find(pst <= max_pst); 

        if strcmp(models{m}(1:3),'FIR') % plot FIR as model fit
            for p = 1:length(pp)
                h(p) = plot(pst(ind),mean(Y(pp{p},ind)),'Color',cp{p},'LineStyle','-','LineWidth',2);
            end
        else % plot FIR as data and store FIR for later scaling of HDM
            fy =[]; for p = 1:length(pp)
                fy(p,:) = mean(Y(pp{p},ind));              
                plot(pst(ind),fy(p,:),'Color',cp{p},'LineStyle',':','LineWidth',2);
            end
            
            Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{m})));
            pst = strvcat(fields(Y)); pst = str2num(pst(:,2:end))/1000;
            Y = struct2array(Y);
            ind = find(pst <= max_pst); 
            
            my = []; for p = 1:length(pp)
                my(p,:) = mean(Y(pp{p},ind));              
            end
             
            minY(r) = min([minY(r) min([fy(:); my(:)])]);
            maxY(r) = max([maxY(r) max([fy(:); my(:)])]);
           
            for p = 1:length(pp)
                plot(pst(ind),my(p,:),'Color',cp{p},'LineStyle','-','LineWidth',2);
            end
        end
        
        if m == 1
            title(roi_names{r})
            if r == 1
                l = legend([h(1) h(2) h(3)],{'Y','M','L'},'Location','NorthEast');
                %set(l,'Position',[0.44 0.77 0.1278 0.1405])
            end
        end
        if m == nmods
            xlabel('PST (s)','FontSize',12);
            set(gca,'XTick',[0:5:max_pst],'FontSize',10); 
        else
            set(gca,'XTick',[]);
        end
        
        if r == 1
            ylabel(sprintf('%s: %% BOLD',models{m}),'FontSize',12)
        end
    end
end

for r = 1:nrois
    for m = 1:nmods
        subplot(nmods,nrois,(m-1)*nmods+r), hold on
        set(gca,'YTick',[-2:0.2:2],'FontSize',10) % assume never more than +/- 2% change!
        axis([0 max_pst minY(r) maxY(r)]) % ensure same scale for all ROIs
    end
end

eval(sprintf('print -dtiff -f%d %s',f4.Number,fullfile(out_dir,'Graphics','HRF_fits_tertiles.tif')))


%% Figure 5 - peak amplitude by age
f5 = figure('OuterPosition',[100 100 1100 1100]);
minY = zeros(1,nmods); maxY = zeros(1,nmods);
for m = 1:nmods
    for r = 1:nrois
        subplot(nmods,nrois,(m-1)*nrois+r)
        yline(0); hold on
        grid on % axis square
    
        if strcmp(models{m}(1:3),'NLF') % Load parameters
            Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_beta.csv.gz',roi_names{r},models{m})));
%            peak_amplitude = Y.amp_scl; % positive  scaling
            peak_amplitude = abs(Y.amp_scl); % positive or negative scaling
            ytitle = 'Amp. Scaling';
            
        else % Calculate from fit
            Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{m})));
            pst = strvcat(fields(Y)); pst = str2num(pst(:,2:end))/1000;
            Y = struct2array(Y);
            ind = find(pst <= max_pst);
            Y = Y(:,ind);
            
%            peak_amplitude = max(Y,[],2);       % Peak (positive only)
            peak_amplitude = max(abs(Y),[],2);   % Peak (positive or negative)
            ytitle = 'Peak Amp. (%)';
        end
                        
        minY(m) = min([minY(m) min(peak_amplitude(:))]);
        maxY(m) = max([maxY(m) max(peak_amplitude(:))]);

        plot(participants.Age, peak_amplitude, '.', 'MarkerSize', 4);
        
        set(gca,'XTick',[18:10:88],'FontSize',10);
        
        [Rval,Pval]=corr(participants.Age, peak_amplitude, 'type','Spearman');
        legend(sprintf('R=%+3.2f, p=%3.2f',Rval,Pval),'FontSize',8,'Location','NorthEast');
        
        if m == 1
            title(roi_names{r})
        elseif m  == nmods           
            xlabel('Age (years)','FontSize',12);
        end
        if r == 1
            ylabel(sprintf('%s: %s',models{m},ytitle),'FontSize',12)
        end
        
        drawnow;
    end
end

for r = 1:nrois
    for m = 1:nmods
        subplot(nmods,nrois,(m-1)*nmods+r)
        if strcmp(models{m},'NLF4')
            axis([min(participants.Age) max(participants.Age) minY(m) maxY(m)]) % ensure same scale for all ROIs, but different for NLF
        else
            others = setdiff([1:nmods],find(strcmp(models,'NLF4')));
            axis([min(participants.Age) max(participants.Age) min(minY(others)) max(maxY(others))]) % ensure same scale for all ROIs and other models
        end
%        axis([min(participants.Age) max(participants.Age) min(minY) max(maxY)]) % ensure same scale for all ROIs
    end
end

eval(sprintf('print -dtiff -f%d %s',f5.Number,fullfile(out_dir,'Graphics','peak_amplitude.tif')))


%% Figure 6 - peak latency by age
f6 = figure('OuterPosition',[100 100 1100 1100]); 
minY = zeros(1,nmods); maxY = zeros(1,nmods);
for m = 1:nmods
    for r = 1:nrois
        subplot(nmods,nrois,(m-1)*nrois+r)
        yline(0); hold on
        grid on % axis square
    
        if strcmp(models{m}(1:3),'NLF') % Load parameters
            Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_beta.csv.gz',roi_names{r},models{m})));
%            peak_latency = Y.lat_scl;
            peak_latency = abs(Y.lat_scl);
            ytitle = 'Latency Scaling';
            
        else % Calculate from fit
            Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{m})));
            pst = strvcat(fields(Y)); pst = str2num(pst(:,2:end))/1000;
            Y = struct2array(Y);
            ind = find(pst <= max_pst);
            Y = Y(:,ind);
            
%            [peak_amplitude,ind] = max(Y,[],2);       % Peak (positive only)
             [peak_amplitude,ind] = max(abs(Y),[],2);   % Peak (positive or negative)
            peak_latency = pst(ind);
            ytitle = 'Peak Latency (s)';
        end
                        
        minY(m) = min([minY(m) min(peak_latency(:))]);
        maxY(m) = max([maxY(m) max(peak_latency(:))]);
        
        plot(participants.Age, peak_latency, '.', 'MarkerSize', 4);
        
        set(gca,'XTick',[18:10:88],'FontSize',10);
        
        [Rval,Pval]=corr(participants.Age, peak_latency, 'type','Spearman');
        legend(sprintf('R=%+3.2f, p=%3.2f',Rval,Pval),'FontSize',8,'Location','NorthEast');
        
        if m == 1
            title(roi_names{r})
        elseif m  == nmods           
            xlabel('Age (years)','FontSize',12);
        end
        if r == 1
            ylabel(sprintf('%s: %s',models{m},ytitle),'FontSize',12)
        end
        
        drawnow;
    end
end

for r = 1:nrois
    for m = 1:nmods
        subplot(nmods,nrois,(m-1)*nmods+r)
        if strcmp(models{m},'NLF4')
            axis([min(participants.Age) max(participants.Age) minY(m) maxY(m)]) % ensure same scale for all ROIs, but different for NLF
        else
            others = setdiff([1:nmods],find(strcmp(models,'NLF4')));
            axis([min(participants.Age) max(participants.Age) min(minY(others)) max(maxY(others))]) % ensure same scale for all ROIs and other models
        end
%        axis([min(participants.Age) max(participants.Age) min(minY) max(maxY)]) % ensure same scale for all ROIs
    end
end

eval(sprintf('print -dtiff -f%d %s',f6.Number,fullfile(out_dir,'Graphics','peak_latency.tif')))


%% Figure 7 - HDM3 parameters by age
m = 4; % Model 4 = HDM
params = {'efficacy','decay','transit'};
labels = {'efficacy (a.u.)','decay (Hz)','transit (Hz)'};
nparams = length(params);

f7 = figure('OuterPosition',[100 100 1100 1100]); 
minY = zeros(1,nparams); maxY = zeros(1,nparams);
for r = 1:nrois    
    Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_beta.csv.gz',roi_names{r},models{m})));
    Y = struct2array(Y);
     
    for p = 1:nparams
        subplot(nparams, nrois, (p-1)*nrois + r)
        grid on % axis square
        
        plot(participants.Age, Y(:,p), '.', 'MarkerSize', 4);
                                
        minY(p) = min([minY(p) min(Y(:,p))]);
        maxY(p) = max([maxY(p) max(Y(:,p))]);

        set(gca,'XTick',[18:10:88],'FontSize',10);
        
        [Rval,Pval]=corr(participants.Age, Y(:,p), 'type','Spearman');
        legend(sprintf('R=%+3.2f, p=%3.2f',Rval,Pval),'FontSize',8,'Location','NorthEast');
        
        if p == 1
            title(roi_names{r})
        elseif p == nparams
            xlabel('Age (years)','FontSize',12);
        end
        if r == 1
            ylabel(sprintf('HDM3: %s',labels{p}),'FontSize',12)
        end
        
        drawnow;
    end
end

for r = 1:nrois
    for p = 1:nparams
        subplot(nparams, nrois, (p-1)*nrois + r)
        axis([min(participants.Age) max(participants.Age) min(minY(p)) max(maxY(p))]) % ensure same scale for all ROIs and other models
     end
end

eval(sprintf('print -dtiff -f%d %s',f7.Number,fullfile(out_dir,'Graphics','HDM3_age.tif')))


% Mediation by vascular factors?
vasc = readtable(fullfile(roi_dir,'vascular_factors.csv'));

% Average HDM parameters across ROIs
mHDM = zeros(nparticipants,3); 
for r = 1:nrois  
    Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_beta.csv.gz',roi_names{r},models{m})));
    labels = fields(Y);
    Y = struct2array(Y);
    mHDM = mHDM + Y/nrois;
end

% Need mediation toolbox: https://github.com/canlab/MediationToolbox
addpath(genpath('/home/rh01/matlab/mediation_toolbox/'))

pval = [];
for l = 1:3
    LVF = sprintf('LVF%d',l);
    vasc_vals = vasc.(LVF);
    vasc_inds = find(~isnan(vasc_vals));
    vasc_ccid = strvcat(vasc.CCID(vasc_inds)); vasc_ccid = str2num(vasc_ccid(:,3:end));
    
    [both_ccid,hdm_both,vasc_both] = intersect(participants.CCID,vasc_ccid);
    length(both_ccid)
    vasc_vals = vasc_vals(vasc_inds(vasc_both));
    %figure,
    for p = 1:3
        %subplot(3,1,p)
        %plot(vasc_vals,mHDM(hdm_both,p),'o'); xlabel(LVF); ylabel(labels{p});
        %    [Rval,Pval]=corr(vasc_vals,mHDM(hdm_both,p),'type','Spearman')
        %    [Rval,Pval]=corr(participants.Age(hdm_both),vasc_vals,'type','Spearman')
        %[Rval,Pval]=partialcorr(vasc_vals,mHDM(hdm_both,p),participants.Age(hdm_both),'type','Spearman');
        %legend(sprintf('R=%+3.2f, p=%3.2f',Rval,Pval),'FontSize',8,'Location','NorthEast');
        [~,mstats] = mediation(zscore(participants.Age(hdm_both)), zscore(mHDM(hdm_both,p)), zscore(vasc_vals), 'stats', 'verbose', 'names', {'Age', params{p}, LVF});
        pval(l,p) = mstats.p(end);
    end
end
pval
[l,p]=find(pval<.05/9)

LVF = sprintf('LVF%d',l);
vasc_vals = vasc.(LVF);
vasc_inds = find(~isnan(vasc_vals));
vasc_ccid = strvcat(vasc.CCID(vasc_inds)); vasc_ccid = str2num(vasc_ccid(:,3:end));
[both_ccid,hdm_both,vasc_both] = intersect(participants.CCID,vasc_ccid);
vasc_vals = vasc_vals(vasc_inds(vasc_both));
mediation(zscore(participants.Age(hdm_both)), zscore(mHDM(hdm_both,p)), zscore(vasc_vals), 'stats', 'verbose', 'names', {'Age', params{p}, LVF});
       
% Mediation by MEG ERFs (better to do MEG power)

meg = readtable(fullfile(roi_dir,'MEG_energy.csv'));

meg_roi_names = meg.Properties.VariableNames(2:end);
pval = [];
for r = 1:nrois      
    Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_beta.csv.gz',roi_names{r},models{m})));
    labels = fields(Y);
    Y = struct2array(Y);
    
    meg_vals = meg.(meg_roi_names{r});
    meg_vals = sqrt(meg_vals);
    meg_inds = find(~isnan(meg_vals));
    meg_ccid = meg.CCID(meg_inds);
    
    [both_ccid,hdm_both,meg_both] = intersect(participants.CCID,meg_ccid);
    %length(both_ccid)
    meg_vals = meg_vals(meg_inds(meg_both));
 
    [Rval,Pval]=corr(meg_vals,participants.Age(hdm_both),'type','Spearman'); fprintf('Age-MEG: %s: R=%4.3f, p=%4.3f\n',meg_roi_names{r},Rval,Pval);
    
    %figure
    for p = 1:3
        %subplot(3,1,p)
        %plot(meg_vals,Y(hdm_both,p),'o'); xlabel(meg_roi_names{r}); ylabel(labels{p});
        %[Rval,Pval]=partialcorr(meg_vals,Y(hdm_both,p),participants.Age(hdm_both),'type','Spearman'); fprintf('MEG-HDM, partial Age: %s: R=%4.3f, p=%4.3f\n',meg_roi_names{r},Rval,Pval);
        %legend(sprintf('R=%+3.2f, p=%3.2f',Rval,Pval),'FontSize',8,'Location','NorthEast'); 
        [~,mstats] = mediation(zscore(participants.Age(hdm_both)), zscore(Y(hdm_both,p)), zscore(meg_vals), 'stats', 'verbose', 'names', {'Age', labels{p}, meg_roi_names{r}});
         pval(r,p) = mstats.p(end);
    end
end
pval


%% Figure 8 - Plot HDM parameters after PEB
params = {'efficacy','decay','transit'};
nparams = length(params);
reord = [3 1 2]; % reorder PEB.Pnames to match L->R order in figure

effs = 'Age'; e = 2;
load(fullfile(hdm_dir,sprintf('BMAs_%d.mat',nparams)));

f8 = figure('OuterPosition',[100 100 1100 600]);
ax = [];
for r = 1:nrois
    ax(r) = subplot(1,nrois,r);
    
    Ep = BMAs{r}.Ep;
    Vp = diag(BMAs{r}.Cp);
    
    Ep = spm_unvec(Ep,PEBs{1}.Ep);
    Vp = spm_unvec(Vp,PEBs{1}.Ep);
    
    Ep = full(Ep(reord,e));
    Vp = full(Vp(reord,e));
    
    % Disregard trivially small parameter values, which we define as those with a rate constant less than 0.001Hz. 
    % This cut-off corresponds to half-life greater than 1000*ln(2)=693s=11.5 minutes, which is too slow to be relevant.    
    Ep(find(abs(Ep)<0.001)) = 0; 
    
    spm_plot_ci(Ep,Vp)
    %    spm_plot_ci(Ep./sqrt(Vp),zeros(length(Ep),1)); ylabel('Zscore')
    
    set(gca,'XTickLabel',params,'FontSize',10,'XTickLabelRotation',90);
    title(roi_names{r});
end
linkaxes(ax,'xy');

eval(sprintf('print -dtiff -f%d %s',f8.Number,fullfile(out_dir,'Graphics',sprintf('HDM%d_Params_%s.tif',nparams,effs))))


%% Figure 9 - cross-validated error 
pnames = cell(1,nmods);
for n = 1:nFIR; pnames{1}{n} = sprintf('t%d',round(tFIR(n)*1000)); end
for n = 1:3; pnames{2}{n} = sprintf('CAN_beta%d',n); end
pnames{3} = {'lat_off','lat_scl','amp_off','amp_scl','R2'};
pnames{4} = {'efficacy','decay','transit'};

for m = 1:nmods
    nparams(m) = length(pnames{m});
end

f9s(1) = figure('OuterPosition',[100 100 1100 600]);
hp(1) = uipanel('Parent',f9s(1),'Position',[0 0 1 1]);
mina = []; maxa = [];
err = nan(nparticipants,nmods); pred_age = nan(nparticipants,nmods);
for m = 1:nmods      
    Y = [];
    for r = 1:nrois
        if strcmp(models{m}(1:3),'FIR')
            rY = spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{m})));
        else
            rY = spm_load(fullfile(roi_dir,sprintf('%s_%s_beta.csv.gz',roi_names{r},models{m})));
        end
        ind = find(ismember(fields(rY),pnames{m})); % exclude R2 from NLF
        rY = struct2array(rY); rY = rY(:,ind);      
        Y = [Y rY];
    end
    
    for s = 1:nparticipants
        sidx = [1:nparticipants]; sidx(s)=[];
        X = [Y(sidx,:) ones(nparticipants-1,1)];
        B = pinv(X)*participants.Age(sidx);
        pred_age(s,m) = [Y(s,:) 1]*B;
        err(s,m) = participants.Age(s) - pred_age(s,m);
    end
    
    subplot(1,nmods,m,'Parent',hp(1)),hold on,
    plot(participants.Age, pred_age(:,m), '.')
    title(models{m})
    
    set(gca,'XTick',[18:35:88],'XTickLabel',[18:35:88],'FontSize',10)
    xlabel(sprintf('%s: Actual Age',models{m}),'FontSize',12)
    if m == 1
        set(gca,'YTick',[18:35:88],'YTickLabel',[18:35:88],'FontSize',10)
        ylabel(sprintf('%s: Predicted Age',models{m}),'FontSize',12)
    else
        set(gca,'YTick',[]);
    end
    axis square
    mina = min([mina; pred_age(:)]); maxa = max([maxa; pred_age(:)]);
    %       axis([min(pred_age(:)) max(pred_age(:)) min(pred_age(:)) max(pred_age(:))]);
    
    R = corr(participants.Age, pred_age(:,m));
    text(10,105,sprintf('R^2=%3.2f',R^2))
end

err = abs(err);
median(err)
p = nan(nmods); 
for m1 = 1:nmods 
    for m2 = (m1+1):nmods
        p(m1,m2) = signrank(err(:,m1),err(:,m2)); 
    end
end
p
%p = p*length(find(~isnan(p)))

f9s(2) = figure('OuterPosition',[100 100 1100 600]);
hp(2) = uipanel('Parent',f9s(2),'Position',[0 0 1 1]);
subplot(1,1,1,'Parent',hp(2)), hold on
line(kron(ones(nparticipants,1),[1:size(err,2)])',abs(err'),'Color',[0.8 0.8 0.8])
boxplot(abs(err))
ylabel('Absolute Error (Years)')
set(gca,'XTickLabel',models)
ylim([0 max(err(:))])
line([0 5],[10 10],'Color',[0 0 0],'LineStyle',':')
for m = 1:nmods
    figure(f9s(1))
    subplot(1,nmods,m,'Parent',hp(1))
    axis([mina maxa mina maxa]);
    
    figure(f9s(2))
    if p(1,m)<.05, text(m-0.03,60,'*','FontSize',18), end
end

f9 = figure('OuterPosition',[100 100 1100 1100]);
npanels = numel(f9s);
hp_sub = nan(1,npanels); 
for idx = 1:npanels
    hp_sub(idx) = copyobj(hp(idx),f9);
    set(hp_sub(idx),'Position',[0,(2-idx)/npanels,1,1/npanels]);
end
whitebg(gcf)

% Cannot work out how to set background to white for this multipanel figure
% so did by hand in graphics editor!

eval(sprintf('print -dtiff -f%d %s',f9.Number,fullfile(out_dir,'Graphics','LOOCV.tif')))
close(f9s(1)), close(f9s(2))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary Figure 1 - effect of changing HDM parameters

hdm_dir = fullfile(out_dir,'HDM_fits');

% Load exemplar model
load(fullfile(hdm_dir,'GCM_HDM6_lAC')); 

M = GCM{1}.M; % M.De same for all subjects, so take first

% Get priors and default values
[pE,pC,D,is_logscale] = spm_hdm_priors_hdm6(1,M);
pE.efficacy = 1; % Needed by spm_bireduce below

params = {'efficacy','decay','transit','alpha','feedback','E0'};
labels = {'efficacy \beta (HDM3)','decay \kappa (HDM3)','transit 1/\tau_h (HDM3)','stiffness \alpha (HDM4)','feedback \gamma (HDM5)','O_2 extraction, E_0 (HDM6)'};
nparam = length(params);

% Priors (plus bit above/below) from HDM.M.De
scales{1} = [0.1 0.2 0.3];
scales{2} = [0.32 0.64 1.28];
scales{3} = [0.51 1.02 2.04];
scales{4} = [0.17 0.33 0.66];
scales{5} = [0.21 0.41 0.82];
scales{6} = [0.20 0.40 0.80];

% Reassign priors
for p = 1:nparam
    D.(params{p}) = scales{p}(2);
end
    
% Legends
for p = [1:3]
    for v=1:length(scales{p}); legends{p}{v} = sprintf('%3.2fHz',scales{p}(v)); end
end
for p = [4:5]
   for v=1:length(scales{p}); legends{p}{v} = sprintf('%3.2f',scales{p}(v)); end
end
for p = 6
    for v=1:length(scales{p}); legends{p}{v} = sprintf('%3.1f%%',100*scales{p}(v)); end
end

% Plot colours
colours = [0 0 0
           0.4 0.4 0.4
           0.7 0.7 0.7;
           ];
ax = [];   
sf1 = figure('OuterPosition',[100 100 1100 1100]);
for p = 1:nparam
    ax(p) = subplot(2,3,p);
    
     if length(scales{p}) == 2
        styles = {':','-'};
    else
        styles = {'-',':','-'};
    end
    
    for j = 1:length(scales{p})
                
        % Set parameter
        de = D;         
        de.(params{p}) = scales{p}(j);
                
        % re-generate kernels
        M.De = de;
        [M0,M1,L1,L2] = spm_bireduce(M,pE);
        dt = M.dt;
        N  = M.N;
        [K0,K1,K2] = spm_kernels(M0,M1,L1,L2,N,dt);
        [H0,H1]    = spm_kernels(M0,M1,M.N,M.dt);    
        
        % Plot BOLD
        t = (1:M.N)*M.dt;    
        plot(t,K1,'LineWidth',3,'LineStyle',styles{j},'Color',colours(j,:)); 
        xlabel('Time (secs)');
        if mod(p-1,3)==0
            ylabel('BOLD'); 
        end
        hold on;            
    end    
    set(gca,'FontSize',12);
    title(labels{p},'FontSize',14);
    legend(legends{p},'FontSize',14);
    linkaxes(ax);
end

eval(sprintf('print -dtiff -f%d %s',sf1.Number,fullfile(out_dir,'Graphics','HDM6_vary_priors.tif')))


%% Supplementary Figure 2 - posterior covariance of 6 parameters

% Needed to reorder columns in covariance matrix to match above
reord = [6 1 3 4 2 5]; 
rpnames = {'\beta';'\kappa';'1/\tau_h';'\alpha';'\gamma';'E_0'}; % shorter versions for axes

sf2 = figure('OuterPosition',[100 100 1100 400]);
ME = zeros(nrois,nparam); MP = ME;
for r = 1:nrois
    load(fullfile(hdm_dir,sprintf('GCM_HDM%d_%s.mat',nparam,roi_names{r})));
    pidx = spm_find_pC(GCM{1});
    
    subplot(1,nrois,r)
    MC = zeros(length(pidx)); 
    for s = 1:length(GCM)
        Cp = full(GCM{s}.Cp(pidx,pidx));
        Cp = Cp(reord,reord);
        CC = spm_cov2corr(Cp);
        MC = MC + CC/length(GCM);
    end
               
    %disp(MC)
    imagesc(MC); axis square; 
    set(gca,'XTick',[1:length(rpnames)],'XTickLabel',rpnames,'XTickLabelRotation',90,'YTick',[1:length(rpnames)],'YTickLabel',rpnames)
    colormap('parula');
    caxis([-0.5 1]); colorbar;
    title(roi_names{r});    
end

eval(sprintf('print -dtiff -f%d %s',sf2.Number,fullfile(out_dir,'Graphics','HDM6_covar.tif')))


%% Supplementary Figure 3 - compare FIR fits with SPM's canonical HRF
%colormap('parula')
sf3 = figure; hold on 
lw = [4 3 2 1 2];
%cl = {'r--','g-.','b:','b-','k-'};
cl = {'r-','g-','b-','b:','k:'};
for r = 1:nrois-1 % to ignore rMC
    FIR = struct2array(spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{1}))));    
%     [U,S,V] = spm_svd(squeeze(FIR));
%     FIRcan(:,r) = full(V(:,1)); 
%     [~,i] = max(abs(FIRcan(:,r))); if sign(FIRcan(i,r))<0; FIRcan(:,r) = -FIRcan(:,r); end
%     FIRcan(:,r) = FIRcan(:,r) - FIRcan(1,r);
%     FIRcan(:,r) = FIRcan(:,r)/max(FIRcan(:,r));
    FIRmean = mean(FIR,1)';
    FIRmean = FIRmean - FIRmean(1);
    FIRmean = FIRmean/max(FIRmean);
%     figure,plot([FIRcan FIRmean/max(FIRmean)]) % nearly identical except rMC!
%     p = plot(tFIR,FIRcan(:,r),cl{r},'LineWidth',lw(r));
    p = plot(tFIR,FIRmean,cl{r},'LineWidth',lw(r));
end
[spmcan,p] = spm_hrf(1); spmcan = spm_hrf(1,p,SPM.xBF.T);
spmcan = spmcan - spmcan(1);
spmcan = spmcan/max(spmcan);
p=plot(tFIR-0.5,spmcan(1:nFIR),cl{5},'LineWidth',lw(5)); % -0.5 because CanHRF is returned from 0s onwards
axis([0 24 -0.2 1])
line([0 24]',[0 0]','Color',[0 0 0],'LineStyle',':');
set(gca,'FontSize',12); set(gca,'YTick',[0])
xlabel('PST (s)'); ylabel('Normalised Amplitude')
legend({roi_names{1:(nrois-1)} 'CanHRF'})

eval(sprintf('print -dtiff -f%d %s',sf3.Number,fullfile(out_dir,'Graphics','mean_HRF.tif')))

%% Produce a better SPM Can HRF (from 2 gammas) and its derivatives?
% fP = []; err = []; flag = [];
% cd(bas_dir) % assumes where can_hrf_fit lives
% for r = 2:nrois %-1 % to ignore rMC
%     FIR = struct2array(spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{1}))));    
%     FIR = mean(FIR,1)';
%     [fP(r,:), err(r), flag(r)] = fminsearch(@(vP) can_hrf_fit(vP,FIR,tFIR),[6 16 1 1 6 0 max(FIR)]);
% end
% 
% mean_FIR = zeros(length(tFIR),1);
% rois = [1:3]; % ignore rMC
% for r = rois
%     FIR = struct2array(spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{1}))));    
%     FIR = mean(FIR,1)';
%     mean_FIR = mean_FIR + FIR/length(rois);
% end
% [fP(end+1,:), err(end+1), flag(end+1)] = fminsearch(@(vP) can_hrf_fit(vP,mean_FIR,tFIR),[6 16 1 1 6 0 max(mean_FIR)]);
% flag
% err
% fP
% 
% [old_can, old_P] = spm_hrf(0.1);
% old_P(1:6)              % 6    16     1     1     6     0 
% new_P = fP(end,1:6)     % 4.47 12.33  0.47  2.95  2.06 -0.01
%
% new_can = spm_hrf(0.1,new_P);
% pst_can = [0:0.1:(length(new_can)-1)*dt]';
% figure,plot(pst_can, new_can/max(new_can))
% hold on,plot(pst_can, old_can/max(old_can),'r')

% % Temporal derivative (same 1s shift as in SPM) 
% bf      = new_can;
% p       = new_P;
% dp      = 1;
% p(6)    = p(6) + dp; % one second shift
% bf(:,2) = (bf(:,1) - spm_hrf(0.1,p))/dp;
% 
% % Dispersion derivative (same dp=0.01 as in SPM, even though p halved)
% p       = new_P;
% dp      = 0.01;
% p(3)    = p(3) + dp; 
% bf(:,3) = (bf(:,1) - spm_hrf(0.1,p))/dp;
% 
% bf = spm_orth(bf);
% figure,plot(pst_can, bf)
% 
% fn_can = fullfile(roi_dir,sprintf('revised_canonical_3bf_across%dROIs.csv',length(rois)));
% fp_can = fopen(fn_can,'w');
% for b = 1:(length(pst_can)-1); fprintf(fp_can,'t%d,',round(1000*pst_can(b))); end; fprintf(fp_can,sprintf('t%d\n',round(1000*pst_can(end))));
% for b = 1:(length(pst_can)-1); fprintf(fp_can,'%6.5f,',bf(b,1)); end; fprintf(fp_can,sprintf('%6.5f\n',bf(end,1)));
% for b = 1:(length(pst_can)-1); fprintf(fp_can,'%6.5f,',bf(b,2)); end; fprintf(fp_can,sprintf('%6.5f\n',bf(end,2)));
% for b = 1:(length(pst_can)-1); fprintf(fp_can,'%6.5f,',bf(b,3)); end; fprintf(fp_can,sprintf('%6.5f\n',bf(end,3)));
% fclose(fp_can);


%% Supplementary Figure 4 - residuals 
sf4b = figure('OuterPosition',[100 100 1100 1100]);  % Supplementary Figure 4
RMSE = [];
labels = models;
labels{3} = 'NLF1'; % Only 1 df when substituted back into 1st-level GLM
for m = 1:nmods
    for r = 1:nrois
        subplot(nmods,nrois,(m-1)*nmods+r)

        Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_RMSE.csv.gz',roi_names{r},models{m})));

        RMSE(:,m,r) = Y;
        
        plot(participants.Age, Y, '.', 'MarkerSize', 4);
        
        set(gca,'XTick',[18:10:88],'FontSize',10);
        %set(gca,'YTick',yvals,'FontSize',10);
        
        [Rval,Pval]=corr(participants.Age,Y,'type','Spearman');
        legend(sprintf('R=%+3.2f, p=%3.2f',Rval,Pval),'FontSize',8,'Location','NorthEast');
        
        if m == 1
            title(roi_names{r})
        elseif m  == nmods           
            xlabel('Age (years)','FontSize',12);
        end
        if r == 1
            ylabel(sprintf('%s: RMSE',labels{m}),'FontSize',12)
        end
        
        drawnow;
    end
end

for m = 1:nmods
    for r = 1:nrois
       subplot(nmods,nrois,(m-1)*nmods+r)
       axis([18 88 0 max(RMSE(:))])
    end
end

eval(sprintf('print -dtiff -f%d %s',sf4b.Number,fullfile(out_dir,'Graphics','residuals_age.tif')))


% Needs RMSE from above and FIR_RMSE from Figure 4
sf4a = figure('OuterPosition',[100 100 1100 1100]);  % Supplementary Figure 4
for r = 1:nrois
    subplot(2,nrois,r)
    e = squeeze(RMSE(:,:,r));
    e = log(e);
    boxplot(e)
    title(roi_names{r})
    set(gca,'XTickLabel',labels)
    t_matrix(e,1);
    if r == 1; ylabel('Log RMSE across scans'); end
    axis([0.5 4.5 min(log(RMSE(:))) max(log(RMSE(:)))])
end

labels = models;
for r = 1:nrois
    subplot(2,nrois,r+nrois)
    e = squeeze(FIR_RMSE(:,:,r));
    e = log(e);
    boxplot(e)
    title(roi_names{r})
    set(gca,'XTickLabel',labels)
    t_matrix(e,1);
    if r == 1; ylabel('Log RMSE across FIR bins'); end
    axis([0.5 4.5 min(log(FIR_RMSE(:))) max(log(FIR_RMSE(:)))])
end

eval(sprintf('print -dtiff -f%d %s',sf4a.Number,fullfile(out_dir,'Graphics','residuals_rmse.tif')))


%% Supplementary Figure S5 - PEB HDM fits for main effect

params = {'efficacy','decay','transit'};
nparams = length(params);

reord = [3 1 2]; % reorder PEB.Pnames to match L->R order in figure

effs = 'Mean'; e = 1; 

load(fullfile(hdm_dir,sprintf('BMAs_%d.mat',nparams)));

sf5 = figure('OuterPosition',[100 100 1100 600]);
ax = [];
for r = 1:nrois
    ax(r) = subplot(1,4,r);
    
    Ep = BMAs{r}.Ep;
    Vp = diag(BMAs{r}.Cp);
    
    Ep = spm_unvec(Ep,PEBs{1}.Ep);
    Vp = spm_unvec(Vp,PEBs{1}.Ep);
    
    Ep = full(Ep(reord,e));
    Vp = full(Vp(reord,e));
    
    % Disregard trivially small parameter values, which we define as those with a rate constant less then 0.001Hz. 
    % This cut-off corresponds to half-life greater than 1000*ln(2)=693s=11.5 minutes, which is too slow to be relevant.    
    Ep(find(abs(Ep)<0.001)) = 0; 

    spm_plot_ci(Ep,Vp);
    %spm_plot_ci(Ep./sqrt(Vp),zeros(length(Ep),1)); ylabel('Zscore') % Z-scores since different scalings?
    
    set(gca,'XTickLabel',params,'FontSize',10,'XTickLabelRotation',90);
    title(roi_names{r});
end
linkaxes(ax,'xy');

eval(sprintf('print -dtiff -f%d %s',sf5.Number,fullfile(out_dir,'Graphics',sprintf('HDM%d_Params_%s.tif',nparam,effs))))



%% Supplementary Figure S6 - LOOCV for Can1-Can3
params = {};
params{1} = {'hrf'};
params{2} = {'hrf (with time derivative)'};
params{3} = {'hrf (with time and dispersion derivatives)'};
nparams = [1 2 3];

label = {};
for m = 1:length(params)
    label{m} = sprintf('Can%d',nparams(m));
end

sf6 = figure('OuterPosition',[100 100 1100 600]); hold on
mina = []; maxa = [];
err = nan(nparticipants,length(params)); pred_age = nan(nparticipants,length(params));
for m = 1:length(params)      
    Y = [];
    for r = 1:nrois
        rY = spm_load(fullfile(roi_dir,sprintf('%s_CAN%d_beta.csv.gz',roi_names{r},nparams(m))));      
        try rY = struct2array(rY); end
        Y = [Y rY];
    end
      
    for s = 1:nparticipants
        sidx = [1:nparticipants]; sidx(s)=[];
        X = [Y(sidx,:) ones(nparticipants-1,1)];
        B = pinv(X)*participants.Age(sidx);
        pred_age(s,m) = [Y(s,:) 1]*B;
        err(s,m) = participants.Age(s) - pred_age(s,m);
    end  
    
%      subplot(1,nrois,r),hold on    
%      err = abs(err);
%      line(kron(ones(nparticipants,1),[1:size(err,2)])',abs(err'),'Color',[0.5 0.5 0.5])
%      boxplot(err); % bar(mean(err))
%      set(gca,'XTickLabel',label,'FontSize',10);
%      
%      if r==1; ylabel('Abs Error (Years)');end
%      % Compare first model with all others
%      p = nan(1,size(err,2)); for m = 2:size(err,2); p(m) = signrank(err(:,1),err(:,m)); end % CARE - tests against first model
%      %p = nan(1,size(err,2)); for m = 2:size(err,2); [~,p(m)] = ttest([err(:,1) - err(:,m)]); end % difference quite Gaussian even if individuals not
%      f = find(p < .05/(length(params)-1)); 
%      MM(:,r) = median(err)';
%      for ff = f; text(ff-0.1,MM(ff,r),'*','FontSize',12,'Color',[0 1 0]); end
%      title(roi_names{r});
end

err = abs(err);
median(err)
p = nan(length(params)); 
for m1 = 1:length(params) 
    for m2 = (m1+1):length(params)
        p(m1,m2) = signrank(err(:,m1),err(:,m2)); 
    end
end
p
%p = p*length(find(~isnan(p)))

line(kron(ones(nparticipants,1),[1:size(err,2)])',abs(err'),'Color',[0.8 0.8 0.8])
boxplot(abs(err))
ylabel('Absolute Error (Years)')
set(gca,'XTickLabel',label)
ylim([0 max(err(:))])
line([0 5],[10 10],'Color',[0 0 0],'LineStyle',':')
axis([mina maxa mina maxa]);
 
for m = 1:length(params)
    if p(1,m)<.05, text(m-0.02,52,'*','FontSize',18), end
end

eval(sprintf('print -dtiff -f%d %s',sf6.Number,fullfile(out_dir,'Graphics','Can_LOOCV.tif')))


%% Supplementary Figure S7 - LOOCV for HDM3-5
params = {};
params{1} = {'efficacy','decay','transit'};
params{2} = {'efficacy','decay','transit','alpha'};
params{3} = {'efficacy','decay','transit','alpha','feedback'};
params{4} = {'efficacy','decay','transit','feedback','alpha','E0'};

label = {};
for m = 1:length(params)
    label{m} = sprintf('HDM%d',length(params{m}));
end

sf7 = figure('OuterPosition',[100 100 1100 600]); hold on
mina = []; maxa = [];
err = nan(nparticipants,nmods); pred_age = nan(nparticipants,nmods);
for m = 1:length(params)      
    Y = [];
    for r = 1:nrois
        rY = spm_load(fullfile(roi_dir,sprintf('%s_HDM%d_beta.csv.gz',roi_names{r},length(params{m}))));      
        rY = struct2array(rY);     
        Y = [Y rY];
    end
      
    for s = 1:nparticipants
        sidx = [1:nparticipants]; sidx(s)=[];
        X = [Y(sidx,:) ones(nparticipants-1,1)];
        B = pinv(X)*participants.Age(sidx);
        pred_age(s,m) = [Y(s,:) 1]*B;
        err(s,m) = participants.Age(s) - pred_age(s,m);
    end  
    
%      subplot(1,nrois,r),hold on    
%      err = abs(err);
%      line(kron(ones(nparticipants,1),[1:size(err,2)])',abs(err'),'Color',[0.5 0.5 0.5])
%      boxplot(err); % bar(mean(err))
%      set(gca,'XTickLabel',label,'FontSize',10);
%      
%      if r==1; ylabel('Abs Error (Years)');end
%      % Compare first model with all others
%      p = nan(1,size(err,2)); for m = 2:size(err,2); p(m) = signrank(err(:,1),err(:,m)); end % CARE - tests against first model
%      %p = nan(1,size(err,2)); for m = 2:size(err,2); [~,p(m)] = ttest([err(:,1) - err(:,m)]); end % difference quite Gaussian even if individuals not
%      f = find(p < .05/(length(params)-1)); 
%      MM(:,r) = median(err)';
%      for ff = f; text(ff-0.1,MM(ff,r),'*','FontSize',12,'Color',[0 1 0]); end
%      title(roi_names{r});
end

err = abs(err);
median(err)
p = nan(nmods); 
for m1 = 1:nmods 
    for m2 = (m1+1):nmods
        p(m1,m2) = signrank(err(:,m1),err(:,m2)); 
    end
end
p
%p = p*length(find(~isnan(p)))

line(kron(ones(nparticipants,1),[1:size(err,2)])',abs(err'),'Color',[0.8 0.8 0.8])
boxplot(abs(err))
ylabel('Absolute Error (Years)')
set(gca,'XTickLabel',label)
ylim([0 max(err(:))])
line([0 5],[10 10],'Color',[0 0 0],'LineStyle',':')
axis([mina maxa mina maxa]);
 
for m = 1:nmods
    if p(1,m)<.05, text(m-0.03,60,'*','FontSize',18), end
end

eval(sprintf('print -dtiff -f%d %s',sf7.Number,fullfile(out_dir,'Graphics','HDM_LOOCV.tif')))



%% Supplementary Figure S8 - Age effects on second-order Volterra kernel, ie nonlinearities as function of poststimulus time?
sf8 = figure('OuterPosition',[100 100 1100 600]);

nparam = 3;
hdm_dir = fullfile(out_dir,'HDM_fits');
zAge = zscore(participants.Age);
mpst = []; T = []; p = []; tp = 1:8:64;
for r = 1:nrois
    load(fullfile(hdm_dir,sprintf('GCM_HDM%d_%s.mat',nparam,roi_names{r})));
    pst = round(1000*[1:GCM{1}.M.N]*GCM{1}.M.dt)/1000;
    
    for pl = 1:2
        Y = [];
        for s = 1:nparticipants
            if pl==1
                Y(s,:,:) = GCM{s}.K2;
            else
                Y(s,:,:) = zAge(s) * GCM{s}.K2;
            end               
        end
        
        mY = squeeze(mean(Y,1));
        %T = mY ./ (squeeze(std(Y)) / sqrt(nparticipants));
        %figure,imagesc(T),colorbar 
        subplot(2,4,(pl-1)*4 + r),hold on, colormap('gray')
        imagesc(mY); axis([1 64 1 64]), axis square
        set(gca,'XTick',tp,'XTickLabel',round(pst(tp)))
        set(gca,'YTick',tp,'YTickLabel',round(pst(tp)))
        title(roi_names{r});
        
        % Statistics on saturation effect (though no multiple comparison correction)
        if pl == 1
            [~,ind] = min(mY(:)); [x,y] = ind2sub(size(mY),ind);
            mpst(r,pl,:) = [x y];
            d = squeeze(Y(:,x,y));
            T(r,pl) = mean(d)/(std(d)/sqrt(size(Y,1)));
            p(r,pl) = 2*tcdf(-abs(T(r,pl)),size(Y,1)-1);
        elseif pl == 2
            [~,ind] = max(mY(:)); [x,y] = ind2sub(size(mY),ind);
            mpst(r,pl,:) = [x y];
            % d = squeeze(Y(:,x,y));
            d = squeeze(Y(:,mpst(r,1,1),mpst(r,1,2))); % Take from maximal mean effect
            T(r,pl) = mean(d)/(std(d)/sqrt(size(Y,1)));
            p(r,pl) = 2*tcdf(-abs(T(r,pl)),size(Y,1)-1);
        end
   end
end

eval(sprintf('print -dtiff -f%d %s',sf8.Number,fullfile(out_dir,'Graphics','Volterra2.tif')))

T, p
mpst = pst(mpst);
squeeze(mpst(:,1,:))
squeeze(mpst(:,2,:))

% Check variance of timeseries explained by K1 and K2 (and how much of that
% variance attributed to age)

K2var = nan(nparticipants,nrois); K2age = nan(1,nrois);
for r = 1:nrois
    load(fullfile(hdm_dir,sprintf('GCM_HDM%d_%s.mat',nparam,roi_names{r})));
    for s = 1:size(GCM,1)
        [~,K2var(s,r)] = kernel_variance(GCM{s});
    end
    [~,~,~,~,stats] = regress(K2var(:,r),[zAge zAge.^2 ones(nparticipants,1)]); 
    K2age(r) = stats(1);
end
median(K2var)
K2age






