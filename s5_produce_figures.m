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

bas_dir = '/imaging/rh01/CamCAN/700/HRF_SMT' % Change to wherever you downloaded/cloned "AgeingHRF" and "HDM-toolbox"

git_dir = fullfile(bas_dir,'AgeingHRF')

spm_dir = '/imaging/local/software/spm_cbu_svn/releases/spm12_latest/' % Your local installation of SPM12
addpath(spm_dir);

addpath(fullfile(git_dir,'Matlab_Utils')) % Some SPM12 updates needed for below (eg to properly turn-off orthogonalisation of basis functions)
spm('Defaults','fMRI')

out_dir = fullfile(bas_dir,'outputs'); % Where results will go

hdm_dir = fullfile(out_dir,'HDM_fits');

try mkdir(fullfile(out_dir,'Graphics')); end % for all figures

participants = spm_load(fullfile(git_dir,'participants.csv'));
nparticipants = length(participants.CCID)
age_range = [min(participants.Age) max(participants.Age)]

roi_names = {'lAC','bVC','lMC','rMC'};
nrois = length(roi_names);
roi_dir = fullfile(git_dir,'ROI_data');

models = {'FIR32','CAN3','NLF4','HDM3'};
nmods  = length(models);


%% Get FIR parameters from example SPM fMRI design matrix
load(fullfile(roi_dir,sprintf('SPM_CC%d_Stim',participants.CCID(1))))
nFIR = SPM.xBF.order;
tFIR = (0.5:1:nFIR) * SPM.xBF.length / SPM.xBF.order;
sFIR = max(SPM.xX.X(:,1)); % Scaling factor for percent signal change (height of FIR bin ~ 16)
nscan = SPM.nscan;

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

% Graphic from Friston paper added manually in powerpoint to bottom of Fig1 to show HDM


%% Figure 2 - MIPs and FIRs from SPM Group Analyses

% Graphic prepared manually from saving component graphics from SPM


%% Figure 3 - HRFs heatmap for each model (rows) and ROI (columns)
s3 = figure('OuterPosition',[100 100 1100 1100]);  % Supplementary Figure 4
ysample = [1:100:nparticipants];
tsample = [0:5:15];
max_pst = 16; smooth_flag = 1;
for m = 1:nmods
    for r = 1:nrois
        subplot(nmods,nrois,(m-1)*nrois+r)

        Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{m})));
        pst = strvcat(fields(Y)); pst = str2num(pst(:,2:end))/1000;
        Y = struct2array(Y);
        ind = find(pst <= max_pst);
        xsample = [];
        for t = 1:length(tsample)
            [~,xsample(t)] = min(abs(pst(ind) - tsample(t)));
        end
        
%         minY = min(min(Y(:,ind)));
%         maxY = max(max(Y(:,ind)));

        if smooth_flag % Smooth the plot across participants (requires Stats Toolbox)
            for t = 1:size(Y,2)
                Y(:,t) = smooth(Y(:,t),5);
            end
        end
        
        imagesc(Y(:,ind));
        
        set(gca,'XTick',xsample,'XTickLabel',tsample,'FontSize',10);
        set(gca,'YTick',ysample,'YTickLabel',participants.Age(ysample),'FontSize',10)
        
        if m == 1
            title(roi_names{r})
        elseif m  == nmods           
            xlabel('PST (s)','FontSize',12);
        end
        if r == 1
            ylabel(sprintf('%s: Age (years)',models{m}),'FontSize',12)
        end
        
%         caxis manual; caxis([minY maxY]); % If want same scale for all ROIs
        colormap('parula');
        colorbar        
        axis square;
        drawnow;
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
for m = 1:nmods
    for r = 1:nrois
        subplot(nmods,nrois,(m-1)*nrois+r)
        yline(0); hold on
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
            for p = 1:length(pp)
                fy{p} = mean(Y(pp{p},ind));
                plot(pst(ind),fy{p},'Color',cp{p},'LineStyle',':','LineWidth',2);
            end
            
            Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{m})));
            pst = strvcat(fields(Y)); pst = str2num(pst(:,2:end))/1000;
            Y = struct2array(Y);
            ind = find(pst <= max_pst);
            
            for p = 1:length(pp)
                yy = mean(Y(pp{p},ind));               
                if strcmp(models{m}(1:3),'HDM') % Scaling lost in HDM kernels
                    yy = yy * (max(fy{p})-min(fy{p})) / (max(yy)-min(yy));
                end
                plot(pst(ind),yy,'Color',cp{p},'LineStyle','-','LineWidth',2);
            end
        end

        % axis([0 max_pst -0.5 1.6]) % hard coding y-axis
        % set(gca,'YTick',[-0.5:0.5:1.5],'FontSize',10)
        set(gca,'XTick',[0:5:max_pst],'FontSize',10); 
        
        if m == 1
            title(roi_names{r})
            if r == 1
                l = legend([h(1) h(2) h(3)],{'Y','M','L'},'Location','NorthEast');
                %set(l,'Position',[0.44 0.77 0.1278 0.1405])
            end
        elseif m == nmods
            xlabel('PST (s)','FontSize',12);
        end
        
        if r == 1
            ylabel(sprintf('%s: %% BOLD',models{m}),'FontSize',12)
        end
        drawnow;
    end
end

eval(sprintf('print -dtiff -f%d %s',f4.Number,fullfile(out_dir,'Graphics','HRF_fits_tertiles.tif')))


%% Figure 5 - peak amplitude by age
f5 = figure('OuterPosition',[100 100 1100 1100]); 
for m = 1:nmods
    for r = 1:nrois
        subplot(nmods,nrois,(m-1)*nrois+r)
        yline(0); hold on
        grid on % axis square
    
        if strcmp(models{m}(1:3),'NLF') % Load parameters
            Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_beta.csv.gz',roi_names{r},models{m})));
            peak_amplitude = Y.amplitude';
            ytitle = 'Amp. Scaling';
        else % Calculate from fit
            Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{m})));
            pst = strvcat(fields(Y)); pst = str2num(pst(:,2:end))/1000;
            Y = struct2array(Y);
            ind = find(pst <= max_pst);
            Y = Y(:,ind);
            
%            [peak_amplitude,ind] = max(Y');       % Peak (positive only)
            [peak_amplitude,ind] = max(abs(Y'));   % Peak (positive or negative)
            ytitle = 'Peak Amp. (%)';
        end
        
        plot(participants.Age, peak_amplitude, '.', 'MarkerSize', 4);
        
        set(gca,'XTick',[18:10:88],'FontSize',10);
        
        [Rval,Pval]=corr(participants.Age, peak_amplitude', 'type','Spearman');
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

eval(sprintf('print -dtiff -f%d %s',f5.Number,fullfile(out_dir,'Graphics','peak_amplitude.tif')))


%% Figure 6 - peak latency by age
f6 = figure('OuterPosition',[100 100 1100 1100]); 
for m = 1:nmods
    for r = 1:nrois
        subplot(nmods,nrois,(m-1)*nrois+r)
        yline(0); hold on
        grid on % axis square
    
        if strcmp(models{m}(1:3),'NLF') % Load parameters
            Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_beta.csv.gz',roi_names{r},models{m})));
            peak_latency = Y.stretch;
            ytitle = 'Latency Scaling';
        else % Calculate from fit
            Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{m})));
            pst = strvcat(fields(Y)); pst = str2num(pst(:,2:end))/1000;
            Y = struct2array(Y);
            ind = find(pst <= max_pst);
            Y = Y(:,ind);
            
%            [peak_amplitude,ind] = max(Y');       % Peak (positive only)
            [peak_amplitude,ind] = max(abs(Y'));   % Peak (positive or negative)
            peak_latency = pst(ind);
            ytitle = 'Peak Latency (s)';
        end
        
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

eval(sprintf('print -dtiff -f%d %s',f6.Number,fullfile(out_dir,'Graphics','peak_latency.tif')))


%% Figure 7 - HDM3 parameters by age
m = 4; % Model 4 = HDM
params = {'efficacy','decay','transit'};
nparams = length(params);

f7 = figure('OuterPosition',[100 100 1100 1100]); 

for r = 1:nrois
    
    Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_beta.csv.gz',roi_names{r},models{m})));
    labels = strvcat(fields(Y));
    Y = struct2array(Y);
     
    for p = 1:nparams
        subplot(nparams, nrois, (p-1)*nrois + r)
        grid on % axis square
        
        plot(participants.Age, Y(:,p), '.', 'MarkerSize', 4);
        
        set(gca,'XTick',[18:10:88],'FontSize',10);
        
        [Rval,Pval]=corr(participants.Age, Y(:,p), 'type','Spearman');
        legend(sprintf('R=%+3.2f, p=%3.2f',Rval,Pval),'FontSize',8,'Location','NorthEast');
        
        if p == 1
            title(roi_names{r})
        elseif p == nparams
            xlabel('Age (years)','FontSize',12);
        end
        if r == 1
            ylabel(sprintf('HDM3: %s',labels(p,:)),'FontSize',12)
        end
        
        drawnow;
    end
end

eval(sprintf('print -dtiff -f%d %s',f7.Number,fullfile(out_dir,'Graphics','HDM3_age.tif')))


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
    
    Ep(find(abs(Ep)<0.001)) = 0; % remove parameters that BMR pruned
    %spm_plot_ci(Ep,Vp);
    spm_plot_ci(Ep./sqrt(Vp),zeros(length(Ep),1)); ylabel('Zscore') % Z-scores since different scalings?
    
    set(gca,'XTickLabel',params,'FontSize',10,'XTickLabelRotation',90);
    title(roi_names{r});
end
linkaxes(ax,'xy');
eval(sprintf('print -dtiff -f%d %s',f8.Number,fullfile(out_dir,'Graphics',sprintf('HDM%d_Params_%s.tif',nparams,effs))))


%% Figure 9 - cross-validated error 
pnames = cell(1,nmods);
for n = 1:nFIR; pnames{1}{n} = sprintf('t%d',round(tFIR(n)*1000)); end
for n = 1:3; pnames{2}{n} = sprintf('CAN_beta%d',n); end
pnames{3} = {'Lat. Off.','Lat. Scl.','Amp. Off.','Amp. Scl.'};
pnames{4} = {'efficacy','decay','transit'};

for m = 1:nmods
    nparams(m) = length(pnames{m});
end

f9 = figure('OuterPosition',[100 100 1100 600]);

MM = []; %AgeR2 = [];
for r = 1:nrois
     err = nan(nparticipants,nmods);   
     for m = 1:nmods    
        if strcmp(models{m}(1:3),'FIR')
            Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{m})));
        else
            Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_beta.csv.gz',roi_names{r},models{m})));
        end
        ind = find(ismember(fields(Y),pnames{m})); % exclude R2 from NLF
        Y = struct2array(Y); Y = Y(:,ind);
        
%        [~,~,~,~,AgeR2(m,r)] = glm(participants.Age,[Y ones(nparticipants,1)],[eye(nparams(m)) zeros(nparams(m),1)]',0);
         
        for s = 1:nparticipants
            sidx = [1:nparticipants]; sidx(s)=[];
            X = [Y(sidx,:) ones(nparticipants-1,1)];
            B = pinv(X)*participants.Age(sidx);
            pred_age = [Y(s,:) 1]*B;
            err(s,m) = participants.Age(s) - pred_age;
        end
     end
     subplot(1,nrois,r),hold on,
     
     err = abs(err);
     line(kron(ones(nparticipants,1),[1:size(err,2)])',abs(err'),'Color',[0.5 0.5 0.5])
     boxplot(err); % bar(mean(err))
     set(gca,'XTickLabel',models,'FontSize',10);
     
     if r==1; ylabel('Abs Error (Years)');end
     %[~,~,p] = t_matrix(AllErr,1);
     % Compare FIR with all others
     p = nan(1,size(err,2)); for m = 2:size(err,2); p(m) = signrank(err(:,1),err(:,m)); end
     %p = nan(1,size(err,2)); for m = 2:size(err,2); [~,p(m)] = ttest([err(:,1) - err(:,m)]); end % difference quite Gaussian even if individuals not
     f = find(p < .05/(nmods-1)); 
     MM(:,r) = median(err)';
     for ff = f; text(ff-0.1,MM(ff,r),'*','FontSize',12,'Color',[0 1 0]); end
     title(roi_names{r});
end
MM
mean(MM,2)
mean(MM(:))

eval(sprintf('print -dtiff -f%d %s',f9.Number,fullfile(out_dir,'Graphics','LOOCV.tif')))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary Figure 1 - effect of changing HDM parameters

hdm_dir = fullfile(out_dir,'HDM_fits');

% Load exemplar model
load(fullfile(hdm_dir,'GCM_HDM6_lAC')); HDM = GCM{1}; clear GCM

M = HDM.M;

% Get priors and default values
[pE,pC,D,is_logscale] = spm_hdm_priors_hdm6(1,M);

params = {'efficacy','decay','transit','feedback','alpha','E0'};
labels = {'efficacy \beta','decay \kappa','transit 1/\tau_h','feedback \gamma','stiffness \alpha','O_2 extraction, E0'};

% Parameter values to set in HDM.M.De
scales{1} = [0.1 0.2 0.3];
scales{2} = [0.32 0.64 1.28];
scales{3} = [0.51 1.02 2.04];
scales{4} = [0.21 0.41 0.82];
scales{5} = [0.17 0.33 0.66];
scales{6} = [0.20 0.40 0.80];

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
for i = 1:length(params)
    ax(i) = subplot(2,3,i);
    
     if length(scales{i}) == 2
        styles = {':','-'};
    else
        styles = {'-',':','-'};
    end
    
    for j = 1:length(scales{i})
                
        de = D;
        P = pE;
        
        % Set parameter
        de.(params{i}) = scales{i}(j);
        
        % Set efficacy      
        P.efficacy  = 1;
        if ~strcmp(params{i},'efficacy')
            de.efficacy = 0.2;
        end
        
        % re-generate kernels
        M.De = de;
        [M0,M1,L1,L2] = spm_bireduce(M,P);
        dt = M.dt;
        N  = M.N;
        [K0,K1,K2] = spm_kernels(M0,M1,L1,L2,N,dt);
        [H0,H1]    = spm_kernels(M0,M1,M.N,M.dt);    
        
        % Plot BOLD
        t = (1:M.N)*M.dt;    
        plot(t,K1,'LineWidth',3,'LineStyle',styles{j},'Color',colours(j,:)); 
        xlabel('Time (secs)');
        if mod(i-1,3)==0
            ylabel('BOLD'); 
        end
        hold on;            
    end    
    set(gca,'FontSize',12);
    title(labels{i},'FontSize',14);
    legend(legends{i},'FontSize',14);
    linkaxes(ax);
end

eval(sprintf('print -dtiff -f%d %s',sf1.Number,fullfile(out_dir,'Graphics','HDM6_vary_priors.tif')))


%% Supplementary Figure 2 - posterior covariance of 6 parameters

% Select 6 param model
nparam = 6
reord = [6 1 2 3 4 5]; rpnames = {'\beta';'\kappa';'\gamma';'1/\tau';'\alpha';'E_0'};

sf2 = figure('OuterPosition',[100 100 1100 400]);
for r = 1:nrois
    load(fullfile(hdm_dir,sprintf('GCM_HDM%d_%s.mat',nparam,roi_names{r})));
    pidx = spm_find_pC(GCM{1});
    
    subplot(1,nrois,r)
    MC = zeros(length(pidx));
    for s = 1:length(GCM)
        CC = spm_cov2corr(full(GCM{s}.Cp));
        CC = CC(pidx,pidx); CC = CC(reord,reord);
        MC = MC + CC/length(GCM);
        %MC = MC + abs(CC(pidx,pidx))/length(GCM);
    end
    %disp(MC)
    imagesc(MC); axis square; 
    set(gca,'XTick',[1:length(rpnames)],'XTickLabel',rpnames,'YTick',[1:length(rpnames)],'YTickLabel',rpnames)
    colormap('parula');
    caxis([-0.5 1]); colorbar;
    title(roi_names{r});    
end

eval(sprintf('print -dtiff -f%d %s',sf2.Number,fullfile(out_dir,'Graphics','HDM6_covar.tif')))


%% Supplementary Figure 3 - compare FIR fits with SPM's canonical HRF
%colormap('parula')
sf3 = figure; hold on 
lw = [2 2 2 2 1];
cl = {'r-','g:','b-.','b--','k-'};
for r = 1:nrois
    FIR = struct2array(spm_load(fullfile(roi_dir,sprintf('%s_%s_fit.csv.gz',roi_names{r},models{1}))));    
    [U,S,V] = spm_svd(squeeze(FIR));
    FIRcan(:,r) = full(V(:,1)); 
    [~,i] = max(abs(FIRcan(:,r))); if sign(FIRcan(i,r))<0; FIRcan(:,r) = -FIRcan(:,r); end
    FIRcan(:,r) = FIRcan(:,r) - FIRcan(1,r);
    FIRcan(:,r) = FIRcan(:,r)/max(FIRcan(:,r));
    FIRmean = mean(FIR,1)';
%    figure,plot([FIRcan FIRmean/max(FIRmean)]) % nearly identical except rMC!
    p = plot(tFIR,FIRcan(:,r),cl{r},'LineWidth',lw(r));
%    set(p,'Color',ones(1,3)*(4-r)/4)
end
spmcan = spm_hrf(1);
%spmcan = spm_hrf(0.5); spmcan=spmcan(2:2:end);
spmcan = spmcan/max(spmcan);
p=plot(tFIR-0.5,spmcan(1:nFIR),cl{5},'LineWidth',lw(5)); % -0.5 because CanHRF is returned from 0s onwards
%set(p,'Color',ones(1,3)/4)
axis([0 24 -0.2 1])
set(gca,'FontSize',12); set(gca,'YTick',[])
xlabel('PST (s)'); ylabel('Normalised Amplitude')
legend({roi_names{:} 'CanHRF'})

eval(sprintf('print -dtiff -f%d %s',sf3.Number,fullfile('Graphics','SVD_HRF.tif')))

% figure,plot(FIRcan(:,1))
% axis([0 32 -0.2 1]); set(gca,'YTick',[])


%% Supplementary Figure 4 - residuals vary with age?
sf4 = figure('OuterPosition',[100 100 1100 1100]);  % Supplementary Figure 4
for m = 1:nmods
    for r = 1:nrois
        subplot(nmods,nrois,(m-1)*nmods+r)

        Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_SSE.csv.gz',roi_names{r},models{m})));
        
        % convert to RMSE
        if strcmp(models{m}(1:3),'NLF')
            Y = sqrt(Y/nFIR); %yvals = [0:0.01:0.07];
        else
            Y = sqrt(Y/nscan); %yvals = [0:0.2:1.2];
        end
        %max(Y)
        yvals = [0:0.2:1.2];
        
        plot(participants.Age, Y, '.', 'MarkerSize', 4);
        
        axis([18 88 0 max(yvals)])
        set(gca,'XTick',[18:10:88],'FontSize',10);
        set(gca,'YTick',yvals,'FontSize',10);
        
        [Rval,Pval]=corr(participants.Age,Y,'type','Spearman');
        legend(sprintf('R=%+3.2f, p=%3.2f',Rval,Pval),'FontSize',8,'Location','NorthEast');
        
        if m == 1
            title(roi_names{r})
        elseif m  == nmods           
            xlabel('Age (years)','FontSize',12);
        end
        if r == 1
            ylabel(sprintf('%s: RMSE',models{m}),'FontSize',12)
        end
        
        drawnow;
    end
end

eval(sprintf('print -dtiff -f%d %s',sf4.Number,fullfile(out_dir,'Graphics','residuals.tif')))



%% Supplementary Figure 5 NLF parameters against age
m = 3;
nparams = 4;
nlabels = {'Lat. Off.','Lat. Scl.','Amp. Off.','Amp. Scl.'}; % just avoiding '_' making subscript in Matlab below
sf5 = figure('OuterPosition',[100 100 1100 1100]); 

for r = 1:nrois
    
    Y = spm_load(fullfile(roi_dir,sprintf('%s_%s_beta.csv.gz',roi_names{r},models{m})));
    labels = strvcat(fields(Y));
    Y = struct2array(Y);
     
    for p = 1:nparams
        subplot(nparams, nrois, (p-1)*nrois + r)
        grid on % axis square
        
        plot(participants.Age, Y(:,p), '.', 'MarkerSize', 4);
        
        set(gca,'XTick',[18:10:88],'FontSize',10);
        
        [Rval,Pval]=corr(participants.Age, Y(:,p), 'type','Spearman');
        legend(sprintf('R=%+3.2f, p=%3.2f',Rval,Pval),'FontSize',8,'Location','NorthEast');
        
        if p == 1
            title(roi_names{r})
        elseif p == nparams
            xlabel('Age (years)','FontSize',12);
        end
        if r == 1
%            ylabel(sprintf('NLF4: %s',labels(p,:)),'FontSize',12)
            ylabel(sprintf('NLF4: %s',nlabels{p}),'FontSize',12)
        end
        
        drawnow;
    end
end

eval(sprintf('print -dtiff -f%d %s',sf5.Number,fullfile(out_dir,'Graphics','NLF4_age.tif')))


%% Supplementary Figure S6 - PEB HDM fits for main effect

params = {'efficacy','decay','transit'};
nparams = length(params);

reord = [3 1 2]; % reorder PEB.Pnames to match L->R order in figure

effs = 'Mean'; e = 1; 

load(fullfile(hdm_dir,sprintf('BMAs_%d.mat',nparams)));

sf6 = figure('OuterPosition',[100 100 1100 600]);
ax = [];
for r = 1:nrois
    ax(r) = subplot(1,4,r);
    
    Ep = BMAs{r}.Ep;
    Vp = diag(BMAs{r}.Cp);
    
    Ep = spm_unvec(Ep,PEBs{1}.Ep);
    Vp = spm_unvec(Vp,PEBs{1}.Ep);
    
    Ep = full(Ep(reord,e));
    Vp = full(Vp(reord,e));
    
    spm_plot_ci(Ep,Vp);
    %spm_plot_ci(Ep./sqrt(Vp),zeros(length(Ep),1)); ylabel('Zscore') % Z-scores since different scalings?
    
    set(gca,'XTickLabel',params,'FontSize',10,'XTickLabelRotation',90);
    title(roi_names{r});
end
linkaxes(ax,'xy');
eval(sprintf('print -dtiff -f%d %s',sf6.Number,fullfile(out_dir,'Graphics',sprintf('HDM%d_Params_%s.tif',nparam,effs))))



%% Supplementary Figure S7 - LOOCV for Can1-Can3
params = {};
params{1} = {'hrf'};
params{2} = {'hrf (with time derivative)'};
params{3} = {'hrf (with time and dispersion derivatives)'};
nparams = [1 2 3];

label = {};
for m = 1:length(params)
    label{m} = sprintf('Can%d',nparams(m));
end

sf7 = figure('OuterPosition',[100 100 1100 600]);

MM = []; %AgeR2 = [];
for r = 1:nrois
     err = nan(nparticipants,length(params));   
     for m = 1:length(params)
         Y = spm_load(fullfile(roi_dir,sprintf('%s_CAN%d_beta.csv.gz',roi_names{r},nparams(m))));      
         
         if m>1; Y = struct2array(Y); end % when only one column, spm_load ignores header!
        
%        [~,~,~,~,AgeR2(p,r)] = glm(participants.Age,[Y ones(nparticipants,1)],[eye(nparams(m)) zeros(nparams(m),1)]',0);
         
        for s = 1:nparticipants
            sidx = [1:nparticipants]; sidx(s)=[];
            X = [Y(sidx,:) ones(nparticipants-1,1)];
            B = pinv(X)*participants.Age(sidx);
            pred_age = [Y(s,:) 1]*B;
            err(s,m) = participants.Age(s) - pred_age;
        end
     end
     subplot(1,nrois,r),hold on
     
     line(kron(ones(nparticipants,1),[1:size(err,2)])',abs(err'),'Color',[0.5 0.5 0.5])
     err = abs(err);
     boxplot(err); % bar(mean(err))
     set(gca,'XTickLabel',label,'FontSize',10);
     
     if r==1; ylabel('Abs Error (Years)');end
     % Compare first model with all others
     p = nan(1,size(err,2)); for m = 2:size(err,2); p(m) = signrank(err(:,m-1),err(:,m)); end  % CARE - pairwise tests now (not against err(:,1))
     %p = nan(1,size(err,2)); for m = 2:size(err,2); [~,p(m)] = ttest([err(:,1) - err(:,m)]); end % difference quite Gaussian even if individuals not
     f = find(p < .05/(length(params)-1)); 
     MM(:,r) = median(err)';
     for ff = f; text(ff-0.1,MM(ff,r),'*','FontSize',12,'Color',[0 1 0]); end
     title(roi_names{r});
end
MM
mean(MM,2)

eval(sprintf('print -dtiff -f%d %s',sf7.Number,fullfile(out_dir,'Graphics','Can_LOOCV.tif')))


%% Supplementary Figure S8 - LOOCV for HDM3-5
params = {};
params{1} = {'efficacy','decay','transit'};
params{2} = {'efficacy','decay','transit','alpha'};
params{3} = {'efficacy','decay','transit','alpha','feedback'};
%params{4} = {'efficacy','decay','transit','feedback','alpha','E0'};

label = {};
for m = 1:length(params)
    label{m} = sprintf('HDM%d',length(params{m}));
end

sf8 = figure('OuterPosition',[100 100 1100 600]);
MM = []; %AgeR2 = [];
for r = 1:nrois
     err = nan(nparticipants,length(params));   
     for m = 1:length(params)
         Y = spm_load(fullfile(roi_dir,sprintf('%s_HDM%d_beta.csv.gz',roi_names{r},length(params{m}))));      
         Y = struct2array(Y);
        
%        [~,~,~,~,AgeR2(p,r)] = glm(participants.Age,[Y ones(nparticipants,1)],[eye(length(params{m})) zeros(length(params{m}),1)]',0);
         
        for s = 1:nparticipants
            sidx = [1:nparticipants]; sidx(s)=[];
            X = [Y(sidx,:) ones(nparticipants-1,1)];
            B = pinv(X)*participants.Age(sidx);
            pred_age = [Y(s,:) 1]*B;
            err(s,m) = participants.Age(s) - pred_age;
        end
     end
     subplot(1,nrois,r),hold on,
     
     err = abs(err);
     line(kron(ones(nparticipants,1),[1:size(err,2)])',abs(err'),'Color',[0.5 0.5 0.5])
     boxplot(err); % bar(mean(err))
     set(gca,'XTickLabel',label,'FontSize',10);
     
     if r==1; ylabel('Abs Error (Years)');end
     % Compare first model with all others
     p = nan(1,size(err,2)); for m = 2:size(err,2); p(m) = signrank(err(:,m-1),err(:,m)); end % CARE - pairwise tests now (not against err(:,1))
     %p = nan(1,size(err,2)); for m = 2:size(err,2); [~,p(m)] = ttest([err(:,1) - err(:,m)]); end % difference quite Gaussian even if individuals not
     f = find(p < .05/(length(params)-1)); 
     MM(:,r) = median(err)';
     for ff = f; text(ff-0.1,MM(ff,r),'*','FontSize',12,'Color',[0 1 0]); end
     title(roi_names{r});
end
MM
mean(MM,2)

eval(sprintf('print -dtiff -f%d %s',sf8.Number,fullfile(out_dir,'Graphics','HDM_LOOCV.tif')))



