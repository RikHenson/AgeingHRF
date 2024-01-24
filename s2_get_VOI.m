% Matlab script for extracting ROI data using SPM VOI functions:
%
%  Download/clone to "bas_dir" below:
%
%       1. "AgeingHRF" from https://github.com/RikHenson/AgeingHRF/
%
%   Link to local SPM12 directory, or install from 
%       https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ 
%
% rik.henson@mrc-cbu.cam.ac.uk, Jan 2023

clear

bas_dir = '/imaging/rh01/CamCAN/700/HRF_SMT/Revision' % Change to your working directory

git_dir = fullfile(bas_dir,'AgeingHRF-main') % https://github.com/RikHenson/AgeingHRF

spm_dir = '/imaging/local/software/spm_cbu_svn/releases/spm12_latest/' % Your local installation of SPM12
addpath(spm_dir);

addpath(fullfile(git_dir,'Matlab_Utils')) % spm_regions needed
spm('Defaults','fMRI')

out_dir = fullfile(bas_dir,'outputs'); % Where results will go (should have been created in s1_fit_SPM_FIR.m)
ana_dir = fullfile(out_dir,'SPM_FIR'); % Where results will go (should have been created in s1_fit_SPM_FIR.m)

participants = spm_load(fullfile(git_dir,'participants.csv'));
nparticipants = length(participants.CCID)

roi_names = {'lAC','bVC','lMC','rMC'};
nrois = length(roi_names);

glms = {'Stim-locked','Resp-locked'};
glm_type = [1 1 2 2]; % 1 = stim-locked, 2 = resp-locked


%% Extract data using SPM VOI, using mask images from s1_fit_SPM_FIR.m

cluster_images{1} = fullfile(ana_dir,'Stim-locked','Group','lAC_Age_F5_280vox.nii'); 
cluster_images{2} = fullfile(ana_dir,'Stim-locked','Group','bVC_Age_F5_182vox.nii'); 
cluster_images{3} = fullfile(ana_dir,'Resp-locked','Group','lMC_Age_F5_54vox.nii'); 
cluster_images{4} = fullfile(ana_dir,'Resp-locked','Group','rMC_Age_F5_88vox.nii'); 

%pval_thrs = []; % no single-subject threshold
pval_thrs = [0.05]; %[0.005 0.05 0.5 1]; 

Nvox = cell(1,nparticipants); 
for r = 1:nrois
    bwd = fullfile(ana_dir,glms{glm_type(r)});
   
    parfor s = 1:nparticipants
        ccid = sprintf('CC%d',participants.CCID(s)); disp(ccid)
        cd(fullfile(bwd,ccid))
        
        Nvox{s}(r) = NaN;
        
        matlabbatch = {};
        matlabbatch{1}.spm.util.voi.spmmat = cellstr('SPM.mat');
        matlabbatch{1}.spm.util.voi.adjust = 1;
        matlabbatch{1}.spm.util.voi.session = 1;
        matlabbatch{1}.spm.util.voi.name = roi_names{r};
        
        matlabbatch{1}.spm.util.voi.roi{1}.mask.image = cellstr(cluster_images{r});
        matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
        
        if isempty(pval_thrs)
            matlabbatch{1}.spm.util.voi.expression = 'i1';
            spm_jobman('run',matlabbatch);
        else
            matlabbatch{1}.spm.util.voi.roi{end+1}.spm.spmmat = {''}; %cellstr(SPMfile);
            matlabbatch{1}.spm.util.voi.roi{end}.spm.contrast = 1;
            matlabbatch{1}.spm.util.voi.roi{end}.spm.conjunction = 1; % not used
            matlabbatch{1}.spm.util.voi.roi{end}.spm.threshdesc = 'none';
            matlabbatch{1}.spm.util.voi.roi{end}.spm.extent = 0;
            matlabbatch{1}.spm.util.voi.roi{end}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
            matlabbatch{1}.spm.util.voi.expression = 'i1&i2';
            ps = 0;
            while isnan(Nvox{s}(r)) & ps < length(pval_thrs) %& Nvox{s}(r) < 3
                ps = ps+1;
                try delete(fullfile(bwd,ccid,sprintf('*%s*',roi_names{r}))); end
                
                matlabbatch{1}.spm.util.voi.roi{end}.spm.thresh = pval_thrs(ps);
                spm_jobman('run',matlabbatch);
                
                try 
                    tmp = rikload(fullfile(bwd,ccid,sprintf('VOI_%s_1',roi_names{r})));
                    Nvox{s}(r) = size(tmp.xY.y,2);
                catch
                    Nvox{s}(r) = 0;
                end
            end
        end
        
        delete(fullfile(bwd,ccid,'VOI_*.nii'))
    end
end
 
%% Exclude participants with no voxels in ROI

Nvox = cat(1,Nvox{:});
[s,r] = find(Nvox==0);
exclude = unique(s);
exclude_ccid = participants.CCID(exclude)
exclude_ages = participants.Age(exclude)

% See why excluded?
% for s = exclude % (s=148, CC221585 has horrible waves!)
%     ccid = sprintf('CC%d',participants.CCID(s)); disp(ccid)
%     cd(fullfile(ana_dir,'Stim-locked',ccid))
%     load SPM
%     figure,imagesc(zscore(SPM.xX.X))
%     spm_movie('load',1,SPM.xY.P,27,1,0);
% end

save(fullfile(out_dir,'Excluded_Participants'),'s','r','exclude','exclude_ccid','exclude_ages','Nvox')

%% Write out valid participants

include = setdiff([1:nparticipants],exclude);
length(include)

fp = fopen(fullfile(git_dir,'participants_include.csv'),'w');
fprintf(fp,'CCID,Age\n');
for s = include
    fprintf(fp,'%d,%d\n',participants.CCID(s),participants.Age(s));
end
fclose(fp);

%% Collect Nvox and Pvar

participants = spm_load(fullfile(git_dir,'participants_include.csv'));
nparticipants = length(participants.CCID)

Nvox = []; Pvar = [];
for r = 1:nrois
    bwd = fullfile(ana_dir,glms{glm_type(r)});
    for s = 1:nparticipants
        ccid = sprintf('CC%d',participants.CCID(s));
        cd(fullfile(bwd,ccid))
        tmp = load(fullfile(bwd,ccid,sprintf('VOI_%s_1',roi_names{r})));
        Nvox(s,r) = size(tmp.xY.y,2);
        Pvar(s,r) = tmp.xY.s(1).^2/sum(tmp.xY.s.^2);
    end
end
[min(Nvox); median(Nvox); max(Nvox)]
[min(Pvar); median(Pvar); max(Pvar)]

figure;
for r=1:nrois
%    subplot(nrois,1,r),plot(participants.Age,Pvar(:,r),'o'),title(roi_names{r}),ylabel('Pvar')
    subplot(nrois,1,r),plot(participants.Age,Nvox(:,r),'o'),title(roi_names{r}),ylabel('Nvox')
%    subplot(nrois,1,r),plot(Nvox(:,r),Pvar(:,r),'o'),title(roi_names{r}),xlabel('Nvox'),ylabel('Pvar')
%%    subplot(1,roi_names,r),plot(participants.Age,SS_pval(:,r),'o'),title(roi_names{r})
    [R,p]=corr(participants.Age,Nvox(:,r)); fprintf('%s: Age-Nvox, R=%3.2f, p=%4.3f\n',roi_names{r},R,p)
end



%% Save BOLD ROI timeseries to CSV for more general access

roi_dir = fullfile(git_dir,'ROI_data') % CSVs (or could use VOI*mat files)
nscans  = 261;

for r = 1:nrois
    glm_dir = glms{glm_type(r)};
    
    fn_bold = fullfile(roi_dir,sprintf('%s_BOLD.csv',roi_names{r}));
    fp_bold = fopen(fn_bold,'w');
    for v = 1:(nscans-1); fprintf(fp_bold,'vol%d,',v); end; fprintf(fp_bold,sprintf('vol%d\n',nscans));
    
    for s = 1:nparticipants 
        name = sprintf('CC%d',participants.CCID(s));
        swd = fullfile(ana_dir,glm_dir,name);

        % Load the prewhitened adjusted timeseries (VOI_XX.mat)
        xY  = load(fullfile(swd,sprintf('VOI_%s_1.mat', roi_names{r})));
        
        for v = 1:(nscans-1); fprintf(fp_bold,'%6.5f,',xY.Y(v)); end; fprintf(fp_bold,'%6.5f\n',xY.Y(nscans));
        fprintf('.')
    end
    fprintf('\n')
    
    fclose(fp_bold);     
    gzip(fn_bold); 
    delete(fn_bold);
end


return
