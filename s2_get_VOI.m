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

bas_dir = '/imaging/rh01/CamCAN/700/HRF_SMT' % Change to your working directory

raw_dir = '/imaging/camcan/cc700/mri/pipeline/release004/data_fMRI_Smooth/aamod_smooth_00001/'; % Change to where you downloaded processed SMT images from CamCAN website above

git_dir = fullfile(bas_dir,'AgeingHRF') % https://github.com/RikHenson/AgeingHRF

spm_dir = '/imaging/local/software/spm_cbu_svn/releases/spm12_latest/' % Your local installation of SPM12
addpath(spm_dir);

addpath(fullfile(git_dir,'Matlab_Utils')) % spm_regions needed
spm('Defaults','fMRI')

out_dir = fullfile(bas_dir,'outputs'); % Where results will go

ana_dir = fullfile(out_dir,'SPM_FIR'); 

participants = spm_load(fullfile(git_dir,'participants.csv'));
nparticipants = length(participants.CCID)

roi_names = {'lAC','bVC','lMC','rMC'};
nrois = length(roi_names);

glms = {'Stim-locked','Resp-locked'};
glm_type = [1 1 2 2]; % 1 = stim-locked, 2 = resp-locked

%% Extract data using SPM VOI
% cluster_images created by hand from viewing "Age" contrast of Group analyses,
% moving to cluster and saving cluster as NII file

cluster_images{1} = fullfile(ana_dir,'Stim-locked','Group','lAC_Age_F12_559vox.nii'); 
cluster_images{2} = fullfile(ana_dir,'Stim-locked','Group','bVC_Age_F12_371vox.nii'); 
cluster_images{3} = fullfile(ana_dir,'Resp-locked','Group','lMC_Age_F12_219vox.nii'); 
cluster_images{4} = fullfile(ana_dir,'Resp-locked','Group','rMC_Age_F12_481vox.nii'); 

pval_thrs = [0.005 0.05 0.5 1]; % [] means no single-subject threshold

SS_pval = cell(1,nparticipants); Nvox = SS_pval;
for r = 1:nrois
    bwd = fullfile(ana_dir,glms{glm_type(r)});
   
    parfor s = 1:nparticipants
        ccid = sprintf('CC%d',participants.CCID(s)); disp(ccid)
        cd(fullfile(bwd,ccid))
        
        Nvox{s}(r) = 0;
        
        matlabbatch = {};
        matlabbatch{1}.spm.util.voi.spmmat = cellstr('SPM.mat');
        matlabbatch{1}.spm.util.voi.adjust = 1;
        matlabbatch{1}.spm.util.voi.session = 1;
        matlabbatch{1}.spm.util.voi.name = roi_names{r};
        matlabbatch{1}.spm.util.voi.roi{1}.mask.image = cellstr(cluster_images{r});
        matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
        matlabbatch{1}.spm.util.voi.expression = 'i1';
        
        if isempty(pval_thrs)
            spm_jobman('run',matlabbatch);
        else
            matlabbatch{1}.spm.util.voi.roi{2}.spm.spmmat = {''}; %cellstr(SPMfile);
            matlabbatch{1}.spm.util.voi.roi{2}.spm.contrast = 1;
            matlabbatch{1}.spm.util.voi.roi{2}.spm.conjunction = 1; % not used
            matlabbatch{1}.spm.util.voi.roi{2}.spm.threshdesc = 'none';
            matlabbatch{1}.spm.util.voi.roi{2}.spm.extent = 0;
            matlabbatch{1}.spm.util.voi.roi{2}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
            matlabbatch{1}.spm.util.voi.expression = 'i1&i2';
            ps = 0;
            while Nvox{s}(r) < 3
                ps = ps+1;
                try delete(fullfile(bwd,ccid,sprintf('*%s*',roi_names{r}))); end
                
                matlabbatch{1}.spm.util.voi.roi{2}.spm.thresh = pval_thrs(ps);
                spm_jobman('run',matlabbatch);
                
                try 
                    tmp = rikload(fullfile(bwd,ccid,sprintf('VOI_%s_1',roi_names{r})));
                    Nvox{s}(r) = size(tmp.xY.y,2);
                end
            end
            SS_pval{s}(r) = pval_thrs(ps);
        end
        
        delete(fullfile(bwd,ccid,'VOI_*.nii'))
    end
end

%% Check VOI p-value and Nvox across age

% SS_pval = cat(1,SS_pval{:});
% figure,hist(SS_pval)

%Nvox    = cat(1,Nvox{:});
Nvox = [];
for r = 1:nrois
    bwd = fullfile(ana_dir,glms{glm_type(r)});
    for s = 1:nparticipants
        ccid = sprintf('CC%d',participants.CCID(s)); 
        cd(fullfile(bwd,ccid))
        tmp = load(fullfile(bwd,ccid,sprintf('VOI_%s_1',roi_names{r})));
        Nvox(s,r) = size(tmp.xY.y,2);
    end
end
[min(Nvox); median(Nvox); max(Nvox)]

figure;
for r=1:nrois
    subplot(nrois,1,r),plot(participants.Age,Nvox(:,r),'o'),title(roi_names{r})
%    subplot(1,roi_names,r),plot(participants.Age,SS_pval(:,r),'o'),title(roi_names{r})
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
