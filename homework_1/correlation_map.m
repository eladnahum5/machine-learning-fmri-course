close all
clear 
clc
%% set paths 
fMRI_folder_path = 'C:\GitHub\machine-learning-fmri-course\Pieman Story\story_selected_reg_wm_csf_hsd_last_half2\'; % add the fmri mat folder path
nii_path = 'C:\GitHub\machine-learning-fmri-course\Pieman Story\VOI  nifti files for HW1-20230329\MNI152_T1_3mm_brain.nii';
VOI_path ='C:\GitHub\machine-learning-fmri-course\homework_1\VOI  nifti files for HW1-20230329\a1_group2_new.voi';

%% find FMRI data and load nifti and VOI file
try
    fmri_data_path = fMRI_folder_path ;
    fmri_files = dir([fmri_data_path '*.mat']); % load fMRI data - mat files
    fmriAVG_data_path = [fMRI_folder_path '\avg\'];    
    fmri_avg_files = dir([fmriAVG_data_path '*.mat']); % load average fMRI data - mat files
    if (length(fmri_avg_files)~=18) 
        disp('Cant load files, please check fMRI folder path ')  
    end 
catch
    disp('Cant load files, please check fMRI folder path ')
end    

try
    ROI_template = xff(nii_path); % load nifti file
catch
    disp('Cant load files, please check nii path ')
end   

try
    voi = xff(VOI_path); % load VOI file 
catch
    disp('Cant load files, please check VOI path ')
end      


%%  Select VOI and change data to 1D vector 
select_voi = 5; % user can select VOI 1-5; 
roinii = voi.CreateMSK(ROI_template,select_voi); %this function create VOI(mask)
roi = single(roinii.VoxelData==1); %Extract data loction
roimask = single(reshape(roi,[(size(roi,1)*size(roi,2)*size(roi,3)),size(roi,4)]));  % reshape the 3d matrix to 1D vector (roimask) 


%% create correlation map 
method = 0 ;% Functional connectivity (FC) = 1, Inter-subject functional correlation (ISFC) = 0 
Nsub = 18;  % number of subjects = 18 
Nsamp =280; % TR
threshold = 6000; % signal threshold

for subject = 1:Nsub
    disp(['Analysis Subject ', num2str(subject)])
    load (fullfile(fmri_data_path, fmri_files(subject).name)); %load fmri data of subject -data_crop
    
    bold_one_temp=data_crop';

    if subject == 1   
        [Nsamp , Nvox] =size(bold_one_temp);
        csub_2 = NaN(Nvox,Nsub);  %avg corr for each subject (initialized to NaN)
    end
    
    mask_single = mean(bold_one_temp) > threshold;  %find bad voxels with low mean
    
    bold_one_temp(:,~mask_single)=NaN;  %set bad voxel as NAN

    ROI=bold_one_temp(:,logical(roimask)); % extract only VOI 
    
    ROI_norm =(ROI - mean(ROI))./std(ROI); % normalized data
    
    B=nanmean(ROI_norm,2); %mean data 
    
    
    if method==1 % FC
        gg_avg = find(mask_single);  %list of good voxels
        Ng_avg = length(gg_avg);   %number of good voxels in the average bold
        B_avg = bold_one_temp(1:Nsamp,gg_avg);  %pull out only the good voxels
    else  %ISFC
        load (fullfile(fmriAVG_data_path, fmri_avg_files(subject).name)); %load fmri avg data of subject - bold_avg
        mask_avg=bold_avg(end,:);    % get an average mask of the z-score average data, from the last row of the bold_avg response
        bold_avg(end,:)=[];
        bold_avg=bold_avg(1:Nsamp,:);

        gg_avg = find(mask_avg);  %list of good voxels
        Ng_avg = length(gg_avg);   %number of good voxels in the average bold
        B_avg = bold_avg(: ,gg_avg);  %pull out only the good voxels      
    end
    
    B = (B - mean(B)) ./(sqrt(Nsamp-1).*std(B)); %remove the mean of each column and convert to z-score units
    
    B_avg = (B_avg - mean(B_avg)) ./((sqrt(Nsamp-1).*std(B_avg))); %remove the mean of each column and convert to z-score units

    cc = (B_avg'*B);
    csub_2(gg_avg,subject) = cc';
    
end 
%% save corrlation map 
method_name = {'ISFC','FC' };

avg_map_dmn=nanmean(csub_2,2);
mean_corr_img = reshape(avg_map_dmn,[61, 73, 61]);
fprintf('Saving maps\n');
nii = xff(nii_path);
nii.VoxelData = mean_corr_img;
nii.VoxelData(isnan(nii.VoxelData)) = 0;
save_name='TaskData_ISFC_A1_GROUP_NEW_VOI_5_A1_L.nii';
nii.SaveAs(save_name);
disp(['Saved ' save_name])

% neuroelf_gui