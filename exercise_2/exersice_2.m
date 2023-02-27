%% Tutorial 2 
close all
clear all
clc

%% load Pie-man story 
[y,Fs] = audioread('exercise_2/Story_Original_MRI.wav');
%sound(y,Fs);% listen your audio input
N = length(y); % sample lenth
slength = N/Fs; % total time span of audio signal
time = linspace(0, N/Fs, N);
signal = mean(y,2);
signal_nor = normalize(signal);

%envelope
% [env_up ,env_down]= envelope(signal,100,'analytic');
[env_up ,env_down]=envelope(signal_nor,Fs/10,'peak');

% Resample 
T1 = resample(env_up,1,Fs);
T = resample(T1,2,3);
T = T(10:289,1);

figure()
hold on
plot(time,signal_nor)
plot ( time,env_up,'-r')
ylabel('Magnitude');
xlabel('Time');
xlim([150 200])
title('Envelope Pieman Audio Signal')
hold off
%% fMRI data - extract average ROI 
roinii1 = xff('exercise_2\MNI152_T1_3mm_brain.nii');
voi = xff('exercise_2\a1_group2_new.voi');   

roinii = voi.CreateMSK(roinii1,1);    % this function create mask:set to 1 the voxels indicies in 3d that match the MNI
roi = single(roinii.VoxelData==1);    % binarize it   
roimask = single(reshape(roi,[(size(roi,1)*size(roi,2)*size(roi,3)),size(roi,4)])); % reshape the 3d matrix to 1D vector (roimask)  
load("exercise_2\A0.mat")
data =bold_avg;
ROI=data(:,logical(roimask));
ROI_nor = (ROI - mean(ROI))./std(ROI);
mean_ROI(:,1) = mean(ROI_nor,2);

% fmri_data_path = 'exercise_2\' ;% add the fmri mat folder path
% fmri_files = dir([fmri_data_path '*.mat']); % load fMRI data - mat filesdict=dir(path);
% for i=1:18
% name=fmri_files(1).name;
% load([fmri_data_path name])
% data =data_crop';
% ROI=data(:,logical(roimask));
% ROI_nor = (ROI - mean(ROI))./std(ROI);
% mean_ROI(:,1) = mean(ROI_nor,2);
% end 

% mean_ROI_2 = mean (mean_ROI,2);

%% cross correlation

[r ,lag]= xcorr(T,mean_ROI(1:end-1),'normalized'); %calculate the correlation 

figure()
plot (lag,r)
xlim([-200 200])
ylabel('Magnitude');
xlabel('Time');
title('Cross Correlation Between BOLD fMRI & Envelope Audio');