%% Importing the data
clear all; clc

nii = xff('exercise_1\sub-025-task-intact2.nii');
voi = xff('exercise_1\a1_group2_new.voi');
data = nii.VoxelData;

%% Get VOIs timecourse
nii_voitc = nii.VOITimeCourse(voi, struct('weight', 3));
nii_voitc_table = cell2table(nii_voitc);

%% Get timecourse average
[B, voin] = nii.VOITimeCourse(voi);
B_table = array2table(B, "VariableNames", voin);
figure(1);
stackedplot(B_table)

%% Calculating correlation
figure(2);
B_corr = corrcoef(B);
imagesc(B_corr);
xticks(1:5)
xticklabels(voin);
yticks(1:5)
yticklabels(voin);
title('Correlation Matrix', "FontSize", 14);
colorbar;