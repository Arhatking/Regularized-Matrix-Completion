clc;clear;close all;
addpath('utils','algorithm','algorithm/nuclear_active','algorithm/Maxide','dataset/MultiLabel');

isSaving = 0;
exp_times = 100;
p_values =[0.1:0.1:0.8];
isDisp = 1;
lists = {'Arts1','Business1','Computers1','Education1','Entertainment1','Health1','Recreation1','Reference1','Science1','Social1','Society1'};
for i = 1:size(lists,2)
    name = lists{1,i};
    disp(['----------------Running experiment on ' name ' dataset: ---------------']);
    tic
    test_mldata(p_values,exp_times,name,isSaving,isDisp);
    toc
    disp(['----------------Experiment on ' name ' dataset is completed: ---------------']);
end



