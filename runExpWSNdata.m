clc;clear;close all;
addpath('libs/PROPACK','utils','algorithm','algorithm/nuclear_active','algorithm/Maxide','dataset/WSN');

isSaving = 1;
exp_times = 500;
p_values =[0.1:0.05:0.5,0.6:0.1:0.8];
isDisp = 1;
recoveryCriterion = 0.01;
disp(['----------------Running experiment on WSN dataset: ---------------']);
tic
test_wsndata('intel_hum',recoveryCriterion,p_values,exp_times,isSaving);
toc
disp(['----------------Experiment on WSN dataset is completed: ---------------']);


