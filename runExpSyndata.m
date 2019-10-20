clear;close all;
addpath('utils','algorithm','algorithm/nuclear_active','algorithm/Maxide');

isSaving = 0;
exp_times = 100;
m=100;n=100;r = 10;
recoveryCriterion = 0.001;
p_values =[ 0.09 0.1 0.12 0.15 0.18 0.21 0.23 0.25 0.28 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8];





tic
delta =0.1;
name = 'noisy_d0.1';
test_syndata(m,n,r,recoveryCriterion,delta,p_values,exp_times,name,isSaving);
toc






