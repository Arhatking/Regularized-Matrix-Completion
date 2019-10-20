
clear;
clc;
addpath('utils','algorithm','algorithm/Maxide','algorithm/nuclear_active','dataset/MultiLabel');

repetition = 1;
percent=0.4;

load('Education1.mat');

Mdata=Mdata(:,1:5000);
B = B(1:5000,:);
[m,n]=size(Mdata);
A = eye(m);

PV =  pinv(B')*B';
PVc =eye(n) - PV;
lambda = 10*sqrt(n);
w = 0.5;
Q = w*PV+PVc;

mu1_space = [-7:5];
mu2_space = [-10:5];
n1 = size(mu1_space,2);
n2 = size(mu2_space,2);

AP_wmc = zeros(n1,n2,repetition);
ERR_wmc = zeros(n1,n2,repetition);
AP_rmc = zeros(n1,n2,repetition);
ERR_rmc = zeros(n1,n2,repetition);

Omega=zeros(m,n);
positive_index=find(Mdata(:)>0);
positive_number=length(positive_index);
positive_random=randperm(positive_number);
positive_select=positive_index(positive_random(1,1:ceil(positive_number*percent)));
Omega(positive_select)=1;
negative_index=find(Mdata(:)<=0);
negative_number=length(negative_index);
negative_random=randperm(negative_number);
negative_select=negative_index(negative_random(1,1:ceil(negative_number*percent)));
Omega(negative_select)=1;
Obs=Mdata.*Omega;

tic
for i = 1:n1 
    mu1 = 2^mu1_space(i);
    parfor j = 1:n2 
        mu2 = 2^mu2_space(j);
        
        L = alm_RMC(Obs, Omega,PVc,lambda,5,mu1,mu2,0);
        AP_rmc(i,j)=AveragePrecision(L',Mdata');
        
        L = alm_WMC(Obs, Omega,Q,5,mu1,mu2,0);
        AP_wmc(i,j)=AveragePrecision(L',Mdata');
        
        fprintf('i = %d, j = %d\n',i,j);
    end   
end
toc

[X,Y]=meshgrid(mu1_space,mu2_space);
figure('numbertitle','off','name','RMC');
surf(X,Y,AP_rmc');
xlabel('\mu_1')
ylabel('\mu_2')

figure('numbertitle','off','name','WMC');
surf(X,Y,AP_wmc');
xlabel('\mu_1')
ylabel('\mu_2')
