close all;
addpath('libs/PROPACK','utils','algorithm','algorithm/nuclear_active','algorithm/Maxide','dataset/WSN');

clims = [0,0.2];
iter_num =100;
mu1 = 0.1;mu2 = 0.1;
p = 0.1;

rand('seed',123);
randn('seed',123);

load('intel_hum.mat');
Rawdata = Mdata;


Mdata=Mdata(:,1:288);
[m,n]=size(Mdata);
A = eye(m);
B = dctmtx(288)';
B = B(:,1:30);

PV =  pinv(B')*B';
PVc =eye(n) - PV;
w = 0.5;
Q = w*PV+PVc;

offset = randsample(288*9,1);
Mdata = Rawdata(:,1+offset:288+offset);
Omega_idx = rand(m,n);
Omega = Omega_idx>1-p;
Obs = Mdata.*Omega;



L = alm_MC(Obs, Omega,0.1,iter_num,0);
err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
fprintf('MC Error = %.3f\n',err);
ErrMtx = abs(Mdata-L)./abs(Mdata);
figure('numbertitle','off','name','MC');
colormap('Jet');
imagesc(ErrMtx,clims); 
set (gcf,'units','pixel','Position',[5,400,200,200]);
set(gca, 'position', [0 0 1 1 ]);axis normal;
AFrame=getframe(gcf);
imwrite(AFrame.cdata,'figures\DrawErrMtx\errmtx_mc.png');

L = alm_WMC(Obs, Omega,Q,iter_num,mu1,mu2,0);
err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
fprintf('WMC Error = %.3f\n',err);
ErrMtx = abs(Mdata-L)./abs(Mdata);
figure('numbertitle','off','name','WMC');
colormap('Jet');
imagesc(ErrMtx,clims); 
set (gcf,'units','pixel','Position',[210,400,200,200]);
set(gca, 'position', [0 0 1 1 ]);axis normal;
AFrame=getframe(gcf);
imwrite(AFrame.cdata,'figures\DrawErrMtx\errmtx_wmc.png');

lambda = 1*sqrt(n);
L = alm_RMC(Obs, Omega,PVc,lambda,iter_num,mu1,mu2,0);
err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
fprintf('RMC Error = %.3f\n',err);            
ErrMtx = abs(Mdata-L)./abs(Mdata);
figure('numbertitle','off','name','RMC');
colormap('Jet');
imagesc(ErrMtx,clims); 
set (gcf,'units','pixel','Position',[420,400,200,200]);
set(gca, 'position', [0 0 1 1 ]);axis normal;
AFrame=getframe(gcf);
imwrite(AFrame.cdata,'figures\DrawErrMtx\errmtx_rmc.png');



[L,~]=Maxide(Obs,find(Omega),A,B,2,iter_num);               
err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
fprintf('IMC Error = %.3f\n',err);
ErrMtx = abs(Mdata-L)./abs(Mdata);
figure('numbertitle','off','name','IMC');
colormap('Jet');
imagesc(ErrMtx,clims); 
set (gcf,'units','pixel','Position',[630,400,200,200]);
set(gca, 'position', [0 0 1 1 ]);axis normal;
AFrame=getframe(gcf);
imwrite(AFrame.cdata,'figures\DrawErrMtx\errmtx_imc.png');

L = DirtyIMC(Obs,find(Omega),A,B,0.1,0.1,10,10);
err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
fprintf('DirtyIMC Error = %.3f\n',err);
ErrMtx = abs(Mdata-L)./abs(Mdata);
figure('numbertitle','off','name','DirtyIMC');
colormap('Jet');
imagesc(ErrMtx,clims); 

set (gcf,'units','pixel','Position',[830,400,200,200]);
set(gca, 'position', [0 0 1 1 ]);axis normal;
AFrame=getframe(gcf);
imwrite(AFrame.cdata,'figures\DrawErrMtx\errmtx_dimc.png');
