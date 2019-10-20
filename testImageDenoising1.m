clc;clear;close all;
addpath('libs/PROPACK','utils','algorithm','algorithm/Maxide','algorithm/nuclear_active');



ImageRaw = imread('dataset/images/10.png');

Mdata = im2double(ImageRaw(:,:,1));
figure('numbertitle','off','name','Orignal Image');
imshow(Mdata)
set (gcf,'units','pixel','Position',[5,400,256,256]);
set(gca, 'position', [0 0 1 1 ]);axis normal;
AFrame=getframe(gcf);
imwrite(AFrame.cdata,'figures\MR_image\sample_rate_15_orig.png');
%close(gcf);

m = size(Mdata,1);
n = size(Mdata,2);

%%


%%
iter_num =100;
sample_rate = 0.15;





Omega_idx = rand(m,n);
Omega = Omega_idx>1 - sample_rate;
Omega_c = ones(m,n) - Omega;
Obs = Mdata.*Omega;

figure('numbertitle','off','name','Undersampled Image');
imshow(Obs);
set (gcf,'units','pixel','Position',[266,400,256,256]);
set(gca, 'position', [0 0 1 1 ]);axis normal;
AFrame=getframe(gcf);
imwrite(AFrame.cdata,'figures\MR_image\sample_rate_15_obsv.png');


mu1 = 0.5;mu2 = 0.01;
delta = 0.00;
r = 50;
[A s B] =lansvd(Mdata,r);
A_est = eye(m);
B_est=B;
Noise = randn(n,r);
B_est = B_est+delta*norm(B_est,'fro')/norm(Noise,'fro')*Noise;
PV = pinv(B_est')*B_est';
PVc =eye(n) - pinv(B_est')*B_est';


L = alm_MC(Obs, Omega, mu1, iter_num,0);
err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
fprintf('MC Error = %.3f\n',err);
figure('numbertitle','off','name','MC');
set(gca,'LooseInset',get(gca,'TightInset'))
imshow(L);
set (gcf,'units','pixel','Position',[266*2,400,256,256]);
set(gca, 'position', [0 0 1 1 ]);axis normal;
AFrame=getframe(gcf);
imwrite(AFrame.cdata,'figures\MR_image\sample_rate_15_mc.png');



Omega_linear = find(Omega);
[L,telapsed_side]=Maxide(Obs,Omega_linear,A_est,B_est,1e-2,iter_num);
err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
fprintf('IMC Error = %.3f\n',err);
figure('numbertitle','off','name','IMC');
imshow(L);
set (gcf,'units','pixel','Position',[5,650,256,256]);
set(gca, 'position', [0 0 1 1 ]);axis normal;
AFrame=getframe(gcf);
imwrite(AFrame.cdata,'figures\MR_image\sample_rate_15_imc.png');



L = DirtyIMC(Obs,Omega_linear,A_est,B_est,0.01,1,10,20);
err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
fprintf('DirtyIMC Error = %.3f\n',err);
figure('numbertitle','off','name','DirtyIMC');
imshow(L);
set (gcf,'units','pixel','Position',[266,650,256,256]);
set(gca, 'position', [0 0 1 1 ]);axis normal;
AFrame=getframe(gcf);
imwrite(AFrame.cdata,'figures\MR_image\sample_rate_15_dirtyIMC.png');



lambda = sqrt(n);
L = alm_RMC(Obs, Omega,PVc,lambda,iter_num,mu1,mu2,0);
err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
fprintf('RMC Error = %.3f\n',err);
figure('numbertitle','off','name','RMC');
imshow(L);
set (gcf,'units','pixel','Position',[266*2,50,256,256]);
set(gca, 'position', [0 0 1 1 ]);axis normal;
AFrame=getframe(gcf);
imwrite(AFrame.cdata,'figures\MR_image\sample_rate_15_rmc.png');

w=0.3;
Q = w*PV+PVc;
L = alm_WMC(Obs, Omega,Q,iter_num,mu1,mu2,0);
err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
fprintf('WMC Error = %.3f\n',err);
figure('numbertitle','off','name','WMC');
imshow(L);
set (gcf,'units','pixel','Position',[266*3,50,256,256]);
set(gca, 'position', [0 0 1 1 ]);axis normal;
AFrame=getframe(gcf);
imwrite(AFrame.cdata,'figures\MR_image\sample_rate_15_wmc.png');

