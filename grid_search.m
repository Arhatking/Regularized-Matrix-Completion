
%%
% rand('seed',1);
% randn('seed',2);

isDisp =0;
iter_num =500;
m=100;n=100;r=10;
sample_rate = 0.5;
recoveryCriterion = 1e-3;
repetition = 20;


A = randn(m,r);
B = randn(n,r);
Mdata = 100*A*B';

A_est = eye(m);
Noise = randn(n,r);
B_est = B;
B_est = B_est +1e-3*norm(B_est,'fro')/norm(Noise,'fro')*Noise;
PV = pinv(B')*B';
PV_est = pinv(B_est')*B_est';
PVc_est = PV_est - eye(n);
theta = subspace(PV, PV_est);
w = sqrt(tan(theta));
Q = w*PV_est+PVc_est;
lambda = 10*sqrt(n);

mu1_space = [0.01 0.1:0.1:2];
mu2_space = [0.001:0.001:0.01];
n1 = size(mu1_space,2);
n2 = size(mu2_space,2);

SUC_wmc = zeros(n1,n2,repetition);
ERR_wmc = zeros(n1,n2,repetition);
SUC_rmc = zeros(n1,n2,repetition);
ERR_rmc = zeros(n1,n2,repetition);

tic
for i = 1:n1        
    for j = 1:n2
        mu1 = mu1_space(i);
        mu2 = mu2_space(j);
        parfor t = 1:repetition                 
            Omega_idx = rand(m,n);
            Omega = Omega_idx>1-sample_rate;
            Obs = Mdata.*Omega;
            try  
                L = alm_RMC(Obs, Omega,PVc_est,lambda,iter_num,mu1,mu2,0);
                err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
                ERR_rmc(i,j,t) = err;
                if err < recoveryCriterion
                    SUC_rmc(i,j,t) = 1;
                end 
            catch ErrorInfo
                disp(ErrorInfo);
            end            
        end
        fprintf('i = %d, j = %d\n',i,j);
    end   
end
toc

err = mean(ERR_rmc,3);
suc = 100*sum(SUC_rmc,3)/repetition;