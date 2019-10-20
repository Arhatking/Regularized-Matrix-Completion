function [ERR_wmc,ERR_rmc,ERR_imc] = test_syndata(m,n,r,recoveryCriterion,delta,p_values,exp_times,name,isSaving)

%%

isDisp =0;
iter_num =100;
mu1 = 0.001;mu2 = 0.001;

p_num = numel(p_values);

SUC_mc = zeros(p_num,exp_times);
ERR_mc = zeros(p_num,exp_times);
SUC_wmc = zeros(p_num,exp_times);
ERR_wmc = zeros(p_num,exp_times);
SUC_rmc = zeros(p_num,exp_times);
ERR_rmc = zeros(p_num,exp_times);
SUC_imc = zeros(p_num,exp_times);
ERR_imc = zeros(p_num,exp_times);
SUC_dirtyImc = zeros(p_num,exp_times);
ERR_dirtyImc = zeros(p_num,exp_times);

A = randn(m,r);
B = randn(n,r);
Mdata = 100*A*B';
Noise = randn(m,n);
Mdata = Mdata +delta*norm(Mdata,'fro')/norm(Noise,'fro')*Noise;


A_est = eye(m);
Noise = randn(n,r);
B_est = B;
B_est = B_est +delta*norm(B_est,'fro')/norm(Noise,'fro')*Noise;

PV = pinv(B')*B';
PV_est = pinv(B_est')*B_est';
PVc_est = PV_est - eye(n);

for i = 1:p_num    
    p = p_values(i);
    
    parfor j = 1:exp_times
        
        
        Omega_idx = rand(m,n);
        Omega = Omega_idx>1-p;
        Obs = Mdata.*Omega;
        
        try
            L = alm_MC(Obs, Omega,0.001,iter_num,0);
            err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
            ERR_mc(i,j) = err;
            if err < recoveryCriterion
                SUC_mc(i,j) = 1;
            end 
        catch ErrorInfo
            disp(ErrorInfo);
        end
        
       %%
        try
            theta = subspace(PV, PV_est);
            w = 0.5;%sqrt(tan(theta));
            Q = w*PV_est+PVc_est;
            L = alm_WMC(Obs, Omega,Q,iter_num,mu1,mu2,0);
            err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
            ERR_wmc(i,j) = err;
            if err < recoveryCriterion
                SUC_wmc(i,j) = 1;
            end 
        catch ErrorInfo
            disp(ErrorInfo);
        end
        
               
       %%
        try
            lambda = 1*sqrt(n);
            L = alm_RMC(Obs, Omega,PVc_est,lambda,iter_num,0.001,0.001,0);
            err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
            ERR_rmc(i,j) = err;
            if err < recoveryCriterion
                SUC_rmc(i,j) = 1;
            end 
        catch ErrorInfo
            disp(ErrorInfo);
        end
        
        try
            [L,~]=Maxide(Obs,find(Omega),A_est,B_est,1e-6,5000);
            err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
            ERR_imc(i,j) = err;
            if err < recoveryCriterion
                SUC_imc(i,j) = 1;
            end 
        catch ErrorInfo
            disp(ErrorInfo);
        end
        
        try
            L = DirtyIMC(Obs,find(Omega),A_est,B_est,0.1,0.1,10,r);
            err = norm((L-Mdata),'fro')/norm(Mdata,'fro');
            ERR_dirtyImc(i,j) = err;
            if err < recoveryCriterion
                SUC_dirtyImc(i,j) = 1;
            end 
        catch ErrorInfo
            disp(ErrorInfo);
        end
            
        fprintf('sample rate = %f, seq = %d \n',p,j);
    end   
end






figure('numbertitle','off','name','SuccessRate_vs_p');
plot(p_values(:),100.0*sum(SUC_dirtyImc,2)/exp_times,'-.k','MarkerSize',8,'linewidth',1);
hold on;
plot(p_values(:),100.0*sum(SUC_imc,2)/exp_times,'-m', 'MarkerSize',8,'linewidth',1);
plot(p_values(:),100.0*sum(SUC_rmc,2)/exp_times,'-xr','MarkerSize',8,'linewidth',1);
plot(p_values(:),100.0*sum(SUC_wmc,2)/exp_times,'--b','MarkerSize',8,'linewidth',1);
plot(p_values(:),100.0*sum(SUC_mc, 2)/exp_times,'-+g','MarkerSize',8,'linewidth',1);
xlim([p_values(1) p_values(end)]);
legend('DirtyIMC','IMC','RMC','WMC','MC','Location','SouthEast');
xlabel({'$p$'},'Interpreter','latex','fontsize',12)
ylabel({'Success Rate(%)'},'fontsize',12)

if isSaving
    %print('-depsc',['figures\syndata_suc_',name,'.eps']);
    print(gcf,'-dpng',['figures\syndata_suc_',name,'.png']);
end

figure('numbertitle','off','name','Err_vs_p');
plot(p_values(:),mean(ERR_dirtyImc,2),'-.k','MarkerSize',8,'linewidth',1);
hold on;
plot(p_values(:),mean(ERR_imc,2),'-m','MarkerSize',8,'linewidth',1);
plot(p_values(:),mean(ERR_rmc,2),'-xr','MarkerSize',8,'linewidth',1);
plot(p_values(:),mean(ERR_wmc,2),'--b','MarkerSize',8,'linewidth',1);
plot(p_values(:),mean(ERR_mc, 2),'-+g','MarkerSize',8,'linewidth',1);
xlim([p_values(1) p_values(end)]);
legend('DirtyIMC','IMC','RMC','WMC','MC','Location','NorthEast');
xlabel({'$p$'},'Interpreter','latex','fontsize',12)
ylabel({'Relative Recovery Error'},'fontsize',12)

if isSaving
    print('-depsc',['figures\syndata_err_',name,'.eps']);
    print(gcf,'-dpng',['figures\syndata_err_',name,'.png']);
end

if isSaving
    save(['savings\syndata_',name,'.mat'],'exp_times','p_values','SUC_mc','SUC_rmc','SUC_wmc','SUC_imc','SUC_dirtyImc','ERR_mc','ERR_rmc','ERR_wmc','ERR_imc','ERR_dirtyImc','iter_num');
end
