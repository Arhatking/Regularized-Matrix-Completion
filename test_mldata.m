function [AP_wmc,AP_rmc,AP_imc,AP_mc,AP_dirtyImc] = test_mldata(p_values,exp_times,name,isSaving,isDisp)

%%


iter_num = 5;
mu1 = 0.1;mu2 = 0.1;

p_num = numel(p_values);

AP_mc = zeros(p_num,exp_times);
AP_wmc = zeros(p_num,exp_times);
AP_rmc = zeros(p_num,exp_times);
AP_imc = zeros(p_num,exp_times);
AP_dirtyImc = zeros(p_num,exp_times);


load([name,'.mat']);

Mdata=Mdata(:,1:5000);
B = B(1:5000,:);

[m,n]=size(Mdata);
A = eye(m);


PV =  pinv(B')*B';
PVc =eye(n) - PV;
lambda_rmc = 1*sqrt(n);
w = 0.1;
Q = w*PV+PVc;

positive_index=find(Mdata(:)>0);
positive_number=length(positive_index);
negative_index=find(Mdata(:)<=0);
negative_number=length(negative_index);

for i = 1:p_num    
    percent = p_values(i);
    parfor j = 1:exp_times
        %Generate the percent% randomly observed entries
        Omega=zeros(m,n);
        positive_random=randperm(positive_number);
        positive_select=positive_index(positive_random(1,1:ceil(positive_number*percent)));
        Omega(positive_select)=1;
        negative_random=randperm(negative_number);
        negative_select=negative_index(negative_random(1,1:ceil(negative_number*percent)));
        Omega(negative_select)=1;
        Omega_linear = find(Omega);
        Obs=Mdata.*Omega;
        
        if isDisp
            fprintf('sample rate = %f, seq = %d \n',percent,j);
        end
       %%
        
        [L,~]=Maxide(Obs,Omega_linear,A,B,4,100);
        AP_imc(i,j)=AveragePrecision(L',Mdata');
      
        L = alm_MC(Obs, Omega,0.1,iter_num,0);
        AP_mc(i,j)=AveragePrecision(L',Mdata');
        
        mu1 = 2^0;mu2 = 2^0;
        L = alm_WMC(Obs, Omega,Q,iter_num,mu1,mu2,0);
        AP_wmc(i,j)=AveragePrecision(L',Mdata');
        
        mu1 = 2^1;mu2 = 2^-7;
        L = alm_RMC(Obs, Omega,PVc,lambda_rmc,iter_num,mu1,mu2,0);
        AP_rmc(i,j)=AveragePrecision(L',Mdata');

               
        L = DirtyIMC(Obs,Omega_linear,A,B,0.1,0.1,10,10);
        AP_dirtyImc(i,j)=AveragePrecision(L',Mdata');
        
       %%
            
        
    end   
end






figure('numbertitle','off','name','AveragePrecision_vs_SampleRate');
plot(p_values(:),mean(AP_dirtyImc,2),'-.k','MarkerSize',8,'linewidth',1);
hold on;
plot(p_values(:),mean(AP_imc,2),'-m', 'MarkerSize',8,'linewidth',1);
plot(p_values(:),mean(AP_rmc,2),'-xr','MarkerSize',8,'linewidth',1);
plot(p_values(:),mean(AP_wmc,2),'--b','MarkerSize',8,'linewidth',1);
plot(p_values(:),mean(AP_mc, 2),'-+g','MarkerSize',8,'linewidth',1);
xlim([p_values(1) p_values(end)]);
ylim([0 1]);
legend('DirtyIMC','IMC','RMC','WMC','MC','Location','SouthEast');
xlabel({'Sampling Rate'},'Interpreter','latex','fontsize',12)
ylabel({'Average Precision'},'fontsize',12)

if isSaving
    %print('-depsc',['figures\syndata_suc_',name,'.eps']);
    print(gcf,'-dpng',['figures\',name,'.png']);
end



if isSaving
    save(['savings\',name,'.mat'],'exp_times','p_values','AP_mc','AP_rmc','AP_wmc','AP_imc','AP_dirtyImc');
end
