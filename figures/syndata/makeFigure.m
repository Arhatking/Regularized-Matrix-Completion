close all
load('../../savings/syndata_n100_r10_d0.mat')
figure('numbertitle','off','name','SuccessRate_vs_p');
plot(p_values(:),100.0*sum(SUC_dirtyImc,2)/exp_times,'-.k','MarkerSize',8,'linewidth',2);
hold on;
plot(p_values(:),100.0*sum(SUC_imc,2)/exp_times,'-m', 'MarkerSize',8,'linewidth',2);
plot(p_values(:),100.0*sum(SUC_rmc,2)/exp_times,'-xr','MarkerSize',8,'linewidth',2);
plot(p_values(:),100.0*sum(SUC_wmc,2)/exp_times,'--b','MarkerSize',8,'linewidth',2);
plot(p_values(:),100.0*sum(SUC_mc, 2)/exp_times,'-+g','MarkerSize',8,'linewidth',2);
xlim([p_values(1) p_values(end)]);
legend('DirtyIMC','IMC','RMC','WMC','MC','Location','NorthEast');
xlabel({'$p$'},'Interpreter','latex','fontsize',12)
ylabel({'Success Rate(%)'},'fontsize',12)
print('-depsc','suc_d0.eps');

clear
load('../../savings/syndata_noisy_d0.1.mat')
figure('numbertitle','off','name','Err_vs_p');
plot(p_values(:),mean(ERR_dirtyImc,2),'-.k','MarkerSize',8,'linewidth',3);
hold on;
plot(p_values(:),mean(ERR_imc,2),'-m','MarkerSize',8,'linewidth',3);
plot(p_values(:),mean(ERR_rmc,2),'-xr','MarkerSize',8,'linewidth',3);
plot(p_values(:),mean(ERR_wmc,2),'--b','MarkerSize',8,'linewidth',3);
plot(p_values(:),mean(ERR_mc, 2),'-+g','MarkerSize',8,'linewidth',3);
xlim([p_values(1) p_values(end)]);
legend('DirtyIMC','IMC','RMC','WMC','MC','Location','NorthEast');
xlabel({'$p$'},'Interpreter','latex','fontsize',12)
ylabel({'Relative Recovery Error'},'fontsize',12)

print('-depsc','err_d10e-2.eps');



clear
load('../../savings/syndata_noisy_d0.15.mat')
figure('numbertitle','off','name','Err_vs_p');
plot(p_values(:),mean(ERR_dirtyImc,2),'-.k','MarkerSize',8,'linewidth',3);
hold on;
plot(p_values(:),mean(ERR_imc,2),'-m','MarkerSize',8,'linewidth',3);
plot(p_values(:),mean(ERR_rmc,2),'-xr','MarkerSize',8,'linewidth',3);
plot(p_values(:),mean(ERR_wmc,2),'--b','MarkerSize',8,'linewidth',3);
plot(p_values(:),mean(ERR_mc, 2),'-+g','MarkerSize',8,'linewidth',3);
xlim([p_values(1) p_values(end)]);
legend('DirtyIMC','IMC','RMC','WMC','MC','Location','NorthEast');
xlabel({'$p$'},'Interpreter','latex','fontsize',12)
ylabel({'Relative Recovery Error'},'fontsize',12)
print('-depsc','err_d15e-2.eps');


clear
load('../../savings/syndata_noisy_d0.05.mat')
figure('numbertitle','off','name','Err_vs_p');
plot(p_values(:),mean(ERR_dirtyImc,2),'-.k','MarkerSize',8,'linewidth',3);
hold on;
plot(p_values(:),mean(ERR_imc,2),'-m','MarkerSize',8,'linewidth',3);
plot(p_values(:),mean(ERR_rmc,2),'-xr','MarkerSize',8,'linewidth',3);
plot(p_values(:),mean(ERR_wmc,2),'--b','MarkerSize',8,'linewidth',3);
plot(p_values(:),mean(ERR_mc, 2),'-+g','MarkerSize',8,'linewidth',3);
xlim([p_values(1) p_values(end)]);
legend('DirtyIMC','IMC','RMC','WMC','MC','Location','NorthEast');
xlabel({'$p$'},'Interpreter','latex','fontsize',12)
ylabel({'Relative Recovery Error'},'fontsize',12)
print('-depsc','err_d5e-2.eps');


clear
load('../../savings/syndata_noisy_d1e-2.mat')
figure('numbertitle','off','name','Err_vs_p');
plot(p_values(:),mean(ERR_dirtyImc,2),'-.k','MarkerSize',8,'linewidth',3);
hold on;
plot(p_values(:),mean(ERR_imc,2),'-m','MarkerSize',8,'linewidth',3);
plot(p_values(:),mean(ERR_rmc,2),'-xr','MarkerSize',8,'linewidth',3);
plot(p_values(:),mean(ERR_wmc,2),'--b','MarkerSize',8,'linewidth',3);
plot(p_values(:),mean(ERR_mc, 2),'-+g','MarkerSize',8,'linewidth',3);
xlim([p_values(1) p_values(end)]);
legend('DirtyIMC','IMC','RMC','WMC','MC','Location','NorthEast');
xlabel({'$p$'},'Interpreter','latex','fontsize',12)
ylabel({'Relative Recovery Error'},'fontsize',12)
print('-depsc','err_d1e-2.eps');