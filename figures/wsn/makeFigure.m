close all
load('../../savings/wsndata_0811.mat')
figure('numbertitle','off','name','Err_vs_p');
plot(p_values(:),mean(ERR_dirtyImc,2),'-.k','MarkerSize',8,'linewidth',2);
hold on;
plot(p_values(:),mean(ERR_imc,2),'-m','MarkerSize',8,'linewidth',2);
plot(p_values(:),mean(ERR_rmc,2),'-xr','MarkerSize',8,'linewidth',2);
plot(p_values(:),mean(ERR_wmc,2),'--b','MarkerSize',8,'linewidth',2);
plot(p_values(:),mean(ERR_mc, 2),'-+g','MarkerSize',8,'linewidth',2);
xlim([p_values(1) p_values(end)]);
legend('DirtyIMC','IMC','RMC','WMC','MC','Location','NorthEast');
xlabel({'$p$'},'Interpreter','latex','fontsize',12)
ylabel({'Recovery Error'},'fontsize',12)

print('-depsc','wsn.eps');





