fileFolder=fullfile('../../savings/yahoo');

dirOutput=dir(fullfile(fileFolder,'*.mat'));

fileNames={dirOutput.name};


for i = fileNames
    name = i{1,1};
    disp(['../../savings/yahoo/' name]);
    load(['../../savings/yahoo/' name]);
    figure('numbertitle','off','name','AveragePrecision_vs_SampleRate');
    plot(p_values(:),mean(AP_dirtyImc,2),'-.k','MarkerSize',8,'linewidth',3);
    hold on;
    plot(p_values(:),mean(AP_imc,2),'-m', 'MarkerSize',8,'linewidth',3);
    plot(p_values(:),mean(AP_rmc,2),'-xr','MarkerSize',8,'linewidth',3);
    plot(p_values(:),mean(AP_wmc,2),'--b','MarkerSize',8,'linewidth',3);
    plot(p_values(:),mean(AP_mc, 2),'-+g','MarkerSize',8,'linewidth',3);
    xlim([p_values(1) p_values(end)]);
    ylim([0 1]);
    legend('DirtyIMC','IMC','RMC','WMC','MC','Location','SouthEast');
    xlabel({'Sampling Rate'},'Interpreter','latex','fontsize',12)
    ylabel({'Average Precision'},'fontsize',12)

    print('-depsc',[name(1:end-4),'.eps']);
    close all;
end