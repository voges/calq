function plot_indv_P(size, P, filter_value, metric, subposition)

markers = {'p','x','o','*','s','d','^','v','>','<','.','h','+'};
colors = colormap((parula(length(size))));
subplot(1,2,subposition);

for i = 1: length(size)
    plot(size(i),P(i),markers{i},'Color', colors(i,:),'MarkerFaceColor',colors(i,:), 'MarkerSize', 14,'LineWidth',2 );
    hold on;
end

grid;
legend('QVZ-T0','QVZ-T1', 'QVZ-T2', 'QVZ-T4', 'QVZ-T8', 'QVZ-T16', 'CALQ', 'Crumble -1','Crumble -9','Quartz');
set(gca,'FontSize',20);

xlhand = get(gca,'xlabel');
set(xlhand,'string','Average bits per quality score','fontsize',20,'Interpreter','Latex')
xlhand = get(gca,'xlabel');
set(xlhand,'fontsize',20);

ylhand = get(gca,'ylabel');
set(ylhand,'string',['Average ' metric ' difference w.r.t. original data'],'fontsize',20,'Interpreter','Latex')
ylhand = get(gca,'ylabel');
set(ylhand,'fontsize',20);

title([metric ' ($\theta=$' filter_value ')'],'fontsize',20,'Interpreter','Latex');

