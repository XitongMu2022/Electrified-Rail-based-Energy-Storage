function PlotResult(RegDeployableRegion,Power_MultiplePeriod,RegDeployingPeriod_MultiplePeriod)
%PLOTRESULT Plot the result of the deployable ERES evaluation and the virtual energy storage aggregation model

figure;
hold on; 
set(gcf, 'Position', [100, 100, 600, 500]);  
set(gca,'position', [0.2 0.2 0.7 0.7]);

vesa = plot(RegDeployableRegion(1,:),RegDeployableRegion(2,:), '-', 'LineWidth',1);
dtp = scatter(RegDeployingPeriod_MultiplePeriod,Power_MultiplePeriod,[], 'filled');

Xlim = RegDeployingPeriod_MultiplePeriod(end);
x_ticks = 0:Xlim/3:Xlim;
xticks(x_ticks);
xlabel('Regulation deploying period (min)')

Ylim = 1.2*Power_MultiplePeriod(1);
ylim([0 Ylim]);
ylabel('Deployable ERES power (MW)')

legend([dtp,vesa],'Deployable ERES evaluation','Aggregation model');
end

