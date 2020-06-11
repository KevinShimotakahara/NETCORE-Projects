
%%
function PlotNetPlan(NetworkPlot)

global eNBs
global UDs
global nD2D
global nUEs

global simParameters

if(NetworkPlot)
    x_enb = [];
    y_enb = [];
    r_enb = [];
    
    for e = 1:length(eNBs)
        x_enb = [x_enb; eNBs(e).nodePosition(1)];
        y_enb = [y_enb; eNBs(e).nodePosition(2)];
        r_enb = [r_enb; eNBs(e).cellRadius];
    end
    
    x_ud = [];
    y_ud = [];
    for u_ = 1:length(UDs)
        x_ud = [x_ud; UDs(u_).nodePosition(1)];
        y_ud = [y_ud; UDs(u_).nodePosition(2)];
    end
    
    figure; hold on; grid on;
%     xlim([0 simParameters.GridSize*1000]);
%     ylim([0 simParameters.GridSize*1000]);

    h1 = scatter(x_enb, y_enb, 'filled', 'd', 'MarkerEdgeColor', 'r',...
        'MarkerFaceColor', 'r');
    
    h2 = scatter(x_ud, y_ud, 'filled', 'x', 'MarkerEdgeColor', 'b',...
        'MarkerFaceColor', 'b');
    
    for k=1:length(x_enb)
        text(x_enb(k), y_enb(k), num2str(eNBs(k).nodeID), 'FontSize', 12, 'FontWeight', 'Bold');
    end
    
    for k=1:length(x_ud)
        text(x_ud(k), y_ud(k), num2str(UDs(k).nodeID), 'FontSize', 12, 'FontWeight', 'Bold');
    end
    
    x = [x_enb];
    y = [y_enb];
    r = [r_enb];
    
    for i = 1:length(x)
        th = 0:pi/50:2*pi;
        xunit = r(i) * cos(th) + x(i);
        yunit = r(i) * sin(th) + y(i);
        plot(xunit, yunit, 'k', 'LineWidth', 1);
    end
    
    legend([h1 h2], {'eNB', 'UDs'}, 'Location','best');
    title({strcat('Nodes deployment; nD2Ds = ', num2str(nD2D),...
        '; nUEs = ', num2str(nUEs)), strcat('eNB radius = ', num2str(eNBs(1).cellRadius))});
    
    ax = gca;
    ax.FontSize = 13;
    ax.FontWeight = 'Bold';
    set(h1, 'LineWidth', 2);
    set(gcf, 'Position', [400, 200, 1000, 700]);
end

end


