clear
close all

set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',10)
set(0,'defaultFigurePosition',[360 198 830 250])

% initialize plot space
t = tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'compact');
for i = 1:4
    nexttile
    if i==4
        set(gca, 'YAxisLocation', 'right')
        ylabel('$\displaystyle\frac{C_L}{C_D}$','FontSize',12,...
            'Interpreter','latex')
        set(get(gca,'ylabel'), 'rotation', 0,...
            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left')
    elseif i~=1
        set(gca, 'YColor', 'none')
    else
        ylabel('$C_L$','FontSize',12,'Interpreter','latex')
        set(get(gca,'ylabel'), 'rotation', 0,...
            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
    end
    grid on
end

files = dir('./AirfoilData/*.txt');
nf = length(files); % number of files
nRe = 3; % number of Reynolds numbers per airfoil
naf = uint16(nf/nRe); % number of airfoils
names = {'NACA 0012H', 'NACA 0015'};
for i = 1:naf
    for j = (1+3*(i-1)):(3*i)
        AFD = readtable(['./AirfoilData/' files(j).name]);
        nexttile(i)
        hold on
        plot(AFD.alpha, AFD.CL)
        title(names(i),'FontWeight','bold')
        nexttile(i+naf)
        hold on
        plot(AFD.alpha, AFD.CL./AFD.CD)
        title(names(i),'FontWeight','bold')
    end
end
xlabel(t,'$\alpha$ (deg)','FontSize',12,'Interpreter','latex')
lgd = legend('Re=100000','Re=200000','Re=500000');
lgd.Orientation = 'horizontal';
lgd.Layout.Tile = 'south';
lgd.FontSize = 9;

saveas(gcf,'AirfoilPlots','epsc')