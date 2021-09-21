% ae460_lab2_Comparison.m ... Airfoil comparisons for Lab 2

% clean up
clear
close all

set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',10)
set(0,'defaultFigurePosition',[360 198 450 350])

%% Import data
% experimental data file
filename = 'catch_22.csv'; % file to read from
% environment
opts = detectImportOptions(filename, 'Range', 'A:G'); % select parameters
opts.DataLines = 6; % data start line
varNames = {'DataPoint','AOA','q_WT_corr',...
    'q_corr_fac','AmbTemp','AmbPress','Re'}; % set valid variable names
param = readtable(filename, opts); % read data using options
param.Properties.VariableNames = varNames; % apply variable names
% airfoil
opts = detectImportOptions(filename, 'Range', 'J4'); % x/c, y/c, pressures
dist = readtable(filename, opts);
% extract variables as arrays
x = dist{1,:}; % x/c
y = dist{2,:}; % y/c
DeltaP = dist{3:end,:}; % differential pressures

% second airfoil read
filename = 'naca0012.csv';
opts = detectImportOptions(filename, 'NumHeaderLines', 1);
opts.VariableNamesLine = 1; % row number which has variable names
foil2 = readtable(filename, opts);

% third airfoil read
filename = 'naca25112-200000.csv';
opts = detectImportOptions(filename, 'NumHeaderLines', 1);
opts.VariableNamesLine = 1;
foil3 = readtable(filename, opts);

% fourth airfoil read
filename = 'n64008a-200000.csv';
opts = detectImportOptions(filename, 'NumHeaderLines', 11);
opts.VariableNamesLine = 11;
foil4 = readtable(filename, opts);

%% Calculate coefficients
% pressure coefficient
C_p = DeltaP ./ param.q_WT_corr;
% normal force coefficient
C_n = sum(0.5*(C_p(:,1:end-1)+C_p(:,2:end)).*(x(2:end)-x(1:end-1)),2);
% axial force coefficient
C_a = sum(0.5*(C_p(:,1:end-1)+C_p(:,2:end)).*(y(1:end-1)-y(2:end)),2);

% leading-edge moment coefficient
% calculate x and y moment arm contributions seperately
% sum the contributions
tmp = C_p.*x;
Xarm = sum(0.5*(tmp(:,1:end-1)+tmp(:,2:end)).*(x(1:end-1)-x(2:end)),2);
tmp = C_p.*y;
Yarm = sum(0.5*(tmp(:,1:end-1)+tmp(:,2:end)).*(y(1:end-1)-y(2:end)),2);
C_mLE = Xarm + Yarm;

% lift coefficient
C_l = C_n.*cosd(param.AOA) - C_a.*sind(param.AOA);
% drag coefficient
C_d = C_n.*sind(param.AOA) + C_a.*cosd(param.AOA);
% quarter-chord moment coefficient
C_mQC = C_mLE + 0.25*C_l;

%% Plotting
% find index of maximum AOA
max_AOA_idx = find(param.AOA == max(param.AOA));

% C_l vs. AOA
figure(3), clf
hold on
plot(param.AOA(1:max_AOA_idx), C_l(1:max_AOA_idx),...
    'LineStyle', 'none', 'Color', [0.1412, 0.5490, 0.0392],...
    'Marker', 'o', 'MarkerSize', 5, 'LineWidth', 0.75)
% plot 2nd, 3rd, 4th Cl data
plot(foil2.Alpha, foil2.Cl,...
    'LineStyle', 'none', 'Color', [0.9063 0.2891 0.1523],...
    'Marker', '^', 'MarkerSize', 5, 'LineWidth', 0.75)
plot(foil3.Alpha, foil3.Cl,...
    'LineStyle', 'none', 'Color', [0.4940 0.1840 0.5560],...
    'Marker', 'h', 'MarkerSize', 5, 'LineWidth', 0.75)
plot(foil4.Alpha, foil4.Cl,...
    'LineStyle', 'none', 'Color', [0.6350 0.0780 0.1840],...
    'Marker', '.', 'MarkerSize', 10, 'LineWidth', 0.75)
hold off
grid on

xlim([param.AOA(1), param.AOA(max_AOA_idx)])
xlabel('\alpha (deg)', 'FontSize', 12)
ylabel('\itC_l', 'FontSize', 12)

legend('Clark Y14', 'NACA 0012', 'NACA 25112', 'NACA 64-008A',...
    'Location', 'SouthEast', 'FontSize', 9)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
    'XMinorTick', 'on', 'YMinorTick', 'on')
set(get(gca,'ylabel'), 'rotation', 0,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
ax = gca;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;

saveas(gcf,'ClComparison','epsc')

% quarter-chord pitching moment vs. AOA
figure(4), clf
hold on
plot(param.AOA, C_mQC, 'LineStyle', 'none', 'Color', [0.1412, 0.5490, 0.0392],...
    'Marker', 'o', 'MarkerSize', 5, 'LineWidth', 0.75)
% plot 2nd, 3rd, 4th Cm data
plot(foil2.Alpha, foil2.Cm,...
    'LineStyle', 'none', 'Color', [0.9063 0.2891 0.1523],...
    'Marker', '^', 'MarkerSize', 3, 'LineWidth', 0.75)
plot(foil3.Alpha, foil3.Cm,...
    'LineStyle', 'none', 'Color', [0.4940 0.1840 0.5560],...
    'Marker', 'h', 'MarkerSize', 3, 'LineWidth', 0.75)
plot(foil4.Alpha, foil4.Cm,...
    'LineStyle', 'none', 'Color', [0.6350 0.0780 0.1840],...
    'Marker', '.', 'MarkerSize', 10, 'LineWidth', 0.75)

hold off
grid on

xlabel('\alpha (deg)', 'FontSize', 12)
ylabel('\itC_{m_{c/4}}', 'FontSize', 12)
xlim([-16 26])
legend('Clark Y14 ','NACA 0012','NACA 25112','NACA 64-008A', 'Location', 'SouthWest', 'FontSize', 9)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
    'XMinorTick', 'on', 'YMinorTick', 'on')
set(get(gca,'ylabel'), 'rotation', 0,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
ax = gca;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;

saveas(gcf,'CmQCComparison','epsc')
