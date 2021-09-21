% ae460_lab2.m ... Code for Lab 2

% clean up
clear
close all

set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',10)
set(0,'defaultFigurePosition',[360 198 450 350])

%% Import data
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


% load comparison airfoil data
% read second airfoil
filename = 'naca0012.csv';
opts = detectImportOptions(filename, 'NumHeaderLines', 1);
opts.VariableNamesLine = 1; % row number which has variable names
foil2 = readtable(filename, opts);

% read third airfoil
filename = 'naca25112-200000.csv';
opts = detectImportOptions(filename, 'NumHeaderLines', 1);
opts.VariableNamesLine = 1;
foil3 = readtable(filename, opts);

% read fourth airfoil
filename = 'n64008a-200000.csv';
opts = detectImportOptions(filename, 'NumHeaderLines', 11);
opts.VariableNamesLine = 11;
foil4 = readtable(filename, opts);

%% Additional parameters
% in imperial units
T0 = 518.67;
R0 = 1716.554;
mu0 = 3.62E-7;
% dynamic viscosity
mu = mu0*(param.AmbTemp/T0).^1.5 .* (T0+198.72)./(param.AmbTemp+198.72);
rho = param.AmbPress*144./param.AmbTemp/R0; % freestream density
V_inf = sqrt(2*param.q_WT_corr*144./rho); % freestream velocity
Mach = V_inf./sqrt(1.4*R0*param.AmbTemp); % Mach number
Re_calc = rho.*V_inf*3.5/12./mu; % calculated Reynolds number

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
% create color, marker, and line style maps
cmap = {[0.1412, 0.5490, 0.0392]; 'blue'; [0.9490, 0.4431, 0]};
mfill = [cmap; strcat(cell(size(cmap)),'none')];
marks = {'^', 'o', 's'};
lines = {'-', '--'};
[ii,jj]=ndgrid(1:numel(marks),1:numel(lines));
% generate marker and line style pairs
style = arrayfun(@(x,y) [marks(y) lines(x)],jj(:),ii(:),'un',0);

% C_p vs. x/c (linear region)
figure(1), clf
hold on
% find index of maximum AOA
max_AOA_idx = find(param.AOA == max(param.AOA));
% grab C_p distributions for increasing AOA
selection = [-4,0,4,8,12];
[indices,~] = find(param.AOA(1:max_AOA_idx) == selection);
wrap = repmat(1:numel(cmap), 1, 2);
for idx=1:length(indices)
    plot(x, C_p(indices(idx),:), 'LineWidth', 0.75,...
        'Color', cmap{wrap(idx)}, 'LineStyle', style{idx}{2},...
        'Marker', style{idx}{1}, 'MarkerSize', 5,...
        'MarkerFaceColor', mfill{idx})
end
hold off
set(gca, 'ydir', 'reverse')
grid on

xlabel('\itx/c', 'FontSize', 12)
ylabel('\itC_p', 'FontSize', 12)

legend('\alpha=' + string(selection) + '^\circ',...
    'Location', 'NorthEast', 'FontSize', 9)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
    'XMinorTick', 'on', 'YMinorTick','on')
set(get(gca,'ylabel'), 'rotation', 0,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
ax = gca;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;

saveas(gcf,'Cp_linear','epsc')

% C_p vs. x/c (near/at stall)
figure(2), clf
hold on
selection = [8,16,18,22,26];
[indices,~] = find(param.AOA(1:max_AOA_idx) == selection);
for idx=1:length(indices)
    plot(x, C_p(indices(idx),:), 'LineWidth', 0.75,...
        'Color', cmap{wrap(idx)}, 'LineStyle', style{idx}{2},...
        'Marker', style{idx}{1}, 'MarkerSize', 5,...
        'MarkerFaceColor', mfill{idx})
end
hold off
set(gca, 'ydir', 'reverse')
grid on

xlabel('\itx/c', 'FontSize', 12)
ylabel('\itC_p', 'FontSize', 12)

legend('\alpha=' + string(selection) + '^\circ',...
    'Location', 'NorthEast', 'FontSize', 9)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
    'XMinorTick', 'on', 'YMinorTick', 'on')
set(get(gca,'ylabel'), 'rotation', 0,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
ax = gca;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;

saveas(gcf,'Cp_stall','epsc')

% C_l vs. AOA
figure(3), clf

lin_bndy = find(param.AOA(1:max_AOA_idx) == 6); % sample linear region
b = [ones(lin_bndy,1) param.AOA(1:lin_bndy)] \ C_l(1:lin_bndy);
b(2) = 2*pi/rad2deg(1); % slope from thin airfoil theory
TAT = b(2)*param.AOA(1:max_AOA_idx) + b(1); % TAT C_l vs. AOA

hold on
plot(param.AOA(1:max_AOA_idx), C_l(1:max_AOA_idx),...
    'LineStyle', 'none', 'Color', [0.1412, 0.5490, 0.0392],...
    'Marker', 'o', 'MarkerSize', 5, 'LineWidth', 0.75)
plot(param.AOA(max_AOA_idx+1:end), C_l(max_AOA_idx+1:end),...
    'LineStyle', 'none', 'Color', 'blue',...
    'Marker', 's', 'MarkerSize', 5, 'LineWidth', 0.75)
plot(param.AOA(1:max_AOA_idx), TAT, 'LineStyle','-',...
    'Color', [0.9490, 0.4431, 0], 'LineWidth', 0.75)
hold off
grid on

v_padding = 0.2*ceil(abs(TAT(1) - C_l(1))*5);
axis([param.AOA(1), param.AOA(max_AOA_idx), TAT(1), max(C_l)+v_padding])
xlabel('\alpha (deg)', 'FontSize', 12)
ylabel('\itC_l', 'FontSize', 12)

legend('{\itC_l}(\alpha\rightarrow\alpha_{max}^-)',...
    '{\itC_l}(\alpha\rightarrow\alpha_{min}^+)',...
    'Thin airfoil theory', 'Location', 'SouthEast', 'FontSize', 9)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
    'XMinorTick', 'on', 'YMinorTick', 'on')
set(get(gca,'ylabel'), 'rotation', 0,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
ax = gca;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;

saveas(gcf,'Cl','epsc')

% quarter-chord pitching moment vs. AOA
figure(4), clf
hold on
plot(param.AOA, C_mQC, 'LineStyle', 'none', 'LineWidth', 0.75,...
    'Color', 'blue', 'Marker', 's', 'MarkerSize', 5)
hold off
grid on

xlabel('\alpha (deg)', 'FontSize', 12)
ylabel('\itC_{m_{c/4}}', 'FontSize', 12)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
    'XMinorTick', 'on', 'YMinorTick', 'on')
set(get(gca,'ylabel'), 'rotation', 0,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
ax = gca;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;

saveas(gcf,'CmQC','epsc')


% ------// Comparison plots //------
Alphas = [param.AOA(:); foil2.Alpha(:); foil3.Alpha(:); foil4.Alpha(:)];

% C_l vs. AOA
figure(5), clf
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
figure(6), clf
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
xlim([min(Alphas), max(Alphas)])
legend('Clark Y14 ','NACA 0012','NACA 25112','NACA 64-008A', 'Location', 'SouthWest', 'FontSize', 9)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
    'XMinorTick', 'on', 'YMinorTick', 'on')
set(get(gca,'ylabel'), 'rotation', 0,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
ax = gca;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;

saveas(gcf,'CmQCComparison','epsc')
