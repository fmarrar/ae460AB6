% ae460_lab2.m ... Code for Lab 2

% clean up
clear
close all

set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',10)
set(0,'defaultFigurePosition',[360 198 460 360])

%% Import data
filename = 'data2.csv'; % file to read from
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

%%
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
Xarm = sum(0.5*(tmp(:,1:end-1)+tmp(:,2:end).*(x(1:end-1)-x(2:end))),2);
tmp = C_p.*y;
Yarm = sum(0.5*(tmp(:,1:end-1)+tmp(:,2:end).*(y(1:end-1)-y(2:end))),2);
C_mLE = Xarm + Yarm;

% lift coefficient
C_l = C_n.*cosd(param.AOA) - C_a.*sind(param.AOA);
% quarter-chord moment coefficient
C_mQC = C_mLE + 0.25*C_l;

%% Plotting
% C_p vs. x/c (linear region)
figure(1), clf
hold on
% find index of maximum AOA
max_AOA_idx = find(param.AOA == max(param.AOA));
% grab C_p distributions for increasing AOA
[indices,~] = find(param.AOA(1:max_AOA_idx) == [-4,0,4,8,12]);
for idx=indices
    plot(x, C_p(idx,:))
end
hold off
set(gca, 'ydir', 'reverse')
grid on

xlabel('\itx/c', 'FontSize', 12)
ylabel('\itC_p', 'FontSize', 12)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
    'XMinorTick', 'on', 'YMinorTick','on')
set(get(gca,'ylabel'), 'rotation', 0,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
ax = gca;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;

% C_p vs. x/c (near/at stall)
figure(2), clf
hold on
[indices,~] = find(param.AOA(1:max_AOA_idx) == [8,16,20,22,26]);
for idx=indices
    plot(x, C_p(idx,:))
end
hold off
set(gca, 'ydir', 'reverse')
grid on

xlabel('\itx/c', 'FontSize', 12)
ylabel('\itC_p', 'FontSize', 12)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
    'XMinorTick', 'on', 'YMinorTick', 'on')
set(get(gca,'ylabel'), 'rotation', 0,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
ax = gca;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;

% C_l vs. AOA
figure(3), clf
lin_bndy = find(param.AOA(1:max_AOA_idx) == 6); % sample linear region
b = [ones(lin_bndy,1) param.AOA(1:lin_bndy)] \ C_l(1:lin_bndy);
b(2) = 2*pi/rad2deg(1); % slope from thin airfoil theory
TAT = b(2)*param.AOA(1:max_AOA_idx) + b(1);
hold on
plot(param.AOA(1:max_AOA_idx), C_l(1:max_AOA_idx),...
    'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 6)
plot(param.AOA(max_AOA_idx+1:end), C_l(max_AOA_idx+1:end),...
    'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 7)
plot(param.AOA(1:max_AOA_idx), TAT)
hold off
grid on
v_padding = abs(TAT(1) - C_l(1));
axis([param.AOA(1), param.AOA(max_AOA_idx), TAT(1), max(C_l)+v_padding])
xlabel('\alpha (deg)', 'FontSize', 12)
ylabel('\itC_l', 'FontSize', 12)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02,.02],...
    'XMinorTick', 'on', 'YMinorTick', 'on')
set(get(gca,'ylabel'), 'rotation', 0,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
ax = gca;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
