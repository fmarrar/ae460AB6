%%
clear
close all

opts = detectImportOptions('data1.csv','NumHeaderLines',2); % number of header lines which are to be ignored
opts.VariableNamesLine = 3; % row number which has variable names
opts.VariableNames = regexprep(opts.VariableNames,'_(\w*)',''); % clear units from variable names
opts.DataLine = 4; % row number from which the actual data starts
tbl = readtable('data1.csv',opts);

%% Rake Calculations
[nrows,ncols] = size(tbl);
TF = contains(tbl.Properties.VariableNames, 'RakeT');
NT = nnz(TF); % number of total pressure rakes
TF = contains(tbl.Properties.VariableNames, 'RakeS');
NS = nnz(TF); % number of static pressure rakes
RakeTAcc = zeros(nrows,1); % accumulating variable
RakeSAcc = zeros(nrows,1);

for idx = 1:NT
    RakeTAcc = RakeTAcc + tbl.(['RakeT' num2str(idx)]);
end

for idx = 1:NS
    RakeSAcc = RakeSAcc + tbl.(['RakeS' num2str(idx)]);
end

P0avg = RakeTAcc/NT; % average total pressures
PSavg = RakeSAcc/NS; % average static pressures
q = P0avg - PSavg; % average dynamic pressures

TPV = zeros(nrows,NT);
for idx = 1:NT
    TPV(:,idx) = 100*(tbl.(['RakeT' num2str(idx)]) - P0avg) ./ q;
end

SPV = zeros(nrows,NS);
for idx = 1:NS
    SPV(:,idx) = 100*(tbl.(['RakeS' num2str(idx)]) - PSavg) ./ q;
end

%% Plotting
figID = 1;
figure(figID), clf
hold on
C = colorSpectrum(NT);
all_marks = {'o','+','*','x','s','d','^','v','>','<','p','h'};
h = struct();
for idx = 1:NT
    h.(['probe_' num2str(idx)]) = line(q, TPV(:,idx), 'Color', C(idx,:));
    set(h.(['probe_' num2str(idx)]), 'LineStyle', 'none', ...
        'Marker', all_marks{mod(idx,12)}, 'MarkerSize', 6)
end
hold off

% Add labels
hXLabel = xlabel('Actual Dynamic Pressure (psi)');
hYLabel = ylabel('Percent Total Pressure Variation');

% Add legend
fields = fieldnames(h);
labels = regexprep(fields,'_',' ');
axes = struct2cell(h);
hLegend = legend([axes{:}], labels{:}, 'Location', 'eastoutside');

% Adjust font
set(gca, 'FontName', 'Helvetica')
% set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hLegend, gca], 'FontSize', 9)
set([hXLabel, hYLabel], 'FontSize', 11)

% Adjust axes properties
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'GridLineStyle','--', ...
    'LineWidth', 1)

set(gcf, 'Position',  [100, 100, 610, 400])

figID = figID + 1;
figure(figID), clf
hold on
C = colorSpectrum(NS);
all_marks = {'o','+','*','x','s','d','^','v','>','<','p','h'};
h = struct();
for idx = 1:NS
    h.(['probe_' num2str(idx)]) = line(q, SPV(:,idx), 'Color', C(idx,:));
    set(h.(['probe_' num2str(idx)]), 'LineStyle', 'none', ...
        'Marker', all_marks{mod(idx,12)}, 'MarkerSize', 6)
end
hold off

% Add labels
hXLabel = xlabel('Actual Dynamic Pressure (psi)');
hYLabel = ylabel('Percent Static Pressure Variation');

% Add legend
fields = fieldnames(h);
labels = regexprep(fields,'_',' ');
axes = struct2cell(h);
hLegend = legend([axes{:}], labels{:}, 'Location', 'eastoutside');

% Adjust font
set(gca, 'FontName', 'Helvetica')
% set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hLegend, gca], 'FontSize', 9)
set([hXLabel, hYLabel], 'FontSize', 11)

% Adjust axes properties
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'GridLineStyle','--', ...
    'LineWidth', 1)

set(gcf, 'Position',  [100, 100, 610, 400])

%% Measurement Differences
% percent differences
OTDP_pd = abs(tbl.OmegaTransmitterDeltaP - q) ./ q * 100;
MDP_pd = abs(tbl.ManometerDeltaP - q) ./ q * 100;

% ---- Plot calibration of OTDP ----

figID = figID + 1;
figure(figID), clf
hold on
hOTDP = line(tbl.OmegaTransmitterDeltaP, q);
set(hOTDP, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 5, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1])
b = [ones(nrows,1) tbl.OmegaTransmitterDeltaP] \ q;
OTDPfit = line(tbl.OmegaTransmitterDeltaP, ...
    b(1)+b(2)*tbl.OmegaTransmitterDeltaP, ...
    'Color', [0.85 0.35 0.01]);
% b1 = tbl.OmegaTransmitterDeltaP \ q;
% OTDPfit = line([0; tbl.OmegaTransmitterDeltaP], b1*[0; tbl.OmegaTransmitterDeltaP], 'Color', 'r');
hold off
axlim = [min(tbl.OmegaTransmitterDeltaP(1), q(1)), ...
    max(tbl.OmegaTransmitterDeltaP(end), q(end))];
axis([axlim, axlim])

hXLabel = xlabel('Omega Transmitter Delta P (psi)');
hYLabel = ylabel('Actual Dynamic Pressure (psi)');

% Add fit equation
hText = text(0.0275, 0.0225, ...
    sprintf('{\\itq = %0.5g + %0.5g\\DeltaP_\\Omega}', b(1), b(2)));
% hText = text(0.0275, 0.0225, ...
%     sprintf('{\\itq = %0.5g\\DeltaP_\\Omega}', b1));

hLegend = legend([hOTDP, OTDPfit], 'Data ({\it\DeltaP_\Omega})', ...
    'Fit ({\it\beta_1x})', 'Location', 'NorthWest');

set(gca, 'FontName', 'Helvetica')
set([hLegend, gca], 'FontSize', 9)
set([hXLabel, hYLabel, hText], 'FontSize', 11)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'XColor', [.3 .3 .3], ...
    'YColor', [.3 .3 .3], 'LineWidth', 1)

set(gcf, 'Position',  [100, 100, 520, 400])

%%
R = 1716; % ft-lbf/slug-R
% 1 psi = 144 psf
rho = tbl.AmbientPress*144 / R ./ tbl.AmbientTemp;
U = sqrt(2*q*144./rho);
Repft = rho .* U ./ tbl.Viscosity;
Mach = U ./ sqrt(1.4*R*tbl.AmbientTemp);

% ---- Plots against RPM ----

figID = figID + 1;
figure(figID), clf
subplot(2,2,1)
line(tbl.MotorSpeed, q, 'LineStyle', 'none', ...
    'Marker', 'o', 'MarkerSize', 5, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1])
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1, 'GridLineStyle', '--')
set(gca, 'FontName', 'Helvetica', 'FontSize', 8)
xlabel(['\fontname{Helvetica}\fontsize{9}Motor Speed (RPM)',newline,...
    '\fontname{Times New Roman}\fontsize{12}\color{black}(a) {\itq} vs. RPM'])
ylabel('Dynamic Pressure (psi)', 'FontSize', 9)

subplot(2,2,2)
line(tbl.MotorSpeed, U, 'LineStyle', 'none', ...
    'Marker', 'o', 'MarkerSize', 5, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1])
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1, 'GridLineStyle','--')
set(gca, 'FontName', 'Helvetica', 'FontSize', 8)
xlabel(['\fontname{Helvetica}\fontsize{9}Motor Speed (RPM)',newline,...
    '\fontname{Times New Roman}\fontsize{12}\color{black}(b) Velocity vs. RPM'])
ylabel('Velocity (ft/s)', 'FontSize', 9)

subplot(2,2,3)
line(tbl.MotorSpeed, Repft, 'LineStyle', 'none', ...
    'Marker', 'o', 'MarkerSize', 5, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1])
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1, 'GridLineStyle','--')
set(gca, 'FontName', 'Helvetica', 'FontSize', 8)
xlabel(['\fontname{Helvetica}\fontsize{9}Motor Speed (RPM)',newline,...
    '\fontname{Times New Roman}\fontsize{12}\color{black}(c) Reynolds number per unit length vs. RPM'])
ylabel('Re/L (ft^{-1})', 'FontSize', 9)

subplot(2,2,4)
line(tbl.MotorSpeed, Mach, 'LineStyle', 'none', ...
    'Marker', 'o', 'MarkerSize', 5, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1])
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1, 'GridLineStyle','--')
set(gca, 'FontName', 'Helvetica', 'FontSize', 8)
xlabel(['\fontname{Helvetica}\fontsize{9}Motor Speed (RPM)',newline,...
    '\fontname{Times New Roman}\fontsize{12}\color{black}(d) Mach number vs. RPM'])
ylabel('Mach Number', 'FontSize', 9)

set(gcf, 'Position',  [100, 100, 610, 420])

%% Velocity Uncertainty
UnPamb = 0.148/2.036;
UnTamb = 1.8;
UnPi = 0.1/100; % pressure probe uncertainty (psi)
UnP0 = sqrt(NT*(UnPi/NT)^2); % total pressure uncertainty
UnPS = sqrt(NS*(UnPi/NS)^2); % static pressure uncertainty

syms RHO P0 PS PAMB TAMB

Uexp   = sqrt(2*(P0 - PS)/RHO); % Bernoulli
rhoexp = PAMB/R/TAMB; % equation of state

% calculate partial derivatives
dUdP0  = double(subs(diff(Uexp, P0), {RHO P0 PS}, {rho P0avg PSavg}));
dUdPS  = double(subs(diff(Uexp, PS), {RHO P0 PS}, {rho P0avg PSavg}));
dUdrho = double(subs(diff(Uexp, RHO), {RHO P0 PS}, {rho P0avg PSavg}));
drhodP = double(subs(diff(rhoexp, PAMB), {PAMB TAMB}, ...
    {tbl.AmbientPress tbl.AmbientTemp}));
drhodT = double(subs(diff(rhoexp, TAMB), {PAMB TAMB}, ...
    {tbl.AmbientPress tbl.AmbientTemp}));

Unrho = sqrt((drhodP.*UnPamb).^2 + (drhodT.*UnTamb).^2);
UnU = sqrt((dUdP0.*UnP0).^2 + (dUdPS.*UnPS).^2 + (dUdrho.*Unrho).^2);
PercentUn = UnU ./ U * 100;
% Section 5.6
tbl5p6 = table(tbl.MotorSpeed, U, UnU, PercentUn);