%% Import data
% clean up
clear
close all

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
C_p = DeltaP ./ param.q_WT_corr;
C_n = sum(0.5*(C_p(:,1:end-1)+C_p(:,2:end)).*(x(2:end)-x(1:end-1)),2);
C_a = sum(0.5*(C_p(:,1:end-1)+C_p(:,2:end)).*(y(1:end-1)-y(2:end)),2);
tmp = C_p.*x;
Xarm = sum(0.5*(tmp(:,1:end-1)+tmp(:,2:end).*(x(1:end-1)-x(2:end))),2);
tmp = C_p.*y;
Yarm = sum(0.5*(tmp(:,1:end-1)+tmp(:,2:end).*(y(1:end-1)-y(2:end))),2);
C_mLE = Xarm + Yarm;
