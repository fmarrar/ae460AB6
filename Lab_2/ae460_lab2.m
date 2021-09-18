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
% C_p vs. x/c
figure(1), clf
grid on
hold on
max_AOA_idx = find(param.AOA == max(param.AOA));
[indices,void] = find(param.AOA(1:max_AOA_idx) == [-4,0,4,8,12]);
for idx=indices
    plot(x,C_p(idx,:))
end
set(gca, 'ydir', 'reverse')