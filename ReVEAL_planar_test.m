% Check to make sure waits are the same
% Check to make sure coil sensitvity maps are the same (fil and averaging)
% Check to make sure normalization is the same

clear
close all;
% clc;

%% Initialize class and import data
op = ReVEAL();

% Set some options
op.options.ReVEALOpts.gamma = 0.4;
op.options.ReVEALOpts.lambda0 = 3;
op.options.ReVEALOpts.sigmaSq = 0.02;
op.options.ReVEALOpts.nit =12;

op.options.ReVEALOpts.uniform_var = 0;
op.options.ReVEALOpts.L1 = 0;
op.options.ReVEALOpts.MAP = 0;
op.options.GAMPOptB.nit = 100;
op.options.GAMPOptX.nit = 100;
op.options.GAMPOptY.nit = 10;
op.options.GAMPOptZ.nit = 10;
op.options.GAMPOptB.verbose = 0;
op.options.GAMPOptX.verbose = 0;
op.options.GAMPOptY.verbose = 0;
op.options.GAMPOptZ.verbose = 0;

% Import data
op.importData('MID28ReVEAL.mat');

% Do some data checks
op.dataChecks();

% Estimate the Sensitivity MAPs
op.estimateSensMaps2D();

op.options.is1Dir = 1;
op.options.isPlanar = 1;
crop = 80;
op.cropData(crop);
op.useSingle();
% op.useGPU();

%% Reconstruct the data
op.ReVEALRecon('x');

cd /rcc4/homes/rich.178/Matlab/4D-Flow-Test
name = 'ReVEAL_LastGood';
op.saveGIFs(name);
cd /rcc4/homes/rich.178/Matlab


%%
cd /eecf/ips/data1/richa/DATA/Flow/MID28
mask = importdata('MID28_mask2.mat');
vel_lsq = importdata('MID28_vel_lsq.mat');
x = importdata('MID28_b_lsq.mat');
x = x/3.2863e-04;

mask = mask(:,crop+1:end-crop,:);
vel_lsq = vel_lsq(:,crop + 1:end-crop,:);
x = x(:,crop +1:end-crop,:);

MSE =20*log10(norm(op.outputs.xHatb(:)-x(:))/norm(x(:)));
fprintf(sprintf('\nMSE = %0.4f, Last Run was -15.3889\n',MSE))

v_hat = EstimFlowParam(mask,op.outputs.thetaX);
v_hat_MAP = EstimFlowParam(mask,op.outputs.thetaMAPX);
v_lsq = EstimFlowParam(mask,vel_lsq);

cd /rcc4/homes/rich.178/Matlab/4D-Flow-Test
f1 = figure(4);
clf
plot(v_hat.mean)
hold on
plot(v_hat_MAP.mean)
plot(v_lsq.mean)
legend('estimate','MAP','lsq')
saveas(f1,'ReVEAL_mean_LastGood.png');

f1 = figure(5);
clf
plot(v_hat.peak)
hold on
plot(v_hat_MAP.peak)
plot(v_lsq.peak)
legend('estimate','MAP','lsq')
saveas(f1,'ReVEAL_max_LastGood.png');
cd ~/Matlab
