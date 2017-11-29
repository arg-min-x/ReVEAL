function estimateSensMaps3D(obj,avg_all)
%ESTIMATESENSITIVITYMAPS Summary of this function goes here
%   Detailed explanation goes here

% fprintf('\n')
% star_display('Estimating Sensitvity Maps',0)
% tic;

if nargin < 2
    avg_all = 0;
end
avg_all
% Create time averaged image
avg_imageb = sum(obj.data.Yb,5);
avg_pattern = sum(obj.data.sampB,4);
avg_pattern(avg_pattern==0) = inf;
avg_imageb = bsxfun(@rdivide,avg_imageb,avg_pattern);
avg_imageb = ifft3_shift(avg_imageb);

avg_imagex = sum(obj.data.Yx,5);
avg_pattern = sum(obj.data.sampX,4);
avg_pattern(avg_pattern==0) = inf;
avg_imagex = bsxfun(@rdivide,avg_imagex,avg_pattern);
avg_imagex = ifft3_shift(avg_imagex);

avg_imagey = sum(obj.data.Yy,5);
avg_pattern = sum(obj.data.sampY,4);
avg_pattern(avg_pattern==0) = inf;
avg_imagey = bsxfun(@rdivide,avg_imagey,avg_pattern);
avg_imagey = ifft3_shift(avg_imagey);

avg_imagez = sum(obj.data.Yz,5);
avg_pattern = sum(obj.data.sampZ,4);
avg_pattern(avg_pattern==0) = inf;
avg_imagez = bsxfun(@rdivide,avg_imagez,avg_pattern);
avg_imagez = ifft3_shift(avg_imagez);

% Remove maxwell phase before averaging
avg_imagex = bsxfun(@times,avg_imagex,conj(obj.data.maxwellCorrX));
avg_imagey = bsxfun(@times,avg_imagey,conj(obj.data.maxwellCorrY));
avg_imagez = bsxfun(@times,avg_imagez,conj(obj.data.maxwellCorrZ));

% Use the average of all 4 maps or just the background image
if avg_all
%     avg_image = 1/4*(abs(avg_imageb) + abs(avg_imagex) + abs(avg_imagey) + abs(avg_imagez))...
%     .*exp( 1j/4*(angle(avg_imageb) + angle(avg_imagex) + angle(avg_imagey) + angle(avg_imagez)) );
    avg_image = (avg_imageb + avg_imagex + avg_imagey + avg_imagez)/4;
else
    avg_image = avg_imageb;
end

% Initialize with SOS time averaged image
obj.data.x0 = sqrt(sum(abs(avg_image).^2,4));
obj.data.x0 = repmat(obj.data.x0,[1,1,1,size(obj.data.sampB,4)]);

% Estimate sensitivity maps
p.mthd = 3; % '1' for espirit, '2' time-average spirit, '3' for walsh
p.reEst = 0; % Res-estimating sensitivities
p.fil = 3;
[obj.data.sensMaps,~] = WalshCoilCombine3D(avg_image,p);

% [obj.data.sensMaps,~] = coilSen3D(obj.data.Yb,obj.data.sampB,p);
% 
% p.mthd   = 1; %'1' espirit, '2' Walsh
% p.fil    = [6,6,1];  % size of kernal for eSPIRiT or size of filter for 'Walsh'; use [6,6,6] for eSPIRiT and [3,3,1] for Walsh
% p.oIter = 12;  % Total outer iteration
% p.oIter_2map = 4; % Outer re-weighting iterations using 2 sensitivity maps (if choose walsh, oIter_2map must be 0)
% p.oIter_1map = p.oIter - p.oIter_2map; % Outer re-weighting iterations using 1 sensitivity maps
% p.eSRT   = [8e-3, 0.95];
% p.avgPhs = 1; % Assign the phase to time-average image. 1: yes, 0: no
% p.ACSsz  = [128,96,30]; % size of the k-space block used for eSPIRiT
% p.ACSco  = [1/sqrt(2),1/sqrt(2)]; % Size of the cut-off filter used in Walsh ???????????????? use [1/2, 1/2 ]
% % p.reEst  = 0; % Res-estimating sentitivities
% [obj.data.sensMaps,~]  = coilSen_3D(obj.data.Yb, obj.data.sampB, p);
% obj.data.sensMaps = obj.data.sensMaps(:,:,:,:,2);
% display('eher')
end

