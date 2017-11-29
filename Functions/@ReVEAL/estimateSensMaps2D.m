function estimateSensMaps2D(obj)
%ESTIMATESENSITIVITYMAPS Summary of this function goes here
%   Detailed explanation goes here

% fprintf('\n')
% star_display('Estimating Sensitvity Maps',0)
% tic;

% For 1 Directional data
if obj.options.is1Dir
    warning('not using xv in sensitivity map estimation')
    % Time Average sampling pattern and k_data
    weights_b = repmat(sum(obj.data.sampB,3),[1,1,size(obj.data.Yb,3)]);
    weights_x = repmat(sum(obj.data.sampX,3),[1,1,size(obj.data.Yx,3)]);
    weights_x(find(weights_x==0)) = Inf;
    weights_b(find(weights_b==0)) = Inf;

    % Normalize k_data by weights and take an inverse Fourier Transform
    time_av_b = ifft2_shift(sum(obj.data.Yb,4)./weights_b);
    time_av_x = ifft2_shift(sum(obj.data.Yx,4)./weights_x);

    % Maxwell Correction
    if ~isempty(obj.data.maxwellCorrX)
        time_av_x = bsxfun(@times,time_av_x,exp(1j*obj.data.maxwellCorrX));
    end
	
    % Combine information from both encodings
    time_av = (time_av_b + time_av_x)/2;
%     time_av = time_av_b;

    % Calculate sensitivity maps
    [x0, walsh_maps] = WalshCoilCombine(time_av,obj.options.sensParamFil);

    % repmat to the size of the dynamic coil data
    obj.data.x0 = repmat(x0,[1,1,size(obj.data.Yb,4)]);
    obj.data.sensMaps = repmat(walsh_maps,[1,1,1,size(obj.data.Yb,4)]);
    
% For 3 Drictional data
else
    % Time Average sampling pattern and k_data
    weights_b = repmat(sum(obj.data.sampB,3),[1,1,size(obj.data.Yb,3)]);
    weights_x = repmat(sum(obj.data.sampX,3),[1,1,size(obj.data.Yx,3)]);
    weights_y = repmat(sum(obj.data.sampY,3),[1,1,size(obj.data.Yy,3)]);
    weights_z = repmat(sum(obj.data.sampZ,3),[1,1,size(obj.data.Yz,3)]);

    weights_b(find(weights_b==0)) = Inf;
    weights_x(find(weights_x==0)) = Inf;
    weights_y(find(weights_y==0)) = Inf;
    weights_z(find(weights_z==0)) = Inf;

    % Normalize k_data by weights and take an inverse Fourier Transform
    time_av_b = ifft2_shift(sum(obj.data.Yb,4)./weights_b);
    time_av_x = ifft2_shift(sum(obj.data.Yx,4)./weights_x);
    time_av_y = ifft2_shift(sum(obj.data.Yy,4)./weights_y);
    time_av_z = ifft2_shift(sum(obj.data.Yz,4)./weights_z);

    % Maxwell Correction
    if ~isempty(obj.data.maxwellCorrX)
        time_av_x = bsxfun(@times,time_av_x,exp(1j*obj.data.maxwellCorrX));
    elseif ~isempty(obj.data.maxwellCorrY)
        time_av_y = bsxfun(@times,time_av_y,exp(1j*obj.data.maxwellCorrY));
    elseif ~isempty(obj.data.maxwellCorrZ)
        time_av_z = bsxfun(@times,time_av_z,exp(1j*obj.data.maxwellCorrZ));
    end
	
    % Combine information from both encodings
    %time_av = (time_av_b + time_av_x + time_av_y + time_av_z)/4;
%     time_av = abs((time_av_b + time_av_x + time_av_y + time_av_z)/4).*exp(1j*angle(time_av_b));
    time_av = time_av_b;

    % Calculate sensitivity maps
    [x0, walsh_maps] = WalshCoilCombine(time_av,obj.options.sensParamFil);

    % repmat to the size of the dynamic coil data
    obj.data.x0 = repmat(x0,[1,1,size(obj.data.Yb,4)]);
    obj.data.sensMaps = repmat(walsh_maps,[1,1,1,size(obj.data.Yb,4)]);
end

end

