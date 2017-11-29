function cropData(obj,crop_val)
%CROPDATA Crops the data in the image domain along the frequence encoded
%dimension to speed up computation.  Call obj.cropData(crop_val); after
%importing the data and estimation the sensitivty maps
%   Detailed explanation goes here

if obj.options.is1Dir && obj.options.isPlanar
    
    % take data to image domain
    obj.data.Yb = ifft2_shift(obj.data.Yb);
    obj.data.Yx = ifft2_shift(obj.data.Yx);

    % Crop along frequency encoded direction
    obj.data.Yb = obj.data.Yb(:,crop_val+1:end-crop_val,:,:);
    obj.data.Yx = obj.data.Yx(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampB = obj.data.sampB(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampX = obj.data.sampX(:,crop_val+1:end-crop_val,:,:);
    obj.data.maxwellCorrX = obj.data.maxwellCorrX(:,crop_val+1:end-crop_val,:,:);
    obj.data.x0 = obj.data.x0(:,crop_val+1:end-crop_val,:,:);
    
    % Crop sensitivity maps if already estimated
    if ~isempty(obj.data.sensMaps)
        obj.data.sensMaps = obj.data.sensMaps(:,crop_val+1:end-crop_val,:,:);
    end

    % convert back to k-space
    obj.data.Yb = fft2_shift(obj.data.Yb);
    obj.data.Yx = fft2_shift(obj.data.Yx);
    
elseif ~obj.options.is1Dir && obj.options.isPlanar
    % take data to image domain
    obj.data.Yb = ifft2_shift(obj.data.Yb);
    obj.data.Yx = ifft2_shift(obj.data.Yx);
    obj.data.Yy = ifft2_shift(obj.data.Yy);
    obj.data.Yz = ifft2_shift(obj.data.Yz);

    % Crop along frequency encoded direction
    obj.data.Yb = obj.data.Yb(:,crop_val+1:end-crop_val,:,:);
    obj.data.Yx = obj.data.Yx(:,crop_val+1:end-crop_val,:,:);
    obj.data.Yy = obj.data.Yy(:,crop_val+1:end-crop_val,:,:);
    obj.data.Yz = obj.data.Yz(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampB = obj.data.sampB(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampX = obj.data.sampX(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampY = obj.data.sampY(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampZ = obj.data.sampZ(:,crop_val+1:end-crop_val,:,:);
    obj.data.maxwellCorrX = obj.data.maxwellCorrX(:,crop_val+1:end-crop_val,:,:);
    obj.data.maxwellCorrY = obj.data.maxwellCorrY(:,crop_val+1:end-crop_val,:,:);
    obj.data.maxwellCorrZ = obj.data.maxwellCorrZ(:,crop_val+1:end-crop_val,:,:);
    obj.data.x0 = obj.data.x0(:,crop_val+1:end-crop_val,:,:);
    
    % Crop sensitivity maps if already estimated
    if ~isempty(obj.data.sensMaps)
        obj.data.sensMaps = obj.data.sensMaps(:,crop_val+1:end-crop_val,:,:);
    end

    % convert back to k-space
    obj.data.Yb = fft2_shift(obj.data.Yb);
    obj.data.Yx = fft2_shift(obj.data.Yx);
    obj.data.Yy = fft2_shift(obj.data.Yy);
    obj.data.Yz = fft2_shift(obj.data.Yz);
    
elseif ~obj.options.is1Dir && ~obj.options.isPlanar
    % take data to image domain
    obj.data.Yb = ifft3_shift(obj.data.Yb);
    obj.data.Yx = ifft3_shift(obj.data.Yx);
    obj.data.Yy = ifft3_shift(obj.data.Yy);
    obj.data.Yz = ifft3_shift(obj.data.Yz);

    % Crop along frequency encoded direction
    obj.data.Yb = obj.data.Yb(:,crop_val+1:end-crop_val,:,:,:);
    obj.data.Yx = obj.data.Yx(:,crop_val+1:end-crop_val,:,:,:);
    obj.data.Yy = obj.data.Yy(:,crop_val+1:end-crop_val,:,:,:);
    obj.data.Yz = obj.data.Yz(:,crop_val+1:end-crop_val,:,:,:);
    obj.data.sampB = obj.data.sampB(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampX = obj.data.sampX(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampY = obj.data.sampY(:,crop_val+1:end-crop_val,:,:);
    obj.data.sampZ = obj.data.sampZ(:,crop_val+1:end-crop_val,:,:);
    obj.data.maxwellCorrX = obj.data.maxwellCorrX(:,crop_val+1:end-crop_val,:,:);
    obj.data.maxwellCorrY = obj.data.maxwellCorrY(:,crop_val+1:end-crop_val,:,:);
    obj.data.maxwellCorrZ = obj.data.maxwellCorrZ(:,crop_val+1:end-crop_val,:,:);
    obj.data.x0 = obj.data.x0(:,crop_val+1:end-crop_val,:,:);
    
    % Crop sensitivity maps if already estimated
    if ~isempty(obj.data.sensMaps)
        obj.data.sensMaps = obj.data.sensMaps(:,crop_val+1:end-crop_val,:,:,:);
    end

    % convert back to k-space
    obj.data.Yb = fft3_shift(obj.data.Yb);
    obj.data.Yx = fft3_shift(obj.data.Yx);
    obj.data.Yy = fft3_shift(obj.data.Yy);
    obj.data.Yz = fft3_shift(obj.data.Yz);
end

end

