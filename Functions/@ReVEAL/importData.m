function importData(obj, fileName)
%IMPORTDATA This function imports the matlab structure containing the PC-MRI
% data
% 	Inputs: fileName - The name of the file that contains the PC-MRI data.  
%			 The file should be "fileName.mat" in the Matlab path
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Last update:  2/18/2015
%**************************************************************************

% Import Data structure
if exist(fileName)
    data = importdata(fileName);

    %% Copy data to object;
    obj.data.Yb = double(data.Yb);
    obj.data.Yx = double(data.Yx);
    obj.data.Yy = double(data.Yy);
    obj.data.Yz = double(data.Yz);
    obj.data.sampB = data.sampB;
    obj.data.sampX = data.sampX;
    obj.data.sampY = data.sampY;
    obj.data.sampZ = data.sampZ;
    obj.data.maxwellCorrX = data.maxwellCorrX;
    obj.data.maxwellCorrY = data.maxwellCorrY;
    obj.data.maxwellCorrZ = data.maxwellCorrZ;
    obj.data.noiseChannel = data.noiseChannel;
elseif exist([fileName(1:end-4),'B.mat'])
    %% Copy data to object;
    obj.data.Yb = importdata([fileName(1:end-4),'B.mat']);
    obj.data.Yx = importdata([fileName(1:end-4),'X.mat']);
    obj.data.Yy = importdata([fileName(1:end-4),'Y.mat']);
    obj.data.Yz = importdata([fileName(1:end-4),'Z.mat']);
    obj.data.sampB = importdata([fileName(1:end-8),'sampB.mat']);
    obj.data.sampX = importdata([fileName(1:end-8),'sampX.mat']);
    obj.data.sampY = importdata([fileName(1:end-8),'sampY.mat']);
    obj.data.sampZ = importdata([fileName(1:end-8),'sampZ.mat']);
    obj.data.maxwellCorrX = [];
    obj.data.maxwellCorrY = [];
    obj.data.maxwellCorrZ = [];
    if exist([fileName(1:end-8),'noiseChannel.mat'])
        obj.data.noiseChannel = importdata([fileName(1:end-8),'noiseChannel.mat']);
    end
%     obj.data.noiseChannel = [];
else
    error('Cannot find Data');
end

% Check for 2D cine data
if ndims(obj.data.Yb) == 4
    
    % 1 Directional Data
    if isempty(obj.data.Yz) && isempty(obj.data.Yy)
        if ~isempty(obj.data.Yb) && ~isempty(obj.data.Yx)
            display('1 Directional data detected');

            if isempty(obj.options.is1Dir)
                obj.options.is1Dir = 1;
            end

            % Get Noise Channel if no noise is provided
            if isempty(obj.options.ReVEALOpts.wvar)
                if isempty(obj.data.noiseChannel)
                    tmp = obj.data.Yb;
                    [~,obj.data.noiseChannel] = coilCombine_old(tmp,obj.options.ReVEALOpts.num_coils);
                    obj.options.ReVEALOpts.wvar = estimNoiseVar(obj.data.noiseChannel,obj.data.sampB)
                else
                     obj.options.ReVEALOpts.wvar = var(obj.data.noiseChannel);
                end
            end    

            % Compress coils jointly across endocings
            kdata_coil = cat(5,obj.data.Yb,obj.data.Yx);
            kdata_coil = coilCombine( kdata_coil, obj.options.ReVEALOpts.num_coils,'2dt');
            obj.data.Yb = kdata_coil(:,:,:,:,1);
            obj.data.Yx = kdata_coil(:,:,:,:,2);

    %         [obj.data.Yb, obj.data.noiseChannel]= coilCombine_old( obj.data.Yb, obj.options.ReVEALOpts.num_coils);
    %         [obj.data.Yx]= coilCombine_old( obj.data.Yx, obj.options.ReVEALOpts.num_coils);
        else
            error('Yb and Yx data are requiured')
        end
    obj.options.is1Dir = 1;
    obj.options.isPlanar = 1;
    % 3 Directional Data
    else
        obj.options.is1Dir = 0;
        obj.options.isPlanar = 1;
        
        % Compress coils jointly across endocings
        kdata_coil = cat(5,obj.data.Yb,obj.data.Yx,obj.data.Yy,obj.data.Yz);
        kdata_coil = coilCombine( kdata_coil, obj.options.ReVEALOpts.num_coils,'2dt');
        obj.data.Yb = kdata_coil(:,:,:,:,1);
        obj.data.Yx = kdata_coil(:,:,:,:,2);  
        obj.data.Yy = kdata_coil(:,:,:,:,3); 
        obj.data.Yz = kdata_coil(:,:,:,:,4); 

        % Get Noise Variance
        if ~isempty(obj.data.noiseChannel)
            obj.options.ReVEALOpts.wvar = var(obj.data.noiseChannel(:));
        else
            error('No Noise Scan provided');
        end
    end

    % Downsample the data
    fun = @(x,y) x.*y;
    obj.data.Yb = bsxfun(fun,obj.data.Yb,permute(obj.data.sampB,[1,2,4,3]));
    obj.data.Yx = bsxfun(fun,obj.data.Yx,permute(obj.data.sampX,[1,2,4,3]));
    obj.data.Yy = bsxfun(fun,obj.data.Yy,permute(obj.data.sampY,[1,2,4,3]));
    obj.data.Yz = bsxfun(fun,obj.data.Yz,permute(obj.data.sampZ,[1,2,4,3]));

    % Time Average sampling pattern and k_data
    weights_b = repmat(sum(cat(3,obj.data.sampB,obj.data.sampX),3),[1,1,size(obj.data.Yb,3)]);
    weights_b(weights_b==0) = Inf;

    % Normalize k_data by weights and take an inverse Fourier Transform
    time_av_b = ifft2_shift(sum(cat(4,obj.data.Yb,obj.data.Yx),4)./weights_b);
    % time_av_b = sos_combine(time_av_b);
    time_av_b = WalshCoilCombine(time_av_b,3);

    % Normalize Data
    x_60 = prctile(abs(time_av_b(:)),60);
    tmp = time_av_b(abs(time_av_b)>=x_60);
    scale = mean(abs(tmp));

    obj.data.x0 = obj.data.x0/scale;
    obj.data.Yb = obj.data.Yb/scale;
    obj.data.Yx = obj.data.Yx/scale;
    obj.data.Yy = obj.data.Yy/scale;
    obj.data.Yz = obj.data.Yz/scale;
    obj.options.ReVEALOpts.wvar = obj.options.ReVEALOpts.wvar/scale^2;

    % Set sigmaSq if none is provided
    if isempty(obj.options.ReVEALOpts.sigmaSq)
        obj.options.ReVEALOpts.sigmaSq = obj.options.ReVEALOpts.wvar/sqrt(size(obj.data.Yb,3));
    end

% 4D data
elseif ndims(obj.data.Yb) == 5
    display('4D detected')
    obj.options.isPlanar = 0;
    obj.options.is1Dir = 0;
    
    % Coil combine if needed
    if size(obj.data.Yb,4) > obj.options.ReVEALOpts.num_coils
        display('coil combining')
        kdata = cat(6,obj.data.Yb,obj.data.Yx,obj.data.Yy,obj.data.Yz);
        kdata = coilCombine(kdata,obj.options.ReVEALOpts.num_coils,'3dt');
        obj.data.Yb = kdata(:,:,:,:,:,1);
        obj.data.Yx = kdata(:,:,:,:,:,2);  
        obj.data.Yy = kdata(:,:,:,:,:,3); 
        obj.data.Yz = kdata(:,:,:,:,:,4);
    end
    
    % Downsample Data
    if (size(obj.data.Yb,1) ~= size(obj.data.sampB,1)) ||...
       (size(obj.data.Yb,2) ~= size(obj.data.sampB,2)) ||...
       (size(obj.data.Yb,3) ~= size(obj.data.sampB,3)) ||...
       (size(obj.data.Yb,5) ~= size(obj.data.sampB,4))
        error('Data and Sampling patterns are in disagreement')
    end
    obj.data.Yb = bsxfun(@times,obj.data.Yb,permute(obj.data.sampB,[1,2,3,5,4]));
    obj.data.Yx = bsxfun(@times,obj.data.Yx,permute(obj.data.sampX,[1,2,3,5,4]));
    obj.data.Yy = bsxfun(@times,obj.data.Yy,permute(obj.data.sampY,[1,2,3,5,4]));
    obj.data.Yz = bsxfun(@times,obj.data.Yz,permute(obj.data.sampZ,[1,2,3,5,4]));
    
    % Estimate Noise Variance
    if ~isempty(obj.data.noiseChannel)
        obj.options.ReVEALOpts.wvar = var(obj.data.noiseChannel(:));
    end
    
    % Normalize the data
    scale = max(abs(obj.data.Yb(:)))/100;
%     obj.data.x0 = obj.data.x0/scale;
    obj.data.Yb = obj.data.Yb/scale;
    obj.data.Yx = obj.data.Yx/scale;
    obj.data.Yy = obj.data.Yy/scale;
    obj.data.Yz = obj.data.Yz/scale;
    obj.options.ReVEALOpts.wvar = obj.options.ReVEALOpts.wvar/scale^2;
    
    % Set sigmaSq if none is provided
    if isempty(obj.options.ReVEALOpts.sigmaSq)
        obj.options.ReVEALOpts.sigmaSq = obj.options.ReVEALOpts.wvar/sqrt(size(obj.data.Yb,3));
    end
    
    %% Setup Maxwell Corrections
    if isempty(obj.data.maxwellCorrX)
        obj.data.maxwellCorrX = ones(size(squeeze(obj.data.sampB(:,:,:,1))));
        obj.data.maxwellCorrY = ones(size(squeeze(obj.data.sampB(:,:,:,1))));
        obj.data.maxwellCorrZ = ones(size(squeeze(obj.data.sampB(:,:,:,1))));
    end
else
    error('Input data should be 4 or five dimensional')  
end

fprintf(sprintf('Estimated Noise Variance = %s\n',num2str(obj.options.ReVEALOpts.wvar)));
fprintf(sprintf('Estimated sigmaSq = %s\n',num2str(obj.options.ReVEALOpts.sigmaSq))); 

end

