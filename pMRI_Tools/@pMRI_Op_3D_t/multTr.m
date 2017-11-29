% Hermitian 3D Operator
function y = multTr(obj,x)

    % Upsample
    
    
    % Put input on GPU
    if strcmpi(obj.compute,'gpu_off')
        x = gpuArray(x);
    end
    
    % Use single
    if strcmpi(obj.precision,'single')
        y = zeros(obj.frame_size(1),obj.frame_size(2),obj.frame_size(3),size(obj.C,4),obj.Q,'single');
    else
        y = zeros(obj.frame_size(1),obj.frame_size(2),obj.frame_size(3),size(obj.C,4),obj.Q);
    end
    
    % Put on GPU
    if strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu')
        y = gpuArray(y);
    end
    
    y(obj.mask_patterns) = x;

    % iFFT
    y = sqrt(size(y,1)*size(y,2)*size(y,3))*ifft(ifft(ifft(y,[],1),[],2),[],3);
    
    % Coil Multiply
    fun = @(x,y) x.*y;
    y = squeeze(sum(bsxfun(fun,y,conj(obj.C)),4));

    % Make a column vector
    y = y(:);
    
    % put array back on cpu
    if strcmpi(obj.compute,'gpu_off')
        y = gather(y);
    end
end