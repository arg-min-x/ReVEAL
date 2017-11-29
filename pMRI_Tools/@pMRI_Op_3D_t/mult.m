% Forward 3D Operator
function y = mult(obj,x)

    % Put input on GPU
    if strcmpi(obj.compute,'gpu_off')
        x = gpuArray(x);
    end
    
    % Reshape
    x = reshape(x,[obj.frame_size(1),obj.frame_size(2),obj.frame_size(3),1,obj.Q]);

    % Image domain multiply and FFT
    x = repmat(x,[1,1,1,size(obj.C,4),1]);

    x = 1/sqrt(size(x,1)*size(x,2)*size(x,3))*fft(fft2(bsxfun(@times,x,obj.C)),[],3);
    
    % Downsample
    y = x(obj.mask_patterns);

    % put array back on cpu
    if strcmpi(obj.compute,'gpu_off')
        y = gather(y);
    end
end