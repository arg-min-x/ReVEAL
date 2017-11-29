function dataChecks( obj )
%DATACHECKS Summary of this function goes here
%   Detailed explanation goes here

% Check if two or four encodings are given
if ~isempty(obj.data.Yb) && ~isempty(obj.data.Yx) && isempty(obj.data.Yy) && isempty(obj.data.Yz)
    display('Detecting 2 directional velocity data');
%     obj.options.is1Dir = true;
elseif ~isempty(obj.data.Yb) && ~isempty(obj.data.Yx) && ~isempty(obj.data.Yy) && ~isempty(obj.data.Yz)
    display('Detecting 3 directional velocity data');
%     obj.options.is1Dir = false;
else
    error(['If using 2 Directional data, supply Yb and Yx.  If using 3',...
        'Directional data suppy Yb, Yx, Yy, and Yz']);
end

% Check to make sure the data sizes make sense
if obj.options.is1Dir
    if ~isempty(find(size(obj.data.Yb) == size(obj.data.Yx)==0))
        error('Yb and Yx are not the same size');
    end
else
    if ~isempty(find(size(obj.data.Yb) ==size(obj.data.Yx)==0))
        error('Yb and Yx are not the same size');
    elseif ~isempty(find(size(obj.data.Yb) ==size(obj.data.Yy)==0))
        error('Yb and Yy are not the same size');
    elseif ~isempty(find(size(obj.data.Yb) ==size(obj.data.Yz)==0))
        error('Yb and Yz are not the same size');
    end
end

if ndims(obj.data.Yb) == 4
    display('Detecting planar data');
    obj.options.isPlanar = true;
elseif ndims(obj.data.Yb) ==5
    display('Detecting volumetric data');
    obj.options.isPlanar = false;
else
    error('Error the data should contain 4 or 5 dimensions'); 
end

if isempty(obj.data.maxwellCorrX)
	warning('No Maxwell correction given for x dimension')
	obj.data.maxwellCorrX = zeros(size(squeeze(obj.data.Yb(:,:,1,1))));
end
if isempty(obj.data.maxwellCorrY)
	warning('No Maxwell correction given for y dimension')
	obj.data.maxwellCorrY = zeros(size(squeeze(obj.data.Yb(:,:,1,1))));
end
if isempty(obj.data.maxwellCorrZ)
	warning('No Maxwell correction given for x dimension')
	obj.data.maxwellCorrZ = zeros(size(squeeze(obj.data.Yb(:,:,1,1))));
end
% if obj.options.Dir2
%     if ~isempty(find(size(squeeze(obj.data.Yb(:) ==size(obj.data.Yx)==0))
%         error('Yb and Yx are not the same size');
%     end
% else
%     if ~isempty(find(size(obj.data.Yb) ==size(obj.data.Yx)==0))
%         error('Yb and Yx are not the same size');
%     elseif ~isempty(find(size(obj.data.Yb) ==size(obj.data.Yy)==0))
%         error('Yb and Yy are not the same size');
%     elseif ~isempty(find(size(obj.data.Yb) ==size(obj.data.Yz)==0))
%         error('Yb and Yz are not the same size');
%     end
% end

end

