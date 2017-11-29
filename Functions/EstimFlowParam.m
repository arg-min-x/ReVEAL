function [v] = EstimFlowParam(mask,velocity)
%ESTIMFLOWPARAM Summary of this function goes here
%   Detailed explanation goes here

% Create a function handle for bsxfun
point_by_point = @(x,y) x.*y;

num_points = round(0.1*length(find(mask(:,:,1)==1)));

% Pull out data from a ROI
v_mask = bsxfun(point_by_point,velocity,mask);
for ind = 1:size(velocity,3)
   tmp = v_mask(:,:,ind);
   tmp = tmp(find(tmp~=0));
   [sorted,inds] = sort(abs(tmp),'descend');
   sorted = sign(tmp(inds)).*sorted;
   [~,peak_ind] = max(abs(tmp));
   if ~isempty(peak_ind)
       if length(sorted>=num_points)
%             v.peak_avg(ind) = mean(sorted(1:num_points));
            v.peak_avg = [];
       else
           v.peak_avg = [];
       end
        v.peak(ind) = tmp(peak_ind);
        v.mean(ind) = mean(tmp);
        v.flow(ind) = sum(tmp);
        if ind>1
            v.sum(ind) = v.sum(ind-1) + sum(tmp);
        else
            v.sum(ind) = sum(tmp);
        end
        v.var(ind) = var(tmp);
   else
        v.peak_avg(ind) = 0;
        v.peak(ind) = 0;
        v.mean(ind) = 0;
        v.flow(ind) = 0;
        v.sum(ind) = 0;
        v.var(ind) = 0;
   end

   v.flow(ind) = sum(tmp);
end

v.total_flow = sum(v.flow);
% v_lsq.total_flow = sum(v_lsq.flow);

end

