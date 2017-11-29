function [ sos ] = sos_combine( x )
%SOS_COMBINE combines multicoil data using Sum of Squares
%   Detailed explanation goes here

sos = squeeze(sqrt(sum(abs(x).^2,3)));
end

