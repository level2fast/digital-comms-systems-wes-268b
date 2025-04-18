function [ r_xy, lags] = c_corr( x, ref )
% C_CORR calculated the shiftd cross-correlation
% y or the ref is assumed to be shorter, conj will be taken of y

[temp, lags] = xcorr(x,ref);

L = numel(x);

r_xy = temp(L:end);

