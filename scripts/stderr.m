function stdErr = stderr(X)
% Calculate the standard error of the mean
% Usage: stdErr = stderr(X)
%
% Used by:
%       /home/Matlab/Adams_Functions/fit_IEI.m
%
% 2017-12-14 Created

stdErr = std(X)./sqrt(length(X));
