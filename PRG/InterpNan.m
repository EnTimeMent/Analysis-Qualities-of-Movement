function [Xi, iNan] = InterpNan(T, X)
% InterpNan replaces Nan in the time series x using spline interpolations
%      X should be a vector (typically a time series)
%      Xi is the series after interpolation (same size as X)
%      iNan is the index of holes in X (nan)

% Version 1 -- D. Mottet -- 2016-05-29

%% Input check
if min(size(X)) > 1
    error('Only vectors can be interpolated')
end
if max(size(X)) < 3
    error('At least 3 elements are needed to interpolate')
end


%% initialisation to 'nothing to do'
Xi   = X ;      % no transformation

%% check for nan to be replaced
iHole = find( isnan(X));        % hole = nan

if ~isempty(iHole)
    iGood = find(~isnan(X));    % not hole = not nan
    % interpolate holes, using splines on the good part
    % and linear extrapolation at the beg-end
    Xinterp = spline(T(iGood), X(iGood), T(iHole));
    % Use Xinterp to fill Xi at the nan
    Xi(iHole) = Xinterp;
end

%% set output 
iNan = iHole;

end
