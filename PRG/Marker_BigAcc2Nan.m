function [Mn, iTooBigAcc] = Marker_BigAcc2Nan(T, M, DoDebug)
% Marker_BigAcc2Nan : replace "big acceleration" with nan = invalid data
% T = time -- time series
% M = marker (X, Y, Z) -- time series
% 
% Mn = marker with nan (where nan = invalid data) 
% iTooBigAcc = index of problems (to be cleaned)
% 
% The methods relies on the idea that acc shoud be lower than a threshold
% Hypothesis 1 :
% -- The dynamics does not exeed typical human movements = 5 g
% Hypothesis 2 :
% -- The acceleration contains "outliers" defined using Tukey's method
% 
% The method also cleans valid solitary points between two areas of nan 

% Version 1 -- D. Mottet -- 2020-07-25
%    adpated from previous script Marker_BigAcc2Nan.sci -- 2018-03-31


%% determine zones where 3D acceleration is problematic

% get dt for the (approximate) time derivative
dt = median(diff(T));   % sampling period

% we get the ***raw*** acceleration using diff
% we want to amplifie all possible problems, to see them
ddX = [0 ; 0 ; diff(diff(M.X))]; % raw acceleration on X 
ddY = [0 ; 0 ; diff(diff(M.Y))]; % raw acceleration on Y
ddZ = [0 ; 0 ; diff(diff(M.Z))]; % raw acceleration on Z
ddXYZ = sqrt(ddX .* ddX + ddY .* ddY + ddZ .* ddZ);     % tangential acc
dT2   = dt*dt;                                          % dT2
a3D   = ddXYZ ./ dT2;                                   % m/s/s
% NB : a3D is the tangential acceleration, always positive
% a3D is computed as d2X/dT2, with correct units for X and T

% Option 1 : set an absolute value, using mechanics based logics
% free arm movement never exeeds 5g (if no contact) 
% free trunk movements should stay below 2g
% ThreshAcc = 20; % 10 = about 1g (g = 9.81 m/s/s) 

% Option 2 : auto threshold definition, from statistics...
% Hypothesis : what is out of "standard" acceleration is an outlier
% -- using Tukey's method (see outliers in wikipedia)
%      * outlier is over median+3*IQR
%      * use 3 * 3 to get a (very) large safety margin
%      * use 3 to stay with the standard outlier definition
% NB : using 2 will make "holes" in the series, but these holes will be
% filled by splines in a second step... Best compromise = trial&error
% ==> Theory usually helps : take 3... unless you have good reason ;-)

ThreshAcc = median(a3D) + 9 * iqr(a3D);

% as we try an automatic threshold, be prudent... limits make sense ! 
if ThreshAcc > 50 
    error ('Data contains very high acceleration : check !!')
end 
if ThreshAcc > 20
    warning('Data contains acceleration higher than 2 G: be careful...') 
end


% inform the user
disp(sprintf('ThreshAcc %0.2f, max = %0.2f', ThreshAcc, max(abs(a3D)) ))


iTooBigAcc = find (a3D > ThreshAcc);    % index of problems
% NB : find provides an index => prepend the result with 'i' for index

%% check valid solitary points between two areas of nan
% Using a high acceleration threshold usually helps keeping as much valid
% data as possible. Yet, some 'close outliers' might not be cancelled out.
% It often help to suppress cases where one single point stays alone
% betwwen two invalid zones...  

% find points that are valid, but alone between 'bad' zones 
iiSinglePoint = find(diff(iTooBigAcc) == 2);

% replace the simgle with the 'best guess' 
%   'best guess' = median of the single and two closest valid points
if ~isempty(iiSinglePoint)
    for i = 1:length(iiSinglePoint)
        % set the index of single point 
        iSingle = iTooBigAcc(iiSinglePoint(i))+1; 
        % find the index of the closest valid point before iSingle
        j = iiSinglePoint(i);
        while j > 1  & (iTooBigAcc(j) == iTooBigAcc(j-1) + 1) 
            j = j-1;
        end
        iBefore = iTooBigAcc(j)-1;
        % find the index of the closest valid point after iSingle
        j = iiSinglePoint(i)+1;
        while j < length(iTooBigAcc)  & (iTooBigAcc(j+1) == iTooBigAcc(j) + 1) 
            j = j+1;
        end
        iAfter  = iTooBigAcc(j)+1;
        % This cancels a local peak, but keeps the data in other cases. 
        if a3D(iSingle) ~= median(a3D([iBefore iSingle iAfter])) 
            iTooBigAcc = [  iTooBigAcc(1:iiSinglePoint(i))...
                            ; iSingle... 
                            ; iTooBigAcc(iiSinglePoint(i)+ 1:end)...
                            ];
            % weadded iSingle into iTooBigAcc => shift after insertion 
            iiSinglePoint(i:end) = iiSinglePoint(i:end) + 1 ; 
        end
    
    end
end

%% inform the user 
MaxContiguousBad = max(SegmentLength(iTooBigAcc));
TotalBad         = length(iTooBigAcc); 
disp(sprintf('BigAcc : %3d samples invalid (%.2f%%), max %d contiguous (%d ms)'...
    , TotalBad  ...
    , 100 * TotalBad / length(M.X)...
    , MaxContiguousBad...
    , round(MaxContiguousBad*dt*1000) )...
    )

if DoDebug
    % one figure for the acceleration in 3D
    figure();clf(); 
    I = (1:length(a3D))';
    semilogy(I, a3D, '-*b')
    hold on
    semilogy(I, ones(size(a3D))*ThreshAcc, '-.k')
    semilogy(I(iTooBigAcc), a3D(iTooBigAcc) , '*r')
    xlabel('Sample (#)')
    ylabel('Raw Acc in approx m/s/s')
    title('DEGUG: findBigAcc')
end


%% set output 
Mn = M; 
% set the "bad" data to "nan" = no useful data
if ~isempty(iTooBigAcc)
    Mn.X(iTooBigAcc) = nan;
    Mn.Y(iTooBigAcc) = nan;
    Mn.Y(iTooBigAcc) = nan;
end