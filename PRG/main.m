% main : entry point
%
%
%   Author(s):
%       D. Mottet, 2020-june-26, Version 0.1

%if ~exist('head_x')
if ~exist('x')
    clear all
    close all
    importCSV;
end

DoDebug = 1;        % 1 = debug actions on
close all           % close all figures

%% import data
if ~exist('rescaled')
    % rescale to standard units
    Pos2meter = 1000;   % likely in mm
    Tim2second = 1;     % already in seconds
    Time   = Time  ./ Tim2second;
    x = x ./ Pos2meter; %%changed from head_x/y/z to plain x/y/z
    y = y ./ Pos2meter;
    z = z ./ Pos2meter;
    rescaled = 1;
end

%% raw data (position, velocity, acceleration) 
% Time 
T = Time;       
dt = median(diff(T));   % sampling period

% raw position 
P.X = x;  %again changed head_x to x 
P.Y = y; 
P.Z = z;

% raw velocity 
%V.Y = gradient(P.Y, dt) ;
%V.Z = gradient(P.Z, dt) ;

% raw acceleration
%A.X = gradient(V.X, dt) ;
%A.Y = gradient(V.Y, dt) ;
%A.Z = gradient(V.Z, dt) ;

%% replace "problems" with nan
[Pn, iTooBigAcc] = Marker_BigAcc2Nan(T, P, DoDebug); 

%% Interpolate nan 
% Pi Position interpolated with splines 
Pi.X = InterpNan(T, Pn.X); 
Pi.Y = InterpNan(T, Pn.Y); 
Pi.Z = InterpNan(T, Pn.Z);

% Vi velocity  
%Vi.X = gradient(Pi.X, dt);
%Vi.Y = gradient(Pi.Y, dt);
%Vi.Z = gradient(Pi.Z, dt);

% Ai acceleration  
%Ai.X = gradient(Vi.X, dt);
%Ai.Y = gradient(Vi.Y, dt);
%Ai.Z = gradient(Vi.Z, dt);

%% get low pass Butterworth filtered data 
SampFreq = 1/dt;  CutFreq = 12; 

Pf.X = LowPassButtDouble (Pi.X, SampFreq, CutFreq);
Pf.Y = LowPassButtDouble (Pi.Y, SampFreq, CutFreq);
Pf.Z = LowPassButtDouble (Pi.Z, SampFreq, CutFreq);

%Vf.X = LowPassButtDouble (Vi.X, SampFreq, CutFreq);
%Vf.Y = LowPassButtDouble (Vi.Y, SampFreq, CutFreq);
%Vf.Z = LowPassButtDouble (Vi.Z, SampFreq, CutFreq);

%Af.X = LowPassButtDouble (Ai.X, SampFreq, CutFreq);
%Af.Y = LowPassButtDouble (Ai.Y, SampFreq, CutFreq);
%Af.Z = LowPassButtDouble (Ai.Z, SampFreq, CutFreq);

%% get SavGol filtered data
%order = 6 ; frameSize = 25; 
%Ps.X = sgolayfilt(Pi.X, order, frameSize);
%Ps.Y = sgolayfilt(Pi.Y, order, frameSize);
%Ps.Z = sgolayfilt(Pi.Z, order, frameSize);

%Vs.X = sgolayfilt(Vi.X, order, frameSize);
%Vs.Y = sgolayfilt(Vi.Y, order, frameSize);
%Vs.Z = sgolayfilt(Vi.Z, order, frameSize);

%As.X = sgolayfilt(Ai.X, order, frameSize);
%As.Y = sgolayfilt(Ai.Y, order, frameSize);
%As.Z = sgolayfilt(Ai.Z, order, frameSize);


%% comparison of raw and nan position time series 
% figure();
% plot(T, [P.X, P.Y, P.Z]);
% hold on
% plot(T, [Pn.X, Pn.Y, Pn.Z], '*')
% xlabel('Time')
% ylabel('Position')
% legend('raw X', 'raw Y', 'raw Z', 'nan X', 'nan Y', 'nan Z' )

%% comparison of raw,  LP and SG acceleration 
%figure(); 
%subplot(3, 1, 1) ; hold on; 
    %plot(T, Ai.X, '-k'); 
    %plot(T, Af.X, '-b', 'linewidth', 2)        
    %plot(T, As.X, '-r', 'linewidth', 2)        
    %ylabel('Acc X (m/s/s)')
%subplot(3, 1, 2) ; hold on;
    %plot(T, Ai.Y, '-k'); 
    %plot(T, Af.Y, '-b', 'linewidth', 2)        
    %plot(T, As.Y, '-r', 'linewidth', 2)
    %ylabel('Acc Y (m/s/s)')

%subplot(3, 1, 3) ; hold on;
    %plot(T, Ai.Z, '-k');     
    %plot(T, Af.Z, '-b', 'linewidth', 2)        
    %plot(T, As.Z, '-r', 'linewidth', 2)
    %ylabel('Acc Z (m/s/s)')
    %xlabel('Time (sec)')
    %legend('raw acc',  sprintf('LP %d Hz', CutFreq), sprintf('SavGol %d %d', order, frameSize ))
%% comparison of the frequency spectrum obtained with LP and SG  
%figure()
%periodogram([A.X, Ai.X, Af.X, As.X])
%legend('Raw', 'cleaned', 'Low pass filter', 'SavGol filter')
%% Save processed position data to txt file
%writetable(struct2table(Pf), 'somefile.txt')

directory = '/Users/olgamatthiopoulou/Desktop/CleanHumanMovementTimeSeries-master/DAT/07-13/t_061/output/B';
fileout = [directory filesep 't_061_B_split_20_out.txt'];

writetable(struct2table(Pf), fileout)

clear 
close all
