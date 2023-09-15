%% Simulating Test Signals for a Radar Receiver
% 

%% Design Specifications

pd = 0.9;            % Probability of detection
pfa = 1e-6;          % Probability of false alarm
max_range = 5000;    % Maximum unambiguous range
range_res = 50;      % Required range resolution
tgt_rcs = 1;         % Required target radar cross section

%% Monostatic Radar System Design
prop_speed = physconst('LightSpeed');   % Propagation speed
pulse_bw = prop_speed/(2*range_res);    % Pulse bandwidth
pulse_width = 1/pulse_bw;               % Pulse width
prf = prop_speed/(2*max_range);         % Pulse repetition frequency
fs = 2*pulse_bw;                        % Sampling rate
waveform = phased.RectangularWaveform(...
    'PulseWidth',1/pulse_bw,...
    'PRF',prf,...
    'SampleRate',fs);

%%


noise_bw = pulse_bw;

receiver = phased.ReceiverPreamp(...
    'Gain',20,...
    'NoiseFigure',0,...
    'SampleRate',fs,...
    'EnableInputPort',true);

%%


snr_db = [-inf, 0, 3, 10, 13];
rocsnr(snr_db,'SignalType','NonfluctuatingNoncoherent');

%%


num_pulse_int = 10;
rocsnr([0 3 5],'SignalType','NonfluctuatingNoncoherent',...
    'NumPulses',num_pulse_int);

snr_min = albersheim(pd, pfa, num_pulse_int)

%%


tx_gain = 20;

fc = 10e9;
lambda = prop_speed/fc;

peak_power = ((4*pi)^3*noisepow(1/pulse_width)*max_range^4*...
    db2pow(snr_min))/(db2pow(2*tx_gain)*tgt_rcs*lambda^2)
%%


%%
% With all this information, we can configure the transmitter.

transmitter = phased.Transmitter(...
    'Gain',tx_gain,...
    'PeakPower',peak_power,...
    'InUseOutputPort',true);


%
% We assume that the antenna is stationary.

antenna = phased.IsotropicAntennaElement(...
    'FrequencyRange',[5e9 15e9]);

sensormotion = phased.Platform(...
    'InitialPosition',[0; 0; 0],...
    'Velocity',[0; 0; 0]);

%%
% With the antenna and the operating frequency, we define both the radiator
% and the collector.

radiator = phased.Radiator(...
    'Sensor',antenna,...
    'OperatingFrequency',fc);

collector = phased.Collector(...
    'Sensor',antenna,...
    'OperatingFrequency',fc);

%%
% This completes the configuration of the radar system. 

%% System Simulation
% *Targets*
%
% To test our radar's ability to detect targets, we must define the targets
% first. Let us assume that there are 3 stationary, non-fluctuating targets
% in space. Their positions and radar cross sections are given below.
tgtpos = [[2024.66;0;0],[3518.63;0;0],[3845.04;0;0]];
tgtvel = [[0;0;0],[0;0;0],[0;0;0]];
tgtmotion = phased.Platform('InitialPosition',tgtpos,'Velocity',tgtvel);

tgtrcs = [1.6 2.2 1.05];
target = phased.RadarTarget('MeanRCS',tgtrcs,'OperatingFrequency',fc);

%% 
% *Propagation Environment*
%
% To simulate the signal, we also need to define the propagation channel
% between the radar system and each target. 
channel = phased.FreeSpace(...
    'SampleRate',fs,...
    'TwoWayPropagation',true,...
    'OperatingFrequency',fc);

%%
% 

fast_time_grid = unigrid(0,1/fs,1/prf,'[)');
slow_time_grid = (0:num_pulse_int-1)/prf;

%%
% The following loop simulates 10 pulses of the receive signal.
%
% We set the seed for the noise generation in the receiver so that we can
% reproduce the same results.

receiver.SeedSource = 'Property';
receiver.Seed = 2007;

% Pre-allocate array for improved processing speed
rxpulses = zeros(numel(fast_time_grid),num_pulse_int);

for m = 1:num_pulse_int
    
    % Update sensor and target positions
    [sensorpos,sensorvel] = sensormotion(1/prf);
    [tgtpos,tgtvel] = tgtmotion(1/prf);

    % Calculate the target angles as seen by the sensor
    [tgtrng,tgtang] = rangeangle(tgtpos,sensorpos);
    
    % Simulate propagation of pulse in direction of targets
    pulse = waveform();
    [txsig,txstatus] = transmitter(pulse);
    txsig = radiator(txsig,tgtang);
    txsig = channel(txsig,sensorpos,tgtpos,sensorvel,tgtvel);
    
    % Reflect pulse off of targets
    tgtsig = target(txsig);
    
    % Receive target returns at sensor
    rxsig = collector(tgtsig,tgtang);
    rxpulses(:,m) = receiver(rxsig,~(txstatus>0));
end

%% Range Detection
% *Detection Threshold*
%
%

npower = noisepow(noise_bw,receiver.NoiseFigure,...
    receiver.ReferenceTemperature);
threshold = npower * db2pow(npwgnthresh(pfa,num_pulse_int,'noncoherent'));

%%
% We plot the first two received pulses with the threshold
num_pulse_plot = 2;
helperRadarPulsePlot(rxpulses,threshold,...
    fast_time_grid,slow_time_grid,num_pulse_plot);

%%


matchingcoeff = getMatchedFilter(waveform);
matchedfilter = phased.MatchedFilter(...
    'Coefficients',matchingcoeff,...
    'GainOutputPort',true);
[rxpulses, mfgain] = matchedfilter(rxpulses);

%%


matchingdelay = size(matchingcoeff,1)-1;
rxpulses = buffer(rxpulses(matchingdelay+1:end),size(rxpulses,1));

%%
% The threshold is then increased by the matched filter processing gain.
threshold = threshold * db2pow(mfgain);

%% 
% The following plot shows the same two pulses after they pass through the
% matched filter.
helperRadarPulsePlot(rxpulses,threshold,...
    fast_time_grid,slow_time_grid,num_pulse_plot);

%%


range_gates = prop_speed*fast_time_grid/2; 

tvg = phased.TimeVaryingGain(...
    'RangeLoss',2*fspl(range_gates,lambda),...
    'ReferenceLoss',2*fspl(max_range,lambda));

rxpulses = tvg(rxpulses);

%%
% Now let's plot the same two pulses after the range normalization 
helperRadarPulsePlot(rxpulses,threshold,...
    fast_time_grid,slow_time_grid,num_pulse_plot);

%%


%%


rxpulses = pulsint(rxpulses,'noncoherent');

helperRadarPulsePlot(rxpulses,threshold,...
    fast_time_grid,slow_time_grid,1);



[~,range_detect] = findpeaks(rxpulses,'MinPeakHeight',sqrt(threshold));

%% 
% The true ranges and the detected ranges of the targets are shown below:

true_range = round(tgtrng)
range_estimates = round(range_gates(range_detect))
figure
scatter(true_range,range_estimates)
%%
scenario = drivingScenario('StopTime',3);
roadcenters = [-35 20 0; -20 -20 0; 0 0 0; 20 20 0; 35 -20 0];
lm = [laneMarking('Solid','Color','w'); ...
    laneMarking('Dashed','Color','y'); ...
    laneMarking('Dashed','Color','y'); ...
    laneMarking('Solid','Color','w')];
ls = lanespec(3,'Marking',lm);
road(scenario,roadcenters,'Lanes',ls);
car = vehicle(scenario, ...
    'ClassID',1, ...
    'Position',[-35 20 0]);
waypoints = [-35 20 0; -20 -20 0; 0 0 0; 20 20 0; 35 -20 0];
smoothTrajectory(car,waypoints);
plot(scenario)
chasePlot(car)
bep = birdsEyePlot('XLim',[-40 40],'YLim',[-30 30]);
olPlotter = outlinePlotter(bep);
lblPlotter = laneBoundaryPlotter(bep,'Color','r','LineStyle','-');
lbrPlotter = laneBoundaryPlotter(bep,'Color','g','LineStyle','-');
rbsEdgePlotter = laneBoundaryPlotter(bep);
legend('off');
while advance(scenario)
    rbs = roadBoundaries(car);
    [position,yaw,length,width,originOffset,color] = targetOutlines(car);
    lb = laneBoundaries(car,'XDistance',0:5:30,'LocationType','Center', ...
        'AllBoundaries',false);
    plotLaneBoundary(rbsEdgePlotter,rbs)
    plotLaneBoundary(lblPlotter,{lb(1).Coordinates})
    plotLaneBoundary(lbrPlotter,{lb(2).Coordinates})
    plotOutline(olPlotter,position,yaw,length,width, ...
        'OriginOffset',originOffset,'Color',color)
end

