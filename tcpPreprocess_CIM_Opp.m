%% Initialization
% Stephen Zhang 2019/07/30
% Edited 230106 Crystian Massengill

% Use previous path if exists
if exist('filepath', 'var')
    if exist('ppCfg', 'var')
        defaultpath = filepath;
        keep defaultpath ppCfg
    else
        defaultpath = filepath;
        keep defaultpath
    end
elseif exist('filepath2', 'var')
    defaultpath = filepath2;
    keep defaultpath
else
    clear
    % common path
    defaultpath = '\\anastasia\data\photometry';
end

% Use UI
hfig = ppCfg_UI(defaultpath);
waitfor(hfig);

% Parse outcome
OPTO_MODE = ppCfg.OPTO_MODE;
PULSE_SIM_MODE = ppCfg.PULSE_SIM_MODE;
MONKEYLOGIC_MODE = false; % manually off for now
data_channel = ppCfg.data_channel;
data_channel2 = ppCfg.data_channel2;
opto_channel = ppCfg.opto_channel;
VF_channel = 13; %This is manual button for vonfrey/other manual stims
ch1_pulse_ind = ppCfg.ch1_pulse_ind;
ch2_pulse_ind = ppCfg.ch2_pulse_ind;
ch1_pulse_thresh = ppCfg.ch1_pulse_thresh;
ch2_pulse_thresh = ppCfg.ch2_pulse_thresh;
filt_stim = ppCfg.filt_stim;
stim_filt_range = ppCfg.stim_filt_range;
use_fnotch_60 = ppCfg.use_fnotch_60;
fnotch_60 = ppCfg.fnotch_60;
blackout_window = ppCfg.blackout_window;
freq = ppCfg.freq;
Ambientpts = ppCfg.Ambientpts;
tone_channel = ppCfg.tone_channel;

%added for ML alignment
Vis_Stim = 10;  %This is the DIO for monkeylogic trial structure

% To manually toggle on/off
% OPTO_MODE = 0;
% VONFREY_MODE = 0;
% MONKEYLOGIC_MODE = 0;
%% IO
% Work out outputpath
[filename, filepath] = uigetfile(fullfile(defaultpath, '*.mat'));
filename_output = [filename(1:end-4), '_preprocessed.mat'];
load(fullfile(filepath, filename), 'data', 'timestamps', 'Fs');

mouse=filename(1:end-21);
date=filename(8:end-14);
run=filename(15:end-10);
%% load BHV2 file, extract data into BHV format
if MONKEYLOGIC_MODE    
    [MLdata,MLConfig,TrialRecord] = mlconcatenate_opto();
    ml = bhv2_to_bhv_RE(MLdata, TrialRecord);
else
    ml = [];
end
%% Basic channel info and point indices
% Gathering pulses

if PULSE_SIM_MODE
    [ch1_pulse, ch2_pulse] = pulsesim(size(data,2), 2500, 9, 10);
else
    % Grab pulse info
    ch1_pulse = data(ch1_pulse_ind,:) > ch1_pulse_thresh;
    ch2_pulse = data(ch2_pulse_ind,:) > ch2_pulse_thresh;
end

% Find pulse points. This step also defines the sampling rate after
% downsampling (which is the rate of pulses)
ch1_data_table = chainfinder(ch1_pulse);
ch2_data_table = chainfinder(ch2_pulse);

% Rearrange data
ch1_data_table(:,3) = ch1_data_table(:,2);
ch1_data_table(:,2) = nan;
ch2_data_table(:,3) = ch2_data_table(:,2);
ch2_data_table(:,2) = nan;

% Equalize the pulse numbers of the two wavelenghts
n_points = min(size(ch1_data_table(:,1),1),size(ch2_data_table(:,2),1)) - 1;
% n_points = min(size(ch2_data_table(:,2),1)) - 1;
% n_points = min(size(ch2_data_table(:,2),1)) - 1;

% Fix pulse 1 if needed
if size(ch1_data_table,1) > n_points
    ch1_data_table = ch1_data_table(1:n_points, :);
end

% Fix pulse 2 if needed
if size(ch2_data_table,1) > n_points
    ch2_data_table = ch2_data_table(1:n_points, :);
end

%% Check if the running file is there
runningfn = sprintf('%srunning.mat', filename(1:end-9));
runningfn_full = fullfile(filepath, runningfn);

if exist(runningfn_full, 'file')
    % Load running data
    running = load(runningfn_full, 'speed');
    
    %Running running sample count
    nrunpulse = size(chainfinder(data(ch2_pulse_ind,:)>0.5),1);
    nrunlength = length(running.speed);
    
    if nrunpulse ~= nrunlength
        % Say something
        fprintf('Running digitization is %0.3f%% off\n', (1 - nrunlength/nrunpulse)*100);
        
        % Upsample running data
        speed_upsampled0 = TDresamp(running.speed', 'resample', nrunpulse/nrunlength * 0.9974);
        speed_upsampled = TDresamp(speed_upsampled0, 'resample',n_points/nrunpulse);
    else
        % Upsample running data
        speed_upsampled = TDresamp(running.speed', 'resample',n_points/length(running.speed));
    end
    
    % Fix the number of points if needed
    if length(speed_upsampled) > n_points
        speed_upsampled = speed_upsampled(1:n_points);
    elseif length(speed_upsampled) < n_points
        speed_upsampled(end:end + n_points - length(speed_upsampled)) = 0;
    end
else
    % Store empty speed matrices
    speedupsampled = [];
end

%% Notch filters
% CIM added accomodation for 2nd channel filtering 230106
if use_fnotch_60
    % Apply notch filter to remove 60 Hz noise
    d_notch = designfilt('bandstopiir','FilterOrder',2, 'HalfPowerFrequency1',...
        fnotch_60(1), 'HalfPowerFrequency2',fnotch_60(2), 'DesignMethod','butter','SampleRate', Fs);
    data_notch_ch1 = filter(d_notch, data(data_channel,:));
    data_notch_ch2 = filter(d_notch, data(data_channel2,:));

else
    data_notch_ch1 = data(data_channel,:);
    data_notch_ch2 = data(data_channel2,:);
end

if OPTO_MODE
    if filt_stim
        d_notch_stim = designfilt('bandstopiir','FilterOrder',2, 'HalfPowerFrequency1',...
            stim_filt_range(1), 'HalfPowerFrequency2',stim_filt_range(2),...
            'DesignMethod','butter','SampleRate', Fs);
        data_notch_ch1 = filter(d_notch_stim, data_notch_ch1);
        data_notch_ch2 = filter(d_notch_stim, data_notch_ch2);
    end
end
% 
%     d2 = designfilt('bandstopiir','FilterOrder',2, ...
%                'HalfPowerFrequency1',0.5,'HalfPowerFrequency2',2, ...
%                'DesignMethod','butter','SampleRate',Fs);
%     data_notch_2 = filter(d2, data_notch_ch1);

    
%     figure;plot(data_notch_ch1);hold on;plot(data_notch_2)
%%
% Apply another filter to filter out stim artifacts if needed
if OPTO_MODE
    if filt_stim
        d_notch_stim = designfilt('bandstopiir','FilterOrder',2, 'HalfPowerFrequency1',...
            stim_filt_range(1), 'HalfPowerFrequency2',stim_filt_range(2),...
            'DesignMethod','butter','SampleRate', Fs);

            d_notch_stim2 = designfilt('bandpassfir','FilterOrder',2, 'CutoffFrequency1',...
            22, 'CutoffFrequency2',24,...
            'SampleRate', Fs);
% 
        d_notch_stim = designfilt('highpassiir','FilterOrder',2, 'PassbandFrequency', 5,...
            'PassbandRipple', 0.5,'SampleRate', Fs);
        fvtool(d_notch_stim)
Example 1:
  %   Design a highpass IIR filter with order 8, passband frequency of 
  %   75 KHz, and a passband ripple of 0.2 dB. Sample rate is 200 KHz.
  %   Visualize the filter response and apply the filter to a vector of
  %   random data. 

  hpFilt = designfilt('highpassiir', 'FilterOrder', 8, ...
           'PassbandFrequency', 75e3, 'PassbandRipple', 0.2,...
           'SampleRate', 200e3);
        
        
        data_notch_2 = filter(d_notch_stim2, data_notch_ch1);
        data_notch_ch2 = filter(d_notch_stim, data_notch_ch2);
    end
end


% Grab data points
% Use median fluorescence during each pulse to calculate fluorescence
% values
%%
for i = 1 : n_points
    % Wavelength 1
    ini_ind = ch1_data_table(i,1) + blackout_window;
    end_ind = ch1_data_table(i,1) + ch1_data_table(i,3) - 1;
    ch1_data_table(i,2) = median(data_notch_ch1(ini_ind:end_ind));
%     
    % Wavelength 2
    ini_ind = ch2_data_table(i,1) + blackout_window;
    end_ind = ch2_data_table(i,1) + ch2_data_table(i,3) - 1;
    ch2_data_table(i,2) = median(data_notch_ch2(ini_ind:end_ind));
end

% figure;
% plot([ch2_data_table(2000:4536,2)])
% ylim([0 6])
%% Ambient-light subtraction
if Ambientpts > 0
    % Initialize matrices
    ch1_amb_table = nan(size(ch1_data_table));
    ch2_amb_table = nan(size(ch2_data_table));
    
    % Get the indices
    ch1_amb_table(:,1) = ch1_data_table(:,1) - Ambientpts;
    ch1_amb_table(:,3) = Ambientpts;
    
    ch2_amb_table(:,1) = ch2_data_table(:,1) - Ambientpts;
    ch2_amb_table(:,3) = Ambientpts;
   
    % Fix out-of-bount indices by taking the next point
    if ch1_amb_table(1,1) < 1
        ch1_amb_table(1,1) = ch1_amb_table(2,1);
    end
    if ch2_amb_table(1,1) < 1
        ch2_amb_table(1,1) = ch2_amb_table(2,1);
    end
    
    % Loop through and take median
    for i = 1 : n_points
        % Wavelength 1
        ini_ind = ch1_amb_table(i,1);
        end_ind = ch1_amb_table(i,1) + ch1_amb_table(i,3) - 1;
        ch1_amb_table(i,2) = median(data_notch(ini_ind:end_ind));

        % Wavelength 2
        ini_ind = ch2_amb_table(i,1);
        end_ind = ch2_amb_table(i,1) + ch2_amb_table(i,3) - 1;
        ch2_amb_table(i,2) = median(data_notch(ini_ind:end_ind));
    end
    
    % Subtract
    ch1_data_table(:,2) = ch1_data_table(:,2) - ch1_amb_table(:,2);
    ch2_data_table(:,2) = ch2_data_table(:,2) - ch2_amb_table(:,2);
end

%% Grab opto pulses
if OPTO_MODE
    % Grab the pulses
    opto_pulse_table = tcpDatasnapper(data(opto_channel,:),data(datachannel,:), 'max', 'pulsetopulse');  
    % Sync the number of pulses
    opto_pulse_table = opto_pulse_table(1 : n_points, :);
else
    opto_pulse_table = [];
end

%% Grab button presses
if VONFREY_MODE
    % Grab the pulses
    button_press_array = detect_button_press(data(VF_channel,:), 125);
    VF_pulse_table = tcpDatasnapper(button_press_array,...
        data(ch2_pulse_ind,:), 'max', 'pulsetopulse');  
    % Sync the number of pulses
    VF_pulse_table = VF_pulse_table(1 : n_points, :);
else
    VF_pulse_table = [];
end

%% Grab ML trial structure pulses
%
if MONKEYLOGIC_MODE
    % Grab the pulses
    ML_pulse_table = tcpDatasnapper(data(Vis_Stim,:),...
        data(ch1_pulse_ind,:), 'max', 'pulsetopulse');  
    
    % Sync the number of pulses
    ML_pulse_table = ML_pulse_table(1 : n_points, :);
    
else
    ML_pulse_table = [];
end


%% Grab tone pulses
if tone_channel < 99
     % Grab the pulses
    tone_pulse_table = tcpDatasnapper(data(:,tone_channel),...
        data(:,ch1_pulse_ind), 'max', 'pulsetopulse');
    
    % Sync the number of pulses
    tone_pulse_table = tone_pulse_table(1 : n_points, :);
else
    tone_pulse_table = [];
end


%% Copy opto table from a different experiment (debug)
%{
[optofn, optofp] = uigetfile(fullfile(filepath, '*.mat'));
loaded = load(fullfile(optofp, optofn), 'opto_pulse_table');
opto_pulse_table = ch1_data_table;
opto_pulse_table(:,2) = 0;
opto_pulse_table(:,3) = loaded.opto_pulse_table(2,3);
ptable = chainfinder(loaded.opto_pulse_table(:,2) > 0.5);
for i = 1 : size(ptable,1)
    istart = ptable(i,1);
    iend = ptable(i,1) + ptable(i,2) - 1;
    opto_pulse_table(istart:iend,2) = 1;
end
%}

%% filter for any movement
% clear smoothed_speed speed_logical ch1_data_speedfilt ch2_data_speedfilt
% 
% %smooth a bit
% smoothed_speed = movmean(speed_upsampled,2);
% 
% lower_threshold = -0.2;
% upper_threshold = 0.2;
% 
% speed_logical = (smoothed_speed < lower_threshold) | (smoothed_speed > upper_threshold);
% 
% onesPositions = find(speed_logical);
% % adding padding bc of slow activity changes
% leftPaddingSize = 10;   % Adjust the left padding size as needed
% rightPaddingSize = 250;  % Adjust the right padding size as needed
% for pos = onesPositions'
%     startLeft = max(1, pos - leftPaddingSize);
%     endLeft = min(length(speed_logical), pos - 1);
%     speed_logical(startLeft:endLeft) = 1;
% 
%     startRight = max(1, pos + 1);
%     endRight = min(length(speed_logical), pos + rightPaddingSize);
%     speed_logical(startRight:endRight) = 1;
% end
% 
% ch1_data_speedfilt=ch1_data_table(:,2);
% ch2_data_speedfilt=ch2_data_table(:,2);
% 
% ch1_data_speedfilt(speed_logical == 1) = NaN;
% ch2_data_speedfilt(speed_logical == 1) = NaN;
% 
% 
% figure('Renderer', 'painters', 'Position', [100 100 1500 600]); hold on
% ax1 = subplot(3,1,1:2);
% plot((1 : n_points)'/freq, ch1_data_speedfilt)
% ylim([1 3])
% ylabel('Signal')
% set(gca,'FontSize',10,'box','off');
% names={'0','10','20','30','40','50'};
% set(gca,'xtick',[0,600,1200,1800,2400,3000],'xticklabel',names);
% xlabel('time (min)')
% line([(11.75*60) (11.75*60)], [1,10],'Color','red');
% line([(22.6*60) (22.6*60)], [1,10],'Color','red');
% 
% 
% % yyaxis right
% % plot((1 : n_points)'/freq, ch2_data_speedfilt)
% ax2 = subplot(3,1,3);
% plot((1 : n_points)'/freq, speed_upsampled)  
% linkaxes([ax1,ax2],'x')
% set(gca,'FontSize',10,'box','off');
% set(gca,'xtick',[0,600,1200,1800,2400,3000],'xticklabel',names);
% line([(11.75*60) (11.75*60)], [-50 50],'Color','red');
% line([(22.6*60) (22.6*60)], [-50 50],'Color','red');
% xlabel('time (min)')
% ylabel('Running')

% data1=ch1_data_table(12000:end,2);
% data2=ch2_data_table(12000:end,2);
% r=corrcoef(data1,data2);
%% Plot raw data
figure('Renderer', 'painters', 'Position', [250 250 750 550]); hold on
set(gcf,'Color','w');
% Plot raw fluorescence data on the left
ax1 = subplot(4,1,1:2);
plot((1 : n_points)'/freq, [ch1_data_table(:,2)]); hold on;
plot((1 : n_points)'/freq, [ch2_data_table(:,2)]);
leg = legend('Left fiber','Right fiber');
ylim([1 5])
ylabel('Raw voltage')
title((strcat(date, '_', mouse, '_', run)),'Interpreter', 'none')   
ax2 = subplot(4,1,3);
plot((1 : n_points)'/freq, [ML_pulse_table(:,2)])  
ylabel('Shock')
xlabel('Time (s)')
ax3 = subplot(4,1,4);
plot((1 : n_points)'/freq, speed_upsampled)  
ylim([-10 30])
ylabel('Speed')
xlabel('Time (s)')
linkaxes([ax1,ax2,ax3],'x');
filen = strcat(date, '_', mouse, '_',run); 
saveas(gcf, [defaultpath, '\', filen, '.png'])
%% Power analysis
% FFT (data, sampling rate, don't plot)
% if OPTO_MODE
%     [Powers, fft_freq] = ft2(ch1_data_table(2:end,2), 50, 0);
% else
%     [Powers, fft_freq] = ft2([ch1_data_table(2:end,2) , ch2_data_table(2:end,2)] , 50, 0);
% end
% 
% % Plot FFT info
% subplot(1,3,2)
% plot(fft_freq, Powers)
% xlim([1, freq/2])
% ylabel('Power')
% xlabel('Frequency')


%% Low pass filter
Ch1_filtered = ch1_data_table(:,2);
Ch2_filtered = ch2_data_table(:,2);
% Interpolate method (last resort)
ch1_data_table = artifact_glm(ch1_data_table, data, 4, 9);
ch2_data_table = artifact_glm(ch2_data_table, data, 4, 9);

% Filter
% Design a filter kernel
d = fdesign.lowpass('Fp,Fst,Ap,Ast',8,10,0.5,40, freq);
Hd = design(d,'equiripple');
% fvtool(Hd)

% Filter data
% % Ch1_filtered = filter(Hd,ch1_data_table(:,2));
Ch2_filtered = filter(Hd,ch2_data_table(:,2));

% Plot filtered fluorescence data on the right
figure(100)
plot((1 : n_points)'/freq, Ch2_filtered)
xlabel('Time (s)')
ylabel('Photodiod voltage (V)')

%% Save
save(fullfile(filepath, filename_output), 'ch1_data_table','ch2_data_table','date','mouse','run',...
    'data', 'freq', 'Fs','n_points','Ch1_filtered','Ch2_filtered', 'ML_pulse_table','PULSE_SIM_MODE', 'OPTO_MODE','MONKEYLOGIC_MODE',...
    'ml','timestamps', 'speed_upsampled','speed_upsampled','opto_pulse_table', 'ppCfg');
