%% Initialization
% Crystian Massengill 04092024
clear all
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
MONKEYLOGIC_MODE = true;
data_channel = ppCfg.data_channel;
data_channel2 = ppCfg.data_channel2;
opto_channel = ppCfg.opto_channel;
ch1_pulse_ind = ppCfg.ch1_pulse_ind;
ch2_pulse_ind = ppCfg.ch2_pulse_ind;
ch1_pulse_thresh = ppCfg.ch1_pulse_thresh;
ch2_pulse_thresh = ppCfg.ch2_pulse_thresh;
filt_stim = ppCfg.filt_stim;
stim_filt_range = ppCfg.stim_filt_range;
use_fnotch_60 = ppCfg.use_fnotch_60;
fnotch_60 = ppCfg.fnotch_60;
blackout_window = ppCfg.blackout_window;
freq = 30;
Ambientpts = ppCfg.Ambientpts;
tone_channel = ppCfg.tone_channel;

%CIM added for ML alignment
Vis_Stim = 10;  %This is the DIO for monkeylogic trial structure

% Window info
prew = 5;
postw = 55;
prew_f = prew * freq;
postw_f = postw * freq;
l = prew_f + postw_f + 1;

%% IO
% Work out outputpath
[filename, filepath] = uigetfile(fullfile(defaultpath, '*.mat'));
filename_output = [filename(1:end-4), '_preprocessed.mat'];
load(fullfile(filepath, filename), 'data', 'timestamps', 'Fs');

mouse=filename(1:end-21);
date=filename(8:end-14);
run=filename(15:end-10);

%% Check if the running file is there
runningfn = sprintf('%srunning.mat', filename(1:end-9));
runningfn_full = fullfile(filepath, runningfn);

if exist(runningfn_full, 'file')
    % Load running data
    running = load(runningfn_full, 'speed');
    
    %Running sample count
    nrunpulse = size(chainfinder(data(7,:)>0.5),1);
    nrunlength = length(running.speed);
else
    % Store empty speed matrices
    speedupsampled = [];
end

%%
opto_pulse_smoothed = smooth(data(6,:),1000)';
opto_pulse = opto_pulse_smoothed > 0.01;
noptopulse = size(opto_pulse,2); 
downfactor = round(noptopulse/nrunpulse);
opto_pulse = downsample(opto_pulse,downfactor);
noptolength = length(opto_pulse);
if nrunlength ~= noptolength 
    opto_pulse = opto_pulse(:,1:nrunpulse);
end
% figure;plot(opto_pulse_smoothed)

opto_ons = chainfinder(opto_pulse > 0.5);

%% Grab the point indices
clear inds_master badtrials
% Indices
inds_master = cell(1,1); %Initialize cell array

inds_master{:,1}=opto_ons(:,1) * [1 1];
inds_master{:,1}(:,1) =  inds_master{:,1}(:,1) - prew_f;
inds_master{:,1}(:,2) =  inds_master{:,1}(:,2) + postw_f;
n_trials = length(inds_master{1,1});

%%
badtrials = ((inds_master{1,1} - prew_f) <= 0) + ((inds_master{1,1} + postw_f) > nrunlength);
inds_master{1,1}(10,:)=[];
% badtrials =[];
n_trials = length(inds_master{1,1});
%% Deal with motion
% Initialize a triggered speed matrix
speedmat = cell(1,1);
for i = 1 : n_trials
    speedmat{1,1}(1:l,i)=zeros(1,1);
end
%populate matrix
for i = 1 : n_trials
   speedmat{1,1}(:,i) = running.speed(inds_master{:,1}(i,1) : inds_master{:,1}(i,2));
end
% 
% Calculate the average triggered results
for i = 1 : n_trials
    speedmat_avg(:,1) = mean(speedmat{1,1},2,'omitnan');
end

%%
y_ax = [-5 15];
f1 = figure(); 
f1.Position = [50 0 1200 800];
set(gcf,'Color','w');
hold on;
seconds = -prew : 1/freq : postw;

    toPlot = {speedmat{:,:}}';
    speedPlot = speedmat_avg';

    ax1=subplot(2,1,1)
    imagesc(seconds, [1 size(toPlot{1,1}', 1)], toPlot{1,1}'); hold on
    line([0 0], y_ax, 'linestyle', '--', 'linewidth', 1.5, 'color', 'k'); hold on
    
    title((strcat(date, '_', mouse,' - Spinal axon opto')),'Interpreter', 'none')   
    ylim([1 n_trials])
    xlim([-prew postw;])
    xlabel('Time (s)')
    ylabel('Trials')
    colormap(flipud(brewermap([],'RdBu')))
    a=colorbar; caxis(y_ax)
    ylabel(a,'cm/s','FontSize',16);
    set(gca, 'fontsize', 14)
    
    ax2=subplot(2,1,2)
    plot(seconds, toPlot{1,1}', 'Color',0.8*[0.75 0.75 0.75]); hold on
    plot(seconds, speedPlot', 'Color','k','LineWidth',4); hold on
    ylim([y_ax]);
    xlim([-prew postw;]);
    line([0 0], y_ax, 'linestyle', '--', 'linewidth', 1.5, 'color', 'k'); hold on
    set(gca, 'fontsize', 14)
    set(gca, 'box','off')
    linkaxes([ax1,ax2],'x')

filename = strcat(date, '_', mouse, '_MeanTrace_Run',run); 
saveas(gcf, [filepath, '\', filename, '.png'])
%% Plot raw data
figure('Renderer', 'painters', 'Position', [100 100 1800 450]); hold on
set(gcf,'Color','w');
plot((1 : nrunpulse)'/freq, [running.speed])  
%     imagesc((1 : nrunpulse)'/freq, [1 size(running.speed, 1)], running.speed); hold on
set(gca,'FontSize',10,'box','off');
% set(gca,'Xticklabel',[])
ylabel('Running')
% % % yyaxis right
% plot((1 : nrunpulse)'/freq, [opto_pulse])   
yticks([0 1])
% ylabel('Ensure Trials')  
xlabel('Time (s)')

filename = strcat(date, '_', mouse, '_RawTrace_Run',run); 
saveas(gcf, [filepath, '\', filename, '.png']);

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
% figure('Renderer', 'painters', 'Position', [100 100 750 450]); hold on
% set(gcf,'Color','w');
% % Plot raw fluorescence data on the left
% % ylim([0.2 1])
% % yyaxis right
% % % plot((1 : n_points)'/freq, [ch2_data_table(:,2)])
% % ylabel('Signal')
% % set(gca,'FontSize',10,'box','off');
% % % set(gca,'Xticklabel',[])
% % % ylim([0.2 1])
% % lgd = legend('Location', 'northeast');
% %     legend('boxoff');
% %     lgd.ItemTokenSize = [30,18];
% % %     if(rem((i/2),1)==0)
% % legend('Ch1:465','Ch2:405');
% % title((strcat(date, '_', mouse)),'Interpreter', 'none')   
% % % line([1190 1190], [1,6],'Color','red');
% 
% 
% % subplot(5,1,3)
% % plot((1 : n_points)'/freq, [opto_pulse_table(:,2)])   
% % set(gca,'FontSize',10,'box','off');
% % set(gca,'Xticklabel',[])
% % ylabel('Opto')
% 
% ax1 = subplot(2,1,1)
% plot((1 : n_points)'/freq, [speed_upsampled])  
% set(gca,'FontSize',10,'box','off');
% % set(gca,'Xticklabel',[])
% ylabel('Running')
% % line([1190 1190], [-20,40],'Color','red');
% % linkaxes([ax1,ax2],'x')
% 
% ax2 = subplot(2,1,2)
% if MONKEYLOGIC_MODE
%     plot((1 : n_points)'/freq, [opto_pulse_table(:,1)])   
% else
% end
% % ylabel('Ensure Trials')  
% xlabel('Time (s)')
% linkaxes([ax2,ax3],'x')
% 
% filename = strcat(date, '_', mouse, '_RawTrace_Run',run); 
% saveas(gcf, [filepath, '\', filename, '.png']);
% 
% 
% %% Power analysis
% % FFT (data, sampling rate, don't plot)
% % if OPTO_MODE
% %     [Powers, fft_freq] = ft2(ch1_data_table(2:end,2), 50, 0);
% % else
% %     [Powers, fft_freq] = ft2([ch1_data_table(2:end,2) , ch2_data_table(2:end,2)] , 50, 0);
% % end
% % 
% % % Plot FFT info
% % subplot(1,3,2)
% % plot(fft_freq, Powers)
% % xlim([1, freq/2])
% % ylabel('Power')
% % xlabel('Frequency')
% 
% 
% %% Low pass filter
% Ch1_filtered = ch1_data_table(:,2);
% Ch2_filtered = ch2_data_table(:,2);
% % Interpolate method (last resort)
% ch1_data_table = artifact_glm(ch1_data_table, data, 4, 9);
% ch2_data_table = artifact_glm(ch2_data_table, data, 4, 9);
% 
% % Filter
% % Design a filter kernel
% d = fdesign.lowpass('Fp,Fst,Ap,Ast',8,10,0.5,40, freq);
% Hd = design(d,'equiripple');
% % fvtool(Hd)
% 
% % Filter data
% Ch1_filtered = filter(Hd,ch1_data_table(:,2));
% Ch2_filtered = filter(Hd,ch2_data_table(:,2));
% 
% % Plot filtered fluorescence data on the right
% % figure(100)
% % plot((1 : n_points)'/freq, Ch1_filtered)
% % xlabel('Time (s)')
% % ylabel('Photodiod voltage (V)')
% 
% %% Save
% save(fullfile(filepath, filename_output), 'ch1_data_table', 'ch2_data_table','date','mouse','run',...
%     'data', 'freq', 'Fs','n_points','Ch1_filtered','Ch2_filtered', 'PULSE_SIM_MODE', 'OPTO_MODE','MONKEYLOGIC_MODE',...
%     'ml','timestamps', 'ML_pulse_table','speed_upsampled','opto_pulse_table', 'ppCfg', 'tone_pulse_table');
