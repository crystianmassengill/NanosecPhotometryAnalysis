%% Initialization
% Stephen Zhang 2019/10/20

% Use previous path if exists
if ~exist('filepath', 'var')
    clear
    % common path
    defaultpath = '\\anastasia\data\photometry';
elseif exist('TrigCfg', 'var')
    defaultpath = filepath;
    keep defaultpath TrigCfg
else
    defaultpath = filepath;
    keep defaultpath
end

% Use UI
hfig = TrigCfg_UI(defaultpath);
waitfor(hfig);

%% IO
% Work out outputpath
[filename, filepath] = uigetfile(fullfile(defaultpath , '*.mat'));
if isempty(TrigCfg.suffix)
    filename_output_triggered = [filename(1:end-4), '_trig.mat'];  
else
    filename_output_triggered = sprintf('%s_trig_%s.mat', filename(1:end-4), TrigCfg.suffix);
end
load(fullfile(filepath, filename), 'date','mouse','run','data', 'freq', 'ch1_data_table','ch2_data_table',...
    'n_points','ml','ML_pulse_table','opto_pulse_table','speed_upsampled');

%% GLM remove channel artifacts
% Issue with small NIDAQ 

% replace with aligned/subtracted/filtered data from tcpAlign
% ch1_data_table(:,2) = signal;
Ch1_filtered = ch1_data_table(:,2);
Ch2_filtered = ch2_data_table(:,2);

if isfield(TrigCfg, 'GLM_artifacts')
    if TrigCfg.GLM_artifacts
        % Interpolate method (last resort)
        ch1_data_table = artifact_glm(ch1_data_table, data, TrigCfg.GLM_ch, 9);

        % Filter
        % Design a filter kernel
        d = fdesign.lowpass('Fp,Fst,Ap,Ast',8,10,0.5,40, freq);
        Hd = design(d,'equiripple');
        % fvtool(Hd)

        % Filter data
        Ch1_filtered = filter(Hd,ch1_data_table(:,2));
    end
else
    TrigCfg.GLM_artifacts = false;
end

figure;plot(Ch1_filtered)
data2use = Ch2_filtered;

%% remove channel artifacts
% Issue with small NIDAQ
if isfield(TrigCfg, 'Remove_artifacts')
    if TrigCfg.Remove_artifacts
        % Interpolate method (last resort)
        datavec_artifactremoved_ch1 = artifact_interpolate(TrigCfg, data, ch1_data_table);
    end
else
    TrigCfg.Remove_artifacts = false;
end

%% Window info
% Window info
prew_f = TrigCfg.prew * freq;
postw_f = TrigCfg.postw * freq;
l = prew_f + postw_f + 1;

% Naming info
mouse = defaultpath(29:end-1);
date = defaultpath(22:end-8);
%% Monkeylogic trials
% Replace with opto-only if needed
useopto = true;

% Grab the trial pulse info and snap it to the photometry pulses
if ~useopto
    trial = ML_pulse_table(:,2);
%     opto = opto_pulse_table(:,2);
else
    trial = opto_pulse_table(:,2);
end

    % For ML
    % Grab trial onsets
%     ML_ons = chainfinder(trial > 0.5);
% %     ML_ons=ML_ons(3:end,:);
%     ml.TimingFileByCond=1;
%     ml.ConditionNumber=1;
    %Add trial markers
%     ML_ons(:,2)=ml.ConditionNumber;
         
%     ML_ons_master = cell(1,length(ml.TimingFileByCond)); %Initialize cell array

    % Collect onsets for each condition into cell array 
%     for Condition = 1:length(ml.TimingFileByCond) 
%         ML_cdn=Condition+zeros(length(ML_ons),1);
%         ML_ons_temp(:,:)= ML_ons(ML_ons(:,2)==ML_cdn);
%         ML_ons_master{:,Condition} = ML_ons_temp;
%         ML_ons_temp=[];
%         n_trials(1,Condition) = length(ML_ons_master{1,Condition});
%     end
%  
%     % Remove trials too close too edges of experiment
%     for Condition = 1:length(ml.TimingFileByCond) 
%         badtrials = ((ML_ons_master{1,Condition} - prew_f) <= 0) + ((ML_ons_master{1,Condition} + postw_f) > n_points);
%         ML_ons_master{1,Condition}(badtrials > 0) = []; 
%         badtrials =[];
%         n_trials(1,Condition) = length(ML_ons_master{1,Condition});
%     end
%     
%     %For opto
    % Find bad opto pulses if needed
    if ~isempty(TrigCfg.minpulsewidth)
        % Get all pulses
        pulseinfo = chainfinder(trial>0.5);
        
        % Bad pulses
        badpulses = pulseinfo(pulseinfo(:,2) < TrigCfg.minpulsewidth, :);
        badpulses(:,2) = badpulses(:,1) + badpulses(:,2) - 1;
        
        % Clean up
        for i = 1 : size(badpulses, 1)
            data(TrigCfg.opto_channel, badpulses(i,1) : badpulses(i,2)) = 0;
        end
    end
    
    % Grab opto onsets
    opto_ons = chainfinder(trial > 0.5);

    % Grab opto inter-stim interval
    opto_isi = diff(opto_ons(:,1));
    opto_isi = [Inf; opto_isi];
    
    % Train lengths
    train_ons = find(opto_isi > TrigCfg.trainlength_threshold * freq);
    tl = opto_ons(train_ons(3)-1) - opto_ons(train_ons(2)) + 2;
    
    % Inter-train interval
    ITI = opto_ons(train_ons(3)) - opto_ons(train_ons(2));
    
    % Determine the actual onsets of trains
    opto_ons = opto_ons(opto_isi > TrigCfg.trainlength_threshold * freq);
    
    % Apply offset in debugging mode
    if TrigCfg.DebugMode
        opto_ons = opto_ons + TrigCfg.opto_on_offset * freq;
    end
    
    % See if any of the pulses is too close to the beginning or the end of the
    % session
    badstims = ((opto_ons - prew_f) <= 0) + ((opto_ons + postw_f) > n_points);
    opto_ons(badstims > 0) = [];
    
    % Number of stims
    n_optostims = length(opto_ons);

%% Flatten data
% Pull data

% Flatten if needed
if TrigCfg.flatten_data
    if TrigCfg.Remove_artifacts
        [data2use, ~, exp_fit, ~] = tcpUIflatten(datavec_artifactremoved, opto);
        data2use_unfilt = datavec_artifactremoved - exp_fit;
    else
        [data2use, ~, exp_fit, ~] = tcpUIflatten(data2use, opto);
        data2use_unfilt = ch1_data_table(:, 2) - exp_fit;
    end
else
    if TrigCfg.Remove_artifacts
        data2use_unfilt = datavec_artifactremoved;
    else
        data2use_unfilt = ch1_data_table(:, 2);
    end
end
figure;plot([data2use, trial])

%% For non-sliding dff

% %Flatten 405 channel
% data1use = Ch1_filtered;
% opto = opto_pulse_table(:,2);
% [data1use, ~, exp_fit, ~] = tcpUIflatten(data1use, opto);
% data1use_unfilt = ch1_data_table(:, 2) - exp_fit;
% 
% data2use = Ch2_filtered;
% Baselines = mean(data2use(1:10000));  %change baseline window here
% 
% speed_binary = logical(speed_upsampled);
% for q = 1:size(data2use,1);
%     df(q,:) = data2use(q,:) - Baselines;
% end
% for q = 1:size(data2use,1);
%     dff(q,:) = df(q,:) / Baselines(:);
% end
% dff = dff*100;
% data1use=data1use*10;
% figure;plot((1 : n_points)'/freq,data1use);
% set(gcf,'color','w')
% hold on; plot((1 : n_points)'/freq, dff);
% ylim([-4 10])
% ylabel('dff (%)')
% yyaxis right
% ylim([-1 20])
% yticks([])
% ylabel('Running bouts')
% hold on; plot((1:n_points)'/freq,speed_binary);
% % set(gca,'Xticklabel',[])
% lgd = legend('Location', 'northwest');
%     legend('boxoff');
%     lgd.ItemTokenSize = [30,18];
% legend('Ch1:rPBN-405','Ch2:lPBN-465');
% names={'0','5','10','15','20','25'};
% set(gca,'xtick',[0,300,600,900,1200,1500],'xticklabel',names);
% xlabel('Time (min)')
% set(gca,'FontSize',18,'box','off');
%% Sliding window dff data
% Dff data if needed
if TrigCfg.dff_data
    % Pull data
    if TrigCfg.Remove_artifacts
        data2use = tcpPercentiledff(datavec_artifactremoved, freq, TrigCfg.dff_win, TrigCfg.dff_prc);
        data2use_unfilt = data2use;
    else
        data2use = tcpPercentiledff(data2use, freq, TrigCfg.dff_win, TrigCfg.dff_prc);
        data2use_unfilt = tcpPercentiledff(ch1_data_table(:, 2), freq, TrigCfg.dff_win, TrigCfg.dff_prc);
    end
    exp_fit = [];
    plot([data2use, trial])
end
figure;plot([data2use, trial])
%% Grab the point indices
% Indices
% Initialize variables
inds_master = zeros(n_optostims, 2); % Initialize indices for the condition
trigmat = zeros(l, n_optostims);     % Triggered matrix for the condition
trigmat_baseline = zeros(1, n_optostims); % Baseline values
trigmat_sub = zeros(l, n_optostims);      % Baseline-subtracted matrix

% Calculate indices for the single condition
inds_master(:,1) = opto_ons - prew_f;           % Start indices
inds_master(:,2) = opto_ons + postw_f;          % End indices

% Populate triggered matrix
for i = 1:n_optostims
    trigmat(:,i) = data2use(inds_master(i,1):inds_master(i,2));
end

% Calculate baseline and subtract from the triggered matrix
for i = 1:n_optostims
    trigmat_baseline(1,i) = mean(trigmat(1:prew_f,i)); % Baseline (pre-event)
    trigmat_sub(:,i) = trigmat(:,i) - trigmat_baseline(1,i); % Baseline subtraction
end

% Calculate average triggered result
trigmat_avg = mean(trigmat_sub, 2, 'omitnan'); % Average across trials
%% Deal with motion
% Check if the running file is there
runningfn = sprintf('%srunning.mat', filename(1:end-22));
runningfn_full = fullfile(filepath, runningfn);

if exist(runningfn_full, 'file')
    % Load running data
    running = load(runningfn_full, 'speed');
    
    % Running running sample count
    nrunpulse = size(chainfinder(data(5,:)>0.5),1);
    nrunlength = length(running.speed);
    if nrunpulse ~= nrunlength
        % Say something
        fprintf('Running digitization is %0.3f%% off\n', (1 - nrunlength/nrunpulse)*100);
        
        % Upsample running data
        speed_upsampled0 = TDresamp(running.speed', 'resample', nrunpulse/nrunlength * 0.9974);
        speed_upsampled = TDresamp(speed_upsampled0, 'resample',...
            n_points/nrunpulse);
    else
        % Upsample running data
        speed_upsampled = TDresamp(running.speed', 'resample',...
            n_points/length(running.speed));
    end
    
    % Fix the number of points if needed
    if length(speed_upsampled) > n_points
        speed_upsampled = speed_upsampled(1:n_points);
    elseif length(speed_upsampled) < n_points
        speed_upsampled(end:end + n_points - length(speed_upsampled)) = 0;
    end
    
    % Initialize a triggered speed matrix
    speedmat = zeros(l, n_optostims);  % Triggered speed matrix for one condition
    speedmat_avg = zeros(l, 1);     % Average triggered speed result for one condition

    % Populate the speed matrix
    for i = 1:n_optostims
        speedmat(:,i) = speed_upsampled(inds_master(i,1):inds_master(i,2));
    end

    % Calculate the average triggered results
    speedmat_avg = mean(speedmat, 2, 'omitnan'); % Average across trials
else
    % Store empty speed matrices
    speedmat = [];
    speedmat_avg = [];
end
%% Deal with licking
% Initialize a triggered lick matrix
% C = [data(2,:)+data(7,:)]; % to deal with switcher freq
% C = [data(2,:)];
% C = C-min(data(2));
% lickvec = tcpDatasnapper(data(TrigCfg.lickch,:)', C', 'max', 'pulsetopulse');
% lickvec = lickvec(:,2);
% lickmat = cell(1,length(ml.TimingFileByCond));
% for Condition = 1:length(ml.TimingFileByCond)
%     for i = 1 : n_trials(1,Condition)
%         lickmat{1,Condition}(1:l,i)=zeros(1,1);
%     end
%     %populate matrix
%     for i = 1 : n_trials(1,Condition)
%        lickmat{1,Condition}(:,i) = lickvec(inds_master{:,Condition}(i,1) : inds_master{:,Condition}(i,2));
%        lickmat{1,Condition}=round(lickmat{1,Condition});
%        lickmat{1,Condition}(:,:)=logical(lickmat{1,Condition});
%     end
% end
% 
% % Calculate the average triggered results
% for Condition = 1:length(ml.TimingFileByCond)
%     for i = 1 : n_trials(1,Condition)
%         lickmat_sum(:,Condition) = sum(lickmat{1,Condition}(:,:),2);
%     end
% end
%%
% Define Y-axis limits
y_ax = [-0.1 0.1];
speed_y_ax = [0 max(speedmat(:)) * 1.1]; % Define Y-axis limits for speed

% Initialize figure
f1 = figure(); 
f1.Position = [50 0 1200 1000]; % Adjust figure height for three rows
set(gcf, 'Color', 'w');
hold on;

% Time axis
seconds = -TrigCfg.prew : 1 / freq : TrigCfg.postw;

% Plot trigmat_sub (heatmap)
subplot(3, 1, 1); % Top subplot
imagesc(seconds, [1 size(trigmat_sub, 2)], trigmat_sub'); % Transpose for correct orientation
hold on;
line([0 0], [0.5 size(trigmat_sub, 2) + 0.5], 'linestyle', '--', 'linewidth', 1.5, 'color', 'k');
title('Heatmap of \DeltaF/F_0');
xlabel('Time (s)');
ylabel('Trials');
colormap(flipud(brewermap([], 'RdBu'))); % Use RdBu colormap
cbar = colorbar; % Add a colorbar
cbar.Ticks = linspace(y_ax(1), y_ax(2), 5); % Set colorbar ticks
cbar.Label.String = '\DeltaF/F_0'; % Set colorbar label
cbar.Label.FontSize = 12; % Set label font size
caxis(y_ax); % Set color axis limits
set(gca, 'fontsize', 14);

% Plot trigmat_avg as a line plot
subplot(3, 1, 2); % Middle subplot
shadedErrorBar(seconds, mean(trigmat_sub, 2), std(trigmat_sub, 0, 2) / sqrt(size(trigmat_sub, 2)), 'lineProps', 'k'); % Mean +/- SEM
hold on;
plot(seconds, trigmat_avg, 'r-', 'LineWidth', 1.5); % Plot average with red line
line([0 0], y_ax, 'linestyle', '--', 'linewidth', 1.5, 'color', 'k'); % Add vertical line at t=0
title('Average \DeltaF/F_0');
xlabel('Time (s)');
ylabel('\DeltaF/F_0');
ylim(y_ax); % Set Y-axis limits
set(gca, 'fontsize', 14);

% Plot speedmat as a line plot
subplot(3, 1, 3); % Bottom subplot
shadedErrorBar(seconds, mean(speedmat, 2), std(speedmat, 0, 2) / sqrt(size(speedmat, 2)), 'lineProps', 'b'); % Mean +/- SEM for speed
hold on;
line([0 0], speed_y_ax, 'linestyle', '--', 'linewidth', 1.5, 'color', 'k'); % Add vertical line at t=0
title('Average Speed');
xlabel('Time (s)');
ylabel('Speed (cm/s)');
ylim(speed_y_ax); % Set Y-axis limits for speed
set(gca, 'fontsize', 14);

% Adjust layout
sgtitle('Triggered Matrices and Speed for One Condition'); % Add a title for the entire figure

% Save figure
% mouse = defaultpath(29:end-8);
% date = defaultpath(36:end-1);
% filen = strcat(date, '_', mouse, '_Trigmat_Plots'); 
% saveas(gcf, fullfile(defaultpath, [filen, '.png']));


%% Save results
if TrigCfg.flatten_data
    save(fullfile(filepath,filename_output_triggered), 'TrigCfg', 'trigmat',...
        'trigmat_sub', 'trigmat_avg',...
        'freq', 'prew_f', 'postw_f', 'l', 'inds_master',...
        'trigmat_avg', 'data2use' , 'data2use_unfilt', 'exp_fit',...
        'speedmat', 'speedmat_avg'  ,'mouse','date');
else
    save(fullfile(filepath,filename_output_triggered), 'TrigCfg', 'trigmat',...
        'trigmat_sub', 'trigmat_avg',...
        'freq', 'prew_f', 'postw_f', 'l', 'inds_master',...
        'trigmat_avg', 'data2use' , 'data2use_unfilt', ...
        'speedmat', 'speedmat_avg','mouse','date');
end