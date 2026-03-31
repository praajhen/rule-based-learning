clear; clc; close all;

%% ================== SUBJECT ==================
subj = 'P02';

%% ================== SETUP ==================
addpath('C:\MyTemp\fieldtrip-20251218');
ft_defaults;

dataset = ['C:\MyTemp\RBLT_project\pilot\RBL_' subj '.edf'];

%% ================== LOAD ==================
hdr   = ft_read_header(dataset);
event = ft_read_event(dataset);
event(arrayfun(@(x) isempty(x.value), event)) = [];

%% ================== EXTRACT CODES ==================
codes   = nan(length(event),1);
samples = nan(length(event),1);

for i = 1:length(event)
    num = regexp(event(i).value,'\d+','match');
    if ~isempty(num)
        codes(i)   = str2double(num{end});
        samples(i) = event(i).sample;
    end
end

valid = ~isnan(codes);
codes   = codes(valid);
samples = samples(valid);

%% ================== PROBABILISTIC STIM ==================
hard_stim = 161:169;
stim_idx = find(ismember(codes, hard_stim));

%% ================== FIND FEEDBACK ==================
fb_codes = [180 181];

fb_samples = nan(length(stim_idx),1);
fb_types   = nan(length(stim_idx),1);

for i = 1:length(stim_idx)
    for j = stim_idx(i)+1 : min(stim_idx(i)+20, length(codes))
        if ismember(codes(j), fb_codes)
            fb_samples(i) = samples(j);
            fb_types(i)   = codes(j);
            break
        end
    end
end

valid = ~isnan(fb_samples);
fb_samples = fb_samples(valid);
fb_types   = fb_types(valid);

fprintf('Correct: %d  | Incorrect: %d\n', ...
    sum(fb_types==180), sum(fb_types==181))

%% ================== TRIALS ==================
fs = hdr.Fs;

prestim  = 0.2;
poststim = 0.8;

trl = [];

for i = 1:length(fb_samples)
    s = fb_samples(i);
    trl = [trl;
        round(s - prestim*fs), ...
        round(s + poststim*fs), ...
        round(-prestim*fs)];
end

%% ================== CONDITIONS ==================
cond = zeros(length(fb_types),1);
cond(fb_types==180) = 1; % correct
cond(fb_types==181) = 2; % incorrect

%% ================== PREPROCESS ==================
badch = {'ECG','EMG','Resp','EDA'};
eeg_channels = setdiff(hdr.label, badch);

cfg = [];
cfg.dataset = dataset;
cfg.trl     = trl;
cfg.channel = eeg_channels;

cfg.demean = 'yes';
cfg.baselinewindow = [-0.2 0];

cfg.hpfilter = 'yes'; cfg.hpfreq = 0.5; cfg.hpfilttype = 'firws';
cfg.lpfilter = 'yes'; cfg.lpfreq = 30;  cfg.lpfilttype = 'firws';

data = ft_preprocessing(cfg);
data.trialinfo = cond;

%% ================== PARIETAL CHANNELS ==================
par_channels = {'62','72','75'};
par_idx = find(ismember(data.label, par_channels));

%% ================== GRAND AVERAGE ==================
correct_trials = find(data.trialinfo==1);
error_trials   = find(data.trialinfo==2);

correct = [];
for i = 1:length(correct_trials)
    correct(:,i) = mean(data.trial{correct_trials(i)}(par_idx,:),1);
end

error = [];
for i = 1:length(error_trials)
    error(:,i) = mean(data.trial{error_trials(i)}(par_idx,:),1);
end

t = data.time{1};

correct_mean = mean(correct,2);
error_mean   = mean(error,2);

correct_sem = std(correct,[],2)/sqrt(size(correct,2));
error_sem   = std(error,[],2)/sqrt(size(error,2));

%% ================== FEEDBACK P3 ==================
p3_idx = t >= 0.30 & t <= 0.50;

P3_correct = mean(correct_mean(p3_idx));
P3_error   = mean(error_mean(p3_idx));

P3_correct_sem = std(mean(correct(p3_idx,:),1)) / sqrt(size(correct,2));
P3_error_sem   = std(mean(error(p3_idx,:),1)) / sqrt(size(error,2));

fprintf('P3 Correct: %.2f\n',P3_correct)
fprintf('P3 Error: %.2f\n',P3_error)

%% ================== FIGURE ==================
figure('Color','w','Position',[200 300 1000 350])

correct_color = [0 0.4 1];
error_color   = [1 0.2 0.2];

%% ================= ERP =================
subplot(1,2,1)
hold on

fill([t fliplr(t)],...
     [correct_mean'-correct_sem' fliplr(correct_mean'+correct_sem')],...
     correct_color,'FaceAlpha',0.25,'EdgeColor','none')

fill([t fliplr(t)],...
     [error_mean'-error_sem' fliplr(error_mean'+error_sem')],...
     error_color,'FaceAlpha',0.25,'EdgeColor','none')

h1 = plot(t,correct_mean,'Color',correct_color,'LineWidth',2);
h2 = plot(t,error_mean,'Color',error_color,'LineWidth',2);

xline(0,'--k')
yline(0,'--k')

legend([h1 h2],{'Correct','Incorrect'})
title('Feedback-locked Parietal ERP')
ylabel('Amplitude (µV)')
xlabel('Time (s)')
box off


%% ================= TRIAL-LEVEL P3 =================
subplot(1,2,2)

% trial-level P3
p3_correct_trials = mean(correct(p3_idx,:),1);
p3_error_trials   = mean(error(p3_idx,:),1);

data  = [p3_correct_trials'; p3_error_trials'];
group = [ones(length(p3_correct_trials),1); ...
         2*ones(length(p3_error_trials),1)];

boxplot(data,group,'Colors','k','Symbol','')

h = findobj(gca,'Tag','Box');

patch(get(h(2),'XData'),get(h(2),'YData'),correct_color,'FaceAlpha',0.35)
patch(get(h(1),'XData'),get(h(1),'YData'),error_color,'FaceAlpha',0.35)

set(gca,'XTickLabel',{'Correct','Incorrect'})
ylabel('P3 amplitude (µV)')
title('Feedback P3')
box off

%% SAVE
exportgraphics(gcf,['Probabilistic_Feedback_ERP_' subj '.png'],'Resolution',300)