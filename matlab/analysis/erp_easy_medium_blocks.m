clear; clc; close all;

%% ================== SETUP ==================
addpath('C:\MyTemp\fieldtrip-20251218');
ft_defaults;

dataset = 'C:\MyTemp\RBLT_project\pilot\RBL_P02.edf';

%% ================== HEADER ==================
hdr = ft_read_header(dataset);
fs  = hdr.Fs;

%% ================== EVENTS ==================
event = ft_read_event(dataset);
event(arrayfun(@(x) isempty(x.value), event)) = [];

codes   = nan(length(event),1);
samples = nan(length(event),1);

for i = 1:length(event)
    num = regexp(event(i).value,'\d+','match');
    if ~isempty(num)
        codes(i)   = str2double(num{end});
        samples(i) = event(i).sample;
    end
end

valid   = ~isnan(codes);
codes   = codes(valid);
samples = samples(valid);

%% ================== STIMULI ==================
stim_codes = [101:109 131:139];
stim_idx   = find(ismember(codes, stim_codes));

prestim  = 0.2;
poststim = 0.8;

trl = [];

for i = 1:length(stim_idx)
    s = samples(stim_idx(i));
    trl = [trl; round(s-prestim*fs) round(s+poststim*fs) round(-prestim*fs)];
end

%% ================== CONDITION ==================
cond = zeros(length(stim_idx),1);

for i = 1:length(stim_idx)
    if ismember(codes(stim_idx(i)),101:109)
        cond(i)=1; % easy
    else
        cond(i)=2; % medium
    end
end

%% ================== CORRECT TRIALS ==================
correct_fb = [30 120 150];
all_fb     = [30 31 120 121 150 151];

is_correct = zeros(length(stim_idx),1);

for i = 1:length(stim_idx)
    for j = stim_idx(i)+1:min(stim_idx(i)+20,length(codes))
        if ismember(codes(j),all_fb)
            if ismember(codes(j),correct_fb)
                is_correct(i)=1;
            end
            break
        end
    end
end

%% ================== ALIGN ==================
min_len = min([length(cond), length(is_correct), size(trl,1)]);
cond = cond(1:min_len);
is_correct = is_correct(1:min_len);
trl = trl(1:min_len,:);

%% ================== PREPROCESS ==================
badch = {'ECG','EMG','Resp','EDA'};
eeg_channels = setdiff(hdr.label,badch);

cfg = [];
cfg.dataset = dataset;
cfg.trl = trl;
cfg.channel = eeg_channels;

cfg.demean = 'yes';
cfg.baselinewindow = [-0.2 0];

cfg.hpfilter = 'yes'; cfg.hpfreq = 0.5;cfg.hpfilttype = 'firws'; 
cfg.lpfilter = 'yes'; cfg.lpfreq = 30; cfg.lpfilttype = 'firws';

data = ft_preprocessing(cfg);

%% ================== SELECT CORRECT ==================
use_trials = find(is_correct==1);

data_sel = data;
data_sel.trial = data.trial(use_trials);
data_sel.time  = data.time(use_trials);
data_sel.trialinfo = cond(use_trials);

%% ================== ELECTRODES ==================
% OCCIPITAL visual N1
n1_channels = {'52','65','70','72','75','83','90','92'};

% PARIETAL P3
p3_channels = {'62','72','75'};

%% ================== WINDOWS ==================
n1_win = [0.12 0.18];
p3_win = [0.25 0.40];

time = data_sel.time{1};

n1_idx = time>=n1_win(1) & time<=n1_win(2);
p3_idx = time>=p3_win(1) & time<=p3_win(2);

n1_idx_chan = find(ismember(data_sel.label,n1_channels));
p3_idx_chan = find(ismember(data_sel.label,p3_channels));

%% ================== TRIAL LEVEL ==================
N1 = nan(length(data_sel.trial),1);
P3 = nan(length(data_sel.trial),1);

for i=1:length(data_sel.trial)

    trial_data = data_sel.trial{i};

    n1_signal = mean(trial_data(n1_idx_chan,:),1);
    p3_signal = mean(trial_data(p3_idx_chan,:),1);

    N1(i) = mean(n1_signal(n1_idx));
    P3(i) = mean(p3_signal(p3_idx));

end

%% ================== FIGURE ==================
figure('Color','w','Position',[200 200 1000 700])

easy_color   = [0 0.4 1];
medium_color = [1 0.2 0.2];

easy_trials   = find(data_sel.trialinfo==1);
medium_trials = find(data_sel.trialinfo==2);

t = data_sel.time{1};

%% ================= OCCIPITAL ERP (N1) =================
subplot(2,2,1)

easy = [];
for i = 1:length(easy_trials)
    easy(:,i) = mean(data_sel.trial{easy_trials(i)}(n1_idx_chan,:),1);
end

medium = [];
for i = 1:length(medium_trials)
    medium(:,i) = mean(data_sel.trial{medium_trials(i)}(n1_idx_chan,:),1);
end

easy_mean = mean(easy,2);
medium_mean = mean(medium,2);

easy_sem = std(easy,[],2)/sqrt(size(easy,2));
medium_sem = std(medium,[],2)/sqrt(size(medium,2));

hold on

fill([t fliplr(t)],...
     [easy_mean'-easy_sem' fliplr(easy_mean'+easy_sem')],...
     easy_color,'FaceAlpha',0.25,'EdgeColor','none')

fill([t fliplr(t)],...
     [medium_mean'-medium_sem' fliplr(medium_mean'+medium_sem')],...
     medium_color,'FaceAlpha',0.25,'EdgeColor','none')

h1 = plot(t,easy_mean,'Color',easy_color,'LineWidth',2);
h2 = plot(t,medium_mean,'Color',medium_color,'LineWidth',2);

xline(0,'--k')
yline(0,'--k')

title('Occipital ERP (N1)')
ylabel('Amplitude (µV)')
legend([h1 h2],{'Easy','Medium'})
box off


%% ================= PARIETAL ERP (P3) =================
subplot(2,2,2)

easy = [];
for i = 1:length(easy_trials)
    easy(:,i) = mean(data_sel.trial{easy_trials(i)}(p3_idx_chan,:),1);
end

medium = [];
for i = 1:length(medium_trials)
    medium(:,i) = mean(data_sel.trial{medium_trials(i)}(p3_idx_chan,:),1);
end

easy_mean = mean(easy,2);
medium_mean = mean(medium,2);

easy_sem = std(easy,[],2)/sqrt(size(easy,2));
medium_sem = std(medium,[],2)/sqrt(size(medium,2));

hold on

fill([t fliplr(t)],...
     [easy_mean'-easy_sem' fliplr(easy_mean'+easy_sem')],...
     easy_color,'FaceAlpha',0.25,'EdgeColor','none')

fill([t fliplr(t)],...
     [medium_mean'-medium_sem' fliplr(medium_mean'+medium_sem')],...
     medium_color,'FaceAlpha',0.25,'EdgeColor','none')

plot(t,easy_mean,'Color',easy_color,'LineWidth',2)
plot(t,medium_mean,'Color',medium_color,'LineWidth',2)

xline(0,'--k')
yline(0,'--k')

title('Parietal ERP (P3)')
box off


%% ================= N1 =================
subplot(2,2,3)

easy = N1(data_sel.trialinfo==1);
medium = N1(data_sel.trialinfo==2);

boxplot([easy;medium],...
        [ones(length(easy),1);2*ones(length(medium),1)],...
        'Colors','k','Symbol','')

h = findobj(gca,'Tag','Box');

patch(get(h(2),'XData'),get(h(2),'YData'),easy_color,'FaceAlpha',0.35)
patch(get(h(1),'XData'),get(h(1),'YData'),medium_color,'FaceAlpha',0.35)

set(gca,'XTickLabel',{'Easy','Medium'})
ylabel('N1 amplitude (µV)')
title('Occipital N1')
box off


%% ================= P3 =================
subplot(2,2,4)

easy = P3(data_sel.trialinfo==1);
medium = P3(data_sel.trialinfo==2);

boxplot([easy;medium],...
        [ones(length(easy),1);2*ones(length(medium),1)],...
        'Colors','k','Symbol','')

h = findobj(gca,'Tag','Box');

patch(get(h(2),'XData'),get(h(2),'YData'),easy_color,'FaceAlpha',0.35)
patch(get(h(1),'XData'),get(h(1),'YData'),medium_color,'FaceAlpha',0.35)

set(gca,'XTickLabel',{'Easy','Medium'})
ylabel('P3 amplitude (µV)')
title('Parietal P3')
box off

%% SAVE
exportgraphics(gcf,'Pilot_ERP_N1_P3.png','Resolution',300)