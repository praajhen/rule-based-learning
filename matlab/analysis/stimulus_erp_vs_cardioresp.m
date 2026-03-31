clear; clc; close all;

subj = 'P02';

addpath('C:\MyTemp\fieldtrip-20251218');
ft_defaults;

dataset = ['C:\MyTemp\RBLT_project\pilot\RBL_' subj '.edf'];

%% ================= LOAD =================
hdr   = ft_read_header(dataset);
event = ft_read_event(dataset);
event(arrayfun(@(x) isempty(x.value), event)) = [];

codes   = nan(length(event),1);
samples = nan(length(event),1);

for i = 1:length(event)
    num = regexp(event(i).value,'\d+','match');
    if ~isempty(num)
        codes(i)=str2double(num{end});
        samples(i)=event(i).sample;
    end
end

valid = ~isnan(codes);
codes = codes(valid);
samples = samples(valid);

fs = hdr.Fs;

%% ================= STIMULUS =================
stim_codes = [101:109 131:139];
stim_idx = find(ismember(codes,stim_codes));

prestim  = 0.2;
poststim = 0.8;

trl = zeros(length(stim_idx),3);

for i=1:length(stim_idx)
    s = samples(stim_idx(i));
    trl(i,1)=round(s-prestim*fs);
    trl(i,2)=round(s+poststim*fs);
    trl(i,3)=round(-prestim*fs);
end

%% ================= LOAD RAW =================
cfg=[];
cfg.dataset=dataset;
raw = ft_preprocessing(cfg);

labels = lower(raw.label);

ecg_chan  = find(contains(labels,'ecg'),1);
resp_chan = find(contains(labels,'resp'),1);

ecg  = raw.trial{1}(ecg_chan,:);
resp = raw.trial{1}(resp_chan,:);

%% ================= ECG =================
ecg_filt = bandpass(ecg,[5 20],fs);

[pks,locs] = findpeaks(ecg_filt,...
    'MinPeakHeight',mean(ecg_filt)+1.2*std(ecg_filt),...
    'MinPeakDistance',0.4*fs);

cardiac = zeros(length(stim_idx),1);

for i=1:length(stim_idx)

    stim_sample = samples(stim_idx(i));

    prev = find(locs<stim_sample,1,'last');
    next = find(locs>stim_sample,1,'first');

    if ~isempty(prev) && ~isempty(next)

        phase=(stim_sample-locs(prev))/(locs(next)-locs(prev));

        if phase<0.5
            cardiac(i)=1; % systole
        else
            cardiac(i)=2; % diastole
        end
    end
end

%% ================= RESP =================
resp_filt = bandpass(resp,[0.1 0.5],fs);
resp_phase = angle(hilbert(resp_filt));

resp_cond = zeros(length(stim_idx),1);

for i=1:length(stim_idx)

    stim_sample = samples(stim_idx(i));

    if stim_sample <= length(resp_phase)

        if resp_phase(stim_sample)>0
            resp_cond(i)=1; % inspiration
        else
            resp_cond(i)=2; % expiration
        end

    end
end

%% ================= EEG =================
badch = {'ECG','EMG','Resp','EDA'};
eeg_channels = setdiff(hdr.label,badch);

cfg=[];
cfg.dataset=dataset;
cfg.trl=trl;
cfg.channel=eeg_channels;

cfg.demean='yes';
cfg.baselinewindow=[-0.2 0];

cfg.hpfilter = 'yes'; cfg.hpfreq = 0.5; cfg.hpfilttype = 'firws';
cfg.lpfilter = 'yes'; cfg.lpfreq = 30;  cfg.lpfilttype = 'firws';

data = ft_preprocessing(cfg);

%% ================= ELECTRODES =================
occ = {'52','65','70','72','75','83','90','92'};
par = {'62','72','75'};

occ_idx = find(ismember(data.label,occ));
par_idx = find(ismember(data.label,par));

t = data.time{1};

n1_win = t>=0.09 & t<=0.16;
p3_win = t>=0.25 & t<=0.40;

%% ================= COLORS =================
blue = [0 0.4 1];
red  = [1 0.2 0.2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ================= CARDIAC ERP ==========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sys_trials = find(cardiac==1);
dia_trials = find(cardiac==2);

sys_occ=[];
dia_occ=[];

sys_par=[];
dia_par=[];

for i=1:length(sys_trials)

    sys_occ(:,i)=mean(data.trial{sys_trials(i)}(occ_idx,:),1);
    sys_par(:,i)=mean(data.trial{sys_trials(i)}(par_idx,:),1);

end

for i=1:length(dia_trials)

    dia_occ(:,i)=mean(data.trial{dia_trials(i)}(occ_idx,:),1);
    dia_par(:,i)=mean(data.trial{dia_trials(i)}(par_idx,:),1);

end

%% means
sys_m = mean(sys_occ,2);
dia_m = mean(dia_occ,2);

sys_sem = std(sys_occ,[],2)/sqrt(size(sys_occ,2));
dia_sem = std(dia_occ,[],2)/sqrt(size(dia_occ,2));

%% amplitudes
N1_sys = mean(sys_occ(n1_win,:),1);
N1_dia = mean(dia_occ(n1_win,:),1);

P3_sys = mean(sys_par(p3_win,:),1);
P3_dia = mean(dia_par(p3_win,:),1);

%% ================= FIGURE CARDIAC =================
figure('Color','w','Position',[200 300 1000 400])

%% OCCIPITAL N1
subplot(1,2,1); hold on

fill([t fliplr(t)],...
 [sys_m'-sys_sem' fliplr(sys_m'+sys_sem')],...
 blue,'FaceAlpha',0.25,'EdgeColor','none')

fill([t fliplr(t)],...
 [dia_m'-dia_sem' fliplr(dia_m'+dia_sem')],...
 red,'FaceAlpha',0.25,'EdgeColor','none')

plot(t,sys_m,'Color',blue,'LineWidth',2)
plot(t,dia_m,'Color',red,'LineWidth',2)

legend({'Systole','Diastole'})
title('Occipital ERP (N1)')
xlabel('Time (s)')
ylabel('Amplitude (µV)')
box off

%% PARIETAL P3
subplot(1,2,2); hold on

sys_par_m = mean(sys_par,2);
dia_par_m = mean(dia_par,2);

sys_par_sem = std(sys_par,[],2)/sqrt(size(sys_par,2));
dia_par_sem = std(dia_par,[],2)/sqrt(size(dia_par,2));

fill([t fliplr(t)],...
 [sys_par_m'-sys_par_sem' fliplr(sys_par_m'+sys_par_sem')],...
 blue,'FaceAlpha',0.25,'EdgeColor','none')

fill([t fliplr(t)],...
 [dia_par_m'-dia_par_sem' fliplr(dia_par_m'+dia_par_sem')],...
 red,'FaceAlpha',0.25,'EdgeColor','none')

plot(t,sys_par_m,'Color',blue,'LineWidth',2)
plot(t,dia_par_m,'Color',red,'LineWidth',2)

legend({'Systole','Diastole'})
title('Parietal ERP (P3)')
xlabel('Time (s)')
ylabel('Amplitude (µV)')
box off

exportgraphics(gcf,...
['ECG_vs_ERP ''.png'],...
'Resolution',300)


%% ================= FIGURE RESP =================
figure('Color','w','Position',[200 300 1000 400])

%% OCCIPITAL N1
subplot(1,2,1); hold on

fill([t fliplr(t)],...
 [insp_m'-insp_sem' fliplr(insp_m'+insp_sem')],...
 blue,'FaceAlpha',0.25,'EdgeColor','none')

fill([t fliplr(t)],...
 [exp_m'-exp_sem' fliplr(exp_m'+exp_sem')],...
 red,'FaceAlpha',0.25,'EdgeColor','none')

plot(t,insp_m,'Color',blue,'LineWidth',2)
plot(t,exp_m,'Color',red,'LineWidth',2)

legend({'Inspiration','Expiration'})
title('Occipital ERP (N1)')
xlabel('Time (s)')
ylabel('Amplitude (µV)')
box off

%% PARIETAL P3
subplot(1,2,2); hold on

insp_par=[];
exp_par=[];

for i=1:length(insp_trials)
    insp_par(:,i)=mean(data.trial{insp_trials(i)}(par_idx,:),1);
end

for i=1:length(exp_trials)
    exp_par(:,i)=mean(data.trial{exp_trials(i)}(par_idx,:),1);
end

insp_par_m = mean(insp_par,2);
exp_par_m  = mean(exp_par,2);

insp_par_sem = std(insp_par,[],2)/sqrt(size(insp_par,2));
exp_par_sem  = std(exp_par,[],2)/sqrt(size(exp_par,2));

fill([t fliplr(t)],...
 [insp_par_m'-insp_par_sem' fliplr(insp_par_m'+insp_par_sem')],...
 blue,'FaceAlpha',0.25,'EdgeColor','none')

fill([t fliplr(t)],...
 [exp_par_m'-exp_par_sem' fliplr(exp_par_m'+exp_par_sem')],...
 red,'FaceAlpha',0.25,'EdgeColor','none')

plot(t,insp_par_m,'Color',blue,'LineWidth',2)
plot(t,exp_par_m,'Color',red,'LineWidth',2)

legend({'Inspiration','Expiration'})
title('Parietal ERP (P3)')
xlabel('Time (s)')
ylabel('Amplitude (µV)')
box off

exportgraphics(gcf,...
['RESP_vs_ERP ''.png'],...
'Resolution',300)
