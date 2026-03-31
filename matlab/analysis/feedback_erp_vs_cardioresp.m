clear; clc; close all;

%% ================= SUBJECT =================
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

%% ================= PROBABILISTIC FEEDBACK =================
hard_stim = 161:169;
stim_idx = find(ismember(codes,hard_stim));

fb_codes = [180 181];

fb_samples = nan(length(stim_idx),1);
fb_types   = nan(length(stim_idx),1);

for i = 1:length(stim_idx)

    for j = stim_idx(i)+1 : min(stim_idx(i)+20,length(codes))

        if ismember(codes(j),fb_codes)

            fb_samples(i) = samples(j);
            fb_types(i)   = codes(j);
            break

        end
    end
end

valid = ~isnan(fb_samples);
fb_samples = fb_samples(valid);
fb_types   = fb_types(valid);

%% ================= LOAD RAW (ECG + RESP) =================
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

cardiac = zeros(length(fb_samples),1);

for i=1:length(fb_samples)

    s = fb_samples(i);

    prev = find(locs<s,1,'last');
    next = find(locs>s,1,'first');

    if ~isempty(prev)&&~isempty(next)

        phase=(s-locs(prev))/(locs(next)-locs(prev));

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

resp_cond = zeros(length(fb_samples),1);

for i=1:length(fb_samples)

    s = fb_samples(i);

    if s <= length(resp_phase)

        if resp_phase(s)>0
            resp_cond(i)=1; % inspiration
        else
            resp_cond(i)=2; % expiration
        end
    end
end

%% ================= TRIALS =================
prestim  = 0.2;
poststim = 0.8;

trl = zeros(length(fb_samples),3);

for i=1:length(fb_samples)

    s = fb_samples(i);

    trl(i,1)=round(s-prestim*fs);
    trl(i,2)=round(s+poststim*fs);
    trl(i,3)=round(-prestim*fs);

end

%% ================= CONDITIONS =================
cond = zeros(length(fb_types),1);
cond(fb_types==180)=1; % correct
cond(fb_types==181)=2; % incorrect

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
data.trialinfo = cond;

%% ================= ELECTRODES =================
fc = {'6','7','31','55','106'};      % FRN
par = {'61','62','67','72','78'};    % P3

fc_idx = find(ismember(data.label,fc));
par_idx = find(ismember(data.label,par));

t = data.time{1};

blue=[0 0.4 1];
red =[1 0.2 0.2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ================= CARDIAC ==============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',[200 300 1000 400])

for phase = 1:2

    if phase==1
        idx = find(cardiac==1);
        title_txt = 'Systole';
    else
        idx = find(cardiac==2);
        title_txt = 'Diastole';
    end

end

corr = find(cond==1 & cardiac==1);
err  = find(cond==2 & cardiac==1);

% correct systole
corr_fc=[];
corr_par=[];
for i=1:length(corr)
corr_fc(:,i)=mean(data.trial{corr(i)}(fc_idx,:),1);
corr_par(:,i)=mean(data.trial{corr(i)}(par_idx,:),1);
end

% incorrect systole
err_fc=[];
err_par=[];
for i=1:length(err)
err_fc(:,i)=mean(data.trial{err(i)}(fc_idx,:),1);
err_par(:,i)=mean(data.trial{err(i)}(par_idx,:),1);
end

corr_fc_m = mean(corr_fc,2);
err_fc_m  = mean(err_fc,2);

corr_fc_sem = std(corr_fc,[],2)/sqrt(size(corr_fc,2));
err_fc_sem  = std(err_fc,[],2)/sqrt(size(err_fc,2));

subplot(1,2,1); hold on

fill([t fliplr(t)],...
 [corr_fc_m'-corr_fc_sem' fliplr(corr_fc_m'+corr_fc_sem')],...
 blue,'FaceAlpha',0.25,'EdgeColor','none')

fill([t fliplr(t)],...
 [err_fc_m'-err_fc_sem' fliplr(err_fc_m'+err_fc_sem')],...
 red,'FaceAlpha',0.25,'EdgeColor','none')

plot(t,corr_fc_m,'Color',blue,'LineWidth',2)
plot(t,err_fc_m,'Color',red,'LineWidth',2)

legend({'Correct','Incorrect'})
title('FRN (systole)')
xlabel('Time (s)')
ylabel('µV')
box off

corr_par_m = mean(corr_par,2);
err_par_m  = mean(err_par,2);

corr_par_sem = std(corr_par,[],2)/sqrt(size(corr_par,2));
err_par_sem  = std(err_par,[],2)/sqrt(size(err_par,2));

subplot(1,2,2); hold on

fill([t fliplr(t)],...
 [corr_par_m'-corr_par_sem' fliplr(corr_par_m'+corr_par_sem')],...
 blue,'FaceAlpha',0.25,'EdgeColor','none')

fill([t fliplr(t)],...
 [err_par_m'-err_par_sem' fliplr(err_par_m'+err_par_sem')],...
 red,'FaceAlpha',0.25,'EdgeColor','none')

plot(t,corr_par_m,'Color',blue,'LineWidth',2)
plot(t,err_par_m,'Color',red,'LineWidth',2)

legend({'Correct','Incorrect'})
title('Feedback P3 (systole)')
xlabel('Time (s)')
ylabel('µV')
box off

exportgraphics(gcf,'feedback_ecg_vs_erp.png','Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ================= RESP =================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',[200 300 1000 400])

corr = find(cond==1 & resp_cond==1);
err  = find(cond==2 & resp_cond==1);

corr_fc=[];
corr_par=[];
for i=1:length(corr)
corr_fc(:,i)=mean(data.trial{corr(i)}(fc_idx,:),1);
corr_par(:,i)=mean(data.trial{corr(i)}(par_idx,:),1);
end

err_fc=[];
err_par=[];
for i=1:length(err)
err_fc(:,i)=mean(data.trial{err(i)}(fc_idx,:),1);
err_par(:,i)=mean(data.trial{err(i)}(par_idx,:),1);
end

corr_fc_m = mean(corr_fc,2);
err_fc_m  = mean(err_fc,2);

corr_fc_sem = std(corr_fc,[],2)/sqrt(size(corr_fc,2));
err_fc_sem  = std(err_fc,[],2)/sqrt(size(err_fc,2));

subplot(1,2,1); hold on

fill([t fliplr(t)],...
 [corr_fc_m'-corr_fc_sem' fliplr(corr_fc_m'+corr_fc_sem')],...
 blue,'FaceAlpha',0.25,'EdgeColor','none')

fill([t fliplr(t)],...
 [err_fc_m'-err_fc_sem' fliplr(err_fc_m'+err_fc_sem')],...
 red,'FaceAlpha',0.25,'EdgeColor','none')

plot(t,corr_fc_m,'Color',blue,'LineWidth',2)
plot(t,err_fc_m,'Color',red,'LineWidth',2)

legend({'Correct','Incorrect'})
title('FRN (inspiration)')
xlabel('Time (s)')
ylabel('µV')
box off

corr_par_m = mean(corr_par,2);
err_par_m  = mean(err_par,2);

corr_par_sem = std(corr_par,[],2)/sqrt(size(corr_par,2));
err_par_sem  = std(err_par,[],2)/sqrt(size(err_par,2));

subplot(1,2,2); hold on

fill([t fliplr(t)],...
 [corr_par_m'-corr_par_sem' fliplr(corr_par_m'+corr_par_sem')],...
 blue,'FaceAlpha',0.25,'EdgeColor','none')

fill([t fliplr(t)],...
 [err_par_m'-err_par_sem' fliplr(err_par_m'+err_par_sem')],...
 red,'FaceAlpha',0.25,'EdgeColor','none')

plot(t,corr_par_m,'Color',blue,'LineWidth',2)
plot(t,err_par_m,'Color',red,'LineWidth',2)

legend({'Correct','Incorrect'})
title('Feedback P3 (inspiration)')
xlabel('Time (s)')
ylabel('µV')
box off

exportgraphics(gcf,'feedback_resp_vs_erp.png','Resolution',300)