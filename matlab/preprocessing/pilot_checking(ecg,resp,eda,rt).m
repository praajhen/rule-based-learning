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

codes = nan(length(event),1);
times = nan(length(event),1);

for i = 1:length(event)
    val = event(i).value;
    num = regexp(val,'\d+','match');
    if ~isempty(num)
        codes(i) = str2double(num{end});
        times(i) = event(i).sample;
    end
end

valid_evt = ~isnan(codes);
codes = codes(valid_evt);
times = times(valid_evt);

%% ================== TRIGGERS ==================
stim_codes = [11:19 101:109 131:139 161:169];
resp_codes = [21:23 111:113 141:143 171:173];
raw_resp   = [1 2 3];
fb_codes   = [30 31 120 121 150 151 180 181];

%% ================== TRIAL EXTRACTION ==================
trial=[]; block=[]; stimulus=[]; response=[]; feedback=[];
RT=[]; stim_time=[]; fb_time=[];

t_idx=1;

for i = 1:length(codes)-3
    
    if ismember(codes(i),stim_codes)
        
        stim = codes(i);
        stim_t = times(i);
        
        resp=NaN; fb=NaN; rt=NaN; fb_t=NaN;
        
        % RT (first button press)
        for j=i+1:min(i+15,length(codes))
            if ismember(codes(j),raw_resp)
                rt = times(j) - stim_t;
                break
            end
        end
        
        % RESPONSE
        for j=i+1:min(i+15,length(codes))
            if ismember(codes(j),resp_codes)
                resp = codes(j);
                break
            end
        end
        
        % FEEDBACK
        for j=i+1:min(i+20,length(codes))
            if ismember(codes(j),fb_codes)
                fb = codes(j);
                fb_t = times(j);
                break
            end
        end
        
        trial(t_idx)=t_idx;
        stimulus(t_idx)=stim;
        response(t_idx)=resp;
        feedback(t_idx)=fb;
        RT(t_idx)=rt;
        stim_time(t_idx)=stim_t;
        fb_time(t_idx)=fb_t;
        
        % BLOCK
        if stim<=19
            block(t_idx)=0;
        elseif stim<=109
            block(t_idx)=1;
        elseif stim<=139
            block(t_idx)=2;
        else
            block(t_idx)=3;
        end
        
        t_idx=t_idx+1;
    end
end

%% ================== TABLE ==================
T = table(trial',block',stimulus',response',feedback',RT',...
    stim_time',fb_time',...
    'VariableNames',{'Trial','Block','Stimulus','Response','Feedback','RT_samples','StimTime','FbTime'});

T.RT = T.RT_samples / fs;
T.StimTimeSec = T.StimTime / fs;
T.FbTimeSec   = T.FbTime / fs;

T.Correct = ismember(T.Feedback,[30 120 150 180]);

% REMOVE INVALID
T = T(~isnan(T.FbTimeSec),:);

%% ================== RT FILTER ==================
fprintf('\nTrials BEFORE RT filter: %d\n', height(T))

T = T(T.RT > 0.2 & T.RT < 1.5, :);

fprintf('Trials AFTER RT filter: %d\n', height(T))

%% ================== LOAD DATA ==================
cfg=[]; cfg.dataset=dataset;
data=ft_preprocessing(cfg);

t = data.time{1};
labels = lower(data.label);

ecg_chan  = find(contains(labels,'ecg') | contains(labels,'ekg'),1);
resp_chan = find(contains(labels,'resp'),1);
eda_chan  = find(contains(labels,'eda') | contains(labels,'gsr'),1);

ecg  = data.trial{1}(ecg_chan,:);
resp = data.trial{1}(resp_chan,:);
eda  = data.trial{1}(eda_chan,:);

%% ================== FILTER ==================
ecg_filt  = bandpass(ecg,[5 20],fs);
resp_filt = bandpass(resp,[0.1 0.5],fs);
eda_smooth = smoothdata(eda,'gaussian',round(fs*0.5));

%% ================== ECG PEAKS ==================
[pks,locs] = findpeaks(ecg_filt,...
    'MinPeakHeight',mean(ecg_filt)+1.2*std(ecg_filt),...
    'MinPeakDistance',0.4*fs);

r_times = t(locs);

%% ================== CARDIAC PHASE (IMPROVED) ==================
CardiacFB = nan(height(T),1);

for i=1:height(T)
    
    fb = T.FbTimeSec(i);
    
    prev_idx = find(r_times < fb,1,'last');
    next_idx = find(r_times > fb,1,'first');
    
    if ~isempty(prev_idx) && ~isempty(next_idx)
        
        t1 = r_times(prev_idx);
        t2 = r_times(next_idx);
        
        phase = (fb - t1) / (t2 - t1); % 0 → 1
        
        if phase < 0.5
            CardiacFB(i) = 1; % systole (early)
        else
            CardiacFB(i) = 2; % diastole (late)
        end
    end
end

T.CardiacFB = CardiacFB;

%% ================== RESP PHASE (FIXED) ==================
resp_phase = angle(hilbert(resp_filt));

RespFB = nan(height(T),1);

for i=1:height(T)
    [~,idx] = min(abs(t - T.FbTimeSec(i)));
    
    if resp_phase(idx) > 0
        RespFB(i) = 1; % inspiration
    else
        RespFB(i) = 2; % expiration
    end
end

T.RespFB = RespFB;
%% ================== EDA ==================
EDA_fb = nan(height(T),1);

for i=1:height(T)
    
    fb = T.FbTimeSec(i);
    
    idx_resp = t>=fb & t<=fb+3;
    idx_base = t>=fb-1 & t<fb;
    
    if sum(idx_resp)>10 && sum(idx_base)>10
        baseline = mean(eda_smooth(idx_base));
        EDA_fb(i) = max(eda_smooth(idx_resp)) - baseline;
    end
end

T.EDA_FB = EDA_fb;

%% ================== SUMMARY ==================
fprintf('\n========== SUMMARY ==========\n')
fprintf('Trials: %d\n',height(T))
fprintf('Accuracy: %.2f\n',mean(T.Correct))
fprintf('Mean RT: %.3f sec\n',mean(T.RT,'omitnan'))

%% ================== FIGURES ==================

% SIGNALS
figure('Name','Signals','Color','w')
subplot(3,1,1); plot(t, ecg_filt); title('ECG')
subplot(3,1,2); plot(t, resp_filt); title('Respiration')
subplot(3,1,3); plot(t, eda_smooth); title('EDA')

% LEARNING
figure('Name','Learning','Color','w')
smoothed = movmean(double(T.Correct),15);
plot(smoothed,'LineWidth',2); hold on
xline(find(diff(T.Block)~=0),'--k')
ylim([0 1]); title('Learning Curve')

% CARDIAC
figure('Name','Cardiac','Color','w')
bar([mean(T.Correct(T.CardiacFB==1)) ...
     mean(T.Correct(T.CardiacFB==2))])
set(gca,'XTickLabel',{'Systole','Diastole'})
ylim([0 1])

% RESP
figure('Name','Respiration','Color','w')
bar([mean(T.Correct(T.RespFB==1)) ...
     mean(T.Correct(T.RespFB==2))])
set(gca,'XTickLabel',{'Inspiration','Expiration'})
ylim([0 1])

% EDA
figure('Name','EDA','Color','w')
valid = ~isnan(T.EDA_FB);
scatter(T.EDA_FB(valid), double(T.Correct(valid)),'filled')
xlabel('EDA'); ylabel('Correct')

figure('Color','w','Position',[200 200 900 600])

easy_color   = [0 0.45 0.9];
medium_color = [0.85 0.33 0.1];
prob_color   = [0.93 0.69 0.13];

%% ================= LEARNING CURVE =================
subplot(2,1,1)

smoothN = 15;

easy_idx   = T.Block==1;
medium_idx = T.Block==2;
prob_idx   = T.Block==3;

easy   = movmean(double(T.Correct(easy_idx)),smoothN);
medium = movmean(double(T.Correct(medium_idx)),smoothN);
prob   = movmean(double(T.Correct(prob_idx)),smoothN);

hold on

plot(easy,'Color',easy_color,'LineWidth',2)

plot(length(easy)+(1:length(medium)),...
     medium,'Color',medium_color,'LineWidth',2)

plot(length(easy)+length(medium)+(1:length(prob)),...
     prob,'Color',prob_color,'LineWidth',2)

xline(length(easy),'--k')
xline(length(easy)+length(medium),'--k')

ylim([0 1])

ylabel('Accuracy')
xlabel('Trial')
title('Learning Curve by Block')

legend({'Easy','Medium','Prob'},'Location','best')

box off


%% ================= ACCURACY PER BLOCK =================
subplot(2,2,3)

plot(acc_block,'o-','LineWidth',2,'MarkerSize',8,...
     'Color',[0.2 0.2 0.2],'MarkerFaceColor',[0.2 0.2 0.2])

ylim([0 1])

set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'Easy','Medium','Prob'})

ylabel('Accuracy')
title('Accuracy per Block')

box off


%% ================= RT PER BLOCK =================
subplot(2,2,4)

plot(rt_block,'o-','LineWidth',2,'MarkerSize',8,...
     'Color',[0.2 0.2 0.2],'MarkerFaceColor',[0.2 0.2 0.2])

set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'Easy','Medium','Prob'})
ylim([0 1])
ylabel('RT (s)')
title('RT per Block')

box off

exportgraphics(gcf,'Pilot_Behaviour.png','Resolution',300)