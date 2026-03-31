clear; clc; close all;

%% ================= SUBJECT =================
subj = 'P02';

addpath('C:\MyTemp\fieldtrip-20251218');
ft_defaults;

dataset = ['C:\MyTemp\RBLT_project\pilot\RBL_' subj '.edf'];

%% ================= LOAD EVENTS =================
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

%% ================= ALL BLOCKS =================
stim_codes = [101:109 131:139 161:169];
stim_idx = find(ismember(codes,stim_codes));

%% ================= BLOCK LABEL =================
block = zeros(length(stim_idx),1);

for i=1:length(stim_idx)

    c = codes(stim_idx(i));

    if ismember(c,101:109)
        block(i)=1; % easy
    elseif ismember(c,131:139)
        block(i)=2; % medium
    elseif ismember(c,161:169)
        block(i)=3; % probabilistic
    end

end

%% ================= TRIALS =================
prestim  = 0.2;
poststim = 0.8;

trl = zeros(length(stim_idx),3);

for i=1:length(stim_idx)
    s = samples(stim_idx(i));
    trl(i,1)=round(s-prestim*fs);
    trl(i,2)=round(s+poststim*fs);
    trl(i,3)=round(-prestim*fs);
end

%% ================= ACCURACY =================
correct_fb = [30 120 150 180];
all_fb     = [30 31 120 121 150 151 180 181];

accuracy = zeros(length(stim_idx),1);

for i = 1:length(stim_idx)
    for j = stim_idx(i)+1:min(stim_idx(i)+20,length(codes))
        if ismember(codes(j),all_fb)
            if ismember(codes(j),correct_fb)
                accuracy(i)=1;
            end
            break
        end
    end
end

%% ================= RAW ECG RESP =================
cfg=[];
cfg.dataset=dataset;
raw = ft_preprocessing(cfg);

labels = lower(raw.label);

ecg_chan  = find(contains(labels,'ecg'),1);
resp_chan = find(contains(labels,'resp'),1);

ecg  = raw.trial{1}(ecg_chan,:);
resp = raw.trial{1}(resp_chan,:);

%% ================= CARDIAC =================
ecg_filt = bandpass(ecg,[5 20],fs);

[pks,locs] = findpeaks(ecg_filt,...
    'MinPeakHeight',mean(ecg_filt)+1.2*std(ecg_filt),...
    'MinPeakDistance',0.4*fs);

cardiac = zeros(length(stim_idx),1);

for i=1:length(stim_idx)

    s = samples(stim_idx(i));

    prev = find(locs<s,1,'last');
    next = find(locs>s,1,'first');

    if ~isempty(prev)&&~isempty(next)

        phase=(s-locs(prev))/(locs(next)-locs(prev));

        if phase<0.5
            cardiac(i)=1;
        else
            cardiac(i)=2;
        end
    end
end

%% ================= RESP =================
resp_filt = bandpass(resp,[0.1 0.5],fs);
resp_phase = angle(hilbert(resp_filt));

resp_cond = zeros(length(stim_idx),1);

for i=1:length(stim_idx)

    s = samples(stim_idx(i));

    if resp_phase(s)>0
        resp_cond(i)=1;
    else
        resp_cond(i)=2;
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
occ = {'70','71','74','75','76'};
par = {'61','62','67','72','78'};

occ_idx = find(ismember(data.label,occ));
par_idx = find(ismember(data.label,par));

t = data.time{1};

n1_win = t>=0.10 & t<=0.18;
p3_win = t>=0.25 & t<=0.40;

%% ================= ERP FEATURES =================
nTrials = length(data.trial);

N1 = nan(nTrials,1);
P3 = nan(nTrials,1);

for i=1:nTrials

    trial = data.trial{i};

    occ_sig = mean(trial(occ_idx,:),1);
    par_sig = mean(trial(par_idx,:),1);

    N1(i) = mean(occ_sig(n1_win));
    P3(i) = mean(par_sig(p3_win));

end

%% ================= ML DATASET =================
X = [N1 P3 cardiac resp_cond];
Y = accuracy;

valid = ~any(isnan(X),2);

X = X(valid,:);
Y = Y(valid);
block = block(valid);
cardiac = cardiac(valid);
resp_cond = resp_cond(valid);

X = zscore(X);

%% ================= MODEL =================
mdl = fitglm(X,Y,'Distribution','binomial');

Y_prob = predict(mdl,X);
Y_pred = Y_prob>0.5;

acc = mean(Y_pred==Y);

fprintf('\nPrediction accuracy = %.2f %%\n',acc*100)

%% ================= CROSS VALIDATION (manual, toolbox-free) =================

K = 5;
N = length(Y);

rng(1);                     % reproducible
perm = randperm(N);         % shuffle indices
foldSize = floor(N/K);

acc_cv = zeros(K,1);

for k = 1:K
    
    if k < K
        testIdx = perm((k-1)*foldSize+1 : k*foldSize);
    else
        testIdx = perm((k-1)*foldSize+1 : end);
    end
    
    trainIdx = setdiff(1:N, testIdx);
    
    Xtrain = X(trainIdx,:);
    Ytrain = Y(trainIdx);
    
    Xtest  = X(testIdx,:);
    Ytest  = Y(testIdx);
    
    mdl_cv = fitglm(Xtrain,Ytrain,...
        'Distribution','binomial',...
        'Link','logit');
    
    prob = predict(mdl_cv,Xtest);
    pred = prob > 0.5;
    
    acc_cv(k) = mean(pred == Ytest);

end

fprintf('Cross-validated accuracy = %.2f %%\n',mean(acc_cv)*100)
%% ================= OPTIMAL =================
optimal = Y_prob>0.6;

%% ================= PHASE =================
acc_sys = mean(Y(cardiac==1));
acc_dia = mean(Y(cardiac==2));

acc_insp = mean(Y(resp_cond==1));
acc_exp  = mean(Y(resp_cond==2));

fprintf('\nCardiac\nSystole %.2f\nDiastole %.2f\n',acc_sys,acc_dia)
fprintf('\nResp\nInspiration %.2f\nExpiration %.2f\n',acc_insp,acc_exp)

%% ================= BLOCK ANALYSIS =================
names = {'Easy','Medium','Probabilistic'};

for b=1:3

    idx = block==b;

    if sum(idx)==0, continue, end

    acc_insp = mean(Y(idx & resp_cond==1));
    acc_exp  = mean(Y(idx & resp_cond==2));

    fprintf('\n%s block\n',names{b})
    fprintf('Inspiration %.2f\n',acc_insp)
    fprintf('Expiration %.2f\n',acc_exp)

end

%% ================= FIGURE 1 =================
% Model validation

figure('Color','w','Position',[200 200 450 350])

data = [Y_prob(Y==1); Y_prob(Y==0)];
group = [ones(sum(Y==1),1); 2*ones(sum(Y==0),1)];

boxplot(data,group,'Colors','k','Symbol','')

h = findobj(gca,'Tag','Box');

patch(get(h(2),'XData'),get(h(2),'YData'),[0 .4 1],'FaceAlpha',0.35)
patch(get(h(1),'XData'),get(h(1),'YData'),[1 .2 .2],'FaceAlpha',0.35)

set(gca,'XTickLabel',{'Correct feedback','Incorrect feedback'})
ylabel('Predicted learning probability')

title('Model predicts learning success')

box off


%% ================= FIGURE 2 =================
% RESPIRATION × BLOCK

means = nan(3,2);
sems  = nan(3,2);

for b = 1:3
    for r = 1:2
        idx = block==b & resp_cond==r;
        vals = Y(idx);
        means(b,r) = mean(vals);
        sems(b,r)  = std(vals)/sqrt(length(vals));
    end
end

figure('Color','w','Position',[200 200 700 400])

x = 1:3;
offset = 0.18;

hold on

errorbar(x-offset,means(:,1),sems(:,1),...
    'o-','LineWidth',2,'Color',[0 .4 1],'MarkerFaceColor',[0 .4 1])

errorbar(x+offset,means(:,2),sems(:,2),...
    'o-','LineWidth',2,'Color',[1 .2 .2],'MarkerFaceColor',[1 .2 .2])

xlim([0.5 3.5])
ylim([0 1])

set(gca,'XTick',[1 2 3])
set(gca,'XTickLabel',{'Easy','Medium','Probabilistic'})

ylabel('Accuracy')

legend({'Inspiration','Expiration'},'Location','northwest')

title('Respiration phase modulates learning across blocks')

box off


%% ================= FIGURE 3 =================
% CARDIAC × BLOCK

means = nan(3,2);
sems  = nan(3,2);

for b = 1:3
    for c = 1:2
        idx = block==b & cardiac==c;
        vals = Y(idx);
        means(b,c) = mean(vals);
        sems(b,c)  = std(vals)/sqrt(length(vals));
    end
end

figure('Color','w','Position',[200 200 700 400])

x = 1:3;
offset = 0.18;

hold on

errorbar(x-offset,means(:,1),sems(:,1),...
    'o-','LineWidth',2,'Color',[0 .4 1],'MarkerFaceColor',[0 .4 1])

errorbar(x+offset,means(:,2),sems(:,2),...
    'o-','LineWidth',2,'Color',[1 .2 .2],'MarkerFaceColor',[1 .2 .2])

xlim([0.5 3.5])
ylim([0 1])

set(gca,'XTick',[1 2 3])
set(gca,'XTickLabel',{'Easy','Medium','Probabilistic'})

ylabel('Accuracy')

legend({'Systole','Diastole'},'Location','northwest')

title('Cardiac phase modulates learning across blocks')

box off

%% ================= OPTIMAL PHASE FROM MODEL =================

figure('Color','w','Position',[200 200 900 350])

%% Cardiac
subplot(1,2,1)

data = [Y_prob(cardiac==1); Y_prob(cardiac==2)];
group = [ones(sum(cardiac==1),1); 2*ones(sum(cardiac==2),1)];

boxplot(data,group,'Colors','k','Symbol','')

h = findobj(gca,'Tag','Box');

patch(get(h(2),'XData'),get(h(2),'YData'),[0 .4 1],'FaceAlpha',0.35)
patch(get(h(1),'XData'),get(h(1),'YData'),[1 .2 .2],'FaceAlpha',0.35)

set(gca,'XTickLabel',{'Systole','Diastole'})
ylabel('Predicted learning probability')

title('Optimal cardiac phase')

box off


%% Respiration
subplot(1,2,2)

data = [Y_prob(resp_cond==1); Y_prob(resp_cond==2)];
group = [ones(sum(resp_cond==1),1); 2*ones(sum(resp_cond==2),1)];

boxplot(data,group,'Colors','k','Symbol','')

h = findobj(gca,'Tag','Box');

patch(get(h(2),'XData'),get(h(2),'YData'),[0 .4 1],'FaceAlpha',0.35)
patch(get(h(1),'XData'),get(h(1),'YData'),[1 .2 .2],'FaceAlpha',0.35)

set(gca,'XTickLabel',{'Inspiration','Expiration'})
ylabel('Predicted learning probability')

title('Optimal respiration phase')

box off


%% ================= BLOCKWISE OPTIMAL (MODEL) =================

names = {'Easy','Medium','Probabilistic'};

resp_mean = nan(3,2);
card_mean = nan(3,2);

for b = 1:3
    
    for r = 1:2
        idx = block==b & resp_cond==r;
        resp_mean(b,r) = mean(Y_prob(idx));
    end
    
    for c = 1:2
        idx = block==b & cardiac==c;
        card_mean(b,c) = mean(Y_prob(idx));
    end
    
end


figure('Color','w','Position',[200 200 1000 400])

%% ================= RESP =================
subplot(1,2,1)

x = 1:3;
offset = 0.15;

hold on

plot(x-offset,resp_mean(:,1),'o-','LineWidth',2,...
    'Color',[0 .4 1],'MarkerFaceColor',[0 .4 1])

plot(x+offset,resp_mean(:,2),'o-','LineWidth',2,...
    'Color',[1 .2 .2],'MarkerFaceColor',[1 .2 .2])

xlim([0.5 3.5])
ylim([0.5 0.7])

set(gca,'XTick',1:3)
set(gca,'XTickLabel',names)

ylabel('Predicted learning probability')

legend({'Inspiration','Expiration'},'Location','northwest')

title('Block-wise optimal respiration phase')

box off


%% ================= CARDIAC =================
subplot(1,2,2)

hold on

plot(x-offset,card_mean(:,1),'o-','LineWidth',2,...
    'Color',[0 .4 1],'MarkerFaceColor',[0 .4 1])

plot(x+offset,card_mean(:,2),'o-','LineWidth',2,...
    'Color',[1 .2 .2],'MarkerFaceColor',[1 .2 .2])

xlim([0.5 3.5])
ylim([0.5 0.7])

set(gca,'XTick',1:3)
set(gca,'XTickLabel',names)

ylabel('Predicted learning probability')

legend({'Systole','Diastole'},'Location','northwest')

title('Block-wise optimal cardiac phase')

box off


%% ================= ADDITIONAL FIGURES =================
% (secondary results)

%% Feature importance
figure('Color','w','Position',[200 200 450 350])

bar(abs(mdl.Coefficients.Estimate(2:end)))
set(gca,'XTickLabel',{'N1','P3','Cardiac','Resp'})
ylabel('Weight')

title('Feature importance')
box off


%% Optimal vs non optimal
figure('Color','w','Position',[200 200 450 350])

bar([mean(Y(optimal)) mean(Y(~optimal))])
set(gca,'XTickLabel',{'Optimal','Non-optimal'})
ylim([0 1])

ylabel('Accuracy')
title('Optimal physiological state')

box off


%% Cardiac box
figure('Color','w','Position',[200 200 450 350])

data = [Y(cardiac==1); Y(cardiac==2)];
group = [ones(sum(cardiac==1),1); 2*ones(sum(cardiac==2),1)];

boxplot(data,group,'Colors','k','Symbol','')

set(gca,'XTickLabel',{'Systole','Diastole'})
ylabel('Accuracy')

title('Cardiac phase effect')

box off


%% Resp box
figure('Color','w','Position',[200 200 450 350])

data = [Y(resp_cond==1); Y(resp_cond==2)];
group = [ones(sum(resp_cond==1),1); 2*ones(sum(resp_cond==2),1)];

boxplot(data,group,'Colors','k','Symbol','')

set(gca,'XTickLabel',{'Inspiration','Expiration'})
ylabel('Accuracy')

title('Respiration phase effect')

box off

