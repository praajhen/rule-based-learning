clear; clc; close all;

addpath('C:\MyTemp\fieldtrip-20251218');
ft_defaults;

%% READ EVENTS
event = ft_read_event('C:\MyTemp\RBLT_project\pilot\RBL_P02.edf');

% Remove rows with empty values
emptyValueIndices = find(arrayfun(@(x) isempty(x.value), event));
event(emptyValueIndices) = [];

%% EXTRACT TRIGGER CODES AND TIMES
codes = zeros(length(event),1);
times = zeros(length(event),1);

for i = 1:length(event)
    
    val = event(i).value;
    num = regexp(val,'\d+','match');
    
    codes(i) = str2double(num{end});
    times(i) = event(i).sample;
    
end

%% DEFINE TRIGGER GROUPS

stim_codes = [11:19 101:109 131:139 161:169];
resp_codes = [21:23 111:113 141:143 171:173];
raw_resp   = [1 2 3];
fb_codes   = [30 31 120 121 150 151 180 181];

%% RECONSTRUCT TRIALS

trial = [];
block = [];
stimulus = [];
response = [];
feedback = [];
RT = [];

t = 1;

for i = 1:length(codes)-3
    
    if ismember(codes(i),stim_codes)
        
        stim = codes(i);
        stim_time = times(i);
        
        resp = NaN;
        fb   = NaN;
        rt   = NaN;
        
        % find raw button for RT
        for j = i+1:min(i+5,length(codes))
            
            if ismember(codes(j),raw_resp)
                
                rt = times(j) - stim_time;
                
            end
            
        end
        
        % find response trigger
        for j = i+1:min(i+6,length(codes))
            
            if ismember(codes(j),resp_codes)
                
                resp = codes(j);
                break
                
            end
            
        end
        
        % find feedback
        for j = i+1:min(i+7,length(codes))
            
            if ismember(codes(j),fb_codes)
                
                fb = codes(j);
                break
                
            end
            
        end
        
        trial(t) = t;
        stimulus(t) = stim;
        response(t) = resp;
        feedback(t) = fb;
        RT(t) = rt;
        
        % identify block
        if stim>=11 && stim<=19
            block(t)=0;
        elseif stim>=101 && stim<=109
            block(t)=1;
        elseif stim>=131 && stim<=139
            block(t)=2;
        elseif stim>=161 && stim<=169
            block(t)=3;
        end
        
        t = t+1;
        
    end
end

%% CREATE TRIAL TABLE

T = table(trial',block',stimulus',response',feedback',RT',...
    'VariableNames',{'Trial','Block','Stimulus','Response','Feedback','RT_samples'});

%% CONVERT RT TO SECONDS

fs = 1000; % change if EDF sampling rate is different
T.RT = T.RT_samples / fs;

%% COMPUTE ACCURACY

correct = ismember(T.Feedback,[30 120 150 180]);
T.Correct = correct;

disp(T)

%% TRIGGER TIMELINE

figure
stem(times,codes,'filled')
xlabel('Sample')
ylabel('Trigger')
title('Trigger Timeline')

%% TRIAL STRUCTURE

figure
plot(T.Stimulus,'o-','LineWidth',1.5)
hold on
plot(T.Response,'o-','LineWidth',1.5)
plot(T.Feedback,'o-','LineWidth',1.5)
legend('Stimulus','Response','Feedback')
xlabel('Trial')
ylabel('Trigger Code')
title('Trial Structure Check')

%% RT DISTRIBUTION

figure
histogram(T.RT,30)
xlabel('Reaction Time (s)')
ylabel('Count')
title('RT Distribution')

%% ACCURACY PER BLOCK

blocks = unique(T.Block);

fprintf('\nAccuracy per block:\n')

for b = blocks'
    
    idx = T.Block==b;
    acc = mean(T.Correct(idx));
    
    fprintf('Block %d Accuracy = %.2f\n',b,acc)
    
end

%% MEAN RT PER BLOCK

fprintf('\nMean RT per block:\n')

for b = blocks'
    
    idx = T.Block==b;
    mrt = mean(T.RT(idx),'omitnan');
    
    fprintf('Block %d RT = %.3f sec\n',b,mrt)
    
end

%% STIMULUS DISTRIBUTION

figure
histogram(T.Stimulus)
xlabel('Stimulus Code')
ylabel('Count')
title('Stimulus Distribution')

%% TRIGGER INTERVAL CHECK

dt = diff(times);

figure
plot(dt)
xlabel('Trigger Index')
ylabel('Sample Interval')
title('Trigger Timing Check')

%% DETECT BAD TRIALS

bad_trials = isnan(T.Response) | isnan(T.Feedback);

fprintf('\nNumber of incomplete trials: %d\n',sum(bad_trials))

if sum(bad_trials)>0
    
    disp('Problematic trials:')
    disp(T(bad_trials,:))
    
end

figure
plot(T.Trial,cumsum(T.Correct)./(1:length(T.Correct))','LineWidth',2)
ylim([0 1])
xlabel('Trial')
ylabel('Cumulative accuracy')
title('Learning curve')