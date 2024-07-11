clc
clear 
close all
%sadaf solgi
%401610045
%% ************************************* Random delay Task Analysis ***************************************************

    fileName = 'RandomDelayTask\withoutPerturbation/HI134_062017_units.mat';
    unitNum  = 25;

 % set parameters for plot
    preCueDur  = 2; % in sec: plot how many seconds before delay  onset. 
    postCueDur = 4; % in sec: plot how many seconds after delay onset.

    
    
    %% laod file & extract info required for plotting
    load(fileName)
        
    % extract info
    trialType      = unit(unitNum).Behavior.Trial_types_of_response_vector; % trial type based on outcome
    % 1: correct lick R trial,  2: correct lick L trial,  3: incorrect lick R trial,  4: incorrect lick L trial,
    % 5: correct lick R trial with early lick,   6: correct lick L trial with early lick, 
    % 7: incorrect lick R trial with early lick, 8: incorrect lick L trial with early lick
    % 9: lick R trial but no response,          10: lick L trial but no response,
    % 11: others(unidentified)
    
    delayDur       = unit(unitNum).Behavior.delay_dur_id; % delay duration type of each trial
       
    photoStimTrial = unit(unitNum).Behavior.stim_trial_vector;  % photo stim trial type
    % 0: no stim,   1: 0.05mW bilateral stim during early dealy,  
    % 2: 0.1mW bilateral stim during early dealy,
    % 3: 0.2mW bilateral stim during early dealy,
    % 4: 0.3mW bilateral stim during early dealy,
    % don't analyze 0.05mW it was too weak to have any behavioral effect

    
    % extract timing info
    T.tAxisPSTH   = -preCueDur-0.1:0.001:postCueDur+0.1;      % T axis for PSTH, shift by 0.1 to remove smoothing artifact at the edge
    T.sampleOnset = [unit(unitNum).Behavior.Sample_start];  % smaple epoch onset
    T.delayOnset  = [unit(unitNum).Behavior.Delay_start];     % delay epoch onset
    T.cueOnset    = [unit(unitNum).Behavior.Cue_start];       % go cue  onset
     
    
    % calculate mean onset time of each epoch for plotting   
    % exclude early lick and non-response trial by "trialType<5"
    T.meanSampleOnset = mean(T.sampleOnset(trialType<5));
    T.meanDelayOnset  = mean(T.delayOnset(trialType<5));
    % calculate cue onset for each delay duration type
    T.delayTypes      = unique(delayDur);
    for p = 1:numel(T.delayTypes)
        T.meanCueOnset(p)    = mean(T.cueOnset(trialType<5 & delayDur==T.delayTypes(p)));
    end
    
    
    
    %% extract spikes  from each trial
    spikeTimes       = unit(unitNum).SpikeTimes; % time bin of spike
    trialIdxOfSpikes = unit(unitNum).Trial_idx_of_spike; % time bin of spike
    
    
    trialRange            = unit(unitNum).Trial_info.Trial_range_to_analyze; % range of trials to analyze
    trialTypeInRange      = trialType(trialRange(1) : trialRange(2));        % trial type in range
    delayDurInRange       = delayDur(trialRange(1) : trialRange(2));         % delay duration of trial in range
    photoStimTrialInRange = photoStimTrial(trialRange(1) : trialRange(2));   % stim trial type in range

    
    PSTH        = nan( trialRange(2)- trialRange(1)+1,numel(T.tAxisPSTH));
    PSTHpreGo   = nan( trialRange(2)- trialRange(1)+1,numel(T.tAxisPSTH));
    
    for tr = trialRange(1) : trialRange(2)

        % extarct spikes of each trial
        SpikesTmp   = spikeTimes(trialIdxOfSpikes == tr) - T.delayOnset(tr); % align to go cue 
        PSTHTmp     = hist(SpikesTmp,T.tAxisPSTH)*1000;
        PSTHSmooth  = smooth(PSTHTmp,101);
        PSTH(tr,:)  = PSTHSmooth;
        
        % extract spikes before go cue
        timePostGoCue = T.tAxisPSTH > (T.cueOnset(tr)-T.delayOnset(tr)) - 0.05; % 50ms before go cue onset (to remove smoothing artifact)
        PSTHSmooth(timePostGoCue) = nan;
        PSTHpreGo(tr,:)  = PSTHSmooth;
    end
    
    
       % plot raster & PSTH of correct non-stim trial (EDF 8k)
    % for the paper, we subsampled trials to plot raster, which we are not
    % doing here
     
    %% correct lick R trial without perturbation    
    subplot(3,1,1);hold on
    trialsToPlot    = find(trialType == 1 & photoStimTrial == 0); % correct lick R trilas
    delayDurToPlot  = delayDurInRange(trialsToPlot);              % delay duration of those trials
    spikeTimesOfTrials = []; trialIndexOfTrials = [];
   
    delayTypes            = unique(delayDurToPlot);
    numOfTrialOfEachDelay = zeros(numel(delayTypes),1);
    spikeCountRateinTrials(tr)=0;
    for p = 1:numel(delayTypes)
        
        trialsWithThisDelay = find(delayDurToPlot == delayTypes(p)); % trial with this delay duration
        
        for q = 1:numel(trialsWithThisDelay)

            tr = trialsToPlot(trialsWithThisDelay(q)); % trial number

            spikeInTrial       = spikeTimes(trialIdxOfSpikes==tr)- T.meanDelayOnset; % spikes in this trial

            spikeTimesOfTrials = [spikeTimesOfTrials;spikeInTrial]; 
            trialIndexOfTrials = [trialIndexOfTrials;ones(size(spikeInTrial))*(sum(numOfTrialOfEachDelay)+q)];
            
           
        end
        
        spikeCountinTrials(tr)= numel(spikeTimesOfTrials);
        spikeCountRateinTrials(tr)=spikeCountRateinTrials(tr)+ spikeCountinTrials(tr)/2.8;

        numOfTrialOfEachDelay(p) = numel(trialsWithThisDelay);
    end

    AvgFrR= spikeCountRateinTrials(tr)/numel(trialsToPlot);
        
    plot(spikeTimesOfTrials, trialIndexOfTrials,'b.')
    formatFigs(T)
    ylim([0.5 max(trialIndexOfTrials)+0.5])
    ylabel CorrectR
    
    %% correct lick L trial without perturbation    
    subplot(3,1,2);hold on
    trialsToPlot    = find(trialType == 2 & photoStimTrial == 0); % correct lick R trilas
    delayDurToPlot  = delayDurInRange(trialsToPlot);              % delay duration of those trials
    spikeTimesOfTrials = []; trialIndexOfTrials = [];
   
    delayTypes            = unique(delayDurToPlot);
    numOfTrialOfEachDelay = zeros(numel(delayTypes),1);
    
    for p = 1:numel(delayTypes)
        
        trialsWithThisDelay = find(delayDurToPlot == delayTypes(p)); % trial with this delay duration
        
        for q = 1:numel(trialsWithThisDelay)

            tr = trialsToPlot(trialsWithThisDelay(q)); % trial number

            spikeInTrial       = spikeTimes(trialIdxOfSpikes==tr)- T.meanDelayOnset; % spikes in this trial

            spikeTimesOfTrials = [spikeTimesOfTrials;spikeInTrial]; 
            trialIndexOfTrials = [trialIndexOfTrials;ones(size(spikeInTrial))*(sum(numOfTrialOfEachDelay)+q)];

        end
        spikeCountinTrials(tr)= numel(spikeTimesOfTrials);
        spikeCountRateinTrials(tr)=spikeCountRateinTrials(tr)+ spikeCountinTrials(tr)/2.8;
        
        numOfTrialOfEachDelay(p) = numel(trialsWithThisDelay);
    end
    AvgFrL= spikeCountRateinTrials(tr)/numel(trialsToPlot);
    plot(spikeTimesOfTrials, trialIndexOfTrials,'r.')
    formatFigs(T)
    ylim([0.5 max(trialIndexOfTrials)+0.5])
    ylabel CorrectL
    
    %% PSTH of correct lick R and L trials   
    subplot(3,1,3);hold on
    plot(T.tAxisPSTH,nanmean(PSTHpreGo(trialTypeInRange==1 & photoStimTrialInRange==0,:)),'b')
    plot(T.tAxisPSTH,nanmean(PSTHpreGo(trialTypeInRange==2 & photoStimTrialInRange==0,:)),'r')
    formatFigs(T)
    ylabel({'Correct non-stim trials';'Spikes per s'})
    hold off 
    
    
    %% incorrect lick R trial without perturbation 
    figure;
    subplot(3,1,1);hold on
    trialsToPlot    = find(trialType == 3 & photoStimTrial == 0); % correct lick R trilas
    delayDurToPlot  = delayDurInRange(trialsToPlot);              % delay duration of those trials
    spikeTimesOfTrials = []; trialIndexOfTrials = [];
   
    delayTypes            = unique(delayDurToPlot);
    numOfTrialOfEachDelay = zeros(numel(delayTypes),1);
    
    for p = 1:numel(delayTypes)
        
        trialsWithThisDelay = find(delayDurToPlot == delayTypes(p)); % trial with this delay duration
        
        for q = 1:numel(trialsWithThisDelay)

            tr = trialsToPlot(trialsWithThisDelay(q)); % trial number

            spikeInTrial       = spikeTimes(trialIdxOfSpikes==tr)- T.meanDelayOnset; % spikes in this trial

            spikeTimesOfTrials = [spikeTimesOfTrials;spikeInTrial]; 
            trialIndexOfTrials = [trialIndexOfTrials;ones(size(spikeInTrial))*(sum(numOfTrialOfEachDelay)+q)];

        end
        
        spikeCountinTrials(tr)= numel(spikeTimesOfTrials);
        spikeCountRateinTrials(tr)=spikeCountRateinTrials(tr)+ spikeCountinTrials(tr)/2.8;
        
        numOfTrialOfEachDelay(p) = numel(trialsWithThisDelay);
    end
    AvgFrNR= spikeCountRateinTrials(tr)/numel(trialsToPlot);
    plot(spikeTimesOfTrials, trialIndexOfTrials,'g.')
    formatFigs(T)
    ylim([0.5 max(trialIndexOfTrials)+0.5])
    ylabel inCorrectR
    
    %% incorrect lick L trial without perturbation    
    subplot(3,1,2);hold on
    trialsToPlot    = find(trialType == 4 & photoStimTrial == 0); % correct lick R trilas
    delayDurToPlot  = delayDurInRange(trialsToPlot);              % delay duration of those trials
    spikeTimesOfTrials = []; trialIndexOfTrials = [];
   
    delayTypes            = unique(delayDurToPlot);
    numOfTrialOfEachDelay = zeros(numel(delayTypes),1);
    
    for p = 1:numel(delayTypes)
        
        trialsWithThisDelay = find(delayDurToPlot == delayTypes(p)); % trial with this delay duration
        
        for q = 1:numel(trialsWithThisDelay)

            tr = trialsToPlot(trialsWithThisDelay(q)); % trial number

            spikeInTrial       = spikeTimes(trialIdxOfSpikes==tr)- T.meanDelayOnset; % spikes in this trial

            spikeTimesOfTrials = [spikeTimesOfTrials;spikeInTrial]; 
            trialIndexOfTrials = [trialIndexOfTrials;ones(size(spikeInTrial))*(sum(numOfTrialOfEachDelay)+q)];

        end
        spikeCountinTrials(tr)= numel(spikeTimesOfTrials);
        spikeCountRateinTrials(tr)=spikeCountRateinTrials(tr)+ spikeCountinTrials(tr)/2.8;
        
        numOfTrialOfEachDelay(p) = numel(trialsWithThisDelay);
    end
    AvgFrNL= spikeCountRateinTrials(tr)/numel(trialsToPlot);
    plot(spikeTimesOfTrials, trialIndexOfTrials,'y.')
    formatFigs(T)
    ylim([0.5 max(trialIndexOfTrials)+0.5])
    ylabel inCorrectL
    
    %% PSTH of correct lick R and L trials   
    subplot(3,1,3);hold on
    plot(T.tAxisPSTH,nanmean(PSTHpreGo(trialTypeInRange==3 & photoStimTrialInRange==0,:)),'g')
    plot(T.tAxisPSTH,nanmean(PSTHpreGo(trialTypeInRange==4 & photoStimTrialInRange==0,:)),'y')
    formatFigs(T)
    ylabel({'inCorrect non-stim trials';'Spikes per s'})

    AvgFr= horzcat(AvgFrR,AvgFrL,AvgFrNR,AvgFrNL);
    stim=1:4;
    [fitresult, gof] = FourierFit(stim, AvgFr)
     
     
function [] = formatFigs(T)
    % format the figure
    % add sample, delay and go cue onset
    
    xlim([min(T.tAxisPSTH)+0.1 max(T.tAxisPSTH)-0.1]);
    yRange = ylim();
    plot([0 0],yRange,'k:')
    plot([T.meanSampleOnset - T.meanDelayOnset  T.meanSampleOnset - T.meanDelayOnset] ,yRange,'k:')
    xlabel('Time from delay epoch osnet (s)')

end