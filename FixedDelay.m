clc
clear 
close all
%sadaf solgi
%401610045
%% ************************************* Fixed delay Task Analysis ***************************************************
fileName=  "C:\Users\Asus\Desktop\neuroModel\data\SiliconProbeData\FixedDelayTask\HI130_051217_units";
unitNum= 32;

% set parameters for plot PSTH
preCueDur  = 4.2; % in sec: plot how many seconds before go cue onset. 
postCueDur = 1.5; % in sec: plot how many seconds after go cue onset.
    

%% laod file & extract info required for plotting
    load("C:\Users\Asus\Desktop\neuroModel\data\SiliconProbeData\FixedDelayTask\HI130_051217_units");
        
    % extract info
    trialType      = unit(unitNum).Behavior.Trial_types_of_response_vector; % trial type based on outcome
    % 1: correct lick R trial,  2: correct lick L trial,  3: incorrect lick R trial,  4: incorrect lick L trial,
    % 5: correct lick R trial with early lick,   6: correct lick L trial with early lick, 
    % 7: incorrect lick R trial with early lick, 8: incorrect lick L trial with early lick
    % 9: lick R trial but no response,          10: lick L trial but no response,
    % 11: others(unidentified)
    
    photoStimTrial = unit(unitNum).Behavior.stim_trial_vector;  % photo stim trial type
    % 0: no stim,   1: 0.05mW bilateral stim during early dealy,  
    % 2: 0.1mW bilateral stim during early dealy,
    % 3: 0.2mW bilateral stim during early dealy,
    % 4: 0.3mW bilateral stim during early dealy,
    % don't analyze 0.05mW it was too weak to have any behavioral effect

    
    % extract timing info
    T.tAxisPSTH   = -preCueDur-0.1:0.001:postCueDur+0.1;      % T axis for PSTH, shift by 0.1 to remove smoothing artifact at the edge
    T.sampleOnset = [unit(unitNum).Behavior.Sample_start];    % smaple epoch onset
    T.delayOnset  = [unit(unitNum).Behavior.Delay_start];     % delay epoch onset
    T.cueOnset    = [unit(unitNum).Behavior.Cue_start];       % go cue  onset
     
    % calculate mean onset time of each epoch for plotting   
    % exclude early lick and non-response trial by "trialType<5"
    T.meanSampleOnset = mean(T.sampleOnset(trialType<5));
    T.meanDelayOnset  = mean(T.delayOnset(trialType<5));
    T.meanCueOnset    = mean(T.cueOnset(trialType<5));
    
    
    
    %% extract spikes  from each trial
    spikeTimes       = unit(unitNum).SpikeTimes; % time bin of spike
    trialIdxOfSpikes = unit(unitNum).Trial_idx_of_spike; % time bin of spike
    
    trialRange            = unit(unitNum).Trial_info.Trial_range_to_analyze; % range of trials to analyze
    trialTypeInRange      = trialType(trialRange(1) : trialRange(2));        % trial type in range
    photoStimTrialInRange = photoStimTrial(trialRange(1) : trialRange(2));   % stim trial type in range
    
    PSTH        = nan( trialRange(2)- trialRange(1)+1,numel(T.tAxisPSTH));
    
    for tr = trialRange(1) : trialRange(2)

        % extarct spikes of each trial
        SpikesTmp   = spikeTimes(trialIdxOfSpikes == tr) - T.cueOnset(tr); % align to go cue 
        PSTHTmp     = hist(SpikesTmp,T.tAxisPSTH)*1000;
        PSTH(tr,:)  = smooth(PSTHTmp,101);
        
    end
    
    %% plot raster & PSTH & tunning curve for Fixed delay task
    %% correct lick R trilas
    trialsToPlot    = find(trialType == 1 & photoStimTrial == 0); 
    spikeTimesOfTrials = []; trialIndexOfTrials = [];
    
    for i = 1:numel(trialsToPlot)
        
        tr = trialsToPlot(i); % trial number
        spikeCountRateinTrials(tr)=0;
        
        spikeInTrial       = spikeTimes(trialIdxOfSpikes==tr)- T.meanCueOnset; % spikes in this trial
        
        spikeTimesOfTrials = [spikeTimesOfTrials;spikeInTrial]; 
        trialIndexOfTrials = [trialIndexOfTrials;ones(size(spikeInTrial))*i];
        
        spikeCountinTrials(tr)= numel(spikeTimesOfTrials);
        spikeCountRateinTrials(tr)=spikeCountRateinTrials(tr)+ spikeCountinTrials(tr)/2.6;
            
    end
    AvgFrR= spikeCountRateinTrials(tr)/numel(trialsToPlot);
    figure;
    subplot(3,1,1);
    plot(spikeTimesOfTrials, trialIndexOfTrials,'b.')
    ylim([0.5 max(trialIndexOfTrials)+0.5])
    xlim([min(T.tAxisPSTH)+0.1 max(T.tAxisPSTH)-0.1])
    ylabel CorrectR
    
    
    %% correct lick L trilas
    trialsToPlot    = find(trialType == 2 & photoStimTrial == 0); 
    spikeTimesOfTrials = []; trialIndexOfTrials = [];
    
    for i = 1:numel(trialsToPlot)
        
        tr = trialsToPlot(i); % trial number
        spikeCountRateinTrials(tr)=0;
        
        spikeInTrial       = spikeTimes(trialIdxOfSpikes==tr)- T.meanCueOnset; % spikes in this trial
        
        spikeTimesOfTrials = [spikeTimesOfTrials;spikeInTrial]; 
        trialIndexOfTrials = [trialIndexOfTrials;ones(size(spikeInTrial))*i];
        
        spikeCountinTrials(tr)= numel(spikeTimesOfTrials);
        spikeCountRateinTrials(tr)=spikeCountRateinTrials(tr)+ spikeCountinTrials(tr)/2.6;
            
    end
    AvgFrL= spikeCountRateinTrials(tr)/numel(trialsToPlot);
    subplot(3,1,2);
    plot(spikeTimesOfTrials, trialIndexOfTrials,'r.')
    ylim([0.5 max(trialIndexOfTrials)+0.5])
    xlim([min(T.tAxisPSTH)+0.1 max(T.tAxisPSTH)-0.1])
    ylabel CorrectL
    
    % PSTH of correct lick R and L trials   
    subplot(3,1,3);hold on
    plot(T.tAxisPSTH,mean(PSTH(trialTypeInRange==1 & photoStimTrialInRange==0,:)),'b')
    plot(T.tAxisPSTH,mean(PSTH(trialTypeInRange==2 & photoStimTrialInRange==0,:)),'r')
    formatFigs(T)
    ylabel({'Correct non-stim trials';'Spikes per s'})
    
    %% incorrect lick R trial
    trialsToPlot    = find(trialType == 3 & photoStimTrial == 0); 
    spikeTimesOfTrials = []; trialIndexOfTrials = [];
    
    for i = 1:numel(trialsToPlot)
        
        tr = trialsToPlot(i); % trial number
        spikeCountRateinTrials(tr)=0;
        
        spikeInTrial       = spikeTimes(trialIdxOfSpikes==tr)- T.meanCueOnset; % spikes in this trial
        
        spikeTimesOfTrials = [spikeTimesOfTrials;spikeInTrial]; 
        trialIndexOfTrials = [trialIndexOfTrials;ones(size(spikeInTrial))*i];
        
        spikeCountinTrials(tr)= numel(spikeTimesOfTrials);
        spikeCountRateinTrials(tr)=spikeCountRateinTrials(tr)+ spikeCountinTrials(tr)/2.6;
            
    end
    AvgFrNR= spikeCountRateinTrials(tr)/numel(trialsToPlot);
    figure;
    subplot(3,1,1);
    plot(spikeTimesOfTrials, trialIndexOfTrials,'g.')
    ylim([0.5 max(trialIndexOfTrials)+0.5])
    xlim([min(T.tAxisPSTH)+0.1 max(T.tAxisPSTH)-0.1])
    ylabel inCorrectR
  
    
    %% incorrect lick L trial
    trialsToPlot    = find(trialType == 4 & photoStimTrial == 0); 
    spikeTimesOfTrials = []; trialIndexOfTrials = [];
    
    for i = 1:numel(trialsToPlot)
        
        tr = trialsToPlot(i); % trial number
        spikeCountRateinTrials(tr)=0;
        
        spikeInTrial       = spikeTimes(trialIdxOfSpikes==tr)- T.meanCueOnset; % spikes in this trial
        
        spikeTimesOfTrials = [spikeTimesOfTrials;spikeInTrial]; 
        trialIndexOfTrials = [trialIndexOfTrials;ones(size(spikeInTrial))*i];
        
        spikeCountinTrials(tr)= numel(spikeTimesOfTrials);
        spikeCountRateinTrials(tr)=spikeCountRateinTrials(tr)+ spikeCountinTrials(tr)/2.6;
            
    end
    AvgFrNL= spikeCountRateinTrials(tr)/numel(trialsToPlot);
    subplot(3,1,2)
    plot(spikeTimesOfTrials, trialIndexOfTrials,'y.')
    ylim([0.5 max(trialIndexOfTrials)+0.5])
    xlim([min(T.tAxisPSTH)+0.1 max(T.tAxisPSTH)-0.1])
    ylabel inCorrectL
    
    % PSTH of incorrect lick R and L trials   
    subplot(3,1,3);hold on
    plot(T.tAxisPSTH,mean(PSTH(trialTypeInRange==3 & photoStimTrialInRange==0,:)),'g')
    plot(T.tAxisPSTH,mean(PSTH(trialTypeInRange==4 & photoStimTrialInRange==0,:)),'y')
    formatFigs(T)
    ylabel({'inCorrect non-stim trials';'Spikes per s'})
    %% Plot Tunning Curve
    AvgFR= horzcat(AvgFrR,AvgFrL,AvgFrNR,AvgFrNL);
    stim= 1:4; 
    [fitresult, gof] = GaussFit(stim, AvgFR)
    
function [] = formatFigs(T)
    % format the figure
    % add sample, delay and go cue onset
    
    xlim([min(T.tAxisPSTH)+0.1 max(T.tAxisPSTH)-0.1]);
    yRange = ylim();
    plot([0 0],yRange,'k:')
    plot([T.meanSampleOnset - T.meanCueOnset  T.meanSampleOnset - T.meanCueOnset]   ,yRange,'k:')
    plot([T.meanDelayOnset - T.meanCueOnset    T.meanDelayOnset - T.meanCueOnset]   ,yRange,'k:')
    xlabel('Time from go cue osnet (s)')

end