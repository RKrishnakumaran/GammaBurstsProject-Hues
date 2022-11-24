cd('./Data_Analysis/');
basePath = '..\GammaHarmonics Segmented Data\'; %pwd;

%%
close all
%%
lenAnalysisProcessing('alpa',basePath, true);
%%
lenAnalysisProcessing('tutu',basePath, true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function lenAnalysisProcessing(subjectName, basePath, selgoodBursts)
% Output folder
saveFolder = fullfile('..', 'savedDataMP_fullLenMP', 'processedData', subjectName);
mkdir(saveFolder);
tag = [saveFolder,'/'];
postfix = '';
if selgoodBursts
    postfix = '_goodBurstsOnly';
end
outputfolder = fullfile(tag,['LengthAnalysis',postfix]); mkdir(outputfolder);
% outputfolder = fullfile(tag,['LengthAnalysis',postfix]); mkdir(outputfolder);
% HighRMSelectrodes
highRMSElectrodesStruct = load(fullfile(basePath, 'Color','analyzeElectrodes',subjectName,'highRMSElectrodes'));        
highRMSElectrodes = highRMSElectrodesStruct.highRMSElectrodes; % 64-M1 and 16-M2 electrodes are common across different hues and achromatic grating
  
% setting up subject-specific variables and figure handles
if strcmp(subjectName,'alpa')
    subjectNum = 1;
%     lenHuesfigure = alpahuefigure;
    % other variables and figures
elseif strcmp(subjectName,'tutu')
    subjectNum = 2;
%     lenHuesfigure = tutuhuefigure;
    % other variables and figures
else
    error('Unknown subject!');
end

% setting up plot parameters and figures
colors = [hsv(36); 0.6*[1,1,1]];

global lenHuesfigure;
global lenElesfigure;

try 
    figure(lenHuesfigure);
catch
    lenHuesfigure = figure('windowstate','maximized','InvertHardcopy',false);
end

% data accumulated in hue loop
medianBL_hues_eles = [];
gammapower_hues_eles = [];
badBurstLens = [];
badBurstDurs = [];
badBurstFreqs = [];
goodBurstFreqs = [];

medianBLs_HUE = {[],[]};
for ihue = 1:37
    %% setting up
    hue = (ihue-1)*10;
    Hue = num2str(hue);
    disp(hue);
    
    if ihue == 37 && strcmp(subjectName, 'alpa')
        expType = 'SForiAchro'; stimType = ''; achroFlag = 1;
    else
        expType = 'Color'; stimType = Hue; achroFlag = 0;
    end
    
    %% loading files
%     if isempty(stimType)
%         stimno = 37;
%     else
%         stimno = str2num(stimType)/10+1;
%     end
    
    % Get Corresponding experiment/data
    subjectID = [subjectName expType];
    if strcmp(subjectID,'alpaColor')
        subjectName = 'alpa';expDate = '301215'; protocolName = 'GRF_001'; % 488: Hue fullscreen
    elseif strcmp(subjectID,'alpaSForiAchro')
        subjectName = 'alpa'; expDate = '301215'; protocolName = 'GRF_005'; % SFOri - alpa Vinay
    elseif strcmp(subjectName,'tutu')
        subjectName = 'tutu'; expDate = '191016'; protocolName = 'GRF_001'; % 111: Hue fullscreen
    end
    
    gridType = 'Microelectrode';
    folderSourceString = fullfile(basePath, expType);
    folderBase = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
    folderLFP = fullfile(folderBase,'segmentedData','LFP');
    % TimeVals, FS and StimPos
    load(fullfile(folderLFP,'lfpInfo'),'timeVals');
    Fs = 1./(timeVals(2)-timeVals(1));
    
    badTrialFile = fullfile(folderBase,'segmentedData','badTrials');
    load(badTrialFile,'badTrials');
    paramCombinationsFile = fullfile(folderBase,'extractedData','parameterCombinations');
    Params = load(paramCombinationsFile);
    
    if strcmp(subjectID,'alpaSForiAchro')
        highRMSElectrodesStruct = load(fullfile(folderSourceString,'analyzeElectrodes', 'alpaLFPElectrodeList.mat'));
        t  = load(fullfile(folderLFP,'lfpInfo.mat'));
        timeVals = t.timeVals;
        
        highRMSElectrodes_hue = highRMSElectrodesStruct.alpaLFPElectrodeList;
        highRMSElectrodes_hue = setdiff(highRMSElectrodes_hue, 4); % 4 is an extra electrode compared to color cases (65 v 64)
        assert(isempty(setdiff(highRMSElectrodes, highRMSElectrodes_hue)),'Check! highRMSElectrodes Mismatch across hues');
        
        cVals = Params.cValsUnique; oVals = Params.oValsUnique; fVals = Params.fValsUnique;
        %         c = find(cVals ==100); % contrast 100
        o = (oVals ==90); %90 degree
        %         f = find(round(fVals) == 2); % SF = 2
        goodPos = Params.parameterCombinations{1,1,1,1,o};
        goodPos = setdiff(goodPos, badTrials);
    else
        highRMSElectrodesStruct = load(fullfile(folderSourceString,'analyzeElectrodes',subjectName,'highRMSElectrodes'));
        highRMSElectrodes_hue = highRMSElectrodesStruct.highRMSElectrodes;
        assert(isempty(setdiff(highRMSElectrodes, highRMSElectrodes_hue)),'Check! highRMSElectrodes Mismatch across hues');
        
        cVals = Params.cValsUnique; oVals = Params.oValsUnique;
        goodPosAll = cell(1,length(oVals));
        c = find(cVals ==100); % contrast 100
        for o = 1:length(oVals)
            goodPosAll{o} = Params.parameterCombinations{1,1,1,1,o,c,1};
            goodPosAll{o} = setdiff(goodPosAll{o},badTrials);
        end
        stimIndex = str2double(stimType)/10 + 1;
        goodPos = goodPosAll{stimIndex};
        % stimIndex - represents Color number, ranging 1 to 36
        % representing hues 0 to 350 with interval size 10
    end
    
    %% MP length analysis
    %             %MP parameters
    %             thresholdFactor=0.5;
    %             maxIteration=50;
    %             adaptiveDictionaryParam=0.9;
    %             dictionarySize=2500000;
    %             % Renaming parameters
    %             gammaFreqRangeHz = gammaRangeHz;
    stimulusPeriodS = [0.1 0.8]; %stimPeriod;
    baselinePeriodS = [-0.5 0]; %baselinePeriod;
    
%     tag = [saveFolder,'/'];
    
    % Output folder
%     outputfolder = fullfile(tag,'LengthAnalysis'); mkdir(outputfolder);
%     outputFolder_hue = fullfile(outputfolder,'Hue-wise',['hue' num2str(ihue)]); mkdir(outputFolder_hue);
    
    %% -
    if ihue == 37
        load(fullfile(saveFolder, [subjectName 'Achro' '.mat']),'electrodeIDs', 'allBurstsAllElectrodes','tBurstsAllElectrodes','fBurstsAllElectrodes','burstTrialAllElectrodes','goodBurstAllElectrodes','numBursts','numBurstsAllElectrodes','medianBurstLength','phaseDiff','gammaFreq', ...
            'diffPower','psdST', 'psdBL', 'goodPos', 'stPos', 'highRMSElectrodes', 'timeVals', 'freqVals', ...
            'subjectName', 'expType', 'stimType', 'expDate', 'protocolName', 'mt', 'Fs', 'gammaRangeHz')
    else
        load(fullfile(saveFolder, [subjectName expType stimType '.mat']), 'electrodeIDs', 'allBurstsAllElectrodes','tBurstsAllElectrodes','fBurstsAllElectrodes','burstTrialAllElectrodes','goodBurstAllElectrodes','numBursts','numBurstsAllElectrodes','medianBurstLength','phaseDiff','gammaFreq', ...
            'diffPower','psdST', 'psdBL', 'goodPos', 'stPos', 'highRMSElectrodes', 'timeVals', 'freqVals', ...
            'subjectName', 'expType', 'stimType', 'expDate', 'protocolName', 'mt', 'Fs', 'gammaRangeHz')
    end
    
    
    %% (Hue,electrode)-wise plots and analysis (plotting done in electrode loop)
    % Accumulate values (to be used in for electrodes loop)
    % scatter-plot burst lengths of different colors in each electrode (1 plot per electrode)
    % burst length
    % TF of avg across trials
    %% Hue-wise plots and analysis
    gammapower_hue = mean(diffPower);
    meanBL_hue_eles = []; semedianBL_hue_eles = [];
    badBurstLen = []; badBurstDur = []; badBurstFreq = []; goodBurstFreq = [];
    for iele = 1:numel(highRMSElectrodes)
        %% only good bursts? or all bursts?
        if selgoodBursts
            
            ind = (tBurstsAllElectrodes{iele}(:)>=0.2 & tBurstsAllElectrodes{iele}(:)<=0.6); %goodBurstAllElectrodes{iele}==1;
            ind = ind & (fBurstsAllElectrodes{iele}(:)>=20 & fBurstsAllElectrodes{iele}(:)<=80);
            ind = ind & (tBurstsAllElectrodes{iele}(:)-allBurstsAllElectrodes{iele}(:)/2>=0 ...
                ...%& tBurstsAllElectrodes{iele}(:)+allBurstsAllElectrodes{iele}(:)/2<=0.8...
                );
%             ind = tBurstsAllElectrodes{iele}(:)>=0.2 & tBurstsAllElectrodes{iele}(:)<=0.6; %goodBurstAllElectrodes{iele}==1;
        else
            ind = ones(size(goodBurstAllElectrodes{iele}))==1;
        end
        burstlens = allBurstsAllElectrodes{iele}(ind) * 1e3; % ms
        bursttrials = burstTrialAllElectrodes{iele}(ind);
        setbursttrials = union([], bursttrials);
        medianBL_trial  = [];
        for bt = setbursttrials(:)' 
            medianBL_trial = [medianBL_trial, median(burstlens(bursttrials==bt))];
        end
        meanBL_hue_eles = [meanBL_hue_eles, mean(medianBL_trial)];
%         semedianBL_hue_eles = [semedianBL_hue_ele, getSEMedian(medianBL_trial)];
        goodburstfreq_vals = fBurstsAllElectrodes{iele}(ind);
        goodBurstFreq = [goodBurstFreq, goodburstfreq_vals(:)' ];
        badburstfreq_vals = fBurstsAllElectrodes{iele}(~ind);
        badBurstFreq = [badBurstFreq, badburstfreq_vals(:)' ];
        badburstlen_vals = allBurstsAllElectrodes{iele}(~ind)*1e3;
        badBurstLen = [badBurstLen, badburstlen_vals(:)' ];
        badburstdur_vals = (1 - tBurstsAllElectrodes{iele}(~ind))*1e3 + badburstlen_vals/2;
        badBurstDur = [badBurstDur, badburstdur_vals(:)' ];
    end
    medianBL_hue = median(meanBL_hue_eles);
%     %% TODO: correct this!!
%     semedianBL_hue = 0.1*median(medianBL_hue_eles); %getSEMedian(medianBL_hue_eles);
    semedianBL_hue = getSEMedian(meanBL_hue_eles);
    % scatter plot _ 1 plot per subject - hue power vs median+-SE gamma
    figure(lenHuesfigure);
    subplot(2,1,subjectNum); hold on;
    medianBLs_HUE{subjectNum} = [medianBLs_HUE{subjectNum}(:)',medianBL_hue];
    scatter(gammapower_hue, medianBL_hue,50, colors(ihue,:),'filled','marker','o');
    xdata = gammapower_hue(:) * [1,1]; ydata = medianBL_hue(:) + semedianBL_hue(:) * [-1, 1];
    %     for xd = 1:size(xdata,1); line(xdata(xd,:), ydata(xd,:), 'color','k' , 'linewidth', 1.5); end
    line(xdata, ydata, 'color','k' , 'linewidth', 1.5);
%     xlim([0 25]); ylim([0 1125]);
    if subjectNum == 2; xlabel('Gamma Power (dB)'); end
    ylabel('Burst length (ms)');
    
    xlim([0 35]); ylim([0 500]);
    % scatter plot _ 1 plot per electrode
    % burst length
    % TF of top 2 trials in top electrode for hue
    
    %% Electrode-wise plots and analysis (plotting done in electrode loop)
    % top 2 trials for top hue for electrode
    % Accumulate values  (to be used in for electrodes loop)
    medianBL_hues_eles = [medianBL_hues_eles, meanBL_hue_eles(:)];
    gammapower_hues_eles = [gammapower_hues_eles; diffPower(:)];
    badBurstLens = [badBurstLens ,badBurstLen(:)' ];
    badBurstDurs = [badBurstDurs ,badBurstDur(:)' ];
    badBurstFreqs = [badBurstFreqs ,badBurstFreq(:)' ];
    goodBurstFreqs = [goodBurstFreqs ,goodBurstFreq(:)' ];
    % scatter plot _ 1 plot per subject - elec power vs median+-SE gamma
    % burst length
end
%% save lenHuesfigure
savebase = fullfile(outputfolder, 'Fig 2_BurstLengths_HueWise');
savefig(lenHuesfigure,[savebase,'.fig'],'compact');
save([savebase,'.mat'],'medianBLs_HUE','meanmedianBLs_HUE','semedianBLs_HUE','-v7.3');
end
