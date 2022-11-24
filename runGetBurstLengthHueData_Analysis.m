cd('./Data_Analysis/');
basePath = '..\GammaHarmonics Segmented Data\'; %pwd; 
alpafigure = []; tutufigure = [];
% For all hues (0:10:350): i = 1 => 0 i = 36 => 350
for i = 1:36
    hue = (i-1)*10;
    Hue = num2str(hue);
    disp(hue);  
%     
    
    plotfigure = tutufigure;
    dataPreProcessing('tutu', 'Color', Hue, basePath, 0,plotfigure);
    
    plotfigure = alpafigure;
    dataPreProcessing('alpa', 'Color', Hue, basePath, 0,plotfigure);
    
    disp(hue);
end
% Achromatic (SForiAchro for alpa and '37th' hue for tutu)
disp('Processing achromatic data')
plotfigure = tutufigure;
dataPreProcessing('tutu', 'Color', '360', basePath, 1, plotfigure);
% 
plotfigure = alpafigure;
dataPreProcessing('alpa', 'SForiAchro', '', basePath, 1, plotfigure);

%%    
function dataPreProcessing(subjectName, expType, stimType, basePath, achroFlag, plotfigure)
    
    if isempty(stimType)
        stimno = 37;
    else
        stimno = str2num(stimType)/10+1;
    end
    
    % Get Corresponding experiment/data
    subjectID = [subjectName expType];
    if strcmp(subjectID,'alpaColor')
        subjectName = 'alpa';expDate = '301215'; protocolName = 'GRF_001'; % 488: Hue fullscreen
    elseif strcmp(subjectID,'tutuColor')
        subjectName = 'tutu'; expDate = '191016'; protocolName = 'GRF_001'; % 111: Hue fullscreen
    elseif strcmp(subjectID,'alpaSForiAchro')
        subjectName = 'alpa'; expDate = '301215'; protocolName = 'GRF_005'; % SFOri - alpa Vinay
    end

    saveFolder = fullfile('..', 'savedDataMP_fullLenMP', 'processedData', subjectName);
    mkdir(saveFolder);
    gridType = 'Microelectrode';
    folderSourceString = fullfile(basePath, expType);
    folderBase = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
    folderLFP = fullfile(folderBase,'segmentedData','LFP');

    stimPeriod = [0.25 0.75];% 500ms
    baselinePeriod = [-0.5 0]; %500 ms

    % TimeVals, FS and StimPos
    load(fullfile(folderLFP,'lfpInfo'),'timeVals');
    Fs = 1./(timeVals(2)-timeVals(1));
    numPoints = round(diff(stimPeriod)*Fs);

    stPos = find(timeVals>=stimPeriod(1),1) + (1:numPoints);
    blPos = find(timeVals>=baselinePeriod(1),1) + (1:numPoints);
    
        
        % Multi-taper parameters
        if ~exist('TW','var'); TW = 2; end
        tw = TW; fmax = 250;
        mt.tapers = [tw (2*tw-1)];
        mt.pad = -1; mt.Fs = Fs;
        mt.fpass = [0 fmax];
        
        %     % load badtrials
        badTrialFile = fullfile(folderBase,'segmentedData','badTrials');
        load(badTrialFile,'badTrials');
        paramCombinationsFile = fullfile(folderBase,'extractedData','parameterCombinations');
        Params = load(paramCombinationsFile);
        
        if strcmp(subjectID,'alpaSForiAchro')
            highRMSElectrodesStruct = load(fullfile(folderSourceString,'analyzeElectrodes', 'alpaLFPElectrodeList.mat'));
            t  = load(fullfile(folderLFP,'lfpInfo.mat'));
            timeVals = t.timeVals;
            
            highRMSElectrodes = highRMSElectrodesStruct.alpaLFPElectrodeList;
            highRMSElectrodes = setdiff(highRMSElectrodes, 4); % 4 is an electrode compared to color cases (65 v 64)
            
            cVals = Params.cValsUnique; oVals = Params.oValsUnique; fVals = Params.fValsUnique;
            o = find(oVals ==90); %90 degree
            goodPos = Params.parameterCombinations{1,1,1,1,o};
            goodPos = setdiff(goodPos, badTrials);
        else
            highRMSElectrodesStruct = load(fullfile(folderSourceString,'analyzeElectrodes',subjectName,'highRMSElectrodes'));
            
            highRMSElectrodes = highRMSElectrodesStruct.highRMSElectrodes;
            
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
        
        % psdST = 3D Array - (:, trialIndex, electrodeIndex) $#
        
        n = 4; % order of butterWorth Filtering of gamma and harmonic signals
        delta = 10; % for filtfilt signal filtering
        gammaRangeHz = [20 80]; % Hz
        maxHarmFreq = 140; %Hz
        GHwidth = 12; % Hz
        
        phaseDiff = [];
        psdST = []; psdBL = [];
        gammaFreq = []; harmonicFreq = []; harmonicFreq_act = [];
        baseCorrLog10PSD = [];
        
        %% Gamma length MP variables
        allBurstsAllElectrodes = cell([length(highRMSElectrodes),1]);
        allBurstzpowersAllElectrodes = cell([length(highRMSElectrodes),1]);
        tBurstsAllElectrodes = cell([length(highRMSElectrodes),1]);
        fBurstsAllElectrodes = cell([length(highRMSElectrodes),1]);
        burstTrialAllElectrodes = cell([length(highRMSElectrodes),1]);
        burstTrialAllElectrodes_allTrials = cell([length(highRMSElectrodes),1]);
        goodBurstAllElectrodes = cell([length(highRMSElectrodes),1]);
        goodBurstAllElectrodes_full = cell([length(highRMSElectrodes),1]);
        numBurstsAllElectrodes = [];
        numBursts=zeros(1,length(highRMSElectrodes));
        medianBurstLength=zeros(1,length(highRMSElectrodes));
        diffPower=[]; zeros()/0;
        oldburstLenghtMP = [];
        
        prevanalogData = [];
        electrodeIDs = highRMSElectrodes;
        for i = 1:length(highRMSElectrodes)
            disp(['elec' num2str(highRMSElectrodes(i))]);
            load(fullfile(folderLFP,['elec' num2str(highRMSElectrodes(i))]),'analogData');
            stLFPP = analogData(goodPos,stPos)';
            
            [psdSTelec,freqVals] = mtspectrumc(stLFPP,mt);
            psdBLelec = mtspectrumc(analogData(goodPos,blPos)',mt);
            psdST = cat(3,psdST,psdSTelec);
            psdBL = cat(3,psdBL,psdBLelec);
            baseCorrLog10PSD = cat(3,baseCorrLog10PSD,log10(psdSTelec)-log10(psdBLelec));
            
            % Find gamma and and following(harmonic) peak
            gammaRangePos = intersect(find(freqVals>=gammaRangeHz(1)),find(freqVals<=gammaRangeHz(2)));
            gammaRange = baseCorrLog10PSD(gammaRangePos,:,i);
            gammaAmpLog = max(gammaRange);
            [gammaPosInd,~] = find(gammaRange==gammaAmpLog);
            gammaPosss = gammaRangePos(gammaPosInd);
            peakGammaFreq = freqVals(gammaPosss);
            
            
            peakHarmonicFreq = 2*peakGammaFreq; % Exact Harmonic (2*fG)
            
            gammaFreq = cat(2,gammaFreq,peakGammaFreq');
            harmonicFreq = cat(2,harmonicFreq,peakHarmonicFreq');
            
            PD = [];
            
            for j =1:length(goodPos)
                [B,A] = butter(n,[peakGammaFreq(j)-delta, peakGammaFreq(j)+delta]/(Fs/2));
                [D,C] = butter(n,[peakHarmonicFreq(j)-delta, peakHarmonicFreq(j)+delta]/(Fs/2));
                
                gammaSignal = filtfilt(B,A,stLFPP(:,j)); %#
                harmonicSignal = filtfilt(D,C,stLFPP(:,j));
                
                % Phasedifference between gamma and harmonic
                G = hilbert(gammaSignal);
                H = hilbert(harmonicSignal);
                PD = cat(2, PD, (angle(H)-2*angle(G)));
            end
            phaseDiff = cat(3, phaseDiff, PD);
            
            %% MP length analysis
            %MP parameters
            thresholdFactor=0.5;
            maxIteration=50;
            adaptiveDictionaryParam=0.9;
            dictionarySize=2500000;
            % Renaming parameters
            gammaFreqRangeHz = gammaRangeHz;
            fullstimulusPeriodS = [0 0.8]; %stimPeriod;
            stimulusPeriodS = [0.1 0.8]; %stimPeriod;
            baselinePeriodS = [-0.5 0]; %baselinePeriod;
            
            runMPAnalysisFlag = true;
                % Save MP results
                tag = [saveFolder,'/'];
                
                % Output folder
                outputFolder = fullfile(tag,'mpAnalysis',['elec' num2str(highRMSElectrodes(i))]);
                mkdir(outputFolder);
                
                % Save gaborInfo as a mat file
                gaborInfoFile = fullfile(outputFolder,['gaborInfo_G30-70_M' num2str(maxIteration) ...
                    '_D' num2str(100*adaptiveDictionaryParam) '_R' num2str(dictionarySize) '_S' num2str(1000*0.25) '-' num2str(1000*0.75) '.mat']);
                
                if ~exist(gaborInfoFile,'file')
                    [~,~,~,gaborInfo,header] = getBurstLengthMP_hue(analogData,[],timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize);
                    save(gaborInfoFile,'gaborInfo','header');
                end
                disp(['Opening saved file ' gaborInfoFile]);
                load(gaborInfoFile);
                
                [burstLengthMP, flist, tlist, alist, ~, ~] = getBurstLengthMP_hue(analogData,goodPos,timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,gaborInfo,header);

            burstTrialSingleElectrode = [];
            burstTrialSingleElectrode_allTrials = [];
            goodBurstSingleElectrode = []; goodBurstSingleElectrode_full = [];
            numBurstsSingleElectrode = [];
            allBurstsSingleElectrode=[];
            allBurstzpowersSingleElectrode = [];
            fBurstsSingleElectrode = [];
            tBurstsSingleElectrode = [];
            numTrials=length(burstLengthMP);
            trialids = 1:size(gaborInfo,1); trialids = trialids(goodPos);
            for iii = 1:numTrials
                burstTrialSingleElectrode = [burstTrialSingleElectrode, iii*ones([1,numel(burstLengthMP{iii})])];
                burstTrialSingleElectrode_allTrials = [burstTrialSingleElectrode_allTrials, trialids(iii)*ones([1,numel(burstLengthMP{iii})])];
                goodBurstSingleElectrode = [goodBurstSingleElectrode, (burstLengthMP{iii}(:)'/2 + tlist{iii}(:)' <= max(stimulusPeriodS)) ...
                    & (tlist{iii}(:)' - burstLengthMP{iii}(:)'/2 >= min(stimulusPeriodS))];
                goodBurstSingleElectrode_full = [goodBurstSingleElectrode_full, (burstLengthMP{iii}(:)'/2 + tlist{iii}(:)' <= max(fullstimulusPeriodS)) ...
                    & (tlist{iii}(:)' - burstLengthMP{iii}(:)'/2 >= min(fullstimulusPeriodS))];
                allBurstsSingleElectrode = [allBurstsSingleElectrode, burstLengthMP{iii}(:)'];
                allBurstzpowersSingleElectrode = [allBurstzpowersSingleElectrode, alist{iii}(:)'];
                fBurstsSingleElectrode = [fBurstsSingleElectrode, flist{iii}(:)'];
                tBurstsSingleElectrode = [tBurstsSingleElectrode, tlist{iii}(:)'];
                numBurstsSingleElectrode = [numBurstsSingleElectrode, numel(burstLengthMP{iii}(:))];
            end
            burstTrialAllElectrodes{i} = burstTrialSingleElectrode;
            burstTrialAllElectrodes_allTrials{i} = burstTrialSingleElectrode_allTrials;
            goodBurstAllElectrodes{i} = goodBurstSingleElectrode;
            goodBurstAllElectrodes_full{i} = goodBurstSingleElectrode_full;
            allBurstsAllElectrodes{i} = allBurstsSingleElectrode;
            allBurstzpowersAllElectrodes{i} = allBurstzpowersSingleElectrode;
            tBurstsAllElectrodes{i} = tBurstsSingleElectrode;
            fBurstsAllElectrodes{i} = fBurstsSingleElectrode;
            numBurstsAllElectrodes= [numBurstsAllElectrodes;numBurstsSingleElectrode(:)'];
            numBursts(1,i)=length(allBurstsSingleElectrode(:))/numTrials;
            medianBurstLength(1,i)=median(allBurstsSingleElectrode);
                
        end
        %% Computing gamma power for plot
        for k = 1:length(highRMSElectrodes)
            powerDB = 10*(log10(mean(psdST(:,:,k),2)) - log10(mean(psdBL(:,:,k),2)));
            
            %         [pk, loc] = findpeaks(powerDB);
            pk = findpeaks(powerDB);
            loc = pk.loc(pk.loc>=1 & pk.loc<=numel(powerDB)); pk = powerDB(loc);
            
            
	    powerDB = powerDB(gammaRangeHz(1)/2+1:gammaRangeHz(2)/2+1);
            pkG = findpeaks(powerDB);
            pkG = powerDB(pkG.loc(pkG.loc>=1 & pkG.loc<=numel(powerDB)));
            
            pGG = max(pkG);
            diffPower = [diffPower, pGG];
        end
        %% - Uncomment this to save
        if achroFlag == 1
            save(fullfile(saveFolder, [subjectName 'Achro' '.mat']),'electrodeIDs', 'allBurstsAllElectrodes','allBurstzpowersAllElectrodes','tBurstsAllElectrodes','fBurstsAllElectrodes','burstTrialAllElectrodes','burstTrialAllElectrodes_allTrials','goodBurstAllElectrodes','goodBurstAllElectrodes_full','numBursts','numBurstsAllElectrodes','medianBurstLength','phaseDiff','gammaFreq', ...
                'diffPower','psdST', 'psdBL', 'goodPos', 'stPos', 'highRMSElectrodes', 'timeVals', 'freqVals', ...
                'subjectName', 'expType', 'stimType', 'expDate', 'protocolName', 'mt', 'Fs', 'gammaRangeHz')
        else
            save(fullfile(saveFolder, [subjectName expType stimType '.mat']), 'electrodeIDs', 'allBurstsAllElectrodes','allBurstzpowersAllElectrodes','tBurstsAllElectrodes','fBurstsAllElectrodes','burstTrialAllElectrodes','burstTrialAllElectrodes_allTrials','goodBurstAllElectrodes','goodBurstAllElectrodes_full','numBursts','numBurstsAllElectrodes','medianBurstLength','phaseDiff','gammaFreq', ...
                'diffPower','psdST', 'psdBL', 'goodPos', 'stPos', 'highRMSElectrodes', 'timeVals', 'freqVals', ...
                'subjectName', 'expType', 'stimType', 'expDate', 'protocolName', 'mt', 'Fs', 'gammaRangeHz')
        end

end
 