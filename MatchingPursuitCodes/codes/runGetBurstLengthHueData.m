basePath = '.\GammaHarmonics Segmented Data\'; %pwd; 

alpafigure = figure; sgtitle('alpa','fontsize',14);

tutufigure = figure; sgtitle('tutu','fontsize',14);        
% For all hues (0:10:350): i = 1 => 0 i = 36 => 350
for i = 1:36
    hue = (i-1)*10;
    Hue = num2str(hue);
    disp(hue);
    
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

% plotfigure = alpafigure;
% dataPreProcessing('alpa', 'SForiAchro', '', basePath, 1, plotfigure);

    
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

    saveFolder = fullfile(pwd, 'savedDataMP', 'processedData', subjectName);
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
    
%     if ~exist(fullfile(saveFolder, [subjectName 'Achro' '.mat']),'file') && ~exist(fullfile(saveFolder, [subjectName expType stimType '.mat']),'file')
        
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
            %         c = find(cVals ==100); % contrast 100
            o = find(oVals ==90); %90 degree
            %         f = find(round(fVals) == 2); % SF = 2
            goodPos = Params.parameterCombinations{1,1,1,1,o};
            goodPos = setdiff(goodPos, badTrials);
        else
            highRMSElectrodesStruct = load(fullfile(folderSourceString,'analyzeElectrodes',subjectName,'highRMSElectrodes'));
            %         t  = load(fullfile(folderName,'lfpInfo.mat'));
            %         timeVals = t.timeVals;
            
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
        gammaRangeHz = [30 70]; % Hz
        maxHarmFreq = 140; %Hz
        GHwidth = 12; % Hz
        
        phaseDiff = [];
        psdST = []; psdBL = [];
        gammaFreq = []; harmonicFreq = []; harmonicFreq_act = [];
        baseCorrLog10PSD = [];
        
        %% Gamma length MP variables
        allBurstsAllElectrodes = cell([length(highRMSElectrodes),1]);
        numBurstsAllElectrodes = [];
        numBursts=zeros(1,length(highRMSElectrodes));
        medianBurstLength=zeros(1,length(highRMSElectrodes));
        diffPower=[]; zeros()/0;
        oldburstLenghtMP = [];
        
        prevanalogData = [];
        for i = 1:length(highRMSElectrodes)
%         for i = 1:3 %length(highRMSElectrodes)
            %         load(fullfile(folderLFP,['elec' num2str(i)]),'analogData');
            disp(['elec' num2str(highRMSElectrodes(i))]);
            load(fullfile(folderLFP,['elec' num2str(highRMSElectrodes(i))]),'analogData');
            stLFPP = analogData(goodPos,stPos)';
            
            [psdSTelec,freqVals] = mtspectrumc(stLFPP,mt);
            psdBLelec = mtspectrumc(analogData(goodPos,blPos)',mt);
            psdST = cat(3,psdST,psdSTelec);
            psdBL = cat(3,psdBL,psdBLelec);
            baseCorrLog10PSD = cat(3,baseCorrLog10PSD,log10(psdSTelec)-log10(psdBLelec));
            
            % Find gamma and and following(harmonic) peak
            gammaRangePos = intersect(find(freqVals>=gammaRangeHz(1)),find(freqVals<gammaRangeHz(2)));
            gammaRange = baseCorrLog10PSD(gammaRangePos,:,i);
            gammaAmpLog = max(gammaRange);
            %         diffPower = [diffPower, 20*gammaAmpLog];
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
            stimulusPeriodS = stimPeriod;
            baselinePeriodS = baselinePeriod;
            
            runMPAnalysisFlag = true;
%             if runMPAnalysisFlag
                % Save MP results
                folderNameMain = fullfile('data',subjectName,gridType,expDate,protocolName);
                tag = [saveFolder,'/'];
                
                % Output folder
                outputFolder = fullfile(tag,'mpAnalysis',['elec' num2str(highRMSElectrodes(i))]);
                mkdir(outputFolder);
                
                % Save gaborInfo as a mat file
                gaborInfoFile = fullfile(outputFolder,['gaborInfo_G' num2str(gammaFreqRangeHz(1)) '-' num2str(gammaFreqRangeHz(2)) '_M' num2str(maxIteration) ...
                    '_D' num2str(100*adaptiveDictionaryParam) '_R' num2str(dictionarySize) '_S' num2str(1000*stimulusPeriodS(1)) '-' num2str(1000*stimulusPeriodS(2)) '.mat']);
                
                if ~exist(gaborInfoFile,'file')
                    [burstLengthMP,~,~,gaborInfo,header] = getBurstLengthMP_hue(analogData,[],timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize);
                    save(gaborInfoFile,'gaborInfo','header');
                end
                disp(['Opening saved file ' gaborInfoFile]);
                load(gaborInfoFile);
                
                burstLengthMP = getBurstLengthMP_hue(analogData,goodPos,timeVals,sqrt(thresholdFactor),0,stimulusPeriodS,baselinePeriodS,gammaFreqRangeHz,maxIteration,adaptiveDictionaryParam,dictionarySize,gaborInfo,header);
%                 oldburstLenghtMP = [oldburstLenghtMP; burstLengthMP];
%             else
%                 burstLengthMP=[];
%             end

            numBurstsSingleElectrode = [];
            allBurstsSingleElectrode=[];
            numTrials=length(burstLengthMP);
            %         for ii=1:numMethods
%             for iii=1:numTrials
%                 allBurstsSingleElectrode=cat(1,allBurstsSingleElectrode,burstLengthMP{iii}(:));
%             end
            for iii = 1:numTrials
                allBurstsSingleElectrode = [allBurstsSingleElectrode, burstLengthMP{iii}(:)'];
                numBurstsSingleElectrode = [numBurstsSingleElectrode, numel(burstLengthMP{iii}(:))];
            end
            allBurstsAllElectrodes{i} = allBurstsSingleElectrode;
            numBurstsAllElectrodes= [numBurstsAllElectrodes;numBurstsSingleElectrode(:)'];
            numBursts(1,i)=length(allBurstsSingleElectrode(:))/numTrials;
            medianBurstLength(1,i)=median(allBurstsSingleElectrode);
            %         end
            %%
            prevanalogData = analogData;
                
        end
        %% Computing gamma power for plot
        for k = 1:length(highRMSElectrodes)
%         for k = 1:3 %length(highRMSElectrodes)
            powerDB = 10*(log10(mean(psdST(:,:,k),2)) - log10(mean(psdBL(:,:,k),2)));
            
            %         [pk, loc] = findpeaks(powerDB);
            pk = findpeaks(powerDB);
            loc = pk.loc(pk.loc>1 & pk.loc<numel(powerDB)); pk = powerDB(loc);
            
            %         pkG = findpeaks(powerDB(gammaRangeHz(1)/2+1:gammaRangeHz(2)/2+1));
            powerDB = powerDB(gammaRangeHz(1)/2+1:gammaRangeHz(2)/2+1);
            pkG = findpeaks(powerDB);
            pkG = powerDB(pkG.loc(pkG.loc>1 & pkG.loc<numel(powerDB)));
            
            pGG = max(pkG);
            %         fG = freqVals(loc((pk == pGG)));
            %         lG = loc((pk == pGG));
            % %%
            %         if k == 27
            %             disp(['k: ',num2str(k)])
            %         end
            %         pkH = pk(intersect(find(loc > lG+GHwidth/2), find(loc < maxHarmFreq/2+1)));
            %         pHH = max(pkH);
            %         fH = freqVals(loc((pk == pHH)));
            
            %         ratioGH(k) = fH/fG;
            %         freqGamma(k) = fG;
            %         freqHarmonic(k) = fH;
            diffPower = [diffPower, pGG];
        end
        %% - Uncomment this to save
        if achroFlag == 1
            save(fullfile(saveFolder, [subjectName 'Achro' '.mat']), 'allBurstsAllElectrodes','numBursts','numBurstsAllElectrodes','medianBurstLength','phaseDiff','gammaFreq', ...
                'diffPower','psdST', 'psdBL', 'goodPos', 'stPos', 'highRMSElectrodes', 'timeVals', 'freqVals', ...
                'subjectName', 'expType', 'stimType', 'expDate', 'protocolName', 'mt', 'Fs', 'gammaRangeHz')
        else
            save(fullfile(saveFolder, [subjectName expType stimType '.mat']), 'allBurstsAllElectrodes','numBursts','numBurstsAllElectrodes','medianBurstLength','phaseDiff','gammaFreq', ...
                'diffPower','psdST', 'psdBL', 'goodPos', 'stPos', 'highRMSElectrodes', 'timeVals', 'freqVals', ...
                'subjectName', 'expType', 'stimType', 'expDate', 'protocolName', 'mt', 'Fs', 'gammaRangeHz')
        end
%     else
%         if achroFlag == 1
%             load(fullfile(saveFolder, [subjectName 'Achro' '.mat']), 'allBurstsAllElectrodes','numBursts','medianBurstLength','phaseDiff','gammaFreq', ...
%                 'diffPower', 'expType', 'stimType');
% %             load(fullfile(saveFolder, [subjectName 'Achro' '.mat']), 'allBurstsAllElectrodes','numBursts','medianBurstLength','phaseDiff','gammaFreq', ...
% %                 'diffPower','psdST', 'psdBL', 'goodPos', 'stPos', 'highRMSElectrodes', 'timeVals', 'freqVals', ...
% %                 'subjectName', 'expType', 'stimType', 'expDate', 'protocolName', 'mt', 'Fs', 'gammaRangeHz')
%         else
%             load(fullfile(saveFolder, [subjectName expType stimType '.mat']), 'allBurstsAllElectrodes','numBursts','medianBurstLength','phaseDiff','gammaFreq', ...
%                 'diffPower', 'expType', 'stimType');
% %             load(fullfile(saveFolder, [subjectName expType stimType '.mat']), 'allBurstsAllElectrodes','numBursts','medianBurstLength','phaseDiff','gammaFreq', ...
% %                 'diffPower','psdST', 'psdBL', 'goodPos', 'stPos', 'highRMSElectrodes', 'timeVals', 'freqVals', ...
% %                 'subjectName', 'expType', 'stimType', 'expDate', 'protocolName', 'mt', 'Fs', 'gammaRangeHz')
%         end
%     end
    %% plot gamma length
    histC = 0.025:0.025:2;
    if exist('plotfigure','var')
        plotcolors = [hsv(36);0.6*[1,1,1]];
        figure(plotfigure);
        subplot(211);
        hold on; disp(num2str(stimno));
%%         scatter(diffPower,medianBurstLength,[],plotcolors(stimno,:),'o','filled'); hold on;
%         scatter(diffPower,medianBurstLength(1:numel(diffPower)),[],plotcolors(stimno,:),'o','filled'); hold on;
        scatter(diffPower,medianBurstLength(1:numel(diffPower)),[],plotcolors(stimno,:),'.'); hold on;
%         axis([0 40 0 0.5]);
        xlim([0, 35]);
        
        subplot(2,1,2);
        xdata = mean(diffPower);
%%         ydata = median(medianBurstLength);
        ydata = median(medianBurstLength(1:numel(diffPower)));
        scatter(xdata, ydata,[],plotcolors(stimno,:),'d','filled'); hold on;
%%         [semedian, ~] = getSEMedian(medianBurstLength);
        [semedian, ~] = getSEMedian(medianBurstLength(1:numel(diffPower)));
        
        line(xdata*[1;1], ydata+semedian*[1,-1],'color','k');
%         subplot(212);
%         histVals = hist(allBurstsAllElectrodes,histC);
%         nHistVals = histVals/sum(histVals);
%         hold on;
%         plot(histC,nHistVals,'color','r'); hold on;
%         disp([expType,stimType,' : Method: ' MP ', median burst length: ' num2str(median(medianBurstLength(i,:))) ', se: ' num2str(getSEMedian(medianBurstLength(i,:)))]);
    end
end
 