cd('./Data_Analysis/');
basePath = '..\GammaHarmonics Segmented Data\'; 
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
outputfolder = fullfile(tag,['TFresults',postfix]); mkdir(outputfolder);
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

% data accumulated in hue loop
medianBL_hues_eles = [];
gammapower_hues_eles = [];
badBurstLens = [];
badBurstDurs = [];
hueTrialIDs = cell([37,1]);

Rpsdlist = cell([37,1]);
fpsd = [];
Gplist = cell([37,1]); 
Gflist = cell([37,1]);
col = [hsv(36);0.6*[1,1,1]];
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
        c = (cVals ==100); % contrast 100
        for o = 1:length(oVals)
            goodPosAll{o} = Params.parameterCombinations{1,1,1,1,o,c,1};
            goodPosAll{o} = setdiff(goodPosAll{o},badTrials);
        end
        stimIndex = str2double(stimType)/10 + 1;
        goodPos = goodPosAll{stimIndex};
        % stimIndex - represents Color number, ranging 1 to 36
        % representing hues 0 to 350 with interval size 10
    end
    hueTrialIDs{ihue} = goodPos;
    %% MP length analysis
    %MP parameters
    thresholdFactor=0.5;
    maxIteration=50;
    adaptiveDictionaryParam=0.9;
    dictionarySize=2500000;
    
    stimulusPeriodS = [0.1 0.7]; %stimPeriod;
    baselinePeriodS = [-0.5 0]; %baselinePeriod;
    
    tag = [saveFolder,'/'];
    
    outputFolder_hue = fullfile(outputfolder,'Hue-wise',['hue' num2str(ihue)]); mkdir(outputFolder_hue);
    
    if ihue == 37
        load(fullfile(saveFolder, [subjectName 'Achro' '.mat']), 'allBurstsAllElectrodes','tBurstsAllElectrodes','fBurstsAllElectrodes','burstTrialAllElectrodes','goodBurstAllElectrodes', ...
            'diffPower','gammaFreq', 'psdST', 'psdBL', 'highRMSElectrodes', 'timeVals', 'freqVals', ...
            'subjectName', 'expType')
    else
        load(fullfile(saveFolder, [subjectName expType stimType '.mat']), 'allBurstsAllElectrodes','tBurstsAllElectrodes','fBurstsAllElectrodes','burstTrialAllElectrodes','goodBurstAllElectrodes', ...
            'diffPower','gammaFreq', 'psdST', 'psdBL', 'highRMSElectrodes', 'timeVals', 'freqVals', ...
            'subjectName', 'expType')
    end
    
    
    %% Hue-wise plots and analysis
    gammapower_hue = mean(diffPower);
    medianBL_hue_eles = []; semedianBL_hue_eles = [];
    badBurstLen = []; badBurstDur = [];
    for iele = numel(highRMSElectrodes):-1:1%1:numel(highRMSElectrodes)
        %% only good bursts? or all bursts?
        if selgoodBursts
%             ind = (tBurstsAllElectrodes{iele}>=stimulusPeriodS(1) & tBurstsAllElectrodes{iele}<=stimulusPeriodS(end)); %goodBurstAllElectrodes{iele}==1;
%             ind = ind & (fBurstsAllElectrodes{iele}>=40 & fBurstsAllElectrodes{iele}<=60);
            ind = (tBurstsAllElectrodes{iele}(:)>=0.2 & tBurstsAllElectrodes{iele}(:)<=0.6); %goodBurstAllElectrodes{iele}==1;
            ind = ind & (fBurstsAllElectrodes{iele}(:)>=30 & fBurstsAllElectrodes{iele}(:)<=70);
            ind = ind & (tBurstsAllElectrodes{iele}(:)-allBurstsAllElectrodes{iele}(:)/2>=0 & tBurstsAllElectrodes{iele}(:)+allBurstsAllElectrodes{iele}(:)/2<=0.8);
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
        medianBL_hue_eles = [medianBL_hue_eles, median(medianBL_trial)];
%         semedianBL_hue_eles = [semedianBL_hue_ele, getSEMedian(medianBL_trial)];
        badburstlen_vals = allBurstsAllElectrodes{iele}(~ind)*1e3;
        badBurstLen = [badBurstLen, badburstlen_vals(:)' ];
        badburstdur_vals = (1 - tBurstsAllElectrodes{iele}(~ind))*1e3 + badburstlen_vals/2;
        badBurstDur = [badBurstDur, badburstdur_vals(:)' ];
    end
    medianBL_hue = median(medianBL_hue_eles);
    Rpsdlist{ihue} = 10*mean(log10(mean(psdST,2))-log10(mean(psdBL,2)),3)';
    fpsd = freqVals;
    %mean(log10(psdST)-log10(psdBL));
    Gplist{ihue} = mean(diffPower,'all');%ean(log10(psdST)-log10(psdBL));
    Gflist{ihue} = mean(gammaFreq,'all');%ean(log10(psdST)-log10(psdBL));
%     %% TODO: correct this!!
%     semedianBL_hue = 0.1*median(medianBL_hue_eles); %getSEMedian(medianBL_hue_eles);
%     semedianBL_hue = getSEMedian(medianBL_hue_eles);
    % scatter plot _ 1 plot per subject - hue power vs median+-SE gamma
%     figure(lenHuesfigure);
%     subplot(2,1,subjectNum); hold on;
    
%     scatter(gammapower_hue, medianBL_hue,50, colors(ihue,:),'filled','marker','o');
%     xdata = gammapower_hue(:) * [1,1]; ydata = medianBL_hue(:) + semedianBL_hue(:) * [-1, 1];
%     %     for xd = 1:size(xdata,1); line(xdata(xd,:), ydata(xd,:), 'color','k' , 'linewidth', 1.5); end
%     line(xdata, ydata, 'color','k' , 'linewidth', 1.5);
%     xlim([0 25]); ylim([0 1125]);
%     if subjectNum == 2; xlabel('Gamma Power (dB)'); end
%     ylabel('Burst length (ms)');
   
    % burst length mean+-se; histogram/scatter -vs freq/starttime; 
    % avg TF across all trials and electrodes hue: consult avging method in
    % hues paper
    %% (Hue,electrode)-wise plots and analysis (plotting done in electrode loop)
    
    % median +- se - burst length in title; TF spectrum: avg of all trials
    % TF of each trial stored as an image in hue/electrode and electrode/hue folders.
    %% Electrode-wise plots and analysis (plotting done in electrode loop)
    % top 2 trials for top hue for electrode
    % Accumulate values  (to be used in for electrodes loop)
    medianBL_hues_eles = [medianBL_hues_eles, medianBL_hue_eles(:)];
    gammapower_hues_eles = [gammapower_hues_eles; diffPower(:)];
    badBurstLens = [badBurstLens ,badBurstLen(:)' ];
    badBurstDurs = [badBurstDurs ,badBurstDur(:)' ];
    % scatter plot _ 1 plot per subject - elec power vs median+-SE gamma
    % burst length
end
%% save lenHuesfigure
savebase = fullfile(outputfolder, 'BurstLength hue_wise');
%% variables accumulating data across electrodes per hue
meanE_hue = cell([37,1]);

%%
for iele = numel(highRMSElectrodes)-(0:7) %1:numel(highRMSElectrodes)
    % output folders
%     outputfolder = fullfile(tag,'LengthAnalysis'); mkdir(outputfolder);
    outputFolder_ele = fullfile(outputfolder,'Electrode-wise',['elec' num2str(highRMSElectrodes(iele))]); mkdir(outputFolder_ele);
    %% Electrode-wise plots and analysis
    gaborFolder = fullfile(tag,'mpAnalysis',['elec' num2str(highRMSElectrodes(iele))]);
    gaborInfoFile = fullfile(gaborFolder,['gaborInfo_G30-70_M' num2str(maxIteration) ...
        '_D' num2str(100*adaptiveDictionaryParam) '_R' num2str(dictionarySize) '_S' num2str(1000*0.25) '-' num2str(1000*0.75) '.mat']);
    load(gaborInfoFile); % gaborInfo, header
    
%     % avg of all trials 
    [~, freqVals] = getEnergyMP3p1(gaborInfo(1,:,:), header(1,:),timeVals);
%     fig_tf_ele   = figure('InvertHardcopy',false);
%     BaseE = mean(meanE_ele(:,timeVals>=baselinePeriodS(1) & timeVals<=baselinePeriodS(2)),1);
%     pcolor(timeVals, freqVals, 10*(log10(meanE_ele)-log10(BaseE))); shading interp; ylim([0 100]); 
%             colormap jet; caxis([-5 5]);
%     % mean +- se - burst length in title
    
    %% (Hue,electrode)-wise plots and analysis (plotting done in electrode loop)
    % Accumulate values (to be used in for electrodes loop)
    % scatter-plot burst lengths of different colors in each electrode (1 plot per electrode)
    % burst length
    % TF of avg across trials
    if isempty(gaborInfo)
        gaborInfo = [];
    end
    if isempty(header)
        header = [];
    end
    %% TF plots
    clear burstTrialAllElectrodes_allTrials burstTrialAllElectrodes
%     parfor (ihue = 1:numel(hueTrialIDs), 6)
    hueIDs = [1, 22, 37];
    subno = 0;
    fignos = [1];
    huefigs = [figure('WindowState','maximized','InvertHardcopy','on')];
    for ihue = hueIDs(:)' %1:numel(hueTrialIDs)
        %% -
        
        if subno==2
            subno = 1;
            fignos = [fignos,fignos(end) + 1];
            huefigs = [huefigs,figure('WindowState','maximized','InvertHardcopy','on')];
        else
            subno = subno + 1;
        end
        
        if ihue == 37
            fv = load(fullfile(saveFolder, [subjectName 'Achro' '.mat']), 'allBurstsAllElectrodes','tBurstsAllElectrodes','fBurstsAllElectrodes','burstTrialAllElectrodes','burstTrialAllElectrodes_allTrials', 'timeVals','goodBurstAllElectrodes')
        else
            stimType = num2str((ihue-1)*10); expType = 'Color';
            fv = load(fullfile(saveFolder, [subjectName expType stimType '.mat']), 'allBurstsAllElectrodes','tBurstsAllElectrodes','fBurstsAllElectrodes','burstTrialAllElectrodes','burstTrialAllElectrodes_allTrials', 'timeVals','goodBurstAllElectrodes' )
        end
        
        [meanE_hue_ele,freqVals] = getEnergyMP3p1(gaborInfo(hueTrialIDs{ihue}(:),:,:), header(hueTrialIDs{ihue}(:),:),fv.timeVals);
%         fig_tf_hue_ele   = figure('WindowState','maximized','InvertHardcopy',false);
        figure(huefigs(end));
        
        subplot(2,10,(subno-1)*5 + [1:3]);
        BaseE = mean(meanE_hue_ele(:,fv.timeVals>=baselinePeriodS(1) & fv.timeVals<=baselinePeriodS(2)),2);
%         pcolor(fv.timeVals, freqVals, 10*(log10(meanE_hue_ele)-log10(BaseE))); shading interp; ylim([20 80]);
%         tsel = fv.timeVals>=0.1 & fv.timeVals<=0.8;
%         fsel = find(freqVals>=20 & freqVals <=75); fsel = fsel(:)';
%         [~,pkfreq_t] = max(log10(meanE_hue_ele(fsel,tsel))-log10(BaseE(fsel)),[],1);
%         pkfreq_t = freqVals(fsel(pkfreq_t));
%         t_pkf = fv.timeVals(tsel);
%         hold on; plot(t_pkf, pkfreq_t,'color',[[1,1,1]*0.7],'LineWidth',2);
%         xlim([-0.2 1]); %0.8]);
%         colormap jet; 
        %% burst marking
        Trials_hue = 1:size(gaborInfo,1); Trials_hue = Trials_hue(hueTrialIDs{ihue});
        brtrids = find(any(fv.burstTrialAllElectrodes{iele}(:) == 1:numel(Trials_hue),1));
        for ii = randi(numel(brtrids))
            i = brtrids(ii)+1;
%             subplot(2,3,i);
%             selTrial = hueTrialIDs{ihue}(2*(i-2)+1);
            selTrial_act = Trials_hue(i-1);
            selTrial = i-1;
            
            [E, freqVals] = getEnergyMP3p1(gaborInfo(selTrial_act,:,:), header(selTrial,:),fv.timeVals);
%             BaseE = mean(E(:,fv.timeVals>=baselinePeriodS(1) & fv.timeVals<=baselinePeriodS(2)),2);
            pcolor(fv.timeVals, freqVals, 10*(log10(E)-log10(BaseE))); shading interp;
            axis tight;
            xlim([-0.2 0.8]); %0.8]);
 
            colormap jet; %caxis([-15 30]);
            hold on;
            burstinds = find((fv.burstTrialAllElectrodes{iele}(:) == selTrial) ...
                & fv.tBurstsAllElectrodes{iele}(:)-fv.allBurstsAllElectrodes{iele}(:)/2>=0 ...
                & fv.tBurstsAllElectrodes{iele}(:)+fv.allBurstsAllElectrodes{iele}(:)/2<=0.8 ...
                & fv.fBurstsAllElectrodes{iele}(:)>=20 ...
                & fv.fBurstsAllElectrodes{iele}(:)<=80); 
                %& fv.goodBurstAllElectrodes{iele});
%             burstinds2 = find(fv.burstTrialAllElectrodes_allTrials{iele} == selTrial_act ...
%                 & fv.goodBurstAllElectrodes{iele});
%             assert(all(burstinds==burstinds2),'TrialID mismatch!')
%             assert(all(burstinds==burstinds2))
            for n = 1:numel(burstinds)
                xdata = fv.tBurstsAllElectrodes{iele}(burstinds(n)) + [-1,1]*fv.allBurstsAllElectrodes{iele}(burstinds(n))/2; 
                ydata = [1,1]*fv.fBurstsAllElectrodes{iele}(burstinds(n));
                
                line(xdata, ydata,'linewidth',20,'color', [[1,0.5,1]*1],'linestyle','-');
            end
        end
%         ylim([20 80]); ylabel('Frequency (Hz)');
        %%
        caxis([0 30]);
        if strcmp(subjectName,'alpa')
            caxis([0 20]);
            if ihue==37
                caxis([0 15]);
            end
        elseif ihue==37
            caxis([0 15]);
%             if strcmp(subjectName,'alpa')
%             caxis([0 5]);
%             if ihue==37
%                 caxis([0 3]);
%             end
%         elseif ihue==37
%             caxis([0 10]);
        end

        c = colorbar('location','east'); 
        c.Label.String = 'dB';  
        c.Position = c.Position + [0.125,-0.0125,0,0];
        c.AxisLocation = 'out';
        xlabel('Time (s)'); 
        if subno == 1 && fignos(end)==1
            ylabel('Frequency (Hz)');
        end
        title(['Hue ',num2str((ihue-1)*10),'\circ'], 'fontweight','bold');
        if ihue==37
            title(['Achromatic grating'], 'fontweight','bold');
        end
        ylim([20 80]);
%         subplot(2,10,(subno-1)*5 + [4]);
% %         Rpsd = mean(log10(psdST))-mean(log10(psdBL)); fpsd = ((1:numel(Rpsd))-1)*2;
%         plot(Rpsdlist{ihue},fpsd,'linewidth',2);
% %         Rpsdlist{ihue} = Rpsdlist{ihue}+Rpsd/numel(highRMSElectrodes);
%         axis tight
%         ylim([20 80]);
        if strcmp(subjectName,'alpa')
            sgtitle(['M1: single trial example'], 'fontweight','bold');
        elseif strcmp(subjectName,'tutu')
            sgtitle(['M2: single trial example'], 'fontweight','bold');
        end 
        
        %% burst length histogram
        allElecallBurstsAllElectrodes=[]
        for iele2 = 1:numel(highRMSElectrodes)
        selTrial = 1:numel(Trials_hue);
        selTrial_act = Trials_hue(selTrial);
        burstinds = find((any(fv.burstTrialAllElectrodes{iele2}(:) == selTrial(:)',2)) ...
                & fv.tBurstsAllElectrodes{iele2}(:) - fv.allBurstsAllElectrodes{iele2}(:)/2>=0 ...
                & fv.tBurstsAllElectrodes{iele2}(:) + fv.allBurstsAllElectrodes{iele2}(:)/2<=0.8 ...
                & fv.tBurstsAllElectrodes{iele2}(:)>=0.2 ...
                & fv.tBurstsAllElectrodes{iele2}(:)<=0.6...
                & fv.fBurstsAllElectrodes{iele2}(:)>=20 ...
                & fv.fBurstsAllElectrodes{iele2}(:)<=80); %& fv.goodBurstAllElectrodes{iele2}(:)); 
            burstinds2 = find(any(fv.burstTrialAllElectrodes_allTrials{iele2}(:) == selTrial_act(:)',2) ...
                & fv.tBurstsAllElectrodes{iele2}(:) - fv.allBurstsAllElectrodes{iele2}(:)/2>=0 ...
                & fv.tBurstsAllElectrodes{iele2}(:) + fv.allBurstsAllElectrodes{iele2}(:)/2<=0.8 ...
                & fv.tBurstsAllElectrodes{iele2}(:)>=0.2 ...
                & fv.tBurstsAllElectrodes{iele2}(:)<=0.6...
                & fv.fBurstsAllElectrodes{iele2}(:)>=20 ...
                & fv.fBurstsAllElectrodes{iele2}(:)<=80); 
%                 & fv.goodBurstAllElectrodes{iele2}(:));
            assert(all(burstinds==burstinds2),'TrialID mismatch!')

         allElecallBurstsAllElectrodes=[allElecallBurstsAllElectrodes,fv.allBurstsAllElectrodes{iele2}(burstinds(:)')*1e3];
         end
         subplot(2,10,10+(subno-1)*5 + [1:4]);
         
         %          histogram(fv.allBurstsAllElectrodes{iele}(burstinds(:)')*1e3,'normalization','pdf','BinWidth',50,'FaceColor',col(ihue,:),'EdgeColor','k');
         histogram(allElecallBurstsAllElectrodes,'normalization','pdf','BinWidth',50,'FaceColor',col(ihue,:),'EdgeColor','k');
         xlabel('Burst Lengths (ms)');
         xlim([0 800]);
          
         if subno == 1 && fignos(end)==1
             ylabel('pdf');
         end
        fig_tf_hue_ele = huefigs(end); figure(fig_tf_hue_ele);
        f = fig_tf_hue_ele; 
        savefig(fig_tf_hue_ele,[fullfile(outputFolder_ele, ['huefig' num2str(fignos(end))]) '.fig'],'compact');
        print(fig_tf_hue_ele,[fullfile(outputFolder_ele, ['huefig' num2str(fignos(end))]) '.png'],'-dpng','-r300');
        outputFolder_hue = fullfile(outputfolder,'Hue-wise',['huefig' num2str(fignos(end))]); mkdir(outputFolder_hue);
        print(fig_tf_hue_ele,[fullfile(outputFolder_hue, ['ele' num2str(highRMSElectrodes(iele))]) '.png'],'-dpng','-r300');
        
%         savefig(fig_tf_hue_ele,[fullfile(outputFolder_ele, ['hue' num2str(ihue)]) '.fig'],'compact');
%         print(fig_tf_hue_ele,[fullfile(outputFolder_ele, ['hue' num2str(ihue)]) '.png'],'-dpng','-r300');
%         outputFolder_hue = fullfile(outputfolder,'Hue-wise',['hue' num2str(ihue)]); mkdir(outputFolder_hue);
%         print(fig_tf_hue_ele,[fullfile(outputFolder_hue, ['ele' num2str(highRMSElectrodes(iele))]) '.png'],'-dpng','-r300');
%         close(fig_tf_hue_ele);
        
        if isempty(meanE_hue{ihue})
            meanE_hue{ihue} = 0;
        end
        meanE_hue{ihue} = meanE_hue{ihue} + meanE_hue_ele/numel(highRMSElectrodes);
    end
    for ifig = 1:numel(huefigs)
        close(huefigs(ifig));
    end
end

%% save lenElecsfigure
savebase = fullfile(outputfolder, 'BurstLength electrode_wise');
% savefig(lenElesfigure,[savebase,'.fig'],'compact');

%% Hue-wise TF plots
hueIDs = [1, 22, 37];
subno = 0;
fignos = [1];
huefigs = [figure('WindowState','maximized','InvertHardcopy','on')];
    col = [hsv(36);0.6*[1,1,1]];    
for ihue = hueIDs(:)' %1:37

    if subno==2
        subno = 1;
        fignos = [fignos,fignos(end) + 1];
        huefigs = [huefigs,figure('WindowState','maximized','InvertHardcopy','on')];
    else
        subno = subno + 1;
    end

    fig_tf_hue = huefigs(end); %figure('InvertHardcopy',false);
    BaseE = mean(meanE_hue{ihue}(:,timeVals>=baselinePeriodS(1) & timeVals<=baselinePeriodS(2)),2);
%     subplot(2,2,subno);
    subplot(2,10,(subno-1)*5 + [1:3]);
    pcolor(timeVals, freqVals, 10*(log10(meanE_hue{ihue})-log10(BaseE))); shading interp; ylim([20 80]); 
    tsel = timeVals>=0.1 & timeVals<=0.9;
        fsel = find(freqVals>=20 & freqVals <=75); fsel = fsel(:)';
        [~,pkfreq_t] = max(log10(meanE_hue{ihue}(fsel,tsel))-log10(BaseE(fsel)),[],1);
        pkfreq_t = freqVals(fsel(pkfreq_t));
%         pkfreq_t = freqVals(pkfreq_t);
        t_pkf = timeVals(tsel);
        hold on; plot(t_pkf, pkfreq_t,'color',[[1,1,1]*0.7],'LineWidth',2);
        xlim([-0.2 0.8]);
        colormap jet; 
        caxis([0 30]);
        if strcmp(subjectName,'alpa')
            caxis([0 15]);
            if ihue==37
                caxis([0 15]);
            end
        elseif ihue==37
            caxis([0 15]);
        end
        c = colorbar('location','east'); 
        c.Label.String = 'dB';  
        c.Position = c.Position + [0.125,-0.0125,0,0];
        c.AxisLocation = 'out';
        xlabel('Time (s)'); 
        if subno == 1 && fignos(end)==1
            ylabel('Frequency (Hz)');
        end
        title(['Hue ',num2str((ihue-1)*10),'\circ'], 'fontweight','bold');
        if ihue==37
            title(['Achromatic grating'], 'fontweight','bold');
        end
        
        
        subplot(2,10,(subno-1)*5 + [4]);
        Rpsd = Rpsdlist{ihue};
        plot(Rpsd,fpsd,'linewidth',2,'Color',col(ihue,:));
        hold on;
        scatter(Gplist{ihue},Gflist{ihue},'kx');
        axis tight
        xl = xlim; xlim(xl(1) + [0 diff(xl)*1.1]);
        ylim([20 80]); xlabel('PSD (dB)'); title(' ','fontweight','bold');
        ax = gca; ax.YAxis.TickLabels={};
% %         ax = subplot(2,10,10 + (subno-1)*5 + [1:4]);
% % %         plot(t_pkf(2:end),-movmean(diff(pkfreq_t)./diff(t_pkf),ceil(0.15/mean(diff(t_pkf)))*[1,1]),'LineWidth',10,'Color',col(ihue,:)); ax.YScale = 'log';
% % %         plot(t_pkf(2:end),-movmean(diff(pkfreq_t)./diff(t_pkf),ceil(0.15/mean(diff(t_pkf)))*[1,1]),'LineWidth',10,'Color',col(ihue,:)); ax.YScale = 'log';
% %         movwin = ceil(0.15/mean(diff(t_pkf)))*[1,1];
% %         plot(t_pkf(3:end),-1e3*...
% %                            (movmean(diff(pkfreq_t(2:end)),movwin)...
% %                            ./diff(t_pkf(2:end)))...
% %             ./(movmean(diff(movmean(diff(pkfreq_t(1:end)),movwin)),movwin)...
% %                            ./(diff(t_pkf(2:end)).^2)),'LineWidth',2,'Color',col(ihue,:)); %ax.YScale = 'log';
% %         xlim([-0.2 0.8]); %ylim([-1,1]*1/0.2);
% %         ylabel('Time Constant (ms)');
        
        subplot(2,10,10+(subno-1)*5 + [1:4]);
         
         %          histogram(fv.allBurstsAllElectrodes{iele}(burstinds(:)')*1e3,'normalization','pdf','BinWidth',50,'FaceColor',col(ihue,:),'EdgeColor','k');
         histogram(allElecallBurstsAllElectrodes,'normalization','pdf','BinWidth',50,'FaceColor',col(ihue,:),'EdgeColor','k');
         xlabel('Burst Durations (ms)');
         xlim([0 800]);
%          if subno==2
%              title('Distribution of Gamma burst durations','fontweight','bold');
%          end
         if subno==1
             ylabel('probability density','fontweight','bold');
         end


        if strcmp(subjectName,'alpa')
            sgtitle(['M1'], 'fontweight','bold');
        elseif strcmp(subjectName,'tutu')
            sgtitle(['M2'], 'fontweight','bold');
        end 
%         xlim([-0.2 0.8]);      
%     colormap jet; caxis([-15 30]);
%     savefig(fig_tf_hue, fullfile(outputFolder_hue,['hue' num2str(ihue) '.fig']),'compact');
%     print(fig_tf_hue, fullfile(outputFolder_hue,['hue' num2str(ihue) '.png']),'-dpng','-r300');
    
    f = fig_tf_hue;
    
    outputFolder_hue = fullfile(outputfolder,'Hue-wise'); mkdir(outputFolder_hue);
    savefig(fig_tf_hue, fullfile(outputFolder_hue,['huefig' num2str(fignos(end)) '.fig']),'compact');
    print(fig_tf_hue, fullfile(outputFolder_hue,['huefig' num2str(fignos(end)) '.png']),'-dpng','-r300');
% %     close(fig_tf_hue);
%     f = fig_tf_hue;
%     postformatFig;
end
end
