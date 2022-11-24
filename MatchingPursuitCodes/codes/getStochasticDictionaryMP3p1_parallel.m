% This program uses MP version 3.1 by Piotr Durka's group. This is the
% original stochastic dictionary code used in their 2001 paper.

function [gaborInfo,header] = getStochasticDictionaryMP3p1_parallel(data,timeVals,maxIteration,adaptiveDictionaryParam,dictionarySize)

if ~exist('maxIteration','var');           maxIteration=50;             end
if ~exist('adaptiveDictionaryParam','var');adaptiveDictionaryParam=0.9; end
if ~exist('dictionarySize','var');      dictionarySize=[];              end

Fs=round(1/(timeVals(2)-timeVals(1)));          % Sampling Frequency
numTrials = size(data,1);
sigLen = size(data,2);

%%%%%%%%%%%%%%%%%%%%%_________MP Parameters _________%%%%%%%%%%%%%%%%%%%
recAc=100; % Pecentage of Maximum reconstructed Energy (Value should be less than 100)

%%%%%%%%%%%%%%%%%%_______ Run MP on single trials ________%%%%%%%%%%%%%%
gaborInfo = zeros(numTrials,maxIteration,7);
header = zeros(numTrials,8);
% paramstruct.data = data; 
% paramstruct.dictionarySize = dictionarySize;
% paramstruct.sigLen = sigLen;
% paramstruct.maxIteration = maxIteration;
% paramstruct.recAc = recAc;
% paramstruct.Fs = Fs;
% paramstruct.adaptiveDictionaryParam = adaptiveDictionaryParam;

for i = 1:numTrials
    disp(['Trial ' num2str(i) ' of ' num2str(numTrials)]);
    iter_ID = i;

    % Save data as ASCII files
    sig=data(i,:)'; %#ok<NASGU>
    save(['sig',num2str(iter_ID),'.txt'],'sig','-ascii');
%     clear sig;
    
    % Write Command File
    fp=fopen(['commands',num2str(iter_ID),'.txt'],'w');
    if ~isempty(dictionarySize)
        fprintf(fp,['reinit -R ' num2str(dictionarySize) ' \n ']);
    end
    fprintf(fp,['set -O ' num2str(sigLen) ' -M ' num2str(maxIteration) ' -E ' num2str(recAc) ' -F ' num2str(Fs) ' -D ' num2str(adaptiveDictionaryParam) ' \n ']);
    fprintf(fp,['loadasc -O sig',num2str(iter_ID),'.txt\n ']);
    fprintf(fp,'mp\n '); % Run MP
    fprintf(fp,['save -S mpresults',num2str(iter_ID),'\n ']);
    fprintf(fp,'exit');        % Exit shell
    fclose(fp);
end

parfor i=1:numTrials
    i, numTrials
    runMPparallel(i);
end
    
parfor i=1:numTrials
    i, numTrials
    [gbinf, hd] = readMPparallel(i);
    gaborInfo(i,:,:) = gbinf;
    header(i,:,:) = hd;
end
end

% function [gaborInfo, header] = runMPparallel(iter_ID)
%     
%     disp('In');
%     iter_ID = iter_ID(1);
% %     paramstruct = i{3};
% %     i = [i{1}, i{2}]; 
%     
%     
%     % Run MP
%     system([which('mp31.exe') ' < commands',num2str(iter_ID),'.txt'],'-echo');
%     
%     % Read Data
%     disp(['Entering Readbook: Trial>',num2str(iter_ID)])
%     [gaborInfo, header]=readbook(['mpresults',num2str(iter_ID)],0);
%     
%     % Delete tmp files
%     delete(['sig',num2str(iter_ID),'.txt']);                         % Deletes the signal file
%     delete(['commands',num2str(iter_ID),'.txt']);                    % Deletes the commamnds file created
%     delete(['mpresults',num2str(iter_ID)]);                       % Deletes mpresults
% end

function runMPparallel(iter_ID)
    
    disp(['In iter #' ,num2str(iter_ID(1))]);
    iter_ID = iter_ID(1);
%     paramstruct = i{3};
%     i = [i{1}, i{2}]; 
    
    
    % Run MP
%     mp31exeloc = which('mp31.exe')
    mp31exeloc = '..\mp31.exe'
    system([mp31exeloc ' < commands',num2str(iter_ID),'.txt'],'-echo');
end

function [gaborInfo, header] = readMPparallel(iter_ID)
    % Read Data
    disp(['Entering Readbook: Trial>',num2str(iter_ID)])
    [gaborInfo, header]=readbook(['mpresults',num2str(iter_ID)],0);
    
    % Delete tmp files
    delete(['sig',num2str(iter_ID),'.txt']);                         % Deletes the signal file
    delete(['commands',num2str(iter_ID),'.txt']);                    % Deletes the commamnds file created
    delete(['mpresults',num2str(iter_ID)]);                       % Deletes mpresults
end