% the whole process of spike sorting
% the basic stratege is to process the data one channel by one channel
%   runstat=spikesort(fileName,varargin)
% 'RunStart': default='start'; {'filt','spike detect', 'spike align',
% 'spike feature', 'cluster'}, or use 'resume'
% 'GoOn': default='on'; {'off'}
function runstat=spikesort(fileName,varargin)
%%%%%%%%%% Init 
% Parameter Default
runOpt='start';
bGoOn=true;
bServer=false;
bResume=false;

% run flag - marking which of processings are needed
pgNum=false(14,1); % default: all stop
% function-number table
funcNumTab=struct('start',1,'filt',2, 'spikeDetect',4, 'spikeFeature',8, 'cluster',11, 'final',13); 
% 888 'spikeAlign',6, 

%%% Option input or struct of data
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'RunStart'
                % see different options
                runOpt=pinfo{parai};
            case 'GoOn'
                if strcmp(pinfo{parai},'off')
                    bGoOn=false;
                else
                    if ~strcmp(pinfo{parai},'on')
                        error('unknown option for ''go_on'' ');
                    end
                end
            case 'Server'
                if pinfo{parai}
                    bServer=true;
                end
            otherwise
                error('unidentified options');
        end
    end
end

%%% Analyze the options
if ischar(runOpt) % 只列单项，说明只需要进行这一步骤
    
elseif iscell(runOpt) % cell中多个项目，表示要炮的步骤
    % ':'代表之前/之后（由符合位置决定）都需要进行。
    
end

%%% Exchange data file - Find at the location of this function ('spikesort')
sspath=which('spikesort');
sspath=fname(sspath); 
sspath=[sspath,'rsdata.mat'];

%%% Determine the procedure used by progress number "pgNum"
switch runOpt
    case 'start' % "start" means begin from start, all to be run
        if bGoOn
            pgNum(:)=true;
        else
            pgNum(funcNumTab.start)=true;
        end        
    % The specified names indicate to run from that step just that  
    case 'filt'
        if bGoOn
            pgNum(funcNumTab.filt:end)=true;
        else
            pgNum(funcNumTab.filt)=true;
        end
    case 'spike detect'
        if bGoOn
            pgNum(funcNumTab.spikeDetect:end)=true;
        else
            pgNum(funcNumTab.spikeDetect)=true;
        end
    case 'spike feature'
        if bGoOn
            pgNum(funcNumTab.spikeFeature:end)=true;
        else
            pgNum(funcNumTab.spikeFeature)=true;
        end
    case 'cluster'
        if bGoOn
            pgNum(funcNumTab.cluster:end)=true;
        else
            pgNum(funcNumTab.cluster)=true;
        end
    case 'final'
        if bGoOn
            pgNum(funcNumTab.final:end)=true;
        else
            pgNum(funcNumTab.final)=true;
        end
    case 'resume'
        bResume=true;
        % Load the exchange data file        
        load(sspath);
        if ~strcmp(info.fileName,fileName)
            disp('warning: the loaded runstat file may not be the right one for the current file');
        end
        if exist('runstat','var')
            switch runstat
                % In resume mode, the stored string indicate that procedure
                % is done, so start from next by "+1" in pgNum
                case 'start'
                    pgNum(:)=true;
                case 'filt'
                    pgNum(funcNumTab.filt+1:end)=true;
                case 'spike detect'
                    pgNum(funcNumTab.spikeDetect+1:end)=true;
                case 'spike feature'
                    pgNum(funcNumTab.spikeFeature+1:end)=true;
                case 'cluster'
                    pgNum(funcNumTab.cluster+1:end)=true;
                case 'final'
                    pgNum(funcNumTab.final+1:end)=true;
                otherwise
                    error('unrecognized run status from loaded file');
            end
        else
            error('can not find the run status variable');
        end
        
    otherwise
        error('invalid run option');
end


%%%%%%%%%%% Main
%%%%%% #1 Initiate the Run status data, Parameter. 
if pgNum(funcNumTab.start)
    %%% Information of data
    if bServer    
        load(fileName,'info'); % now we get "info" variable
        load(fileName,'srate'); % now we get "srate" variable
        info.srate=srate;
    else
        datafile=matfile([fileName,'.mat']);
        info=datafile.info;
        info.srate=datafile.srate;
    end
    info.fileName=fileName;
    
    % check completeness of fields
    tp=isfield(info,{'chAmt','chID','ptsAmt','srate'});
    if ~multilogic(tp,'and',2)
        error('the information for processing is not complete');
    end

    %%% Also init the parameters
    paras=struct('bDirectMem',bServer);

    %%% Also the run state
    if ~bResume
        save(sspath,'info','paras');
        % status need to consider the pgNum
        if pgNum(1)    
            runstat='start';
            save(sspath,'-append','runstat');
        end
    end
end

%%%%%% #2 Filtering and removal of bad channels
fnF=[fileName,'_f'];
% load(sspath);
if pgNum(funcNumTab.filt)    
%     [paras.fb,paras.fa]=butter(4,10000/(srate/2),'low');
    [paras.fb,paras.fa]=butter(4,[200,10000]/(info.srate/2));
    paras.bRmc=true; % whether to remove big noise channel
    paras.rmChThres=50; % (mV)
    save(sspath,'-append','paras');
    
    filerun('filt',fnF,fileName);
    
    % Update run status
    runstat='filt';
    save(sspath,'-append','runstat');
    disp('filtering done');
end

%%%%%% #4 Spike detect
fnS=[fileName,'_s'];
if pgNum(funcNumTab.spikeDetect)
    load(sspath,'paras');
    paras.bNegaSep=false; % whether to separate the positive and negative peaks to separate channels
    save(sspath,'-append','paras');
    
    filerun('spike detect',fnS,fnF);
    
    % Update run status
    runstat='detect';
    save(sspath,'-append','runstat');
    disp('detection done');
end

%%%%%% #8 Feature extract
fnFEA=[fileName,'_fea'];
if pgNum(funcNumTab.spikeFeature)
    filerun('spike feature',fnFEA,fnS);
    
    % Update run status
    runstat='spike feature';
    save(sspath,'-append','runstat');
    disp('feature extract done');
end

%%%%%% #11 Clustering    
fnC=[fileName,'_c'];
if pgNum(funcNumTab.cluster)
    filerun('cluster',fnC,fnFEA,fnS);
    
    % Update run status
    runstat='cluster';
    save(sspath,'-append','runstat');
    disp('clustering done');
end

%%%%%% #13 cleaning up
if pgNum(funcNumTab.final)
    delete(sspath);
end

disp('run finish, all done');