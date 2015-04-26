% the whole process of spike sorting
% the basic stratege is to process the data one channel by one channel
%   runstat=spikesort(fileName,varargin)
% fileName does not have postfix
% 'RunStep': default=':'; {'filt','spike detect', 'spike align',
% 'spike feature', 'cluster'}, or use 'resume'
% 'GoOn': default='on'; {'off'}
function runstat=spikesort(fileName,varargin)
%%%%%%%%%% Parameter Default
%%% 运行相关
runOpt=':';
bUserSpecifiedMem=false;

% Process control usage
% run flag - marking which of processings are needed
runFlag=false(11,1); % default: all stop
% function-number table
funcNumTab=struct('start',1,'filt',2,'filtLFP',11,'spikeDetect',4, 'spikeAlign',5, ...
    'spikeFeature',6, 'cluster',8, 'final',10); 
% 888 'spikeAlign',6, 

%%% 处理步骤的?
paras=struct();
% filtering
paras.spkFiltFreq=[200,5000];
paras.bRmc=false; % whether to remove big noise channel
paras.rmChThres=5/1000; % (V)

% LFP filtering
paras.lfpFiltFreq=[4,300];
paras.lfpSrate=1000;

% Spike detect
paras.movThres=6;
paras.bPNSep=false; % whether to separate the positive and negative peaks to separate channels

% Spike align
paras.alignWin=[-1,1.2]; % (ms)
paras.bAlignSmooth=true;

% Spike feature
paras.feaDim=3;

% Spike cluster
cluMinSpkThres=6; % cluster number minimum (/min)

% Names of files used
fnRaw=[fileName,'.mcd'];
fnF=[fileName,'_f'];
fnFL=[fileName,'_fl'];
fnS=[fileName,'_s'];
fnFEA=[fileName,'_fea'];
fnC=[fileName,'_c'];


%%%%%%%%%%%% User input 
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'RunStep'
                % see different options
                runOpt=pinfo{parai};            
            case 'direct mem'
                bUserSpecifiedMem=true;
                bDirectMem=pinfo{parai};
            otherwise
                error('unidentified options');
        end
    end
end

%%% Process of stages

if ischar(runOpt) % 只列单项
    runOpt={runOpt};
elseif iscell(runOpt) % cell中多个项目，表示要进行的步骤  <<<<
else
    error('invalid run step option');
end

sna=length(runOpt);
% 输入只有一项的情况下
if sna==1 
    if strcmp(runOpt{1},':') % run all steps
        runFlag(:)=true;
    else
        runFlag(rs2pgnum(runOpt{1},funcNumTab))=true;
    end

% 输入不止一项的情况下     
else
    for sni=1:sna
        if strcmp(runOpt{sni},':') % 包含通配符
            if sni==1
                runFlag(1:rs2pgnum(runOpt{sni+1},funcNumTab))=true;
            elseif sni==sna
                runFlag(rs2pgnum(runOpt{sni-1},funcNumTab):rs2pgnum('final',funcNumTab))=true;
            else
                runFlag(rs2pgnum(runOpt{sni-1},funcNumTab):rs2pgnum(runOpt{sni+1},funcNumTab))=true;
            end
        else
            runFlag(rs2pgnum(runOpt{sni},funcNumTab))=true;
            % * 这里标记有些重复，不过为了代码简洁，就这样了
        end
    end
end

%%% Process of others

%%% Exchange data file - exist at the location of this file
% Computer system
computerType=computer();
if strcmp(computerType(1:5),'PCWIN')
    mark='\';
else
    mark='/';
end

sspath=which([fileName,'.mcd']);
if isempty(sspath)
    sspath=[pwd(),mark];
    fprintf('did not find %s.mcd, put exchange file at: %s\n',fileName,sspath);
else
    sspath=fname(sspath);
end
sspath=[sspath,'rsdata.mat'];

%%% Determine the file reading method <<<< change to memory determined
%%% method.
if ~bUserSpecifiedMem % if not specified by user    
    switch computerType
        case {'PCWIN64','PCWIN'}
            bDirectMem=false;
        otherwise % Linux
            bDirectMem=true;
    end
end

% process fileName
[~,fileName,~]=fname(fileName);


%%% #1 Initiate the Run status data, Parameter. 这是每次运行必须的
    % Information of data    
%     if exist([fnS,'.mat'],'file') % load processed information from _s.mat file
%         load(fnS,'info');
%     else % directely from .mcd file
        info=getMCDinfo([fileName,'.mcd'],bDirectMem);
        info.chAmt=info.rawchAmt;
%     end
    
    % save the source file name that data is from
    info.fileName=fileName; 
    
%     % Check completeness of key fields
%     tp=isfield(info,{'chAmt','chID','ptsAmt','srate'});
%     if ~multilogic(tp,'and',2)
%         error('the information for processing is not complete');
%     end

    %%% save the above 2 + runFlag
    paras.bDirectMem=bDirectMem;
    save(sspath,'info','paras','runFlag');



%%%%%%%%%%% Main
% 输出数据对应文件名称
% outfile1=FnF; outfile2=FnS; outfile3=FnC;

%%%%%% #2 Filtering and removal of bad channels
if runFlag(funcNumTab.filt)   
    runstat='filt';    
    % Get filter  
    [paras.fb,paras.fa]=butter(4,paras.spkFiltFreq/(info.srate/2));    
    
    % Update run status and parameters & Execute
%     save(sspath,'-append','paras','runstat');
    
    disp('filt >>>');
    % Set up input
    [nsresult, infile1] = ns_OpenFile(fnRaw);
    if nsresult==-1,     error('open .mcd file error'); end    
    % Set up output
    if bDirectMem
        outfile1=struct('X',[],'T',[]);
    else
        outfile1=matfile(fnF,'Writable',true);
    end
    % add time points data
    outfile1.T=(1:info.ptsAmt)'/info.srate; % measured in (s)    
    
    %%% Do the filtering
%     filerun('filt',fnF,fileName);
    outfile1.X=zeros(info.ptsAmt,info.chAmt);
    rmMark=false(info.chAmt,1); % marking bad channels
    for chi=1:info.chAmt        
        [~,~,temp]=ns_GetAnalogData(infile1,info.dataChEntity(chi),1,info.ptsAmt);
%         temp=double(temp); % <<< double()?
        temp=filtfilt(paras.fb,paras.fa,temp);

        % When the noise level (measured by STD) overpass a threshold
        if paras.bRmc && std(temp)>paras.rmChThres 
            rmMark(chi)=true;
            fprintf('X');
        else
            outfile1.X(:,chi)=temp;
            fprintf('|');
        end
    end        
    fprintf('\n');    
    
    % Remove extra space of X chLabel if any
    if sum(rmMark)>0
        outfile1.X(:,rmMark)=[];
        
        % Update info
        info.badChannel=rmMark;       
        temp=info.chLabel;
        temp(rmMark)=[];
        info.chLabel=temp;
        info.rawChAmt=info.chAmt;
        info.chAmt=info.chAmt-sum(rmMark);
    end
    info.rawchLabel=info.chLabel;
    
    ns_CloseFile(infile1);
%     save(sspath,'-append','info');
    disp('filtering done'); % 结束后报喜
    %%% #S Save the filtering output data if runs in direct memory mode
    if bDirectMem
        disp('saving filtered X ...');
        outfile1.info=info;
        save(fnF,'-v7.3','-struct','outfile1');
        disp('saving done.');
    end
end


%%%%%% #4 Spike detect
if runFlag(funcNumTab.spikeDetect)
    runstat='spike detect';    
    % Update run status and parameters
%     save(sspath,'-append','runstat');    
    
    % Set up input - if no filtered data (outfile1) in memory, try load it. In this case,
    % it is input.
    if ~exist('outfile1','var')
        fprintf('loading raw signal from %s\n',fnF);
        if bDirectMem
            outfile1=load(fnF);
        else
            outfile1=matfile(fnF);
        end
    end
    % Set up output
    if bDirectMem
        outfile2=struct();
    else
        outfile2=matfile(fnS,'Writable',true);
    end
    
    %%% Do the detection
    disp('spike detect >>>');
    SD=cell(info.chAmt,1);
    SQ=cell(info.chAmt,1);
    SA=cell(info.chAmt,1);
    for chi=1:info.rawchAmt
        [SD{chi},SQ{chi},SA{chi}]=spike_detect(outfile1.X(:,chi),info.srate,'movThres',paras.movThres);        
        fprintf('|');
    end
    fprintf('\n');    
    outfile2.SD=SD; outfile2.SQ=SQ; outfile2.SA=SA;
    outfile2.info=info;
    
    if bDirectMem
        save(fnS,'-v7.3','-struct','outfile2');        
    end
    disp('detection done');
    
    
    %%% Splite the positive and negative spikes 
    if paras.bPNSep
        disp('split N/P spikes >>>');
        newSD=cell(0,1);
        newSA=cell(0,1);
        newSQ=cell(0,1);
        newchID=zeros(0,1);
        newchLabel=cell(0,1);

        chcount=0;
        for chi=1:info.rawchAmt
            PN=SQ{chi}(:,1); % get posi/nega mark
            pnum=sum(PN); % number of positive peaks
            if pnum>0
                chcount=chcount+1;

                % assign to the newSD
                newSD{chcount,1}=SD{chi}(PN);
                % assign to new SQ
                newSQ{chcount,1}=SQ{chi}(PN,:);
                % assign to new sA
                newSA{chcount,1}=SA{chi}(PN,:);
                % assign new chID
                newchID(chcount,1)=chi;
                % assign new chLabel
                newchLabel{chcount,1}=info.rawchLabel{chi};
            end

            if pnum<length(SD{chi}) % number of negative peaks > 0
                chcount=chcount+1;

                PN=~PN;
                % assign to the newSD
                newSD{chcount,1}=SD{chi}(PN);
                % assign to new SQ
                newSQ{chcount,1}=SQ{chi}(PN,:);
                % assign to new sA
                newSA{chcount,1}=SA{chi}(PN,:);
                % assign new chID
                newchID(chcount,1)=chi;
                % assign new chLabel
                newchLabel{chcount,1}=info.rawchLabel{chi};
            end
        end
        
        SD=newSD;
        outfile2.SD=newSD;
        outfile2.SA=newSA;
        outfile2.SQ=newSQ;
        outfile2.info.chID=newchID;
        outfile2.info.rawchLabel=info.chLabel;
        outfile2.info.chLabel=newchLabel;
        outfile2.info.chAmt=chcount;
        
        if bDirectMem
            save(fnS,'-v7.3','-struct','outfile2');        
        end
        disp('splited');
        clear newSD newSQ newSA newchID newchLabel
    end
    
end


%%%%%%% Spike alignment is separated here.
if runFlag(funcNumTab.spikeAlign)    
    %%% load
    if ~exist('outfile1','var')
        fprintf('loading raw signal from %s\n',fnF);
        if bDirectMem
            outfile1=load(fnF);
        else
            outfile1=matfile(fnF);
        end
    end
    if ~exist('outfile2','var')
        if bDirectMem
            outfile2=load(fnS);
        else
            outfile2=matfile(fnS);
        end        
    end    
    
    disp('spike alignment >>>');
    if bDirectMem
        [A,rmlist]=spike_align(outfile1.X,outfile2.SD,info.srate,'chAssign',outfile2.info.chID,'window',paras.alignWin,'bSmooth',paras.bAlignSmooth);
    else % 由于matfile方法处理2D矩阵的限制？
        [A,rmlist]=spike_align(outfile1,outfile2.SD,info.srate,'chAssign',outfile2.info.chID,'window',paras.alignWin,'bSmooth',paras.bAlignSmooth);
    end
    if ~isempty(rmlist)
        SD=outfile2.SD;
        % * must do it in reversed order
        for k=size(rmlist,1):-1:1
            SD{rmlist(k,1)}(rmlist(k,2))=[];
        end
        if bDirectMem
            save(fnS,'-v7.3','-append','SD');
        else
            outfile2.SD=SD;
        end
    end
    
    if bDirectMem
        save(fnS,'-v7.3','-append','A');
    end
    disp('data alignment done');
end


%%%%%% #3 Filtering for LFP
if runFlag(funcNumTab.filtLFP)
    runstat='filtLFP';
    % Get filter
    [paras.lfblp,paras.lfalp]=butter(4,paras.lfpFiltFreq(2)/(info.srate/2),'low');
    [paras.lfbhp,paras.lfahp]=butter(4,paras.lfpFiltFreq(1)/(info.srate/2),'high');
    
    paras.lfpDSratio=paras.lfpSrate/info.srate;
    
    save(sspath,'-append','paras','runstat');
    
    filerun('filtLFP',fnFL,fileName);
    disp('LFP filtering done');
end


%%%%%% #8 Feature extract
if runFlag(funcNumTab.spikeFeature)
    % Set input
    if ~exist('A','var')
        fprintf('loading aligned spike data from %s\n',fnS);
        load(fnS,'A');
    end
    
%     filerun('spike feature',fnFEA,fnS);
    disp('spike feature >>>');
    SF=cell(info.chAmt,1);
    for chi=1:info.chAmt
        if isempty(A{chi,1})
            SF{chi,1}=[];
            fprintf('X');
        else
            SF{chi,1}=spike_feature(A{chi,1},'dim',paras.feaDim);
            fprintf('|');
        end        
    end
    fprintf('\n');
    
    % Save output
    save(fnC,'-v7.3','SF');
    disp('feature extract done');
end


%%%%%% #11 Clustering    
if runFlag(funcNumTab.cluster)    
%     filerun('cluster',fnC,fnFEA,fnS,T);
    disp('spike clustering >>>');
    
    % Set input
    load(fnS,'SD');
    if ~exist('SF','var'), load(fnC,'SF'); end
    if ~exist('A','var'), load(fnS,'A'); end
    
    %%% Set output
    if exist('outfile1','var')
        time=[outfile1.T(1,1), outfile1.T(info.ptsAmt,1)];
    else
        load(fnF,'T');
        time=[T(1,1), T(info.ptsAmt,1)];
    end
    totalTime=(time(2)-time(1)); % total time span of data        
    % cluster number minimum is based on spikes/min, so the actual number
    % should be calculated by time length
    cluMinSpkThres=cluMinSpkThres*(totalTime/60);
    
    % Channel SI: cluster identity of each spike in raw channel.
    CSI=cell(info.chAmt,1); % * note there may be 0 in CSI - those not assigned to any cluster.       
    % newchID: from which raw channel (electrode) the neuron is recorded.
    newchID=zeros(0,1);
    % NSD is neuron spk data (index format)
    newSD=cell(0,1);

    %%% Clutering each channel. 
    chcount=0; % number of effective neurons - 数量不能预知，所以累计计数。
    for chi=1:info.chAmt
        tp=length(SD{chi});
        if tp<=cluMinSpkThres % spike数量过少的通道直接复制原值（不删除为好）。
            chcount=chcount+1;
            newSD{chcount,1}=SD{chi};            
            newchID(chcount,1)=info.chID(chi);
            CSI{chi}=zeros(tp,1);
        else
            [CSI{chi},lbinfo]=spike_cluster(SF{chi,1},'kmeans',cluMinSpkThres,A{chi,1},info.srate);
            if lbinfo.cAmt>0 % * it is possible that non is left after clustering. 
                % 聚类结果的各个cluster依次添加
                for k=1:lbinfo.cAmt
                    chcount=chcount+1;
                    % sort the index in ascending order, since lbinfo.ids{} does not guarantee monotonic increasing.
                    temp=sort(lbinfo.ids{k},'ascend'); 
                    newSD{chcount,1}=SD{chi}(temp);
                    newchID(chcount,1)=info.chID(chi);
                end
            end
        end
        fprintf('|');
    end
    fprintf('\n');
    
    
    %%% 还原的原始的通道数量下的格式
    % * 此步骤的作用主要在于可视化分类结果：因为在detect and clutering之间的
    % 步骤可能发生过通道数量的变化和spike的删减，之后的spkalignplot()要正确显示
    % 很麻烦，所以生成此数据。 *
    % CSI to raw channels
    newCSI=cell(info.rawchAmt,1);
    rawSD=cell(info.rawchAmt,1);
    for k=1:info.chAmt
        % find raw channel index of the neuron
        chi=info.chID(k);
        % Insert to new CSI
        if isempty(newCSI{chi})
            newCSI{chi}=CSI{k};
            rawSD{chi}=SD{k};
        else
            % find biggest number of existing neuron id
            nn=max(newCSI{chi});
            % except 0, add the nn as base
            I=(CSI{k}~=0);
            tp=zeros(length(CSI{k}),1); tp(I)=nn+CSI{k}(I);
            newCSI{chi}=[newCSI{chi};tp];
            rawSD{chi}=[rawSD{chi};SD{k}];
        end
    end
    % sort to time order (different cluster mixed in order)
    for k=1:info.rawchAmt
        if ~isempty(newCSI{k})
            [rawSD{k},I]=sort(rawSD{k},'ascend');
            newCSI{k}=newCSI{k}(I);
        end
    end
%     CSI=newCSI;
    
    %%% Update channel ID & chAmt
    info.chID=newchID;
    info.chAmt=chcount;
    SD=newSD;
    
    % Save output
    save(fnC,'-v7.3','-append','time','rawSD','SD','CSI','newCSI','info');
    disp('clustering done');
end


%%%%%% #13 cleaning up
delete(sspath);
disp('run finish, all done');
end


%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
% To identify the runFlag of run step option in string form
function N=rs2pgnum(rsstr,funcNumTab)
switch rsstr
    case 'start'
        N=funcNumTab.start;
    case 'filt'
        N=funcNumTab.filt;
    case 'filtLFP'
        N=funcNumTab.filtLFP;
    case 'spike detect'                
        N=funcNumTab.spikeDetect;
    case 'spike align'
        N=funcNumTab.spikeAlign;
    case 'spike feature'
        N=funcNumTab.spikeFeature;
    case 'cluster'
        N=funcNumTab.cluster;
    case 'final'
        N=funcNumTab.final;
    otherwise
        error('invalid run option');
end

end