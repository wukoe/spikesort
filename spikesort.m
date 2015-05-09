% the whole process of spike sorting
% the basic stratege is to process the data one channel by one channel
%   runstat=spikesort(fileName,varargin)
% fileName does not have postfix
% 'RunStep': default=':'; {'filt','spike detect', 'spike align',
% 'spike feature', 'cluster'}, or use 'resume'
% 'GoOn': default='on'; {'off'}
function spikesort(fileIn,varargin)
%%%%%%%%%% Parameter Default
%%% 运行相关
% bFullProc=true;
runOpt=':';
bUserSpecifiedMem=false;
bNewOut=false;

% Process control usage
% run flag - marking which of processings are needed
runFlag=false(11,1); % default: all stop
% function-number table
funcNumTab=struct('start',1,'filt',2,'spikeDetect',4,'merge',5,'spikeAlign',6, ...
     'cluster',8,'final',10); % 'spikeFeature',6,

%%% 处理步骤的?
paras=struct();
% filtering
paras.spkFiltFreq=[200,5000];
paras.bRmc=false; % whether to remove big noise channel
paras.rmChThres=5/1000; % (V)
paras.bSaveFiltSignal=true;

% Spike detect
paras.movThres=6.5;
paras.bNP=true;
paras.bPP=true;

% Spike align
paras.alignWin=[-0.6,0.7]; % [-0.8 1.1] (ms)
paras.bAlignSmooth=false;
paras.bExactST=false;

% Spike cluster
paras.bPNSep=true; % whether to separate the positive and negative peaks to separate channels first.
% *this option is related to feature selection and clustering.
% Test found that it's bad to set it as false!!
paras.bSpkShapeOL=true; % whether to filter out spike shape outlier.
paras.bSpkFeatureOL=true; % whether to filter out spike feature outlier.
paras.cluMinSpkThres=30; % cluster number minimum (用绝对数量而非spiking rate来定义)
paras.feaDim=3;% Spike feature
paras.bDrawForCluster=true;
paras.drawNum=2000;
paras.cluMethod='kmeans';
paras.bMatchNoisySpk=true;
paras.SNratioThres_noise=4;
paras.bTimeShiftMerge=true;
paras.lagWin=[-0.5,0.5]; % shift merging's maximum range (ms)
paras.ccThres=0.8;%0.91; %<<< cluster mean curve correlation threshold 0.85
paras.SNratioThres_cm=12; % 15 was tested.
paras.bIndivMove=false; % move each spike with individual distance or not.
%
paras.STactThre=0.2; % (spikes/sec)

%%%%%%%%%%%% User input 
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'RunStep'
                % see different options
                runOpt=pinfo{parai};
%                 bFullProc=false;
            case 'direct mem'
                bUserSpecifiedMem=true;
                bDirectMem=pinfo{parai};
            case 'new output'
                bNewOut=true;
                fileOut=pinfo{parai};
            otherwise
                error('unidentified options');
        end
    end
end

%%% Process of stages
if ischar(runOpt) % 只列单项
    runOpt={runOpt};
elseif iscell(runOpt) % cell中多个项目，表示要进行的步骤 
else
    error('invalid run step option');
end

sna=length(runOpt);
% 输入只有一项的情况下
if sna==1 
    if strcmp(runOpt{1},':') % run all steps
        runFlag(:)=true;
    else
        runFlag(rs2pgnum(runOpt{1}))=true;
    end

% 输入不止一项的情况下     
else
    for sni=1:sna
        if strcmp(runOpt{sni},':') % 包含通配符
            if sni==1
                runFlag(1:rs2pgnum(runOpt{sni+1}))=true;
            elseif sni==sna
                runFlag(rs2pgnum(runOpt{sni-1}):rs2pgnum('final'))=true;
            else
                runFlag(rs2pgnum(runOpt{sni-1}):rs2pgnum(runOpt{sni+1}))=true;
            end
        else
            runFlag(rs2pgnum(runOpt{sni}))=true;
            % * 这里标记有些重复，不过为了代码简洁，就这样了
        end
    end
end

% Computer system
computerType=computer();
% Determine the file reading method <<<< change to memory determined
% method.
if ~bUserSpecifiedMem % if not specified by user    
    switch computerType
        case {'PCWIN64','PCWIN'}
            bDirectMem=false;
        otherwise % Linux
            bDirectMem=true;
    end
end

% In&Out fileNames
[inPath,fileIn,~]=fname(fileIn);
fileIn=[inPath,fileIn];
if bNewOut
    [outPath,fileOut,~]=fname(fileOut);
    fileOut=[outPath,fileOut];
else
    fileOut=[inPath,fileIn];
end

%%% Names of files used
fnRaw=[fileIn,'.mcd'];
fnF=[fileOut,'_f'];
% fnFL=[fileName,'_fl'];
fnS=[fileOut,'_s'];
% fnFEA=[fileName,'_fea'];
fnC=[fileOut,'_c'];
fnST=[fileOut,'_st'];
disp(fnRaw);


%%% #1 Initiate the Run status data, Parameter. 这是每次运行必须的
% Get info. struct data.
% * 如果要运行filt或spike detect步骤，则从原始文件中获得信息；如果是其他步骤，则查看fnS文件是否存在，若是，则从中加载现成的info
if ~runFlag(funcNumTab.filt) && ~runFlag(funcNumTab.spikeDetect) && exist([fnS,'.mat'],'file') % load processed information from _s.mat file
    load(fnS,'info');
else % directely from .mcd file
    info=getMCDinfo(fnRaw,bDirectMem);    
    % 删除多余的事件通道 <<< 目前没有更好策略的情况下暂定。
    if info.rawchAmt>120
        info.rawchAmt=120;
        disp('extra electrode removed:');
        disp(info.elecLabel(121:end)');
        info.elecLabel=info.elecLabel(1:120);
        info.chID=info.chID(1:120);
    end
    
    info.chAmt=info.rawchAmt;
end
paras.bDirectMem=bDirectMem;

% save the source file name that data is from
info.fileIn=fnRaw; 

%     % Check completeness of key fields
%     tp=isfield(info,{'chAmt','chID','ptsAmt','srate'});
%     if ~multilogic(tp,'and',2)
%         error('the information for processing is not complete');
%     end

%%% Process of others
paras.lagWin=floor(paras.lagWin/1000*info.srate);


%%%%%%%%%%%% Main
% 输出数据对应文件名称
% outfile1=FnF; outfile2=FnS; outfile3=FnC;
%%%%%% #2 Filtering and removal of bad channels
if runFlag(funcNumTab.filt)    
    % Set up input
    [nsresult, infile1] = ns_OpenFile(fnRaw);
    if nsresult==-1,     error('open .mcd file error'); end
    
    disp('filt >>>');
    % Get filter
    [paras.fb,paras.fa]=butter(4,paras.spkFiltFreq/(info.srate/2));
    
    % add time points data
    FtT=(1:info.ptsAmt)'/info.srate; % measured in (s)
    
    %%% Do the filtering
    FtX=zeros(info.ptsAmt,info.chAmt);
    rmMark=false(info.chAmt,1); % marking bad channels
    for chi=1:info.chAmt  %%<<< parfor注意
        [~,~,temp]=ns_GetAnalogData(infile1,info.dataChEntity(chi),1,info.ptsAmt);
        temp=filtfilt(paras.fb,paras.fa,temp);
        
        % When the noise level (measured by STD) overpass a threshold
        if paras.bRmc && std(temp)>paras.rmChThres
            rmMark(chi)=true;
            fprintf('X');
        else
            FtX(:,chi)=temp;
            fprintf('|');
        end
    end
    fprintf('\n');
    
    % Remove noisy channel of X if any
    if sum(rmMark)>0
        % FtX(:,rmMark)=[];        
        % Update info
        info.badChannel=rmMark;
        % info.chAmt=info.chAmt-sum(rmMark);
    end
    
    ns_CloseFile(infile1); 
    disp('filtering done'); % 结束后报喜
    
    %%% #S Save the filtering output data if runs in direct memory mode    
    if paras.bSaveFiltSignal && bDirectMem
        disp('saving filtered X ...');        
        save(fnF,'-v7.3','FtX','FtT','info'); %
        disp('saving done.');
    end
end


%%%%% #4 Spike detect
if runFlag(funcNumTab.spikeDetect)
    % Set up input - if no filtered data (outfile1) in memory, try load it.
    loadoutfile1();
    load(fnF,'info'); %info=outfile1.info;
    % Set up output
    outfile2=struct();    
    
    disp('spike detect >>>');    
    rawSD=cell(info.rawchAmt,1);
    rawSQ=rawSD; rawSA=rawSD; rawSW=rawSD;
    parfor chi=1:info.rawchAmt
        [rawSD{chi},rawSQ{chi},rawSA{chi},rawSW{chi}]=spike_detect(FtX(:,chi),info.srate,paras);
        fprintf('|');
    end
    fprintf('\n');
    
    % Save raw SD
    outfile2.rawSD=rawSD; outfile2.rawSQ=rawSQ; outfile2.rawSA=rawSA; outfile2.rawSW=rawSW;    
    % Save SD and information
    info.chID=(1:info.rawchAmt)';
    outfile2.info=info;  
    outfile2.paras=paras;
    % Write to disk
    save(fnS,'-struct','outfile2'); %,'-v7.3'
    disp('detection done');
end


%%%%% Spike merging 
% *!!! (Merging is put here only because it uses biggest amplitude as
% standard to keep spikes in a cluster. -this can do the merging when
% unsorted data is used.)
if runFlag(funcNumTab.merge)
    if ~exist('outfile2','var')
        outfile2=load(fnS);
    end    
    
    disp('spike merge >>>');
    [outfile2.SD,rmlist]=spikemerge(outfile2.rawSD,info,outfile2.rawSA);
    % delete extra SA too
    outfile2.SA=outfile2.rawSA;
    for chi=1:outfile2.info.rawchAmt
        I=rmlist{chi};
        outfile2.SA{chi}(I)=[];
    end
    % Save
    save(fnS,'-struct','outfile2');
    disp('merging done');
    
    %%% Save the ST form of spike_detect too. 
    loadoutfile1();
    if paras.bExactST % Get Precise time by splining spike peaks.        
        SDT=exactST(FtX,outfile2.SD,FtT,info.srate,info.chID);
    else
        SDT=idx2time(outfile2.SD,FtT);
    end
    SDTchID=info.chID;
    save(fnST,'SDT','SDTchID');
    disp('detected SD converted to ST');
end


%%%%% Spike alignment is separated here.
if runFlag(funcNumTab.spikeAlign)    
    % Set up input
    loadoutfile1();
    if ~exist('outfile2','var')
        outfile2=load(fnS);
    end    
    if ~isfield(outfile2,'SD')
        outfile2.SD=outfile2.rawSD;
        outfile2.SA=outfile2.rawSA;
    end
    
    disp('spike alignment >>>');
    %%% Remove low activity channels and save results as SD
    sa=cellstat(outfile2.SD,'length');
    I=(sa>=5);
    SD=outfile2.SD(I); SA=outfile2.SA(I); % * SQ not need to be saved for further use.
    info=outfile2.info;
    info.chAmt=sum(I);
    info.chID=info.chID(I);
    
    %%% Spike alignment
    [A,rmlist,SO]=spike_align(FtX,SD,info.srate,'chAssign',info.chID,'window',paras.alignWin,'bSmooth',paras.bAlignSmooth);
    % 消除窗口无法覆盖的边缘spike，同步更新SD,SA,SQ等信息。
    if ~isempty(rmlist)
        for k=size(rmlist,1):-1:1 % * must do it in reversed order
            SD{rmlist(k,1)}(rmlist(k,2))=[];
            SA{rmlist(k,1)}(rmlist(k,2))=[];
        end
        save(fnS,'-v7.3','-append','SD','SA');
    end
    % Update with new data
    outfile2.A=A;
    outfile2.SD=SD; outfile2.SA=SA;
    % Get spike morphology data length (of each sample).
    info.spklen=SO.spklen; info.spkprew=SO.preww; info.spkpostw=SO.postww;
    outfile2.info=info;
    
    %%% Save results
    save(fnS,'-v7.3','-struct','outfile2');    
    % a final check
    assert(isequal(cellstat(SD,'length'),cellstat(A,'size',2)),'SD and A number not equal');
    disp('data alignment done');
end


%%%%% Spike clustering
if runFlag(funcNumTab.cluster)
    % Set up input
    loadoutfile1();
    if ~exist('outfile2','var')
        outfile2=load(fnS);        
    end
    info=outfile2.info;
    
    disp('spike cluster >>>');
    CST=cell(info.rawchAmt,1);
    reconSD=CST;    reconA=CST;
    for k=1:info.chAmt
        fprintf('|');
        chi=info.chID(k);
        [CST{chi},reconSD{chi},reconA{chi}]=chcluster(FtX(:,chi),outfile2.SD{k},outfile2.A{k},info,paras); % outfile2.SA{k},
    end
    fprintf('\n');
    
    %%% 生成新的NSD（每个分开的"神经元"算一个通道）and new chID and chAmt.
    assert(isequal(cellstat(reconSD,'length'),cellstat(CST,'length')),'reconSD and CSI length not match after clustering');
%     NSD=cell(0,1);
%     chID=zeros(0,1);
%     chcount=0;
%     for chi=1:info.rawchAmt
%         spkclu=reabylb(CST{chi});
%         % only choose non-noise channel
%         nzI=find(spkclu.types>0);  na=length(nzI);
%         if na>0
%             temp=cell(na,1);
%             for k=1:length(nzI)
%                 I=spkclu.ids{nzI(k)};
%                 temp{k}=reconSD{chi}(I);
%             end
%             NSD(chcount+1:chcount+na,1)=temp;
%             chID(chcount+1:chcount+na,1)=chi;
%             chcount=chcount+na;
%         end
%     end
%     info.chID=chID;
%     info.chAmt=chcount;
    [NSD,info.chID]=getNSD(reconSD,CST,info);
    info.chAmt=length(info.chID);

%     % Add amplitude of NSD (NSA), for use in spike merge.
%     NSA=cell(info.chAmt,1);
%     for chi=1:info.chAmt
%         NSA{chi}=FtX(NSD{chi},info.chID(chi));
%     end
    
    % save
    save(fnC,'CST','NSD','info','reconSD','reconA','paras');%,'NSA');
    disp('clustering done');
    
    %%% Save the ST form of spike_detect too.
    % First delete low-activity channels
    temp=cellstat(NSD,'length');
    thr=paras.STactThre*info.TimeSpan;
    I=(temp>thr);
    NSD=NSD(I);
    info.chID=info.chID(I);
    info.chAmt=length(info.chID);
    
    % Remaining converted to ST
    if paras.bExactST % Get Precise time by splining spike peaks.
        ST=exactST(FtX,NSD,FtT,info.srate,info.chID);
    else
        ST=idx2time(NSD,FtT);
    end
    STchID=info.chID;
    save(fnST,'-append','ST','STchID','info');
    disp('sorted SD converted to ST');
end

% Cleaning up
% delete(sspath);
disp('run finish, all done');
%%%%%%%%%%%% End of main.


%%%%%%%%%%%%%%% Load variable sub-functions
    function loadoutfile1()
        if ~exist('FtX','var')
            fprintf('loading raw signal from %s\n',fnF);
            outfile1=load(fnF);
            FtX=outfile1.FtX;
            FtT=outfile1.FtT;
        end
    end


%%%%%%%%%%%%%%% Assisting
% To identify the runFlag of run step option in string form
    function N=rs2pgnum(rsstr)
        switch rsstr
            case 'start'
                N=funcNumTab.start;
            case 'filt'
                N=funcNumTab.filt;
            case 'spike detect'
                N=funcNumTab.spikeDetect;
            case 'spike align'
                N=funcNumTab.spikeAlign;
            case 'spike cluster'
                N=funcNumTab.cluster;
            case 'spike merge'
                N=funcNumTab.merge;
            case 'final'
                N=funcNumTab.final;
            otherwise
                error('invalid run option');
        end        
    end

end