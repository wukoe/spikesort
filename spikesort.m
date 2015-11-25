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
funcNumTab=struct('start',1,'filt',2,'detect',4,'align',6, ...
     'cluster',9,'merge',8,'toST',10,'final',12); % 'spikeFeature',6,

%%% 处理步骤的?
paras=struct();
% filtering
paras.spkFiltFreq=[200,5000];
paras.bRmc=false; % whether to remove big noise channel
paras.rmChThres=5/1000; % (V)
paras.bSaveFiltSignal=true;

% Spike detect
paras.movThres=6.5; %<<<
paras.bNP=true;
paras.bPP=true;

% Spike align
paras.alignWin=[-0.6,0.7]; % [-0.8 1.1] (ms)
paras.bAlignSmooth=false;
paras.bExactST=true;

% Spike merge
paras.bClusterBeforeMerge=true;

% Spike cluster
paras.bPNSep=true; % whether to separate the positive and negative peaks to separate channels first.
% *this option is related to feature selection and clustering.
% Test found that it's bad to set it as false!!
paras.bSpkShapeOL=false; % whether to filter out spike shape outlier.
paras.bSpkFeatureOL=true; % whether to filter out spike feature outlier.
paras.cluMinSpkThres=30; % cluster number minimum (用绝对数量而非spiking rate来定义)
paras.feaDim=3;% Spike feature
paras.bDrawForCluster=true;
paras.drawNum=1200;
paras.cluMethod='kmeans';
paras.bFineCluster=true;
paras.bMatchNoisySpk=true;
paras.bAlwaysMatch=false; % whether always match noise to at least one template (not allow CSI==0)
paras.SNratioThres_noise=4;
paras.bTimeShiftMerge=true;
paras.lagWin=[-0.5,0.5]; % shift merging's maximum range (ms)
paras.ccThres=0.8;%0.91; %<<< cluster mean curve correlation threshold 0.85
paras.SNratioThres_cm=12; % 15 was tested.
paras.bIndivMove=false; % move each spike with individual distance or not.

% To ST
paras.STactThre=0.2; % (spikes/sec)


%%%%%%%%%%%% User input 
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'run step'
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
    error('run-step option must be string or cell of strings.');
end

sna=length(runOpt);
% 输入只有一项的情况下
if sna==1
    if strcmp(runOpt{1},':') % run all steps
        runFlag(:)=true;
    else
        runFlag(funcNumTab.(runOpt{1}))=true;
    end

% 输入不止一项的情况下
else
    for sni=1:sna
        if strcmp(runOpt{sni},':') % 包含通配符
            if sni==1
                runFlag(1:funcNumTab.(runOpt{sni+1}))=true;
            elseif sni==sna
                runFlag(funcNumTab.(runOpt{sni-1}):funcNumTab.final)=true;
            else
                runFlag(funcNumTab.(runOpt{sni-1}):funcNumTab.(runOpt{sni+1}))=true;
            end
        else
            runFlag(funcNumTab.(runOpt{sni}))=true;
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
if ~runFlag(funcNumTab.filt) && exist([fnF,'.mat'],'file') % load processed information from _f.mat file
    load(fnF,'info');
elseif ~runFlag(funcNumTab.detect) && exist([fnS,'.mat'],'file') % load processed information from _s.mat file
    load(fnS,'info');
else % directely from .mcd file
    info=getMCDinfo(fnRaw,bDirectMem);    
    % 删除多余的事件通道 <<< 目前没有更好策略的情况下暂定。
    if info.rawchAmt>120
        info.rawchAmt=120;
        disp('extra electrode removed:');
        disp(info.elecLabel(121:end)');
        info.elecLabel=info.elecLabel(1:120);
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
    T=(1:info.ptsAmt)'/info.srate; % measured in (s)
    
    %%% Do the filtering
    X=zeros(info.ptsAmt,info.chAmt);
    rmMark=false(info.chAmt,1); % marking bad channels
    for chi=1:info.chAmt  %%<<< parfor注意
        [~,~,temp]=ns_GetAnalogData(infile1,info.dataChEntity(chi),1,info.ptsAmt);
        temp=filtfilt(paras.fb,paras.fa,temp);
        
        % When the noise level (measured by STD) overpass a threshold
        if paras.bRmc && std(temp)>paras.rmChThres
            rmMark(chi)=true;
            fprintf('X');
        else
            X(:,chi)=temp;
            fprintf('|');
        end
    end
    fprintf('\n');
    
    % Remove noisy channel of X if any
    if sum(rmMark)>0
        % X(:,rmMark)=[];        
        % Update info
        info.badChannel=rmMark;
        % info.chAmt=info.chAmt-sum(rmMark);
    end
    
    ns_CloseFile(infile1); 
    disp('filtering done'); % 结束后报喜
    
    %%% #S Save the filtering output data if runs in direct memory mode    
    if paras.bSaveFiltSignal && bDirectMem
        disp('saving filtered X ...');        
        save(fnF,'-v7.3','X','T','info'); %
        disp('saving done.');
    end
end


%%%%% #4 Spike detect
if runFlag(funcNumTab.detect)
    % Set up input - if no filtered data (outfile1) in memory, try load it.
    if ~exist('X','var')
        fprintf('loading raw signal from %s\n',fnF);
        load(fnF,'X');
        load(fnF,'T');
    end
    load(fnF,'info'); %info=outfile1.info;
    % Set up output
    outfile2=struct();    
    
    disp('spike detect >>>');    
    rawSD=cell(info.rawchAmt,1);
    rawSQ=rawSD; rawSA=rawSD; rawSW=rawSD;
    sigSD=zeros(info.rawchAmt,1); %<<<
    parfor chi=1:info.rawchAmt
        [rawSD{chi},rawSA{chi},rawSW{chi},rawSQ{chi},tp]=spike_detect(X(:,chi),info.srate,paras);
        sigSD(chi)=tp.SSD;
        fprintf('|');
    end
    fprintf('\n');
    
    % Save raw SD and information
    outfile2.rawSD=rawSD; outfile2.rawSA=rawSA; outfile2.rawSW=rawSW; outfile2.rawSQ=rawSQ;
    outfile2.sigSD=sigSD;
    outfile2.info=info;
    outfile2.paras=paras;
    % Write to disk
    save(fnS,'-struct','outfile2'); %,'-v7.3'
    disp('detection done');
end


%%%%% This section is preparation in form of alignment cutting.
% the immediate cutting is in "cluster" stage.
if runFlag(funcNumTab.align)    
    % Set up input
    if ~exist('X','var')
        fprintf('loading raw signal from %s\n',fnF);
        load(fnF,'X');
        load(fnF,'T');
    end
    if ~exist('outfile2','var')
        outfile2=load(fnS);
    end
    info=outfile2.info;
    
    disp('Get window length by spike alignment test >>>');
    chinfo=struct('rawchAmt',info.rawchAmt);
    
    % Remove low activity channels and save results as SD
    thr=paras.STactThre*info.TimeSpan;
    sa=cellstat(outfile2.rawSD,'length');
    I=(sa>=thr); sa=sa(I);
    SD=outfile2.rawSD(I); NSA0=outfile2.rawSA(I); % * SQ not need to be saved for further use.    
    chinfo.AlichAmt=sum(I);
    temp=(1:chinfo.rawchAmt)';    chinfo.AlichID=temp(I);
    
    % Get window Spike alignment (纯用于获取窗口宽度，选最少spike的通道即可)
    [~,idx]=min(sa);
    [~,~,SO]=spike_align(X(:,1),SD(idx),info.srate,'window',paras.alignWin);
    % 消除窗口无法覆盖的边缘spike，同步更新SD,SA等信息。    
    for chi=1:chinfo.AlichAmt
        % * must remove large indexed SD first.
        rmlist=(SD{chi}+SO.postww>info.ptsAmt);
        SD{chi}(rmlist)=[];
        NSA0{chi}(rmlist)=[];
        rmlist=(SD{chi}<=SO.preww);
        SD{chi}(rmlist)=[];
        NSA0{chi}(rmlist)=[];
    end
    
    % Update with new data
    outfile2.SD=SD; outfile2.SA=NSA0;
    info.spklen=SO.spklen; info.spkprew=SO.preww; info.spkpostw=SO.postww;
    outfile2.info=info;
    outfile2.chinfo=chinfo;
    
    % Save results
    save(fnS,'-v7.3','-struct','outfile2');    
    disp('data alignment done');
end


%%%%%%% Spike merging 
if runFlag(funcNumTab.merge)
    if ~exist('outfile2','var')
        outfile2=load(fnS);
    end
    info=outfile2.info;
    chinfo=outfile2.chinfo;
    
    %%% in this 1st clustering, do not separate noise spikes in final results,
    % keep all noises and use spike merge to see if can delete them.
    if paras.bClusterBeforeMerge
        saveopt=paras.bAlwaysMatch;
        paras.bAlwaysMatch=true; % match all spikes.
        
        % Set up input
        if ~exist('X','var')
            fprintf('loading raw signal from %s\n',fnF);
            load(fnF,'X');
            load(fnF,'T');
        end
        % Set up output
        outfile3=struct();
        
        %%% Clustering of each channel
        disp('spike pre-clustering for merging >>>');
        CST=cell(chinfo.AlichAmt,1);
        SD=outfile2.SD; tx=X(:,chinfo.AlichID); % special for parfor parallel
        parfor k=1:chinfo.AlichAmt
            [CST{k},SD{k}]=chcluster(tx(:,k),SD{k},info,paras);
            fprintf('|');
        end
        fprintf('\n');
        
        % Corrected SD.
        outfile3.SD0=SD; % their chID is AlichID.
%         outfile3.CST0=CST; no necessary <<<
%         outfile3.SA0=cell(chinfo.AlichAmt,1);
%         for chi=1:chinfo.AlichAmt
%             outfile3.SA0{chi}=X(SD{chi},chinfo.AlichID(chi));
%         end
        
        %%% 生成新的NSD（每个分开的"unit"算一个通道）and new chID and chAmt.
        [outfile3.NSD0,tp]=getNSD(SD,CST);
        chinfo.NSD0chID=chinfo.AlichID(tp);
        chinfo.NSD0chAmt=length(chinfo.NSD0chID);
        
        % Add amplitude of NSD (NSA), for use in spike merge.
        outfile3.NSA0=cell(chinfo.NSD0chAmt,1);
        for chi=1:chinfo.NSD0chAmt
            outfile3.NSA0{chi}=X(outfile3.NSD0{chi},chinfo.NSD0chID(chi));
        end
        outfile3.info=info;
        outfile3.chinfo=chinfo;
        outfile3.paras=paras;
        
        % save
%         save(fnC,'-struct','outfile3');
        paras.bAlwaysMatch=saveopt;
        disp('pre-clustering done');
    end
    

    %%%%%%%%%%%%
    % *!!! (Merging is put here only because it uses biggest amplitude as
    % standard to keep spikes in a cluster. -this can do the merging when
    % unsorted data is used.) 
    if ~exist('outfile3','var')
        outfile3=load(fnC);
        chinfo=outfile3.chinfo;
    end
    if ~exist('X','var')
        fprintf('loading raw signal from %s\n',fnF);
        load(fnF,'X');
        load(fnF,'T');
    end
    
    disp('spike merge >>>');
%     % Get correct spike peak location by time (if necessary)
%     [A,~,O]=spike_align(X,outfile3.NSD0,info.srate,'chAssign',chinfo.NSD0chID,'window',[-1,1]);
%     checkstat=spikecs_check(A,info);
%     % 由于spike数量少而stype==0的通道进行location重置。    
%     I=(checkstat.stype==0);
%     checkstat.ploc(I)=O.preww+1;        
%     % Find channel whose location is changed.
%     tloc=checkstat.ploc-(O.preww+1);
%     tprange=(0.1/1000)*(info.srate);
%     changeI=(abs(tloc)>=tprange);
%     for chi=1:chinfo.NSD0chAmt
%         if changeI(chi)
%             temp=A{chi};
%             I=checkstat.ploc(chi)-tprange:checkstat.ploc(chi)+tprange;
%             I=I(I>=1);
%             I=I(I<=O.preww+O.postww+1);
%             temp=temp(I,:);
%             [~,idx]=max(abs(temp));
%             % index to start of A
%             idx=I(idx);
%             % distance to original position
%             dt=idx-(O.preww+1);
%             
%             % Change the NSDO accordingly.
%             outfile3.NSD0{chi}=outfile3.NSD0{chi}+dt;
%         end
%     end

    % Accurate time by up-sampling.
    [NST0,NSA0]=exactST(X,outfile3.NSD0,T,info.srate,chinfo.NSD0chID);
    % Merge function
    [~,outfile3.NSD0csmark,outfile3.NSD0csrmlb]=spikemerge(NST0,NSA0,info,'time span',info.TimeSpan);
    NSD=outfile3.NSD0;
    for chi=1:chinfo.NSD0chAmt
        NSD{chi}(outfile3.NSD0csrmlb{chi}~=0)=[];
    end
    
    %%% 需要将NSD重新组合成reconSD以利于2nd clustering 或 display.
    reconSD0=cell(chinfo.rawchAmt,1);
    reconD0CSrmlb=cell(chinfo.rawchAmt,1);% 通过重组NSD0到reconSD0得到reconSD0数据对应的rmlb.
    for chi=1:chinfo.rawchAmt
        I=find(chinfo.NSD0chID==chi);
        if ~isempty(I)
            temp=[]; temp2=[]; temp3=[];
            for k=1:length(I)
                temp=[temp;NSD{I(k)}];
                temp3=[temp3;outfile3.NSD0{I(k)}];
                temp2=[temp2;outfile3.NSD0csrmlb{I(k)}];
            end
            [reconSD0{chi}]=sort(temp(:),'ascend');
            [~,idx]=sort(temp3,'ascend');
            reconD0CSrmlb{chi}=temp2(idx);
        end
    end
    % Delete any very sparse spike channels
    sa=cellstat(reconSD0,'length');
    I=(sa>10);
    reconSD0=reconSD0(I);
    temp=(1:chinfo.rawchAmt)';
    chinfo.MergechID=temp(I);    chinfo.MergechAmt=length(chinfo.MergechID);
    
    % Save
    outfile3.reconSD0=reconSD0;
    outfile3.chinfo=chinfo;
    outfile3.reconD0CSrmlb=reconD0CSrmlb;
    save(fnC,'-struct','outfile3');
    NST0chID=chinfo.NSD0chID;
    save(fnST,'NST0','NSA0','NST0chID','info');
    disp('merging done');
end


%%%%%%% Spike clustering
% Second round of clustering - the noisy spikes can be deleted in final of this round.
if runFlag(funcNumTab.cluster)
    % Set up input
    if ~exist('X','var')
        fprintf('loading raw signal from %s\n',fnF);
        load(fnF,'X');
        load(fnF,'T');
    end
    if ~exist('outfile3','var')
        outfile3=load(fnC);
    end    
    info=outfile3.info;
    chinfo=outfile3.chinfo;
    
    disp('2nd spike cluster >>>');
    CST=cell(chinfo.MergechAmt,1);    reconSD=CST; reconAL=CST;
    % special for parfor parallel
    SD=outfile3.reconSD0; tx=X(:,chinfo.MergechID);
    parfor k=1:chinfo.MergechAmt        
        [CST{k},reconSD{k},reconAL{k}]=chcluster(tx(:,k),SD{k},info,paras);
        fprintf('|');
    end    
    fprintf('\n');
    temp=cell(chinfo.rawchAmt,1); temp(chinfo.MergechID)=CST; outfile3.CST=temp;
    temp=cell(chinfo.rawchAmt,1); temp(chinfo.MergechID)=reconSD; outfile3.reconSD=temp;
    temp=cell(chinfo.rawchAmt,1); temp(chinfo.MergechID)=reconAL; outfile3.reconAL=temp;
    % Get reconSA
    reconSA=cell(chinfo.rawchAmt,1);
    for k=1:chinfo.MergechAmt
        chi=chinfo.MergechID(k);
        reconSA{chi}=X(reconSD{k},chi);
    end
    outfile3.reconSA=reconSA;
    
    %%% 生成新的NSD（每个分开的"神经元"算一个通道）and new chID and chAmt.
    assert(isequal(cellstat(reconSD,'length'),cellstat(CST,'length')),'reconSD and CSI length not match after clustering');
    [outfile3.NSD,tp]=getNSD(reconSD,CST);
    chinfo.CluchID=chinfo.MergechID(tp);
    chinfo.CluchAmt=length(chinfo.CluchID);
    outfile3.chinfo=chinfo;
    
    % save
    save(fnC,'-struct','outfile3');
    disp('clustering done');
end


%%%%%%%
if runFlag(funcNumTab.toST)
    if ~exist('X','var')
        fprintf('loading raw signal from %s\n',fnF);
        load(fnF,'X');
        load(fnF,'T');
    end
    if ~exist('outfile2','var')
        outfile2=load(fnS);
    end
    if ~exist('outfile3','var')
        outfile3=load(fnC);
    end
    info=outfile3.info;
    chinfo=outfile3.chinfo;
    
    %%% Save the ST form of spike_detect.
    rawSDTchID=1:chinfo.rawchAmt;
    if paras.bExactST % Get Precise time by splining spike peaks.
        [rawSDT,rawSDA]=exactST(X,outfile2.rawSD,T,info.srate,rawSDTchID);
        save(fnST,'-append','rawSDT','rawSDA');
    else
        rawSDT=idx2time(outfile2.rawSD,T);
        save(fnST,'-append','rawSDT');
    end
    disp('detected SD converted to ST');    
    
    %%% Save the ST form of spike_cluster.
    if ~exist('outfile3','var')
        outfile3=load(fnC); % loading right "info" is important.
    end    
    % First delete low-activity channels
    temp=cellstat(outfile3.NSD,'length');
    thr=paras.STactThre*info.TimeSpan;
    I=(temp>thr);
    NSD=outfile3.NSD(I);
    NSTchID=chinfo.CluchID(I);    
    % Remaining converted to ST
    if paras.bExactST % Get Precise time by splining spike peaks.
        [NST,NSA]=exactST(X,NSD,T,info.srate,NSTchID);
    else
        stAmt=length(NSTchID);
        NST=idx2time(NSD,T);
        NSA=cell(stAmt,1);
        for chi=1:stAmt
            NSA{chi}=X(NSD{chi},NSTchID(chi));
        end
    end
    save(fnST,'-append','NST','NSTchID','info','NSA');
    disp('sorted SD converted to ST');
    
    %%% Save reconstructed ST (after final clustering).
    reconSD=outfile3.reconSD;
    temp=cellstat(reconSD,'length');
    I=(temp>thr);
    reconSD=reconSD(I);
    reconSTchID=rawSDTchID(I);
    % Convert to ST.
    if paras.bExactST
        [reconST]=exactST(X,reconSD,T,info.srate,reconSTchID);
    else
        reconST=idx2time(reconSD,T);
    end
    save(fnST,'-append','reconST','reconSTchID');
    disp('reconstructed SD converted to ST');
end

% Cleaning up
% delete(sspath);
disp('run finish, all done');

end %%%%%%% End of main.