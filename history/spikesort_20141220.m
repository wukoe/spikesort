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
paras.alignWin=[-0.8,1.1]; % (ms)
paras.bAlignSmooth=false;

% Spike feature
paras.feaDim=3;

% Spike cluster
paras.cluMinSpkThres=6; % cluster number minimum (/min)
paras.bTimeShiftMerge=true; 
paras.lagWin=[-0.8,0.8]; % shift merging's window (ms)
paras.ccThres=0.8; %<<< cluster mean curve correlation threshold 0.85
paras.SNratioThres=10;

%%% Names of files used
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
if ~runFlag(funcNumTab.filt) && ~runFlag(funcNumTab.spikeDetect) && exist([fnS,'.mat'],'file') % load processed information from _s.mat file
    load(fnS,'info');
else % directely from .mcd file
    info=getMCDinfo([fileName,'.mcd'],bDirectMem);
    info.chAmt=info.rawchAmt;
    info.chLabel=info.rawchLabel;
end

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


%%% Process of others
paras.lagWin=floor(paras.lagWin/1000*info.srate);


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
    % 消除窗口无法覆盖的边缘spike，同步更新SD,SA,SQ等信息。
    if ~isempty(rmlist)
        SD=outfile2.SD; SA=outfile2.SA; SQ=outfile2.SQ;
        % * must do it in reversed order
        for k=size(rmlist,1):-1:1
            SD{rmlist(k,1)}(rmlist(k,2))=[];
            SA{rmlist(k,1)}(rmlist(k,2))=[];
            SQ{rmlist(k,1)}(rmlist(k,2))=[];
        end
        if bDirectMem
            save(fnS,'-v7.3','-append','SD','SA','SQ');
        else
            outfile2.SD=SD;
            outfile2.SA=SA;
            outfile2.SQ=SQ;
        end
    end
    
    if bDirectMem
        save(fnS,'-v7.3','-append','A');
        save(fnS,'-v7.3','-append','paras');
    end
    disp('data alignment done');
end


%%%%%% #3 Filtering for LFP
if false; %runFlag(funcNumTab.filtLFP)
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
    disp('spike feature >>>');
    % Set input
    if ~exist('A','var')
        fprintf('loading aligned spike data from %s\n',fnS);
        load(fnS,'A');
    end
    if ~exist('SA','var')
        load(fnS,'SA');
    end
    
    SF=cell(info.chAmt,1);
    for chi=1:info.chAmt
        if isempty(A{chi,1})
            SF{chi,1}=[];
            fprintf('X');
        else
            SF{chi,1}=spike_feature(A{chi,1},'dim',paras.feaDim);
            % 添加Spike amplitude 信息！
            SF{chi,1}=[SF{chi,1},SA{chi}];
            SF{chi,1}=zscore(SF{chi,1});
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
    
    %%% Proc
    totalTime=(time(2)-time(1)); % total time span of data        
    % cluster number minimum is based on spikes/min, so the actual number
    % should be calculated by time length
    paras.cluMinSpkThres=paras.cluMinSpkThres*(totalTime/60);
    
    % Channel SI: cluster identity of each spike in raw channel.
    CSI=cell(info.chAmt,1); % * note there may be 0 in CSI - those not assigned to any cluster.       
    % newchID: from which raw channel (electrode) the neuron is recorded.
    newchID=zeros(0,1);
    % NSD is neuron spk data (index format)
    newSD=cell(0,1);

    %%% Clutering each channel. 
    chcount=0; % number of effective neurons - 数量不能预知，所以累计计数。
    for chi=1:info.chAmt
        tp=size(SF{chi},1);
        if tp<=paras.cluMinSpkThres % spike数量过少的通道直接复制原值（不删除为好）。
            chcount=chcount+1;
            newSD{chcount,1}=SD{chi};            
            newchID(chcount,1)=info.chID(chi);
            CSI{chi}=zeros(tp,1);
        else
            [CSI{chi},lbinfo]=spike_cluster(SF{chi,1},'kmeans',paras.cluMinSpkThres);
            
            %%%%%%%%%%% Shift the signal and get distance measure - could be merged after that.
            if paras.bTimeShiftMerge
            spkData=A{chi,1};
            %%% Recursive process to combine the pair of smallest
            while lbinfo.cAmt>1
                %%% Get mean curve of each cluster
                cm=zeros(size(spkData,1),lbinfo.cAmt);
                for k=1:lbinfo.cAmt
                    cm(:,k)=mean(spkData(:,lbinfo.ids{k}),2);
                end

                %%% 计算每一对是否应该合并 - 找到最接近的那一对。
                sameSpkMark=false(lbinfo.cAmt); % mark whether 2 clusters (nominated by row and column) should merge.
                cc=zeros(lbinfo.cAmt); % corr coef
                ML=zeros(lbinfo.cAmt); % lag value
                for m=1:lbinfo.cAmt-1
                    for n=m+1:lbinfo.cAmt
                        %%% Check the shape similarity by CC (use the mean curve of 2 clusters)
                        % move along time to get max CC and min distance
                        [tp1,tp2]=lagcorr([cm(:,m),cm(:,n)],paras.lagWin);
                        cc(m,n)=tp1(1,2); ML(m,n)=tp2(1,2);
                        if cc(m,n)>paras.ccThres
                            %%% Also check absolute difference between curves - this requirement
                            % should also be fullfilled.
                            % * This is because 2 types with similar shape but different
                            % amplitude is also different

                            % Align 2 curves according to ML
                            if ML(m,n)>0 % m behind n
                                tp1=cm(1:ML(m,n),m); tp2=cm(end-ML(m,n)+1:end,n);
                            else
                                tp1=cm(end-ML(m,n)+1:end,m); tp2=cm(1:ML(m,n),n);
                            end
                            % overlap part length
                            ol=length(tp1);
                            % average difference (of time)
                            S=mean(abs(tp1-tp2));
                            % amplitude of signals (max amplitude of two curves) -
                            % measure "signal' in SNR
                            tp=[tp1,tp2];
                            D=max(max(tp)-min(tp));                

                            if S<=D/paras.SNratioThres
                                sameSpkMark(m,n)=true;
                            end
                        end
                    end
                end

                % 若没有任何可以合并的，放弃这一通道
                if sum(sum(sameSpkMark))==0
                    break
                end

                %%% Merge the "same spike" pair with highest cc value
                cc(~sameSpkMark)=0; %<<<<
                [~,idx]=matmax(cc,'triu');

                % 以cluster成员数量多者为时间轴对齐标准
                num=lbinfo.typeAmt(idx);
                [~,movidx]=min(num); % movidx是被移动者的（在idx中的）序号，值只有1、2两个选项
                movidx=idx(movidx); % 变成具体要移动哪个cluster。
                movpts=abs(ML(idx(1),idx(2)));
                
                % 获得被移动者需要移动的距离（每个spike单独进行，cluster的距离作为限制搜索范围的参考）
                if (ML(idx(1),idx(2))<0 && movidx==1) || (ML(idx(1),idx(2))>0 && movidx==2)%被移动者向前
                    temp=spkData(:,lbinfo.ids{movidx});
                    spkData(:,lbinfo.ids{movidx})=[temp(movpts+1:end,:); zeros(movpts,lbinfo.typeAmt(movidx))];
                else
                    temp=spkData(:,lbinfo.ids{movidx});
                    spkData(:,lbinfo.ids{movidx})=[zeros(movpts,lbinfo.typeAmt(movidx)); temp(1:end-movpts,:)];
                end
%                 
%                 if ML(idx(1),idx(2))<0 % row proceed column
%                     if movidx==1 % row不必动，move column upstream
%                         temp=spkData(:,lbinfo.ids{idx(2)});
%                         spkData(:,lbinfo.ids{idx(2)})=[temp(movpts+1:end,:); zeros(movpts,lbinfo.typeAmt(idx(2)))];
%                     else % move row down
%                         temp=spkData(:,lbinfo.ids{idx(1)});
%                         spkData(:,lbinfo.ids{idx(1)})=[zeros(movpts,lbinfo.typeAmt(idx(1))); temp(1:end-movpts,:)];
%                     end
% 
%                 elseif ML(idx(1),idx(2))>0 % row behind column
%                     if movidx==1 % move column down
%                         temp=spkData(:,lbinfo.ids{idx(2)});
%                         spkData(:,lbinfo.ids{idx(2)})=[zeros(movpts,lbinfo.typeAmt(idx(2))); temp(1:end-movpts,:)];
%                     else % move row up
%                         temp=spkData(:,lbinfo.ids{idx(1)});
%                         spkData(:,lbinfo.ids{idx(1)})=[temp(movpts+1:end,:); zeros(movpts,lbinfo.typeAmt(idx(1)))];
%                     end
%                 end

                % 合并
                if movidx==1
                    lbinfo.typeAmt(idx(1))=lbinfo.typeAmt(idx(1))+lbinfo.typeAmt(idx(2)); lbinfo.typeAmt(idx(2))=[];
                    lbinfo.ids{idx(1)}=[lbinfo.ids{idx(1)};lbinfo.ids{idx(2)}];  lbinfo.ids(idx(2))=[];
                    lbinfo.types(idx(1))=lbinfo.types(idx(1))+lbinfo.types(idx(2)); lbinfo.types(idx(2))=[];
                else
                    lbinfo.typeAmt(idx(2))=lbinfo.typeAmt(idx(2))+lbinfo.typeAmt(idx(1)); lbinfo.typeAmt(idx(1))=[];
                    lbinfo.ids{idx(2)}=[lbinfo.ids{idx(2)};lbinfo.ids{idx(1)}]; lbinfo.ids(idx(1))=[];
                    lbinfo.types(idx(2))=lbinfo.types(idx(2))+lbinfo.types(idx(1)); lbinfo.types(idx(1))=[];
                end
                lbinfo.cAmt=lbinfo.cAmt-1;
            end
            end
            
            %%% % * it is possible that non is left after clustering. 
            if lbinfo.cAmt>0 
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
    save(fnC,'-v7.3','-append','paras');
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
    case 'spike cluster'
        N=funcNumTab.cluster;
    case 'final'
        N=funcNumTab.final;
    otherwise
        error('invalid run option');
end

end