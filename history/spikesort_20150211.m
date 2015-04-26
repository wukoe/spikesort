% the whole process of spike sorting
% the basic stratege is to process the data one channel by one channel
%   runstat=spikesort(fileName,varargin)
% fileName does not have postfix
% 'RunStep': default=':'; {'filt','spike detect', 'spike align',
% 'spike feature', 'cluster'}, or use 'resume'
% 'GoOn': default='on'; {'off'}
function spikesort(fileName,varargin)
%%%%%%%%%% Parameter Default
%%% 运行相关
% bFullProc=true;
runOpt=':';
bUserSpecifiedMem=false;

% Process control usage
% run flag - marking which of processings are needed
runFlag=false(11,1); % default: all stop
% function-number table
funcNumTab=struct('start',1,'filt',2,'spikeDetect',4, 'spikeAlign',5, ...
     'cluster',8, 'final',10); % 'spikeFeature',6,

%%% 处理步骤的?
paras=struct();
% filtering
paras.spkFiltFreq=[200,5000];
paras.bRmc=false; % whether to remove big noise channel
paras.rmChThres=5/1000; % (V)

% Spike detect
paras.movThres=6;
paras.bPNSep=true; % whether to separate the positive and negative peaks to separate channels first.
% *this option is related to feature selection.

% Spike align
paras.alignWin=[-0.6,1]; % [-0.8 1.1] (ms)
paras.bAlignSmooth=false;

% Spike cluster
paras.cluMinSpkThres=12; % cluster number minimum (/min)
paras.feaDim=3;% Spike feature
paras.bDrawForCluster=false;
paras.drawNum=4000;
paras.bTimeShiftMerge=true;
paras.lagWin=[-0.5,0.5]; % shift merging's maximum range (ms)
paras.ccThres=0.8;%0.91; %<<< cluster mean curve correlation threshold 0.85
paras.SNratioThres=8; % 17; % 17 is tested.
paras.bIndivMove=false; % move each spike with individual distance or not.

%%% Names of files used
fnRaw=[fileName,'.mcd'];
fnF=[fileName,'_f'];
% fnFL=[fileName,'_fl'];
fnS=[fileName,'_s'];
% fnFEA=[fileName,'_fea'];
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
%                 bFullProc=false;
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

% process fileName
[~,fileName,~]=fname(fileName);


%%% #1 Initiate the Run status data, Parameter. 这是每次运行必须的
% Get info. struct data.
% * 如果要运行filt或spike detect步骤，则从原始文件中获得信息；如果是其他步骤，则查看fnS文件是否存在，若是，则从中加载现成的info
if ~runFlag(funcNumTab.filt) && ~runFlag(funcNumTab.spikeDetect) && exist([fnS,'.mat'],'file') % load processed information from _s.mat file
    load(fnS,'info');
else % directely from .mcd file
    info=getMCDinfo([fileName,'.mcd'],bDirectMem);
    info.chAmt=info.rawchAmt;
end
paras.bDirectMem=bDirectMem;

% save the source file name that data is from
info.fileName=fileName; 

%     % Check completeness of key fields
%     tp=isfield(info,{'chAmt','chID','ptsAmt','srate'});
%     if ~multilogic(tp,'and',2)
%         error('the information for processing is not complete');
%     end

%%% Process of others
paras.lagWin=floor(paras.lagWin/1000*info.srate);

totalTime=info.TimeSpan; 
% cluster number minimum is based on spikes/min, so the actual number should consider time length
paras.cluMinSpkThres=paras.cluMinSpkThres*(totalTime/60);


%%%%%%%%%%%% Main
% 输出数据对应文件名称
% outfile1=FnF; outfile2=FnS; outfile3=FnC;
if runFlag(funcNumTab.filt)    
    % Set up input
    [nsresult, infile1] = ns_OpenFile(fnRaw);
    if nsresult==-1,     error('open .mcd file error'); end
    % Set up output
    if bDirectMem
        outfile1=struct('X',[],'T',[]);
    else
        outfile1=matfile(fnF,'Writable',true);
    end
    
    run_Filt();
end
if runFlag(funcNumTab.spikeDetect)
    % Set up input - if no filtered data (outfile1) in memory, try load it.
    loadoutfile1();
    % Set up output
    if bDirectMem
        outfile2=struct();
    else
        outfile2=matfile(fnS,'Writable',true);
    end
    
    run_Detect();
end
if runFlag(funcNumTab.spikeAlign)
    % Set up input
    loadoutfile1(); 
    if ~exist('outfile2','var')
        if bDirectMem
            outfile2=load(fnS);
        else
            outfile2=matfile(fnS);
        end
    end
    
    run_Align();
end

if runFlag(funcNumTab.cluster)
    % Set up input
%     if ~exist('SF','var'), load(fnS,'info','A','SA','SD'); end % to produce SF
    if ~exist('outfile2','var')
        outfile2=load(fnS);
        info=outfile2.info;
    end
    if ~exist('T','var'), load(fnF,'T'); end
        
    run_Cluster();
end

%%%%%%%%%%%%%% Whether to enter the spike time-shift merge stage
if runFlag(funcNumTab.cluster) && paras.bTimeShiftMerge
    loadoutfile1();
    load(fnC,'info','reconSD','CSI');
    % Update to get the new alignment data
    A=spike_align(outfile1.X,reconSD,info.srate,'window',paras.alignWin,'bSmooth',paras.bAlignSmooth);
    % * 用了reconSD, 就不需要'chAssign'选项了。
    
    newSD=reconSD;
    for chi=1:info.rawchAmt
        spkclu=reabylb(CSI{chi});
        % 排除噪声选项？ <<<<<<
        idx=find(spkclu.types==0);
        if ~isempty(idx)
            spkclu.types(idx)=[];
            spkclu.typeAmt(idx)=[];
            spkclu.cAmt=spkclu.cAmt-1;
            spkclu.ids(idx)=[];
        end
        
        %%% Recursive process to combine the pair of smallest distance
        if spkclu.cAmt>1
            spkData=A{chi};
            %%% Get mean curve of each cluster
            cm=zeros(info.spklen,spkclu.cAmt);
            for k=1:spkclu.cAmt
                cm(:,k)=mean(spkData(:,spkclu.ids{k}),2);
            end
            
            %%% 计算每一对是否应该合并 - 找到最接近的那一对。
            sameSpkMark=false(spkclu.cAmt); % mark whether 2 clusters (nominated by row and column) should merge.
            cc=zeros(spkclu.cAmt); % corr coef
            D=cc;
            ML=zeros(spkclu.cAmt); % lag value
            for m=1:spkclu.cAmt-1
                for n=m+1:spkclu.cAmt
                    %%% Check the whether 2 shapes should merge (use the mean curve of 2 clusters)
                    % First Move along time to get max CC and related distance
                    [tp1,tp2]=lagcorr([cm(:,m),cm(:,n)],paras.lagWin);
                    cc(m,n)=tp1(1,2); ML(m,n)=tp2(1,2);
                    if cc(m,n)>paras.ccThres                        
                        % Also check absolute difference between curves (* this requirement
                        % should also be fullfilled, because 2 types with similar shape but different
                        % amplitude is also different).                        
                        % Align 2 curves according to ML.
                        if ML(m,n)>0 % m behind n
                            tp1=cm(ML(m,n)+1:end,m); tp2=cm(1:end-ML(m,n),n);
                        else
                            tp1=cm(1:end+ML(m,n),m); tp2=cm(1-ML(m,n):end,n);
                        end                        
                        % average difference (of time)
                        D(m,n)=mean(abs(tp1-tp2));
                        % amplitude of signals (max amplitude of two curves) -
                        % measure "signal' in SNR
                        tp=[tp1,tp2];
                        S=max(max(tp)-min(tp));
                        
                        if D(m,n)<=S/paras.SNratioThres
                            sameSpkMark(m,n)=true;
                        end
                    end
                end
            end
            ML=ML-ML'; % fill the bottom triangle of lag matrix
            
            % 若有任何可以合并的，则进行处理
            if sum(sum(sameSpkMark))>0
                % Get the direction of merging (which cluster merged into
                % which)  * method: cluster with more member wins. (others
                % to consider: less within-cluster variance wins)
                MD=false(spkclu.cAmt);
                for m=1:spkclu.cAmt-1
                    for n=m+1:spkclu.cAmt
                        if sameSpkMark(m,n)
                            if spkclu.typeAmt(m)>spkclu.typeAmt(n)
                                MD(n,m)=1; % n->m
                            else
                                MD(m,n)=1;
                            end
                        end
                    end
                end
                
                % Collapse single source (e.g., if 1->2, 2->3, then get
                % 1->3 instead of 1->2).
                for m=1:spkclu.cAmt
                    idx=find(MD(m,:));
                    while ~isempty(idx)
                        if length(idx)>1
                            disp('warning for multi-target'); 
                            %idx=idx(1); 
                            % max number/best fit
                            temp=D(m,idx);
                            [~,tp]=min(temp);
                            idx=idx(tp);
                        end
                        % move to new target
                        MD(m,:)=0; MD(m,idx)=1; 
                        idx=find(MD(idx,:));
                    end
                end
                
                %%% Merge the "same spike" pair
                for m=1:spkclu.cAmt
                    movidx=find(MD(m,:));
                    if ~isempty(movidx)
                        movpts=ML(m,movidx); %movpts P/N property decide whether the data move forward/backward.
                        % Get more accurate move distance by matching the
                        % peaks. <<<<<<<<<<<
                        peakloc=info.spkprew+1; % location of peak of spike data
                        tp=peakloc+movpts;
                        [~,idx]=max(abs(cm(tp-2:tp+2,m)));
                        if idx ~= 3 % which match the number 4 used
                            movpts=movpts+(idx-3);
                        end                            
                        
                        % 获得被移动者需要移动的距离
                        if paras.bIndivMove %每个spike单独进行，cluster的距离作为限制搜索范围的参考
                            for k=1:spkclu.typeAmt(m)
                                % Decide the move length of this spike
                                % 1/N use max correlation location
                                sp=spkData(:,spkclu.ids{m}(k));
                                [tp,mp]=lagcorr([cm(:,movidx),sp],[-movpts-2,-movpts+2]);
                                tp=tp(1,2); %mp=mp(1,2);
                                % 2/N make peaks match.
                                
                                % e/N
                                
                                if true % tp>=paras.ccThres
                                % Make re-location
                                temp=newSD{chi}(spkclu.ids{m}(k));
                                temp=temp-mp; % * note this is different from batch method below
                                newSD{chi}(spkclu.ids{m}(k))=temp;
                                end
                            end
                            
                        else % move in a batch.
                            newSD{chi}(spkclu.ids{m})=newSD{chi}(spkclu.ids{m})+movpts;
                        end
%                         CSI{chi}(spkclu.ids{m})=spkclu.types(movidx);                        
                    end
                end
            end
            
        end
    end
    
    %%% Update those unsorted information as well.
    SD=newSD;
    info.chID=1:info.rawchAmt;
    SA=cell(info.rawchAmt,1);
    for chi=1:info.rawchAmt
        SA{chi}=outfile1.X(SD{chi},chi);
    end
    save(fnS, '-append','info','SD', 'SA');
    
    %%% Do clustering again
    outfile2=load(fnS);
    info=outfile2.info;
    run_Align();
    info=outfile2.info;
    run_Cluster();
end

% Cleaning up
% delete(sspath);
disp('run finish, all done');
%%%%%%%%%%%% End of main.


%%%%%%%%%%%%%% Separate stages
%%%%%% #2 Filtering and removal of bad channels
    function run_Filt()
        disp('filt >>>');
        % Get filter
        [paras.fb,paras.fa]=butter(4,paras.spkFiltFreq/(info.srate/2));
        
        % add time points data
        outfile1.T=(1:info.ptsAmt)'/info.srate; % measured in (s)
        
        %%% Do the filtering
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
        
        % Remove noisy channel of X if any
        if sum(rmMark)>0
%             outfile1.X(:,rmMark)=[];
            
            % Update info
            info.badChannel=rmMark;
%             info.chAmt=info.chAmt-sum(rmMark);
        end
        
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


%%%%% #4 Spike detect
    function run_Detect()        
        disp('spike detect >>>');
        info=outfile1.info;
        rawSD=cell(info.rawchAmt,1);
        rawSQ=cell(info.rawchAmt,1);
        rawSA=cell(info.rawchAmt,1);        
        for k=1:info.rawchAmt
            chi=info.chID(k);
            [rawSD{chi},rawSQ{chi},rawSA{chi}]=spike_detect(outfile1.X(:,chi),info.srate,'movThres',paras.movThres);
            fprintf('|');
        end
        fprintf('\n');
        % Save raw SD
        outfile2.rawSD=rawSD; outfile2.rawSQ=rawSQ; outfile2.rawSA=rawSA;
        
        % Save SD for next stage
        info.chAmt=info.rawchAmt;
        info.chID=(1:info.rawchAmt)';
        outfile2.info=info;
        outfile2.SD=rawSD;
        outfile2.SA=rawSA; % * SQ not need to be saved for further use.
        
        if bDirectMem
            save(fnS,'-v7.3','-struct','outfile2');
        end
        disp('detection done');        
    end


%%%%% Spike alignment is separated here.
    function run_Align()
        bST=false;
        
        %%% Remove low activity channels and save results as SD
        sa=cellstat(outfile2.SD,'length');
        I=(sa>=5);
        SD=outfile2.SD(I); SA=outfile2.SA(I); % * SQ not need to be saved for further use.
        info=outfile2.info;
        info.chAmt=sum(I);
        info.chID=info.chID(I);        
        
        
        %%% Splite the positive and negative spikes
        if paras.bPNSep
            disp('split N/P spikes >>>');
            newSD=cell(0);
            newSA=cell(0);
            newchID=zeros(0,1);
            
            chcount=0;
            for chi=1:info.chAmt
                PM=(SA{chi}>0); % get posi/nega mark
                pnum=sum(PM); % number of positive peaks
                if pnum>0
                    chcount=chcount+1;
                    
                    % assign to the newSD
                    newSD{chcount,1}=SD{chi}(PM);
                    % assign to new sA
                    newSA{chcount,1}=SA{chi}(PM,:);
                    % assign new chID
                    newchID(chcount,1)=info.chID(chi);
                end
                
                if pnum<length(SD{chi}) % number of negative peaks > 0
                    chcount=chcount+1;
                    
                    PM=~PM;
                    % assign to the newSD
                    newSD{chcount,1}=SD{chi}(PM);
                    % assign to new sA
                    newSA{chcount,1}=SA{chi}(PM,:);
                    % assign new chID
                    newchID(chcount,1)=info.chID(chi);
                end
            end
            
            info.chID=newchID; info.chAmt=chcount;            
            outfile2.SD=newSD;
            outfile2.SA=newSA;
            disp('splited');
        end
        
        outfile2.info=info;        
        if bDirectMem
            save(fnS,'-v7.3','-struct','outfile2');
        end
        
        %%%
        disp('spike alignment >>>');
        if bDirectMem
            [A,rmlist,SO]=spike_align(outfile1.X,outfile2.SD,info.srate,'chAssign',info.chID,'window',paras.alignWin,'bSmooth',paras.bAlignSmooth);
        else % 由于matfile方法处理2D矩阵的限制？
            [A,rmlist]=spike_align(outfile1,outfile2.SD,info.srate,'chAssign',info.chID,'window',paras.alignWin,'bSmooth',paras.bAlignSmooth);
        end
        outfile2.A=A;
        % 消除窗口无法覆盖的边缘spike，同步更新SD,SA,SQ等信息。
        if ~isempty(rmlist)
            SD=outfile2.SD; SA=outfile2.SA;
            % * must do it in reversed order
            for k=size(rmlist,1):-1:1
                SD{rmlist(k,1)}(rmlist(k,2))=[];
                SA{rmlist(k,1)}(rmlist(k,2))=[];
            end
            save(fnS,'-v7.3','-append','SD','SA');
        end
        
        % Get spike morphology data length.
        tp=cellstat(outfile2.SD,'length');
        idx=find(tp>0,1);
        info.spklen=SO.spklen; info.spkprew=SO.preww; info.spkpostw=SO.postww;
        outfile2.info=info;
        
        if bDirectMem
            if bST
                % Get Precise time by splining spike peaks
                ST=exactST(outfile1.X,outfile2.SD,outfile1.T,info.srate,info.chID);
                save(fnS,'-v7.3','-append','A','info','ST');
            else
                save(fnS,'-v7.3','-append','A','info');
            end
        end
        disp('data alignment done');
    end


%%%%%%% #11 Clustering
    function run_Cluster(varargin)
        disp('spike clustering >>>');        
              
        % check
        assert(length(outfile2.SD)==length(outfile2.A),'SD and A channel number not match');
        
        %%% Remove low activity channels
        % * It must be done because feature selection and clustering both
        % can not work well with few samples.
        sAmt=cellstat(outfile2.A,'size',2);
        
        %%% Clutering each channel.
        SF=cell(info.chAmt,1);
        % Channel SI: cluster identity of each spike in each channel.
        CSI=cell(info.chAmt,1); % * note there may be 0 in CSI - detected as noise, not assigned to any cluster.
        
        for chi=1:info.chAmt
            % spike数量过少的通道直接复制原值（不删除为好，因为shift merge等需要完整的信息），否则进行正常聚类。
            if sAmt(chi)<paras.cluMinSpkThres %通道spike数量过少
                CSI{chi}=ones(sAmt(chi),1); % 标记为1而非0
            else
                Ach=outfile2.A{chi}; 
                SAch=outfile2.SA{chi};
                
                if paras.bDrawForCluster && sAmt(chi)>paras.drawNum
                    flagDoDraw=true;
                else
                    flagDoDraw=false;
                end
                %%% Make draw from samples in case there're too many spikes
                if flagDoDraw
                    drawI=randi(sAmt(chi),paras.drawNum,1);
                    restI=true(sAmt(chi),1); restI(drawI)=false;
                    restAch=Ach(:,restI);
                    Ach=Ach(:,drawI); 
                    SAch=SAch(drawI);
                end
                
                %%% Feature extraction and selection
                % Exclude spike shape outlier
                spknol=~outlier_detect(Ach'); % spike non-ol
                temp=Ach(:,spknol);
                
                % Feature extraction
                temp=spike_feature(temp,'dim',paras.feaDim);
                % 添加Spike amplitude 信息！
                temp=[temp,SAch(spknol)];
                % Normalize
                SF{chi,1}=zscore(temp);
                
                %%% Clustering
                % Exclude SF outlier
                % * 注意如果outlier_detect放在取feature之后，则在进行time-shift
                % merge时应该对outlier同样进行分析并考虑重新加入下一轮feature selection，因为有可能feature上的少数派可能在实际上是位移造成的。
                sfol=outlier_detect(SF{chi,1});
                temp=SF{chi}(~sfol,:);
                
                % Clustering
                tpCSI=spike_cluster(temp,'kmeans');
                
                %%%
                % Make spike morphology templates out of clusters
                if flagDoDraw
                    tplb=reabylb(tpCSI);
                    % find non-0 types
                    idx=find(tplb.types==0,1);
                    if ~isempty(idx)
                        tplb.types(idx)=[];
                        tplb.typeAmt(idx)=[];
                        tplb.ids(idx)=[];
                        tplb.cAmt=tplb.cAmt-1;
                    end
                    if tplb.cAmt>0
                        spktemplate=zeros(info.spklen,tplb.cAmt);
                        for k=1:tplb.cAmt
                            spktemplate(:,k)=mean(Ach(:,tplb.ids{k}),2);
                        end
                    end
                end               
                
                % Fit to template
                if flagDoDraw && tplb.cAmt>0
                    rAmt=size(restAch,2);
                    restCSI=zeros(rAmt,1);
                    for k=1:rAmt
                        % Distance to all templates.
                        D=zeros(tplb.cAmt,1); S=D;
                        for m=1:tplb.cAmt
                            % average difference (as difference)
                            D(m)=mean(abs(restAch(:,k)-spktemplate(:,m)));
                            % amplitude of signals (max amplitude of two curves) -
                            % measure "signal' in SNR
                            tp=[restAch(:,k),spktemplate(:,m)];
                            S(m)=max(max(tp)-min(tp));
                        end
                        
                        % Find Most qualified template.
                        R=D./S; % difference/signal amplitude ratio
                        [tp,idx]=min(R);
                        if tp<1/paras.SNratioThres
                            restCSI(k)=tplb.types(idx);
                        else
                            restCSI(k)=0;
                        end
                    end
                end
                
                %%%                
                
                % SF outlier
                temp=zeros(size(SF{chi},1),1);
                temp(~sfol)=tpCSI;
                
                CSI{chi}=zeros(sAmt(chi),1);
                if flagDoDraw
                    % spike outlier
                    temp2=zeros(length(drawI),1);
                    temp2(spknol)=temp;
                    % Insert the non-outliers to results with original index.                
                    CSI{chi}(drawI)=temp2;
                    CSI{chi}(restI)=restCSI;
                else
                    CSI{chi}(spknol)=temp;
                end
            end            
            fprintf('|');
        end
        fprintf('\n');
        % *至此，通道数(SD length)=CSI length,依然是info.chAmt,没变化。
        assert(isequal(cellstat(outfile2.SD,'length'),cellstat(CSI,'length')),'SD and CSI number not equal');
        
        %%% Combine different channels from the same electrode(raw channel)
        % together - only using the class labels.
        % * 此步骤的另一主要作用在于可视化分类结果：因为在detect and clutering之间的
        % 步骤可能发生过通道数量的变化和spike的删减，之后的spkalignplot()要正确显示
        % 很麻烦，所以生成此数据。   
        newCSI=cell(info.rawchAmt,1);
        reconSD=cell(info.rawchAmt,1);
        for chi=1:info.chAmt % ！注意必须是chAmt而非rawchAmt
            % 首先各个通道的CSI,比如是[1,5]的要压缩为[1,2].
            spkclu=reabylb(CSI{chi});
            nzI=find(spkclu.types>0); %所有非0的类别
            if ~isempty(nzI)                
                for k=1:length(nzI)
                    CSI{chi}(spkclu.ids{nzI(k)})=k;
                end
            end
                    
            % Find electrode(raw channel) index of the new channel
            rawchi=info.chID(chi);
            % If no raw channel data has been put into this new channel, new
            % data is inserted; else new data is added to existing data.
            if isempty(newCSI{rawchi})
                newCSI{rawchi}=CSI{chi};
                reconSD{rawchi}=outfile2.SD{chi};
            else
                % find biggest number of existing neuron id
                nn=max(newCSI{rawchi});
                % except 0, add the nn as base
                I=(CSI{chi}~=0);
                tp=zeros(length(CSI{chi}),1); tp(I)=nn+CSI{chi}(I);
                newCSI{rawchi}=[newCSI{rawchi};tp];
                reconSD{rawchi}=[reconSD{rawchi};outfile2.SD{chi}];
            end
        end
        % Sort reconSD to time order, same order apply to newCSI. (经过此步骤，different cluster mixed together)
        for chi=1:info.rawchAmt
            if ~isempty(newCSI{chi})
                [reconSD{chi},I]=sort(reconSD{chi},'ascend');
                newCSI{chi}=newCSI{chi}(I);
            end
        end        
        
        %%% Update
        CSI=newCSI;
        
        % 生成新的SD（每个分开的神经元算一个通道）
        newSD=cell(0); newchID=zeros(0,1);
        chcount=0;
        for chi=1:info.rawchAmt
            spkclu=reabylb(newCSI{chi});
            % only choose non-noise channel
            nzI=find(spkclu.types>0); 
            for k=1:length(nzI)
                chcount=chcount+1;
                I=spkclu.ids{nzI(k)};
                newSD{chcount}=reconSD{chi}(I);
                newchID(chcount)=chi;
            end
        end
        SD=newSD;        
        % Update channel ID & chAmt
        info.chID=newchID;
        info.chAmt=length(SD);
        
        % Save output
        save(fnC,'-v7.3','SF','info');% ,'-append'
        save(fnC,'-v7.3','-append','info','paras','reconSD','SD','CSI');
        disp('clustering done');
    end


%%%%%%%%%%%%%%% Load variable sub-functions
    function loadoutfile1()
        if ~exist('outfile1','var')
            fprintf('loading raw signal from %s\n',fnF);
            if bDirectMem
                outfile1=load(fnF);
            else
                outfile1=matfile(fnF);
            end            
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
%             case 'spike feature'
%                 N=funcNumTab.spikeFeature;
            case 'spike cluster'
                N=funcNumTab.cluster;
            case 'final'
                N=funcNumTab.final;
            otherwise
                error('invalid run option');
        end        
    end

end