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
paras.alignWin=[-0.6,0.7]; % [-0.8 1.1] (ms)
paras.bAlignSmooth=false;
paras.bExactST=false;

% Spike cluster
paras.cluMinSpkThres=12; % cluster number minimum (/min)
paras.feaDim=3;% Spike feature
paras.bDrawForCluster=true;
paras.drawNum=2000;
paras.bMatchNoisySpk=true;
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
%%%%%% #2 Filtering and removal of bad channels
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
if runFlag(funcNumTab.spikeDetect)
    % Set up input - if no filtered data (outfile1) in memory, try load it.
    loadoutfile1();
    % Set up output
    if bDirectMem
        outfile2=struct();
    else
        outfile2=matfile(fnS,'Writable',true);
    end
    
    
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
    
    outfile2.paras=paras;    
    if bDirectMem
        save(fnS,'-v7.3','-struct','outfile2');
    end
    disp('detection done');
end


%%%%% Spike alignment is separated here.
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
    
    disp('spike alignment >>>');
    %%% Remove low activity channels and save results as SD
    sa=cellstat(outfile2.rawSD,'length');
    I=(sa>=5);
    SD=outfile2.rawSD(I); SA=outfile2.rawSA(I); % * SQ not need to be saved for further use.
    info=outfile2.info;
    info.chAmt=sum(I);
    temp=(1:info.rawchAmt)'; info.chID=temp(I);
    
    %%% Spike alignment
    [A,rmlist,SO]=spike_align(outfile1.X,SD,info.srate,'chAssign',info.chID,'window',paras.alignWin,'bSmooth',paras.bAlignSmooth);
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
    if bDirectMem
        save(fnS,'-v7.3','-struct','outfile2');
        % Get Precise time by splining spike peaks if required.
        if paras.bExactST
            ST=exactST(outfile1.X,outfile2.SD,outfile1.T,info.srate,info.chID);
            save(fnS,'-v7.3','-append','ST');
        end
    end
    % a final check
    assert(isequal(cellstat(SD,'length'),cellstat(SA,'length')),'SD and SA number not equal');
    assert(isequal(cellstat(SD,'length'),cellstat(A,'size',2)),'SD and A number not equal');
    disp('data alignment done');
end

if runFlag(funcNumTab.cluster)
    % Set up input
    loadoutfile1();
    if ~exist('outfile2','var')
        outfile2=load(fnS);        
    end
    info=outfile2.info;
        
%     run_Cluster();
    % 2/N
    CST=cell(info.rawchAmt,1);
    reconSD=CST;    reconA=CST;
    for k=1:info.chAmt
        fprintf('|');
        chi=info.chID(k);
        [CST{chi},reconSD{chi},reconA{chi}]=chcluster(outfile1.X(:,chi),outfile2.SD{k},outfile2.SA{k},outfile2.A{k},info,paras);
    end
    fprintf('\n');
    
    %%% 生成新的NSD（每个分开的"神经元"算一个通道）and new chID and chAmt.
    assert(isequal(cellstat(reconSD,'length'),cellstat(CST,'length')),'reconSD and CSI length not match after clustering');
    NSD=cell(0,1);
    chID=zeros(0,1);
    chcount=0;
    for chi=1:info.rawchAmt
        spkclu=reabylb(CST{chi});
        % only choose non-noise channel
        nzI=find(spkclu.types>0);  na=length(nzI);
        if na>0
            temp=cell(na,1);
            for k=1:length(nzI)
                I=spkclu.ids{nzI(k)};
                temp{k}=reconSD{chi}(I);
            end
            NSD(chcount+1:chcount+na,1)=temp;
            chID(chcount+1:chcount+na,1)=chi;
            chcount=chcount+na;
        end
    end
    info.chID=chID;
    info.chAmt=chcount;
    
    %%% save
%     save(fnS,'-v7.3','-append','SD');
    save(fnC,'CST','NSD','info','reconSD','reconA','paras');
end

% Cleaning up
% delete(sspath);
disp('run finish, all done');
%%%%%%%%%%%% End of main.


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

% obsolete %%%%%%%%%%%%%% Whether to enter the spike time-shift merge stage
% if runFlag(funcNumTab.cluster) && paras.bTimeShiftMerge
%     loadoutfile1();
%     load(fnC,'info','reconSD','CSI');
%     % Update to get the new alignment data
%     A=spike_align(outfile1.X,reconSD,info.srate,'window',paras.alignWin,'bSmooth',paras.bAlignSmooth);
%     % * 用了reconSD, 就不需要'chAssign'选项了。
%     
%     newSD=reconSD;
%     for chi=1:info.rawchAmt
%         spkclu=reabylb(CSI{chi});
%         % 排除噪声选项？ <<<<<<
%         idx=find(spkclu.types==0);
%         if ~isempty(idx)
%             spkclu.types(idx)=[];
%             spkclu.typeAmt(idx)=[];
%             spkclu.cAmt=spkclu.cAmt-1;
%             spkclu.ids(idx)=[];
%         end
%         
%         %%% Recursive process to combine the pair of smallest distance
%         if spkclu.cAmt>1
%             spkData=A{chi};
%             %%% Get mean curve of each cluster
%             cm=zeros(info.spklen,spkclu.cAmt);
%             for k=1:spkclu.cAmt
%                 cm(:,k)=mean(spkData(:,spkclu.ids{k}),2);
%             end
%             
%             %%% 计算每一对是否应该合并 - 找到最接近的那一对。
%             sameSpkMark=false(spkclu.cAmt); % mark whether 2 clusters (nominated by row and column) should merge.
%             cc=zeros(spkclu.cAmt); % corr coef
%             D=cc;
%             ML=zeros(spkclu.cAmt); % lag value
%             for m=1:spkclu.cAmt-1
%                 for n=m+1:spkclu.cAmt
%                     %%% Check the whether 2 shapes should merge (use the mean curve of 2 clusters)
%                     % First Move along time to get max CC and related distance
%                     [tp1,tp2]=lagcorr([cm(:,m),cm(:,n)],paras.lagWin);
%                     cc(m,n)=tp1(1,2); ML(m,n)=tp2(1,2);
%                     if cc(m,n)>paras.ccThres                        
%                         % Also check absolute difference between curves (* this requirement
%                         % should also be fullfilled, because 2 types with similar shape but different
%                         % amplitude is also different).                        
%                         % Align 2 curves according to ML.
%                         if ML(m,n)>0 % m behind n
%                             tp1=cm(ML(m,n)+1:end,m); tp2=cm(1:end-ML(m,n),n);
%                         else
%                             tp1=cm(1:end+ML(m,n),m); tp2=cm(1-ML(m,n):end,n);
%                         end                        
%                         % average difference (of time)
%                         D(m,n)=mean(abs(tp1-tp2));
%                         % amplitude of signals (max amplitude of two curves) -
%                         % measure "signal' in SNR
%                         tp=[tp1,tp2];
%                         S=max(max(tp)-min(tp));
%                         
%                         if D(m,n)<=S/paras.SNratioThres
%                             sameSpkMark(m,n)=true;
%                         end
%                     end
%                 end
%             end
%             ML=ML-ML'; % fill the bottom triangle of lag matrix
%             
%             % 若有任何可以合并的，则进行处理
%             if sum(sum(sameSpkMark))>0
%                 % Get the direction of merging (which cluster merged into
%                 % which)  * method: cluster with more member wins. (others
%                 % to consider: less within-cluster variance wins)
%                 MD=false(spkclu.cAmt);
%                 for m=1:spkclu.cAmt-1
%                     for n=m+1:spkclu.cAmt
%                         if sameSpkMark(m,n)
%                             if spkclu.typeAmt(m)>spkclu.typeAmt(n)
%                                 MD(n,m)=1; % n->m
%                             else
%                                 MD(m,n)=1;
%                             end
%                         end
%                     end
%                 end
%                 
%                 % Collapse single source (e.g., if 1->2, 2->3, then get
%                 % 1->3 instead of 1->2).
%                 for m=1:spkclu.cAmt
%                     idx=find(MD(m,:));
%                     while ~isempty(idx)
%                         if length(idx)>1
%                             disp('warning for multi-target'); 
%                             %idx=idx(1); 
%                             % max number/best fit
%                             temp=D(m,idx);
%                             [~,tp]=min(temp);
%                             idx=idx(tp);
%                         end
%                         % move to new target
%                         MD(m,:)=0; MD(m,idx)=1; 
%                         idx=find(MD(idx,:));
%                     end
%                 end
%                 
%                 %%% Merge the "same spike" pair
%                 for m=1:spkclu.cAmt
%                     movidx=find(MD(m,:));
%                     if ~isempty(movidx)
%                         movpts=ML(m,movidx); %movpts P/N property decide whether the data move forward/backward.
%                         % Get more accurate move distance by matching the
%                         % peaks. <<<<<<<<<<<
%                         peakloc=info.spkprew+1; % location of peak of spike data
%                         tp=peakloc+movpts;
%                         [~,idx]=max(abs(cm(tp-2:tp+2,m)));
%                         if idx ~= 3 % which match the number 4 used
%                             movpts=movpts+(idx-3);
%                         end                            
%                         
%                         % 获得被移动者需要移动的距离
%                         if paras.bIndivMove %每个spike单独进行，cluster的距离作为限制搜索范围的参考
%                             for k=1:spkclu.typeAmt(m)
%                                 % Decide the move length of this spike
%                                 % 1/N use max correlation location
%                                 sp=spkData(:,spkclu.ids{m}(k));
%                                 [tp,mp]=lagcorr([cm(:,movidx),sp],[-movpts-2,-movpts+2]);
%                                 tp=tp(1,2); %mp=mp(1,2);
%                                 % 2/N make peaks match.
%                                 
%                                 % e/N
%                                 
%                                 if true % tp>=paras.ccThres
%                                 % Make re-location
%                                 temp=newSD{chi}(spkclu.ids{m}(k));
%                                 temp=temp-mp; % * note this is different from batch method below
%                                 newSD{chi}(spkclu.ids{m}(k))=temp;
%                                 end
%                             end
%                             
%                         else % move in a batch.
%                             newSD{chi}(spkclu.ids{m})=newSD{chi}(spkclu.ids{m})+movpts;
%                         end
% %                         CSI{chi}(spkclu.ids{m})=spkclu.types(movidx);                        
%                     end
%                 end
%             end
%             
%         end
%     end
%     
%     %%% Update those unsorted information as well.
%     SD=newSD;
%     info.chID=1:info.rawchAmt;
%     SA=cell(info.rawchAmt,1);
%     for chi=1:info.rawchAmt
%         SA{chi}=outfile1.X(SD{chi},chi);
%     end
%     save(fnS, '-append','info','SD', 'SA');
%     
%     %%% Do clustering again
%     outfile2=load(fnS);
%     info=outfile2.info;
%     run_Align();
%     info=outfile2.info;
%     run_Cluster();
% end