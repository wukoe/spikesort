% the whole process of spike sorting
% the basic stratege is to process the data one channel by one channel
%   runstat=spikesort(fileName,varargin)
% fileName does not have postfix
% 'RunStep': default=':'; {'filt','spike detect', 'spike align',
% 'spike feature', 'cluster'}, or use 'resume'
% 'GoOn': default='on'; {'off'}
function spikesort(fileName,varargin)
%%%%%%%%%% Parameter Default
%%% �������
% bFullProc=true;
runOpt=':';
bUserSpecifiedMem=false;

% Process control usage
% run flag - marking which of processings are needed
runFlag=false(11,1); % default: all stop
% function-number table
funcNumTab=struct('start',1,'filt',2,'spikeDetect',4, 'spikeAlign',5, ...
     'cluster',8, 'final',10); % 'spikeFeature',6,

%%% �������?
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
if ischar(runOpt) % ֻ�е���
    runOpt={runOpt};
elseif iscell(runOpt) % cell�ж����Ŀ����ʾҪ���еĲ��� 
else
    error('invalid run step option');
end

sna=length(runOpt);
% ����ֻ��һ��������
if sna==1 
    if strcmp(runOpt{1},':') % run all steps
        runFlag(:)=true;
    else
        runFlag(rs2pgnum(runOpt{1}))=true;
    end

% ���벻ֹһ��������     
else
    for sni=1:sna
        if strcmp(runOpt{sni},':') % ����ͨ���
            if sni==1
                runFlag(1:rs2pgnum(runOpt{sni+1}))=true;
            elseif sni==sna
                runFlag(rs2pgnum(runOpt{sni-1}):rs2pgnum('final'))=true;
            else
                runFlag(rs2pgnum(runOpt{sni-1}):rs2pgnum(runOpt{sni+1}))=true;
            end
        else
            runFlag(rs2pgnum(runOpt{sni}))=true;
            % * ��������Щ�ظ�������Ϊ�˴����࣬��������
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


%%% #1 Initiate the Run status data, Parameter. ����ÿ�����б����
% Get info. struct data.
% * ���Ҫ����filt��spike detect���裬���ԭʼ�ļ��л����Ϣ��������������裬��鿴fnS�ļ��Ƿ���ڣ����ǣ�����м����ֳɵ�info
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
% ������ݶ�Ӧ�ļ�����
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
    % * ����reconSD, �Ͳ���Ҫ'chAssign'ѡ���ˡ�
    
    newSD=reconSD;
    for chi=1:info.rawchAmt
        spkclu=reabylb(CSI{chi});
        % �ų�����ѡ� <<<<<<
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
            
            %%% ����ÿһ���Ƿ�Ӧ�úϲ� - �ҵ���ӽ�����һ�ԡ�
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
            
            % �����κο��Ժϲ��ģ�����д���
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
                        
                        % ��ñ��ƶ�����Ҫ�ƶ��ľ���
                        if paras.bIndivMove %ÿ��spike�������У�cluster�ľ�����Ϊ����������Χ�Ĳο�
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
        disp('filtering done'); % ������ϲ
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
        else % ����matfile��������2D��������ƣ�
            [A,rmlist]=spike_align(outfile1,outfile2.SD,info.srate,'chAssign',info.chID,'window',paras.alignWin,'bSmooth',paras.bAlignSmooth);
        end
        outfile2.A=A;
        % ���������޷����ǵı�Եspike��ͬ������SD,SA,SQ����Ϣ��
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
            % spike�������ٵ�ͨ��ֱ�Ӹ���ԭֵ����ɾ��Ϊ�ã���Ϊshift merge����Ҫ��������Ϣ������������������ࡣ
            if sAmt(chi)<paras.cluMinSpkThres %ͨ��spike��������
                CSI{chi}=ones(sAmt(chi),1); % ���Ϊ1����0
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
                % ���Spike amplitude ��Ϣ��
                temp=[temp,SAch(spknol)];
                % Normalize
                SF{chi,1}=zscore(temp);
                
                %%% Clustering
                % Exclude SF outlier
                % * ע�����outlier_detect����ȡfeature֮�����ڽ���time-shift
                % mergeʱӦ�ö�outlierͬ�����з������������¼�����һ��feature selection����Ϊ�п���feature�ϵ������ɿ�����ʵ������λ����ɵġ�
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
        % *���ˣ�ͨ����(SD length)=CSI length,��Ȼ��info.chAmt,û�仯��
        assert(isequal(cellstat(outfile2.SD,'length'),cellstat(CSI,'length')),'SD and CSI number not equal');
        
        %%% Combine different channels from the same electrode(raw channel)
        % together - only using the class labels.
        % * �˲������һ��Ҫ�������ڿ��ӻ�����������Ϊ��detect and clutering֮���
        % ������ܷ�����ͨ�������ı仯��spike��ɾ����֮���spkalignplot()Ҫ��ȷ��ʾ
        % ���鷳���������ɴ����ݡ�   
        newCSI=cell(info.rawchAmt,1);
        reconSD=cell(info.rawchAmt,1);
        for chi=1:info.chAmt % ��ע�������chAmt����rawchAmt
            % ���ȸ���ͨ����CSI,������[1,5]��Ҫѹ��Ϊ[1,2].
            spkclu=reabylb(CSI{chi});
            nzI=find(spkclu.types>0); %���з�0�����
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
        % Sort reconSD to time order, same order apply to newCSI. (�����˲��裬different cluster mixed together)
        for chi=1:info.rawchAmt
            if ~isempty(newCSI{chi})
                [reconSD{chi},I]=sort(reconSD{chi},'ascend');
                newCSI{chi}=newCSI{chi}(I);
            end
        end        
        
        %%% Update
        CSI=newCSI;
        
        % �����µ�SD��ÿ���ֿ�����Ԫ��һ��ͨ����
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