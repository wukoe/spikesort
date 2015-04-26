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
    'spikeFeature',6, 'cluster',8, 'final',10); 
% 888 'spikeAlign',6, 

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
paras.alignWin=[-0.8,1.1]; % (ms)
paras.bAlignSmooth=false;

% Spike feature
paras.feaDim=3;

% Spike cluster
paras.cluMinSpkThres=6; % cluster number minimum (/min)
paras.bTimeShiftMerge=false; 
paras.lagWin=[-0.8,0.8]; % shift merging's window (ms)
paras.ccThres=0.8; %<<< cluster mean curve correlation threshold 0.85
paras.SNratioThres=10;

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
        runFlag(rs2pgnum(runOpt{1},funcNumTab))=true;
    end

% ���벻ֹһ��������     
else
    for sni=1:sna
        if strcmp(runOpt{sni},':') % ����ͨ���
            if sni==1
                runFlag(1:rs2pgnum(runOpt{sni+1},funcNumTab))=true;
            elseif sni==sna
                runFlag(rs2pgnum(runOpt{sni-1},funcNumTab):rs2pgnum('final',funcNumTab))=true;
            else
                runFlag(rs2pgnum(runOpt{sni-1},funcNumTab):rs2pgnum(runOpt{sni+1},funcNumTab))=true;
            end
        else
            runFlag(rs2pgnum(runOpt{sni},funcNumTab))=true;
            % * ��������Щ�ظ�������Ϊ�˴����࣬��������
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


%%% #1 Initiate the Run status data, Parameter. ����ÿ�����б����
% Get info. struct data.
% * ���Ҫ����filt��spike detect���裬���ԭʼ�ļ��л����Ϣ��������������裬��鿴fnS�ļ��Ƿ���ڣ����ǣ�����м����ֳɵ�info
if ~runFlag(funcNumTab.filt) && ~runFlag(funcNumTab.spikeDetect) && exist([fnS,'.mat'],'file') % load processed information from _s.mat file
    load(fnS,'info');
else % directely from .mcd file
    info=getMCDinfo([fileName,'.mcd'],bDirectMem);
    info.chAmt=info.rawchAmt;
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
% save(sspath,'info','paras','runFlag'); obsolete

%%% Process of others
paras.lagWin=floor(paras.lagWin/1000*info.srate);


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
    info=outfile1.info;
    if ~exist('outfile2','var')
        if bDirectMem
            outfile2=load(fnS);
        else
            outfile2=matfile(fnS);
        end
        info=outfile2.info;
    end
        
    run_Align();
end
if runFlag(funcNumTab.spikeFeature)
    if ~exist('A','var'), load(fnS,'info'); end
    run_Feature();
end
if runFlag(funcNumTab.cluster)
    if ~exist('SF','var'), load(fnC,'info'); end
    run_Cluster();
end

%%%%%%%%%%%%%% Whether to enter the spike time-shift merge stage
if paras.bTimeShiftMerge
end


% cleaning up
delete(sspath);
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
        rawSD=cell(info.chAmt,1);
        rawSQ=cell(info.chAmt,1);
        rawSA=cell(info.chAmt,1);
        for chi=1:info.rawchAmt
            [rawSD{chi},rawSQ{chi},rawSA{chi}]=spike_detect(outfile1.X(:,chi),info.srate,'movThres',paras.movThres);
            fprintf('|');
        end
        fprintf('\n');
        outfile2.rawSD=rawSD; outfile2.rawSQ=rawSQ; outfile2.rawSA=rawSA;
        
        % Remove low activity channels and save results as SD
        sa=cellstat(rawSD,'length');
        I=sa>=5;
        SD=rawSD(I); SQ=rawSQ(I); SA=rawSA(I);
        info.chAmt=sum(I);
        chID=1:info.rawchAmt; 
        info.chID=chID(I);        
        
        outfile2.SA=SA; % * SQ not need to be saved for further use.
        outfile2.info=info;
        
        if bDirectMem
            save(fnS,'-v7.3','-struct','outfile2');
        end
        disp('detection done');
        
        %%% Splite the positive and negative spikes
        if paras.bPNSep
            disp('split N/P spikes >>>');
            newSD=cell(0);
            newSA=cell(0);
            newchID=zeros(0,1);
            
            chcount=0;
            for chi=1:info.chAmt
                PM=SQ{chi}(:,1); % get posi/nega mark
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
            outfile2.info=info;          
            if bDirectMem
                save(fnS,'-v7.3','-struct','outfile2');
            end
            disp('splited');
        end        
    end


%%%%% Spike alignment is separated here.
    function run_Align()         
        disp('spike alignment >>>');
        if bDirectMem
            [A,rmlist,ST]=spike_align(outfile1.X,outfile2.SD,info.srate,'chAssign',outfile2.info.chID,'window',paras.alignWin,'bSmooth',paras.bAlignSmooth,'T',outfile1.T);
        else % ����matfile��������2D��������ƣ�
            [A,rmlist]=spike_align(outfile1,outfile2.SD,info.srate,'chAssign',outfile2.info.chID,'window',paras.alignWin,'bSmooth',paras.bAlignSmooth);
        end
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
        
        if bDirectMem
            save(fnS,'-v7.3','-append','A','paras','ST');
        end
        disp('data alignment done');
    end


%%%%% #8 Feature extract
    function run_Feature(varargin)
        % Set input
        if ~exist('A','var'), load(fnS,'A','SA'); end
    
        if nargin>0
            flagChSpec=true;
            chSpec=varargin{1};
        else
            flagChSpec=false;
            disp('spike feature >>>');
        end     
        
        if flagChSpec
            if ~exist('SF','var'), load(fnC,'SF'); end
            temp=spike_feature(A{chSpec,1},'dim',paras.feaDim);
            temp=[temp,SA{chSpec}];
            SF{chSpec,1}=temp;
            
            % Save output
            save(fnC,'-v7.3','-append','SF');
        else
            SF=cell(info.chAmt,1);
            rmlist=[];
            for chi=1:info.chAmt
                if isempty(A{chi,1})
                    rmlist=[rmlist,chi];
                    fprintf('X');
                else
                    temp=spike_feature(A{chi,1},'dim',paras.feaDim);
                    % ���Spike amplitude ��Ϣ��
                    temp=[temp,SA{chi}];
                    SF{chi,1}=zscore(temp);
                    fprintf('|');
                end
            end
            fprintf('\n');
            
            if ~isempty(rmlist)
                SF(rmlist)=[];
                info.chID(rmlist)=[];
                info.chAmt=info.chAmt-length(rmlist);                
            end
            
            % Save output
            save(fnC,'-v7.3','SF','info');% ,'-append'
            disp('feature extract done');
        end
    end


%%%%%%% #11 Clustering
    function run_Cluster(varargin)
        disp('spike clustering >>>');        
        % Set up input
        if ~exist('SD','var'), load(fnS,'SD','A'); end 
        if ~exist('SF','var'), load(fnC,'SF'); end
        if ~exist('T','var'), load(fnF,'T'); end
        
        %%% Proc
        totalTime=info.TimeSpan; 
        % cluster number minimum is based on spikes/min, so the actual number
        % should be calculated by time length
        paras.cluMinSpkThres=paras.cluMinSpkThres*(totalTime/60);
        
        % Channel SI: cluster identity of each spike in one channel.
        CSI=cell(info.chAmt,1); % * note there may be 0 in CSI - those not assigned to any cluster.
        % newchID: from which raw channel (electrode) the neuron is recorded.
%         % NSD is neuron spk data (index format)
%         newSD=cell(0,1);
        
        %%% Clutering each channel.
        chcount=0; % number of effective neurons - ��������Ԥ֪�������ۼƼ�����
        for chi=1:info.chAmt
            sAmt=size(SF{chi},1);
            % spike�������ٵ�ͨ��ֱ�Ӹ���ԭֵ����ɾ��Ϊ�ã���Ϊshift merge����Ҫ��������Ϣ������������������ࡣ
            if sAmt<=paras.cluMinSpkThres %ͨ��spike��������
                chcount=chcount+1;
%                 newSD{chcount,1}=SD{chi};
                CSI{chi}=zeros(sAmt,1);
            else
                % Exclude outlier first
                % * ע�����outlier_detect����ȡfeature֮�����ڽ���time-shift
                % mergeʱӦ�ö�outlierͬ�����з������������¼�����һ��feature selection����Ϊ�п���feature�ϵ������ɿ�����ʵ������λ����ɵġ�
                temp=SF{chi,1};
                ol=outlier_detect(temp);
                temp(ol,:)=[];
                
                % Clustering
                tpCSI=spike_cluster(temp,'kmeans');
                                
                % Insert the non-outliers to results with original index.
                CSI{chi}=zeros(sAmt,1);
                CSI{chi}(~ol)=tpCSI;
            end            
            fprintf('|');
        end
        fprintf('\n');
        % *���ˣ�ͨ����(SD length)=CSI length,��Ȼ��info.chAmt,û�仯��        
        
        %%% Combine different channels from the same electrode(raw channel)
        % together - only using the class labels.
        % * �˲������һ��Ҫ�������ڿ��ӻ�����������Ϊ��detect and clutering֮���
        % ������ܷ�����ͨ�������ı仯��spike��ɾ����֮���spkalignplot()Ҫ��ȷ��ʾ
        % ���鷳���������ɴ����ݡ�   
        newCSI=cell(info.rawchAmt,1);
        reconSD=cell(info.rawchAmt,1);
        for chi=1:info.chAmt % ��ע�������chAmt����rawchAmt
            % ���ȸ���ͨ����CSI,������[1,5]��Ҫѹ��Ϊ[1,2].
            lbinfo=reabylb(CSI{chi});
            nzI=find(lbinfo.types>0); %���з�0�����
            if ~isempty(nzI)                
                for k=1:length(nzI)
                    CSI{chi}(lbinfo.ids{nzI(k)})=k;
                end
            end
                    
            % Find electrode(raw channel) index of the new channel
            rawchi=info.chID(chi);
            % If no raw channel data has been put into this new channel, new
            % data is inserted; else new data is added to existing data.
            if isempty(newCSI{rawchi})
                newCSI{rawchi}=CSI{chi};
                reconSD{rawchi}=SD{chi};
            else
                % find biggest number of existing neuron id
                nn=max(newCSI{rawchi});
                % except 0, add the nn as base
                I=(CSI{chi}~=0);
                tp=zeros(length(CSI{chi}),1); tp(I)=nn+CSI{chi}(I);
                newCSI{rawchi}=[newCSI{rawchi};tp];
                reconSD{rawchi}=[reconSD{rawchi};SD{chi}];
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
            lbinfo=reabylb(newCSI{chi});
            % only choose non-noise channel
            nzI=find(lbinfo.types>0); 
            for k=1:length(nzI)
                chcount=chcount+1;
                I=lbinfo.ids{nzI(k)};
                newSD{chcount}=reconSD{chi}(I);
                newchID(chcount)=chi;
            end
        end
        SD=newSD;        
        % Update channel ID & chAmt
        info.chID=newchID;
        info.chAmt=length(SD);
        
        % Save output
        save(fnC,'-v7.3','-append','info','paras','reconSD','SD','CSI');
        disp('clustering done');
    end


%%%%%%%%%%%%%%% One channel goes through all stages.
% <<<������
    function run_channel()
        disp('all stage >>>');
        SD=cell(info.chAmt,1);
        SQ=cell(info.chAmt,1);
        SA=cell(info.chAmt,1);
        SDT=cell(info.chAmt,1);
        
        for chi=1:info.rawchAmt
            [SD{chi},SQ{chi},SA{chi}]=spike_detect(outfile1.X(:,chi),info.srate,'movThres',paras.movThres);
            
            % <<< output is cell, attention.
            [A,rmlist,SDT(chi)]=spike_align(outfile1.X,outfile2.SD,info.srate,'chAssign',outfile2.info.chID,'window',paras.alignWin,'bSmooth',paras.bAlignSmooth,'T',outfile1.T);
            
            if paras.bPNSep
            end
            
            temp=spike_feature(A{chSpec,1},'dim',paras.feaDim);
            temp=[temp,SA{chSpec}];
            SF{chSpec,1}=temp;
            
            sAmt=size(SF{chi},1);
            % spike�������ٵ�ͨ��ֱ�Ӹ���ԭֵ����ɾ��Ϊ�ã���Ϊshift merge����Ҫ��������Ϣ������������������ࡣ
            if sAmt<=paras.cluMinSpkThres %ͨ��spike��������
                chcount=chcount+1;
                newSD{chcount,1}=SD{chi};
                newchID(chcount,1)=info.chID(chi);
                CSI{chi}=zeros(sAmt,1);
            else
                %%% Exclude outlier first
                % * ע�����outlier_detect����ȡfeature֮�����ڽ���time-shift
                % mergeʱӦ�ö�outlierͬ�����з������������¼�����һ��feature selection����Ϊ�п���feature�ϵ������ɿ�����ʵ������λ����ɵġ�
                temp=SF{chi,1};
                ol=outlier_detect(temp);
                temp(ol,:)=[];
                % Clustering
                [tpCSI,lbinfo]=spike_cluster(temp,'kmeans');
                % Insert the outliers back to results with CSI=0
                CSI{chi}=zeros(sAmt,1);
                CSI{chi}(ol)=0;
                CSI{chi}(~ol)=tpCSI;
                idx=find(lbinfo.types==0); % find the noise type
                % If no outlier in original clustering, append this new
                % type; else add to existing noise type.
                if isempty(idx)
                    lbinfo.cAmt=lbinfo.cAmt+1;
                    lbinfo.types(lbinfo.cAmt)=0;
                    lbinfo.typeAmt(lbinfo.cAmt)=sum(ol);
                    lbinfo.ids{lbinfo.cAmt}=find(ol);
                else
                    lbinfo.typeAmt(idx)=lbinfo.typeAmt(idx)+sum(ol);
                    lbinfo.ids{idx}=[lbinfo.ids{idx};find(ol)];
                end
                
                %%%%%%%%%%% Shift the signal and get distance measure
                % * - should be done when all channels of one electrode are
                % put together, since the previous posi/nega separation
                % have already split electrode.
                if paras.bTimeShiftMerge
                    spkData=A{chi,1};
                    %%% Recursive process to combine the pair of smallest
                    while lbinfo.cAmt>1
                        %%% Get mean curve of each cluster
                        cm=zeros(size(spkData,1),lbinfo.cAmt);
                        for k=1:lbinfo.cAmt
                            cm(:,k)=mean(spkData(:,lbinfo.ids{k}),2);
                        end
                        
                        %%% ����ÿһ���Ƿ�Ӧ�úϲ� - �ҵ���ӽ�����һ�ԡ�
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
                        
                        % ��û���κο��Ժϲ��ģ�������һͨ��
                        if sum(sum(sameSpkMark))==0
                            break
                        end
                        
                        %%% Merge the "same spike" pair with highest cc value
                        cc(~sameSpkMark)=0; %<<<<
                        [~,idx]=matmax(cc,'triu');
                        
                        % ��cluster��Ա��������Ϊʱ��������׼
                        num=lbinfo.typeAmt(idx);
                        [~,movidx]=min(num); % movidx�Ǳ��ƶ��ߵģ���idx�еģ���ţ�ֵֻ��1��2����ѡ��
                        movidx=idx(movidx); % ��ɾ���Ҫ�ƶ��ĸ�cluster��
                        movpts=abs(ML(idx(1),idx(2)));
                        
                        % ��ñ��ƶ�����Ҫ�ƶ��ľ��루ÿ��spike�������У�cluster�ľ�����Ϊ����������Χ�Ĳο���
                        if (ML(idx(1),idx(2))<0 && movidx==1) || (ML(idx(1),idx(2))>0 && movidx==2)%���ƶ�����ǰ
                            temp=spkData(:,lbinfo.ids{movidx});
                            spkData(:,lbinfo.ids{movidx})=[temp(movpts+1:end,:); zeros(movpts,lbinfo.typeAmt(movidx))];
                        else
                            temp=spkData(:,lbinfo.ids{movidx});
                            spkData(:,lbinfo.ids{movidx})=[zeros(movpts,lbinfo.typeAmt(movidx)); temp(1:end-movpts,:)];
                        end
                        
                        % �ϲ�
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
                
                fprintf('|');
            end
            fprintf('\n');
        end
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
    function N=rs2pgnum(rsstr,funcNumTab)
        switch rsstr
            case 'start'
                N=funcNumTab.start;
            case 'filt'
                N=funcNumTab.filt;
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

end