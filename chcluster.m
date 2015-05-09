% process the several steps of spikesort. (align, split, clustering,...)
%   [CSI,SD,AL]=chcluster(X,SD,AL,info,paras)
function [CSI,SD,AL]=chcluster(X,SD,AL,info,paras) % make back later.
roundnum=4;
rawSD=SD;
for round=1:roundnum
    SA=X(SD);
    chAmt=1;    
%%% Splite the positive and negative spikes
if paras.bPNSep
    newSD=cell(0);
    newSA=cell(0);
    newAL=cell(0);
    
    chcount=0;
    PM=(SA>0); % get posi/nega mark
    pnum=sum(PM); % number of positive peaks
    % The positive peaks.
    if pnum>0
        chcount=chcount+1;
        newSD{chcount,1}=SD(PM);
        newSA{chcount,1}=SA(PM);
        newAL{chcount,1}=AL(:,PM);
    end
    % The negative peaks.
    PM=~PM; % reverse index
    if pnum<length(SD) % number of negative peaks > 0
        chcount=chcount+1;        
        newSD{chcount,1}=SD(PM);
        newSA{chcount,1}=SA(PM);
        newAL{chcount,1}=AL(:,PM);
    end
    
    SA=newSA; 
    chAmt=chcount;
else
    SA={SA};
    newAL={AL};
end

% %%% temporary for debug (only for the SO.spklen).
% [AL,rmlist,SO]=spike_align(X,SD,info.srate,'chAssign',ones(chAmt,1),'window',paras.alignWin);
% % Get spike morphology data length.
% info.spklen=SO.spklen; 


%%%%%%%%%%%%%%% Clutering each channel.
% if chAmt <<<
sAmt=cellstat(SA,'length');

% Channel SI: cluster identity of each spike in each channel.
CSI=cell(chAmt,1); % * note there may be 0 in CSI - detected as noise, not assigned to any cluster.
parfor chi=1:chAmt
    Ach=newAL{chi};
    SAch=SA{chi};

    %%% Make draw from samples in case there're too many spikes.
    % Decide whether to make draw based on both option and spike number
    if paras.bDrawForCluster && sAmt(chi)>paras.drawNum
        flagDoDraw=true;
    else
        flagDoDraw=false;
    end
    % Make draw - spike separated to do clustering + left.
    if flagDoDraw
        temp=randperm(sAmt(chi)); % ! Do not use randi(), it's sampling with replacement
        temp=temp(1:paras.drawNum);
        drawI=false(sAmt(chi),1); drawI(temp)=true;
%         leftAch=Ach(:,restI);
        Ach=Ach(:,drawI);
        SAch=SAch(drawI);
    end

    %%% Exclude spike shape outlier
    if paras.bSpkShapeOL
        spknol=~outlier_detect(Ach',1.2);% spike non-outlier
    else
        spknol=true(size(Ach,2),1);
    end
    
    % "clean" spike数量过少的通道直接复制原值（不删除为好，因为shift merge等需要完整的信息），否则进行正常聚类。
    % * It must be done because feature selection and clustering both can not work well with few samples.
    if sum(spknol)<paras.cluMinSpkThres %30 
        CSI{chi}=zeros(size(Ach,2),1); 
        CSI{chi}(spknol)=1; % outlier spike标记为0
    else        
        %%% Feature extraction and selection        
        temp=Ach(:,spknol);        
        % Feature extraction
        temp=spike_feature(temp,'dim',paras.feaDim);
        % 添加Spike amplitude 信息！
        temp=[temp,SAch(spknol)];
        % Normalize
        sf=zscore(temp);
        
        %%% Clustering
        % Exclude SF outlier
        % * 注意如果outlier_detect放在取feature之后，则在进行time-shift
        % merge时应该对outlier同样进行分析并考虑重新加入下一轮feature selection，因为有可能feature上的少数派可能在实际上是位移造成的。
        if paras.bSpkFeatureOL
            sfnol=~outlier_detect(sf,1.1); % SF non-outlier
        else
            sfnol=true(size(sf,1),1);
        end
        temp=sf(sfnol,:);        
        % Clustering
        tpCSI=spike_cluster(temp,paras.cluMethod);
        
        % * 下面这个部分现在的策略是将其置为0，后面和噪声spike共同决定归属。 
        % 另外一个方法是：维持此处独立的template matching, 但把前面抽取的步骤放到outlier detech之后。

        %%% Add the cluster back to data along outliers.        
        % Put back to temp alongside SF outlier
        temp=zeros(size(sf,1),1);
        temp(sfnol)=tpCSI;
                
        CSI{chi}=zeros(sAmt(chi),1);
        if flagDoDraw
            % spike outlier
            temp2=zeros(paras.drawNum,1);
            temp2(spknol)=temp;
            % Insert the non-outliers to results with original index.
            CSI{chi}(drawI)=temp2;
        else
            CSI{chi}(spknol)=temp;
        end
    end
end

%%% Combine positive and negative channels together as one - reorganize the class labels.
if chAmt>1
    % 首先各个通道的CSI要最小化,比如是[1,5]的要压缩为[1,2].
    for chi=1:chAmt
        spkclu=reabylb(CSI{chi});
        nzI=find(spkclu.types>0); %所有非0的类别
        if ~isempty(nzI)
            for k=1:length(nzI)
                CSI{chi}(spkclu.ids{nzI(k)})=k;
            end
        end
    end    
    
    % find biggest number of existing neuron id
    nn=max(CSI{1});
    % except 0, add the nn as base
    I=(CSI{2}~=0);
    CSI{2}(I)=nn+CSI{2}(I);
    newCSI=[CSI{1}; CSI{2}];
    
    % Sort CSI to time order (经过此步骤，different cluster mixed together).
    newSD=[newSD{1}; newSD{2}];
    [~,I]=sort(newSD,'ascend');
    CSI=newCSI(I);
else
    CSI=CSI{1};
end

%%% Match 1. noisy spikes & 2. left aside in spike draw; to the templates (choice 1, the other is at end of whole program) <<<
% Separate noise and clean spikes.
spkclu=reabylb(CSI);
nidx=find(spkclu.types==0,1);
if ~isempty(nidx)
    % Record noisy spike indexes.
    noiseI=spkclu.ids{nidx};
    % 排除噪声选项
    spkclu.types(nidx)=[];
    spkclu.typeAmt(nidx)=[];
    spkclu.cAmt=spkclu.cAmt-1;
    spkclu.ids(nidx)=[];
end

if spkclu.cAmt>0
    % Get mean curve (template) of each "clean" cluster.
    cm=zeros(info.spklen,spkclu.cAmt);
    for k=1:spkclu.cAmt
        cm(:,k)=mean(AL(:,spkclu.ids{k}),2);
    end
    
    if paras.bMatchNoisySpk && ~isempty(nidx)
        % 获得每个noisy Spike的belonging cluster    
        na=length(noiseI);
        if na>0 % spkclu.cAmt>1代表还存在非noisy的模板可供选择。
            for ni=1:na
                sp=AL(:,noiseI(ni));
                % Distance to each template
                D=zeros(spkclu.cAmt,1); S=D;
                for k=1:spkclu.cAmt
                    D(k)=mean(abs(cm(:,k)-sp));
                    tp=[cm(:,k),sp];
                    S(k)=max(max(tp)-min(tp));
                end            
                % Assign cluster if fullfil the criteria.
                R=D./S; % difference/signal amplitude ratio
                [tp,nidx]=min(R);
                if tp<1/paras.SNratioThres_noise
                    CSI(noiseI(ni))=spkclu.types(nidx);
                end
            end
        end
    end
end


%%%%%%%%%%%%%% Whether to enter the spike time-shift merge stage
if paras.bTimeShiftMerge 
    newSD=SD;
    spkclu=reabylb(CSI);
    
    % 排除噪声选项
    idx=find(spkclu.types==0,1);
    if ~isempty(idx)
        spkclu.types(idx)=[];
        spkclu.typeAmt(idx)=[];
        spkclu.cAmt=spkclu.cAmt-1;
        spkclu.ids(idx)=[];
    end

    %%% Recursive process to combine the pair of smallest distance
    if spkclu.cAmt>1
        % * Get mean curve of each cluster (<<<不要在这里重新计算cm,因为之前已经把部分噪声spike加进来了）

        %%% 根据cm计算每一对是否应该合并 - 找到最接近的那一对。
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
                    
                    if D(m,n)/S<=1/paras.SNratioThres_cm
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

            %%% Merge the pair supposed to be "same spike"
            for m=1:spkclu.cAmt
                movidx=find(MD(m,:));
                if ~isempty(movidx)
                    movpts=ML(m,movidx); %movpts P/N property decide whether the data move forward/backward.
%                     % 之前尝试过Get more accurate move distance by matching the
%                     % peaks.(结果：No! matching peaks 不是最好的方法，因为即使是mean curve，peak往往也有对齐问题。)

                    % 被移动者移动
                    if paras.bIndivMove %每个spike单独进行，cluster的距离作为限制搜索范围的参考
                        for k=1:spkclu.typeAmt(m)
                            % Calculate the difference
                            sp=AL(:,spkclu.ids{m}(k));
                            
                            if movpts>0 % m behind movidx??? <<<
                                tp1=sp(movpts+1:end,m); tp2=cm(1:end-movpts,movidx);
                            else
                                tp1=sp(1:end+movpts,m); tp2=cm(1-movpts:end,movidx);
                            end
                            % average difference (of time)
                            D=mean(abs(tp1-tp2));
                            % amplitude of signals (max amplitude of two curves) -
                            % measure "signal' in SNR
                            tp=[tp1,tp2];
                            S=max(max(tp)-min(tp));
                             
                            if D<=S/paras.SNratioThres                                
                                % Make re-location
                                temp=newSD(spkclu.ids{m}(k));
                                temp=temp+movpts; % * note this is different from batch method below
                                newSD(spkclu.ids{m}(k))=temp;
                            end
                        end

                    else % 整个类一起移动，移动同样的距离.
                        newSD(spkclu.ids{m})=newSD(spkclu.ids{m})+movpts;
                    end
                    
                    % UPdate CSI
                    CSI(spkclu.ids{m})=spkclu.types(movidx);
                end
            end
        end
    end

    %%% Update those unsorted spike timing information as well.
    SD=newSD;
    if round<roundnum % >>> chi<chAmt ??
        % Update to get the new alignment data
        AL=spike_align(X,{SD},info.srate,'window',paras.alignWin,'bSmooth',paras.bAlignSmooth); % * 用了reconSD, 就不需要'chAssign'选项了。
        AL=AL{1};
    end
end


% * --> Go loop back to do clustering again 
% End the loop if: 1. only 1 type left; 2. no relocation requested (case of
% merging on the same location should be taken into account for re-doing
% the process).
if spkclu.cAmt<2 || ~paras.bTimeShiftMerge || sum(abs(ML(sameSpkMark)))==0
    break
end
end
fprintf('%d',round);
assert(length(SD)==length(rawSD),'SD spike number changed');

%%% "Squeeze" the CSI number.
spkclu=reabylb(CSI);
nzI=find(spkclu.types>0); %所有非0的类别
if ~isempty(nzI)
    for k=1:length(nzI)
        CSI(spkclu.ids{nzI(k)})=k;
    end
end

%%%%%%%%%%%%%%% Match noisy spikes to the templates (choice 2)
end % func end