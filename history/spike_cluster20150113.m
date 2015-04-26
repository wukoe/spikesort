% spike clustering
%   [SI,lb]=spike_cluster(SF,method,totalTime)
function [SI,lb]=spike_cluster(SF,method)
bRemoveLargeNoiseCluster=false; % Remove clusters with very large within-cluster variance (!但问题在于，peak很大的类别也会被处理掉）
% 用于cluster合并的参数
bFDRmerge=true;
mergeFDRthres=4; % <<< 经测试，15是比较好的值。 4 for NP not separated.

% Processing
sAmt=size(SF,1);

%%%%%%%%%%%%% clustering by different methods
switch method
    case 'kmeans'
        %%%%%%%%%%%%% Appended K-means
        % This is a method based on kmeans, or you can view as an"shell"
        % overlay on kmeans function.        
        cluDistThres=2; % maximum within-cluster variance (fold to standard - which will present later)        
        
        %%% Start from initial number of clusters
        icn=5; % initial cluster number        
        icn=beinrange(icn,size(SF,1),'high'); % restrict it to numbers of samples in total
        
        if icn<=4 % no need to do clustering            
            withinCluVar=0;
            SI=ones(icn,1);
            lb=reabylb(SI);
            return %
        else
            %%% K-means
            % ??? for some unknown reason, the matlab kmean failed
            % sometimes with same data, indicating problem that relate to random
            % initiation. Temporary solution: repetative try
            failtime=0;
            retry=true;
            while retry && failtime<=5
                try                 
                    [SI,C]=kmeans(SF,icn);
                    retry=false;
                catch EM
                    failtime=failtime+1;            
                end
            end
            if retry % if the error flag is still there
                rethrow(EM);
            end
            
            lb=reabylb(SI);

            % Get within-cluster variance - the mean distance to centroid
            % (for bRemoveLargeNoiseCluster option)
            if bRemoveLargeNoiseCluster
                withinCluVar=zeros(icn,1);
                for k=1:icn
                    I=(SI==k);
                    temp=SF(I,:);
                    withinCluVar(k)=mean(juli(temp,C(k,:)));
                end
            end
        end
        
        %%%%%%%%%%%%%% Merge those should be together & Remove "outlier" clusters
        % Remove those with very few samples.         
        rmlist=find(lb.typeAmt<=2);
        % Remove clusters with very large within-cluster variance (!但问题在于，peak很大的类别也会被处理掉）
        if bRemoveLargeNoiseCluster
            % Get the standard for comparison - the group with majority of members
            [~,idx]=max(lb.typeAmt);
            cluDist=withinCluVar(idx); % this is standard of variance: the cluster with most samples
            for k=1:icn            
                if withinCluVar(k)>cluDist*cluDistThres
                    rmlist=[rmlist;k];
                end
            end
        end
        % Remove
        if ~isempty(rmlist)
            fprintf(length(rmlist));
            lb.types(rmlist)=[];
            lb.cAmt=length(lb.types);
            lb.typeAmt(rmlist)=[];
            lb.ids(rmlist)=[];
        end
        
        %%% Combine clusters with low separability 
        % - measured by the multi-dimension FDR        
        if bFDRmerge
%             fNoMoreReduct=false;
            while lb.cAmt>1 % && ~fNoMoreReduct                
                % Calculate FDR between pairs of clusters
                FV=zeros(lb.cAmt);
                for m=1:lb.cAmt-1
                    for n=m+1:lb.cAmt
                        FV(m,n)=abs(fdr({SF(lb.ids{m},:),SF(lb.ids{n},:)}));
                    end
                end
                
                % Find the pair with minimum FDR, check, and merge.
                [tp,idx]=matmin(FV,'triu');                
                if tp<mergeFDRthres
                    % Merge cluster: n->m
                    m=idx(1); n=idx(2);                  
                    lb.types(n)=[];
                    lb.cAmt=lb.cAmt-1;
                    lb.typeAmt(m)=lb.typeAmt(m)+lb.typeAmt(n);
                    lb.typeAmt(n)=[];
                    lb.ids{m}=[lb.ids{m};lb.ids{n}];
                    lb.ids(n)=[];
                else
                    break % else no merger happened, end the process.
                end
            end
        end
        
        % ! No want to do this now, that neurons should be remained, deleting is supposed to be at analysis stage, -just
        % as those low firing channel is remained 
%         %%% Remove those types with under threshould number
%         rmlist=find(lb.typeAmt<cluNumThres);
%         lb.types(rmlist)=[];        
%         lb.typeAmt(rmlist)=[];
%         lb.ids(rmlist)=[];
%         lb.cAmt=length(lb.types);
        
    case 'spc'
        % Super-paramagnetic clustering (SPC)
        [SI,res]=spc(SF);
        lb=reabylb(SI);
end

%%% Sort the clusters in the order of original existing time
for k=1:lb.cAmt
    lb.ids{k}=sort(lb.ids{k},'ascend');
end


%%% Convert the lb back to SI
% SI is initiated as 0 - here 0 denote samples not included in any cluster
SI=zeros(sAmt,1);
for ci=1:lb.cAmt
    SI(lb.ids{ci})=lb.types(ci);
end

end

%%%%%%%% obsolete
% 分析是否存在的一个spike两峰都检出的情况
% 这个不适合在spike detect步骤进行，因为可能会造成和两个独立spike的情况无法区分（*引入deep learning 的必要性）
% 而分类后统计上大量样本的存在使区分这两种状况变得容易。
%%% Check if any two types have a fixed time lag that is small enough 
% Get variance and mean of difference pair of types (note that the lag uses
% the index ID difference instead of time difference! - this will be more accurate)
% %%% Time
%         if n~=m
%             % For each spk in n, find the one in m that is just 1 post upstream
%             % to it. this is the most possible scenario, since it is very
%             % unlikely there is a another spike between 2 peaks of a single spike.
%             existMark=false(lb.typeAmt(n),1);
%             for k=1:lb.typeAmt(n)
%                 [near1up,nearidx]=find(lb.ids{m}==(lb.ids{n}(k)-1),1,'first'); % (lb.ids{n}(k)-1) is one pose upstream
%                 lag=spkTime(lb.ids{m}(nearidx))-spkTime(lb.ids{n}(k));
%                 if ~isempty(near1up) && lag<3 % index condition && lag should < 3ms
%                     existMark(k)=true;
%                 end
%             end
%             
%             if sum(existMark)/lb.typeAmt(n)>0.80  % threshold: at least 95% should follow spk in group m
%                 effectiveSpkID=lb.ids{n}(existMark);
%                 
%                 % Further check if the time difference is stable                
%                 spkTimeDiff=spkTime(effectiveSpkID)-spkTime(effectiveSpkID-1);
%                 if std(spkTimeDiff)/mean(spkTimeDiff)<0.5 % threshold: the coef of variance 
%                     checkCV=true;
%                 else
%                     checkCV=false;
%                 end
%                 
%                 % Further check the shape similarity (use the mean curve of 2 clusters)
%                 cmm=mean(spkData(:,effectiveSpkID-1),2);
%                 cmn=mean(spkData(:,effectiveSpkID),2);
%                 % convolute (cnm only move up in time) to get min distance
%                 convLen=floor(size(spkData,1)/2); % at most half the full length window (15ms/2)
%                 cctemp=zeros(convLen,1);
%                 for k=1:convLen
%                     temp=corrcoef(cmm(k:end),cmn(1:end-k+1));
%                     cctemp(k)=temp(1,2);
%                 end
%                 minDist=1-max(cctemp);
%                 
%                 if minDist<0.1 % threshold: distance
%                     checkSimilar=true;
%                 else
%                     checkSimilar=false;
%                 end
%                 
%                 % If both checked, merge mark +1
%                 if checkCV && checkSimilar
%                     sameSpkMark(m,n)=true;
%                 end
%             end
%             
%         end