% spike clustering
%   [SI,lb]=spike_cluster(SF,method)
% method: 'kmeans', 'spc', 'graphcomm'
function SI=spike_cluster(SF,method)
bRemoveLargeNoiseCluster=false; % Remove clusters with very large within-cluster variance (!但问题在于，peak很大的类别也会被处理掉）
% 用于cluster合并的参数
bFDRmerge=true;
mergeFDRthres=4; % <<< 经测试，4是比较好的值。 曾经用15，但那个导致分太粗。一般希望在有shift-merge的情况下初始的还是分细些比较好。

% Processing
sAmt=size(SF,1);


%%%%%%%%%%%%% clustering by different methods
switch method
    %%%%%%%%%%%%% Appended K-means
    % This is a method based on kmeans, or you can view as an"shell"
    % overlay on kmeans function.
    case 'kmeans'                
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
%             fprintf('%d',length(rmlist));
            lb.types(rmlist)=[];
            lb.cAmt=length(lb.types);
            lb.typeAmt(rmlist)=[];
            lb.ids(rmlist)=[];
        end
        
        %%% Combine clusters with low separability 
        % - measured by the multi-dimension FDR        
        if bFDRmerge
            while lb.cAmt>1
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
    
    %%%%%%%%%%%%%%% Super-paramagnetic clustering (SPC)
    case 'spc'        
        [SI,res]=spc(SF);
        
    %%%%%%%%%%%%%%%% Based on network community detection
    % This could give number of clusters automatically.
    case 'graphcomm'
        SI=netcluster(SF); 
end

end