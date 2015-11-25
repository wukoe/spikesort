% spike clustering
%   [SI,lb]=spike_cluster(SF,method)
% method: 'kmeans', 'spc', 'graphcomm','optics'
function SI=spike_cluster(SF,method)

switch method
    %%% Appended K-means
    % This is a method based on kmeans, or you can view as an"shell"
    % overlay on kmeans function.
    case 'kmeans'
        opt=struct();
        opt.bRemoveLargeNoiseCluster=false; % Remove clusters with very large within-cluster variance (!但问题在于，peak很大的类别也会被处理掉）
        % 用于cluster合并的参数
        opt.mergeFDRthres=4; % <<< 经测试，4是比较好的值。 曾经用15，但那个导致分太粗。一般希望在有shift-merge的情况下初始的还是分细些比较好。
        opt.cluDistThres=2; % maximum within-cluster variance (fold to standard - which will present later)
        opt.icn=5; % initial cluster number
        opt.icn=beinrange(opt.icn,size(SF,1),'high'); % restrict it to numbers of samples in total
        
        if opt.icn<=4 % no need to do clustering            
            SI=ones(opt.icn,1);
            return %
        else
            SI=kmmerge(SF,opt);
        end
    
    %%% Super-paramagnetic clustering (SPC)
    case 'spc'        
        [SI,res]=spc(SF);
        
    %%% Based on network community detection
    % This could give number of clusters automatically.
    case 'graphcomm'
        SI=netcluster(SF,2); 
        
    case 'pointks'
        SI=pointKS(SF',false); % bSpk=false
        
    %%% OPTICS-based method
    case 'optics'
        minpts=10;
        peakwin=20;
        SI=cluster_optics(SF,minpts,peakwin);
end

end