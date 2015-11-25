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
        opt.bRemoveLargeNoiseCluster=false; % Remove clusters with very large within-cluster variance (!���������ڣ�peak�ܴ�����Ҳ�ᱻ�������
        % ����cluster�ϲ��Ĳ���
        opt.mergeFDRthres=4; % <<< �����ԣ�4�ǱȽϺõ�ֵ�� ������15�����Ǹ����·�̫�֡�һ��ϣ������shift-merge������³�ʼ�Ļ��Ƿ�ϸЩ�ȽϺá�
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