% analyze the info from spikemerge.m
%   [ksnum,rnum,ksrnum]=spikemerge_see(X,SD,info)
function [ksnum,rnum,ksrnum]=spikemerge_see(X,SD,info)
% Find the space range of each cluster (- the real "range" by 查看所有通道中出现aligned
% spike的情况。如果两个cluster其实是同一个，这个方法可以更准确地把这个找出来。
%%% Setting
xcha=120;
Rthr=2; % <<< partially tested
KSthr=0.15; % < based on data with outlier removal.
 
ksnum=zeros(xcha,1);
rnum=ksnum;
ksrnum=ksnum;
for chi=1:xcha
    if length(SD{chi})<6
        rnum(chi)=0;
        ksnum(chi)=0;
        ksrnum(chi)=0;
    else
        % Align the signal in all X channels to representative SD of cluster.
        sd=cell(xcha,1);
        for k=1:xcha
            sd{k}=SD{chi};
        end
        A=spike_align(X,sd,info.srate,'window',[-1,1]);
        aLen=size(A{1},1);
        
        parfor k=1:xcha
            ol=outlier_detect(A{k}');
            A{k}=A{k}(:,~ol);
        end
        
        % 看A中各个通道的叠加的有无和类别数。
        R=ratioPeakNoise(A);
        
        kst=zeros(xcha,1);
        for chm=1:xcha
            al=A{chm};
            p=zeros(aLen,1);
            parfor k=1:aLen
                p(k)=test_ks(al(k,:));
            end
            
            kst(chm)=max(p);
        end
        
        % ID of channels over uni-model threshold
        IKS=kst>KSthr;
        
        % ID of channels over exist-nonexist threshold
        IR=R>Rthr;
        
        rnum(chi)=sum(IR);
        ksnum(chi)=sum(IKS);
        ksrnum(chi)=sum(IR&IKS);
        fprintf('|');
    end
end
end % of main


function R=ratioPeakNoise(A)
aa=cellstat(A,'size',2);
idx=find(aa>0,1); alen=size(A{idx},1);
xcha=length(A);
mc=zeros(alen,xcha); v=zeros(1,xcha);
I=find(aa>=5);
for k=1:length(I)
    chi=I(k);
    mc(:,chi)=mean(A{chi},2); % mean of all curves.
    v(chi)=max(std(A{chi},[],2)); % standard curves.
end    
ap=max(mc)-min(mc); % amplitude of mean curves
% Ratio
R=zeros(xcha,1);
R(I)=(ap(I)./v(I))';
end

% seqAmt=length(seqlist);
% sd=cell(xcha,1);
% Fch=false(xcha,seqAmt);
% for seqi=1:seqAmt
%     % Align the signal in all X channels to representative SD of cluster.
%     chi=seqlist{seqi}(idxS(seqi));
%     temp=SD{chi}(seqCS{seqi}(:,idxS(seqi)));
%     for chi=1:xcha
%         sd{chi}=temp;
%     end
%     A=spike_align(X,sd,info.srate,'window',[-1,1]);
%     
%     %%% Get ratio of spike peaks to noise
%     R=ratioPeakNoise(A);
%     
%     %%% whether there's one more more population in the morphology
%     % <<<
%     
%     % related channels
%     Fch(:,seqi)=R;     
% end
% I=(Fch>Rthr);