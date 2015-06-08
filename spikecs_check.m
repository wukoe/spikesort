% analyze the info from spikemerge() and spikecs().
% purpose: whether there's false positive channels in the CS seqs. Checking
% on raw data is the main method.
%   [stat,IKS1,IKS2,IR]=spikecs_check(X,SD,info)
% SD should be marker spikes only.
function [stat,IR1,IR2,IKS]=spikecs_check(X,SD,info)
% Find the space range of each cluster (- the real "range" by 查看所有通道中出现aligned
% spike的情况。如果两个cluster其实是同一个，这个方法可以更准确地把这个找出来。
%%% Setting
xcha=120;
Rthr1=2; % SNR thres: normal standard. (tested)
Rthr2=4; % higher standard.
KSthr=0.15; % < based on data with outlier removed.

%
sda=length(SD);

IKS=false(xcha,sda);
IR1=IKS; IR2=IKS;
ksnum=zeros(xcha,1);
rnum=ksnum; ksrnum=ksnum;
for chi=1:sda
    if length(SD{chi})<0.1*info.TimeSpan
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
        
        ksd=zeros(xcha,1);
        parfor chm=1:xcha
            al=A{chm};
            d=zeros(aLen,1);
            for k=1:aLen
                d(k)=test_ks(al(k,:));
            end            
            ksd(chm)=max(d);
        end
        
        % ID of channels over uni-model threshold
        IKS(:,chi)=ksd>KSthr;
        
        % ID of channels over exist-nonexist threshold
        IR1(:,chi)=R>Rthr1;
        IR2(:,chi)=R>Rthr2;
        
        rnum(chi)=sum(IR1(:,chi));
        ksnum(chi)=sum(IKS(:,chi));
        ksrnum(chi)=sum(IR1(:,chi)&IKS(:,chi));
        fprintf('|');
    end
end

stat.rnum=rnum;
stat.ksnum=ksnum;
stat.ksrnum=ksrnum;
end % of main


%%%%%%%%%%%%%%%
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