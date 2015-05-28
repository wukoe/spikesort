% analyze the info from spikemerge() and spikecs().
% purpose: whether there's false positive channels in the CS seqs. Checking
% on raw data is the main method.
%   [stat,IKS1,IKS2,IR]=spikecs_check(X,SD,info)
function [stat,IKS1,IKS2,IR]=spikecs_check(X,SD,info)
% Find the space range of each cluster (- the real "range" by �鿴����ͨ���г���aligned
% spike��������������cluster��ʵ��ͬһ��������������Ը�׼ȷ�ذ�����ҳ�����
%%% Setting
xcha=120;
Rthr=2; % SNR thres <<< partially tested
KSthr1=0.2; % < based on data with outlier removed.
KSthr2=0.15;

sda=length(SD);

IKS1=false(xcha,sda);
IKS2=IKS1; IR=IKS1;
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
        
        % ��A�и���ͨ���ĵ��ӵ����޺��������
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
        IKS1(:,chi)=ksd>KSthr1;
        IKS2(:,chi)=ksd>KSthr2;
        
        % ID of channels over exist-nonexist threshold
        IR(:,chi)=R>Rthr;
        
        rnum(chi)=sum(IR(:,chi));
        ksnum(chi)=sum(IKS1(:,chi));
        ksrnum(chi)=sum(IR(:,chi)&IKS1(:,chi));
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