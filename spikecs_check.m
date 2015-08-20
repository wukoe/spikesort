% Check whether one CS exist at specific location by checking raw signal(����ͨ���г���aligned
% spike�����)���������cluster��ʵ��ͬһ��������������Ը�׼ȷ��ȷ����һ�㡣
% analyze the info from spikemerge() and spikecs().
% purpose: whether there's false positive channels in the CS seqs. Checking
% on raw data is the desired method.
%   [stat,IR1,IR2,ID,IKS]=spikecs_check(X,SD,info)
% apply each individual member of SD on all X channels.
%   spikecs_check(X,SD,info,IW)
% IW is cell same size as SD, IW{chi} determine which of spikes of SD{chi}
% should be used to align each channel of X. (in case this option is not
% used, SD should be marker spikes only.) 
% stat:rnum1,rnum2,ksnum,ksrnum.
function [stat,IR1,IR2,O]=spikecs_check(A,info,varargin)
% Setting
spkThr=0.1*info.TimeSpan; % spike number threshold.
Rthr1=2.5; % SNR thres: normal standard.
% * thr=2 tested when the amplitude is measured by max(c)-min(c). (c is
% mean curve), not sure how much difference would make to measure
% abs(max(c)).
Rthr2=5; % higher standard.
KSthr=0.15; % < based on data with outlier removed.

% User
flagKS=false;
flagFoClu=false;
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'KS'
                flagKS=true;
            case 'clu' % ���align�Ƿ���ڶ���.
                flagFoClu=true;
            otherwise
                error('unidentified options');
        end
    end
end

% Proc
cha=length(A);

%%%
IR1=false(cha,1); IR2=IR1;
if flagKS, IKS=false(cha,1); end
if flagFoClu, IC=zeros(cha,1); end
stype=zeros(cha,1); ploc=stype; pamp=stype; ssd=stype;
% flagPeakLocChange=false(cha,1);
for chi=1:cha
    spka=size(A{chi},2);
    if spka>=spkThr
        aLen=size(A{chi},1);
        
        % Spike type identification, and adjust peak accordingly.
        [stype(chi),tloc,tamp,tssd]=spktype(A{chi});
        if stype(chi)==0
            pamp(chi)=0;
        else
            ploc(chi)=tloc; pamp(chi)=tamp; ssd(chi)=tssd;
        end
        
        %%% Check A�и���ͨ�����Ӻ�peak�����ޡ�
        R=abs(pamp(chi))./ssd(chi); %>>>spkSNR(A{chi});
        % ID of channels over exist-nonexist threshold
        IR1(chi)=R>Rthr1;
        IR2(chi)=R>Rthr2;
        
        %%% Check�Ƿ���ڶ����������
        if flagKS
            d=zeros(aLen,1);
            for k=1:aLen
                d(k)=test_ks(A{chi}(k,:));
            end
            ksd=max(d);
            % ID of channels over uni-model threshold
            IKS(chi)=ksd>KSthr;
        end
        
        %%% by sorting.
        if flagFoClu
            % Feature
            temp=spike_feature(A{chi},'dim',3);
            % Cluster
            temp=zscore(temp);
            CSI=spike_cluster(temp,'kmeans');
            % Count the number
            tp=reabylb(CSI);
            I=(tp.typeAmt>=1/6*info.TimeSpan);
            IC(chi)=sum(I);
        end        
    end
end

% Number statistics and Output.
stat=struct('stype',stype,'ploc',ploc,'pamp',pamp,'ssd',ssd);

stat.rnum1=sum(IR1);
stat.rnum2=sum(IR2);
O=struct();
if flagKS
    stat.ksnum=sum(IKS);
    stat.ksrnum=sum(IR1&IKS);
    O.IKS=IKS;
end
if flagFoClu
    O.IC=IC;
end

end % of main