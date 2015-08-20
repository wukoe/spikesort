% Check whether one CS exist at specific location by checking raw signal(所有通道中出现aligned
% spike的情况)。如果两个cluster其实是同一个，这个方法可以更准确地确认这一点。
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
            case 'clu' % 检测align是否存在多类.
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
        
        %%% Check A中各个通道叠加后peak的有无。
        R=abs(pamp(chi))./ssd(chi); %>>>spkSNR(A{chi});
        % ID of channels over exist-nonexist threshold
        IR1(chi)=R>Rthr1;
        IR2(chi)=R>Rthr2;
        
        %%% Check是否存在多类别的情况。
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