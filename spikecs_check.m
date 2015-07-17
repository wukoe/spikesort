% Find the space range of each CS by check raw signal(所有通道中出现aligned
% spike的情况)。如果两个cluster其实是同一个，这个方法可以更准确地确认这一点。
% analyze the info from spikemerge() and spikecs().
% purpose: whether there's false positive channels in the CS seqs. Checking
% on raw data is the main method.
%   [stat,IR1,IR2,ID,IKS]=spikecs_check(X,SD,info)
% Use all spikes in SD. (in this case, SD should be marker spikes only.)
%   spikecs_check(X,SD,info,IW)
% IW is cell same size as SD, IW{chi} determine which of spikes of SD{chi}
% should be used to align each channel of X.
% stat:rnum1,rnum2,ksnum,ksrnum.
function [stat,IR1,IR2,ID,O,A]=spikecs_check(X,SD,info,varargin)
% Setting
spkThr=0.1*info.TimeSpan;
Rthr1=2; % SNR thres: normal standard.
% * thr=2 tested when the amplitude is measured by max(c)-min(c). (c is
% mean curve), not sure how much difference would make to measure
% abs(max(c)).
Rthr2=4; % higher standard.
KSthr=0.15; % < based on data with outlier removed.
% threshold for time jitter (by dt STD)
DTJthr=0.15;%(ms) = 3 points at 20000HZ SR.
% threshold for lagging time.
DLthr=1;%(ms)
twin=[-1,1];

% User
flagIW=false;
flagKS=false;
flagFoClu=false;
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'IW'
                IW=pinfo{parai};
                flagIW=true;
            case 'KS'
                flagKS=true;
            case 'fo clu'
                flagFoClu=true;
            otherwise
                error('unidentified options');
        end
    end
end

% Proc
DTJthr=DTJthr*info.srate/1000;
DLthr=DLthr*info.srate/1000;
xcha=size(X,2);
sda=length(SD);


%%%%%%%%%%%%
IR1=false(xcha,sda); IR2=IR1;
IKS=false(xcha,sda);
ID=false(xcha,sda);
rnum1=zeros(xcha,1); rnum2=rnum1;
ksnum=zeros(xcha,1); ksrnum=ksnum;
mt=cell(sda,1);
for chi=1:sda
    if length(SD{chi})<spkThr
        rnum1(chi)=0;
        ksnum(chi)=0;
        ksrnum(chi)=0;
    else
        if flagIW
            % Check the number of IW channels and X channels if match
            assert(size(IW{chi},2)==xcha,'IW and X channel numbers not match.');
        end
        
        % Align the signal in all X channels to one specific channel of SD.
        tsd=cell(xcha,1);
        for xi=1:xcha
            if flagIW
                tsd{xi}=SD{chi}(IW{chi}(:,xi));
            else
                tsd{xi}=SD{chi};
            end
        end
        [A,~,O]=spike_align(X,tsd,info.srate,'window',twin);
        baseidx=O.preww+1;
        aLen=O.spklen;
%         % Remove outlier to make results more accurate. <<<
%         parfor xi=1:xcha
%             ol=outlier_detect(A{xi}');
%             A{xi}=A{xi}(:,~ol);
%         end
        
        %%% Check A中各个通道的叠加的有无。
        R=ratioPeakNoise(A);
        % ID of channels over exist-nonexist threshold
        IR1(:,chi)=R>Rthr1;
        IR2(:,chi)=R>Rthr2;
        
        %%% Check是否存在较大的时间滞后差。record peak time at the same time.
        pt=zeros(xcha,1);
        for xi=1:xcha
            if IR1(xi,chi)
                % find peak location
                cm=mean(A{xi},2);
                [~,idx]=max(abs(cm)); % <<< only hunt for max peak for EPSP.
                pt(xi)=idx-baseidx;
                if pt(xi)>=DLthr
                    ID(xi,chi)=true;
                end                
            end
        end
        % transform it to time (s) relative to the markers spikes (SD)
        temp=pt/info.srate;
        mt{chi}=temp(IR1(:,chi));        
        
        %%% Check是否存在多类别的情况。
        if flagKS
            ksd=zeros(xcha,1);
            parfor xi=1:xcha
                al=A{xi};
                d=zeros(aLen,1);
                for k=1:aLen
                    d(k)=test_ks(al(k,:));
                end
                ksd(xi)=max(d);
            end
            % ID of channels over uni-model threshold
            IKS(:,chi)=ksd>KSthr;
        end
        
        %%% by sorting.
        if flagFoClu
            IC=zeros(xcha,1);
            for xi=1:xcha
                al=A{xi};
                anum=size(al,2);
                % Feature
                temp=spike_feature(al,'dim',3);
                % Cluster
                temp=zscore(temp);
                CSI=spike_cluster(temp,'kmeans');
                % Count the number 
                tp=reabylb(CSI);
                I=(tp.typeAmt>=1/6*info.TimeSpan);
                IC(xi)=sum(I);
            end
        end        
        
        %%% number statistics.
        rnum1(chi)=sum(IR1(:,chi));
        rnum2(chi)=sum(IR2(:,chi));
        ksnum(chi)=sum(IKS(:,chi));
        ksrnum(chi)=sum(IR1(:,chi)&IKS(:,chi));
        fprintf('|');
    end
end

stat.rnum1=rnum1;
stat.rnum2=rnum2;
% stat.meanTime=mt;

O=struct();
if flagKS
    stat.ksnum=ksnum;
    stat.ksrnum=ksrnum;
    O.IKS=IKS;
end
if flagFoClu
    O.IC=IC;
end
stat.meanTime=mt;
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
    v(chi)=max(std(A{chi},[],2)); % STD of curves. choose max of all points
end    
% ap=max(mc)-min(mc); % amplitude of mean curves
ap=max(abs(mc));

% Ratio
R=zeros(xcha,1);
R(I)=(ap(I)./v(I))';
end