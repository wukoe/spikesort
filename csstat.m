%用于检测conduction cluster 的各种特征。
% summary and results in sm:
% chnum: 出现的channle数量
% chmax: 单个channel参与seq数量的最大值
% ccr1: CC ratio of all spikes (in 1st channel)
% ccrall: CC ratio of all spikes (mean in all channel)
% ndbig: (neighbor distance big) seq中相邻两者的距离过大的数量
function res=csstat(stat,seq,seqT,varargin)
seqAmt=length(seq);
chAmt=length(stat.chfr);

res=struct('seqAmt',seqAmt);

%% path time 
% time of traveling the whole path
PT=cell(seqAmt,1);
res.PTmean=zeros(seqAmt,1);
for seqi=1:seqAmt
    PT{seqi}=zeros(stat.seqcount(seqi),1);
    for k=1:stat.seqcount(seqi)
        tp=seqT{seqi}(k,:);
        PT{seqi}(k)=tp(end)-tp(1);
    end
    res.PTmean(seqi)=mean(PT{seqi});
end
res.PT=PT;

%% path length (Distance of the path of deduced propagation)
if ~isempty(varargin)
    CL=varargin{1};
    res.seqPL=cell(seqAmt,1);
    res.seqPLmax=zeros(seqAmt,1);
    for k=1:seqAmt
        sn=stat.seqlen(k);
        res.seqPL{k}=zeros(sn-1,1);
        for m=1:sn-1
            res.seqPL{k}(m)=electrodedist(CL,seq{k}(m),seq{k}(m+1));
        end
    end    
    
    % Length of whole path
    res.PL=zeros(seqAmt,1);
    for k=1:seqAmt
        res.PL(k)=sum(res.seqPL{k});
    end
    
    % Maximum path neighboring distance
    res.seqPLmax(seqi)=max(res.seqPL{seqi});
end

%% CV of PT
% 每个seq成对元素（不只是相邻的）间,以及全长的time difference的CV
res.totalPTCV=zeros(seqAmt,1);
res.PTCVmax=zeros(seqAmt,1);
res.PTCVmin=zeros(seqAmt,1);
res.PTCVmedian=zeros(seqAmt,1);
for k=1:seqAmt
    TDV=zeros(stat.seqlen(k));
    for m=1:stat.seqlen(k)-1
        for n=m+1:stat.seqlen(k)
            tp=seqT{k}(:,n)-seqT{k}(:,m);
            TDV(m,n)=std(tp)/mean(tp);            
        end
    end
    
    % CV of total path time
    res.totalPTCV(k)=TDV(1,stat.seqlen(k));
    % maximum of all path length
    res.PTCVmax(k)=matmax(TDV);
    res.PTCVmin(k)=matmin(TDV,'triu');
end

%% repeat seq 现象涉及channel数量，以及是否有某seq中ch重复出现的情况
chcount=zeros(chAmt,1); % 每个channel出现在seq中的次数
flagRep=false;
res.doubleChNum=0;
for seqi=1:seqAmt    
    for k=1:stat.seqlen(seqi)
        % participating channel ID
        chi=seq{seqi}(k);
        chcount(chi)=chcount(chi)+1;
        
        if sum(seq{seqi}==chi)>1
            flagRep=true;
        end
    end
end
if flagRep
    res.doubleChNum=res.doubleChNum+1;
end
res.chCount=chcount;


%% CC在参与的各个channel的所有spike中所占的比例。
res.ccratio=cell(seqAmt,1);
for seqi=1:seqAmt
    res.ccratio{seqi}=zeros(stat.seqlen(seqi),1);
    for k=1:stat.seqlen(seqi)
        % participating channel ID
        chi=seq{seqi}(k);
        
        res.ccratio{seqi}(k)=stat.seqcount(seqi)/stat.chspknum(chi);
    end
end

%%
%%%%%%%%%% 总结成单个表征值的部分
sm=struct();
sm.seqCmean=mean(stat.seqcount);
sm.PLmean=mean(res.PL);
sm.PTmean=mean(res.PTmean);
sm.PTmedian=median(res.PTmean);

I=~isnan(res.totalPTCV); 
sm.totalPTCVmean=mean(res.totalPTCV(I));
sm.totalPTCVmedian=median(res.totalPTCV(I));
I=~isnan(res.PTCVmax); sm.PTCVmaxmean=mean(res.PTCVmax(I));
I=~isnan(res.PTCVmin); sm.PTCVminmean=mean(res.PTCVmin(I));

% seq涉及的channel数量
sm.chnum=sum(res.chCount>0);
% 单个ch最多参与了多少个seq
sm.chcmax=max(res.chCount);

% seq中相邻两者的距离过大的数量
if ~isempty(varargin)
    thres=sqrt(10);
    mark=(res.seqPLmax>thres);
    sm.ndbig=sum(mark);
end

% mean of CC ratio of first channel and all channels
temp=zeros(seqAmt,1);
temp2=temp;
for seqi=1:seqAmt
    temp(seqi)=res.ccratio{seqi}(1);
    temp2(seqi)=mean(res.ccratio{seqi});
end
sm.ccr1=mean(temp);
sm.ccrall=mean(temp2);
sm.ccrmax=max(temp2);


res.sm=sm;

end
