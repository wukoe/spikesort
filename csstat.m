%���ڼ��conduction cluster �ĸ���������
%   res=csstat(seqlist,seqMark,seqType,seqCS,stat,chID,CL)
% res contains
% 
function res=csstat(seqlist,seqMark,seqType,seqCS,stat,chID,CL)
ied=100; % inter-electrode distance.(um)

%%% pre-Proc.
% Number of CS and channel number in each CS
seqAmt=length(seqlist);
seqlen=cellstat(seqlist,'length'); % length of whole sequence
axonlen=zeros(seqAmt,1); % length of axon (suppose spike with negative peak)
axonseq=cell(seqAmt,1);
for seqi=1:seqAmt
    axonlen(seqi)=sum(seqType{seqi}==-1);
    axonseq{seqi}=seqlist{seqi}(seqType{seqi}==-1);
end

res=struct('seqAmt',seqAmt);
res.seqlen=seqlen; res.axonlen=axonlen;

%%% * �ڽ�������path time, path time���������²��Լ����ٶȣ�ÿ��ͨ�����ҵ������������һ��ͨ����
% ֻ�������߼��ʱ���;��룬�Ӷ��ó��ٶȡ�
% channel -> electrode location
tss=cell(seqAmt,1);
for seqi=1:seqAmt
    tss{seqi}=chID(axonseq{seqi});
end
% Path Distance.
[PL,pairI]=csdist(tss,CL,3);
% unit to (m)
PL=PL*ied/1e6;
res.PL=PL;

% Path Time (sum of between nearest pairs)
PT=zeros(seqAmt,1); 
for seqi=1:seqAmt
    if axonlen(seqi)<2
        PT(seqi)=NaN;
    else
        % restrict to axons.
        cstime=stat.meanTime{seqi}(seqType{seqi}==-1);
        segtime=zeros(axonlen(seqi),1);
        tp=cstime(pairI{seqi});
        for k=1:axonlen(seqi)
            segtime(k)=abs(tp(k,1)-tp(k,2));
        end
        PT(seqi)=sum(segtime);
    end
end
res.PT=PT; 

%%% Speed of conduction
I=isnan(PT);
temp=zeros(seqAmt,1);
temp(~I)=PL(~I)./PT(~I);
temp(I)=NaN;
res.apspeed=temp;

%%% CS�ڲ���ĸ���channel������spike����ռ�ı�����
% first marker's idx in each sequence.
markI=zeros(seqAmt,1);
for seqi=1:seqAmt
    markI(seqi)=find(seqlist{seqi}==seqMark(seqi));
end

cspMarker=cell(seqAmt,1);
for seqi=1:seqAmt
    tp=seqCS{seqi}{markI(seqi)};
    csnum=cellstat(seqCS{seqi},'sum');
    cspMarker{seqi}=csnum'/sum(tp);
end
res.cspMarker=cspMarker;

% % seq���������ߵľ�����������
% thres=sqrt(10);
% mark=(res.seqPLmax>thres);
% sm.ndbig=sum(mark);

end
