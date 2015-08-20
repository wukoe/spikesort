%用于检测conduction cluster 的各种特征。
%   res=csstat(seqlist,seqMark,seqType,seqCS,stat,chID,CL)
% res contains
% 
function res=csstat(seqlist,seqMark,seqType,seqCS,stat,chID,CL)
ied=100; % inter-electrode distance.(um)


%%% Number of CS and channel number in each CS
seqAmt=length(seqlist);
seqlen=cellstat(seqlist,'length'); % length of whole sequence
axonlen=cell(seqAmt,1); % length of axon (suppose spike with negative peak)
axonseq=cell(seqAmt,1);
for seqi=1:seqAmt
    axonlen{seqi}=sum(seqType{seqi}==-1);
    axonseq{seqi}=seqlist{seqi}(seqType{seqi}==-1);
end

res=struct('seqAmt',seqAmt);
res.seqlen=seqlen; res.axonlen=axonlen;

%%% path time (only at axon segment)
PT=cell(seqAmt,1); htT=zeros(seqAmt,1);
for seqi=1:seqAmt
    tp=stat.meanTime{seqi}(seqType{seqi}==-1);
    PT{seqi}=diff(tp);
    % time of traveling the whole path
    if isempty(PT{seqi})
        htT(seqi)=NaN;
    else
        htT(seqi)=tp(end)-tp(1);
    end
end
res.PT=PT; res.htT=htT;

%%% path length (axon only)
% channel -> electrode location
tss=cell(seqAmt,1);
for seqi=1:seqAmt
    tss{seqi}=chID(axonseq{seqi});
end
%
[~,PL,PLseq]=csdist(tss,CL);
% unit to (m)
PL=PL*ied/1e6;
res.PL=PL;
res.PLseq=PLseq;

%%% Speed of conduction
res.apspeed=PL./htT;

%%% CS在参与的各个channel的所有spike中所占的比例。
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

% % seq中相邻两者的距离过大的数量
% thres=sqrt(10);
% mark=(res.seqPLmax>thres);
% sm.ndbig=sum(mark);

end
