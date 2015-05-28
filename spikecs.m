% Locate conduction signal in the original spike train data.
% Use results from spikemerge().
function [seq,seqMark,seqCS,stat]=spikecs(ST,marklb,info)
chAmt=length(ST);
ext=0.001;
minThres=4/60*info.TimeSpan;

%%% Find the trace by the marked "top" units
stat=struct();
% Get number of "marker" spikes of each unit.
markSDnum=zeros(chAmt,1);
for chi=1:chAmt
    markSDnum(chi)=sum(marklb{chi}==1);
end
% Filt out units with enough marker.
markI=find(markSDnum>0); markAmt=length(markI);

% End game if no marker neuron is found.
if markAmt==0
    seq=[];    seqMark=[];    seqCS=[];
    stat.marklb=marklb;
    return
end

%%% Use the markers to re-scan the SD.
nsd2=ST(markI);
marklb=marklb(markI); stat.markSDnum=markSDnum(markI);
seqMark=markI;% marker neuron ID.

% Take mark neurons in turn to scan.
seq=cell(markAmt,1);
seqCS=cell(markAmt,1);
stat.meanTime=seq;
stat.fonum=seq;
rml=[];
for mi=1:markAmt
    % 当前marker neuron spike里面 marklb==1者选出，并就这些进行全ST扫描
    % * ST里面自身也要包括在内 （只指自身那些marklb==1的spike）
    temp=nsd2{mi}(marklb{mi}==1);
    [Ichi,Isd,rnum]=cspair(temp,ST,ext);
    
    % 所有follower通道里面只保留足够重复数量的通道    
    foI=find(rnum>=minThres); foAmt=length(foI);
    Ichi=Ichi(:,foI); Isd=Isd(foI);
    if foAmt<=1, rml=[rml,mi]; end % 如果没有follower,这个通道不存在seq,取消。
    % seq channels composition
    tp=(1:chAmt); seq{mi}=tp(foI);
    stat.fonum{mi}=rnum(foI); % number of followers in each channel
    
    % Average time
    mt=zeros(foAmt,1);
    for k=1:foAmt
        mt(k)=mean(ST{foI(k)}(Isd{k}));
    end
    % sort by time
    [mt,I]=sort(mt,'ascend');
    seq{mi}=seq{mi}(I);
    seqCS{mi}=Isd(I);
    stat.fonum{mi}=stat.fonum{mi}(I);
    stat.meanTime{mi}=mt;
end
if ~isempty(rml)
    seq(rml)=[];
    seqMark(rml)=[];
    seqCS(rml)=[];
    stat.markSDnum(rml)=[];
    stat.meanTime(rml)=[];
    stat.fonum(rml)=[];
    marklb(rml)=[];
end
stat.marklb=marklb;

end % of main