% Locate conduction signal in the original spike train data.
% Find the trace by the marked "top" units, Use results from spikemerge().
%   [seq,seqMark,seqType,seqCS,stat]=spikecs('spktrain',ST,marklb,info)
%   [seq,seqMark,seqType,seqCS,stat]=spikecs('rawsig',ST,marklb,info,X,chID)
% seqType is property of each member in a CS. 1=marker; 2=core member;
% 3x=edge member - 31=low amp with no large time diff, 32=low amp with
% large time delay (could be the signal of post-synaptic potential?)
% In method=='rawsig', function need to call spikecs_check().
function [seq,seqMark,seqType,seqCS,stat]=spikecs(method,SD,ST,marklb,info,X,chID)
ext=0.0008;%(s)
param=struct();
param.bAutoThres=true;
param.minThres=4/60*info.TimeSpan;
param.numPthr=1e-11;%<<<<
param.timeSpan=info.TimeSpan;
param.DTJthr=0.1/1000;%(ms/1000) = 2 points at 20000HZ SR. 

%%% Process
% Method type and Data format
if strcmp(method,'spktrain') % find CS by scanning spike train data
    method=1;
elseif strcmp(method,'rawsig') % find CS by scanning raw signal
    method=2;
else
    error('invalid method');
end
% more info
chAmt=length(SD);
sAmt=cellstat(SD,'length');
assert(isequal(sAmt,cellstat(marklb,'length')),'SD, marklb spike num');
assert(isequal(sAmt,cellstat(ST,'length')),'SD, ST spike num');
stat=struct();

% Get number of "marker" spikes of each unit.
markSDnum=zeros(chAmt,1);
for chi=1:chAmt
    markSDnum(chi)=sum(marklb{chi}==1);
end
% Filt units with enough marker. (will be the number of CS)
markI=find(markSDnum>0); markAmt=length(markI);
% End game if no marker neuron is found.
if markAmt==0
    seq=[];    seqMark=[];    seqCS=[];
    return
end


%%%%%%%%%%%%
% Keep only useful marker channels
marklb=marklb(markI); stat.markSDnum=markSDnum(markI);
seqMark=markI;% marker neuron ID.
% Get ST of markers only.
mst=cell(markAmt,1); msd=mst;
for mi=1:markAmt
    mst{mi}=ST{markI(mi)}(marklb{mi}==1);
    msd{mi}=SD{markI(mi)}(marklb{mi}==1);
end

seq=cell(markAmt,1);
seqCS=cell(markAmt,1);
stat.meanTime=seq; stat.fospknum=seq;
if method==1
    %%% Use the markers to re-scan the ST.
    for mi=1:markAmt
        % 当前marker neuron spike里面 marklb==1者选出，并就这些进行全ST扫描
        % * ST里面自身也要包括在内 （只指自身那些marklb==1的spike）
        [Ichi,Isd,rnum,foI]=cspair(mst{mi},ST,ext,param);
        foAmt=length(foI);
        if foAmt<=1
            stat.fospknum{mi}=rnum;
            continue
        end
        
        %%% Add spikecs_check        
        tx=X(:,chID(foI));
        [~,IR1,IR2]=spikecs_check(tx,msd(mi),info,{Ichi});
        %        
        foI=foI(IR1);
        foAmt=sum(IR1);
        rnum=rnum(IR1);        Ichi=Ichi(:,IR1);        Isd=Isd(IR1);
        if foAmt<=1
            stat.fospknum{mi}=rnum;
            continue
        end
        
        % seq channels composition
        tp=(1:chAmt); seq{mi}=tp(foI);
        stat.fospknum{mi}=rnum; % number of followers in each channel

        % Sort by time.
        % * time measured by difference to "base time" - marker's time.
        mt=zeros(foAmt,1); %Average time
        for k=1:foAmt
            dt=ST{foI(k)}(Isd{k}) - mst{mi}(Ichi(:,k));
            mt(k)=mean(dt);
        end        
        [mt,I]=sort(mt,'ascend');
        seq{mi}=seq{mi}(I);
        seqCS{mi}=Isd(I);
        stat.fospknum{mi}=stat.fospknum{mi}(I);
        stat.meanTime{mi}=mt;
    end    
    seqType=[];
    
    % 如果没有达标的follower,不存在以这个通道为marker的seq,取消。
    temp=cellstat(stat.fospknum,'length');
    rml=(temp<=1);
    if sum(rml)>0
        seq(rml)=[];
        seqMark(rml)=[];
        seqCS(rml)=[];
        stat.markSDnum(rml)=[];
        stat.meanTime(rml)=[];
        stat.fospknum(rml)=[];
        marklb(rml)=[];
    end
    
else
    %%% Use the markers to re-scan the raw signal.
    xcha=size(X,2);
    seqType=cell(markAmt,1);
    seqMark=chID(seqMark);
    
    [sttp,IR1,IR2,ID]=spikecs_check(X,msd,info); % IR2 is of higher standard.    
    for mi=1:markAmt
        temp=zeros(xcha,1);
        temp(IR1(:,mi))=31;
        temp(IR2(:,mi))=21; % core member
        I=IR1(:,mi)&ID(:,mi); temp(I)=32;
        I=IR2(:,mi)&ID(:,mi); temp(I)=22;
        temp(seqMark(mi))=1;
        
        seqType{mi}=temp(temp~=0);
        
        seq{mi}=find(IR1(:,mi));
        
        % sort by time
        mt=sttp.meanTime{mi}(IR1(:,mi));
        [mt,I]=sort(mt,'ascend');
        seq{mi}=seq{mi}(I);
        stat.meanTime{mi}=mt;
    end
    seqCS=[]; % <<< no way to do the CS right now.
    
    % 如果没有达标的follower,不存在以这个通道为marker的seq,取消。
    temp=sum(IR1);
    rml=(temp<=1);
    if sum(rml)>0
        seq(rml)=[];
        seqMark(rml)=[];
        seqType(rml)=[];
        marklb(rml)=[];
    end
end

stat.marklb=marklb;
end % of main