% Locate conduction signal in the original spike train data.
% Find the trace by the marked "top" units, Use results from spikemerge().
%   [seq,seqMark,seqType,seqCS,stat]=spikecs('spktrain',SD,ST,marklb,info,X,chID)
%   [seq,seqMark,seqType,seqCS,stat]=spikecs('rawsig',SD,ST,marklb,info,X,chID)
% seqType is property of each member in a CS. 1=marker; 2=core member;
% 3x=edge member - 31=low amp with no large time diff, 32=low amp with
% large time delay (could be the signal of post-synaptic potential?)
% In method=='rawsig', function need to call spikecs_check().
function [seqlist,seqMark,seqType,seqCS,stat,SD,ST]=spikecs(method,SD,ST,marklb,info,X,chID)
% 
param=struct();
param.ext=0.001;%(s)
param.bAutoThres=true;
param.minThres=4/60*info.TimeSpan;
param.numPthr=1e-10;% here the thres should be more strict than spikemerge.m
param.timeSpan=info.TimeSpan;
param.DTJthr=0.15/1000;%(ms/1000) = 2 points at 20000HZ SR. 

twin=[-1,1];
flagSeparateIchi=false;

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
% xcha=size(X,2);
sAmt=cellstat(SD,'length');
assert(isequal(sAmt,cellstat(marklb,'length')),'SD, marklb spike num');
assert(isequal(sAmt,cellstat(ST,'length')),'SD, ST spike num');

% Get number of "marker" spikes of each unit.
markSDnum=zeros(chAmt,1);
for chi=1:chAmt
    markSDnum(chi)=sum(marklb{chi}==1);
end
% Filt units with enough marker. (will be the number of CS)
markI=find(markSDnum>0); markAmt=length(markI);
% End game if no marker neuron is found.
if markAmt==0
    seqlist=[];    seqMark=[];    seqCS=[];
    return
end


%%%%%%%%%%%%
stat=struct();
% Keep only useful marker channels
marklb=marklb(markI); stat.markSDnum=markSDnum(markI);
seqMark=markI;% marker neuron ID.
% Get ST of markers only.
mst=cell(markAmt,1); msd=mst;
for mi=1:markAmt
    mst{mi}=ST{markI(mi)}(marklb{mi}==1); % <<<到底需要不需要marklb
    msd{mi}=SD{markI(mi)}(marklb{mi}==1);
end
sAmt=sAmt(markI);

seqlist=cell(markAmt,1);
seqCS=cell(markAmt,1); seqType=cell(markAmt,1);
stat.meanTime=seqlist; stat.fospknum=seqlist; 
stat.IW=seqlist; stat.spkAmp=seqlist; stat.spkSSD=seqlist;
% flagAnyPeakLocChange=false;
if method==1 % Use the markers to re-scan the ST.
    for mi=1:markAmt
        %%% Scan the ST.
        % 当前marker neuron spike里面 marklb==1者选出，并就这些进行全ST扫描
        % * ST里面自身也要包括在内 （只指自身那些marklb==1的spike）
        [Ichi,Isd,rnum,foI]=cspair(mst{mi},ST,sAmt(mi),param);
        foAmt=length(foI);
        if foAmt<=1
            stat.fospknum{mi}=rnum;
            continue
        end
        unionIchi=any(Ichi,2);        
        % From follower channel (foI) to follower electrode.
        foElect=chID(foI);
        
        %%% Add spikecs_check
        % Align the signal in all X channels to one specific channel of SD.        
        tsd=cell(foAmt,1);
        for xi=1:foAmt
            if flagSeparateIchi
                tsd{xi}=msd{mi}(Ichi(:,xi));
            else
                tsd{xi}=msd{mi}(unionIchi);
            end
        end        
        % Align signal with up-sampling.
        [A,~,O]=spike_align(X(:,foElect),tsd,info.srate,'window',twin);
        % * can not use 'up sample'option, the adjust should only happen on
        % marker channel, others should follow that instead of having their
        % own adjusted realignment.
        UA=cell(foAmt,1);
        idx=(foI==seqMark(mi));
        [UA{idx},~,sttp]=spkUpSamp(A{idx},'realign',[],'center',O.preww+1);
        UAspklen=size(UA{idx},1);
        usSrate=(UAspklen-1)/(size(A{idx},1)-1);
        % then the others
        idx=find(~idx);
        for k=1:foAmt-1
            UA{idx(k)}=spkUpSamp(A{idx(k)},'realign',sttp.newCentIdx);
        end
        
        % Check morphology of alignments.
        [checkstat,IR1]=spikecs_check(UA,info);
        % Filtering according to IR1.
        foI=foI(IR1); foElect=foElect(IR1);
        foAmt=checkstat.rnum1;
        rnum=rnum(IR1);        Ichi=Ichi(:,IR1);        Isd=Isd(IR1);
        ploc=checkstat.ploc(IR1); pamp=checkstat.pamp(IR1);   
        stype=checkstat.stype(IR1); ssd=checkstat.ssd(IR1);
        
        if foAmt<=1
            stat.fospknum{mi}=rnum;
            continue
        end
        
        % Get time by aligned signal.
        idx=(foI==seqMark(mi));
        baseloc=ploc(idx);
        mt=(ploc-baseloc)/(info.srate*usSrate); % (s)
        
        % 当出现两者相同的通道时，再选择amp较大的那个。
        lb=reabylb(foElect);
        for k=1:lb.cAmt
            if lb.typeAmt(k)>1
                tp=abs(pamp(lb.ids{k}));
                [~,idx]=min(tp);
                pamp(lb.ids{k}(idx))=0;%通过将较小者置0的方法
            end
        end
        I=(pamp==0);
        foI(I)=[]; rnum(I)=[]; Ichi(:,I)=[]; Isd(I)=[];
        pamp(I)=[]; stype(I)=[]; ssd(I)=[];
        mt(I)=[];
        
        
%         % 另一种方法：measured by time difference to "base time" - marker's time.
%         mt2=zeros(foAmt,1); %Average time
%         for k=1:foAmt
%             dt=ST{foI(k)}(Isd{k}) - mst{mi}(Ichi(:,k));
%             mt2(k)=mean(dt);
%         end 

        % Sort by time.
        [mt,I]=sort(mt,'ascend');
        % fill output
        seqlist{mi}=foI(I);% seq channels composition
        seqCS{mi}=Isd(I);
        stat.fospknum{mi}=rnum(I);% number of followers in each channel
        stat.meanTime{mi}=mt;
        stat.IW{mi}=Ichi(:,I);
        stat.spkAmp{mi}=pamp(I);
        stat.spkSSD{mi}=ssd(I);
        seqType{mi}=stype(I);
        
%         % when location shift happen in any channel.
%         if flagChPLChange
%             flagAnyPeakLocChange=true;
%         end
    end
    
    % 如果没有达标的follower,不存在以这个通道为marker的seq,取消。
    temp=cellstat(stat.fospknum,'length');
    rml=(temp<=1);
    if sum(rml)>0
        seqlist(rml)=[];
        seqMark(rml)=[];
        seqCS(rml)=[];
        seqType(rml)=[];
        stat.markSDnum(rml)=[];
        stat.meanTime(rml)=[];
        stat.fospknum(rml)=[];
        stat.IW(rml)=[];
        stat.spkAmp(rml)=[];
        stat.spkSSD(rml)=[];
        marklb(rml)=[];
    end
    
else % Use the markers to re-scan the raw signal.
    xcha=size(X,2);
    seqType=cell(markAmt,1);
    seqMark=chID(seqMark);
    
    [stat,IR1,IR2,ID]=spikecs_check(X,msd,info); % IR2 is of higher standard.    
    for mi=1:markAmt
        temp=zeros(xcha,1);
        temp(IR1(:,mi))=31;
        temp(IR2(:,mi))=21; % core member
        I=IR1(:,mi)&ID(:,mi); temp(I)=32;
        I=IR2(:,mi)&ID(:,mi); temp(I)=22;
        temp(seqMark(mi))=1;
        
        seqType{mi}=temp(temp~=0);
        
        seqlist{mi}=find(IR1(:,mi));
        
        % sort by time
        mt=stat.meanTime{mi}(IR1(:,mi));
        [mt,I]=sort(mt,'ascend');
        seqlist{mi}=seqlist{mi}(I);
        stat.meanTime{mi}=mt;
    end
    seqCS=[]; % <<< no way to do the CS right now.
    
    % 如果没有达标的follower,不存在以这个通道为marker的seq,取消。
    temp=sum(IR1);
    rml=(temp<=1);
    if sum(rml)>0
        seqlist(rml)=[];
        seqMark(rml)=[];
        seqType(rml)=[];
        marklb(rml)=[];
    end
end

%%% Output
stat.marklb=marklb;
% stat.flagAnyPeakLocChange=flagAnyPeakLocChange;

end % of main