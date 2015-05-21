function [seq,seqMark,seqCS,stat]=spikecs(SD,marklb,rmlb)

%%% Find the trace by the marked "top" units
stat=struct();
% Get number of "marker" spikes of each unit.
markSDnum=zeros(chAmt,1);
for chi=1:chAmt
    markSDnum(chi)=sum(marklb{chi}==1);
end
% Filt out units with enough marker.
markI=find(markSDnum>0); markAmt=length(markI);

if markAmt>0
    nsd2=SD(markI);
    marklb=marklb(markI); stat.markSDnum=markSDnum(markI);
    seqMark=markI;
    
    % Use the markers to re-scan the SD
    seq=cell(markAmt,1);
    seqCS=cell(markAmt,1);
    stat.meanTime=seq;
    stat.fonum=seq;
    rml=[];
    for m=1:markAmt        
        % 挑选出marklb==1者并就这些进行全SD扫描
        temp=nsd2{m}(marklb{m}==1);        
        [Ichi,Isd,rnum]=cspair(temp,SD,ext);
%         % 消除自身
%         Ichi(:,chi)=false;
%         Isd{chi}=[];
%         rnum(chi)=0;
        
        % 只保留足够重复数量的通道
        % * 包括了自身,并且只有自身那些marklb==1的spike
        foI=find(rnum>=minThres); foAmt=length(foI); % foI的取值范围是全部120通道
        if foAmt<=1, rml=[rml,m]; end
        Ichi=Ichi(:,foI); Isd=Isd(foI); 
        % seq channels
        tp=(1:chAmt); seq{m}=tp(foI);
        stat.fonum{m}=rnum(foI); % followers in each channel
        
        % Average time
        mt=zeros(foAmt,1);
        for k=1:foAmt
            mt(k)=mean(SD{foI(k)}(Isd{k}));
        end
        % sort by time
        [mt,I]=sort(mt,'ascend');
        seq{m}=seq{m}(I);
        seqCS{m}=Isd(I);
        stat.fonum{m}=stat.fonum{m}(I);
        stat.meanTime{m}=mt;
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
else
    seq=[];
    seqMark=[];
    seqCS=[];
end

stat.marklb=marklb;
end % of main

% %%%
% [stat,seqlist,seqCS]=spikerepeat(SD,info.TimeSpan,'ext',ext,'length limit',totaltimelimit);
% seqAmt=length(seqlist);
% 
% %%% Removal procedure
% % Initiate rmlist, each spike has individual marker. 0=not removed,
% % other=index of channel that this spike is removed according to.
% rmlist=cell(chAmt,1);
% for chi=1:chAmt
%     rmlist{chi}=false(sAmt(chi),1);
% %     rmlist{chi}=zeros(sAmt(chi),1);
% end
% % Remove
% idxS=zeros(seqAmt,1);
% newSD=SD;
% for seqi=1:seqAmt
%     seq=seqlist{seqi};
%     seqlen=stat.seqlen(seqi);
%     if optSelect==1 % keep channel with maximum total firing.
%         % percentage of repeats to each unit's total firing number
%         P=stat.seqcount(seqi)./sAmt(seq);
%         [~,idx]=max(P);
%     else % keep channel with maximum spike amplitude.
%         % Get amplitude of spikes in each channel (use only spikes participating current cluster,
%         % remove possible artifact)
%         sAmp=zeros(seqlen,1);
%         for k=1:seqlen
%             I=seqCS{seqi}(:,k);
%             temp=abs(SA{seq(k)}(I));
%             I=outlier_detect(temp);
%             temp=temp(~I);
%             sAmp(k)=mean(temp);
%         end
%         [~,idx]=max(sAmp);
%     end
%     
%     % Keep only "keep" channel, delete others. (mark in "remove list")
%     for k=1:seqlen
%         if k~=idx
%             I=seqCS{seqi}(:,k);
%             rmlist{seq(k)}(I)=true;
% %             rmlist{seq(k)}(I)=seq(idx); % This can tell one spike is deleted by which channel's "main" spikes
%         end
%     end
%     
%     % Save selected
%     idxS(seqi)=idx;
% end
% % Do deletion
% for chi=1:chAmt
%     I=rmlist{chi};
% %     I=(rmlist{chi}>0);
%     newSD{chi}(I)=[];
% end