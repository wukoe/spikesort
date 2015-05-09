% [SD,rmlist,seqlist,seqCS,stat,idxS]=spikemerge(SD,info,SA)
% 
function [newSD,rmlist,seqlist,seqCS,stat,idxS]=spikemerge2(SD,info,SA)
% standard to find cluster sequence.
ext=0.00025; %0.25ms
totaltimelimit=0.002; % full length
optSelect=2; % way to determine the representative of the cluster

% Change time to fit SD
ext=ceil(ext*info.srate);
totaltimelimit=ceil(totaltimelimit*info.srate);

chAmt=length(SD);
% Get number of firing in each channel
sAmt=cellstat(SD,'length');
% Check number of SD SA
assert(isequal(sAmt,cellstat(SA,'length')),'merging: SD and SA number not equal');

%%%
% Find all tight-following spikes (under ext) in the channel pairs
wd=[SD{chi}-ext,SD{chi}+ext];


end % of main

function [Iw,Isd]=cspair(wd,SD)
wdlen=size(wd,1);
cha=length(SD);
Iw=false(wdlen,cha);
for chi=1:cha
    [BI,BA]=binid(SD{chi},wd);
    Iw(:,chi)=(BA>0);
    Isd=(BI>0);
    num=sum(Iw);
end
end

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