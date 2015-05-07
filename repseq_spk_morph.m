% check the aligned spike form of repeat spike sequence.
%   repseq_spk_morph(X,SD,SDchID,seq,seqCS)
% seq and seqCS is of one specific sequence. e.g. seq=[3,10],
% seqCS=[seqcount X 2]
% * seq and seqCS. should be from ST converted from same SD here. 
function dx=repseq_spk_morph(X,SD,SDchID,seq,seqCS)
prewin=10; postwin=20;

seqlen=length(seq);
sc=size(seqCS,1);
dx=zeros(prewin+postwin+1,seqlen,sc);
for si=1:sc %stat.seqcount(seqidx)
    sd=SD{seq(1)}(seqCS(si,1));
    dx(:,:,si)=X(sd-prewin:sd+postwin,SDchID(seq));
end
for k=1:seqlen
    subplot(seqlen,1,k);
    plot(squeeze(dx(:,k,:)),'b');
    ylabel(sprintf('ch:%d',SDchID(seq(k))));
end

% display the seq repeat number, and total spike number of each
% participating channel
subplot(seqlen,1,1);
tp=num2str(cellstat(SD(seq),'length')');
str=sprintf('seq count: %d, SD spike num: %s',sc,tp);
title(str);
end