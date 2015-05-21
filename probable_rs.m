% 获得特定数量的repeated sequence的概率 - 一次性给出各个数量的概率。
%   P=probable_rs(totaltime,spkAmt,intv,N)
% totaltime: (s)
% spkAmt is number of total spikes in participating channels. length(spkAmt)=2
% means seq length of 2, etc.
% intv: interval length (s)
% N is specified number of repeats.
function P=probable_rs(totaltime,spkAmt,intv,N)
seqlen=length(spkAmt); % = number of channels participating.

if seqlen==2
    P=prob_sl2(totaltime,spkAmt,intv,N);
else
    ER=0:50; % the last one should be 1
    ERlen=length(ER);
    
    P=prob_sl2(totaltime,spkAmt,intv,ER);
    for si=3:seqlen
        EP=zeros(ERlen,ERlen);
        for ei=1:ERlen
            % each row is for fixed number of spikes, each column for fixed
            % number of repeats.
            tp=prob_sl2(totaltime,[ER(ei),spkAmt(si)],intv,ER);
            EP(ei,:)=tp*P(ei);
        end
        % expected number
        P=sum(EP);
    end
    P=P(N+1);
end

end % main


%%%%%%%%%%%%%% 
% Probability of seq length of 2.
function p=prob_sl2(totaltime,spkAmt,intv,N)
nl=length(N);
% Time length occupied by spikes of first channel.
ot=spkAmt(1)*intv;

% Expected number of the repeat
rT=ot/totaltime * spkAmt(2);
% probability based on Poisson.
p=zeros(1,nl);
for k=1:nl
    if N(k)>spkAmt(1)
        p(k)=0;
    else
        if N>50 % Use Stirling's approaximation
            p(k)=exp(-rT)*(exp(1)*rT/N(k))^N(k)/sqrt(2*pi*N(k));
        else
            p(k)=rT^N(k)*exp(-rT)/prod(1:N(k));
        end
    end
end

end