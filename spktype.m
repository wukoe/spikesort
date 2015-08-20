% Analyze spike types, and detect curve peak. 
%   [stype,ploc,pamp,ssd]=spktype(ALD)
% A= multi align matrix.
% pamp= amplitude with polarity (non abs).
function [stype,ploc,pamp,ssd,flagPeakLocChange]=spktype(ALD)
mcthr=1; % Threshold (fold of STD) of mean curve for determining peaks.

dlen=size(ALD,1);

%%% 
% Mean curve.
mc=mean(ALD,2);
% Get pseudo-STD of each point - outlier-proof.
ssd=zeros(dlen,1);
for k=1:dlen
    tp=ALD(k,:)-mc(k);
    ssd(k)=median(abs(tp))/0.6745;
end
ssd=mean(ssd);

% Thres update
mcthr=ssd*mcthr;

% Get all over-thres segments
D=(mc>mcthr);
psm=continuous_segment(D);
D=(mc<-mcthr);
nsm=continuous_segment(D);
sm=[psm;nsm];
[~,I]=sort(sm(:,1));
sm=sm(I,:);
sma=size(sm,1);

% Amplitude of each segment
amp=zeros(sma,1);
loc=amp;
for smi=1:sma
    temp=mc(sm(smi,1):sm(smi,2));
    [~,loc(smi)]=max(abs(temp));
    amp(smi)=temp(loc(smi));
    loc(smi)=loc(smi)+sm(smi,1)-1;
end
[ampabs,I]=sort(abs(amp),'descend');
amp=amp(I);
loc=loc(I);

% Select by rule.
flagPeakLocChange=false;
% if sma==0
%     stype=0;
% elseif sma==1
%     if amp(1)>0
%         stype=1;
%     else
%         stype=-1;
%     end
% else
%     if amp(1)<0
%         if amp(2)>0 && ampabs(2)>ampabs(1)*0.9
%             stype=1;
%             % in this case, the peak location of spike is changed, notify
%             % others about it.
%             flagPeakLocChange=true; 
%         else
%             stype=-1;
%         end
%     else
%         stype=1;
%     end
% end
if sma==0
    stype=0;
else
    stype=sign(amp(1));
end

if stype==1
    I=find(ampabs>0);
    [pamp,idx]=max(amp(I));
    idx=I(idx);
elseif stype==-1
    pamp=amp(1);
    idx=1;
else
    pamp=[]; idx=[];
end
ploc=loc(idx);

end