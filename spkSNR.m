% Signal-noise ratio of spikes. "signal" measured by peak amplitude.
%   [R,ap,ssd]=spkSNR(A)
% A= aligned spikes. Can be matrix of aligments or cell of multiple matrix.
% R= SNR; ap= absolute peak amplitude.
% spkSNR only calculate for aligements with >=5 spikes (unfit is set 0).
function [R,ap,ssd]=spkSNR(A)
if ~iscell(A)
    A={A};
end    
aa=cellstat(A,'size',2);
cha=length(A);

I=find(aa>=5);
ap=zeros(cha,1); v=ap;
for m=1:length(I)    
    chi=I(m);
    
    % mean of all curves.
    tpA=A{chi}';
    mc=mean(tpA); 
    [~,idx]=max(abs(mc));
    ap(chi)=mc(idx);
    % ap=max(mc)-min(mc);
    
    % pseudo-STD of all curves.    
    tpA=bsxfun(@minus,tpA,mc);
    ssd=median(abs(tpA))/0.6745;
    v(chi)=mean(ssd); % choose mean of all points
%     v(chi)=mean(std(tpA));
end    

% Ratio
R=zeros(cha,1);
R(I)=(abs(ap(I))./v(I))';
end