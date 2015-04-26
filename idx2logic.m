% SD=idx2logic(SI)
% SD=idx2logic(SI,dlen)
function SD=idx2logic(SI,varargin)

cha=length(SI);

% max of idx
maxidx=max(cellstat(SI,'max'));
if nargin==2
    ptsAmt=varargin{1};
    if ptsAmt<maxidx
        error('the length can not hold biggest index number');
    end    
else
    ptsAmt=maxidx;
end

%
SD=false(ptsAmt,cha);
for chi=1:cha    
    SD(SI{chi},chi)=true;
end