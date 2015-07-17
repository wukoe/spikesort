% SD=idx2logic(SI)
% SD=idx2logic(SI,dlen)
% SD=idx2logic(SI,dlen,interval) interval in (pts)
function SD=idx2logic(SI,varargin)
cha=length(SI);

maxidx=max(cellstat(SI,'max')); % max index number of all channels
if nargin>=2
    if isempty(varargin{1})
        ptsAmt=maxidx;
    else
        ptsAmt=varargin{1};
        if ptsAmt<maxidx
            error('the length can not hold biggest index number');
        end    
    end
end
if nargin==3
    interval=varargin{2};
    ptsAmt=ceil(ptsAmt/interval);
else
    interval=20;
end

%
SD=false(ptsAmt,cha);
for chi=1:cha
    I=ceil(SI{chi}/interval);
    SD(I,chi)=true;
end

end