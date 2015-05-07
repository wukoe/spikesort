% SD=idx2logic(SI)
% SD=idx2logic(SI,dlen)
% SD=idx2logic(SI,dlen,interval)
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
    interval=1;
end

%
SD=false(ptsAmt,cha);
for chi=1:cha
    tp=ceil(SI{chi}/interval);
    SD(tp,chi)=true;
end

end