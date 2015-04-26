%   SI=logic2idx(SD) 
% output the index
%   ST=logic2idx(SD,T) 
% This is actually "logic2time". It give time of occurance instead of
% index.
% SI is a cell array, one for each channel
function SI=logic2idx(SD,varargin)
chAmt=size(SD,2);
if nargin==2
    T=varargin{1};
    if length(T)~=size(SD,1)
        error('SD and T length mismatch');
    end
    useT=true;
else
    useT=false;
end

SI=cell(chAmt,1);
if useT
    for chi=1:chAmt
        SI{chi}=T(SD(:,chi));
    end
else
    for chi=1:chAmt
        SI{chi}=find(SD(:,chi));
    end
end
