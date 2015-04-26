% Y=sdmerge(x1,x2,x3,...)
function Y=sdmerge(varargin)
sa=nargin;
if sa<=1
    Y=varargin{1};
    return
end

Y=varargin{1};
if iscell(Y)
    cha=length(Y);
    for si=2:sa
        for chi=1:cha
            Y{chi}=[Y{chi};varargin{si}{chi}];
        end
        if ~isempty(Y{chi})
        Y{chi}=sort(Y{chi});
        end
    end
else
    for si=2:sa
        Y=[Y;varargin{si}];
    end
    Y=sort(Y);
end

end