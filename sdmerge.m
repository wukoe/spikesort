% Merge together multiple SD data sets, and order train by time.
%   Y=sdmerge(x1,x2,x3,...)
function Y=sdmerge(varargin)
sa=nargin;
if sa<=1
    Y=varargin{1};
    return
end

% Initialize
Y=varargin{1};
if iscell(Y)
    cha=length(Y);
    % Append
    for si=2:sa
        for chi=1:cha
            Y{chi}=[Y{chi};varargin{si}{chi}];
        end        
    end
    % Sort to right time order
    for chi=1:cha
        if ~isempty(Y{chi})
            Y{chi}=sort(Y{chi},'ascend');
        end
    end
else
    for si=2:sa
        Y=[Y;varargin{si}];
    end
    Y=sort(Y);
end

end