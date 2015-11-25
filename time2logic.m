% Convert spike train from index to logical format
%   SD=time2logic(ST,'time',T)
%   parameter {'bin size',bs} specify size (s);
function SD=time2logic(ST,varargin)
bTime=false;
binSize=0.01; % (s)

% User input
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'time'
                bTime=true;
                T=pinfo{parai};
            case 'bin size'
                binSize=pinfo{parai};
            otherwise
                error('unidentified options');
        end
    end
end


if bTime
    SI=time2idx(ST,T);
    SD=idx2logic(SI,length(T));
else
    chAmt=length(ST);
%     ts=min(cellstat(ST,'min'));
    te=max(cellstat(ST,'max'));
    
    % Determine the bins
    bAmt=ceil(te/binSize);

    % Transfer SD to bin data
    SD=false(bAmt,chAmt);
    for chi=1:chAmt
        I=ceil(ST{chi}/binSize);
        SD(I,chi)=true;
        fprintf('|');
    end
    fprintf('\n');
end

end