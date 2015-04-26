% Convert spike train from index to logical format
%   SD=time2logic(ST,'time',T)
%   parameter {'bin size',bs} specify size; {'bin',bm} specify bins.
function SD=time2logic(ST,varargin)
bTime=false;
binSize=0.01; % (s)
bAutoBin=true;

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
            case 'bin'
                bAutoBin=false;
                binm=pinfo{parai};
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
    % Determine the bins
    if bAutoBin
        ts=min(cellstat(ST,'min'));
        te=max(cellstat(ST,'max'));
        binm=(ts:binSize:te)';
        binm=[binm(1:end-1),binm(2:end)];
    end
    bAmt=size(binm,1);

    % Transfer SD to bin data
    BD=false(bAmt,chAmt);
    for chi=1:chAmt
        [~,BA]=binid(BD{chi},binm);
        BD(:,chi)=(BA>0);
    end
end