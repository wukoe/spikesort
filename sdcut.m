% to cut a specified interval out from (time) indexed spike data
%   Y=sdcut(X,[segS,segE]) get data in [segS,segE].
%   Y=sdcut(X,[segS,segE],'del') delete data in [segS,segE] (get data other than [segS,segE]).
% both X and Y are cell structure, each for one channel. segS and segE are
% with unit of (s).
function X=sdcut(X,iv,varargin)
cha=length(X);

flagDel=false;
if nargin==3
    if strcmp(varargin{1},'del')
        flagDel=true;
    else
        error('invalid option');
    end
end

if flagDel
    for chi=1:cha
        I=(X{chi}>=iv(1)) & (X{chi}<=iv(2));
        X{chi}(I)=[];
    end
else
    for chi=1:cha
        I=(X{chi}>=iv(1)) & (X{chi}<=iv(2));
        X{chi}=X{chi}(I);
    end
end

end