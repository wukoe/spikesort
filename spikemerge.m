%   [SD,rmlist,seqlist,seqCS,stat,idxS]=spikemerge(SD,info,SA)
% ! currently, this must be used on sorted data.
function [SD,rmlist]=spikemerge(SD,info,SA)
% standard to find cluster sequence.
ext=0.0005; %0.25ms
% totaltimelimit=0.002; % full length
minThres=4/60*info.TimeSpan;

% Change time to fit SD
ext=ceil(ext*info.srate);
% totaltimelimit=ceil(totaltimelimit*info.srate);

chAmt=length(SD);
% Get number of firing in each chAmtnnel
sAmt=cellstat(SD,'length');
% Check number of SD SA
assert(isequal(sAmt,cellstat(SA,'length')),'merging: SD and SA number not equal');

%%% Removal procedure
rmlist=cell(chAmt,1);
for chi=1:chAmt
    rmlist{chi}=false(sAmt(chi),1);
end
% 
for chi=1:chAmt
    if ~isempty(SD{chi})
        wd=[SD{chi}-ext,SD{chi}+ext];
        for m=1:chAmt
            if m~=chi && ~isempty(SD{m})
                % Find all tight-following spikes (under ext) in the other channel
                [BI,BA]=binid(SD{m},wd);
                Ichi=(BA>0);
                Im=(BI>0);
                num=sum(Ichi);
                % Whether number of repeats over threshold
                if num>minThres
                    % Find mean amplitude of chi and m
                    ampchi=abs(mean(SA{chi}(Ichi)));
                    ampm=abs(mean(SA{m}(Im)));
                    if ampchi>ampm
                        SD{m}(Im)=[];
                        rmlist{m}(Im)=true;
                    else
                        SD{chi}(Ichi)=[];
                        wd(Ichi,:)=[];
                        rmlist{chi}(Ichi)=true;
                    end
                end
            end
        end
    end
    fprintf('|');
end
fprintf('\n');

end