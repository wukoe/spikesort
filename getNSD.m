% [NSD,chID]=getNSD(reconSD,CST)
% Noise type is removed.
function [NSD,chID]=getNSD(reconSD,CST)
NSD=cell(0,1);
chID=zeros(0,1);
chcount=0;
for chi=1:length(reconSD)
    spkclu=reabylb(CST{chi});
    % only choose non-noise channel
    nzI=find(spkclu.types>0);  na=length(nzI);
    if na>0
        temp=cell(na,1);
        for k=1:length(nzI)
            I=spkclu.ids{nzI(k)};
            temp{k}=reconSD{chi}(I);
        end
        NSD(chcount+1:chcount+na,1)=temp;
        chID(chcount+1:chcount+na,1)=chi;
        chcount=chcount+na;
    end
end

end