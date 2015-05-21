% Look for the closely-intervaled spike pairs between specified channels
% 相对当前通道（由marksd代表，化为wd的segment形式），其他所有通道的follower的出现位置。
%   [Iw,Isd,num]=cspair(marksd,ST,ext)
function [Iw,Isd,num]=cspair(marksd,ST,ext)
marklen=length(marksd);
cha=length(ST);

wd=[marksd-ext,marksd+ext];
Iw=false(marklen,cha);
Isd=cell(cha,1);
for chi=1:cha    
    [BI,BA]=binid(ST{chi},wd);
    Iw(:,chi)=(BA>0);
    Isd{chi}=(BI>0);
end
num=sum(Iw);

end