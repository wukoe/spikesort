%   SI=time2idx(ST,T)
% ST and T both in (s)
function SI=time2idx(ST,T)
cha=length(ST);

SI=cell(cha,1);
for chi=1:cha
    SI{chi}=closest(ST{chi},T);
    fprintf('|');
end
fprintf('\n');

end