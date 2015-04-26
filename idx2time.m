% transfer the SD from index of sampling point to exact time of spike
%   ST=idx2time(SI,T)
function ST=idx2time(SI,T)
if iscell(SI)
    cha=length(SI);
    ST=cell(cha,1);
    for chi=1:cha
        ST{chi}=T(SI{chi});
    end
else
    ST=T(SI);
end

end