% transfer the SD from index of sampling point to exact time of spike
%   ST=idx2time(SI,T)
%   ST=idx2time(SI,sampleRate)
function ST=idx2time(SI,T)
if length(T)==1 % the input is sampling rate rather than Time
    ml=max(cellstat(SI,'max'));
    sr=T;
    T=(1:ml)'/sr;
end

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