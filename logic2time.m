% please use logic2idx(SD,T) for this function
function ST=logic2time(SL,T)
SI=logic2idx(SL);
ST=idx2time(SI,T);
end