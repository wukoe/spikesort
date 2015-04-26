% mannual separate spike cluster by amplitude.
% Specify the list here: [ch,loc,thres]
movl=[];

% Start move
for mi=1:size(movl,1)
    cst=CST{movl(mi,1)};
    a=A{movl(mi,1)};
    clu=reabylb(cst);
    maxlb=max(clu.types);
    for k=clu.cAmt
        if clu.types(k)>0
            temp=a(movl(mi,2),clu.ids{k}); % amplitude of one spike cluster at specified location
            I=(temp>movl(mi,3));
            if sum(I)>0 && sum(I)<clu.typeAmt(k)
                cst(clu.ids{k}(I))=maxlb+1;
                maxlb=maxlb+1;
            end
        end
    end
    CST{movl(mi,1)}=cst;
end