% % mannual separate spike cluster by amplitude.
% % Specify the list here: [ch,unit,loc,thres]
% movl=[21,3,17,-2e-5; 31,3,17,-3e-5; 37,3,13,-2e-6];
% 
% % Start move
% for mi=1:size(movl,1)
%     cst=CST{movl(mi,1)};
%     a=A{movl(mi,1)};    
%     doall=isnan(movl(mi,2)); % whether apply separation on all units.
%     
%     clu=reabylb(cst);
%     % exclude noise
%     idx=find(clu.types==0);
%     clu.types(idx)=[];
%     clu.typeAmt(idx)=[];
%     clu.ids(idx)=[];
%     clu.cAmt=clu.cAmt-1;
%     
%     maxlb=max(clu.types);
%     for k=1:clu.cAmt
%         if doall || (k==movl(mi,2))
%             if clu.types(k)>0
%                 temp=a(movl(mi,3),clu.ids{k}); % amplitude of one spike cluster at specified location
%                 I=(temp>movl(mi,4));
%                 if sum(I)>0 && sum(I)<clu.typeAmt(k)
%                     cst(clu.ids{k}(I))=maxlb+1;
%                     maxlb=maxlb+1;
%                 end
%             end
%         end
%     end
%     CST{movl(mi,1)}=cst;
% end

%%
fnC='05_25_15_I9123_0002_c';
outfile3=load(fnC);
chinfo=outfile3.chinfo;

[outfile3.NSD,tp]=getNSD(outfile3.reconSD,CST);
chinfo.CluchID=tp;
chinfo.CluchAmt=length(chinfo.CluchID);
outfile3.chinfo=chinfo;

% save
save(fnC,'-struct','outfile3');