% A=spkalplot_cs(X,seq,seqMark,seqCS,SD,chID,CL,elecLabel)
function varargout=spkalplot_cs(X,seq,seqMark,seqCS,SD,chID,CL,elecLabel)
% <<<< the signal is not aligned to mark but to their own signals.
srate=20000;
pagemaxch=7;
pagemaxseq=8;
bSameScale=true;
pathLayoutOpt=1;

% Proc
seqAmt=length(seq);
seqlen=cellstat(seq,'length');
chMax=max(seqlen);% max number of channels used.

%%%
AD=cell(seqAmt,1);
seqcount=cell(seqAmt,1);
for seqi=1:seqAmt
    tsd=cell(seqlen(seqi),1);
    % Use the marker spikes to align all channels.    
    tp=SD{seqMark(seqi)};
    idx=find(seq{seqi}==seqMark(seqi));
    tp=tp(seqCS{seqi}{idx});
    for k=1:seqlen(seqi)
        tsd{k}=tp;
    end
    AD{seqi}=spike_align(X(:,chID(seq{seqi})),tsd,srate,'window',[-1,1]);
    seqcount{seqi}=cellstat(seqCS{seqi},'sum');
end
alignlen=size(AD{1}{1},1);

%%% Initialize GUI controls
% * location start from left bottom of page
vpos=0;
hpos=20;
Fig=figure('Position',[100 10 1000 900]);%,'ResizeFcn',@resizegui); 
HpreB=uicontrol('Style','pushbutton','String','Up','Position',[hpos,vpos,50,25],'Callback',@PreButton_CB);
HpageEdit=uicontrol('Style','edit','Position',[hpos+50,vpos,30,25],'Callback',@PageEdit_CB);
HnextB=uicontrol('Style','pushbutton','String','Down','Position',[hpos+80,vpos,50,25],'Callback',@NextButton_CB);
HpreseqB=uicontrol('Style','pushbutton','String','<--','Position',[hpos+140,vpos,50,25],'Callback',@PreSeqButton_CB);
HnextseqB=uicontrol('Style','pushbutton','String','-->','Position',[hpos+190,vpos,50,25],'Callback',@NextSeqButton_CB);
HselectEdit=uicontrol('Style','edit','Position',[hpos+250,vpos,100,25],'Callback',@SelectEdit_CB);
HeqchCheck=uicontrol('Style','checkbox','String','same order','Position',[hpos+360,vpos,100,20],'Callback',@OrderCk_CB);
% HcurvePop=uicontrol('Style','popup','String',{'samples','mean'},'Position',[hpos+470,vpos,100,25],'Callback',@CurvePop_CB);
HtextDisp=uicontrol('Style','text','Position',[hpos+470,vpos,200,25]);


% Control status.
flagEqChOrder=false;
ckminVal=get(HeqchCheck,'Min');
set(HeqchCheck,'Value',ckminVal);
chAmt=chMax;
% number of pages.
pgAmt=[ceil(seqAmt/pagemaxseq),ceil(chAmt/pagemaxch)]; %[hor, ver]

% Draw 1st page
% Union of channels used in all seqs
chColle=[];
for k=1:seqAmt
    chColle=union(chColle,seq{k});
end

selectSeqAmt=seqAmt;
selectSeqI=(1:seqAmt);
pgi=[1,1]; % seq,ch
drawpage(pgi);

%%% Output
if nargout==1
    varargout{1}=AD;
end


%%%%%%%%%%%%% Callback functions for the controls.
%%% pre button
    function PreButton_CB(hObject,eventdata)
       if pgi(2)<=1
           pgi(2)=1;
       else
           % Plot
            pgi(2)=pgi(2)-1;
           drawpage(pgi);
       end
    end
    
%%% next button
    function NextButton_CB(hObject,eventdata)
       if pgi(2)>=pgAmt(2)
           pgi(2)=pgAmt(2);
       else
           % Plot
           pgi(2)=pgi(2)+1;
           drawpage(pgi);
       end
    end

%%% up button
    function PreSeqButton_CB(hObject,eventdata)
       if pgi(1)<=1
           pgi(1)=1;
       else
           % Plot
            pgi(1)=pgi(1)-1;
           drawpage(pgi);
       end
    end
    
%%% next button
    function NextSeqButton_CB(hObject,eventdata)
       if pgi(1)>=pgAmt(1)
           pgi(1)=pgAmt(1);
       else
           % Plot
           pgi(1)=pgi(1)+1;
           drawpage(pgi);
       end
    end

%%% page number text input
    function PageEdit_CB(hObject,eventdata)
        tp=get(HpageEdit,'String');
        pgi(2)=str2double(tp);
        if pgi(2)<1, pgi(2)=1; end
        if pgi(2)>pgAmt(2), pgi(2)=pgAmt(2); end
        drawpage(pgi);
    end

%%% seq number text input
    function SelectEdit_CB(hObject,eventdata)
        tp=get(HselectEdit,'String');
        if isempty(tp) || strcmp(tp,'*')
            selectSeqI=(1:seqAmt);
        else
            selectSeqI=str2num(tp);
        end
        % number of pages.
        selectSeqAmt=length(selectSeqI);
        % Union of channels used in all seqs
        chColle=[];
        for k=1:selectSeqAmt
            chColle=union(chColle,seq{selectSeqI(k)});
        end

        pgAmt(1)=ceil(selectSeqAmt/pagemaxseq); 
        pgi=[1,1];
        drawpage(pgi);
    end

%%% check box equal channel order option
    function OrderCk_CB(hObject,eventdata)
        if get(HeqchCheck,'Value')==get(HeqchCheck,'Min')
            flagEqChOrder=false;
            chAmt=chMax;            
            pgAmt(2)=ceil(chAmt/pagemaxch); % number of pages.%[hor, ver]
        else
            flagEqChOrder=true;
            chAmt=length(chColle);
            pgAmt(2)=ceil(chAmt/pagemaxch); % number of pages.%[hor, ver]
        end
        pgi(2)=1; % only change channel page but not seq page.
        drawpage(pgi);
    end

%%%
    function drawpage(pgi)
        set(Fig,'Name',sprintf('%d seqs %d channels, page%d-%d',selectSeqAmt,chAmt,pgi(1),pgi(2)));
        set(HpageEdit,'String',pgi(2));
        % Clean the page
        subplot(1,1,1);
        % label the processing
        set(HtextDisp,'String','apply re-draw');
        
        pgseqI=(pgi(1)-1)*pagemaxseq+1:min(pgi(1)*pagemaxseq,selectSeqAmt);
        pgchI=(pgi(2)-1)*pagemaxch+1:min(pgi(2)*pagemaxch,chAmt);
        if bSameScale
            tp=[0,0];
            for m=1:length(pgseqI)
                seqi=selectSeqI(pgseqI(m));
                for k=1:seqlen(seqi)
                    temp=mean(AD{seqi}{k},2);
                    tp(1)=min(tp(1),min(temp));
                    tp(2)=max(tp(2),max(temp));
                end
            end
            axisScale=[0,alignlen,tp];
        end
        
        %%%
        for coli=1:length(pgseqI)
            seqi=selectSeqI(pgseqI(coli));
            
            %%% seq MEA layout
            subplot(pagemaxch+1,pagemaxseq,coli);
            % Different sequence path layout options.
            if pathLayoutOpt==0 % basic: "marker" in special color.
                tp=zeros(120,1); tp(chID(seq{seqi}))=0.8; % the participating channels in red.
                tp(chID(seqMark(seqi)))=0.6; % the marker use a different color.
            else % color by the sequence time.
                tp=zeros(120,1);
                tp(chID(seq{seqi}))=linspace(0.2,0.8,seqlen(seqi));
                %<<< change the color difference to match time difference later.
            end
            chlayout(tp,CL,'clim',[0,1]);
            axis off
            % seq infor.
            title(sprintf('#%d(%dch)',seqi,seqlen(seqi)));
            
            %%% spikes of each seq.
            if flagEqChOrder % keep same order of channels
                for rowi=1:length(pgchI)
                    subplot(pagemaxch+1,pagemaxseq,rowi*(pagemaxseq)+coli);
                    chi=chColle(pgchI(rowi));
                    
                    idx=find(seq{seqi}==chi);
                    if ~isempty(idx)
                        temp=AD{seqi}{idx};
                        temp=temp-mean(mean(temp));
                        if chi==seqMark(seqi)
                            plot(temp,'b');
                        else
                            plot(temp,'k');
                        end
                        
                        spknum=seqcount{seqi}(idx);
                    else
                        spknum=0;
                    end
                    
                    if bSameScale
                        axis(axisScale);
                    end
                    axis off
                    % chi info + related spike number.
                    title(sprintf('%s(%d)',elecLabel{chID(chi)},spknum));
                end
                
            else % each seq with its own time order of channels.
                dispchI=pgchI(pgchI<=seqlen(seqi));
                for rowi=1:length(dispchI)
                    subplot(pagemaxch+1,pagemaxseq,rowi*(pagemaxseq)+coli);
                    chi=dispchI(rowi);
                    
                    temp=AD{seqi}{chi};
                    temp=temp-mean(mean(temp));
                    if seq{seqi}(chi)==seqMark(seqi)
                        plot(temp,'b');
                    else
                        plot(temp,'k');
                    end
                    
                    if bSameScale
                        axis(axisScale);
                    end
                    axis off
                    % chi info.
                    title(sprintf('%s(%d)',elecLabel{chID(seq{seqi}(chi))},seqcount{seqi}(chi)));
                end
            end
        end
        set(HtextDisp,'String','');
    end

end