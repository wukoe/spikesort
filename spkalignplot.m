% Program for making the plot of spikes of one neuron/channel aligned to time of spikes(peak). 
% Use button and text to navigate.
% Accept 2 types of data:
%   spkalignplot(fnX,SD,srate) uses full signal X and index format data spike SD.
% X can be either data in memory or MatFile handle. In this case the actual
% number of channels shown depends on SD.
%   spkalignplot(A) uses aligned data
%   A=spkalignplot() export the aligned data.
% Other varargin:
% 'NCI', (neuron cluster identity) - add clustering information using color;
% 'window',[-pre,post] - specify the time window (ms) for spikes.
% 'limit',0 - no limit for the number of spikes of each plot. non-0 specify the number (limited to
% 100 by default).
% 'elecLabel',L : display electrode names in L.
function varargout=spkalignplot(X,varargin)
%%% Parameters default
% window for spike waveform
pww=[-0.6,0.7]; % pre-event & post-event window width (ms)
bCluster=false;
bShow0=1; % whether to display "outliers" (those clustered id as 0)
bNumLimit=true; % limit the maximum number of curves available in each plot
limitNum=100;

%%% Determine input datatypes (3 different types).
if iscell(X) % 直接输入波形(variable A in spike sorting)的情况
    dataType=1;
    upara=varargin;    
    bXaxis=false;    
elseif isa(X,'matlab.io.MatFile') % 输入文件句柄+SD
    dataType=2;
    SD=varargin{1}; srate=varargin{2};
    if nargin>3
        upara=varargin(3:end);
    else
        upara=[];
    end  
    bXaxis=true;
elseif isnumeric(X) % Input of X in memory + SD
    dataType=3;
    SD=varargin{1}; srate=varargin{2};
    if nargin>3
        upara=varargin(3:end);
    else
        upara=[];
    end
    bXaxis=true;    
else
    error('invalid data type');
end

%%% User input
bChAssign=false;
bUseLabel=false;
if ~isempty(upara)
    [pname,pinfo]=paramoption(upara{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
%             case 'select'
%                 seleI=pinfo{parai};
%                 bChSelect=true;
            case 'NCI' % neuron cluster identity
                CST=pinfo{parai};
                bCluster=true;
            case 'chAssign'
                chAssign=pinfo{parai};
                bChAssign=true;
            case 'window'
                pww=pinfo{parai};
                bXaxis=true;
            case 'limit'
                if pinfo{parai}==0
                    bNumLimit=false;
                else
                    bNumLimit=true;
                    limitNum=pinfo{parai};
                end
            case 'elecLabel'
                eleLabel=pinfo{parai};
                bUseLabel=true;
            otherwise
                error('unidentified options');
        end
    end
end    

%%% Process
% Get Channel number (=total number to display)
if dataType==1
    chAmt=length(X);
else % for both 2 and 3
    chAmt=length(SD);
    % If ch num in SD is different from X, 'chAssign' is mandatory.
    if chAmt~=size(X,2) && ~bChAssign
        disp('warning: SD and X channel number mismatch! consider "chAssign"');
    end
    if bChAssign && length(chAssign)~=chAmt
        error('SD channel number not match to number of chAssign items');
    end
end
if ~bChAssign
    chAssign=1:chAmt;
end

% Check numbers
if bCluster    
    if dataType==1
        sAmt=cellstat(X,'size',2);
        assert(isequal(cellstat(CST,'length'),sAmt),'number of spikes not match');
    else
        sAmt=cellstat(SD,'length');
        assert(isequal(cellstat(CST,'length'),sAmt),'number of spikes not match');
    end
end


%%% Get aligned segment data if input is not A data.
if dataType~=1
    X=spike_align(X,SD,srate,'window',pww,'chAssign',chAssign);%,'select',seleI);
end

%%% Ploting
sAmt=cellstat(X,'size',2); 
% Find the ms unit of x axis
% * 找到第一个存在的alignment为止，否则（如果全为空）采用预定点数。
idx=find(sAmt>0,1);
if isempty(idx)
    pntAmt=10;
else
    pntAmt=size(X{idx},1);
end
if bXaxis
    ax=linspace(pww(1),pww(2),pntAmt);
else
    ax=1:pntAmt;
end

% If need clustered spikes, define curve color.
if bCluster
    coltab=['b','r','g','y','m','c'];
    colnum=6;
end

%%% Initialize GUI controls
pgAmt=ceil(chAmt/16);
pgCh=cutseg([1,chAmt],16);

Fig=figure('Position',[100 0 1170 790],'ResizeFcn',@resizegui);
HpreB=uicontrol('Style','pushbutton','String','<--','Position',[50,800,50,25],...
    'Callback',@PreButton_CB);
HnextB=uicontrol('Style','pushbutton','String','-->','Position',[160,800,50,25],...
    'Callback',@NextButton_CB);
Hedit=uicontrol('Style','edit','Position',[110,800,40,25],...
    'Callback',@EditBox_CB);
HcurvePop=uicontrol('Style','popup','String',{'samples','mean'},'Position',[220,800,100,25],...
    'Callback',@CurvePop_CB);

set(Fig,'Name',sprintf('%d channels in total, %d pages',chAmt,pgAmt));
align(Fig,'Center','None');

% Draw 1st page
pgi=1;
curve=1; % 1=some of spikes, 2=all<<, 3=mean curves
drawpage(pgi);

if nargout>0
    varargout{1}=X;
end
%%%%%%%%% End of main


%%%%%%%%%%%%% Callback functions for the controls.
%%% pre button
    function PreButton_CB(hObject,eventdata)
       if pgi<=1
           pgi=1;
       else
           % Plot
            pgi=pgi-1;
           drawpage(pgi);
       end
    end
    
%%% next button
    function NextButton_CB(hObject,eventdata)
       if pgi>=pgAmt
           pgi=pgAmt;
       else
           % Plot
           pgi=pgi+1;
           drawpage(pgi);
       end
    end
    
%%% page number text input
    function EditBox_CB(hObject,eventdata)
        tp=get(Hedit,'String');
        pgi=str2double(tp);
        if pgi<1, pgi=1; end
        if pgi>pgAmt, pgi=pgAmt; end
        drawpage(pgi);
    end

%%% change the way curves are displayed
    function CurvePop_CB(hObject,eventdata)
        tp=get(HcurvePop,'Value');
        if tp~=curve
            curve=tp;
            drawpage(pgi);
        end
    end

%%% resize function
    function resizegui(hObject,eventdata)
        loc=get(Fig,'Position'); loc=loc(4);
        loc=loc-35;
        set(HpreB,'Position',[50,loc,50,25]);
        set(HnextB,'Position',[160,loc,50,25]);
        set(Hedit,'Position',[110,loc,40,25]);
        set(HcurvePop,'Position',[220,loc,40,25]);
    end


%%%%%%%%%%%%% Draw functions for one page.
    function drawpage(pgi)
        set(Hedit,'String',pgi);
        tp=pgCh(pgi,2)-pgCh(pgi,1)+1; % number of subplot in this page.
        [rn,cn]=subpnum(tp); % get column and row number for best subplot layout
        chI=(pgCh(pgi,1):pgCh(pgi,2));
        for drawi=1:length(chI);
            subplot(rn,cn,drawi,'align');
            cla;            hold on;
            % Find real channel index
            chi=chI(drawi);
            
            if sAmt(chi)>0
                %%% If need clustered spikes, determine color for each spike
                if bCluster
                    tp=length(CST{chi});
                    if tp~=sAmt(chi)
                        fprintf('SI length %d and spike number %d not match at %d\n',tp,sAmt(chi),chi);
                    end
                    
                    %%% Order the clusters. target: the un-grouped is black, others
                    % ordered in descending order of cluster size.
                    lb=reabylb(CST{chi});                    
                    % noisy (0) type
                    ni=find(lb.types==0); 
                    % number of noise spikes and index
                    if isempty(ni)
                        nn=0;
                        nidx=[];
                    else
                        nn=lb.typeAmt(ni);
                        nidx=lb.ids{ni};
                    end
                    % Remove noise type from lb:
                    if ~isempty(ni)
                        % remove the 0 type from list
                        lb.cAmt=lb.cAmt-1;
                        lb.types(ni)=[];
                        lb.typeAmt(ni)=[];
                        lb.ids(ni)=[];
                    end
                    
                    % Others
                    cmark=cell(lb.cAmt,1);
                    for k=1:lb.cAmt
                        cmark{k}=coltab(mod(k-1,colnum)+1);
                    end
                end
                
                %%% Decide which curves are to be shown, when total number
                %%% is limited.
                if bNumLimit
                    plotNum=min(limitNum,sAmt(chi));
                    % If clustered, ensure each cluster get a minimum
                    % number of curves.
                    if bCluster
                        tp=min(nn,floor(plotNum/lb.cAmt)); % number to show for each type on average
                        dispI=nidx(1:tp); % first include the noise spikes
                        for k=1:lb.cAmt % then append to valide spikes (ensure everyone have enough number counted in)
                            tp=min(lb.typeAmt(k),floor(plotNum/lb.cAmt));
                            dispI=[dispI;lb.ids{k}(1:tp)];
                        end
                    end
                else
                    plotNum=sAmt(chi);
                    dispI=1:plotNum;
                end
                
                %%%
                if bCluster
                    % Get data for display
                    temp=CST{chi}(dispI);
                    
                    % First draw the noises
                    I=temp==0; 
                    if bShow0 && sum(I)>0
                        plot(ax,X{chi}(:,I),'k');
                    end
                    % Effective ones depend on the form of display
                    if curve==1                        
                        for k=1:lb.cAmt
                            I=temp==lb.types(k);
                            if sum(I)>0
                                plot(ax,X{chi}(:,dispI(I)),cmark{k});
                            end
                        end
                    elseif curve==2
                        for k=1:lb.cAmt
                            mc=mean(X{chi}(:,lb.ids{k}),2);
                            plot(ax,mc,cmark{k},'LineWidth',3);
                        end
                    end
                else
                    if curve==1
                        plot(ax,X{chi}(:,1:plotNum),'b');
                    elseif curve==2
                        plot(ax,mean(X{chi}(:,1:plotNum),2),'b');
                    end
                end
            end
            
            %%%
            if bUseLabel
                tp=eleLabel{chAssign(chi)};
            else
                tp=['#',num2str(chAssign(chi))];
            end            
            str=sprintf('%s: %d',tp,sAmt(chi));
            if sAmt(chi)>0 && bCluster
                str=[str, sprintf('(%d)-%d',lb.cAmt,nn)];
            end
            if bNumLimit && sAmt(chi)>limitNum
                % when displayed number is limited, title in red.
                title(str,'Color','r');
            else
                title(str,'Color','k');
            end
            if bXaxis, xlabel('(ms)'); end
            hold off; axis tight;
        end
    end

end % end of main

%%%  obsolete
%                     for si=1:plotNum
%                         % when 1: show all (include outlier); 2: not show outlier and the current
%                         % neuron is not outlier.
%                         if bShow0 || (~bShow0 && ~strcmp(sColor{si},'k'))
%                             plot(ax,X{chi}(:,si),sColor{si});
%                         end
%                     end