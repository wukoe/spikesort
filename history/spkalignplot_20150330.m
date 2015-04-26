% Program for making the plot of spikes of one neuron/channel aligned to time of spikes(peak). 
% Use button and text to navigate.
% Accept 2 types of data:
%   spkalignplot(fnX,SD,srate) uses full signal X and index format data spike SD.
% X can be either data in memory or MatFile handle. In this case the actual
% number of channels shown depends on SD.
%   spkalignplot(A) uses aligned data
% 
%   A=spkalignplot() export the aligned data.
% Other varargin:
% 'NCI', (neuron cluster identity) - add clustering information using color;
% 'select',I - plot only selected channels, chI is all about SD, not raw
% channel indexes.
% 'window',w - specify the time window (ms) for spikes.
% 'limit',0 - no limit for the number of spikes of each plot. non-0 specify the number (limited to
% 100 by default).
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
bChSelect=false;
if ~isempty(upara)
    [pname,pinfo]=paramoption(upara{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'select'
                seleI=pinfo{parai};
                bChSelect=true;
            case 'NCI' % neuron cluster identity
                CST=pinfo{parai};
                bCluster=true;
            case 'chAssign'
                chAssign=pinfo{parai};
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
            otherwise
                error('unidentified options');
        end
    end
end    

%%% Process
% Channel number
if dataType==1
    chAmt=length(X);
else % for both 2 and 3
    chAmt=length(SD);
    % If ch num in SD is different from X, 'chAssign' is mandatory.
    if chAmt~=size(X,2)
        
    end
end

% Whether channels are specified.
if bChSelect
    seleChAmt=length(seleI);
    % Check the chI does not go across chAmt
    if size(seleI,2)>size(seleI,1)
        seleI=seleI';
    end

    I=(seleI<1) | (seleI>chAmt);
    if any(I)
        error('invalid in chI');
    end
else
    seleChAmt=chAmt;
    seleI=1:chAmt;
end

if ~exist('chAssign','var')
    chAssign=seleI;
end


%%%%%%%%%%%%%% Get aligned segment data.
if dataType==1
    if bChSelect
        X=X(seleI);
    end
else
    if bChSelect
        X=spike_align(X,SD(seleI),srate,'window',pww,'chAssign',chAssign);%,'select',seleI); % here X is file handle input
    else
        X=spike_align(X,SD,srate,'window',pww);
    end
end

if bCluster && bChSelect
    CST=CST(seleI);
end


%%%%%%%%%%%%%%% Ploting
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
pgAmt=ceil(seleChAmt/16);
pgCh=cutseg([1,seleChAmt],16);

Fig=figure('Position',[50 50 1120 840],'ResizeFcn',@resizegui);
HpreB=uicontrol('Style','pushbutton','String','<--','Position',[50,800,50,25],...
    'Callback',@PreButton_CB);
HnextB=uicontrol('Style','pushbutton','String','-->','Position',[160,800,50,25],...
    'Callback',@NextButton_CB);
Hedit=uicontrol('Style','edit','Position',[110,800,40,25],...
    'Callback',@EditBox_CB);

set(Fig,'Name',sprintf('%d channels in total, %d pages',seleChAmt,pgAmt));
align(Fig,'Center','None');

% Draw 1st page
pgi=1;
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
    
%%% text input
    function EditBox_CB(hObject,eventdata)
        tp=get(Hedit,'String');
        pgi=str2double(tp);
        if pgi<1, pgi=1; end
        if pgi>pgAmt, pgi=pgAmt; end
        drawpage(pgi);
    end

%%% resize function
    function resizegui(hObject,eventdata)
        loc=get(Fig,'Position'); loc=loc(4);
        loc=loc-35;
        set(HpreB,'Position',[50,loc,50,25]);
        set(HnextB,'Position',[160,loc,50,25]);
        set(Hedit,'Position',[110,loc,40,25]);
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
%                     sColor=cell(sAmt(chi),1);
                    tp=length(CST{chi});
                    if tp~=sAmt(chi)
                        fprintf('SI length %d and spike number %d not match at %d\n',tp,sAmt(chi),chi);
                    end
                    
                    %%% Order the clusters. target: the un-grouped is black, others
                    % ordered in descending order of cluster size.
                    lb=reabylb(CST{chi});
                    
                    % 0 type
                    cidx=find(lb.types==0);
                    if ~isempty(cidx)
%                         sColor(lb.ids{cidx})={'k'};
                        % remove the 0 type from list
                        lb.cAmt=lb.cAmt-1;
                        lb.types(cidx)=[];
                        lb.typeAmt(cidx)=[];
                        lb.ids(cidx)=[];
                    end
                    
                    % others
%                     [~,lbOrder]=sort(lb.typeAmt,'descend');
                    cmark=cell(lb.cAmt,1);
                    for k=1:lb.cAmt
                        cmark{k}=coltab(mod(k-1,colnum)+1);
%                         cidx=lbOrder(k);
%                         sColor(lb.ids{cidx})={cmark{k}};
                    end
                end
                
                %%% Maximum number limit
                if bNumLimit
                    plotNum=min(limitNum,sAmt(chi));
                else
                    plotNum=sAmt(chi);
                end
                
                %%%
                if bCluster
                    %1/2
%                     for si=1:plotNum
%                         % when 1: show all (include outlier); 2: not show outlier and the current
%                         % neuron is not outlier.
%                         if bShow0 || (~bShow0 && ~strcmp(sColor{si},'k'))
%                             plot(ax,X{chi}(:,si),sColor{si});
%                         end
%                     end
                    %2/2                    
                    temp=CST{chi}(1:plotNum);
                    I=temp==0; 
                    if bShow0 && sum(I)>0
                        plot(ax,X{chi}(:,I),'k');
                    end
                    for k=1:lb.cAmt
                        I=temp==lb.types(k);
                        if sum(I)>0
                            plot(ax,X{chi}(:,I),cmark{k});
                        end
                    end                    
                    %e/2
                else
                    plot(ax,X{chi}(:,1:plotNum),'b');
                end
            end
            
            %%%
            str=['#',num2str(seleI(chi)),': ',num2str(sAmt(chi))];
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

end