%  this program is for making the overlapping plot of spikes of one neuron
%  (channel) aligned to the peaks of spikes.
% Two types of data:
%   spkalignplot(fnX,SD,srate) uses the full signal X (by the corresponding X) and index format data SD.
%   spkalignplot(A) uses aligned data
%   A=spkalignplot(fnX,SD,srate) export the aligned data.
% Other variable parameters: 'NCI',nci - add clustering information by color;
% 'chI',I - plot only selected channels.
% 'window',w - specify the time window (ms) for spikes.
% 'limit',0 - no limit for the number of channels of each plot. (limited to
% 100 by default).
function varargout=spkalignplot(X,varargin)
%%% Parameters
% window for spike waveform
pww=[-1,1.5]; % pre-event & post-event window width (ms)

bCluster=false;
bShow0=1; % whether to display "outliers" (those clustered id as 0)
bNumLimit=true; % limit the maximum number of curves available in each plot
limitNum=100;

%%% Determine two different datatypes.
if iscell(X) % 直接输入波形的情况
    dataType=1;
    upara=varargin;    
    bXaxis=false;
    
elseif isa(X,'matlab.io.MatFile') % 输入文件句柄+SD，在这里截出分段
    dataType=2;
    SD=varargin{1}; srate=varargin{2};
    if nargin>3
        upara=varargin(3:end);
    else
        upara=[];
    end  
    bXaxis=true;

elseif isnumeric(X)
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
bChSpecify=false;
if ~isempty(upara)
    [pname,pinfo]=paramoption(upara{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'chI'
                chI=pinfo{parai};
                bChSpecify=true;
            case 'NCI' % neuron cluster identify
                CST=pinfo{parai};
                bCluster=true;
            case 'chAssign'
                chAssign=pinfo{parai};
            case 'window'
                pww=pinfo{parai};
                bXaxis=true;
            case 'limit'
                bNumLimit=pinfo{parai}==true;
            otherwise
                error('unidentified options');
        end
    end
end    

%%% Process
if dataType==1
    chAmt=length(X);
else % for both 2 and 3
    chAmt=length(SD);
end

% Whether channels are specified.
if bChSpecify
    dispAmt=length(chI);
    % Check the chI does not go across chAmt
    if bChSpecify
        if size(chI,2)>size(chI,1)
            chI=chI';
        end

        I=(chI<1) | (chI>chAmt);
        if any(I)
            error('invalid in chI');
        end
    end
else
    dispAmt=chAmt;
    chI=1:chAmt;
end

if ~exist('chAssign','var')
    chAssign=chI;
end


%%%%%%%%%%%%%% When the input is the file handle+SD, get aligned data first.
if dataType==1
    if bChSpecify
        X=X(chI);
    end
else
    if bChSpecify
        X=spike_align(X,SD(chI),srate,'window',pww,'chAssign',chAssign,'chI',chI); % here X is file handle input
    else
        X=spike_align(X,SD,srate,'window',pww);
    end
end

if bCluster && bChSpecify
    CST=CST(chI);
end


%%%%%%%%%%%%%%% Ploting
sAmt=zeros(1,dispAmt);
for chi=1:dispAmt
    sAmt(chi)=size(X{chi},2);
end

% Find the ms unit of x axis
% * 找到第一个长度不为0的Align为止，否则（如果全为空）采用预定点数。
for chi=1:dispAmt
    pntAmt=size(X{chi},1);
    if pntAmt>0, break; end
end
if pntAmt==0, pntAmt=10; end

if bXaxis
    ax=linspace(pww(1),pww(2),pntAmt);
else
    ax=1:pntAmt;
end
        
% If need clustered spikes, add color
if bCluster
    coltab=['b','r','g','y','m','c'];
    colnum=6;
end

%%%
Fig=figure('Name',num2str(dispAmt));
HpreB=uicontrol('Style','pushbutton','String','<--','Position',[315,180,70,25],...
    'Callback',@PreButton_CB);
HnextB=uicontrol('Style','pushbutton','String','-->','Position',[400,180,70,25],...
    'Callback',@NextButton_CB);
Htext=uicontrol('Style','text','Position',[315,180,70,25],...
    'Callback',@TextBox_CB);

[rn,cn]=subpnum(dispAmt); % get column and row number for best subplot layout
for chi=1:dispAmt
    subplot(rn,cn,chi,'align');
    cla;            hold on;
    
    if sAmt(chi)>0
        %%% If need clustered spikes, determine color for each spike
        if bCluster
            sColor=cell(sAmt(chi),1);
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
                sColor(lb.ids{cidx})={'k'};
                % remove the 0 type from list
                lb.cAmt=lb.cAmt-1;
                lb.types(cidx)=[];
                lb.typeAmt(cidx)=[];
                lb.ids(cidx)=[];
            end

            % others
            [~,lbOrder]=sort(lb.typeAmt,'descend');
            for k=1:lb.cAmt
                cidx=lbOrder(k);
                cmark=coltab(mod(k-1,colnum)+1);
                sColor(lb.ids{cidx})={cmark};
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
            for si=1:plotNum
                % when 1: show all (include outlier); 2: not show outlier and the current
                % neuron is not outlier.
                if bShow0 || (~bShow0 && ~strcmp(sColor{si},'k'))
                    plot(ax,X{chi}(:,si),sColor{si});
                end
            end
        else
            for si=1:plotNum
                plot(ax,X{chi}(:,si));
            end
        end
    end
    
    %%%
    str=['#',num2str(chI(chi)),': ',num2str(sAmt(chi))];
    if bNumLimit && sAmt(chi)>limitNum
        % when displayed number is limited, title in red.
        title(str,'Color','r');
    else
        title(str);
    end
    if bXaxis, xlabel('(ms)'); end
    hold off; axis tight; 
end

if nargout==1
    varargout{1}=X;
end


%%%%%%%%%%%%% Callback functions for the controls.
    %%% pre button
    function PreButton_CB(source,data)
    set(Fig,'Name','pre');
    end
    
    %%% next button
    function NextButton_CB(source,data)
    set(Fig,'Name','new');
    end
    
    %%% text input
    function TextBox_CB()
    end

end