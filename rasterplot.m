% raster plot of spike trains
%  rasterplot(ST)   ST is time index format or 
%  rasterplot(SD)   SD is logical format
%  rasterplot(SD,'time',[start_time, end_time]) will add time axes
%  ('group',CG) add color to neurons belonging to different channels
%  defined by CG.
function rasterplot(SD,varargin)
%%% User input
% Differentiate 2 types of data format
% pre-fetch the chAmt
if iscell(SD)
    % In case of time index format    
    bIdx=true;
    chAmt=length(SD);    
else
    % In case of logical indication format    
    bIdx=false;
    if ~islogical(SD)
        error('not the logical form expected');
    end
    [ptsAmt,chAmt]=size(SD);
end

% Specified options
% Default
bTime=false; % be time info.
bGroup=false; % whether group the channels 
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'time'
                bTime=true;
                time=pinfo{parai};
                
            case 'group'
                bGroup=true;
                gp=pinfo{parai}; % group information is here
                if length(gp)~=chAmt, error('grouping info length not match to channels'); end 
                
                % initiate color list for draw;
                coltab=['b','r','g'];%,'y','m','c'];
                colnum=3;
                
            case 'dt'
                dispTime=pinfo{parai};
                
            otherwise
                error('unidentified options');
        end
    end
end

%%% Process
if ~bTime
    time=[min(cellstat(SD,'min')),max(cellstat(SD,'max'))];
end
totalTime=time(end)-time(1);
if ~exist('dispTime','var')
    if totalTime<=30
        dispTime=totalTime;
    elseif totalTime>30 && totalTime<=300
        dispTime=30;
    else
        dispTime=totalTime/10;
    end
end

if ~bIdx
    if bTime
        srate=ptsAmt/(totalTime/1000);
    else
        srate=1;
    end
end

%%% When grouping is needed, color differently.
% Assign colors
colstr=cell(chAmt,1);
if bGroup% if grouping is needed
    gpcount=1; 
    lastChGp=gp(1); % group identity of last channel, initiated as the 1st one.
    for chi=1:chAmt
        if gp(chi)~=lastChGp
            gpcount=gpcount+1;
            lastChGp=gp(chi);
        end
        idx=mod(gpcount-1,colnum)+1;
        colstr{chi}=[coltab(idx),'-'];
    end
else % all black
    for chi=1:chAmt
        colstr{chi}='k-';
    end
end

%%% Initialize GUI controls
pgAmt=ceil(totalTime/dispTime);
pgTime=cutseg([time(1),time(end)],dispTime,'common');

Fig=figure(); %('Position',[50 50 1120 840]);
HpreB=uicontrol('Style','pushbutton','String','<--','Position',[10,0,50,25],...
    'Callback',@PreButton_CB);
HnextB=uicontrol('Style','pushbutton','String','-->','Position',[100,0,50,25],...
    'Callback',@NextButton_CB);
HpageE=uicontrol('Style','edit','Position',[60,0,40,25],'Callback',@PageEdit_CB);
HrangeE=uicontrol('Style','edit','Position',[160,0,40,25],'Callback',@RangeEdit_CB);

set(Fig,'Name',sprintf('%0.1g seconds in total, %d pages',totalTime,pgAmt));
set(HrangeE,'String',dispTime);
align(Fig,'Center','None');

% Draw 1st page
pgi=1;
drawpage(pgi);
title('raster plot');
if bTime, xlabel('(s)'); end
ylabel('channels');
        
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
    function PageEdit_CB(hObject,eventdata)
        tp=get(HpageE,'String');
        pgi=str2double(tp);
        if pgi<1, pgi=1; end
        if pgi>pgAmt, pgi=pgAmt; end
        drawpage(pgi);
    end

% %%% resize function
    function RangeEdit_CB(hObject,eventdata)
        tp=get(HrangeE,'String'); dispTime=str2double(tp);
        pgAmt=ceil(totalTime/dispTime);
        pgTime=cutseg([time(1),time(end)],dispTime,'common');
        pgi=1;
        drawpage(pgi);
    end


%%%%%%%%%%%%% Draw functions for one page.
    function drawpage(pgi)
        set(HpageE,'String',pgi);
        SDtp=sdcut(SD,pgTime(pgi,:));
        
        cla;
        hold on
        if bIdx
            for chi=1:chAmt
                I=SDtp{chi};
                spikeNum=length(I);
                if size(I,2)<spikeNum, I=I'; end
                
                % use chi as draw height of each neuron (1st to last: top 2 bottom)
                plot([I;I],[chi+zeros(1,spikeNum);chi+ones(1,spikeNum)],colstr{chi});
            end
            axis([pgTime(pgi,1),pgTime(pgi,2),0,chAmt+1]);
            
        else
            for chi=1:chAmt
                I=find(SDtp(:,chi))';
                % transfer into time based on sampling rate
                I=I/srate;
                spikeNum=length(I);
                % use chi as draw height of each neuron
                plot([I;I],[chi+zeros(1,spikeNum);chi+ones(1,spikeNum)],colstr{chi});
            end
        end
        
        hold off
    end

end