%  this program is for making the overlapping plot of spikes of one neuron
%  (channel) aligned to the peaks of spikes.
%   spikealign(X,SD,srate,'plot') plot all channels
%   spikealign(X,SD,srate,'plot',SI) plot with clustering information
%   A=spikealign() return the window data
function A=spike_align(X,SD,srate,varargin)
% parameter default
preww=2.5; % pre-event window width (ms)
postww=5; % post-event

% options
% yalign in the code % at maximum (for positive spike) or minimum (for negative spike)
% =0 for align at positive-negative cross 0 point.

bPlot=false;
bCluster=false;

%%%%%%%%%%%%%%%%%% processing of data
[pntAmt,chAmt]=size(X);
if size(SD,2)~=chAmt,    error('size'); end

% whether select specific channels for display
nvarin=length(varargin);
if nvarin>0
    for k=1:nvarin
        if ischar(varargin{k})
            if strcmpi(varargin{k},'plot')
                bPlot=true;
            else                
                error('unknown input type: %s\n',varargin{k});
            end
        else
            CST=varargin{k};
            bCluster=true;
        end
    end    
end

%
sAmt=sum(SD);

% turn ww in (ms) to ww in (pts)
preww=fix(preww/1000*srate);
postww=fix(postww/1000*srate);

% output init
A=cell(chAmt,1);
for chi=1:chAmt
    A{chi}=zeros(preww+postww+1,sAmt(chi));
end

%%% If need clustered spikes, add color
if bCluster
    coltab=['b','r','g','y','m','c'];
    colnum=6;
end


%%%%%%%%%%%%%%%%%%%%%%%%
if bPlot
    clf; 
    [rn,cn]=subpnum(chAmt); % get column and row number for best subplot layout
end

for chi=1:chAmt    
    SI=find(SD(:,chi));
    if sAmt(chi)>0
        if bPlot
            subplot(rn,cn,chi,'align');
            cla;            hold on;
            
            %%% If need clustered spikes, determine color for each spike
            if bCluster
                sColor=cell(sAmt(chi),1);
                if length(CST{chi})~=sAmt(chi)
                    error('SI length and spike number match not');
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
        end
        
        %%%% If need to insert the maximum number limit, it is here
        rmlist=[];
        for si=1:sAmt(chi)
            % data segment index
            temp=beinrange([SI(si)-preww,SI(si)+postww],1,pntAmt);          
            dslow=temp(1); dshigh=temp(2);
            
%             % 1/N align at maximum
%             temp=X(dslow:dshigh,chi);
%             temp=temp-X(I(si),chi);

            % 2/N align at mean
            temp=X(dslow:dshigh,chi);
            temp=temp-mean(temp);
            
%             % 3/N align at mean-crossing point
%             temp=temp-mean(temp);
%             temp=(temp<0);
%             idx=find(temp(I(si)-dslow:end),1);
%             dslow=I(si)+idx-fix(ww/2);
%             dshigh=I(si)+idx+fix(ww/2);
%             temp=X(dslow:dshigh,chi);
            
%             % 4/N
%             vmax=X(I(si),chi);
%             [vmin,idx]=min(X(I(si):I(si)+120,chi));
%             idx=I(si)+round(idx/2);
%             temp=X(idx-ww:idx+ww,chi);
%             temp=temp-(vmax+vmin)/2;
            
            % save data in output
            if length(temp)<preww+postww+1
                % discard the spike at ends of data segment that are not as long
                rmlist=[rmlist,si];
                
            %%%%%%%%%% If want to discard more stuff, the code can be here
            %%%%%%%%%% with extra "elseif" command.
            else
                A{chi}(:,si)=temp;
            end
            
            if bPlot
                % find the ms unit of x axis
                ax=dslow-SI(si):dshigh-SI(si);
                ax=ax/srate*1000;
                
                if bCluster
                    plot(ax,temp,sColor{si});
                else
                    plot(ax,temp);
                end
            end
        end
        A{chi}(:,rmlist)=[];
        
        if bPlot
            title(['#',num2str(chi),': ',num2str(sAmt(chi))]);
            hold off; axis tight; 
        end
    end
end