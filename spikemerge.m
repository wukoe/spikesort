% Delete the duplicates from conduction signal
%   [newST,marklb,rmlb]=spikemerge(ST,SA,info)
% options:
% 'ext'=0.0005 (s)
% 'auto thres'=true
function [ST,marklb,rmlb]=spikemerge(ST,SA,info,varargin)
% Standard to find cluster sequence.
ext=0.001; % ! 1ms is necessary, testing shows some delay can >0.7ms.
% optSelect=2; % way to determine the representative of the cluster
param=struct();
param.bAutoThres=true; % whether use algorithm-determined adaptive threshold for pair number.
param.numPthr=1e-8; % for with auto threshold
minThres=4/60*info.TimeSpan; % static threshold if chosen to use.
flagUserTimeSpan=false;
bRmIntv=false; % whether remove the large "empty" segments of spike train
% spike time difference jitter limit.
param.DTJthr=0.15/1000;%(ms/1000) = 3 points at 20000HZ SR. <<<
param.bCScheck=false;

% User
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'ext'
                ext=pinfo{parai};
            case 'time span'
                param.timeSpan=pinfo{parai};
                flagUserTimeSpan=true;
            case 'auto thres'
                param.bAutoThres=pinfo{parai};
            otherwise
                error('unidentified options');
        end
    end
end

% Process
chAmt=length(ST);
% Get number of firing in each channel
sAmt=cellstat(ST,'length');
% Check number of SD SA
assert(isequal(sAmt,cellstat(SA,'length')),'merging: SD and SA number not equal');

if ~flagUserTimeSpan
    if bRmIntv
        % First remove the large "empty" segments of spike train
        % (most typical case: the inter-burst interval)
        for chi=1:chAmt
            % <<<
        end
    else
        param.timeSpan=max(cellstat(ST,'max'))-min(cellstat(ST,'min'));
    end
end


%%%%%%%%%%%%% Identify the spikes involved in conduction signal.
% Initialize the markers (Spike remain/delete, etc).
% cslb=cell(chAmt,1);% 被代表的通道号
marklb=cell(chAmt,1);% 各个spike是否被选为代表的标记
rmlb=cell(chAmt,1); % 各个spike是否被删除
for chi=1:chAmt
%     cslb{chi}=zeros(sAmt(chi),1);
    marklb{chi}=zeros(sAmt(chi),1);
    rmlb{chi}=false(sAmt(chi),1);
end
for chi=1:chAmt
    if ~param.bAutoThres
        if length(ST{chi})<minThres
            continue
        end
    else
        if length(ST{chi})<3
            continue
        end
    end
    
    %%% Find all tight-following spikes (under ext) in the channel pairs
    [Ichi,Isd,~,foI]=cspair(ST{chi},ST(chi+1:end),ext,sAmt(chi),param); % 不需要对全体进行扫描，只需要对后面的单元进行即可。因为与前面的已经进行过了。
    
    % 再次，若没有任何满足条件的conduction通道，放弃这个channel(chi).
    foAmt=length(foI);
    if foAmt==0 
        continue
    end
    
    %%% Add spikecs_check <<<,<<<<
    if param.bCScheck
        tx=X(:,chID(foI));
        tp=any(Ichi,2);
        tpI=false(size(Ichi));
        for k=1:size(Ichi,2)
            tpI(:,k)=tp;
        end
        [~,IR1,IR2]=spikecs_check(tx,msd(mi),info,'IW',{tpI});
        %
        foI=foI(IR1);        foAmt=sum(IR1);
        rnum=rnum(IR1);        Ichi=Ichi(:,IR1);        Isd=Isd(IR1);
        if foAmt<=1
            stat.fospknum{mi}=rnum;
            continue
        end
    end
    
    %%% Determine which to delete of each pair.
    % * Comparing rule:
    % if the amplitude difference is larger than 1.25 times of the
    % other, select by larger amplitude; else if the amplitude
    % difference is not large enough, select by more total spike number
    % (also only suitable for sorted data).
    % * 若完全根据chi中被关联的spike群体中已分配/未分配的比例决定是否将这些找到的视为单独的群体，
    % 面临detection缺陷的挑战。
    
    % 对于每个通道的“following”spike，找各自的幅值。
    ampchi=zeros(foAmt,1);
    ampsd=ampchi;
    for fi=1:foAmt
        % itself
        ampchi(fi)=mean(abs(SA{chi}(Ichi(:,fi))));
        % conducters
        ampsd(fi)=mean(abs(SA{foI(fi)+chi}(Isd{fi})));
    end
    
    % Apply the comparing rule.
    for fi=1:foAmt
%         if ampchi(fi)>ampsd(fi)*1.25 || ampsd(fi)>ampchi(fi)*1.25  %<<<<
            if ampchi(fi)>ampsd(fi)
                flagKeepchi=true;
            else
                flagKeepchi=false;
            end
%         else
%             %
%             if sAmt(chi)>sAmt(foI(fi)+chi)
%                 flagKeepchi=true;
%             else
%                 flagKeepchi=false;
%             end
%         end
        
        if  flagKeepchi %则保留chi
            % mark the "marker" spikes
            marklb{chi}(Ichi(:,fi))=1;
            % 标记要删除的sd(fi)
            rmlb{foI(fi)+chi}(Isd{fi})=true;
        else
            marklb{foI(fi)+chi}(Isd{fi})=1;
            rmlb{chi}(Ichi(:,fi))=true;
        end
    end
    fprintf('|');
end
fprintf('\n');

% Delete from the spike train & label of representative identity marker
for chi=1:chAmt
    ST{chi}(rmlb{chi})=[];
    marklb{chi}(rmlb{chi})=0;
end

end % of main