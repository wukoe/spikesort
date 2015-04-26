% Find repeat pattern
%  [allelist,allecount]=findrepeat(S,lenOpt,varargin)
% options part ---
% 'ms': mini support
% 'mc': mini confidence
function [allelist,allecount]=findrepeat(S,lenOpt,varargin)
Lthres=30;
mcThres=0.9;
singleThres=0.05;

% User input
if isempty(lenOpt)
    bSpecificLen=false;
else
    bSpecificLen=true;
end
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'ms' % minimum support
                Lthres=pinfo{parai};
            case 'mc' % minimum confidence
                mcThres=pinfo{parai};
            case 'single'
                singleThres=pinfo{parai};
            otherwise
                error('unidentified options');
        end
    end
end

% Proc
sAmt=length(S); % seq number
sLen=cellstat(S,'length'); % length of each seq

%%% Get list of type of elements (=L1)
elist=[];  ecount=[];
for si=1:sAmt
    for m=1:sLen(si)
        idx=find(elist==S{si}(m),1);
        if isempty(idx)
            elist=[elist;S{si}(m)];
            ecount=[ecount;1];
        else
            ecount(idx)=ecount(idx)+1;
        end
    end
end
eAmt=length(elist);
if bSpecificLen && lenOpt==1
    allelist=elist; allecount=ecount; 
    return; 
end

%%% Get L2 seqs
% Use table, to count occurance (in input seqs) of all possible combination of 2 elements
R=zeros(eAmt,eAmt);
for si=1:sAmt
    sq=S{si}; sql=length(sq);
    for m=1:sql-1
        % transform the element of m to elist index
        mi=find(elist==sq(m),1);
        for n=m+1:sql
            % transform the element of n to elist index
            ni=find(elist==sq(n),1);
            R(mi,ni)=R(mi,ni)+1;
        end
    end
end

% Get effective combinations. - update items list, now list contains L2
% items
% Filter R for combination with enough frequency - by number here.
[r,c,~]=find(R>Lthres);
elAmt=length(r); % elemebt
R2=zeros(elAmt,2);
ecount=zeros(elAmt,1);
for m=1:elAmt
    R2(m,:)=[elist(r(m)),elist(c(m))];
    ecount(m)=R(r(m),c(m));
end
elist=R2; clear('R','R2');
if bSpecificLen && lenOpt==2    
    allelist=elist; allecount=ecount; 
    return; 
end

%%%%%%%%%%% Extending! with Apriori rule
eAmt=size(elist,1);
temp=cell(eAmt,1);
for k=1:eAmt
    temp{k}=elist(k,:);
end
temp(ecount<Lthres)=[];
elist=temp;

% items of all length.
allelist=cell(2,1); % make allelist{1} empty for L=1.
allecount=cell(2,1);

L=2; % length of member in pool.
flagContinue=true;
while flagContinue
    %%% Find L+1 candidates by L-joining
    % * L-joining的方法是保证有L-1个相同的情况下往一头增加一位。
    eAmt=length(elist);
    
    newelist=cell(0); newecount=zeros(0,1);
    listcount=0;
    for m=1:eAmt
        tar=elist{m}(2:end);%原序列的2:end element作为target
        for n=1:eAmt
            % if 1st in n = last in m
            if elist{n}(1:end-1)==tar
                % Determine one candidate
                candy=[elist{m},elist{n}(end)];
                % Check if candidate is already in new list <<< necessary?
                
                % 使用Apriori法则，Check if every subset is a member of L -
                % 目标长为L+1，只需找长L的subset即可：共L+1个，只用挨个去掉单个元素即可。
                % 又由于已经有了去掉第一和最后元素的，所以只需要找去掉2~L+1-1的即可。
                flagNotAp=false; % flag: not qualified Apriori property.
                for k=2:L
                    tp=candy; tp(k)=[];
                    if isempty(cellfind(elist,tp))
                        flagNotAp=true;
                        break;
                    end
                end
                if flagNotAp % move to next candidate in n-loop
                    break; 
                end
                
                % Find how many times candidate occurs in raw seqs
                I=scan_seq_pattern(S,candy);
                N=sum(I);
                
                % Check if N is above "minimum support threshold".                
                if N>=Lthres % <<<<<<<<
                    % put it into new list
                    listcount=listcount+1;
                    newelist{listcount}=candy;
                    newecount(listcount,1)=N;
                end
            end
        end
    end
    
    %%%
    % Save this Level of items
    allelist{L}=elist;
    allecount{L}=ecount;
    
    % Judge whether to end loop based on number of items
    if listcount==0
        % If no L+1 item is found, end the search.
        flagContinue=false;
    else        
        if listcount<=2
            % If less than 3 L+1 items, also end the search, since even if L==2,
            % if less than 3 in L=3, there's no possibility any L=4 items could be
            % made.
            % Save next Level of items,and exit.
            allelist{L+1}=newelist;
            allecount{L+1}=newecount;
            flagContinue=false;
        else            
            % Update
            elist=newelist;
            ecount=newecount;
        end
        L=L+1;
    end
end

%%%%%%%%%% Association stage
% we use association as method to collapse subset item to larger item.
% * 从L到L+1的聚合按这个规则：假设L的m为L+1中n1,n2两者的子集，并且将进行合并。
% n1,n2任何一者的数量必须达到m的最低比例（5%）才能参与分配。
% 若所有合格参与者的数量的和达到m的高位-minimum confidence threshold(90%), 则m
% 被n1,n2取代——自己被删除，剩余的计数（10%或更少)被累计入n1,n2中,分配比例按照n1,n2的原数量计算；
% 若达不到MCT,则保留m独立存在，但同时要将L+1,n的数量从m中减去！。
% n1,n2,n3...数量更多时类推。
% 注意！！ 这个顺序必须从高到低进行，即若最长为L=5的话，第一进行L=4->5,再L=3->4...
% 因为若非最后一层，完全可能出现"n1+n2+...>m"的情况！
if L>2
    for li=L-1:-1:2
        LayerItemAmt=length(allelist{li});
        %%% 分配网络
        WR=cell(LayerItemAmt,1); %remove-删除L的item
        WS=cell(LayerItemAmt,1); %subtract-只减少L item的count数    
        for m=1:LayerItemAmt
            % Identify all items in L+1 which m is a subset of.
            I=scan_seq_pattern(allelist{li+1},allelist{li}{m});
            WR{m}=find(I);            
            % Get repeat number of m and all ni
            LPOItemRepeat=allecount{li+1}(WR{m});% 参与分配者的数量
            LItemRepeat=allecount{li}(m);
            % By these number, Filter out qualified items in L+1
            I=(LPOItemRepeat>=LItemRepeat*singleThres);
            WR{m}=WR{m}(I);
            
            % If sum of all qualified is over minimum confidence threshold,
            % mark as remaining WR(erase m). Otherwise mark as WS(subtract n from m).
            LPOItemRepeat=LPOItemRepeat(I);
            if sum(LPOItemRepeat)<LItemRepeat*mcThres
                WS{m}=WR{m};
                WR{m}=[];
            end
        end
        
        %%% Process the L items
        rmlist=false(LayerItemAmt,1);
        Nadd=zeros(length(allelist{li+1}),1);
        for m=1:LayerItemAmt
            tar=WR{m};            
            if ~isempty(tar)
                LItemRepeat=allecount{li}(m);                
                LPOItemRepeat=allecount{li+1}(tar);
                %
                if LItemRepeat>sum(LPOItemRepeat)
                    N=round((LItemRepeat-sum(LPOItemRepeat)) * LPOItemRepeat/sum(LPOItemRepeat));
                    Nadd(tar)=Nadd(tar)+N;
                end
                rmlist(m)=true;
            elseif ~isempty(WS{m})
                LPOItemRepeat=allecount{li+1}(WS{m});
                allecount{li}(m)=allecount{li}(m)-sum(LPOItemRepeat);
            end            
        end        
        allecount{li+1}=allecount{li+1}+Nadd;
        
        % Delete items in L specified by WR
        allelist{li}(rmlist)=[];
        allecount{li}(rmlist)=[];
        % Combine L+1 to L
        allelist{li}={allelist{li}{:},allelist{li+1}{:}};
        allecount{li}=[allecount{li};allecount{li+1}];
        allelist(li+1)=[]; allecount(li+1)=[];
    end
    
    rmlist=(allecount{2}<Lthres);
    allelist{2}(rmlist)=[]; allecount{2}(rmlist)=[];
    
%     %%% 重新按长度分配层级。
%     sl=cellstat(allelist{2},'length');
%     temp=cell(L,1); temp2=temp;
%     for li=2:L
%         I=(sl==li);
%         temp{li}=allelist{2}(I);
%         temp2{li}=allecount{2}(I);
%     end
%     allelist=temp; allecount=temp2;
    
%     %%% 每一层的整合分配决定
%     for li=2:L-1
%         WR{li}=cell(LayerItemAmt(li),1); WS{li}=cell(LayerItemAmt(li),1);
%         for m=1:LayerItemAmt(li)
%             % Identify all items in L+1 which m is a subset of.
%             I=scan_seq_pattern(allelist{li+1},allelist{li}{m});
%             WR{li}{m}=find(I);            
%             % From these, Filter out qualified items in L+1
%             LPOItemRepeat=allecount{li+1}(WR{li}{m});
%             LItemRepeat=allecount{li}(m);
%             I=(LPOItemRepeat>=LItemRepeat*singleThres);
%             WR{li}{m}=WR{li}{m}(I);
%             
%             % If sum of all qualified is over minimum confidence threshold,
%             % mark as remaining WR(erase m). Otherwise mark as WS(subtract n from m).
%             LPOItemRepeat=LPOItemRepeat(I);
%             if sum(LPOItemRepeat)<LItemRepeat*mcThres
%                 WS{li}{m}=WR{li}{m};
%                 WR{li}(m)=[];
%                 LayerItemAmt(li)= LayerItemAmt(li)-1;
%             end            
%         end
%     end
%     
%     %%% 实施整合
%     for li=L-1:-1:2
%         rmlist=[];
%         for m=1:LayerItemAmt(li)
%             tar=WR{li}{m};
%             if ~isempty(tar)
%                 LItemRepeat=allecount{li}(m);
%                 % 参与分配者的数量
%                 LPOItemRepeat=allecount{li+1}(tar);
%                 % 
%                 if LItemRepeat>sum(LPOItemRepeat)
%                 N=(LItemRepeat-sum(LPOItemRepeat)) * LPOItemRepeat/sum(LPOItemRepeat);
%                 allecount{li+1}(tar)=allecount{li+1}(tar)+N;
%                 end
%                 %
%                 rmlist=[rmlist;m];
%             elseif ~isempty(WS{li}{m})
%                 tar=WS{li}{m};
%                 LPOItemRepeat=allecount{li+1}(tar);
%                 allecount{li}(m)=allecount{li}(m)-sum(LPOItemRepeat);
%             end
%         end
%         
%         % Remove
%         allelist{li}(rmlist)=[];
%         allecount{li}(rmlist)=[];
%     end
end


%%%
if bSpecificLen
    if lenOpt>L
        allelist=[]; allecount=[];
    else
        allelist=allelist{lenOpt}; allecount=allecount{lenOpt};
    end
% else
%     allelist=allelist{2}; allecount=allecount{2};
end

end