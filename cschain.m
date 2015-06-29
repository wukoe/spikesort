% See also cspair
function [Isd]=cschain(SD,ext)%,bOrder)
cha=length(SD);
Isd=cell(cha,1);

% First find cs pairs between 1 and 2
[Iw,I]=cspair(SD{1},SD(2),ext);
Isd{1}=find(Iw);
Isd{2}=find(I{1});

if cha>2
    % For the next node, find cs pair with all existing channels and link
    % it with one having maximum number of pairs with it.
    for chi=3:cha
        % Analyze only linked spikes
        temp=cell(chi-1,1);
        for k=1:chi-1
            temp{k}=SD{k}(Isd{k});
        end
        [Iw,Iexist,num]=cspair(SD{chi},temp,ext);
        [~,idx]=max(num);
        
        Isd{chi}=find(Iw(:,idx));
        
        % trim the existing channels\
        seleI=Iexist{idx};
        for k=1:chi-1
            Isd{k}=Isd{k}(seleI);
        end
    end    
end

end