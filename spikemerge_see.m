% analyze the info from spikemerge.m
function [Fch]=spikemerge_see(seqlist,seqCS,idxS,X,SD,info)
% Find the space range of each cluster (- the real "range" by 查看所有通道中出现aligned
% spike的情况。如果两个cluster其实是同一个，这个方法可以更准确地把这个找出来。
%%% Setting
xcha=120;
Rthr=2; % partially tested


end

% seqAmt=length(seqlist);
% sd=cell(xcha,1);
% Fch=false(xcha,seqAmt);
% for seqi=1:seqAmt
%     % Align the signal in all X channels to representative SD of cluster.
%     chi=seqlist{seqi}(idxS(seqi));
%     temp=SD{chi}(seqCS{seqi}(:,idxS(seqi)));
%     for chi=1:xcha
%         sd{chi}=temp;
%     end
%     A=spike_align(X,sd,info.srate,'window',[-1 1.5]);
%     
%     % Get ratio of spike peaks to noise
%     alen=size(A{1},1);
%     mc=zeros(alen,xcha); v=zeros(1,xcha);
%     for chi=1:xcha
%         mc(:,chi)=mean(A{chi},2); % mean of all curves.
%         v(chi)=max(std(A{chi},[],2)); % standard curves.
%     end    
%     ap=max(mc)-min(mc); % amplitude of mean curves
%     % Ratio
%     R=ap./v;
%     
%     %%% related channels
%     Fch(:,seqi)=R; 
% end