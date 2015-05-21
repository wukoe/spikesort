% Calculates the spike features - one channel
%   spkfea = spike_feature(A,varargin)
% A: each column as a spike, spikes stack up in row direction
function [spkfea,O]= spike_feature(spikes,varargin)
% Default
runMode=1; % 1=learning features, 2=use features

scales = 5;
outDim = 3; % 10 for wavelet, pca will change it to 3
spkn=size(spikes,2);

feature='dwt';
wlname='db2'; 
% use only 'db2' for db family
% 'coif2' for double negative peaks
% or 'bior1.3'
% otherwise in CWT: gaus

% User input
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'mode'
                runMode=pinfo{parai};
            case 'pararm'
                param=pinfo{parai};
            case 'feature'
                feature=pinfo{parai};
            case 'wltype'
                wlname=pinfo{parai};
            case 'dim'
                outDim=pinfo{parai};
            otherwise
                error('unidentified options');
        end
    end
end
if ~ismember(runMode,[1,2])
    error('runMode opt invalid');
end

%%% CALCULATES FEATURES
switch feature
    case 'dwt'
        [c,~]=wavedec(spikes(:,1),scales,wlname);
        slen = length(c);
        cc=zeros(slen,spkn); % discrete wavelet decomp coeffcients
        cc(:,1)=c;
        for spi=2:spkn
            [c,~]=wavedec(spikes(:,spi),scales,wlname);
%             cc(:,spi)=c(1:slen);
            cc(:,spi)=c;
        end
        
        % KS test for coefficient selection
        if runMode==1
            sd=zeros(slen,1);
            for pti=1:slen
                thr_dist = std(cc(pti,:)) * 3;
                thr_dist_min = mean(cc(pti,:)) - thr_dist;
                thr_dist_max = mean(cc(pti,:)) + thr_dist;
                % find the under-threshold
                aux = cc(pti,(cc(pti,:)>thr_dist_min & cc(pti,:)<thr_dist_max));
                
                if length(aux) > 10; % as standard to rule out the single model point
                    sd(pti)=test_ks(aux);
                else
                    sd(pti)=0;
                end
            end
            [~,idx]=sort(sd,'descend');
        else
            idx=param.selectFeaI;
        end
        spkfea=cc(idx(1:outDim),:)';
        
        O.selectFeaI=idx;
        
    case 'cwt'
        % Continuous wavelet coefficient method
        
    case 'pca'
        [~,cc,~] = princomp(spikes');        
        spkfea=cc(:,1:outDim);
    otherwise
        error('invalid feature option');
end

end