% Restore the raw SD in _s.mat file.
% For processing missing variables due to sorting update in old files.
%   procoldfile(fn)
function procoldfile(fn)
% Load
fnS=[fn,'_s'];
load(fnS);

% Restore the raw SD in _s.mat file.
SD=rawSD; SA=rawSA;
info.chAmt=info.rawchAmt;
info.chID=(1:info.rawchAmt)';

% Save
save(fnS,'-append','SD','SA','info');
end