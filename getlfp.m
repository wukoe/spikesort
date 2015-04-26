% Get LFP from .mcd file

function getlfp(fileName)
bServer=true;
paras=struct();
paras.lfpFiltFreq=[4,300];
paras.lfpSrate=1000;

%%%
dinfo=getMCDinfo([fileName,'.mcd'],bServer);

% Get filter
[lfblp,lfalp]=butter(4,paras.lfpFiltFreq(2)/(dinfo.srate/2),'low');
% [lfbhp,lfahp]=butter(4,paras.lfpFiltFreq(1)/(dinfo.srate/2),'high');
% Process
paras.lfpDSratio=paras.lfpSrate/dinfo.srate;
T=(1:dinfo.ptsAmt)'/dinfo.srate; % measured in (s)    

%%%
% 先对T进行降采样，所得降采样索引在后面对数据进行时重复使用，免去调用ds()
[T,dsIndex]=ds(T,paras.lfpDSratio);
newpts=length(dsIndex); % length after DS

% Set up input
[nsresult, infile1] = ns_OpenFile([fileName,'.mcd']);
if nsresult==-1,     error('open .mcd file error'); end    

%%% Do the filtering
X=zeros(newpts,dinfo.rawchAmt);
for chi=1:dinfo.rawchAmt        
    [~,~,temp]=ns_GetAnalogData(infile1,dinfo.dataChEntity(chi),1,dinfo.ptsAmt);
    temp=filtfilt(lfblp,lfalp,temp);
%     temp=filtfilt(lfbhp,lfahp,temp);
    X(:,chi)=temp(dsIndex);
    fprintf('|');
end
fprintf('\n');

%%% Save LFP
save([fileName,'_lfp'],'X','T','paras');

end