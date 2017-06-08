clear all
close all
codeFolder = '/Users/kohler/code';
addpath(genpath([codeFolder,'/git/schlegel/matlab_lib']));
addpath(genpath([codeFolder,'/git/mrC']));
addpath(genpath([codeFolder,'/git/MRI/matlab']));

loadData = true;
loadRingData = false;
loadFoveaData = true;
plotROIs = false;
%plotSNR = false;
suffix = 'vr.sc.dt';

topFolder = '/Volumes/Denali_MRI/kohler/fMRI_EXP/MULTIFOVEA';%'/Users/kohler/Dropbox/PRESENTATIONS/2016/FENS_2016'
figFolder = '/Volumes/Denali_MRI/kohler/fMRI_EXP/MULTIFOVEA/figures';

subjFolders = subfolders([topFolder,'/201*'],1);

benoitFolder = '/Volumes/Denali_4D2/kohler/fMRI_EXP/CS_DISP';
subjFolders = [subjFolders;subfolders([benoitFolder,'/201*'],1)];
benoitSubj = cell2mat(cellfun(@(x) ~isempty(strfind(x,'CS_DISP')),subjFolders,'uni',false));

subCount = [length(subjFolders),sum(~benoitSubj),sum(benoitSubj)]; % all, mine, benoit

saveName = [figFolder,'/',num2str(subCount(1)),'sub_',suffix,'.mat'];
ringSave = [figFolder,'/ring_',num2str(subCount(1)),'sub_',suffix,'.mat'];
foveaSave = [figFolder,'/fovea_',num2str(subCount(1)),'sub_',suffix,'.mat'];

nCycles = 10;
nHarm = 5;
nTR = 120;
whichSNR = 2; % which SNR to use [1 = noisefloor within condition, 2=noisefloor across conditions]
whichHarm = 1; % which harmonic to use, usually 1
setenv('DYLD_LIBRARY_PATH','')

roiSelection = {'V1' 'V2' 'V3' 'hV4' 'VO1' 'VO2','PHC1','PHC2',...
    'TO1','LO1','LO2','V3A','V3B','IPS0','IPS1','IPS2','IPS3','SPL1'}; % exclude: 'IPS4','IPS5'

wangROIs = {'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'}; % Wang ROIs, in order

myConds = {'cont','mofo','disp'};
benoitConds = {'C1vsC2','C2vsC3','C1vsC3'};

if ~loadData
    roiFiles = '*wangatlas_al_mask.nii.gz';
    roiLabel = [cellfun(@(x) [x,'-L'],wangROIs,'uni',false),cellfun(@(x) [x,'-R'],wangROIs,'uni',false)];
    for s=1:subCount(1)
        curROI = subfiles([subjFolders{s},'/ROIs/',roiFiles],1);
        disp(['Reading ',subjFolders{s},' ...']);
        if ~benoitSubj(s)
            curConds = myConds;
        else
            curConds = benoitConds;
        end
        curROI = subfiles([subjFolders{s},'/ROIs/',roiFiles],1);
        for c = 1:length(curConds)
            curFiles = subfiles([subjFolders{s},'/run*',curConds{c},'*',suffix,'.nii.gz'],1);
            if curFiles{1}
                [roiData.(curConds{c})(s,:)] = mriRoiFFT(curFiles,curROI,roiLabel,nCycles,roiSelection);
            else
            end
        end
    end
    save(saveName,'roiData','-v7.3');
else
    load(saveName,'roiData')
end

%% RING ROIS
if ~loadRingData
    ringIncr = 0.25;
    ringSize = .5;
    ringMax = 6;
    ringList = (ringSize/2):ringIncr:ringMax; % list of centers
    wangFiles = '*wangatlas_al_mask.nii.gz';
    ringFiles = '*eccen.V1-V3model_al_mask.nii.gz';
    roiLabel = [cellfun(@(x) [x,'-L'],wangROIs,'uni',false),cellfun(@(x) [x,'-R'],wangROIs,'uni',false)];
    evcIdx{1} = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'V1')),roiLabel,'uni',false)));
    evcIdx{2} = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'V2')),roiLabel,'uni',false)));
    evcIdx{3} = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'V3')),roiLabel,'uni',false)));
    abIdx = find(cell2mat(cellfun(@(x) sum([isempty(strfind(x,'V3A')),isempty(strfind(x,'V3B'))],2)==2,roiLabel,'uni',false)));
    evcIdx{3} = evcIdx{3}(ismember(evcIdx{3},abIdx));
    evcNames = {'V1','V2','V3'};
    
    for s = 1:subCount(1);
        if ~benoitSubj(s)
            curConds = myConds;
        else
            curConds = benoitConds;
        end
        disp(['Reading ring: ',subjFolders{s},' ...']);
        wangROI = subfiles([subjFolders{s},'/ROIs/',wangFiles],1);
        ringROI = subfiles([subjFolders{s},'/ROIs/',ringFiles],1);
        dataFiles = cellfun(@(x) subfiles([subjFolders{s},'/run*',x,'*',suffix,'.nii.gz'],1),curConds,'uni',false);
        wangData = mriReadBrainData(wangROI);
        ringData = mriReadBrainData(ringROI);
        for c=1:length(dataFiles)
            if dataFiles{c}{1}
                funcData = mean(mriReadBrainData(dataFiles{c}),5);
                
                for r = 0:length(ringList)
                    if r == 0
                        curMin = 0;
                        curMax = 6.25;
                    else
                        curMin = ringList(r)-(ringSize/2);
                        curMax = ringList(r)+(ringSize/2);
                    end
                    curNames = cellfun(@(x) [x,'ring',num2str(curMin),'-',num2str(curMax)],evcNames,'uni',false);
                    curRing = (ringData>curMin).*(ringData<=curMax);
                    curMask = cellfun(@(x) repmat(curRing.*ismember(wangData,x),1,1,1,nTR),evcIdx,'uni',false);
                    [tmpResults(:,r+1)] = cell2mat(arrayfun(@(x) ...
                        mriFFT(reshape(funcData(curMask{x}==1),length(find(curMask{x}==1))/nTR,nTR)',nCycles,nHarm,curNames{x}),...
                        1:3,'uni',false));
                end
                ringRoiData.(curConds{c})(s,:) = reshape(tmpResults,1,[]);
                clear tmpResults;
            else
            end 
        end
    end
    save(ringSave,'ringRoiData','-v7.3');
else
    load(ringSave,'ringRoiData')
end

%% FOVEA ROIS

foveaTypes = {'inner','border','outer'};

if ~loadFoveaData
    wangFiles = '*wangatlas_al_mask.nii.gz';
    foveaFiles = '*fovea_3mm_al_mask.nii.gz';
    roiLabel = [cellfun(@(x) [x,'-L'],wangROIs,'uni',false),cellfun(@(x) [x,'-R'],wangROIs,'uni',false)];
    
    evcNames = {'V1','V2','V3','V3A','V3B','hV4','VO1','LO1','LO2','TO1'};
    for z=1:length(evcNames)
        evcIdx{z} = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,evcNames{z})),roiLabel,'uni',false)));
    end
    abIdx = find(cell2mat(cellfun(@(x) sum([isempty(strfind(x,'V3A')),isempty(strfind(x,'V3B'))],2)==2,roiLabel,'uni',false)));
    evcIdx{3} = evcIdx{3}(ismember(evcIdx{3},abIdx));
    
    for s = 1:subCount(1);
        if ~benoitSubj(s)
            curConds = myConds;
        else
            curConds = benoitConds;
        end
        disp(['Reading fovea: ',subjFolders{s},' ...']);
        wangROI = subfiles([subjFolders{s},'/ROIs/',wangFiles],1);
        foveaROI = subfiles([subjFolders{s},'/ROIs/',foveaFiles],1);
        dataFiles = cellfun(@(x) subfiles([subjFolders{s},'/run*',x,'*',suffix,'.nii.gz'],1),curConds,'uni',false);
        wangData = mriReadBrainData(wangROI);
        foveaData = double(mriReadBrainData(foveaROI));
        for c=1:length(dataFiles)
            if dataFiles{c}{1}
                funcData = mean(mriReadBrainData(dataFiles{c}),5); % average over runs
                fovMask{1} = foveaData(:,:,:,2); % inner
                fovMask{2} = foveaData(:,:,:,4); % border
                fovMask{3} = sum(foveaData(:,:,:,[1,2]),4)==0; % everything else
                for f=1:length(foveaTypes)
                    curNames = cellfun(@(x) [x,'_',foveaTypes{f}],evcNames,'uni',false);
                    curFovea = fovMask{f};
                    curMask = cellfun(@(x) repmat(curFovea.*ismember(wangData,x),1,1,1,nTR),evcIdx,'uni',false);  
                    tmpResults(:,f) = cell2mat(arrayfun(@(x) ...
                        mriFFT(reshape(funcData(curMask{x}==1),length(find(curMask{x}==1))/nTR,nTR)',nCycles,nHarm,curNames{x}),...
                        1:length(evcNames),'uni',false));
                end
                foveaRoiData.(curConds{c})(s,:) = reshape(tmpResults,[],1);
                clear tmpResults;
            else
            end 
        end
    end
    save(foveaSave,'foveaRoiData','-v7.3');
else
    load(foveaSave,'foveaRoiData')
end


%% RING TEST ANALYSIS
% testSubj = 1;
% ringFiles = '*_ringAll_al_mask.nii.gz';
% ringROI = subfiles([subjFolders{testSubj},'/ROIs/',ringFiles],1);
% roiFiles = '*wangatlas_al_mask.nii.gz';
% allROI = subfiles([subjFolders{testSubj},'/ROIs/',roiFiles],1);
% 
% tmp = niftiRead(allROI{:});
% roiData = tmp.data;
% roiDims = size(roiData);
% tmp = niftiRead(ringROI{:});
% ringData = tmp.data;
% ringDims = size(ringData);
% 
% roiLabel = [cellfun(@(x) [x,'-L'],wangROIs,'uni',false),cellfun(@(x) [x,'-R'],wangROIs,'uni',false)];
% ringList = {'0-1.5','1.5-2.5','2.5-3.5','3.5-4.5','4.5-5.5','5.5-6.5','6.5-7.5','7.5-8.5','8.5-9.5'};
% evcIdx = cell2mat(cellfun(@(x) sum([~isempty(strfind(x,'V1')),~isempty(strfind(x,'V2')),~isempty(strfind(x,'V3'))],2),roiLabel,'uni',false));
% abIdx = cell2mat(cellfun(@(x) sum([isempty(strfind(x,'V3A')),isempty(strfind(x,'V3B'))],2)==2,roiLabel,'uni',false));
% evcIdx = logical(abIdx.*evcIdx);
% tempLabels = roiLabel(evcIdx);
% ringNames=arrayfun(@(x) sprintf('ring%d',x),1:length(ringList),'uni',false); % ring ROIs
% [ii,jj]=ndgrid(1:numel(tempLabels),1:numel(ringNames));
% ringLabel=arrayfun(@(x,y) [ringNames{x},'_',tempLabels{y}],jj(:),ii(:),'uni',false)';
% 
% for t = 1:50
%     cmpIdx = cell2mat(cellfun(@(x) ~isempty(strfind(x,roiLabel{t})),ringLabel,'uni',false));
% 
%     ringMask = ismember(ringData,find(cmpIdx));
%     roiMask = ismember(roiData,t);
% 
%     ringSize(t) = numel(find(ringMask==1));
%     roiSize(t) = numel(find(roiMask==1));
%     ringOverlap(t) = numel(find(ismember(find(roiMask==1),find(ringMask==1))))./numel(find(roiMask==1));
% end
% figure;
% subplot(1,2,1);imshow(roiData(:,:,30)>0)
% subplot(1,2,2);imshow(ringData(:,:,30)>0);



% ringList = cat(1,{CONTring(1,:).name});
% roiList = cat(1,{CONTdata(1,:).name});
% testROI = 4;
% testSubj = 1;
% cmpIdx = cell2mat(cellfun(@(x) ~isempty(strfind(x,roiList{testROI})),ringList,'uni',false));
% cmpData = cat(1,CONTring(testSubj,cmpIdx).rawData);
% cmpROI = mean(mean(cmpData,3));
% wholeROI = mean(mean(CONTdata(testSubj,testROI).rawData,3));

%% COMBINE ROIs and COUNT VOXELS
allConds = [myConds,benoitConds];
for c=1:length(allConds)
     roiData.(allConds{c}) = cat(2,roiData.(allConds{c}),ringRoiData.(allConds{c}),foveaRoiData.(allConds{c}));
end
%get fovea indices
roiList = cat(1,{roiData.cont(1,:).name}); % all ROIs, including fovea
innerIdx = find(cell2mat(arrayfun(@(x) ~isempty(strfind(roiData.cont(1,x).name,'inner')),1:length(roiList),'uni',false)));
outerIdx = find(cell2mat(arrayfun(@(x) ~isempty(strfind(roiData.cont(1,x).name,'outer')),1:length(roiList),'uni',false)));
borderIdx = find(cell2mat(arrayfun(@(x) ~isempty(strfind(roiData.cont(1,x).name,'border')),1:length(roiList),'uni',false)));

for s=1:subCount(1)
    if ~benoitSubj(s)
        curConds = myConds;
    else
        curConds = benoitConds;
    end
    voxCount(s,1:length(roiList)) = NaN;
    for c=1:length(curConds)
        subjROIs = cat(1,{roiData.(curConds{c})(s,:).name});
        runCount(c,s) = size(roiData.(curConds{c})(s,1).rawData,3); % ROI does not matter, just grab first
        runCount(c,s) = runCount(c,s).*~isempty(roiData.(curConds{c})(s,1).rawData);
        if c==1
            voxCount(s,1:length(subjROIs)) = cell2mat(arrayfun(@(x) size(roiData.(curConds{c})(s,x).rawData,2), 1:length(subjROIs),'uni',false));
        else
        end
    end
    for r=1:length(innerIdx)
        foveaCount(s,r,:) = [voxCount(s,innerIdx(r)),voxCount(s,borderIdx(r)),voxCount(s,outerIdx(r))];
    end
end
missingROIs = roiList(nansum(voxCount<5,1)>2); % names of ROIs where more than 2 subjects have fewer than 5 voxels

% combine DIPS and Benoit DISP (C2vsC3):
roiData.alldisp = cat(1,roiData.disp(:,1:size(roiData.C2vsC3,2)),roiData.C2vsC3(benoitSubj,:));
allConds{end+1} = 'alldisp';


%% COMPUTE INCOHERENT SNR
% fitComplex = [];
% for r=1:length(roiSelection)
%     incohIdx = ismember(roiList,roiSelection{r});
%     for s=1:subCount(2)
%         clear runMean;
%         runMean(:,:,1) = mean(roiData.cont(s,incohIdx).rawData,3);
%         runMean(:,:,2) = mean(roiData.mofo(s,incohIdx).rawData,3);
%         if ~isnan(roiData.disp(s,r).harmonics)
%             runMean(:,:,3) = mean(roiData.disp(s,incohIdx).rawData,3);
%         else
%             runMean(:,:,3) = nan(size(runMean(:,:,2)));
%         end
%         nTR = size(runMean,1);
%         maxCycles = round(nTR/2);
%         fftComplex = 2.*fft(runMean,[],1) ./ nTR;
%         harmComplex = fftComplex(2:maxCycles+1,:,:); % don't include DC
%         lst2 = false(size(1:maxCycles));
%         lst2([nCycles-1,nCycles-2,nCycles+1,nCycles+2])=true;
%         voxAmp = squeeze(abs(harmComplex(nCycles,:,:)));
%         voxSNR = squeeze(abs(harmComplex(nCycles,:,:))./nanmean(abs(harmComplex(lst2,:,:)),1));
%         ampIncoh(s,:,r) = nanmean(voxAmp,1);
%         snrIncoh(s,:,r) = nanmean(voxSNR,1);
%         % grab contrast values for fitting
%         if r == 1
%             fitComplex = [fitComplex;squeeze(harmComplex(nCycles,:,1))'];
%         else
%         end
%     end
% end
% % 
% snrIncohMean = squeeze(nanmean(snrIncoh,1));
% snrIncohStderr = squeeze(nanstd(snrIncoh,0,1))./sqrt(subCount(2)');
% ampIncohMean = squeeze(nanmean(ampIncoh,1));
% ampIncohStderr = squeeze(nanstd(ampIncoh,0,1))./sqrt(subCount(2)');
% % figure;
% errorbar(1:length(roiSelection),snrIncohMean(1,:),snrIncohStderr(1,:),'r');
% hold on;
% errorbar(1:length(roiSelection),snrIncohMean(2,:),snrIncohStderr(2,:),'b');
% errorbar(1:length(roiSelection),snrIncohMean(3,:),snrIncohStderr(3,:),'g');
% hold off

%% PROCESS DATA
close all;
for r=1:length(roiList)
    fprintf('processing %s\n',roiList{r});
    for c=1:length(allConds)
        for s=1:subCount(1);
            if (c<4 && benoitSubj(s)) || (r>size(roiData.(allConds{c})(s,:),2))
                % if my conditions and Benoit subject, or if data is empty, simply replace with NaNs 
                addNans = true;
            else
                if isempty(roiData.(allConds{c})(s,r).harmonics) || any(isnan(roiData.(allConds{c})(s,r).harmonics))
                    addNans = true;
                else
                    addNans = false;
                end
            end
            if benoitSubj(s) 
                repData = roiData.C1vsC2(s,1); % replacement data, used for size
            else
                repData = roiData.cont(s,1);
            end
            if addNans
                
                roiData.(allConds{c})(s,r).name = NaN;
                roiData.(allConds{c})(s,r).harmonics = NaN(size(repData.harmonics));
                roiData.(allConds{c})(s,r).meanCycle = NaN(size(repData.meanCycle));
                roiData.(allConds{c})(s,r).zScore = NaN(size(repData.zScore));
                roiData.(allConds{c})(s,r).SNR = NaN(size(repData.SNR));
                roiData.(allConds{c})(s,r).phase = NaN(size(repData.phase));
                roiData.(allConds{c})(s,r).rawData = NaN(size(repData.rawData,1),1,size(repData.rawData,3)); % since this is fake data, make number of voxels 1.
                roiData.(allConds{c})(s,r).realSignal = NaN(size(repData.realSignal));
                roiData.(allConds{c})(s,r).imagSignal = NaN(size(repData.imagSignal));
                roiData.(allConds{c})(s,r).realNoise = NaN(size(repData.realNoise));
                roiData.(allConds{c})(s,r).imagNoise = NaN(size(repData.imagNoise));
            else
            end
            % incoherent SNR field does not exist, so always add NaNs
            roiData.(allConds{c})(s,r).incohSNR = NaN(size(repData.SNR));
            clear repData;
        end
        
        % find those subjects whose sum(meanCycle) is not NAN
        nanIdx{r,c} = cell2mat(arrayfun(@(x)~isnan(sum(x.meanCycle)),roiData.(allConds{c})(:,r),'uni',false )); 
        numSubs(r,c) = length(find(nanIdx{r,c} == 1));
        % concatenate each zScore by row, create a subjectXzScore matrix. (zScore is a 1x5 vector)
        allZ{r}(:,:,c) = cat(2,roiData.(allConds{c})(:,r).zScore); 
        % calculate the column mean of zScore matrix
        allMeanZ(c,r,:)=nanmean(cat(2,roiData.(allConds{c})(:,r).zScore));
        % and standard deviation
        allStdevZ(c,r,:)=nanstd(cat(2,roiData.(allConds{c})(:,r).zScore))/sqrt(numSubs(r,c));
    
        % COMPLEX
        allRealSignal(:,:,c,r) = cat(2,roiData.(allConds{c})(:,r).realSignal)';
        allImagSignal(:,:,c,r) = cat(2,roiData.(allConds{c})(:,r).imagSignal)';  
        allRealNoise(:,:,c,r) = cat(2,roiData.(allConds{c})(:,r).realNoise)';
        allImagNoise(:,:,c,r) = cat(2,roiData.(allConds{c})(:,r).imagNoise)';
 
        % ADJUST PHASE phaseS
        for n=1:2
            if n==1
                tempComplex = complex(allRealSignal(:,:,c,r),allImagSignal(:,:,c,r));
            else
                tempComplex = complex(allRealNoise(:,:,c,r),allImagNoise(:,:,c,r));
            end
            tempMag = abs(tempComplex);
            tempAngle = angle(tempComplex);
            if c>2 % if disp or Benoit data, rotate by pi (180 deg.)
                tempAngle = tempAngle+pi;     % condition 3 is 180 degrees out of phase
            else
            end
            if c<4 % phase angle adjusment due to mux, only on my data
                tempAngle = tempAngle-pi/3;                 % pi/3 radians = 60 degrees = 4 seconds (pre-TR 16, rather than 12, seconds)
            else
            end
            if c>6 % if combined disp data
                tempAngle(~benoitSubj,:) = tempAngle(~benoitSubj,:)-pi/3;
            else
            end
            tempAngle = phaseMod(tempAngle); % adjust to sine base
            tempComplex = tempMag.*exp(1i*tempAngle); % convert back to complex numbers
            if n==1
                allRealSignal(:,:,c,r) = real(tempComplex);
                allImagSignal(:,:,c,r) = imag(tempComplex);
            else
                allRealNoise(:,:,c,r) = real(tempComplex);
                allImagNoise(:,:,c,r) = imag(tempComplex);
            end
            clear temp*
        end
        
        % COMPUTE VECTOR MEAN AND ERRORS
        if numSubs(r,c) > 1 % if at least two subjects have data
            tempVecErr = fitErrorEllipse([allRealSignal(nanIdx{r,c},whichHarm,c,r),allImagSignal(nanIdx{r,c},whichHarm,c,r)],'SEM',false);
            vecErrLow(c,r) = tempVecErr(1);
            vecErrHigh(c,r) = tempVecErr(2);
            meanComplex(c,r) = complex(mean(allRealSignal(nanIdx{r,c},whichHarm,c,r),1),mean(allImagSignal(nanIdx{r,c},whichHarm,c,r),1));
            vecMeanAmp(c,r) = abs(meanComplex(c,r));
            vecMeanPhase(c,r) = phase(meanComplex(c,r));
            if strcmp(roiList{r},'LO1') && c==7
                disp('LO1');
            else
            end
            tempResults = tSquaredFourierCoefs([allRealSignal(nanIdx{r,c},whichHarm,c,r),allImagSignal(nanIdx{r,c},whichHarm,c,r)],'testMu',[0,0],'alphaVal',0.05);
            vecMeanP(c,r) = tempResults.pVal;
            if vecMeanPhase(c,r)< 0 
                vecMeanPhase(c,r) = 2*pi+vecMeanPhase(c,r);
            else
            end
        else
            vecErrLow(c,r) = NaN;
            vecErrHigh(c,r) = NaN;
            meanComplex(c,r) = NaN;
            vecMeanAmp(c,r) = NaN;
            vecMeanPhase(c,r) = NaN;
            vecMeanP(c,r) = NaN;
        end
        
        % COMPUTE MEAN CYCLE
        % get average time courses across ROI for each subject
        % average across the runs, then average across voxels, results is a 120X#ofsubj matrix (mean for each time point) 
        allTime(:,:,c,r) = cell2mat(arrayfun(@(x) double(nanmean(nanmean(roiData.(allConds{c})(x,r).rawData,3),2)),1:subCount(1),'uni',false)); 
        aveTime(:,c,r) = nanmean(allTime(:,:,c,r),2);  % for each ROI (indexed by r), average across subject  

        % recompute average cycle
        for s=1:subCount(1)
            tempTime = allTime(:,s,c,r);
            if c<3 % if mofo or cont conditions 
                tempTime = [nan(2,1);tempTime;nan(10,1)]; % add 2 nans in the beginning, 10 in the end
                tempMean = nanmean(reshape(tempTime,size(tempTime,1)/(nCycles+1),nCycles+1),2); % take "fake cycle" into account
            elseif c== 3 % if disp condition
                tempTime = [nan(8,1);tempTime;nan(4,1)]; % add 8 nans in the beginning, 4 in the end
                tempMean = nanmean(reshape(tempTime,size(tempTime,1)/(nCycles+1),nCycles+1),2); % take "fake cycle" into account
            elseif ( c == 7 && ~benoitSubj(s) ) % if all disp, and my subjects
                tempTime = [nan(8,1);tempTime;nan(4,1)]; % add 8 nans in the beginning, 4 in the end
                tempMean = nanmean(reshape(tempTime,size(tempTime,1)/(nCycles+1),nCycles+1),2); % take "fake cycle" into account
            else
                % benoit data have no mux, but are off by half a cycle
                tempTime = [nan(6,1);tempTime;nan(6,1)]; % add 12 nans in the beginning, 0 in the end
                tempMean = nanmean(reshape(tempTime,size(tempTime,1)/(nCycles+1),nCycles+1),2); % take "fake cycle" into account
            end
            roiData.(allConds{c})(s,r).meanCycle = tempMean - (max(tempMean)+min(tempMean))./2;
        end

        allMeanCycle(c,:,r) = nanmean(cat(2,roiData.(allConds{c})(:,r).meanCycle),2);
        allStdevCycle(c,:,r) = nanstd(cat(2,roiData.(allConds{c})(:,r).meanCycle),0,2)/sqrt(numSubs(r,c));

        % COMPUTE HARMONICS
        nScans = size(roiData.(allConds{c})(1,1).rawData,1); % same for all ROIs and subjects
        maxCycles = round(nScans/2);
        xHarm= 1:maxCycles;
        allMeanHarm(c,:,r) = nanmean(cat(2,roiData.(allConds{c})(:,r).harmonics));
        allStdevHarm(c,:,r) = nanstd(cat(2,roiData.(allConds{c})(:,r).harmonics))/sqrt(numSubs(r,c));
        
        % FIX SNR
        if c==length(allConds) && whichSNR == 2 % if last condition and we are using across-condition noise floor
        % generate new SNR value and replace old one
            nHarm = length(roiData.cont(1,1).SNR); % number of harmonics in data
            for s=1:subCount(1)
                if ~benoitSubj(s)
                    curConds = [myConds,'alldisp'];
                else
                    curConds = [benoitConds,'alldisp'];
                end
                 
                for h=1:nHarm
                    %make the two bands around the target cycle non-zero in order to calculate noise
                    lst = false(size(roiData.cont(1,1).harmonics));
                    lst([nCycles*h-1,nCycles*h-2,nCycles*h+1,nCycles*h+2])=true;                     
                    % get sidebands from all conditions
                    tempNoise = cell2mat(arrayfun(@(x) double(roiData.(curConds{x})(s,r).harmonics(lst)'),1:length(curConds),'uni',false));
                    tempNoise = nanmean(nanmean(tempNoise)); % average over sidebands and conditions
                    
                    % compute noise values for incoherent SNR
                    runMean = arrayfun(@(x) nanmean(roiData.(curConds{x})(s,r).rawData,3),1:length(curConds),'uni',false);
                    fftComplex = cellfun(@(x) 2.*fft(x,[],1) ./ nTR,runMean,'uni',false);
                    fftComplex = cellfun(@(x) x(2:maxCycles+1,:),fftComplex,'uni',false);
                    tempNanIdx = cell2mat(cellfun(@(x) ~isnan(x(1,1)), fftComplex,'uni',false));
                    % average the amplitude across side bands and conditions
                    voxNoise = nanmean(cell2mat(cellfun(@(x) nanmean(abs(x(lst,:)),1),fftComplex(tempNanIdx),'uni',false)'),1);
                    
                    for cCur=1:length(curConds)
                        if ~isnan(roiData.(curConds{cCur})(s,r).SNR(1)) % if subject has data
                            roiData.(curConds{cCur})(s,r).SNR(h,1) = roiData.(curConds{cCur})(s,r).harmonics(nCycles*h) / tempNoise;
                            voxSNR = abs(fftComplex{cCur}(nCycles*h,:))./voxNoise;
                            roiData.(curConds{cCur})(s,r).incohSNR(h,1) = nanmean(voxSNR,2);
                        else
                        end
                    end
                end
            end
        else
        end
    end
    % combine SNR values
    for c = 1:length(allConds)
        allSNR(:,:,c,r) = cat(2,roiData.(allConds{c})(:,r).SNR)';
        allIncohSNR(:,:,c,r) = cat(2,roiData.(allConds{c})(:,r).incohSNR)';
    end
end
       

%% SET UP COLORS AND BASIC PLOTTING VARIABLES
% condition colors
condColors = [1 0 0; 0 0 1; 0 1 0; 1 0.83 0 ; 0 1 0; 1 0.1034 0.72]; 


% roi colors
colormap('jet');
tempColors = colormap;
tempIdx = round(linspace(0,length(tempColors),length(roiSelection)+1));
if isempty(which('flip'))  
    roiColors = tempColors([tempIdx(3:2:end),fliplr(tempIdx(2:2:end))],:);
else
    roiColors = tempColors([tempIdx(3:2:end),flip(tempIdx(2:2:end))],:);
end    
fontSize = 12;
gcaOpts = {'tickdir','out','box','off','fontsize',fontSize,'fontname','Arial','linewidth',1};

condNames = myConds;
condNames{1} = 'contrast';
condNames{2} = 'motion';
condNames{3} = 'disparity';
condNames{4} = 'in-phase/uncorr'; % 'C1vsC2' - note that A and B block are swapped
condNames{5} = 'anti-phase/in-phase'; % 'C2vsC3' - note that A and B block are swapped
condNames{6} = 'anti-phase/uncorr'; % 'C2vsC3' - note that A and B block are swapped
condNames{7} = 'alldisp';

%% PLOT ROI FIGURES
close all;
if plotROIs
    roiFigWidth=[16,6];
    for p = 1:2
        condIdx = (1:length(myConds))+3.*(p-1);
        for r=1:length(roiSelection)
            roiIdx = ismember(roiList,roiSelection{r});
            if length(nanIdx{roiIdx,condIdx(1)})>8 % if enough subjects have this ROI
                figure;
                % PLOT MEAN CYCLE
                subplot(2,3,[1,4])
                hold on;
                % label figure with ROI name
                text(1,0.6,roiSelection{r},'fontname','Arial','fontsize',fontSize,'horizontalalignment','left');
                % make condition boxes
                paH(1) = patch([0.1,0.1,6,6],[-.5,-.6,-.6,-.5],[153 153 153]./255);
                paH(2) = patch([6.0,6.0,12.1,12.1],[-.5,-.6,-.6,-.5],[51 51 51]./255);
                set(paH(:),'EdgeColor','none')
                text(3,-.55,'A Block','horizontalAlignment','center','FontSize',fontSize,'fontname','Arial')
                text(9,-.55,'B Block','horizontalAlignment','center','FontSize',fontSize,'fontname','Arial','Color', [1 1 1]);    
                for c=1:length(condIdx)
                    pH(condIdx(c))=plot(.5:11.5,allMeanCycle(condIdx(c),:,roiIdx),'LineWidth',2,'color',condColors(condIdx(c),:));
                    errorb(.5:11.5,allMeanCycle(condIdx(c),:,roiIdx),allStdevCycle(condIdx(c),:,roiIdx),'color',condColors(condIdx(c),:))
                end
                cycleLim = [-.6,.6];
                ylabel(gca,'% signal change','FontSize',fontSize)
                xlabel(gca,'secs','FontSize',fontSize);
                text(28,min(cycleLim)+diff(cycleLim)*.35,['n=',num2str(cellfun(@(x) length(x),nanIdx(roiIdx,condIdx)),' %0.0f    |')],'FontSize',fontSize,'fontname','Arial','horizontalalignment','right');
                text(28,min(cycleLim)+diff(cycleLim)*.25,['vecAmp=',num2str(vecMeanAmp(condIdx,roiIdx)',' %0.2f |')],'FontSize',fontSize,'fontname','Arial','horizontalalignment','right');
                text(28,min(cycleLim)+diff(cycleLim)*.15,['delay(s)=',num2str(vecMeanPhase(condIdx,roiIdx)'./(2*pi).*24,' %0.1f |')],'FontSize',fontSize,'fontname','Arial','horizontalalignment','right');       
                text(28,min(cycleLim),sprintf('mean vox=%0.1f',mean(voxCount(condIdx,roiIdx))),'FontSize',fontSize,'fontname','Arial','horizontalalignment','right');
                set(gca,'xtick',0:2:12,'xticklabel',0:4:24,'ytick',min(cycleLim):.2:max(cycleLim),...
                        'ylim',[min(cycleLim),max(cycleLim)],'xlim',[0,12],gcaOpts{:},'ticklength',[.03,.03]); 
                hold off;    
                % PLOT HARMONICS
                hIdx = [2,3,6];
                for c=1:length(condIdx)
                    subplot(2,3,hIdx(c))
                    plot(xHarm,allMeanHarm(condIdx(c),1:maxCycles,roiIdx),'-','color',[.5 .5 .5],'LineWidth',2);
                    hold on
                    if nCycles>1
                        pH = plot(xHarm(nCycles-1:nCycles+1),allMeanHarm(condIdx(c),nCycles-1:nCycles+1,roiIdx),'color',condColors(condIdx(c),:),'LineWidth',2);
                        lH = legend(pH,condNames{condIdx(c)},'location','northeast');
                        legend(gca,'boxoff');
                        lPos = get(lH,'position');
                        lPos(2) = lPos(2)+0.005;
                        set(lH,'position',lPos);

                    else
                        plot(xHarm,allMeanHarm(condIdx(c),:,r),'k','LineWidth',2)
                    end
                    ymax = .6;
                    ylim([0,ymax]);
                    xlim([.5 60.5]);
                    hold off
                    % Ticks
                    xtick=nCycles:nCycles:(maxCycles+1);
                    set(gca,'xtick',xtick,'ytick',0:.2:ymax,gcaOpts{:},'ticklength',[.04,.04])
                end
                % PRINT FIGURE
                set(gcf, 'units', 'centimeters'); % make figure size units centimeters
                oldPos = get(gcf, 'pos');
                newPos = oldPos;
                newPos(3) = roiFigWidth(1);
                newPos(4) = roiFigWidth(2);
                set(gcf, 'pos',newPos);
                drawnow;
                tightfig;
                if p==1
                    figName = [figFolder,'/',roiSelection{r},'.',suffix];
                else
                    figName = [figFolder,'/benoit_',roiSelection{r},'.',suffix];
                end
                export_fig([figName,'.pdf'],'-pdf','-transparent',gcf),
            else
            end
        end
    end
else
end
close all;
%% COMBINE DATA ACROSS D/V AND L/R
% rIdx = 1:length(roiSelection); % swap out MT and V1 to use MT as the baseline
% 
% for r=1:length(roiSelection)
%     curIdx{1} = find(ismember(roiList,roiSelection{r}));
%     curIdx{2} = find(ismember(roiList,[roiSelection{r},'-L']));
%     curIdx{3} = find(ismember(roiList,[roiSelection{r},'-R']));
%     for z=1:3
%         dataSNR{z}(:,:,:,r)=allSNR(:,:,:,curIdx{z});
%         dataIncohSNR{z}(:,:,:,r)=allIncohSNR(:,:,:,curIdx{z});
%         partVecAmp{z}(r,:) = vecMeanAmp(:,curIdx{z});
%         partVecReal{z}(r,:) = real(meanComplex(:,curIdx{z}));
%         partVecImag{z}(r,:) = imag(meanComplex(:,curIdx{z}));
%         partVecPhase{z}(r,:) = vecMeanPhase(:,curIdx{z});
%         partVecLow{z}(r,:) = vecErrLow(:,curIdx{z});
%         partVecHigh{z}(r,:) = vecErrHigh(:,curIdx{z});
%         partVecP{z}(r,:) = vecMeanP(:,curIdx{z});
%         numSubs(r,:,z) = arrayfun(@(x) length(nanIdx{curIdx{z},x}),1:size(dataSNR{1},3));
%         partIncohSNR
%     end
% end
% partSNRmean = cellfun(@(x) permute(squeeze(nanmean(x,1)),[3,2,1]),dataSNR,'uni',0);
% partSNRmean = cellfun(@(x) permute(squeeze(nanmean(x,1)),[3,2,1]),dataIncohSNR,'uni',0);
% 
% 
% partSNRstd = cellfun(@(x) permute(squeeze(nanstd(x,0,1)),[3,2,1]),dataSNR,'uni',0);
% for z=1:3
%     partSNRerr{z} = partSNRstd{z}./sqrt(repmat(numSubs(:,:,1),1,1,5));
% end

%% PHASE PLOT
% figure;
% roiIdx = cell2mat(cellfun(@(x) ismember(roiList,x),roiSelection,'uni',false));
% phaseColors = [eye(3);[0,1,1]];
% phaseTitle = roiSelection(roiIdx);
% 
% plot(1:length(roiSelection),vecMeanPhase(:,roiIdx)./(pi*2)*24);

%% SNR PLOT
% close all
% vecPlot = true;
% lrFig = figure; bothFig = figure; EVCfig = figure;
% orderIdx(1,:) = cell2mat(cellfun(@(x) sum(strcmp(x,{'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' 'LO1' 'LO2'})),roiSelection,'uni',false))*2;
% orderIdx(2,:) = cell2mat(cellfun(@(x) sum(strcmp(x,{'TO1' 'TO2' 'V3A' 'V3B'})),roiSelection,'uni',false))*3;
% orderIdx(3,:) = cell2mat(cellfun(@(x) sum(strcmp(x,{'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4','IPS5' 'SPL1' 'FEF'})),roiSelection,'uni',false))*4;
% orderIdx(4,:) = ~sum(orderIdx,1);
% orderIdx = sum(orderIdx,1);
% if plotSNR
%     titleStr = {'Early Visual Areas','Dorsal Stream','Ventral Stream','Parietal'};
%     xAxis = 1:3;
%     evcLabels = {'LH dorsal','RH dorsal','LH ventral','RH ventral'};
%     
%     if whichHarm == 1
%         legendPos = 'NorthEast';
%         labelPos = min(xAxis);
%     else
%         legendPos = 'NorthEast';
%         labelPos = max(xAxis);
%     end
%     
%     for plotType = 1:7
%         if plotType<4
%             numROIs = length(roiSelection);
%         else
%             numROIs = 3;
%         end
%         for r=1:numROIs
%             plotCount = orderIdx(r);
%             switch plotType
%                 case 1 % plot grand mean
%                     figure(bothFig); subplot(1,4,plotCount); hold on;
%                     plotLabel = ' ';
%                     plotSize = [1.5,.75];
%                     figName = [figFolder,'/bothSNR_harm',num2str(whichHarm),'.',suffix];
%                 case 2 % plot left and right
%                     figure(lrFig); 
%                     subplot(2,4,plotCount); hold on;
%                     plotLabel = 'LH';
%                     plotSize = [1.5,1.5];
%                     figName = [figFolder,'/hemiSNR_harm',num2str(whichHarm),'.',suffix];
%                 case 3
%                     figure(lrFig); 
%                     subplot(2,4,plotCount+4); hold on;
%                     plotLabel = 'RH';
%                     plotSize = [1.5,1.5];
%                     figName = [figFolder,'/hemiSNR_harm',num2str(whichHarm),'.',suffix];
%                 otherwise
%                     figure(EVCfig);
%                     subplot(2,2,plotType-3); hold on;
%                     plotLabel = evcLabels(plotType-3);
%                     plotSize = [1.5,1.5];
%                     figName = [figFolder,'/evcSNR_harm',num2str(whichHarm),'.',suffix];
%             end
%             if vecPlot
%                 snrPlot(r) = plot(xAxis,partVecAmp{plotType}(r,:),'o-','LineWidth',2,'color',roiColors(r,:));
%                 set(snrPlot(r), 'MarkerFaceColor', get(snrPlot(r), 'Color'));
%                 lErr = abs(partVecAmp{plotType}(r,:)-partVecLow{plotType}(r,:));
%                 uErr = abs(partVecAmp{plotType}(r,:)-partVecHigh{plotType}(r,:));
%                 errorbar(xAxis,partVecAmp{plotType}(r,:),lErr,uErr,'LineWidth',2,'color',roiColors(r,:));
%                 yMax = 0.6; yUnit = 0.2; yMin = -0.05;
%             else
%                 snrPlot(r)=plot(xAxis,partSNRmean{plotType}(r,:,whichHarm),'o-','LineWidth',2,'color',roiColors(r,:));
%                 set(snrPlot(r), 'MarkerFaceColor', get(snrPlot(r), 'Color'));
%                 errorb(xAxis,partSNRmean{plotType}(r,:,whichHarm),partSNRstd{plotType}(r,:,whichHarm)./sqrt(numSubs(r,1)),'LineWidth',2,'color',roiColors(r,:));
%                 yMax = 9; yUnit = 1;yMin = 0;
%             end
%             text(labelPos,9.5,plotLabel,'HorizontalAlignment','center','VerticalAlignment','top','fontname','Arial');
%             text(median(xAxis),0.5, titleStr{plotCount},'HorizontalAlignment','center','VerticalAlignment','top','fontname','Arial');
%             set(gca,'xtick',1:3,'xticklabel',conditions,'ytick',-1:yUnit:yMax,gcaOpts{:});
%             xlim([.5 4.5]);
%             ylim([yMin,yMax]);
%             if whichSNR==1
%                 ylabel('SNR','fontname','Arial'); % SNR: 4 sidebands, within-condition noise'
%             else
%                 ylabel('SNR','fontname','Arial'); % 'SNR: 4 sidebands, across-condition noise'
%             end
%             xlabel('wallpaper group','fontname','Arial');
%             if r==length(partSNRmean{plotType}(:,1,whichHarm))
%                 if plotType <4
%                     for z=1:4
%                         legend(snrPlot(orderIdx==z),roiSelection(orderIdx==z),'location',legendPos,'fontsize',12)
%                     end
%                 else
%                     legend(snrPlot(1:3),roiSelection(1:3),'location',legendPos,'fontsize',12)
%                 end
%             else
%             end
%             hold off
%         end
%         clear snrPlot;
%         if sum(plotType == [1,3,7])
%             pos = get(gcf, 'Position');
%             pos(3) = pos(3)*plotSize(1); % Select the height of the figure in [cm]
%             pos(4) = pos(4)*plotSize(2);
%             set(gcf, 'Position', pos);
%             export_fig([figName,'.pdf'],'-pdf','-transparent',gcf);
%         else
%         end
%     end
% else
% end

%% CONVERT RING DATA
clear ring*; 
pCutOff = 0.05;
evcROIs = {'V1','V2','V3'};
ringList = unique(cellfun(@(x) x(7:end),roiList,'uni',false));
ringList = ringList(cellfun(@(x) ~isempty(strfind(x,'-')), ringList));
for r = 1:length(ringList);
    ringSNR(:,:,:,r)=cell2mat(cellfun(@(x) allSNR(:,whichHarm,:,ismember(roiList,[x,'ring',ringList{r}])),evcROIs,'uni',false));
    % ringSNR dimensions are subj x ROI x condition x ring
    ringReal(:,:,:,r)=cell2mat(cellfun(@(x) allRealSignal(:,whichHarm,:,ismember(roiList,[x,'ring',ringList{r}])),evcROIs,'uni',false));
    ringImag(:,:,:,r)=cell2mat(cellfun(@(x) allImagSignal(:,whichHarm,:,ismember(roiList,[x,'ring',ringList{r}])),evcROIs,'uni',false));
    
    ringVecMean(r,:,:)=cell2mat(cellfun(@(x) vecMeanAmp(:,ismember(roiList,[x,'ring',ringList{r}])),evcROIs,'uni',false));
    ringVecPhase(r,:,:) = cell2mat(cellfun(@(x) vecMeanPhase(:,ismember(roiList,[x,'ring',ringList{r}])),evcROIs,'uni',false));
    ringVecHigh(r,:,:) = cell2mat(cellfun(@(x) vecErrHigh(:,ismember(roiList,[x,'ring',ringList{r}])),evcROIs,'uni',false));
    ringVecLow(r,:,:) = cell2mat(cellfun(@(x) vecErrLow(:,ismember(roiList,[x,'ring',ringList{r}])),evcROIs,'uni',false));
    ringVecP(r,:,:) = cell2mat(cellfun(@(x) vecMeanP(:,ismember(roiList,[x,'ring',ringList{r}])),evcROIs,'uni',false));
    ringComplex(r,:,:) = cell2mat(cellfun(@(x) meanComplex(:,ismember(roiList,[x,'ring',ringList{r}])),evcROIs,'uni',false));
    
    tempCycleMean = cellfun(@(x) allMeanCycle(:,:,ismember(roiList,[x,'ring',ringList{r}])),evcROIs,'uni',false);
    ringCycleMean(:,:,:,r) = reshape(cell2mat(tempCycleMean),length(condNames),12,length(evcROIs));
    tempCycleErr = cellfun(@(x) allStdevCycle(:,:,ismember(roiList,[x,'ring',ringList{r}])),evcROIs,'uni',false);
    ringCycleErr(:,:,:,r) = reshape(cell2mat(tempCycleErr),length(condNames),12,length(evcROIs));
    for z=1:length(evcROIs)
        ringNum(r,:,z) = sum(~isnan(ringSNR(:,z,:,r)));
    end
end
ringSNRMean = permute(squeeze(nanmean(ringSNR,1)),[3,2,1]); % ring x condition x ROI
ringSNRerr = permute(squeeze(nanstd(ringSNR,0,1)),[3,2,1])./sqrt(ringNum); % ring x condition x ROI

%% CONVERT FOVEA DATA
clear fovea*; 
evcROIs = {'V1','V2','V3','hV4'};
foveaList = {'inner','border','outer'};
for f = 1:length(foveaList);
    for r=1:length(evcROIs)
        if isempty(allSNR(:,whichHarm,:,ismember(roiList,[evcROIs{r},'_',foveaList{f}])))
            foveaSNR(:,r,:,f)=NaN(subCount,3);
            foveaVecMean(f,:,r) = NaN;
        else
            foveaSNR(:,r,:,f) = allSNR(:,whichHarm,:,ismember(roiList,[evcROIs{r},'_',foveaList{f}]));
            foveaReal(:,r,:,f)=allRealSignal(:,whichHarm,:,ismember(roiList,[evcROIs{r},'_',foveaList{f}]));
            foveaImag(:,r,:,f)=allImagSignal(:,whichHarm,:,ismember(roiList,[evcROIs{r},'_',foveaList{f}]));
            foveaVecMean(f,:,r)=vecMeanAmp(:,ismember(roiList,[evcROIs{r},'_',foveaList{f}]));
            foveaVecPhase(f,:,r)=vecMeanPhase(:,ismember(roiList,[evcROIs{r},'_',foveaList{f}]));
            foveaVecHigh(f,:,r)=vecErrHigh(:,ismember(roiList,[evcROIs{r},'_',foveaList{f}]));
            foveaVecLow(f,:,r)=vecErrLow(:,ismember(roiList,[evcROIs{r},'_',foveaList{f}]));
            foveaVecP(f,:,r)=vecMeanP(:,ismember(roiList,[evcROIs{r},'_',foveaList{f}]));
            foveaComplex(f,:,r)=meanComplex(:,ismember(roiList,[evcROIs{r},'_',foveaList{f}]));

            foveaCycleMean(:,:,r,f) = allMeanCycle(:,:,ismember(roiList,[evcROIs{r},'_',foveaList{f}]));
            foveaCycleErr(:,:,r,f) = allStdevCycle(:,:,ismember(roiList,[evcROIs{r},'_',foveaList{f}]));
            foveaNum(r,:,f) = sum(~isnan(foveaSNR(:,r,:,f)));
        end
    end
end
foveaSNRMean = permute(squeeze(nanmean(foveaSNR,1)),[3,2,1]); % fovea x condition x ROI
foveaSNRerr = permute(squeeze(nanstd(foveaSNR,0,1))./sqrt(foveaNum),[3,2,1]); % fovea x condition x ROI

%% FIND PHASE SIGN
fitSeperate = false;
for c=1:length(condNames)
    % fit on all data across ROIs
    if fitSeperate
        fitVals = ringComplex(:,c,:);
    else
        % use values from V1, across all subjects, for contrast only
        %fitVals = fitComplex;
        %fitIdx = ismember(roiList,'V1');
        %fitVals = complex(allRealSignal(:,whichHarm,1,fitIdx),allImagSignal(:,whichHarm,1,fitIdx));
        %fitVals = fitVals(~isnan(fitVals));
        fitVals = ringComplex(:,1,:); % contrast only
    end
    fitVals = fitVals(:);
    lineFit(:,c) = polyfitZero(real(fitVals(:)),imag(fitVals(:)),1);
    linePerp(c) = 1/(-lineFit(1,c));
    for z=1:3 % loop over bilateral, left and right hemisphere
        % do whole ROIs
        switch z
            case 1
                % bilateral
                curIdx = find(ismember(roiList,roiSelection));
            case 2
                curIdx = find(ismember(roiList,cellfun(@(x) [x,'-L'],roiSelection,'uni',false)));
            case 3
                curIdx = find(ismember(roiList,cellfun(@(x) [x,'-R'],roiSelection,'uni',false)));
        end
        
        % compute aggregate values for sign, amplitude, error
        roiAggregate(z,c).numSubs = numSubs(curIdx,c);
        roiAggregate(z,c).Vec.sign = sign( imag(meanComplex(c,curIdx)) - linePerp(c)*real(meanComplex(c,curIdx)) );
        roiAggregate(z,c).Vec.ampMean  = vecMeanAmp(c,curIdx);
        roiAggregate(z,c).Vec.errLow  = vecErrLow(c,curIdx);
        roiAggregate(z,c).Vec.errHigh  = vecErrHigh(c,curIdx);
        roiAggregate(z,c).Vec.Pval = vecMeanP(c,curIdx);
        roiAggregate(z,c).Vec.Phase  = vecMeanPhase(c,curIdx);
        
        roiAggregate(z,c).SNRcoh.Mean = nanmean(squeeze(allSNR(:,whichHarm,c,curIdx)));
        roiAggregate(z,c).SNRcoh.stdErr = nanstd(squeeze(allSNR(:,whichHarm,c,curIdx)),0,1)./numSubs(curIdx,c)';
        roiAggregate(z,c).SNRincoh.Mean = nanmean(squeeze(allIncohSNR(:,whichHarm,c,curIdx)));
        roiAggregate(z,c).SNRincoh.stdErr = nanstd(squeeze(allIncohSNR(:,whichHarm,c,curIdx)),0,1)./numSubs(curIdx,c)';
        
        [~,roiAggregate(z,c).SNRcoh.Pval,~,tempStruct] = ttest(squeeze(allSNR(:,whichHarm,c,curIdx)),1,'dim',1,'tail','right');
        roiAggregate(z,c).SNRcoh.Tval = tempStruct.tstat;
        [~,roiAggregate(z,c).SNRincoh.Pval,~,tempStruct] = ttest(squeeze(allIncohSNR(:,whichHarm,c,curIdx)),1,'dim',1,'tail','right');
        roiAggregate(z,c).SNRincoh.Tval = tempStruct.tstat;
    end
    for r=1:length(evcROIs)
        % note: use ring fits for fovea sign
        foveaSign(:,c,r) = sign(imag(foveaComplex(:,c,r)) - linePerp(c)*real(foveaComplex(:,c,r)));
        if r<4
            % compute signs based on fits
            ringSign(:,c,r) = sign(imag(ringComplex(:,c,r)) - linePerp(c)*real(ringComplex(:,c,r)));
        else
        end
    end
end

%% SIMPLER PLOT
close all;
orderIdx(1,:) = cell2mat(cellfun(@(x) sum(strcmp(x,{'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2'})),roiSelection,'uni',false))*2;
orderIdx(2,:) = cell2mat(cellfun(@(x) sum(strcmp(x,{'LO1' 'LO2' 'TO1' 'TO2' 'V3A' 'V3B'})),roiSelection,'uni',false))*3;
orderIdx(3,:) = cell2mat(cellfun(@(x) sum(strcmp(x,{'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4','IPS5' 'SPL1' 'FEF'})),roiSelection,'uni',false))*4;
orderIdx(4,:) = ~sum(orderIdx,1);
orderIdx = sum(orderIdx,1);
missingIdx = roiAggregate(1,4).numSubs<8; % missing ROIs in Benoit subs
expIdx{1} = orderIdx;
expIdx{2} = orderIdx;
expIdx{2}(missingIdx>0)=0; 
plotVec = true;
plotSign = true;
simpleSelection = roiSelection;
pCutOff = 0.05;
offSet = [0 .15 -.15];
splitPoints{1} = [3.75,9.25,14.75];
splitPoints{2} = [3.75,6.25,11.75];
legOffset = [.05,.075,-.1,0];

for p = 1:2
    simpFig(p) = figure;
    condIdx = (1:length(myConds))+3.*(p-1);
    xAxis = 0; xTicks =[]; xOrder=[]; yMax = 0.6; yMin=-0.1; yUnit = 0.1;
    sigPos = [-.03 -.055 -.08];
    for z=1:max(expIdx{p})
        if sum(expIdx{p}==z)>0
            figure(simpFig(p));
            hold on
            if xAxis == 0
                xAxis = 1:length(find(expIdx{p}==z));
            else
                xAxis = (max(xAxis)+.5)+(1:length(find(expIdx{p}==z)));
            end
            if z==min(expIdx{p}(expIdx{p}>0))
                plot(linspace(0,200,10),zeros(10,1),'linewidth',2,'color',[.5 .5 .5]);
                sigH = plot(-1000,-1000,'ksq','linewidth',2);
                set(sigH, 'MarkerFaceColor', get(sigH, 'Color'));
                plot(ones(10,1)*splitPoints{p},repmat(linspace(-10,10,10),length(splitPoints{p}),1)','linewidth',2,'color',[.5 .5 .5]);
            else
            end
            for c = 1:length(condIdx)
                if plotVec
                    if plotSign
                        
                        plotVals = roiAggregate(1,condIdx(c)).Vec.ampMean(expIdx{p}==z).* roiAggregate(1,condIdx(c)).Vec.sign(expIdx{p}==z);
                        sigPos = [.4 .36 .32];
                        yMax = 0.4; yMin=-0.6; yUnit = 0.1;
                        legLoc = 'southeast';
                    else
                        plotVals = roiAggregate(1,condIdx(c)).Vec.ampMean(expIdx{p}==z);
                        legLoc = 'northeast';

                    end
                    lErr = abs(roiAggregate(1,condIdx(c)).Vec.ampMean(expIdx{p}==z)-roiAggregate(1,condIdx(c)).Vec.errLow(expIdx{p}==z));
                    uErr = abs(roiAggregate(1,condIdx(c)).Vec.ampMean(expIdx{p}==z)-roiAggregate(1,condIdx(c)).Vec.errHigh(expIdx{p}==z));
                    simpPlot(z,c) = plot(xAxis+offSet(c),plotVals,'o-','color',condColors(condIdx(c),:),'LineWidth',2);
                    set(simpPlot(z,c), 'MarkerFaceColor', get(simpPlot(z,c), 'Color'),'markerSize',4);
                    e_H(c,z) = errorbar(xAxis+offSet(c),plotVals,lErr,uErr,'LineWidth',2,'color',condColors(condIdx(c),:));
                    curSig = double(roiAggregate(1,condIdx(c)).Vec.Pval(expIdx{p}==z)<pCutOff);
                    curSig(curSig==0) = NaN; % use NaNs to generate gaps in the plot
                else
                    simpPlot(z,c) = plot(xAxis+offSet(c),partSNRmean{1}(expIdx{p}==z,condIdx(c)),'o-','color',condColors(condIdx(c),:),'LineWidth',2);
                    set(simpPlot(z,c), 'MarkerFaceColor', get(simpPlot(z,c), 'Color'),'markerSize',4);
                    e_H(c,z) = errorbar(xAxis+offSet(c),partSNRmean{1}(expIdx{p}==z,condIdx(c)),partSNRerr{1}(expIdx{p}==z,condIdx(c))./2,'LineWidth',2,'color',condColors(condIdx(c),:)); % divide by two because of errorbar.m
                    %curSig = double(partVecP{1}(expIdx{p}==z,c)<pCutOff);
                    %curSig(curSig==0) = NaN; % use NaNs to generate gaps in the plot
                    %tempH = plot(xAxis,curSig*sigPos(c),[lColor{c},'sq'],'linewidth',2);
                    %set(tempH, 'MarkerFaceColor', get(tempH, 'Color'),'markerSize',5);
                end
                tempH = plot(xAxis,curSig*sigPos(c),'sq','color',condColors(condIdx(c),:),'linewidth',2);
                set(tempH, 'MarkerFaceColor', get(tempH, 'Color'),'markerSize',5);
            end
            xTicks = [xTicks,xAxis];
            xOrder = [xOrder,find(expIdx{p}==z)];
            if z==max(expIdx{p})
                set(gca,'xtick',xTicks,'xticklabel',roiSelection(xOrder),'ytick',yMin:yUnit:yMax,gcaOpts{:});
                drawnow;
                rotateXLabels(gca,45);
                xlim([min(xTicks)-0.5,max(xTicks)+0.5]);
                ylim([yMin,yMax]);
                lH = legend([simpPlot(z,:),sigH],[condNames(condIdx),{['p<',num2str(pCutOff)]}],'location',legLoc);
                lPos = get(lH,'position');
                lPos(1) = lPos(1)+legOffset(1);
                lPos(2) = lPos(2)+legOffset(2);
                lPos(3) = lPos(3)+legOffset(3);
                lPos(4) = lPos(4)+legOffset(4);
                set(lH,'position',lPos);
                legend boxoff
                ylabel('Signed Vector Average Amplitude');
                arrayfun(@(x) errorbar_tick(x,50),e_H(e_H>0));
            else
            end
            hold off
        else
        end
    end
    figWidth=[16,8];
    set(gcf, 'units', 'centimeters'); % make figure size units centimeters
    oldPos = get(gcf, 'pos');
    newPos = oldPos;
    newPos(3) = figWidth(1);
    newPos(4) = figWidth(2);
    set(gcf, 'pos',newPos);
    if p ==1
        figName = [figFolder,'/simplePlot'];
    else
        figName = [figFolder,'/simplePlot_benoit'];
    end
    %export_fig([figName,'.pdf'],'-pdf','-transparent',simpFig(p));
    hold off
end


%% MAKE RING PLOT
close all
sigColors = rgb2hsv(condColors);
sigColors(:,2) = sigColors(:,2)*.2;
sigColors = hsv2rgb(sigColors);


% compute ring values
ringSet = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'V1ring')),roiList,'uni',false)));
loVal = roiList{ringSet(1)}(7:strfind(roiList{ringSet(1)},'-')-1);
hiVal = roiList{ringSet(1)}(strfind(roiList{ringSet(1)},'-')+1:end);
hiVal2 = roiList{ringSet(2)}(strfind(roiList{ringSet(2)},'-')+1:end);
% assign values
ringSize = str2double(hiVal)-str2double(loVal); 
ringIncr = str2double(hiVal2)-str2double(hiVal);
ringMax = str2double(roiList{ringSet(end)}(strfind(roiList{ringSet(end)},'-')+1:end));
ringList = (ringSize/2):ringIncr:(ringMax-ringSize/2);
figure;
subFigLabel = {'A','B','C','D'};
lStyle = {'o','sq'};
numRingPlots = 3; % V1, V2, V3
numConds = 3;

xMin = 0; xMax = 6.25;
plotCount = 0;
for c=1:(length(condNames)-1)
    if c==1
        yMax = 2;
        yMin = -3;
        yUnit = 1;
        figure(1);
        figName = [figFolder,'/ringROI'];
    elseif c<4
        yMax = 0.4;
        yMin = -0.4;
        yUnit = 0.2;
    else
        if c==4
            figure(2);
            figName = [figFolder,'/ringROI_benoit'];
            plotCount = 0;
        else
        end
        yMax = 0.6;
        yMin = -0.4;
        yUnit = 0.2;
    end
    sigPos = yMin+(yMax-yMin)*.05;
    plotCount = plotCount+1;
    for r=1:numRingPlots
        %if c==3; c=7; else end
        subplot(numConds+1,numRingPlots,r+(plotCount-1)*numRingPlots);
        hold on;
        plot(ringList,zeros(size(ringList)),'color',[.5 .5 .5],'linewidth',2);
        plot(ones(1,10)*2,linspace(-10,10,10),'color',[.5 .5 .5],'linewidth',2);
        ringH(plotCount,r)=plot(ringList,ringVecMean(:,c,r).*ringSign(:,c,r),'o', 'color',condColors(c,:),'linewidth',2,'markersize',5);
        plot(ringList,ringVecHigh(:,c,r).*ringSign(:,c,r),'-', 'color',condColors(c,:),'linewidth',.5);
        plot(ringList,ringVecLow(:,c,r).*ringSign(:,c,r),'-', 'color',condColors(c,:),'linewidth',.5);
        curSig = ringVecP(:,c,r)<pCutOff;
        sigRegions = bwlabel(curSig);
        
        for z=1:max(sigRegions)
            minIdx = find(sigRegions == z,1,'first');
            maxIdx = find(sigRegions == z,1,'last');
            xPatch = [ringList(minIdx)-(ringSize/4),ringList(minIdx)-(ringSize/4),ringList(maxIdx)+(ringSize/4),ringList(maxIdx)+(ringSize/4)];
            yPatch = [yMin,yMax,yMax,yMin];
            patchH = patch(xPatch,yPatch,sigColors(c,:),'edgecolor','none');
            uistack(patchH, 'bottom')
        end

                
        
        %curSig(curSig==0) = NaN; % use NaNs to generate gaps in the plot
        %tempH = plot(ringList,curSig*sigPos,'o', 'color',condColors(c,:),'linewidth',10);
        set(gca,'xtick',0:1:7,'ytick',yMin:yUnit:yMax,gcaOpts{:},'ticklength',[.025,.025]); 
        ylim([yMin,yMax]); xlim([xMin,xMax]);
        if r==numRingPlots && plotCount == 3
            lH = legend(ringH(:,numRingPlots),condNames(c-2:c),'location','northeast');
            lPos = get(lH, 'position');
            lPos(1:2) = lPos(1:2)+0.09;
            set(lH, 'position',lPos);
        elseif r==1
            text(-1,yMax*1.3,subFigLabel{plotCount},'fontsize',24,'fontname','arial');
            if plotCount == 3
                text(-1,yMin+yMin*.9,subFigLabel{4},'fontsize',24,'fontname','arial');
            else
            end
        else
        end
        text(.25,yMax*1.1,evcROIs{r},'fontsize',18,'fontname','arial');
        if r==1 && plotCount==2
           ylabel('Signed Vector Average Amplitude','fontsize',18,'fontname','arial'); 
        elseif r==2 && plotCount==3
           xlabel('Eccentricity at Ring Center (º/vis. angle)','fontsize',18,'fontname','arial'); 
        else
        end
        hold off
        % fit plots
        subplot(numConds+1,numRingPlots,plotCount+numConds*numRingPlots);
        hold on;
        if r==1 
            %plot(linspace(yMin,yMax),lineFit(1,plotCount).*linspace(yMin,yMax),'color',[.5 .5 .5],'linewidth',2);
            plot(linspace(yMin,yMax),linePerp(plotCount).*linspace(yMin,yMax),'--','color',[.5 .5 .5],'linewidth',2);
        else
        end
        for ri = 1:size(ringComplex,1)
            if ringSign(ri,c,r)==1
                plot(squeeze(real(ringComplex(ri,c,r))),squeeze(imag(ringComplex(ri,c,r))),lStyle{1},'color',condColors(c,:))
            else
                plot(squeeze(real(ringComplex(ri,c,r))),squeeze(imag(ringComplex(ri,c,r))),lStyle{2},'color',condColors(c,:))
            end
        end
        xlim([yMin,yMax]);ylim([yMin,yMax])
        xlabel('real part'); ylabel('imaginary part'); 
        set(gca,'xtick',yMin:yUnit:yMax,'ytick',yMin:yUnit:yMax,gcaOpts{:},'ticklength',[.025,.025]);
        axis square 
    end
    if c==3 || c==6 || c==7
        figWidth=[30,30];
        set(gcf, 'units', 'centimeters'); % make figure size units centimeters
        oldPos = get(gcf, 'pos');
        newPos = oldPos;
        newPos(3) = figWidth(1);
        newPos(4) = figWidth(2);
        set(gcf, 'pos',newPos);
        export_fig([figName,'.pdf'],'-pdf','-transparent',gcf);
        hold off
    else
    end
end


%% MAKE FOVEA PLOT
close all
yMax = [2,.6,.3,.3,.3,.3];
yMin = [-1.5,-.2,-.1,-.1,-.1,-.1];
yUnit = [.5,.2,.1,.1,.1,.1];
sigPos = yMax-(yMax-yMin)*.05;
numConds = length(condNames)-1;
for r=1:length(evcROIs)
    for c=1:numConds
         if c==1
            figure(1);
            figName = [figFolder,'/foveaROI'];
            plotCount = 0;
         elseif c==4
            figure(2);
            figName = [figFolder,'/foveaROI_benoit'];
            plotCount = 0;
         end
        plotCount = plotCount+1;
        subplot(3,length(evcROIs),r+(plotCount-1)*length(evcROIs));
        roiToPlot = find(~isnan(foveaVecMean(:,c,r)));
        hold on;
        plot(0:length(foveaList)+1,zeros(1,length(foveaList)+2),'k','linewidth',2);
        tempErrLow = (foveaVecMean(roiToPlot,c,r)-foveaVecLow(roiToPlot,c,r)).*foveaSign(roiToPlot,c,r);
        tempErrHigh = (foveaVecHigh(roiToPlot,c,r)-foveaVecMean(roiToPlot,c,r)).*foveaSign(roiToPlot,c,r);
        ha = errorbar(roiToPlot,foveaVecMean(roiToPlot,c,r).*foveaSign(roiToPlot,c,r),tempErrLow,tempErrHigh,'-','color',condColors(c,:),'linewidth',2);
        hb = get(ha,'children');  
        Xdata = get(hb(2),'Xdata');
        temp = 4:3:length(Xdata);
        temp(3:3:end) = [];
        % xleft and xright contain the indices of the left and right
        %  endpoints of the horizontal lines
        xleft = temp; xright = temp+1; 
        % Increase line length by 0.2 units
        Xdata(xleft) = Xdata(xleft) - .05;
        Xdata(xright) = Xdata(xright) + .05;
        set(hb(2),'Xdata',Xdata)
        curSig = double(foveaVecP(roiToPlot,c,r)<pCutOff);
        curSig(curSig~=1) = NaN; % use NaNs to generate gaps in the plot
        tempH = plot(roiToPlot,curSig*sigPos(c),'sq','color',condColors(c,:),'linewidth',2);
        set(tempH, 'MarkerFaceColor', get(tempH, 'Color'),'markerSize',5);
        if plotCount==1 
            title(evcROIs{r},'fontsize',fontSize,'fontname','Arial');
        else
        end
        if plotCount==3
            set(gca,'xtick',1:length(foveaList),'xticklabel',{'inner','border','outer'},'ytick',yMin(c):yUnit(c):yMax(c),gcaOpts{:},'ticklength',[.05,.05]);
        else
            if plotCount==2 && r == 1
                ylabel('Signed Vector Average Amplitude','fontsize',fontSize,'fontname','Arial');
            else
            end
            set(gca,'xtick',4,'ytick',yMin(c):yUnit(c):yMax(c),gcaOpts{:},'ticklength',[.05,.05]); 
        end
        ylim([yMin(c),yMax(c)]);
        xlim([.5,3.5]);
        hold off
    end
end
figWidth=[length(evcROIs)*6,12];
set(gcf, 'units', 'centimeters'); % make figure size units centimeters
oldPos = get(gcf, 'pos');
newPos = oldPos;
newPos(3) = figWidth(1);
newPos(4) = figWidth(2);
set(gcf, 'pos',newPos);
export_fig([figName,'.pdf'],'-pdf','-transparent',gcf);

    
