clear all
close all
codeFolder = '/Volumes/Denali_4D2/kohler/code';
addpath(genpath([codeFolder,'/matlab/others/afni']));
addpath(genpath([codeFolder,'/matlab/self']));
addpath(genpath([codeFolder,'/matlab/others/fitting']));
addpath(genpath([codeFolder,'/matlab/others/Gerhard']));
addpath(genpath([codeFolder,'/matlab/others/surfing']));



readAFNI = true;
smoothData = true;
runSmooth = true;
topFolder = '/Volumes/Denali_MRI/kohler/fMRI_EXP/MULTIFOVEA';
setenv('DYLD_LIBRARY_PATH','');

hemi = {'lh','rh'};

%if smoothData
%    suffix = '.3fwhm';
%    if runSmooth
        % 
        %oldSuffix = ['_runSurf.vr.sc.dt',suffix];
        suffix = '_vr.sc.dt.3fwhm.runSmooth';
%    else
%    end
%else
%    suffix = '';
%end
subjFolders = subfolders([topFolder,'/201*'],1);
for h = 1:2
    if readAFNI
        clear inputFile;
        for s =1:length(subjFolders);
            inputFile{1} = sprintf('%s/SURF/%s.std.141.CONT_signal%s.niml.dset',subjFolders{s},hemi{h},suffix);
            inputFile{2} = sprintf('%s/SURF/%s.std.141.MOFO_signal%s.niml.dset',subjFolders{s},hemi{h},suffix);
            if exist(sprintf('%s/SURF/%s.std.141.DISP_signal%s.niml.dset',subjFolders{s},hemi{h},suffix),'file');
                inputFile{3} = sprintf('%s/SURF/%s.std.141.DISP_signal%s.niml.dset',subjFolders{s},hemi{h},suffix);
                %oldFile = sprintf('%s/SURF/%s.std.141.DISP_signal%s.niml.dset',subjFolders{s},hemi{h},oldSuffix);
                %if exist(oldFile,'file'); movefile(oldFile,inputFile{3}); else end
            else
            end
            % move old files
            %oldFile = sprintf('%s/SURF/%s.std.141.CONT_signal%s.niml.dset',subjFolders{s},hemi{h},oldSuffix);
            %if exist(oldFile,'file'); movefile(oldFile,inputFile{1}); else end
            %oldFile = sprintf('%s/SURF/%s.std.141.MOFO_signal%s.niml.dset',subjFolders{s},hemi{h},oldSuffix);
            %if exist(oldFile,'file'); movefile(oldFile,inputFile{2}); else end
            
            M=cellfun(@(x) afni_niml_readsimple(x),inputFile,'uni',false);
            numBricks = size(M{1}.data,2);
            %snrM = afni_niml_readsimple(sprintf('%s/SURF/%s.std.141.MULTIFOVEA.SNR%s.niml.dset',subjFolders{s},hemi{h},suffix));
            for z=1:length(M)
                numBricks = size(M{z}.data,2);
                brainData(:,1:numBricks,z,s) = (M{z}.data);
                %brainData(:,numBricks+1,z,s) = snrM.data(:,z);
            end
            if length(M) < 3
                brainData(:,:,3,s) = NaN(size(brainData(:,:,1,s)));
            else
            end
            %if s==length(subjFolders)
            %    origM = M{1};
            %else
            %end        
        end
        if h == 1
            save(sprintf('/Volumes/Denali_MRI/kohler/fMRI_EXP/MULTIFOVEA/scripts/leftSurfData%s.mat',suffix),'brainData');
        else
            save(sprintf('/Volumes/Denali_MRI/kohler/fMRI_EXP/MULTIFOVEA/scripts/rightSurfData%s.mat',suffix),'brainData');
        end
    else
        if h == 1
            load(sprintf('/Volumes/Denali_MRI/kohler/fMRI_EXP/MULTIFOVEA/scripts/leftSurfData%s.mat',suffix),'brainData');
        else
            load(sprintf('/Volumes/Denali_MRI/kohler/fMRI_EXP/MULTIFOVEA/scripts/rightSurfData%s.mat',suffix),'brainData');
        end
    end
    groupNames = {'CONT','MOFO','DISP'};
    for c=1:length(groupNames)
        disp(['populating dataset:',hemi{h},'-',groupNames{c},' starting ...']);
        include = find(~isnan(brainData(1,5,c,:)));
        snrTemp = squeeze(brainData(:,2,c,include))-1;  % snr is the second value, subtract one to compare to zerocd 
        realTemp = squeeze(brainData(:,5,c,include)); % real part
        imagTemp = squeeze(brainData(:,6,c,include)); % imaginary part
        
        % FIX PHASE
        tempComplex = complex(realTemp,imagTemp);
        tempMag = abs(tempComplex);
        tempAngle = angle(tempComplex);
        if c==3
            tempAngle = tempAngle+pi;     % condition 3 is 180 degrees out of phase
        else
        end
        tempAngle = tempAngle-pi/3;                 % pi/3 radians = 60 degrees = 4 seconds (pre-TR 16, rather than 12, seconds)
        tempAngle = pkPhaseMod(tempAngle); % adjust to sine base
        tempAngle(tempAngle<0) = 2*pi+tempAngle(tempAngle<0); % take care of wraparound
        
        tempComplex = tempMag.*exp(1i*tempAngle); % convert back to complex numbers
        realTemp = real(tempComplex);
        imagTemp = imag(tempComplex);
        
        % find sign
        if c==1 % fit to CONT condition
            fitVals = complex(mean(realTemp,2),mean(imagTemp,2)); % fit based on means across subjects
        else
        end
        fitVals = fitVals(:);
        lineFit = polyfitZero(real(fitVals(:)),imag(fitVals(:)),1);
        linePerp = 1/(-lineFit(1));
        signTemp = sign(mean(imagTemp,2) - linePerp*mean(realTemp,2));
        
        % replot individual subject data
        subCount = 0;
        for s =1:length(subjFolders)
            if  ismember(s,include) % if subject has data
                subCount = subCount+1;
                subOutData{s}(:,1+(c-1)*3) = abs(complex(realTemp(:,subCount),imagTemp(:,subCount))); % amp
                subOutLabel{s}(1+(c-1)*3) = {[groupNames{c},': AMP']};
                subOutData{s}(:,2+(c-1)*3) = tempAngle(:,subCount); % phase
                subOutLabel{s}(2+(c-1)*3) = {[groupNames{c},': PHA']};
                subOutData{s}(:,3+(c-1)*3) = snrTemp(:,subCount); % snr
                subOutLabel{s}(3+(c-1)*3) = {[groupNames{c},': SNR']};
            else
            end
            if c==length(groupNames) % save data
                subId = subjFolders{s}(max(strfind(subjFolders{s},'_'))+1:end);
                if strcmp(subId,'nl-0033');
                    subId = 'skeri0003';
                else
                end
                subOutName{s} = sprintf('%s/group_surf/%s.std.141.%s.COMBINED%s.niml.dset',topFolder,hemi{h},subId,suffix);
                subS = struct();
                subS.data = subOutData{s};
                subS.labels = subOutLabel{s};
                subS.stats = repmat({'none'},1,length(subOutLabel{s}));
                afni_niml_writesimple(subOutName{s},subS);
            else
            end
        end

        % populate data set
        readyData(:,1) = abs(complex(mean(realTemp,2),mean(imagTemp,2))); % vector mean
        newLabels(1) = {[groupNames{c},': AmpVecMean']};
        newStats(1) = {'none'};
        for v=1:size(realTemp,1)
            if sum(sum([realTemp(v,:)',imagTemp(v,:)']<=0))<(length(subjFolders)*2)
                %ampBounds(v,:) =fitErrorEllipse([realTemp(v,:)',imagTemp(v,:)'],'95CI',false);
                vecStats(v) = tSquaredFourierCoefs([realTemp(v,:)',imagTemp(v,:)'],[],0.05);
                vecPvals(v) = cat(2,vecStats(v).pVal);
                vecTsquared(v) = cat(2,vecStats(v).tSqrd);
            else
                %ampBounds(v,:) = [0,0];
                vecPvals(v) = 1;
            end
        end
        
        tempPhase = arrayfun(@(x) phase(complex(mean(realTemp(x,:),2),mean(imagTemp(x,:),2))),1:size(realTemp,1)); % mean phase
        tempPhase(tempPhase<0) = 2*pi+tempPhase(tempPhase<0); % take care of wraparound
        readyData(:,2) = tempPhase;
        newLabels(2) = {[groupNames{c},': PhaVecMean']};
        newStats(2) = {'none'};
        
        readyData(:,3) = (1-vecPvals)*1000;
        newLabels(3) = {[groupNames{c},': vecPvals']};
        newStats(3) = {'none'};

        
        [H,snrP,CI,STATS] = ttest(snrTemp,0,'alpha',0.05,'dim',2,'tail','right');
        readyData(:,4) = nanmean(snrTemp,2);
        newLabels(4) = {[groupNames{c},': SNR']};
        newStats(4) = {'none'};
        
        readyData(:,5) = STATS.tstat;
        newLabels(5) = {[groupNames{c},': Tscore']};
        newStats(5) = {sprintf('Ttest(%d)',length(subjFolders))};
        
        readyData(:,6) = signTemp;
        newLabels(6) = {[groupNames{c},': SignVec']};

        disp(['populating dataset:',hemi{h},'-',groupNames{c},' done!']);
        
        newM = struct();
        newM.data = readyData;
        newM.labels = newLabels;
        newM.stats = newStats;
        outName = sprintf('%s/group_surf/%s.%s.std.141.vectorResults%s_new.niml.dset',topFolder,hemi{h},groupNames{c},suffix);
        afni_niml_writesimple(outName,newM); 
        clear readyData;
    end
end