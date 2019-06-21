clear all
close all
setenv('DYLD_LIBRARY_PATH','')

doSmooth = true;
experiment = 3;

switch experiment
    case 1
        topDir = '/Volumes/Denali_MRI/kohler/fMRI_EXP/MULTIFOVEA';
        conditions = {'cont','mofo','disp'};
        suffix = 'vr.sc.dt';
        if doSmooth
            surfSuffix = [suffix,'.3fwhm'];
        else
            surfSuffix = suffix;
        end
    case 3
        topDir = '/Volumes/Denali_MRI/kohler/fMRI_EXP/MULTIFOVEA_ATT';
        conditions = {'motAtt','dispAtt'};
        suffix = 'vr.sc.dt_al';
        if doSmooth
            surfSuffix = [suffix,'_3fwhm'];
        else
            surfSuffix = suffix;
        end
    otherwise
        msg = sprintf('\nunknown experiment %0.0d\n',experiment);
        error(mgs);
end


subjFolders = subfolders([topDir,'/201*'],1);
for s =1:length(subjFolders);
    for c = 1:length(conditions)
        volInput = subfiles(...
            sprintf('%s/run*%s*%s.nii.gz',subjFolders{s},conditions{c},suffix),1);
        outSurf = [subjFolders{s},'/SURF'];
        volOutput = sprintf('%s/%s_signal.nii.gz',outSurf,conditions{c});
        if volInput{1}(1)~=0 && ~exist(volOutput,'file')
            mriBrainFFT(volInput,volOutput,10);
        else
        end
        runSurfDir = [subjFolders{s},'/SURF/run_surf'];
        if ~exist(outSurf,'dir')
            mkdir(outSurf);
        elseif ~exist(runSurfDir,'dir')
            mkdir(runSurfDir);
        end
        hemi = {'lh','rh'};
        for h = 1:length(hemi)
            surfInput = subfiles(sprintf('%s/%s.std.141.run*%s*%s.niml.dset',runSurfDir,hemi{h},conditions{c},surfSuffix),1);
            surfOutput = sprintf('%s/%s.std.141.%s_signal_%s.runSmooth.niml.dset',outSurf,hemi{h},conditions{c},surfSuffix);
            if surfInput{1}(1)~=0 && ~exist(surfOutput,'file')
                mriBrainFFT(surfInput,surfOutput,10);
            else
            end
        end
    end
end