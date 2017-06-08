clear all
close all

topDir = '/Volumes/Denali_MRI/kohler/fMRI_EXP/MULTIFOVEA';
suffix = 'vr.sc.dt';
doSmooth = true;

if doSmooth
    surfSuffix = [suffix,'.3fwhm'];
else
    surfSuffix = suffix;
end

subjFolders = subfolders([topDir,'/201*'],1);
for s =1:length(subjFolders);
    CONTinput = subfiles([subjFolders{s},'/run*cont*',suffix,'.nii.gz'],1);
    MOFOinput = subfiles([subjFolders{s},'/run*mofo*',suffix,'.nii.gz'],1);
    DISPinput = subfiles([subjFolders{s},'/run*disp*',suffix,'.nii.gz'],1);
    outSurf = [subjFolders{s},'/SURF'];
    runSurfDir = [subjFolders{s},'/SURF/run_surf'];
    if ~exist(outSurf,'dir')
        mkdir(outSurf);
    elseif ~exist(runSurfDir,'dir')
        mkdir(runSurfDir);
    end
    CONToutput = [outSurf,'/CONT_signal.nii.gz'];
    MOFOoutput = [outSurf,'/MOFO_signal.nii.gz'];
    DISPoutput = [outSurf,'/DISP_signal.nii.gz'];
    if ~exist(CONToutput,'file')
        VOL_FFT(CONTinput,CONToutput,10);
    else
    end
    if ~exist(MOFOoutput,'file')
        VOL_FFT(MOFOinput,MOFOoutput,10);
    else
    end
    if DISPinput{1}(1)~=0 && ~exist(DISPoutput,'file')
           VOL_FFT(DISPinput,DISPoutput,10);
    else
    end
%     
%     topSurf = subfiles(sprintf('%s/*.niml.dset*',subjFolders{s}),1);
%     if topSurf{1}~=0 % % if surf files in the top level folder, move to SURF folder
%         for t=1:length(topSurf);
%             movefile(topSurf{t},[runSurfDir,'/.']);
%         end
%     else
%     end
    hemi = {'lh','rh'};
    for h = 1:length(hemi)
        if ~exist(sprintf('%s/%s.std.141.CONT_signal_%s.runSmooth.niml.dset',outSurf,hemi{h},surfSuffix),'file') && ...
           ~exist(sprintf('%s/%s.std.141.MOFO_signal_%s.runSmooth.niml.dset',outSurf,hemi{h},surfSuffix),'file');
                CONTsurfinput = subfiles(sprintf('%s/%s.std.141.run*cont*%s.niml.dset',runSurfDir,hemi{h},surfSuffix),1);
                MOFOsurfinput = subfiles(sprintf('%s/%s.std.141.run*mofo*%s.niml.dset',runSurfDir,hemi{h},surfSuffix),1);
                CONTsurfoutput = sprintf('%s/%s.std.141.CONT_signal_%s.runSmooth.niml.dset',outSurf,hemi{h},surfSuffix);
                MOFOsurfoutput = sprintf('%s/%s.std.141.MOFO_signal_%s.runSmooth.niml.dset',outSurf,hemi{h},surfSuffix);
                SURF_FFT(CONTsurfinput,CONTsurfoutput,10);
                SURF_FFT(MOFOsurfinput,MOFOsurfoutput,10);
        else
        end
        DISPsurfinput = subfiles(sprintf('%s/%s.std.141.run*disp*%s.niml.dset',runSurfDir,hemi{h},surfSuffix),1);
        if DISPsurfinput{1}(1)~=0 && ~exist(sprintf('%s/%s.std.141.DISP_signal_%s.runSmooth.niml.dset',outSurf,hemi{h},surfSuffix),'file')
            DISPsurfoutput = sprintf('%s/%s.std.141.DISP_signal_%s.runSmooth.niml.dset',outSurf,hemi{h},surfSuffix);
            SURF_FFT(DISPsurfinput,DISPsurfoutput,10);
        else
        end
    end
end