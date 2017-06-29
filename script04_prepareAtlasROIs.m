function roiReport = script04_prepareAtlasROIs(varargin)    
    % Description:	main function for analysis of multicue daata
    % 
    % Syntax:	script04_prepareAtlasROIs(<options>)
    % <options>
    %   doExp    - scalar indication whether to run Exp. 1 or [2]
    %    
    %   doKGS    - logical indicating whether to include KGS ROIs
    %                   (true), or not([false])
    %
    %   saveROIs - logical indicating whether to save the ROIs (true)
    %                   or just run diagnostics ([false])
    %
    %   hemiOverlap - logical indicating whether to 
    %                 allow L/R overlap ([true]) or not (false)
    %
    %   doMask - logical indicating whether to 
    %                 create new mask ([true]) or not (false)
    %
    %   clipFrac    - scalar indicating the clipping fraction for 3dAutoMask [0.4]   

    %% ADD PATHS   
    codeFolder = '/Users/kohler/code';
    addpath(genpath([codeFolder,'/git/schlegel/matlab_lib']));
    addpath(genpath([codeFolder,'/git/mrC']));
    addpath(genpath([codeFolder,'/git/MRI/matlab']))
    setenv('DYLD_LIBRARY_PATH','')
    
    %% RUN FUNCTION
    opt	= ParseArgs(varargin,...
            'doExp', 1, ...
            'doKGS', false, ...
            'saveROIs', true, ...
            'hemiOverlap', true, ...
            'clipFrac',0.5, ...
            'doMask',true ...
            );
    
    if opt.doExp == 2
        disp('running Experiment 2 (collected by B. Cottereau)');
        topFolder = '/Volumes/Denali_4D2/kohler/fMRI_EXP/CS_DISP';
        subjFolders = subfolders([topFolder,'/2012*'],1);
    elseif opt.doExp == 1
        disp('running Experiment 1 (collected by P.J. Kohler)');
        topFolder = '/Volumes/Denali_MRI/kohler/fMRI_EXP/MULTIFOVEA';
        subjFolders = subfolders([topFolder,'/201*'],1);
    else
        error('only two experiments!')
    end
    
    atlasName{1}='wangatlas';
    if opt.doKGS
        disp('including KGS ROIs');
        atlasName{2}='KGSatlas';
        roiNames{2} = {'IOG' 'OTS' 'mFUS' 'pFUS' 'PPA' 'VWFA1' 'VWFA2'}; % KGS ROIs, in order
    else
    end

    roiNames{1} = {'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
        'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
        'IPS5' 'SPL1' 'FEF'}; % WANG ROIs, in order

    for s = 1:length(subjFolders)
        if isempty(strfind(subjFolders{s},'nl-0033'))
            subName = subjFolders{s}(max(strfind(subjFolders{1},'_'))+1:end);
        else
            subName = 'skeri0003';
        end
        if opt.doMask
            rmCall = sprintf('rm %s/ROIs/%s_mask.nii.gz',subjFolders{s},subName);
            maskCall = sprintf('3dAutomask -prefix %s/ROIs/%s_mask.nii.gz -clfrac %.1f %s/ref_epi.ts.do.vr.nii.gz',subjFolders{s},subName,opt.clipFrac,subjFolders{s});
            system(join({rmCall,maskCall},';'));
        else
        end
        
        for a=1:length(atlasName)
            % determine how much is lost in the masking
            nii = NIfTI.Read([subjFolders{s},'/ROIs/rh.',atlasName{a},'_rs.nii.gz']);
            rightData = nii.data;
            nii = NIfTI.Read([subjFolders{s},'/ROIs/lh.',atlasName{a},'_rs.nii.gz']);
            leftData = nii.data;
            nii = NIfTI.Read([subjFolders{s},'/ROIs/',subName,'_mask.nii.gz']);
            maskData = nii.data;

            if a==1
                finalNii = nii;
            else
            end
            finalROIs = zeros(size(leftData));
            
            for z=1:length(roiNames{a})
                numVox(z,1) = length(find(maskData==1 & leftData==z)); % number of voxels in final ROI
                numVox(z,2) = length(find(maskData==1 & rightData==z)); % number of voxels in final ROI
                percVox(z,1) = 1-((length(find(leftData==z))-numVox(z,1))./length(find(leftData==z)));
                percVox(z,2) = 1-((length(find(rightData==z))-numVox(z,2))./length(find(rightData==z)));
                percVox(isnan(percVox)) = 0;
                %percVox(percVox>1) = 1;
                roiReport.numOutside{a}(z,s) = ( length(find(leftData==z))-numVox(z,1) ) + ( length(find(rightData==z))-numVox(z,2) );
                if mean(percVox(z,:),2)>(3/4) && sum(numVox(z,:),2)>30 % if less than 25% was loss across both ROIs                
                    finalROIs( maskData==1 & leftData==z ) = z;
                    finalROIs( maskData==1 & rightData==z ) = z+length(roiNames{a});
                    if ~opt.hemiOverlap
                        % take out L/R overlap
                        finalROIs( maskData==1 & rightData==z & leftData==z ) = 0;
                    else
                    end
                    roiReport.included{a}(z,s) = 1;
                    roiReport.hemiOverlap{a}(z,s) =  length( find ( maskData==1 & rightData==z & leftData==z ) ); 
                else
                    roiReport.included{a}(z,s) = 0;
                    roiReport.hemiOverlap{a}(z,s) = NaN;
                end
            end
            if a==1
                wangROIs = finalROIs;
            else
                kgsROIs = finalROIs;
                kgsROIs(kgsROIs>0) = kgsROIs(kgsROIs>0)+length(roiNames{1}); % add 50 to KGS rois
            end
            % save individual sets
            if opt.saveROIs 
                finalNii.data = finalROIs;
                outputName = [subjFolders{s},'/ROIs/',subName,'_',atlasName{a},'_al_mask.nii.gz'];
                if exist(outputName,'file')
                    delete(outputName)
                else
                end
                NIfTI.Write(finalNii,outputName);
            else
            end
            clear left*
            clear right*
            roiReport.num{a}(:,s) = sum(numVox,2);
            roiReport.perc{a}(:,s) = mean(percVox,2);
            roiReport.names{a} = roiNames{a};
            clear percVox; clear numVox;
        end
    end
end
            
        
        

        
    
    %tempData = niftiRead('/Volumes/Denali_4D2/kohler/fMRI_EXP/SYM_4GROUPS/20140825_nl-0025'/'nl-0014_wangatlas_al_mask.nii.gz');