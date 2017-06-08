clear all
close all

topFolder = '/Volumes/Denali_MRI/kohler/fMRI_EXP/MULTIFOVEA';
subjFolders = subfolders([topFolder,'/201*'],1);
atlasName{1}='wangatlas';
%atlasName{2}='KGSatlas';

roiNames{1} = {'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'}; % WANG ROIs, in order
roiNames{2} = {'IOG' 'OTS' 'mFUS' 'pFUS' 'PPA' 'VWFA1' 'VWFA2'}; % KGS ROIs, in order

for s=1:length(subjFolders)
    if isempty(strfind(subjFolders{s},'nl-0033'))
        subName = subjFolders{s}(max(strfind(subjFolders{1},'_'))+1:end);
    else
        subName = 'skeri0003';
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
        if s==7
            disp('hello');
        else
        end
        for z=1:length(roiNames{a})
            numVox(z,1) = length(find(maskData==1 & leftData==z)); % number of voxels in final ROI
            numVox(z,2) = length(find(maskData==1 & rightData==z)); % number of voxels in final ROI
            percVox(z,1) = 1-((length(find(leftData==z))-numVox(z,1))./length(find(leftData==z)));
            percVox(z,2) = 1-((length(find(rightData==z))-numVox(z,2))./length(find(rightData==z)));
            percVox(isnan(percVox)) = 0;
            %percVox(percVox>1) = 1;
            if mean(percVox(z,:),2)>(3/4) && sum(numVox(z,:),2)>30 % if less than 25% was loss across both ROIs                
                finalROIs(maskData==1 & leftData==z) = z;
                finalROIs(maskData==1 & rightData==z) = z+length(roiNames{a});
                included{a}(z,s) = 1;
            else
                included{a}(z,s) = 0;
            end
        end
        finalROIs(find(leftData >0 & rightData>0)) = 0; % take out L/R overlap
        overlap(s,a) = length(find(leftData >0 & rightData>0));
        if a==1
            wangROIs = finalROIs;
        else
            kgsROIs = finalROIs;
            kgsROIs(kgsROIs>0) = kgsROIs(kgsROIs>0)+length(roiNames{1}); % add 50 to KGS rois
        end
        % save individual sets
        finalNii.data = finalROIs;
        outputName = [subjFolders{s},'/ROIs/',subName,'_',atlasName{a},'_al_mask.nii.gz'];
        if exist(outputName,'file')
            delete(outputName)
        else
        end
        NIfTI.Write(finalNii,outputName);
        clear left*
        clear right*
        allNum{a}(:,s) = sum(numVox,2);
        allPerc{a}(:,s) = mean(percVox,2);
        clear percVox; clear numVox;
    end
    
    % combine ROIs
%     comboROIs = zeros(size(finalROIs));
%     comboROIs = wangROIs+kgsROIs;
%     comboROIs(find(kgsROIs >0 & wangROIs>0)) = 0;
%     setOverlap(s) = length(find(kgsROIs >0 & wangROIs>0));
%     
%     finalNii.data = comboROIs;
%     outputName = [subjFolders{s},'/ROIs/',subName,'_wangKGS_al_mask.nii.gz'];
%     if exist(outputName,'file')
%         delete(outputName)
%     else
%     end
%     NIfTI.Write(finalNii,outputName);
end
   
            
        
        

        
    
    %tempData = niftiRead('/Volumes/Denali_4D2/kohler/fMRI_EXP/SYM_4GROUPS/20140825_nl-0025'/'nl-0014_wangatlas_al_mask.nii.gz');