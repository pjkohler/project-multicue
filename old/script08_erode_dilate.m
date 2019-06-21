% input: roi_ids is an Nx1 vector (where N is the number of vertices)
%        so that roi_ids(k)=j means that node k falls in ROI j.
%        If j==0, then node k is not part of any ROI
%
% (the output from afni_niml_read_simpleroi can be used to set roi_ids)

% Define the surface used to measure distances.
% This surface should be anatomical (not inflated or flat) for more
% accurate measurement of distances between nodes.
% Here a node-wise average of the coordinates of the pial and white matter 
% surface is used; slightly less ideal would be to use the pial or white
% matter surface directly.

clear all;
close all;

topFolder = '/Volumes/Denali_MRI/kohler/fMRI_EXP/MULTIFOVEA';
surfFolder = '/Volumes/svndl/anatomy/FREESURFER_SUBS';
setenv('DYLD_LIBRARY_PATH','');
subjFolders = subfolders([topFolder,'/201608*'],1);

hemi = {'lh','rh'};

% the erosion radius is in the same units as the vertex coordinates; 
% typically this is millimeters
erode_radius=3;

% for the metric to measure distances along the curved surface:
% - 'geodesic' is most accurate but requires the fast marching toolbox, and
%   in some cases, a working mex installation
% - 'dijkstra' is a good approximation for geodesic and does not require a
%   working mex installation
distance_metric='geodesic';

for s=1:length(subjFolders)
    subID = subjFolders{s}(max(strfind(subjFolders{s},'_'))+1:end);
    if strcmp(subID,'nl-0033')
        subID = 'skeri0003';
    else
    end
    for h = 1:length(hemi)
        surf_fn = [surfFolder,'/',subID,'_fs4/SUMA/',hemi{h},'.inflated.asc'];
        [v,f]=surfing_read(surf_fn);

        roi_file = subfiles([subjFolders{s},'/SURF/',hemi{h},'*fovea.niml.roi'],1);
        roiStrct = afni_niml_readsimpleroi(roi_file{1});
        roi_idx = roiStrct{1}.region{1};
        roi_ids = zeros(size(v,1),1);
        roi_ids(roi_idx)=1;

        nv=size(v,1);
        if nv~=numel(roi_ids);
            error('roi_ids must be of size %d x 1', nv);
        end

        non_zero=roi_ids(roi_ids>0);
        fprintf('Before erosion: %d nodes in %d ROIs\n',...
                    numel(non_zero), numel(unique(non_zero)));

        % find nodes on the border of each cluster
        nbrs=surfing_surface_nbrs(f);
        cl=surfing_clusterize(roi_ids,nbrs);
        cl_ids=surfing_cluster_borders(cl,nbrs);
        on_border_node_msk=cl_ids<0;
        on_border_node_ids=find(on_border_node_msk);

        n2f=surfing_nodeidxs2faceidxs(f');

        eroded_roi_ids=roi_ids;
        dilated_roi_ids=roi_ids;
        ring_roi_ids=zeros(size(v,1),1);
        for node_id=on_border_node_ids(:)'
            nbr_ids=surfing_circleROI(v,f,node_id,erode_radius,...
                                        distance_metric,n2f);
            eroded_roi_ids(nbr_ids)=0;
            dilated_roi_ids(nbr_ids)=1;
            ring_roi_ids(nbr_ids)=1;
        end

        eroded_non_zero=eroded_roi_ids(eroded_roi_ids>0);
        dilated_non_zero=dilated_roi_ids(dilated_roi_ids>0);
        ring_non_zero=ring_roi_ids(ring_roi_ids>0);

        fprintf('After erosion: %d nodes in %d ROIs\n',...
            numel(eroded_non_zero), numel(unique(eroded_non_zero)));
        fprintf('After dilation: %d nodes in %d ROIs\n',...
            numel(dilated_non_zero), numel(unique(dilated_non_zero)));
        fprintf('After ring: %d nodes in %d ROIs\n',...
            numel(ring_non_zero), numel(unique(ring_non_zero)));

        ring_idx = find(ring_roi_ids>0);
        dilated_idx = find(dilated_roi_ids>0);
        eroded_idx = find(eroded_roi_ids>0);

        S=struct();
        S.data=[roi_ids,eroded_roi_ids,dilated_roi_ids,ring_roi_ids];
        S.labels={'original','eroded','dilated','ring'};
        afni_niml_writesimple([subjFolders{s},'/SURF/',hemi{h},'.fovea_',num2str(erode_radius),'mm.niml.dset'],S)
    end
end

