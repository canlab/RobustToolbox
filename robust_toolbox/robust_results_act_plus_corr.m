% res = robust_results_act_plus_corr(u1, u2, k1, k2, [covfield], [['mask', mask], ['overlay', overlay], ['writeimgs', 0|1]])
%
% gets activation blobs that overlap with at least one covariate blob, and
% vice versa
% generates output orthviews and tables
%
% tor wager

function res = robust_results_act_plus_corr(u1, u2, k1, k2, covfield, varargin)
    writeimgs = 0;
    overlay = [];

    if(~exist('covfield', 'var') || isempty(covfield)), covfield = 'cov1'; end
    for i=1:2:length(varargin)
        switch(varargin{i})
            case 'mask'
                mask = varargin{i+1};
            case 'overlay'
                overlay = varargin{i+1};
            case 'writeimgs'
                writeimgs = varargin{i+1};
        end
    end

    % -------------------------------------------------------------------------
    % Threshold activation and covariate images
    % -------------------------------------------------------------------------

    % first image: Activation
    [res.act.clpos, res.act.clneg, res.act.dat, volInfo] = robust_results_threshold(u1, k1, 'mask', mask, 'overlay', overlay);

    % second image: Correlation
    [res.(covfield).clpos, res.(covfield).clneg, res.(covfield).dat] = robust_results_threshold(u2, k2, covfield, 'mask', mask, 'overlay', overlay);

    % -------------------------------------------------------------------------
    % Find which clusters have any overlapping values
    % -------------------------------------------------------------------------

    % intersection of both
    [res.act.cloverlap, res.(covfield).cloverlap, res.act.overlapdat, res.(covfield).overlapdat] = ...
        iimg_cluster_intersect(res.act.dat, res.(covfield).dat, volInfo);

    % save max stats in each set of clusters for later printout in table
    [clindx, nvox] = iimg_cluster_index(res.act.overlapdat, volInfo.xyzlist');
    for i = 1:length(res.act.cloverlap)
        %%%get touching clusters for this overlap area
        whcl = (clindx == i);

        res.act.cloverlap(i).act_tmax = get_maxstat_from_dat(whcl, res.act.dat, volInfo);
        res.act.cloverlap(i).cov_tmax = get_maxstat_from_dat(whcl, res.(covfield).dat, volInfo);
    end

    [clindx, nvox] = iimg_cluster_index(res.(covfield).overlapdat, volInfo.xyzlist');
    for i = 1:length(res.(covfield).cloverlap)
        %%%get touching clusters for this overlap area
        whcl = (clindx == i);

        res.(covfield).cloverlap(i).act_tmax = get_maxstat_from_dat(whcl, res.act.dat, volInfo);
        res.(covfield).cloverlap(i).cov_tmax = get_maxstat_from_dat(whcl, res.(covfield).dat, volInfo);
    end

    %cluster_orthviews(cl1, 'bivalent');
    %cluster_orthviews(cl2, 'bivalent');

    % print tables
    fprintf(1, '\n\nActivated clusters overlapping with correlated clusters\n')
    fprintf(1, 'Height threshold is p < %3.6f, Extent is %3.0f contiguous voxels\n', u1, k1)
    fprintf(1, '---------------------------------------------------------------------\n\n')

    cluster_table(res.act.cloverlap, 0, 0, 'act_tmax', 'cov_tmax');

    fprintf(1, '\n\nCorrelated clusters overlapping with correlated clusters\n')
    fprintf(1, 'Height threshold is p < %3.6f, Extent is %3.0f contiguous voxels\n', u2, k2)
    fprintf(1, '---------------------------------------------------------------------\n\n')

    cluster_table(res.(covfield).cloverlap, 0, 0, 'act_tmax', 'cov_tmax');

    fprintf(1, '\n\n')

    % -------------------------------------------------------------------------
    % Find overlapping voxels
    % -------------------------------------------------------------------------

    if writeimgs
        iname = ['act_u' num2str(u1) 'k' num2str(k1) '_' covfield '_u' num2str(u2) 'k' num2str(k2) '_overlap.img'];
        [int_dat] = iimg_intersection(res.act.overlapdat, res.(covfield).overlapdat, volInfo, 'name', iname, 'posneg');
    else
        [int_dat] = iimg_intersection(res.act.overlapdat, res.(covfield).overlapdat, volInfo, 'posneg');
    end

    % get clusters
    res.cl_pospos = iimg_indx2clusters(int_dat(:, 1), volInfo);
    res.cl_posneg = iimg_indx2clusters(int_dat(:, 2), volInfo);
    res.cl_negpos = iimg_indx2clusters(int_dat(:, 3), volInfo);
    res.cl_negneg = iimg_indx2clusters(int_dat(:, 4), volInfo);
    
    
    cluster_orthviews(res.cl_pospos, {[1 0 0]}, 'overlay', overlay);
    cluster_orthviews(res.cl_posneg, {[1 .5 0]}, 'add');
    cluster_orthviews(res.cl_negpos, {[0 .5 1]}, 'add');
    cluster_orthviews(res.cl_negneg, {[0 0 1]}, 'add');
    
    combined_cls = {res.cl_pospos, res.cl_posneg, res.cl_negpos, res.cl_negneg};
    cl_colors = {[1 0 0] [1 .5 0] [0 .5 1] [0 0 1]};
    whempty = cellfun(@isempty, combined_cls);
    combined_cls(whempty) = [];
    cl_colors(whempty) = [];
    montage_clusters(overlay, combined_cls{:}, cl_colors);
    montage_clusters_medial(overlay, combined_cls{:}, cl_colors);

    res.volInfo = volInfo;

    iname = ['act_u' num2str(u1) 'k' num2str(k1) '_' covfield '_u' num2str(u2) 'k' num2str(k2) '_overlap'];
    fprintf(1, '\nDone: %s\n', iname);
    fprintf(1, 'res.act contains activation clusters and data vectors for + and - voxels.\n');
    fprintf(1, 'res.%s contains covariate clusters and data vectors for + and - voxels.\n', covfield);
    fprintf(1, 'res.act.cloverlap contains activation clusters overlapping with covariate.\n');
    fprintf(1, 'res.%s.cloverlap contains covariate clusters overlapping with activation.\n', covfield);
    fprintf(1, 'res.volInfo contains mask volume information\n');
    fprintf(1, '\n\n')
end



function maxstat = get_maxstat(cl)
    % max absolute value
    if size(cl.Z, 2) ~= size(cl.XYZmm, 2)
        cl.Z = cl.Z';
    end
    [maxabs, whmax] = max(abs(cl.Z(1, :)));
    maxstat = cl.Z(1, whmax(1));
end


function maxstat = get_maxstat_from_dat(whcl, dat, volInfo)
    [tmp, tmp, tmp, touching_this_cl] = ...
        iimg_cluster_intersect(whcl, dat, volInfo);
    touching_this_cl( touching_this_cl == 0) = [];
    [maxabs, whmax] = max(abs(touching_this_cl));
    maxstat = touching_this_cl(whmax(1));
end

