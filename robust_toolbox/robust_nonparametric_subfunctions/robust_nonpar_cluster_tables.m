function robust_nonpar_cluster_tables(R,doextended,wh_regressor,cl_all,cl_sig,cl_pos,cl_neg,rfieldname,pfieldname,uthr,dat,cl_extent_pos,cl_extent_neg)
% robust_nonpar_cluster_tables(R,doextended,wh_regressor,cl_all,cl_sig,cl_pos,cl_neg,rfieldname,pfieldname,uthr,dat,cl_extent_pos,cl_extent_neg)
% 
%
% tor wager
% used in robust_nonparam_displayregions
% run that high-level function to configure necessary inputs
%
% dat is matrix of image at different thresholds (one thr.image per column)
% [dat,dat_pos,dat_neg,A] = robust_nonparam_orthviews(R,wh_regressor,overlay,showrois,A);
% robust_nonpar_cluster_tables(R,1,1,A.cl_all,A.cl_sig,A.cl_pos,A.cl_neg,A.rfieldname,A.pfieldname,.05,dat,A.cl_uncor_pos,A.cl_uncor_neg);

numc = R.volInfo.c;         % number of clusters
%n = R.volInfo.nvox;
%[k,v] = size(R.correct.t);

if ~isempty(cl_sig)

    % Tables: ROIs with sig values in them
    % ---------------------------------------------------------------------

    fprintf(1,'\n* Note: maxstat is %s','abs(t threshold)')

    fprintf(1,'\n* Significant ROIs                         ')
    fprintf(1,'\n* ========================================================================\n')
    cluster_table(cl_all,1,0,'threshold');

    % ---------------------------------------------------------------------


    % Tables: Significant contig regions within ROIs
    % ---------------------------------------------------------------------

    fprintf(1,'\n* Significant regions within ROIs                         ')
    fprintf(1,'\n* ========================================================================\n')
    if doextended
        fprintf(1,'\n> Positive                         ')
        extended_table(cl_pos,R,rfieldname,pfieldname)
        fprintf(1,'\n> Negative                         ')
        extended_table(cl_neg,R,rfieldname,pfieldname)
    else
        fprintf(1,'\n> Positive                         ')
        cluster_table(cl_pos,1,0,'threshold','from_cluster');
        fprintf(1,'\n> Negative                         ')
        cluster_table(cl_neg,1,0,'threshold','from_cluster');
    end

    % ---------------------------------------------------------------------


    % Tables: Extent at uncorrected thresholds
    % ---------------------------------------------------------------------

    fprintf(1,'\nExtent at uncorrected thresholds (numbers are voxels corrected or at this uncorr. thresh. or higher)\n')
    
    % add number significant at uncorrected thresholds and print table
    nthresh = length(uthr);

    ndat = size(dat,2);
    if ndat == nthresh + 1
        % we have significant 'corrected' in 1st col. of dat
        addcol = 1;
    elseif ndat == ndat
        % we have uncorrected thresh only in dat
        addcol = 0;
    else
        error('number of thresholds does not match size of dat input.')
    end
        
    for j = 1:nthresh
        dats = sum(dat(:,1:j+addcol)>0,2);   % signficant uncor. voxels at this threshold

        for c = 1:numc

            sigvox = R.svc.tsig{c}(:,wh_regressor);

            % get sig voxels (dats) adjacent to at least one
            dats_pruned = iimg_cluster_prune(dats,sigvox,R.volInfo);

            nsig(c,j) = sum(dats_pruned);  % num sig at this threshold in or adjacent

        end

    end

    % header
    fprintf(1,'\nCluster\t')
    for nt = 1:nthresh
        fprintf(1,'%s %3.3f\t','N @ ', uthr(nt))
    end
    fprintf(1,'\n');

    % body
    for c = 1:numc
        if any(nsig(c,:))
            fprintf(1,'Cl. %3.0f\t', c)
            fprintf(1,repmat('%3.0f\t',1,nthresh),nsig(c,:))
            fprintf(1,'\n');
        end
    end
    fprintf(1,'\n');

    % ---------------------------------------------------------------------


    fprintf(1,'\n* Significant regions : uncorrected p < %3.3f within ROIs  ',uthr(end))
    fprintf(1,'\n* ========================================================================\n')
    fprintf(1,'\n> Positive                         ')
    cluster_table(cl_extent_pos,1,0);
    fprintf(1,'\n> Negative                         ')
    cluster_table(cl_extent_neg,1,0);

end

return



% -------------------------------------------------------------------
% TABLE: extended output table for a set of clusters
% -------------------------------------------------------------------

function extended_table(cl_sig,R,rfieldname,pfieldname)
%cluster_table(cl_sig,1,0,'p','threshold','from_cluster');
str = 'cluster_table(cl_sig,1,0,''p'',''threshold'',''from_cluster''';
for myregressor = 1:size(R.X,2)
    str = [str ',''' rfieldname{myregressor} ''',''' pfieldname{myregressor} ''''];
end
str = [str ');'];
eval(str)

return

