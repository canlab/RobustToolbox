function [allcl] = robseed_results(pthresh,kthresh,df,varargin)
% function [allcl] = robseed_results(pthresh,kthresh,df,['seeds',wh seeds],['overlay',overlayimg],['imagenames',extract_data_img_names])
%
% can extract data, given image names from the corresponding files
%
% Example:
% 5 vox, 15 df (17 subjects-2 params), seed 1 only, Extract data
% cl = robseed_results(.005,5,15,1,EXPT.SNPM.P{1}); 
%
% to plot correlated cluster 2 against seed 1:
% load ../seed_clusters.mat
% tor_fig(1,2);montage_clusters_maxslice([],clcorr(2),{'r'});   % or: canlab_montage_clusters(clcorr(2), 'colors', {'r'})
% subplot(1,2,2); plot_correlation_samefig(EXPT.seeds(:,1),cl(2).timeseries,[],'ko',0,1);
% xlabel('L Amygdala seed region (Look - Reapp Neg)')
%
% seedcl = robseed_results(.005,3,19,'seeds',1,'overlay','../mean_funct.img');

allcl = [];
names = {'rob_tmap_0002.img' 'rob_tmap_0003.img' 'rob_tmap_0004.img' 'rob_tmap_0005.img'};
datanames = [];
ovl = [];

for i = 1:length(varargin)
    if isstr(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'overlay',ovl = varargin{i+1}; varargin{i+1} = [];
            case 'seeds',whseeds = varargin{i+1};
            case 'imagenames',datanames = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if ~exist('whseeds','var'), whseeds = 1:length(names); end

disp('Robseed_results.m')
disp('Running for each of these results images: ')
disp(char(names{whseeds}))
disp(' ')


[dum,dirnm] = fileparts(pwd);

% % whseeds = 1:length(names);
% % if length(varargin) > 0, whseeds = varargin{1}; end
% % 
% % datanames = [];
% % if length(varargin) > 1, datanames = varargin{2}; end
% % 
% % ovl = [];
% % if length(varargin) > 2, ovl = varargin{3}; end


for i = whseeds
    
    figname = ['Contrast ' dirnm ', Seed ' num2str(i) '.'];
    
    if exist(names{i}, 'file') == 2
        disp(names{i}); disp('_________________________________');
        
        % Threshold the robust t-map: |t| > tinv(1 - pthresh, df), extent >= kthresh
        % (same rule as the legacy threshold_imgs(..., 'both') call)
        tcrit = tinv(1 - pthresh, df);
        t = statistic_image('image_names', names{i}, 'type', 'T', 'dfe', df);
        t = threshold(t, [-tcrit tcrit], 'raw-outside', 'k', kthresh);
        r = region(t);
        cl = region2struct(r);   % legacy clusters structure for tor_extract_rois and downstream code
        
        % Axial + midline sagittal montage (replaces montage_clusters and
        % montage_clusters_medial)
        canlab_montage_clusters(cl, 'colors', {'r'}, 'overlay', ovl, 'montagetype', 'compact', 'figname', figname, 'noverbose');
        disp(' ');
        
        
        % extract data, if names entered
        if ~isempty(datanames), cl = tor_extract_rois(datanames,cl); end
            
        allcl{i} = cl;
    end
    
    
end

if isempty(allcl), disp('No seeds found. (Are you in a robseed directory?) Exiting.'); end

return
