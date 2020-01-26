function robust_results_make_invp_maps(varargin)
% robust_results_make_invp_maps([options])
% DESCRIPTION
%   run in robfit directory
%   takes rob_p_* maps and makes invpmaps/*_invp maps
%   if EXPT.SNPM.connames and/or EXPT.covnames exist, will use in naming
%       otherwise: invpmaps/lcon<i>_hcon<j>_invp.nii
%   converts NaNs to 0s in output maps
% OPTIONS
%   'space', <img>
%       resample inverse p maps to space of <img>
%       make <img> blank string for no resampling
%       DEFAULT: avg152T1.nii
%   'd', <robfitdir>
%       DEFAULT: pwd
%
% NOTES
%   to view, try the shell function ~ruzicl/scripts/invpview (wrapper for FSLview)
%
% Luka Ruzic (2012)
%

RESAMPLE = which('avg152T1.nii');
WD = pwd;

i=1;
while i<=numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'space'
                i=i+1;
                RESAMPLE = varargin{i};
            case 'd'
                i=i+1;
                WD = varargin{i};             
            otherwise
                error(sprintf('Unrecognized argument: %s',varargin{i}))                
        end
        
    else
        disp(varargin{i})
        error('ABOVE ARGUMENT UNRECOGNIZED')
    end
    i=i+1;
end


cd(WD)

robdirs = filenames('robust[0-9]*','absolute');
if isempty(robdirs), error('No directories: robust[0-9]*'); end

dout = fullfile(pwd,'invpmaps');
if ~exist(dout,'dir'), mkdir(dout); end

lconnames = {};
hconnames = {};
try
    load('EXPT.mat')
    
    if isfield(EXPT.SNPM,'connames')
        lconnames = cellstr(EXPT.SNPM.connames);
        lconnames = regexprep(lconnames,' ','');
        lconnames = regexprep(lconnames,'[^A-Za-z0-9_+.-]','');
        
        fprintf('Lower level contrast names:\n')
        disp(lconnames)
    end
    
    if isfield(EXPT,'covnames')
        hconnames = cellstr(EXPT.covnames);
        hconnames = regexprep(hconnames,' ','');
        hconnames = regexprep(hconnames,'[^A-Za-z0-9_+.-]','');
        
        fprintf('Higher level covariate names:\n')
        disp(hconnames)
    end
end

fprintf('... CREATING inverse p maps\n')
fout={};
for i = 1:numel(robdirs)    
    l = regexprep(robdirs{i},'.*/robust','');
    if isempty(lconnames)
        lconname = sprintf('lcon%s',l);
    else
        lconname = lconnames{str2num(l)};
    end
    
    pmaps = filenames(fullfile(robdirs{i},'rob_p_*.img'),'absolute');
    for j = 1:numel(pmaps)
        [dummy f e] = fileparts(pmaps{j});
        
        h = regexprep(f,'rob_p_','');
        if isempty(hconnames)            
            hconname = sprintf('hcon%s',h);
        else
            hconname = hconnames{str2num(h)};
        end

        
        fout{end+1} = fullfile(dout,[lconname '_' hconname '_invp.nii']);
        
        fprintf('%-7s %s\n','pmap',pmaps{j},'invpmap',fout{end});
        
%         [pvol pdat] = iimg_read_img(pmaps{j});
%         [tvol tdat] = iimg_read_img(regexprep(pmaps{j},'_p_','_tmap_'));
        pvol = spm_vol(pmaps{j});
        pdat = spm_read_vols(pvol);
        tvol = spm_vol(regexprep(pmaps{j},'_p_','_tmap_'));
        tdat = spm_read_vols(tvol);
        
        % invert
        pdat = 1 - pdat;
        
        % add valence
        pdat(tdat<0) = pdat(tdat<0) * -1;
        
        % write out file
        pvol.fname = fout{end};
        pvol.descrip = 'IRLS robust regression inverse p values';
        spm_write_vol(pvol,pdat);
%         iimg_write_images(pdat,pvol,fout{end});
        
        if ~isempty(RESAMPLE)
            fprintf('resampling to space of: %s\n',RESAMPLE);
            scn_map_image(fout{end},RESAMPLE,'write',fout{end});
        end        
    end
    fprintf('\n')
end
fprintf('... REMOVING NaNs\n')
nans2zeros(fout);

end

