function [o2, t] = robust_results_threshold_2014(varargin)
% Makes maps of results from a CANlab Robust Regression Toolbox directory
% IN DEVELOPMENT - VALID BUT NOT COMPLETE IN TERMS OF OPTIONS, ETC.
% Other robust_results methods are also valid.
%
% Usage:
% -------------------------------------------------------------------------
% o2 = robust_results_threshold_2014([optional inputs])
%
% Author and copyright information:
% -------------------------------------------------------------------------
%     Copyright (C) 2014  Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Inputs:
% -------------------------------------------------------------------------
% p            Followed by (uncorrected) p-value threshold
% k            Followed by cluster extent threshold (k contiguous voxels)
%
% Outputs:
% -------------------------------------------------------------------------
% o2           fmridisplay object with handles/info about display elements.
% t            A series of statistic_image objects for all t-maps in dir
%
% Examples:
% -------------------------------------------------------------------------
%
% give examples here
%
% See also:
% * list other functions related to this one, and alternatives*

% Programmers' notes:
% List dates and changes here, and author of changes

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

% Defaults
% -----------------------------------
pthr = .005;
kthr = 1;


% optional inputs with default values
% -----------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            
            case {'p', 'pthr'}, pthr = varargin{i+1}; varargin{i+1} = [];
            case {'k', 'kthr'}, kthr = varargin{i+1}; varargin{i+1} = [];
                    
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% -------------------------------------------------------------------------
% LOAD SETUP AND ALL MAPS
% -------------------------------------------------------------------------

[dd, dirname] = fileparts(pwd);

setupfile = fullfile(pwd, 'SETUP.mat');

if ~exist(setupfile, 'file')
    error('You must be in a Robust Toolbox results directory. SETUP.mat is missing.')
else   
    load(setupfile)
end

names = {'Intercept (group activation)'};

tnames{1} = fullfile(pwd, 'rob_tmap_0001.img');   %dir(fullfile(pwd, 'rob_tmap*img'));
pnames{1} = fullfile(pwd, 'rob_p_0001.img');   %dir(fullfile(pwd, 'rob_tmap*img'));

indx = 2;
while exist(tnames{end}, 'file')
    tnames{end + 1} = fullfile(pwd, sprintf('rob_tmap_%04d.img', indx));
    pnames{end + 1} = fullfile(pwd, sprintf('rob_p_%04d.img', indx));
    names{end + 1} = sprintf('Cov %d', indx-1);
  indx = indx + 1;
end
tnames(end) = [];
pnames(end) = [];
names(end) = [];

if isempty(tnames), error('You must be in a Robust Toolbox results directory. No valid t-maps.'), end

n = length(tnames);

for i = 1:n
    t{i} = statistic_image('image_names', tnames{i}, 'type', 'T', 'dfe', SETUP.dfe);
    p{i} = statistic_image('image_names', pnames{i});
    
    if size(t{i}.dat, 1) ~= size(p{i}.dat, 1), error('T and p images do not match in size/num voxels.'); end
    
    % replace p-values with robust regression ones - correct p-vals
    
    t{i}.p = p{i}.dat;
end

clear p

% -------------------------------------------------------------------------
% THRESHOLD AND DISPLAY ALL MAPS
% -------------------------------------------------------------------------

o2 = [];

for i = 1:n
t{i} = threshold(t{i}, pthr, 'unc', 'k', kthr);
o2 = canlab_results_fmridisplay(t{i}, o2, 'nooutline', 'addmontages');

title([dirname ' ' names{i}], 'FontSize', 18)
end


end % function
