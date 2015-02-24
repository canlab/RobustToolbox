function varargout = robust_htw_results_plots(meth, varargin)
%
% This function is a specialized function for the scnlab htw estimation
% tools
% For use with robust regression directories specifically
% IN this framework, first level models must be run with SPM 5/8, and 
% apply_derivative_boost.m must be used to reconstruct htw (height, time to
% peak, amplitude images) in each first-level directory.
% This function should actually work for any robust results directory with
% SPM models with basis functions at the first level, though.
% It does not use robust reg weights (yet), unfortunately.
%
% This function makes plots of reconstructed HRFs across subjects for significant
% regions stored in clusters .mat files
%
% This function has several methods, and can be run in the following
% sequence of steps:
%
%cd('robust0001'); % change to a robust results directory
% robust_htw_results_plots('setup');  % Gather image names and append to SETUP.mat
% robust_htw_results_plots('define'); % Define conditions to plot/average over and append to  SETUP.mat
% clpos_data = robust_htw_results_plots('extract', clpos_data);  % extract data from images
% robust_htw_results_plots('plot', clpos_data);  % Plot and save pngs


    load SETUP
    imgs = SETUP.files;
    n = size(imgs, 1);
    
    switch meth
        case 'setup'
        % ===================================================
% Get list of files with image names

    for i = 1:n

        dd = fileparts(SETUP.files(i, :));

        subjdirs{i} = dd;
        
        imglistnames{i} = fullfile(dd, 'db_amplitude_names.mat');

        if ~exist(imglistnames{i}, 'file')
            disp(['Warning!  Cannot find ' imglistnames{i}])
        end

        spmmatnames{i} = fullfile(dd, 'SPM.mat');

        if ~exist(imglistnames{i}, 'file')
            disp(['Warning!  Cannot find ' spmmatnames{i}])
        end
    end

    SETUP.HTWinfo.subjdirs = subjdirs;
    SETUP.HTWinfo.imglistnames = imglistnames;
    SETUP.HTWinfo.spmmatnames = spmmatnames;

    %% Load sample SPM mat and get trial types
    for i = 1:n

        [SETUP.HTWinfo.trialtypes{i}, SETUP.HTWinfo.nsess{i}, SETUP.HTWinfo.ntrialtypes{i}, SETUP.HTWinfo.bf{i}, SETUP.HTWinfo.betas_of_interest{i}, SETUP.HTWinfo.beta_names{i}] ...
            = get_spm_htw_info(spmmatnames{i});

    end

    save SETUP -append SETUP
    %%

        case 'define'  % define trial types to average and plot
        % ===================================================
        disp(SETUP.HTWinfo.trialtypes{1})
        emptytostop = 1;
        conditions = zeros(SETUP.HTWinfo.ntrialtypes{1}, 1);
        
        indx = 1;
        gook = 1;
        while gook
            disp(['Enter vector of ones and zeros for condition ' num2str(indx) 'in brackets (e.g., [1 0 0 1 0 0]), or return to exit']);
            disp('Entries of 1 on the same row will be averaged over in plots.');
            
            v = input(': ');
            if isempty(v)
            gook = 0;
            elseif length(v) ~= SETUP.HTWinfo.ntrialtypes{1}, disp('Invalid length.'); 
            else
                disp('Accepted.');
                conditions(:, indx) = v ./ sum(abs(v));
                indx = indx + 1;
            end
        end
        
        disp('Entered conditions (weights for weighted sum in each condition to plot): '); 
        SETUP.HTWinfo.conditions = conditions;
        print_matrix(conditions);
        
        for i = 1:size(conditions, 2)
            SETUP.HTWinfo.condnames{i} = input(['Enter string name for contrast ' num2str(i)], 's');
        end
        
        disp(' '); 
        disp('Adding contrast vectors to SETUP.HTWinfo and saving'); 
        for i = 1:n
            SETUP.HTWinfo.contrastvectors{i} = repmat(conditions, SETUP.HTWinfo.nsess{i}, 1); 
        end
        
        save SETUP -append SETUP
    
        
        case 'extract'  % extract data from clusters
        % ===================================================
        cl = varargin{1};
        
        fprintf('Extracting betas and reshaping: Subject %03d', 0);
        
        for i = 1:n
            fprintf('\b\b\b%03d', i);
            
            v = SETUP.HTWinfo.contrastvectors{i};
            cli = tor_extract_rois(SETUP.HTWinfo.beta_names{i}, cl);

            for cc = 1:length(cl)
                cl(cc).HTW.betas{i} = cli(cc).timeseries;

                % Parse by basis function
                nbf = size(SETUP.HTWinfo.bf{i}, 2);
                cl(cc).HTW.betas_by_bf{i} = cell(1, nbf);
                for bb = 1:nbf
                    cl(cc).HTW.betas_by_bf{i}{bb} = cl(cc).HTW.betas{i}(bb:nbf:end);
                end

                % Average : basis function x condition (rows x cols)
                for bb = 1:nbf , bf_x_condavg(:, bb) = (cl(cc).HTW.betas_by_bf{i}{bb}' * v)'; end
                cl(cc).HTW.bf_x_condavg{i} = bf_x_condavg;

                % Fits
                cl(cc).HTW.fit{i} = SETUP.HTWinfo.bf{i} * cl(cc).HTW.bf_x_condavg{i};
                cl(cc).HTW.condnames = SETUP.HTWinfo.condnames;

                for jj = 1:size(v, 2) % for each condition
                    cl(cc).HTW.group_fit{jj}(i, :) = cl(cc).HTW.fit{i}(:, jj)';
                end

               
                
            end % cc, cluster loop
        end % i, subject loop
        
        for i = 1:n
            for cc = 1:length(cl)
                v = SETUP.HTWinfo.contrastvectors{i}; % same for all subjects

                % Weights
                for jj = 1:size(v, 2) % for each condition
                    [avg,t,p,se, w] = robust_mean(mean(cl(cc).HTW.group_fit{jj}, 2));
                    w(isnan(w)) = 0;
                    w = w ./ sum(w);
                    cl(cc).HTW.weights{jj} = w;
                end

                % Weighted averages
                for jj = 1:size(v, 2) % for each condition
                   % cl(cc).HTW.weighted_avg{jj} = cl(cc).HTW.group_fit{jj} * cl(cc).HTW.weights{jj};
                end
            end
        end
            
        varargout{1} = cl;

        case 'plot'  % extract data from clusters
            % ===================================================
            dirname = input('Enter dir name to save pngs: ', 's');
            mkdir(fullfile(pwd, dirname))

            cl = varargin{1};

            colors = {'b' 'g' 'r'};

            for cc = 1:length(cl)
                f1 = create_figure('HRF plot');
                for jj = 1:length(cl(cc).HTW.group_fit)
                    hh(jj) = tor_fill_steplot(cl(cc).HTW.group_fit{jj}, colors(jj), 0);
                end

                legend(hh, cl(cc).HTW.condnames)

                rlabel = ['Region' num2str(cc)];
                savename = fullfile(pwd, dirname, [rlabel '.png']);

                title([dirname ' ' rlabel]);
                xlabel('Time');
                ylabel('BOLD units');

                saveas(f1, savename);

            end
            
        otherwise error('Unknown method input string.')

    end % case


end % main function


function [trialtypes, nsess, ntrialtypes, bf, betas_of_interest, beta_names] =  get_spm_htw_info(spmmatname)

    load(spmmatname)
    trialtypes = [];
    nsess = length(SPM.Sess);
    ntrialtypes = length(SPM.Sess(1).U);

    for i = 1:ntrialtypes
    trialtypes = strvcat(trialtypes,  SPM.Sess(1).U(i).name{1});
    end

    disp(trialtypes)

    % Basis set
    bf = SPM.xBF.bf;
    nbf = size(bf, 2);

    % Betas of interest
    betas_of_interest = [];

    for i = 1:nsess
        %trialtypes = strvcat(trialtypes,  SPM.Sess(1).U(i).name{1});
        ntrialtypes = length(SPM.Sess(i).U);

        betas_of_interest = [betas_of_interest SPM.Sess(i).col(1:ntrialtypes * nbf)];
    end

    create_figure; imagesc(zscore(SPM.xX.X(:, betas_of_interest))); set(gca, 'YDir', 'reverse'); axis tight; colormap gray; drawnow
    title('All columns of interest -- verify that these correspond to conditions of interest')

    % Get list of beta names
    % -------------------------------------------------
    [dd, ff, ee] = fileparts(spmmatname);


    beta_names = cell(length(betas_of_interest), 1);
    for i = 1:length(betas_of_interest)
        beta_names{i} = sprintf('%s%sbeta_%04d.img', dd, filesep, betas_of_interest(i));
    end

    beta_names = char(beta_names{:});

end % get_spm_htw_info

%%

