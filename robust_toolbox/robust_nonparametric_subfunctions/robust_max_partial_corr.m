function [rrob,prob,Ymaxeffect,whY,Xadj,Yadj] = robust_max_partial_corr(X,Y,doplots,doprint)
% function [rrob,prob,Ymaxeffect,whY,Xadj,Yadj] = robust_max_partial_corr(X,Y,doplots,doprint)
%
% tor wager, aug. 06
%
% finds maximum partial corr. for each column of X in a set of data columns
% Y.  Plots are optional (default = no plots, no stat printout to screen).
%
% Example:
% Y = cl_pos(1).CONTRAST.all_data;
% X = R.X;  %3 columns, 3rd is intercept
% [rrob,prob,Ymaxeffect,whY,Xadj,Yadj] = robust_max_partial_corr(X,Y,1);

warning('off','stats:statrobustfit:IterationLimit');

if nargin < 3, doplots = 0; end
if nargin < 4, doprint = 0; end

% run robust reg. on each col of Y
[b,t,p,sig,f,fp,fsig,stat] = robust_reg_matrix(X,Y,1);

% get voxels x regs matrix of max effect
whmat = abs(t') == repmat(max(abs(t'),[],1),size(t,2),1);

% equivalent:
% whmat = (p') == repmat(min((p'),[],1),size(p,2),1);

% whY is which Y vector shows max effect
[whY,Yindx] = find(whmat);
if size(whY,2) < length(whY), whY = whY'; end
if size(Yindx,2) < length(Yindx),Yindx = Yindx'; end

% Obs. x regressors matrix of Y data with max effect for each reg.
Ymaxeffect = Y(:,whY);

% partial correlations
% rrob and prob show max partial corr for each regressor
for i = Yindx
    [Xadj(:,i),Yadj(:,i),rrob(i),prob(i)] = partialcor(X,Ymaxeffect(:,i),i,doprint);
end


% plotting stuff
% ---------------------------------------
if doplots
    wh_intercept = find_intercept(X);
    k = length(Yindx);
    tor_fig(1,k)

    for i = Yindx
        subplot(1,k,i)
        if i == wh_intercept
            barplot_columns(Yadj(:,i),[],[],'nofig','dorob');
        else
            plot_correlation(X,Ymaxeffect(:,i),'col',i,'robust');
        end

        ylabel(['Y column ' num2str(whY(i))])
        xlabel(['Regressor ' num2str(i)]);
    end

end


warning('on','stats:statrobustfit:IterationLimit');

return


% -------------------------------------------------------------------
% COMPUTE: misc functions
% -------------------------------------------------------------------

% -------------------------------------------------------------------

function val = find_intercept(X)

val = find(~any(diff(X)));

if isempty(val), error('X should contain an intercept.'); end
return


