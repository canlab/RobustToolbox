function [br,tr,pr,se] = robust_ttest2(d,grp)
% [mean,t,p,se] = robust_ttest2(d,grp)
%
% takes a data matrix d and returns a robust 2-sample t-test
% on each column, as estimated with IRLS
%
% grp is a model matrix
% grouping vector, should be [1 -1] contrast
%
% also t and p values for Ho: u = 0 on each column.
%
% tor wager



for i = 1:size(d,2)
    
    try
        [br(i),stats] = robustfit(ones(size(d,1),1),d(:,i),'bisquare',[],'off');
        tr(i) = stats.t;
        pr(i) = stats.p;
        se(i) = stats.se;
    catch
        tr(i) = NaN;
        pr(i) = NaN;
        se(i) = NaN;
        br(i) = NaN;
    end
        
end