function [br,tr,pr,se, w] = robust_mean(d)
% [mean,t,p,se, w] = robust_mean(d)
%
% takes a data matrix d and returns the robust mean
% of each column, as estimated with IRLS
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
        
        w(:, i) = stats.w;
    catch
        tr(i) = NaN;
        pr(i) = NaN;
        se(i) = NaN;
        br(i) = NaN;
    end
        
end