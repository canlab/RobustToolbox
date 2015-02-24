function [t,tr,tvdiff,pdiff] = rob_vs_ols_matrix(d)
% [t,tr,zdiff,pdiff] = rob_vs_ols_matrix(d)
%
% takes a data matrix d and does t-tests on the columns, robust and OLS
% plots output.

b = mean(d);
t = mean(d) ./ ste(d);
p = 2*(1 - tcdf(t,size(d,1)-1));

for i = 1:size(d,2)
    
    z = (d(:,i) - mean(d(:,i))) ./ std(d(:,i)); wh = find(abs(z) >= 3);
    tmp = d(:,i); tmp(wh) = [];
    
    % trimmed univariate outliers - OLS
    br2(i) = mean(tmp); tr2(i) = mean(tmp) ./ ste(tmp); pr2(i) = 2*(1 - tcdf(tr2(i),size(tmp,1)-1));
    
    % Least Trimmed Squares (univariate)
    %[res]=fastltsm_noplot(ones(size(d,1),1),d(:,i),struct('intercept',0,'alpha',.8));
    %br2(i) = res.coefficients; tr2(i) = res.coefficients ./ (res.scale ./ sqrt(size(d,1)));
    %pr2(i) = 2*(1 - tcdf(tr2(i),size(d,1)-1));
    
    
    % IRLS
    
    [br(i),stats] = robustfit(ones(size(d,1),1),d(:,i),'bisquare',[],'off');
    tr(i) = stats.t;
    pr(i) = stats.p;
    
    
end

df = (size(d,1) - 1);
sigma = sqrt(  (1/df) + (1/df)  );
z1 = spm_t2z(t,df);
z2 = spm_t2z(tr,df);        % spm_t2z progressively underestimates z for large values (> 5, starts at >3)
zdiff = abs((z1 - z2) ./ sigma);
pdiff = 2 * (1 - normcdf(zdiff));   % 2-tailed
% which is correct?

[mn,v] = tstat(df);
tdiff = abs(t - tr);
sigma = sqrt(  (v/df) + (v/df)  );  % something is wrong with this one, but I'm not sure what
tdiff = tdiff ./ sigma;
pdiff = 2 * (1 - normcdf(tdiff));   % 2-tailed

tv = t ./ sqrt(v); trv = tr ./ sqrt(v);
tvdiff = abs(tv - trv); % z-scores by dividing by t-distribution variance
pdiff = 2 * (1 - normcdf(tvdiff));   % 2-tailed

tv = t ./ sqrt(v); tr2v = tr2 ./ sqrt(v);
tv2diff = abs(tv - tr2v); % z-scores by dividing by t-distribution variance
pdiff2 = 2 * (1 - normcdf(tv2diff));   % 2-tailed


figure('Color','w'); set(gca,'FontSize',18); hold on;

bar([t' tr' tr2']); colormap gray
%tor_bar_steplot(t,ste(t),{'k'},-.13)
%tor_bar_steplot(rmean(:,2)',rste(:,2)',{'k'},.13)
%set(gca,'YLim',[.5 1])
title('T-values for OLS and robust methods')
ylabel('t-score')
xlabel('Task Condition')
set(gca,'XLim',[0 size(d,2)+1],'XTick',1:size(d,2),'XTickLabel',1:size(d,2))

for i = 1:size(d,2), 
    if pdiff(i) < .001, 
        text(i-.25,max([t(i) tr(i)]) + .1,'***','FontSize',18);,
    elseif pdiff(i) < .01, 
        text(i-.15,max([t(i) tr(i)]) + .1,'**','FontSize',18);,
    elseif pdiff(i) < .05, 
        text(i,max([t(i) tr(i)]) + .1,'*','FontSize',18);,
    end
end

for i = 1:size(d,2), 
    if pdiff2(i) < .001, 
        text(i-.1,max([t(i) tr2(i)]) + .1,'***','FontSize',18);,
    elseif pdiff2(i) < .01, 
        text(i+.1,max([t(i) tr2(i)]) + .1,'**','FontSize',18);,
    elseif pdiff2(i) < .05, 
        text(i+.2,max([t(i) tr2(i)]) + .1,'*','FontSize',18);,
    end
end
    
legend({'OLS' 'Robust IRLS' 'Univar'})

return



