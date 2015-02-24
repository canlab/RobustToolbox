function rob_vs_ols_corr1(it,ns)
%rob_vs_ols_corr1(it,ns)
% it = pop size of realizations
% ns = [10 20 30] vector of group sample sizes

%it = 5000;  % population size

nind = 1;
for n = ns %[10 20 30 40]   % number of subjects
    
    g = randn(n,it*2);
    r = []; nout = []; rr = [];
    
    for i = 1:2:size(g,2)
        
        tmp = corrcoef(g(:,i:i+1));
        r(end+1) = tmp(1,2);
        
        
        [res]=fastmcd_noplot([g(:,i) g(:,i+1)]);
        % remove n most extreme outliers and recompute correlation
        
        wh = res.flag==0; nout(end+1) = sum(res.flag==0);
        tmp = g(:,i:i+1); tmp(wh,:) = []; 
        tmp = corrcoef(tmp);
        rr(end+1) = tmp(1,2);
        
        if mod(i,100)==0, fprintf(1,'.'),end
        
    end
    
    figure('Color','w');
    [h1,x]=hist(r,50); h2 = hist(rr,x); bar(h1,'b'); hold on; bar(h2,'r');alpha(.5)
    R{nind} = r; RR{nind} = rr; 
    
    % fpr: false positive rates
    
    [rci,sig]=r2z(r,n,.05);
    fpr(nind) = sum(sig) ./ it;
    
    for i = 1:it
        [tmp,sig(i)] = r2z(rr(i),n-nout(i),.05);
    end
    sig(isnan(sig)) = [];
    fprr(nind) = sum(sig)./ length(sig);
    
    % number of outliers removed
    meanout(nind) = mean(nout);
    
    % fnr: false negative rates
    % not done yet
    
    % can do signal detection - Zhr - Zfar
    
    nind = nind + 1;
    drawnow
end


save robust_ols_corr_output

figure;plot(fpr,'bo-','LineWidth',2)
hold on;plot(fprr,'rs-','LineWidth',2)
legend({'OLS' 'Robust MCD'})
set(gca,'XTick',1:length(ns),'XTickLabel',ns)
saveas(gcf,'rob_vs_ols','fig')

keyboard

return

