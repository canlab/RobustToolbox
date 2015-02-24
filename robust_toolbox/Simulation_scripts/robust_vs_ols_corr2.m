function rob_vs_ols_corr1(it,ns)
%rob_vs_ols_corr1(it,ns)
% it = pop size of realizations
% ns = [10 20 30] vector of group sample sizes

%it = 5000;  % population size

    
    
nind = 1;
for n = ns %[10 20 30 40]   % number of subjects
for rep = 1:10
    
    g = randn(n,it*2);
    r = []; nout = []; rr = []; r2 = []; r3 = []; t2 = []; t3 = []; t = []; p = []; rt = [];
    rp = []; p2 = []; p3 = []; rt2 =[]; rp2 = []; rt3 =[]; rp3 = [];Mp = [];Mp2 = []; Mp3 = [];
    Mr = []; Mr2 = []; Mr3 = []; Mt = []; Mt2 = []; Mt3 = [];
    
    if n == 30 & rep == 1, figure('Color','w');,end
    
    for i = 1:2:size(g,2)
        
        % null hypothesis
        tmp = corrcoef(g(:,i:i+1));
        r(end+1) = tmp(1,2);
        [B,dev,stat]=glmfit(g(:,i),g(:,i+1));
        t(end+1) = stat.t(2);
        p(end+1) = stat.p(2);
        
        % first is intercept, 2nd is correlation
        [B,stat] = robustfit(g(:,i),g(:,i+1));
        rt(end+1) = stat.t(2);
        rp(end+1) = stat.p(2);
        
        % robust MCD
        [res]=fastmcd_noplot([g(:,i) g(:,i+1)]);
        wh = res.flag==0; nout(end+1) = sum(res.flag==0);
        tmp = g(:,i:i+1); tmp(wh,:) = []; 
        tmp2 = corrcoef(tmp); Mr(end+1) = tmp2(1,2);
        [B,dev,stat]=glmfit(tmp(:,1),tmp(:,2));
        Mt(end+1) = stat.t(2);
        Mp(end+1) = stat.p(2);
        
        
        % alternative - d' = .5?;
        tmp = g(:,i) + 2 * g(:,i+1); tmp2 = corrcoef(tmp,g(:,i));
        r2(end+1) = tmp2(1,2);
        [B,dev,stat]=glmfit(g(:,i),tmp);
        t2(end+1) = stat.t(2);
        p2(end+1) = stat.p(2);
        
        [B,stat] = robustfit(g(:,i),tmp);
        rt2(end+1) = stat.t(2);
        rp2(end+1) = stat.p(2);
        
        [res]=fastmcd_noplot([g(:,i) tmp]);
        wh = res.flag==0; nout(end+1) = sum(res.flag==0);
        tmp3 = [g(:,i) tmp]; tmp3(wh,:) = []; 
        tmp2 = corrcoef(tmp3); Mr2(end+1) = tmp2(1,2);
        [B,dev,stat]=glmfit(tmp3(:,1),tmp3(:,2));
        Mt2(end+1) = stat.t(2);
        Mp2(end+1) = stat.p(2);
        
        
        % alternative with 2 outliers
        tmp(1:2) = tmp(1:2) + 10 * randn(2,1);  tmp2 = corrcoef(tmp,g(:,i));
        r3(end+1) = tmp2(1,2);
        [B,dev,stat]=glmfit(g(:,i),tmp);
        t3(end+1) = stat.t(2);
        p3(end+1) = stat.p(2);
        
        [B,stat] = robustfit(g(:,i),tmp);
        rt3(end+1) = stat.t(2);
        rp3(end+1) = stat.p(2);
        
         [res]=fastmcd_noplot([g(:,i) tmp]);
        wh = res.flag==0; nout(end+1) = sum(res.flag==0);
        tmp3 = [g(:,i) tmp]; tmp3(wh,:) = []; 
        tmp2 = corrcoef(tmp3); Mr3(end+1) = tmp2(1,2);
        clear stat
        [B,dev,stat]=glmfit(tmp3(:,1),tmp3(:,2));
        Mt3(end+1) = stat.t(2);
        Mp3(end+1) = stat.p(2);
        
        if mod(i,100)==0, fprintf(1,'.'),end
        
    end
    if n == 30 & rep == 1,
        figure('Color','w')
    subplot(1,3,1),tmp3 = [];tmp3{1} = [];
    [h1,x]=hist(t,50); h2 = hist(rt,x); bar(h1,'b'); hold on; bar(h2,'r');alpha(.5)
    tmp = 1:10:length(x);, for ii = 1:length(tmp),tmp3{ii}=sprintf('%3.2f',x(tmp(ii)));,end
    set(gca,'XTick',tmp,'XTickLabel',tmp3,'FontSize',16)
    T{nind} = t; RT{nind} = rt; legend({'OLS' 'Robust IRLS' 'Robust MCD'})
    P{nind} = p; RP{nind} = rp;
    title(['No effect: Ho, n = ' num2str(n)])
    
    subplot(1,3,2)
    [h1,x]=hist(t2,50); h2 = hist(rt2,x); bar(h1,'b'); hold on; bar(h2,'r');alpha(.5)
    for ii = 1:length(tmp),tmp3{ii}=sprintf('%3.2f',x(tmp(ii)));,end
    set(gca,'XTick',tmp,'XTickLabel',tmp3,'FontSize',16)
    T2{nind} = t2; RT2{nind} = rt2; %legend({'OLS' 'Robust IRLS'})
    P2{nind} = p2; RP2{nind} = rp2;
    title(['True effect: H1'])
        
    subplot(1,3,3)
    [h1,x]=hist(t3,50); h2 = hist(rt3,x); bar(h1,'b'); hold on; bar(h2,'r');alpha(.5)
    tmp = get(gca,'XTick'); tmp = 1:10:50; 
    for ii = 1:length(tmp),tmp3{ii}=sprintf('%3.2f',x(tmp(ii)));,end
    set(gca,'XTick',tmp,'XTickLabel',tmp3,'FontSize',16)
    T3{nind} = t3; RT3{nind} = rt3; %legend({'OLS' 'Robust IRLS'})
    P3{nind} = p3; RP3{nind} = rp3;
    title(['H1 with 2 outliers'])
    drawnow
    end
    % fpr: false positive rates and false neg rates
    
    fpr(rep,nind) = sum(p < .05) ./ it;
    fprr(rep,nind) = sum(rp < .05)./ (it);
    fprM(rep,nind) = sum(Mp < .05)./ (it);
    
    fpr2(rep,nind) = sum(p2 < .05) ./ it;
    fprr2(rep,nind) = sum(rp2 < .05)./ (it);
    fprM2(rep,nind) = sum(Mp2 < .05)./ (it);
    
    fpr3(rep,nind) = sum(p3 < .05) ./ it;
    fprr3(rep,nind) = sum(rp3 < .05)./ (it);
    fprM3(rep,nind) = sum(Mp3 < .05)./ (it);
    
    % number of outliers removed
    %meanout(nind) = mean(nout);
    
    % fnr: false negative rates
    % not done yet
    
    % can do signal detection - Zhr - Zfar
    
    end % reps

    nind = nind + 1;
    drawnow
end

save robust_ols_corr_output2

figure('Color','w'); subplot(1,3,1);set(gca,'FontSize',16)
plot(mean(fpr),'bo-','LineWidth',2)
hold on;plot(mean(fprr),'rs-','LineWidth',2)
legend({'OLS' 'Robust IRLS' 'Robust MCD'})
set(gca,'XTick',1:length(ns),'XTickLabel',ns)
plot(mean(fprM),'g^-','LineWidth',2)
tor_bar_steplot(mean(fpr),std(fpr)./sqrt(10),{'b'})
tor_bar_steplot(mean(fprr),std(fprr)./sqrt(10),{'r'})


tor_bar_steplot(mean(fprM),std(fprM)./sqrt(10),{'g'})

title('False positive rate','FontSize',16)

hh(1) = subplot(1,3,2);set(gca,'FontSize',16)
plot(mean(fpr2),'bo-','LineWidth',2)
hold on;plot(mean(fprr2),'rs-','LineWidth',2)
set(gca,'XTick',1:length(ns),'XTickLabel',ns)
tor_bar_steplot(mean(fpr2),std(fpr)./sqrt(10),{'b'})
tor_bar_steplot(mean(fprr2),std(fprr)./sqrt(10),{'r'})

plot(mean(fprM2),'g^-','LineWidth',2)
tor_bar_steplot(mean(fprM2),std(fprM2)./sqrt(10),{'g'})

title(['Power, r = ' sprintf('%3.2f',mean(r2))],'FontSize',16)

xlabel('Number of subjects (n)')

hh(2) = subplot(1,3,3);set(gca,'FontSize',16)
plot(mean(fpr3),'bo-','LineWidth',2)
hold on;plot(mean(fprr3),'rs-','LineWidth',2)
plot(mean(fprM3),'g^-','LineWidth',2)
set(gca,'XTick',1:length(ns),'XTickLabel',ns)
tor_bar_steplot(mean(fpr3),std(fpr)./sqrt(10),{'b'})
tor_bar_steplot(mean(fprr3),std(fprr)./sqrt(10),{'r'})


tor_bar_steplot(mean(fprM3),std(fprM3)./sqrt(10),{'g'})

title(['Power: r = ' sprintf('%3.2f',mean(r3)) ' (2 outliers)'],'FontSize',16)
equalize_axes(hh)

saveas(gcf,'rob_vs_ols2','fig')

return

