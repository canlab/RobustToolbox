n = 20;
iter = 500;

% UNDER THE NULL HYPOTHESIS

clear b, clear t, clear br, clear tr

t1 = clock;
for i = 1:iter

    x = randn(n,1); y = randn(size(x)); X = [x ones(size(x))];
    %[bb,stats]=robustfit(x,y); 
    [res,raw] = fastltsm_noplot(x,y); bb = res.coefficients;
    br(i) = bb(2); tr(i) = 0; %tr(i)=stats.t(2);, 
    [bb,bint,r,rint,stats]=regress(y,X);,b(i)=bb(1); t(i)=stats(2).^.5 * sign(bb(1));
    
end
fprintf(1,['Time elapsed for %3.0f iterations is %4.2f s\n'],iter,etime(clock,t1))

[h,xxb] = hist([b br],20);
[h,xxt] = hist([t tr],20);
xxb = [xxb(1)-xxb(end):mean(diff(xxb)):xxb(1)-mean(diff(xxb)) xxb xxb(end):mean(diff(xxb)):max(xxb).*2];
xxt = [xxt(1)-xxt(end):mean(diff(xxt)):xxt(1)-mean(diff(xxt)) xxt xxt(end):mean(diff(xxt)):max(xxt).*2];

figure('Color','w'); subplot 221; [h] = hist(b,xxb); plot(xxb,h./sum(h)),title('OLS betas')

subplot 222; [h] = hist(br,xxb); plot(xxb,h./sum(h)),title('Robust IRLS betas')

subplot 223; [h] = hist(t,xxt); plot(xxt,h./sum(h)),title('OLS t-scores')

subplot 224; [h] = hist(tr,xxt); plot(xxt,h./sum(h)),title('Robust IRLS t-scores')

diary robust_vs_ols_test1.txt
fprintf(1,'Normal Ho Data, simulation of n = %3.0f\n----------------------------------------------------------------\n',n)
fprintf(1,'Bias under : %3.4f for OLS and %3.4f for IRLS\n',mean(b)-1,mean(br)-1)
fprintf(1,'SE under : %3.4f for OLS and %3.4f for IRLS: %3.2f%% improvement in beta stability for IRLS\n',std(b),std(br),100*(std(b)./std(br)))
fprintf(1,'Mean t-score under : %3.4f for OLS and %3.4f for IRLS: %3.2f%% improvement for IRLS\n',mean(t),mean(tr),100*(mean(t)./mean(tr)))
fprintf(1,'SE of t-scores under : %3.4f for OLS and %3.4f for IRLS: %3.2f%% improvement for IRLS\n',std(t),std(tr),100*(std(t)./std(tr)))
diary off


for i = 1:iter

    x = randn(n,1); y = randn(size(x)); X = [x ones(size(x))];
    
    % add an outlier
    y(1) = y(1) + randn(1) * 10;
    mystd(i) = max((y - mean(y)) ./ std(y(2:end)));
    
    [bb,stats]=robustfit(x,y); br(i) = bb(2); tr(i)=stats.t(2);, 
    [res,raw] = fastltsm_noplot(x,y); bb = res.coefficients;
    br(i) = bb(2); tr(i) = 0; %tr(i)=stats.t(2);, 
    [bb,bint,r,rint,stats]=regress(y,X);,b(i)=bb(1); t(i)=stats(2).^.5 * sign(bb(1));
    
end


hh(1)=subplot(2,2,1); hold on; [h] = hist(b,xxb); plot(xxb,h./sum(h),'r'),title('OLS betas')
legend({'Normal Ho data' 'Ho with 1 outlier'})

hh(2)=subplot(2,2,2); hold on; [h] = hist(br,xxb); plot(xxb,h./sum(h),'r'),title('Robust IRLS betas')

hh(3)=subplot(2,2,3); hold on; [h] = hist(t,xxt); plot(xxt,h./sum(h),'r'),title('OLS t-scores')

hh(4)=subplot(2,2,4); hold on; [h] = hist(tr,xxt); plot(xxt,h./sum(h),'r'),title('Robust IRLS t-scores')

equalize_axes(hh(1:2)); equalize_axes(hh(3:4)); 
legend({'Normal Ho data' 'Ho with 1 outlier'})

diary on
fprintf(1,'Normal Ho Data with Outlier at %3.2f st. deviations, simulation of n = %3.0f\n----------------------------------------------------------------\n',mean(mystd),n)
fprintf(1,'Bias under : %3.4f for OLS and %3.4f for IRLS\n',mean(b)-1,mean(br)-1)
fprintf(1,'SE under : %3.4f for OLS and %3.4f for IRLS: %3.2f%% improvement in beta stability for IRLS\n',std(b),std(br),100*(std(b)./std(br)))
fprintf(1,'Mean t-score under : %3.4f for OLS and %3.4f for IRLS: %3.2f%% improvement for IRLS\n',mean(t),mean(tr),100*(mean(t)./mean(tr)))
fprintf(1,'SE of t-scores under : %3.4f for OLS and %3.4f for IRLS: %3.2f%% improvement for IRLS\n',std(t),std(tr),100*(std(t)./std(tr)))
diary off





% UNDER THE ALTERNATIVE HYPOTHESIS, WITH SNR = 1

clear b, clear t, clear br, clear br2, clear tr

for i = 1:iter

    x = randn(n,1); y = x + randn(size(x)); X = [x ones(size(x))];
    [bb,stats]=robustfit(x,y); br(i) = bb(2); tr(i)=stats.t(2);, 
    [bb,bint,r,rint,stats]=regress(y,X);,b(i)=bb(1); t(i)=stats(2).^.5 * sign(bb(1));
    
end

[h,xxb] = hist([b br],20);
[h,xxt] = hist([t tr],20);
xxb = [xxb(1)-xxb(end):mean(diff(xxb)):xxb(1)-mean(diff(xxb)) xxb xxb(end):mean(diff(xxb)):max(xxb).*2];
xxt = [xxt(1)-xxt(end):mean(diff(xxt)):xxt(1)-mean(diff(xxt)) xxt xxt(end):mean(diff(xxt)):max(xxt).*2];

figure('Color','w'); subplot 221; [h] = hist(b,xxb); plot(xxb,h./sum(h)),title('OLS betas')

subplot 222; [h] = hist(br,xxb); plot(xxb,h./sum(h)),title('Robust IRLS betas')

subplot 223; [h] = hist(t,xxt); plot(xxt,h./sum(h)),title('OLS t-scores')

subplot 224; [h] = hist(tr,xxt); plot(xxt,h./sum(h)),title('Robust IRLS t-scores')

diary on
fprintf(1,'Normal Ha Data, SNR = 1, simulation of n = %3.0f\n----------------------------------------------------------------\n',n)
fprintf(1,'Bias under Ha: %3.4f for OLS and %3.4f for IRLS\n',mean(b)-1,mean(br)-1)
fprintf(1,'SE under Ha: %3.4f for OLS and %3.4f for IRLS: %3.2f%% improvement in beta stability for IRLS\n',std(b),std(br),100*(std(b)./std(br)))
fprintf(1,'Mean t-score under Ha: %3.4f for OLS and %3.4f for IRLS: %3.2f%% improvement for IRLS\n',mean(t),mean(tr),100*(mean(t)./mean(tr)))
fprintf(1,'SE of t-scores under Ha: %3.4f for OLS and %3.4f for IRLS: %3.2f%% improvement for IRLS\n',std(t),std(tr),100*(std(t)./std(tr)))
diary off

% this is the critical one for testing lts

for i = 1:iter

    x = randn(n,1); y = x + randn(size(x)); X = [x ones(size(x))];
    
    % add an outlier
    y(1) = y(1) + randn(1) * 10;
    mystd(i) = max((y - mean(y)) ./ std(y(2:end)));
    
    [bb,stats]=robustfit(x,y); br(i) = bb(2); tr(i)=stats.t(2);, 
    [res,raw] = fastltsm_noplot(x,y); bb = res.coefficients;
    br2(i) = bb(1); br2int(i) = bb(2);
    [bb,bint,r,rint,stats]=regress(y,X);,b(i)=bb(1); t(i)=stats(2).^.5 * sign(bb(1));
    
end

hh(1)=subplot(1,3,1); hold on; [h] = hist(b,xxb); plot(xxb,h./sum(h),'r'),title('OLS betas')
legend({'Normal Ha data' 'Ha with 1 outlier'})

hh(2)=subplot(1,3,2); hold on; [h] = hist(br,xxb); plot(xxb,h./sum(h),'r'),title('Robust IRLS betas')

hh(3)=subplot(1,3,3); hold on; [h] = hist(br2,xxb); plot(xxb,h./sum(h),'r'),title('Robust LTS betas')

%hh(4)=subplot(2,2,4); hold on; [h] = hist(tr,xxt); plot(xxt,h./sum(h),'r'),title('Robust IRLS t-scores')

equalize_axes(hh(1:3)); %equalize_axes(hh(3:4)); 
legend({'Normal Ha data' 'Ha with 1 outlier'})

diary on
fprintf(1,'Normal Ha Data (SNR = 1) with Outlier at %3.2f st. deviations, simulation of n = %3.0f\n----------------------------------------------------------------\n',mean(mystd),n)
fprintf(1,'Bias under Ha: %3.4f for OLS and %3.4f for IRLS\n',mean(b)-1,mean(br)-1)
fprintf(1,'SE under Ha: %3.4f for OLS and %3.4f for IRLS: %3.2f%% improvement in beta stability for IRLS\n',std(b),std(br),100*(std(b)./std(br)))
fprintf(1,'Mean t-score under Ha: %3.4f for OLS and %3.4f for IRLS: %3.2f%% improvement for IRLS\n',mean(t),mean(tr),100*(mean(t)./mean(tr)))
fprintf(1,'SE of t-scores under Ha: %3.4f for OLS and %3.4f for IRLS: %3.2f%% improvement for IRLS\n',std(t),std(tr),100*(std(t)./std(tr)))
diary off
