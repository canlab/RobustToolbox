% go to results directory for robust000*
df = input('Enter degrees of freedom (n images - k parameters): ');

p = dir('*0001.img'); p = str2mat(p.name)

[P2,P,s,sn] = threshold_imgs(p([5 8],:),tinv(1-.005,df),0,'both');
[p2,p1] = threshold_imgs('irls-ols_z_0001.img',norminv(.95),0,'both');  % this is one tailed...

%h = image_scatterplot(str2mat(P,p1),'avgvs3');
%xlabel('Average OLS and Robust t-value'), ylabel('Z-score of Robust - OLS difference')

figure('Color','w')
compare_filtered_t([],P2(1,:),P2(2,:),p2)

disp('OLS (left), Robust (right), Difference (bottom left)')

