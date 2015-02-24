function t_image_ztest(P,df)
%
% do Z test comparing reliability of two t images
%
% t_image_ztest(P,df)
% e.g., t_image_ztest([],5)
%
% P is string matrix of 2 t-image file names, or empty to choose with GUI
% df is degrees of freedom going into each t-image (must be same), e.g, df
% = 9 for 10 subject one-sample t-test.
% 
% Warning:
% The z-test here is not truly a test on normally distributed data
% because it divides t-scores by their distribution standard deviation,
% computes differences between standardized t-scores, and then dividides by
% sqrt(2) to correct for number of tests.  
% Simulation shows that this test is too conservative with small df (e.g.,
% p-values are too high by a factor of 1.7 with 10 df, but this drops as the
% threshold gets more extreme).  P-values are accurate with df = 40, but at
% the tails (high threshold, Z = 3) the test may become slightly anticonservative.
%


if isempty(P), P = spm_get(2,'*img','Choose two t-images');,end

V = spm_vol(P); v = spm_read_vols(V);
v(any(v == 0 | isnan(v),4)) = NaN;

 % * calculate Z-test for comparison
 % ----------------------------------------------------
 
 % SIMULATION: shows this test is too conservative
 % [mn,vv] = tstat(df); 
 %x = trnd(df,2,5000)'; x = x ./ vv;
%zdiff = x(:,2) - x(:,1); zdiff = zdiff ./ sqrt(2); 
%[xax,h] = hist(zdiff,50); figure;plot(h,xax)
%xn = randn(5000,1);
%[xnx] = hist(xn,h); hold on; plot(h,xnx,'r')
%figure;
%h = normplot(zdiff);
% SIMULATION SUGGESTS that P-values are too high by a factor of 1.7 with 10
% df
% for z= 1.6:.01:3
%nor(end+1)=sum(xn>z) ./ length(xn);
%zd(end+1)=sum(zdiff>z) ./ length(zdiff);
%end
%figure;plot(nor,'kx');hold on; plot(zd,'ro')
%plot(nor./zd,'gs')
 
        [mn,vv] = tstat(df); 
        
        zdiff = squeeze(v(:,:,:,1)) ./ sqrt(vv) - squeeze(v(:,:,:,2)) ./ sqrt(vv); % z-scores by dividing by t-distribution variance
        zdiff = zdiff ./ sqrt(2);                                % diff between z-scores is distributed with var=sum of var(z1) + var(z2)
                                                            % ...thus,sigma
                                                            % (z1 - z2) = sqrt(2) 
                                                            % e.g., - http://www.mathpages.com/home/kmath046.htm
        pdiff = 2 * (1 - normcdf(abs(zdiff)));              % 2-tailed
         
        
        zthr = zdiff; zthr(abs(zdiff) < 1.96) = NaN;
        
        V = V(1);
    V.fname = fullfile(pwd,'z_difference.img');
    V.descrip = [P(1,:) P(2,:)];
    spm_write_vol(V,zdiff);
    
    V.fname = fullfile(pwd,'p_difference.img');
    V.descrip = [P(1,:) P(2,:)];
    spm_write_vol(V,pdiff);
    
        V.fname = fullfile(pwd,'z_diff_thresholded.img');
    V.descrip = [P(1,:) P(2,:)];
    spm_write_vol(V,zdiff);
    
    try, 
        spm_image('init','z_diff_thresholded.img'); colormap jet
        %cl = mask2clusters('z_diff_thresholded.img');
        %montage_clusters([],cl,{'r'})
        %  icbm_localize(cl)
      catch
      end
      
      
      
      return
      
      
      