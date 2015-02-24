function robust_results2(EXPT,u,k)
%function robust_results2(EXPT,u,k)
% u and k are uncorrected threshold or 'FDR' and extent threshold
% FDR IS NOW VALID ONLY FOR INCREASES

warning off

% convert u to t-threshold
df = length(EXPT.subjects)-length(dir('rob_p_0*.img'));
fprintf(1,'robust_results2 thinks there are %3.0f degrees of freedom.',df);

uo = u;

% if uncorrected threshold
if ~isstr(u)
    t = tinv(1-u,df);
    str = sprintf('Height thresh: t = %3.2f (%3.0f Ss, %3.0f df @ p < %3.4f, extent = %3.0f\n',t,length(EXPT.subjects),df,uo,k);
    disp(str)
else
    uo = 0.05;
end

ovl = [];
if isfield(EXPT,'overlay');
    disp(['Using overlay image: ' EXPT.overlay]);
    ovl = EXPT.overlay;
end


% intercept results
% -------------------------------
fprintf(1,'\n---------------------------------------------\nOverall pos activation\n---------------------------------------------\n');

%FDR threshold
if isstr(u)
    if strcmp(u,'FDR')
        [u2,Ps,Ts] = spm_uc_FDR(.05,[1 df],'T',1,spm_vol('rob_tmap_0001.img'),0);
        t = u2;
        str = sprintf('Height thresh FDR-corr: t = %3.2f (%3.0f Ss, %3.0f df @ p < %3.4f, extent = %3.0f\n',t,length(EXPT.subjects),df,uo,k);
        disp(str)
    end
end

[P2,P,sigmat] = threshold_imgs('rob_tmap_0001.img',t,[k],['pos']);
d=fileparts(P2); [d,f]=fileparts(d);
cl = mask2clusters(P2);
eval(['save ' f '_intercept_pos cl'])

fprintf(1,'\n---------------------------------------------\nOverall neg activation\n---------------------------------------------\n')
[P2,P,sigmat] = threshold_imgs('rob_tmap_0001.img',t,[k],['neg']);
cl = mask2clusters(P2);
eval(['save ' f '_intercept_neg cl'])



if exist('rob_tmap_0002.img') == 2
    
    fprintf(1,'\n---------------------------------------------\nCovariate pos effect\n---------------------------------------------\n')
    %FDR threshold
    if isstr(u)
        if strcmp(u,'FDR')
            [u2,Ps,Ts] = spm_uc_FDR(.05,[1 df],'T',1,spm_vol('rob_tmap_0002.img'),0);
            t = u2;
            str = sprintf('Height thresh FDR-corr: t = %3.2f (%3.0f Ss, %3.0f df @ p < %3.4f, extent = %3.0f\n',t,length(EXPT.subjects),df,uo,k);
            disp(str)
        end
    end

    [P2,P,sigmat] = threshold_imgs('rob_tmap_0002.img',t,[k],['pos']);
    cl = mask2clusters(P2);
    eval(['save ' f '_0002_pos cl'])
    
    fprintf(1,'\n---------------------------------------------\nCovariate neg effect\n---------------------------------------------\n')
    [P2,P,sigmat] = threshold_imgs('rob_tmap_0002.img',t,[k],['neg']);
    cl = mask2clusters(P2);
    eval(['save ' f '_0002_neg cl'])
end

if exist('rob_tmap_0003.img') == 2
    
    fprintf(1,'\n---------------------------------------------\nCovariate 2 pos effect\n---------------------------------------------\n')
    %FDR threshold
    if isstr(u)
        if strcmp(u,'FDR')
            [u2,Ps,Ts] = spm_uc_FDR(.05,[1 df],'T',1,spm_vol('rob_tmap_0003.img'),0);
            t = u2;
            str = sprintf('Height thresh FDR-corr: t = %3.2f (%3.0f Ss, %3.0f df @ p < %3.4f, extent = %3.0f\n',t,length(EXPT.subjects),df,uo,k);
            disp(str)            
        end
    end
    
    [P2,P,sigmat] = threshold_imgs('rob_tmap_0003.img',t,[k],['pos']);
    cl = mask2clusters(P2);
    eval(['save ' f '_0003_pos cl'])
    
    fprintf(1,'\n---------------------------------------------\nCovariate 2 neg effect\n---------------------------------------------\n')
    [P2,P,sigmat] = threshold_imgs('rob_tmap_0003.img',t,[k],['neg']);
    cl = mask2clusters(P2);
    eval(['save ' f '_0003_neg cl'])
end
    
if exist('rob_tmap_0004.img') == 2
    
    fprintf(1,'\n---------------------------------------------\nCovariate 3 pos effect\n---------------------------------------------\n')
    %FDR threshold
    if isstr(u)
        if strcmp(u,'FDR')
            [u2,Ps,Ts] = spm_uc_FDR(.05,[1 df],'T',1,spm_vol('rob_tmap_0004.img'),0);
            t = u2;
            str = sprintf('Height thresh FDR-corr: t = %3.2f (%3.0f Ss, %3.0f df @ p < %3.4f, extent = %3.0f\n',t,length(EXPT.subjects),df,uo,k);
            disp(str)            
        end
    end
    
    [P2,P,sigmat] = threshold_imgs('rob_tmap_0004.img',t,[k],['pos']);
    cl = mask2clusters(P2);
    eval(['save ' f '_0004_pos cl'])
    
    fprintf(1,'\n---------------------------------------------\nCovariate 3 neg effect\n---------------------------------------------\n')
    [P2,P,sigmat] = threshold_imgs('rob_tmap_0004.img',t,[k],['neg']);
    cl = mask2clusters(P2);
    eval(['save ' f '_0004_neg cl'])
end

    
if exist('rob_tmap_0005.img') == 2
    
    fprintf(1,'\n---------------------------------------------\nCovariate 4 pos effect\n---------------------------------------------\n')
    %FDR threshold
    if isstr(u)
        if strcmp(u,'FDR')
            [u2,Ps,Ts] = spm_uc_FDR(.05,[1 df],'T',1,spm_vol('rob_tmap_0005.img'),0);
            t = u2;
            str = sprintf('Height thresh FDR-corr: t = %3.2f (%3.0f Ss, %3.0f df @ p < %3.4f, extent = %3.0f\n',t,length(EXPT.subjects),df,uo,k);
            disp(str)            
        end
    end
    
    [P2,P,sigmat] = threshold_imgs('rob_tmap_0005.img',t,[k],['pos']);
    cl = mask2clusters(P2);
    eval(['save ' f '_0005_pos cl'])
    
    fprintf(1,'\n---------------------------------------------\nCovariate 4 neg effect\n---------------------------------------------\n')
    [P2,P,sigmat] = threshold_imgs('rob_tmap_0005.img',t,[k],['neg']);
    cl = mask2clusters(P2);
    eval(['save ' f '_0005_neg cl'])
end




fprintf(1,'\n---------------------------------------------\nDisplaying Intercept\n---------------------------------------------\n')
cl = multi_threshold('rob_tmap_0001.img','T',df,ovl);
%save cl_intercept cl

input('Press return to continue')

if exist('rob_tmap_0002.img') == 2
    fprintf(1,'\n---------------------------------------------\nDisplaying Cov 1\n---------------------------------------------\n')
    cl = multi_threshold('rob_tmap_0002.img','T',df,ovl);
    input('Press return to continue')
    %save cl_0002 cl
end

if exist('rob_tmap_0003.img') == 2
    fprintf(1,'\n---------------------------------------------\nDisplaying Cov 2\n---------------------------------------------\n')
    cl = multi_threshold('rob_tmap_0003.img','T',df,ovl);
    input('Press return to continue')
end

if exist('rob_tmap_0004.img') == 2
    fprintf(1,'\n---------------------------------------------\nDisplaying Cov 3\n---------------------------------------------\n')
    cl = multi_threshold('rob_tmap_0004.img','T',df,ovl);
    input('Press return to continue')
end

if exist('rob_tmap_0005.img') == 2
    fprintf(1,'\n---------------------------------------------\nDisplaying Cov 4\n---------------------------------------------\n')
    cl = multi_threshold('rob_tmap_0005.img','T',df,ovl);
    input('Press return to continue')
end

warning on

return
