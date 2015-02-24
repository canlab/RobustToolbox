function robust_correlation_working

% robust correlation
rr = robustfit(scale(x),scale(y)); rr = rr(2);

nullrr = zeros(1000,1);
y = scale(y);

warning off
for i = 1:1000
    xi = getRandom(y);
    rr = robustfit(xi,y);
    nullrr(i) = rr(2);
    
end

warning on

