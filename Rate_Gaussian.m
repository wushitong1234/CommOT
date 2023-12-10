% the rate distortion function R(D)
function [b_rate] = Rate_Gaussian(Dx,sigma)
if Dx <= sigma^2
    b_rate = 0.5*log(sigma^2/Dx);
else
    b_rate = 0;
end
end
