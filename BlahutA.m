function [obj,Dtx,Errs] = BlahutA(p,Nu,la)
% setup the parametr
[ms,ns] = size(Nu); 
r = ones(ns,1)/ns;
ones_p = ones(ms,1);
ones_r = ones(ns,1);
W = ones(ms,ns);
iter = 200;
tol_BA = 1e-15;

% the main iteration
for i = 1:iter
    % update r
    r = W'*p;
    % updata w
    W = exp(-la*Nu)*diag(r);
    regs = W*ones_r;
    pp = ones_p./regs;
    W = diag(pp)*W;
    % save the previous data
    % the error 
    Err(1) = norm(sum(r)-1);
    Err(2) = norm(W'*p-r);
    Err(3) = norm(W*ones_r-ones_p);
    Errs = max(Err);
    if Errs < tol_BA
       break;
    end
%     if isnan(Err_sink(i,2)) || isinf(Err_sink(i,2))
%         break;
%     end
end
% W/r
Wr = W*diag(ones_r./r);
% the output
mat_obj = diag(p)*(W.*log(Wr));
obj = ones_p'*mat_obj*ones_r;
mat_Dtx = diag(p)*(W.*Nu);
Dtx = ones_p'*mat_Dtx*ones_r;
end