function [obj,Wr,a,b,D_tx,Err_sink,Err] = sinkhorn_admm(p,Nu,Dtx)
% perform the sinkhorn 
niter_sinkhorn = 3000;
tol_newton = 1e-15;
tol_did = 1e-15;
tol_sink = 1e-15;

% function G and G'
Ga = @(ax,bx,lax,rx) (ax.*p)'*(Nu.*exp(-lax*Nu))*(bx.*rx)-Dtx;
Ga_gd = @(ax,bx,lax,rx) -(ax.*p)'*((Nu.*Nu).*exp(-lax*Nu))*(bx.*rx);

% set up for sinkhorn
[ms,ns] = size(Nu); r = ones(ns,1)/ns;
ones_p = ones(ms,1); ones_r = ones(ns,1);
a = ones(ms,1); % a is \phi
b = ones(ns,1); % b is \psi
la = 0; Err_sink = [];

% function F and F'
Fa = @(Wrx,betx,etx) ones_r'*((Wrx'*p)./(etx*ones_r-betx))-1;
% Fa_gd = @(Wrx,betx,etx) -(ones_r'*((Wrx'*p)./((etx*ones_r-betx).^2)));

% numerical settings
Ker = exp(-la*Nu);
Ker(Ker<1e-30)=1e-30; % Safe

% the main iteration
for i = 1:niter_sinkhorn
    
    % perform the sinkhorn
    for ij = 1:1
        b = ones_r ./ (Ker'*(a.*p)); % important for numerical settings
        a = ones_p ./ (Ker*(b.*r));
        % b = ones_r ./ (Ker'*(a.*p));
    end
    
    % compute lambda by newton
    if Ga(a,b,0,r) >= 0
        la = newton_la(a,b,la,r,tol_newton,Ga,Ga_gd);
    else
        la = 0;
    end
    
    % numerical settings
    Ker = exp(-la*Nu);
    Ker(Ker<1e-30)=1e-30; % Safe 
    % 可能取的太小了！！！ 
    
    % set Wr and beta
    Wr = diag(a)*Ker*diag(b.*r);
    beta = -log(b)-0.5*ones_r;
    % compute eta by newton
    [eta,psn] = newton_et(Wr,beta,tol_did,Fa);
    
    % update variable r first
    r = ((b.*r).*(Ker'*(a.*p)))./(eta*ones_r-beta);
    regs = (a.*p)'*Ker*(b.*r);
    
    % compute errors after all perfomance; the KKT
    Err_sink(i,1) = norm(b.*(Ker'*(a.*p))-ones_r);
    Err_sink(i,2) = norm(a.*(Ker*(b.*r))-ones_p); 
    Err_sink(i,3) = norm(la*Ga(a,b,la,r)/(la+1));
    Err_sink(i,4) = norm(regs-1);
    Err_sink(i,5) = norm(Ga(a,b,la,r));
    Err_sink(i,6) = eta;
    Err_sink(i,7) = Fa(Wr,beta,eta);
    Err_sink(i,8) = psn;
    if max(Err_sink(i,1:4)) < tol_sink
        break;
    end
    % for NAN data
    if isnan(Err_sink(i,2)) || isinf(Err_sink(i,2))
        break;
    end
end

% we must make this again verse the Safe
Ker = exp(-la*Nu);

% output the result
W = diag(a)*Ker*diag(b);
Wr = diag(a)*Ker*diag(b.*r); % the output 
Wpr = diag(a.*p)*Ker*diag(b.*r);
obj = ones_p'*(log(W).*Wpr)*ones_r;
D_tx = ones_p'*(Wpr.*Nu)*ones_r;

Errs(1) = norm(sum(r)-1);
Errs(2) = norm(Wr'*p-r);
Errs(3) = norm(Wr*ones_r-ones_p);
Err = max(Errs);
end

% Newton method for lambda
function [lax] = newton_la(ax,bx,lax,rx,tol,Ga,Ga_gd)
for i = 1:15
    lax = lax - Ga(ax,bx,lax,rx)/(Ga_gd(ax,bx,lax,rx));
    if norm(Ga(ax,bx,lax,rx)) < tol
        break;
    end
end
end

% search the root by divide the interval
function [etx,p] = newton_et(Wrx,betx,tol,Fa)
a = max(betx)+1e-30; 
b = max(betx)+length(betx); 
p = -1;
etx = a;
% p<12000 is good parameter;
while((Fa(Wrx,betx,a)*Fa(Wrx,betx,b)<=0 && p<12000) && (abs(a-b)>tol || norm(Fa(Wrx,betx,etx))>tol))
    etx = (a+b)/2;
    if Fa(Wrx,betx,etx)*Fa(Wrx,betx,b) <= 0
        a=etx;
        p = p+1;
    else
        b=etx;
        p = p+1;
    end
end
end