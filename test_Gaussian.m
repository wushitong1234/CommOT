clc;clear

%% set conition
% set conition
delta = 1;
int_x = 10; int_y = 10;
X = -int_x:delta:int_x;
mu = 0; sigma = 1; % when sigam = 1 is good
for i = 1:length(X)
    Px(i,1) = normcdf(X(i)+delta/2,mu,sigma)-normcdf(X(i)-delta/2,mu,sigma);
end
Px_sum = sum(Px);
Px = Px/Px_sum;
Y = -int_y:delta:int_y; Y = Y';
m = length(X); n = length(Y);
% set distance matrix
temp_x = repmat(X,n,1);
temp_y = repmat(Y,1,m);
temp_D = (temp_x-temp_y)';
temp_D = temp_D.^2;
Nu = temp_D;


%% run cvx and sinkhorn
upper_Dtx = 5;

% % the sinkhorn time
ans_time_sink = tic;
% set the parameter Dx for the rate distortion constriants
Dtx = 0:0.1:upper_Dtx;
for i = 1:length(Dtx)
    [obj_sink(i),w{i},phi{i},psi{i},Dtx_sink(i),Err_sink{i},Err_sinkr{i}] = sinkhorn_admm(Px,Nu,Dtx(i));
end
ans_time_sink = toc(ans_time_sink);

% the Blahut-Arimoto Algorithm
ans_time_BA = tic;
lam = 0:0.01:10;
for i = 1:length(lam)
    [obj_BA(i),Dtx_BA(i),Err_BA(i)] = BlahutA(Px,Nu,lam(i));
end
ans_time_BA = toc(ans_time_BA);
Dtx_BA(Dtx_BA>upper_Dtx)=upper_Dtx;


% the exact result
for i = 1:length(Dtx)
    obj_Rate(i) = Rate_Gaussian(Dtx(i),sigma);
end

% plot figure
plot(Dtx,obj_sink,'r*');
% plot(D_sink,obj_sink,'r.');
hold on
% plot(D_sink,obj_sink,'g--');
plot(Dtx_BA,obj_BA,'g.');
hold on
% plot(Dtx_BA,obj_BA,'b.');
plot(Dtx,obj_Rate,'b-');


