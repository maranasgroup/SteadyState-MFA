function [x,actcon] = initialize(emod, randscale)

% Setting scale to variation using random function
% If do not set scaling for randomization, set scaling to
% 100 mmolglucose/gDW/h
if nargin < 2
    randscale = 100;
end

N = emod.vardata.N;
nx = size(N,2);
x0 = zeros(nx,1);
A = [-N;N];
vb = emod.vardata.vb;
b = [-vb(:,1)-1e-7;vb(:,2)];
x0 = x0 + randscale*rand(size(x0));
xb = vb(emod.vardata.vfree,:);
opt = optimset('Display','off','MaxFunEvals',10000000,'TolCon',0);%,'ConstraintTolerance',-1e-7);
x = fmincon(@(x) 1,x0,A,b,[],[],[],[],[],opt);
if any(A*x > b)
    v = N*x;
    v = min(max(v,(vb(:,1)+1e-7)),(vb(:,2)-1e-7));
    options = optimset('display','Iter','maxiter',10000,'Algorithm','interior-point');
    [x,~,~,exitflag] = lsqlin(N,v,A,b,[],[],xb(:,1),xb(:,2),x,options);
end
actcon = find(A*x>=b);
nh = sum([emod.data.nh]);
x = [x;ones(nh,1)];
