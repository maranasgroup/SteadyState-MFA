function [res, foptCell, resCell] = flxestimate_proper(emod, repeat, randseed)
% Flux estimation function
% Sarat's flxestimate edited by Hoang Dinh

if nargin < 3
    rng('shuffle')
else
    rng(randseed)
end    

if nargin < 2
    repeat = 100;
end

N = emod.vardata.N;
nu = emod.vardata.nu;
nh = sum([emod.data.nh]);
A = [N;-N];
A = blkdiag(A,eye(nh));
b = emod.vardata.vb;
b(:,2) = -1*b(:,2);
b = b(:);
b = [b;1e-7*ones(nh,1)];

% Copy from initialize.m
nx = size(N,2);
A_init = [-N;N];
vb_init = emod.vardata.vb;
b_init = [-vb(:,1)-1e-7;vb(:,2)];
xb_init = vb_init(emod.vardata.vfree,:);

% Matrix containing randomize x-deviations
xdev_mat = rand(nx, repeat);

% Start iteration
for i = 1:repeat
    
    % Initialize: copy from initialize.m
    x0 = zeros(nx,1);
    x0 = x0 + 100*xdev_mat(:,i);
    
    opt = optimset('Display','off','MaxFunEvals',10000000,'TolCon',0);%,'ConstraintTolerance',-1e-7);
    x = fmincon(@(x) 1,x0,A_init,b_init,[],[],[],[],[],opt);
    if any(A_init*x > b_init)
        v = N*x;
        v = min(max(v,(vb_init(:,1)+1e-7)),(vb_init(:,2)-1e-7));
        options = optimset('display','Iter','maxiter',10000,'Algorithm','interior-point');
        [x,~,~,exitflag] = lsqlin(N,v,A_init,b_init,[],[],xb_init(:,1),xb_init(:,2),x,options);
    end
    actcon = find(A_init*x>=b);
    nh = sum([emod.data.nh]);
    x = [x;ones(nh,1)];
    
    % Initial run: step size = emod.options.dfbase (set in defopt)
    [x_init,f_init,~,actcon] = lsqsolve(x,emod,A,b,actcon);
        
    % Extended run: step size = MATLAB's limit step size (eps)
    emod.options.dfbase = eps;
    iter = 1; fail = true;
    x = xopt; f = fopt; ac = actcon;
    while fail || iter <= 5
        [x,f,fail,ac] = lsqsolve(x,emod,A,b,ac);
        if f<fopt
            xopt = x;
            fopt = f;
            actcon = ac;
        end
        iter = iter+1;
    end
    
end


res = compileresult(xopt,emod);
end

