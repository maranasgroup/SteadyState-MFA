function [res, foptCell] = flxestimate_fast(emod, randscale)
% Flux estimation function
% Sarat's flxestimate
% First, run non-linear SSR minimization with step size 1e-6
%   with randomized start, run for number of times equal
%   to the user's settings (multistart in defopt.m)
% Second, best solutions among multistart runs are selected
%   for subsequent run with step size of MATLAB's limit (eps)
%   to check for convergence

% Setting scale to variation using random function
% If do not set scaling for randomization, set scaling to
% 100 mmolglucose/gDW/h
if nargin < 2
    randscale = 100;
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
nms = emod.options.multistart;
fopt = Inf;
foptCell = zeros(nms,1);
for i = 1:nms
    [x,actcon] = initialize(emod,randscale);
    [xf,f,~,actcon] = lsqsolve(x,emod,A,b,actcon);
    foptCell(i) = f;
    if f<fopt
        fopt = f;
        xopt = xf;
    end
end

emod.options.dfbase = eps; % Use smaller step size here
iter = 1;
fail = true;
x = xopt;
f = fopt;
ac = actcon;
while fail || iter <= 5
    [x,f,fail,ac] = lsqsolve(x,emod,A,b,ac);
    if f<fopt
        xopt = x;
        fopt = f;
        actcon = ac;
    end
    iter = iter+1;
end
res = compileresult(xopt,emod);
end

