function [res] = flxestimate(model)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N = model.vardata.N;
nu = model.vardata.nu;
nh = sum([model.data.nh]);
A = [N;-N];
A = blkdiag(A,eye(nh));
b = model.vardata.vb;
b(:,2) = -1*b(:,2);
b = b(:);
b = [b;1e-7*ones(nh,1)];
nms = model.options.multistart;
fopt = Inf;
for i = 1:nms
    [x,actcon] = initialize(model);
    [xf,f,~,actcon] = lsqsolve(x,model,A,b,actcon);
    if f<fopt
        fopt = f;
        xopt = xf;
    end
end

model.options.dfbase = eps;
iter = 1;
fail = true;
x = xopt;
f = fopt;
ac = actcon;
while fail || iter <= 5
    [x,f,fail,ac] = lsqsolve(x,model,A,b,ac);
    if f<fopt
        xopt = x;
        fopt = f;
        actcon = ac;
    end
    iter = iter+1;
end
res = compileresult(xopt,model);
end

