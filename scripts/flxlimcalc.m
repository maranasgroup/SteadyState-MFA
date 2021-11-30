function [range,xbest,fbest] = flxlimcalc(x,A,b,actcon,Aeq,model,bounds,r,W,J)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

output = model.options.output_display;
if output
   % display output 
    
    
end
f0 = r'*W*r;
fbest = f0;
xbest = x;
xmax = bounds(2);
xmin = bounds(1);

% calculating upper bound
[maxx,xb,fb] = limflx(x,A,b,actcon,model,Aeq,xmax,r,W,J);
if fb<fbest
    fbest = fb;
    xbest = xb;
end
range(2) = maxx;

% calculating lower bound

[minx,xb,fb] = limflx(x,A,b,actcon,model,-Aeq,-xmin,r,W,J);
if fb<fbest
    fbest = fb;
    xbest = xb;
end
range(1) = -minx;


end


%Helper function for tracing bounds of X
function [xlim,xbest,fbest] = limflx(x,A,b,actcon,model,Aeq,xbound,r,W,J)

% some preliminary calculations
output = model.options.output_display;
dFmax = chi2inv(model.options.conf_lvl,1);
nsteps = 5; % minimum steps required to reach bound
dFstep = dFmax/nsteps;
dFacc = 1e-4;
x0 = x;
f0 = r'*W*r;
g0 = J'*W*r;
H0 = J'*W*J;
r0 = r;
J0 = J;
fbest = f0;
xbest = x0;
dxmin = 1e-7;   %minimum change in desired parameter
hmin = 1e-5;
fmax = f0+dFmax;
f = f0;
xcurr = Aeq*x;
xsc = max(x,1e-5);
df = eps*max(diag(H0).*xsc.^2);
[Q,R,Z,actcon] = initlqz(zeros(0,length(x)),A,-Aeq',actcon);
alf = 2;


while f<fmax && xcurr <= (xbound-dxmin)
    
    delx = 0;
    pdiff = 0;
    citer = 0;
    h = 1;
    nostep = false;
    f0 = r0'*W*r0;
    xlim0 = Aeq*x0;
    fail = false;
    
    %search direction of constrained step
    [dx,stepfail] = searchdir(Z,H0,-Aeq',df,xsc);
    if ~stepfail
        % step length determined by solving quadratic equation
        DM = df*diag(1./xsc.^2);
        bx = 2*dx'*g0;   %linear term coefficient
        c = dFstep;     % constant term coefficient
        a = abs(dx'*((H0+H0)/2 + DM)*dx);
        disc = (bx^2 + 4*a*c);
        hex = (-bx+sqrt(disc))/(2*a);
        dx = dx*hex;
        
        % performing step
        [h,blkcon] = steplength(x,dx,actcon,A,b);
        x = x0+h*dx;
        
        %predictor calculations
        rp = r0 + J0*h*dx;
        fp = rp'*W*rp;
        delx = Aeq*(x-x0);
        xcurr = xcurr + delx;
        
        % corrector calculations
        if h >= hmin || xcurr >= (xbound - dxmin)
            [r,W,J] = simlabdist(x,model);
            fact = r'*W*r;
            pdiff = abs(fp-fact);
            
            % step acceptance criteria
            if (pdiff <= dFstep && delx > 0) || (pdiff > dFstep && delx < dxmin)
                
                % acceptable step
                % update constraints
                if h < 1
                    [Q,R,Z,actcon] = include_constraints(Q,R,actcon,A,blkcon);
                else
                    [Q,R,Z,actcon] = binding_constraints(Q,R,actcon,-Aeq');
                end
                
                % apply corrections
                if pdiff < dFacc
                    fx = fact;
                else
                    act1 = find(A*x<=b);
                    model.options.output_display = false;
                    [x,fx,fail,act1,r,J,citer] = lsqsolve(x,model,A,b,act1,Aeq);
                    model.options.output_display = true;
                end
                
            else
                % unacceptable step
                fail = true;
            end
            
        else
            % too small step size
            [Q,R,Z,actcon] = include_constraints(Q,R,actcon,A,blkcon);
            nostep = true;
        end
    end
    
    % Take appropriate action
    if fail
        
        % increase damping factor
        df = df*alf;
        alf = alf*2;
    
    elseif nostep
        % no action
    
    else
        f = fx;
        
        if f < fbest
            opt = optimset('Display','off','MaxFunEvals',10000000,'TolCon',0);%,'ConstraintTolerance',-1e-7);
            x1 = fmincon(@(x) 1,x,-A,-b,[],[],[],[],[],opt);
            [rx,W,Jx] = simlabdist(x1,model);
            fz = rx'*W*rx;
            if fz < fbest
                fbest = fz;
                xbest = x1;
            end
        end
        
        % display output
        if output
            % set up output display
        end
        
        alf = 2;
        if h == 1
            Rz = dFstep/pdiff;
            df = df*max(1-(2*Rz-1)^3,1/3);
        end
        
        x0 = x;
        r0 = r;
        J0 = J;
        xsc = max(x0,1e-5);
        g0 = J0'*W*r0;
        H0 = J0'*W*J0;
        
        if citer>0
            actcon = find(A*x0 <= b);
            [Q,R,Z,actcon] = initlqz(zeros(0,length(x)),A,-Aeq',actcon);
        end
    end
        
end

f = r0'*W*r0;

%linear interpolation to compute xlim
if xcurr <= (xbound - dxmin)
    xlim = xlim0 + (((xcurr-xlim0)/(f-f0))*(fmax-f0));
else
    xlim = xcurr;
end



end


function [Q,R,Z,actcon] = initlqz(Aeq,A,g,actcon)
% INITQRZ performs a QR decomposition

[neq,nx] = size(Aeq);
nactineq = length(actcon);

% Handling dependent equality constraints
% Equality constraints will appear during confidence interval determination
% in which case a double projection must be performed. The first projection
% is on the space of variables remaining after elimintation of equality
% constraint. The second projection is on the set of active and *binding*
% constraints. Binding constraints will be identified using their
% corresponding Lagrange Multipliers

if neq == 0 % For the simple parameter estimation problem
    Z = eye(nx,nx);
    Q = Z;
    R = zeros(nx,0);

else        % For the parameter sensitivity and confidence interval problem
    r = rank(Aeq);
    if r < neq
        [~,~,P] = qr(Aeq'); % The MATLAB function QR returns an optional output P, a matrix of zeros and ones indicating dependent constraints
        i = any(P(:,1:r),2);% Identify dependent constraints
        Aeq = Aeq(i,:);
    end
    [Q,R] = qr(Aeq');
    Z = Q(:,size(R,2)+1:nx);% Matrix Q2 (see function description)
end
% Handling the active equalities within inequality constraints
if nactineq > 0 %some parameters are hitting their bounds
    % First projection to account for supplied equality constraints
    actA_p1 = A(actcon,:)*Z;
    %tol = eps;
    
    % Handling Rank-deficiency within active inequalities
    r = rank(actA_p1);
    if r < nactineq
        [~,~,P] = qr(actA_p1');
        i = any(P(:,1:r),2);
        %actA_p1 = actA_p1(i,:);
        actcon = actcon(i);
        nactineq = sum(i);
    end
    eq_con = [A(actcon,:),Aeq]; %all active equality constraints
    [Q,R] = qr(eq_con');
    [Q,R,Z,actcon] = binding_constraints(Q,R,actcon,g);
    % Eliminating non-binding constraints
    nactineq = nactineq + 1;
    while nactineq > length(actcon)
        nactineq = length(actcon);
        [Q,R,Z,actcon] = binding_constraints(Q,R,actcon,g);
    end
end
Z = Q(:,size(R,2)+1:nx);
end

function [dx,stepfail] = searchdir(Z,H,g,df,x)
% SEARCHDIR computes the Newtonian step direction by solving the linear
% system of equations: H*dx = -g. Appropriate null-space projection must be
% performed to account for active inequality constraints. The Hessian is
% scaled using a damping factor related to DF. This contributed to a
% steepest-descent-like step during the earlier stages of the optimization
% procedure and a Newtonian step closer to the optimal solution.
% Singularity of the Hessian is also accounted for using appropriate
% modifications on the Cholesky decomposition of H

% Damping the Hessian
DM = df*diag(1./x.^2);
H = (H+H')/2 + DM;

% derivative projection
H = Z'*H*Z;
G = Z'*g;

% Check whether Hessian is positive definite
[R,p] = chol(H);
warning('off','all');
% compute search direction
if p > 0
    dx = -Z*real(H\G);
    if g'*dx > 0
        dx = -dx;
    end
else
    dx = -Z*(R\(R'\G));
end

if ~any(dx)
    stepfail = true;
else
    stepfail = false;
end
end

function [h,blkcon] = steplength(x0,dx,actcon,A,b)
% STEPLENGTH identifies the maximum step length H that can be taken along
% the search direction DX before a constraint from A*X >= B is violated.
% The indices of binding constraints are contained within ACTCON and A*DX
% for these constraints is always zero. See Documentation for mathematical
% details

% mask of binding constraints
nc = size(A,1);
cons = false(nc,1);
cons(actcon) = true;

% Setting a tolerance to prevent division by zero
db = A*dx;
tol = (1e-10)*norm(db);

% Finding encountered constraints and corresponding allowed step lengths
blkcon = find(db < -tol & ~cons);
if isempty(blkcon)
    h = 1;
    blkcon = [];
else
    db = db(blkcon);
    slack = min(0,(b(blkcon)-(A(blkcon,:)*x0)));
    [h,ndx] = min(slack./db);
    if h<=1
        h = min(h,1);
        blkcon = blkcon(ndx);
    else
        h = 1;
        blkcon = [];
    end
end
end

function [Q,R,Z,actcon] = include_constraints(Q,R,actcon,A,blkcon)
% INCLUDE_CONSTRAINTS appends the newly encountered constraints contained
% within BLKCON to the existing active constraints ACTCON and updates the
% QR decomposition and the intermediate null-space projection matrix Z. The
% required rows are extracted from A.

% Where to insert new constraints (See documentation of QRINSERT)
j = length(actcon)+1;
nx = length(Q);

% Updating matrices
if any(blkcon)
    [Q1,R1] = qrinsert(Q,R,j,A(blkcon,:)');
    if rank(R1) == size(R1,2)   %full rank, no dependent constraints
        actcon = [actcon;blkcon];
        Q = Q1;
        R = R1;
    else                        % dependent constraints introduced
        actcon = [actcon;blkcon];
        Aeq = A(actcon,:);
        r = rank(Aeq');
        [~,~,P] = qr(Aeq');
        i = any(P(:,1:r),2);
        actcon = actcon(i);
        Aeq = A(actcon,:);
        [Q,R] = qr(Aeq');
    end
end
Z = Q(:,size(R,2)+1:nx);

end


function [Q,R,Z,actcon] = binding_constraints(Q,R,actcon,g)
% BINDING_CONSTRAINTS updates the decomposition matrices Q and R from the
% QR decomposition and generates an apropriate projection matrix Z to
% project the solution space onto the null-space of binding (in)equality
% constraints. Binding constraints have a negative Lagrange Multiplier
% which can be computed by solving the equation dL(x,lambda)/dx = 0 where
% L(x, lambda) is the Lagrangian of the objective function

% Storing dimensions for ease
nx = size(Q,1);
nact = length(actcon);

if nact > 0
    
    % computing Lagrange Multipliers
    lambda = R\(Q'*g);
    lambda = lambda(1:nact);
    lneg = lambda < 0;
    if any(lneg)
        nbcind = find(actcon == min(actcon(lneg)));
        [Q,R] = qrdelete(Q,R,nbcind);
        actcon(nbcind) = [];
    end
end
Z = Q(:,size(R,2)+1:nx);% Projection matrix
end