function [ xopt,fopt,fail,actcon,ropt,drdxopt,iter ] = lsqsolve( x,model,A,b,actcon,Aeq )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%initialize technical parameters for least-squares solver
iter = 0;       %iteration counter
%h = 1;          %default step length. Will be modified in each iteration
done = false;   %optimization termination indicator
fail = true;    %optimization success indicator
first = true;   %the first iteration initializes QR decomposition
                %every subsequent iteration involves optimality checks

% Levenberg-Marquardt algorithm parameters
nx = length(x);
tolx = 0.01; % Termination tolerance on x
tolg = 0.1; % Termination tolerance on gradient
alf = 2;                % parameter to update damping factor
itermax = 500*nx;
dfmax = 1e50;           %maximum allowed damping of Hessian. Exceeding this value terminates the least-squares solver and reports failure
dfbase = model.options.dfbase;
hmin = 1e-5;            %Required minimum step length. h < hmin updates active constraints
dispop = model.options.output_display;
%additional optional inputs to integrate range estimation
if nargin < 6
    Aeq = zeros(0,nx);  %equality constraints to account for parameters hitting the bounds
end
if nargin < 7
    [r,W,drdx] = simlabdist(x,model); % initial function evaluation
end


%computing initial function value
f = r'*W*r;
if isnan(f)
    msg = 'Initial function evaluation failed. Please check network decomposition';
    done = true;
end

%store initial values
f0 = f;
r0 = r;
drdx0 = drdx;
x0 = x;
xsc = max(x,1e-5);

%output display initalization
if dispop
    fprintf('   Iter           SSRES          dSSRES         st-size              DF\n')
    fprintf('%7.0f    %12.4g\n',iter,f0);
    fmt = ('%7.0f    %12.4g    %12.4g    %12.4f    %12.4g\n');
end

%main optimization loop
%Note that the minimum allowed function value is 0 due to squared
%deviations
while ~done && f0 > 0
    
    %derivatives
    g = drdx'*W*r; % Approximate gradient
    H = drdx'*W*drdx; % Approximate Hessian
    
    % Before computing next step, generate LQ decompositions if first
    % iteration or check for termination when iter > 1. iter is not used as
    % a conditional variable here as it is only updated in every SUCCESSFUL
    % iteration. A failed first iteration will cause the program to loop
    % around the first iteration and repeatedly perform LQ decomposition
    
    if first
        
        %store initial derivative information
        g0 = g;
        H0 = H;
        
        %initial damping factor
        %df = max(0.001*max(diag(H).*xsc.^2),1e6);
        df = dfbase*max(diag(H).*xsc.^2);
        % Initial search direction
        [Q,R,Z,actcon] = initlqz(Aeq,A,g0,actcon);  %initialize LQ decomposition and obtain the projection matrix Z to handle active constraints
        [dx,stepfail] = searchdir(Z,H,g0,df,xsc);  %compute feasible search direction. Stepfail occurs when dx is 0 indicating that the initial point is optimal
        first = false;
        
    else
        
        %display output
        if dispop
            if h > hmin
                dfx = g'*(h*dx);
                fprintf(fmt,iter,f,dfx,h,df)
            end
        end
        
        % Check step acceptance
        if f < f0   % accept step
            
            % update damping factor (Madsen 2002)
            alf = 2; % reset update factor
            if h == 1
                Rz = (f0-f)/max(f0-f1,eps);
                df = df*max(1-(2*Rz-1)^3,1/3);
            end
            
            %updating values
            f0 = f;
            r0 = r;
            drdx0 = drdx;
            g0 = g;
            H0 = H;
            x0 = x;
            xsc = max(x0,1e-5);
            
            %updating active constraints
            if h < 1    % New active constraints to be included
                [Q,R,Z,actcon] = include_constraints(Q,R,actcon,A,blkcon);
            else
                [Q,R,Z,actcon] = binding_constraints(Q,R,actcon,g0);
            end
            
        else    % Step not accepted: Update active constraints if step size is too small or reduce stesize by increasing damping
            
            if h < hmin     % too small step
                [Q,R,Z,actcon] = include_constraints(Q,R,actcon,A,blkcon);
            else            % reduce steo by updating damping factor
                df = df*alf;
                alf = alf*2;
            end
        end
            
        % Compute search direction
        [dx,stepfail] = searchdir(Z,H0,g0,df,xsc);

        % Termination Criteria (Madsen 2002)
        % 1. No feasible search direction and (a) no active constraints
        % or (b) with active constraints: Optimum found
        % 2. Converged to non-optimal point due to :
            % a. too small step size caused by excessive damping df > dfmax
            % b. Maximum iterations exceeded
        TermX = max(abs(dx./x0)) < tolx;
        TermG = max(abs(g0.*x0)) < tolg;
        if ~stepfail && isempty(blkcon) && (TermX || TermG)
            done = true;
            if dispop
            if any(actcon)
                disp('Optimal solution found with active constraints');
            else
                disp('Unconstrained optimal solution found');
            end
            end
            fail = false;
        elseif df > dfmax
            done = true;
            disp('Converged to non-optimal point due to excessive Hessian Damping');
        elseif iter > itermax
            done = true;
            disp('Maximum iterations exceeded');
        end
        
    end     % Ending differences between first and subsequent iterations
    
    % End if done
    if done
        break
    end
    
    % Performing step
    [h,blkcon] = steplength(x0,dx,actcon,A,b);
    
    if h > hmin && ~stepfail
        x = x0 + h*dx;
        [r,~,drdx] = simlabdist(x,model);
        f = r'*W*r;
        iter = iter+1;
        r1 = r0 + drdx0*dx*h;
        f1 = r1'*W*r1;
    else
        r = r0;
        f = f0;
        drdx = drdx0;
    end
    
end

%outputs

%disp(char('',msg)); %termination message
xopt = x0;
fopt = f0;
ropt = r0;
drdxopt = drdx0;

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
    eq_con = [A(actcon,:);Aeq]; %all active equality constraints
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
    slack = (b(blkcon)-(A(blkcon,:)*x0));
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


    
    
    
    
    
    
    