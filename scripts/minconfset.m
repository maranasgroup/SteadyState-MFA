function [vconf] = minconfset(emod)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

vrev = emod.vardata.vrev;
vfwd = emod.vardata.vfwd;
%rangecompute = false(length(vrev),1);
nfix = false(length(vrev),1);
for i = 1:length(emod.data)
    nfix(emod.data(i).flxind) = true;
end


% EMU cut

dxdv = cell2mat([emod.mdvsim.dAdv;emod.mdvsim.dBdv]);
dx1dv = dxdv;
dx1dv(:,vrev) = 0;
k = any(dx1dv,2);
dx1dv = dx1dv(k,:);
dx1dv = abs(dx1dv);
dx1dv = spones(dx1dv);
dx1dv = unique(dx1dv,'rows');
nv = any(dx1dv,1)';


% Coupled flux cut
params.outputflag = 0;
params.DualReductions = 0;
%srcv = false(size(nv));
Sf = emod.vardata.S;
Sf1 = spones(abs(Sf));
Sf1 = sum(Sf1,2);
k = Sf1<=1;
Sf = Sf(k,:);
bdry_v = any(Sf,1)';
Sf(Sf>0) = 0;
srcv = any(Sf,1);
%nt = sum(srcv);
N = emod.vardata.N;
N1 = N;
N1(vfwd,:) = N1(vfwd,:)-N1(vrev,:);
%N1(vrev,:) = [];

%setting up gurobi model for flux coupling analysis
%I = zeros(length(srcv),sum(srcv));
%I(srcv,:) = eye(sum(srcv));
%A = [N];%,-I];
A = emod.vardata.S_bal;
b = zeros(length(A(:,1)),1);
I = zeros(length(b),sum(srcv));
%I(srcv,:) = eye(sum(srcv));
A = [A,-I];
I1 = zeros(sum(srcv),length(srcv));
I1(:,srcv) = eye(sum(srcv));
I1 = [I1,-eye(sum(srcv))];
A = [A;I1];
b = zeros(length(A(:,1)),1);
sense = '';
sense(1:length(b)) = '=';
%lb = 1e-8*ones(length(N(1,:)),1);%+sum(srcv),1);
%ub = 1e12*ones(length(N(1,:)),1);%+sum(srcv),1);
lb = 1e-8*ones(length(A(1,:)),1);
ub = 1e6*ones(length(A(1,:)),1);
bdry_c = bdry_v;
bdry_c(nfix) = false;
bdry_ind = find(bdry_c);

%flx_ind = find(nv);
N1x = N1;
ncpl = false(size(nv));
N1 = eye(length(lb));%,zeros(size(I))];
N1(vfwd,:) = N1(vfwd,:)-N1(vrev,:);
%N1 = [full(N1),zeros(length(srcv),sum(srcv))];
nind = find(nfix);
for i = 1:length(nind)
    gmod.A = sparse(A);
    gmod.rhs = b;
    gmod.sense = sense;
    gmod.lb = lb;
    gmod.ub = ub;
    %gmod.A = [gmod.A;sparse(A(nind(i),:))];
    gmod.A = [gmod.A;sparse(N1(nind(i),:))];
    gmod.rhs = [gmod.rhs;1];
    gmod.sense = [gmod.sense,'='];
    gmod.vtype(1:length(lb)) = 'C';
    for j = length(bdry_ind):-1:1
        gmod.obj = N1(bdry_ind(j),:)';
        gmod.modelsense = 'max';
        r = gurobi(gmod,params);
        if ismember({r.status},{'OPTIMAL'})
            r1 = r.objval;
        else
            r1 = 1e12;
        end
        gmod.modelsense = 'min';
        r = gurobi(gmod,params);
        if ismember({r.status},{'OPTIMAL'})
            r2 = r.objval;
        else
            r2 = 0;
        end
        %r2 = r.objval;
        if abs(r1-r2)<1e-5
            bdry_ind(j) = [];
            %ncpl(flx_ind(j))=true;
        end
    end
end

bdry_v = nfix;
bdry_v(bdry_ind) = true;
bdry_ind = find(bdry_v);
nv(bdry_ind) = false;
flx_ind = find(nv);

for i = 1:length(bdry_ind)
    gmod.A = sparse(A);
    gmod.rhs = b;
    gmod.sense = sense;
    gmod.lb = lb;
    gmod.ub = ub;
    %gmod.A = [gmod.A;sparse(A(bdry_ind(i),:))];
    gmod.A = [gmod.A;sparse(N1(bdry_ind(i),:))];
    gmod.rhs = [gmod.rhs;1];
    gmod.sense = [gmod.sense,'='];
    for j = length(flx_ind):-1:1
        gmod.obj = N1(flx_ind(j),:)';
        gmod.modelsense = 'max';
        r = gurobi(gmod,params);
        if ismember({r.status},{'OPTIMAL'})
            r1 = r.objval;
        else
            r1 = 1e12;
        end
        gmod.modelsense = 'min';
        r = gurobi(gmod,params);
        if ismember({r.status},{'OPTIMAL'})
            r2 = r.objval;
        else
            r2 = 0;
        end
        if abs(r1-r2)<1e-5
            ncpl(flx_ind(j))=true;
            flx_ind(j) = [];
        end
    end
end

%nv(ncpl) = false;

% balanced metabolite cut
Sf = emod.vardata.S;
Sf(:,vrev) = 0;
C = cell2mat(emod.mdvsim.dcxdc);
nc = any(C);
Sf = Sf(nc,:);
Sf = full(spones(Sf));
Sf1 = sum(Sf,2);
k = Sf1<2;
Sf(k,:) = [];
nb = false(size(Sf(:,1)));
nch = or(nv,nfix);
for i = 1:length(nb)
    nbx = logical(Sf(i,:));
    if sum(nch(nbx))==sum(nbx)
        nb(i) = true;
    end
end
Sf1 = sum(Sf,2);
k = Sf1==2;
nb = or(nb,k);
Sf = Sf(nb,:);
Sf(:,ncpl) = 0;
k = any(Sf,2);
Sf = Sf(k,:);
inds = [];
for i = length(Sf(:,1)):-1:1
    if any(Sf(i,:))
        f = find(Sf(i,:));
        Sf(:,f(1)) = 0;
        inds = [inds,f(1)];
    end
end
ncpl(inds) = true;

    


% grouped flux cut
dx1dv(:,ncpl) = 0;
k = any(dx1dv,2);
dx1dv = dx1dv(k,:);
d1 = sum(dx1dv,2);
k = d1<=1;
X1 = dx1dv(k,:);
X2 = dx1dv(~k,:);
nf = any(X1,1);
X2(:,nf) = 0;
k = any(X2,2);
X2 = X2(k,:);
d1 = sum(X2,2);
k = d1<=1;
while any(k)
    
    nf = or(nf,any(X2(k,:),1));
    X2(:,nf) = 0;
    k = any(X2,2);
    X2 = X2(k,:);
    d1 = sum(X2,2);
    k = d1<=1;
end
nf = nf';
nf(bdry_v) = true;
nf(nfix) = true;
nf(ncpl) = false;
X2 = unique(X2,'rows');
nf_find = find(nf);
vconf = zeros(length(nf));
for i = 1:length(nf_find)
    vconf(i,nf_find(i)) = 1;
end
vconf = [vconf;X2];


end

