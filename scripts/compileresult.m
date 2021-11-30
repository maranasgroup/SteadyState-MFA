function res = compileresult(xopt,model)

alpha = model.options.conf_lvl;
[r,W,~,X] = simlabdist(xopt,model);
res = struct('DOF',[],'frange',[],'fmin',[],'fluxes',[],'residuals',[]);
reinit_data = struct('u',[],'h',[]);
flx = struct('id', '','type','','val',[],'vLB',[],'vUB',[]);
fres = struct('id','','expt','','data',[],'val',[],'WRES',[],'SSRES',[]);   %flux residual
mres = struct('id','','expt','','time',inf,'data',[],'val',[],'WRES',[],'SSRES',[]);   %MDV residual
%cres = struct('id','','expt','','data',[],'val',[],'WRES',[],'SSRES',[]);   %pool size residual

% statistical stuff
dof = length(r)-length(xopt);
if dof<=0
    frange = [0,0];
else
    p = (1-alpha)/2;
    frange = [chi2inv(p,dof),chi2inv(1-p,dof)];
end
fmin = r'*W*r;

res.fmin = fmin;
res.DOF = dof;
res.frange = frange;


% reinitialization info
nu = model.vardata.nu;
u = xopt(1:nu);
h = xopt(nu+1:end);
reinit_data.u = u;
reinit_data.h = h;

res.reinit_data = reinit_data;

% reporting fluxes
N = model.vardata.N;
v = N*u;
vf = model.vardata.vfwd;
vr = model.vardata.vrev;
v(vf) = v(vf)-v(vr);
nv = length(v);
fluxes(1:nv) = deal(flx);
for i = 1:nv
    fluxes(i).id = model.vardata.flxdata(i).name;
    if vr(i)
        typ = 'exc';
    else
        typ = 'net';
    end
    fluxes(i).type = typ;
    fluxes(i).val = v(i);
    fluxes(i).vUB = 10000000;
    if vf(i)
        fluxes(i).vLB = -10000000;
    else
        fluxes(i).vLB = 0;
    end
end
res.fluxes = fluxes;

%residuals
nflx = 0;
nmdv = 0;
%nc = 0;
for i = 1:length(model.data)
nflx = nflx + length(model.data(i).flxval);
nmdv = nmdv + length(model.data(i).msval);
%nc = nc + length(model.data(i).cval);
end
sf = 0;
sm = 0;
%sc = 0;
dflx(1:nflx) = deal(fres);
dmdv(1:nmdv) = deal(mres);
%dc(1:nc) = deal(cres);
for i = 1:length(model.data)
    % lack of fit for fluxes
    for j = 1:length(model.data(i).flxval)
        sf = sf+1;
        dflx(sf).id = fluxes(model.data(i).flxind(j)).id;
        dflx(sf).expt = model.data(i).exptname;
        dflx(sf).data = model.data(i).flxval(j);
        dflx(sf).val = v(model.data(i).flxind(j));
        std = model.data(i).flxwt(j);
        dflx(sf).WRES = (dflx(sf).val-dflx(sf).data)/std;
        dflx(sf).SSRES = (dflx(sf).WRES).^2;
    end
    
    %pool sizes Lack of fit
    
    % MDV lack of fit
    for j = 1:length(model.data(i).msval)
        sm = sm+1;
        dmdv(sm).id = model.data(i).msid{j};
        dmdv(sm).expt = model.data(i).exptname;
        dmdv(sm).data = model.data(i).msval{j};
        j1 = model.data(i).msind{j}(2);
        j2 = model.data(i).msind{j}(3);
        dmdv(sm).val = h(sm)*conv(X{j1,i}(j2,:),model.data(i).mscorr{j});
        std = model.data(i).mswt{j};
        dmdv(sm).WRES = (dmdv(sm).val - dmdv(sm).data)./std;
        dmdv(sm).SSRES = dmdv(sm).WRES*dmdv(sm).WRES';
    end
end
res.residuals.flxfit = dflx;
res.residuals.mdvfit = dmdv;

end

    
    

