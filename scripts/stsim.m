function [res,W,J,X1,Y] = stsim(x,model)
%STSIM simulates steady-state labeling distributions of EMUs using the free
%fluxes contained in X and the model structure contained within MODEL. the
%output RES is the difference between the predicted labeling distributions
%and the experimentally measured values contained within MODEL. J is the
%sensitivity of the residuals contained in RES to free fluxes in X.

%initialize sizes of variables
nu = model.vardata.nu;
nh = sum([model.data.nh]);
%np = length(x);
nemu = model.mdvsim.nemu;
nip = model.mdvsim.nip;
nb = length(nemu);
u = x(1:nu);
h = x(nu+1:nu+nh);
nh = [model.data.nh];
simsens = nargout>2;
nxp = length(model.data);

% computing matrices A and B (Refer to Gopalakrishnan and Maranas, Metab
% Eng, 2018 for notations).

A = model.mdvsim.A;
B = model.mdvsim.B;
dAdu = model.mdvsim.dAdu;
dBdu = model.mdvsim.dBdu;
Y = model.mdvsim.Y;
X = model.mdvsim.X;
cY = model.mdvsim.cY;
if simsens
    sA = model.mdvsim.sA;
    sB = model.mdvsim.sB;
    sX = model.mdvsim.sX;
    sY = model.mdvsim.sY;
end
%simulating labeling distributions

for i = 1:nb
    nX = nemu(i);
    nY = nip(i);
    A{i} = A{i} + reshape(dAdu{i}*u,nX,nX);
    B{i} = B{i} + reshape(dBdu{i}*u,nX,nY);
    [Q,R] = qr(A{i});
    for xpt = 1:nxp
        X{i,xpt} = -R\(Q'*B{i})*Y{i,xpt};
        si = length(X{i,xpt}(1,:));
        if simsens
            D = sB{i}*Y{i,xpt} + sA{i}*X{i,xpt};
            D = reshape(D,nX,(si)*nu);
            sX{i,xpt} = -R\Q'*(D + B{i}*sY{i,xpt});
        end
        if i<nb
            if simsens
                [Y{i+1,xpt},sY{i+1,xpt}] = emuconv(X(:,xpt),i+1,cY{i+1},Y{i+1,xpt},sX(:,xpt),sY{i+1,xpt});
            else
                Y{i+1,xpt} = emuconv(X(:,xpt),i+1,cY{i+1},Y{i+1,xpt});
            end
        end
    end
end
X1 = X;

% organizing data for easy collection:
for xpt = 1:nxp
    for i = 1:nb
        nX = nemu(i);
        nmdv = size(X{i,xpt},2);
        X{i,xpt} = mat2cell(X{i,xpt},ones(1,nX),nmdv);
        X{i,xpt} = X{i,xpt}';
        %X{i} = cell2mat(X{i});
        if simsens
            sX{i,xpt} = sX{i,xpt}';
            sX{i,xpt} = reshape(sX{i,xpt},nu,nX*nmdv);
            sX{i,xpt} = mat2cell(sX{i,xpt},nu,nmdv*ones(1,nX));
            sX{i,xpt} = [sX{i,xpt}(:)]';
            %sX{i} = cell2mat(sX{i});
        end
    end
end
%X = cell2mat(X');
%if simsens
%    sX = cell2mat(sX');
%end

%collecting data and computing residuals
%Xh = model.Eh;
%hmult = Xh*h;
res = zeros(0,1);
std = res;
if simsens
    J = zeros(0,length(x));
end
ctr = 0;
for xpt = 1:nxp
    nvmeas = length(model.data(xpt).flxind);
    vpred = model.vardata.N*u;
    vpred = vpred(model.data(xpt).flxind);
    vmeas = model.data(xpt).flxval;
    if simsens
        dvpred = [model.vardata.N(model.data(xpt).flxind,:),zeros(nvmeas,sum(nh))];
    end
    verr = model.data(xpt).flxwt;
    xmeas = model.data(xpt).msval';
    xpred = cell(size(xmeas));
    dxpred = xpred;
    for i = 1:length(xmeas)
        ctr = ctr+1;
        msind = model.data(xpt).msind{i};
        xpred{i} = X{msind(2),xpt}{msind(3)};
        xpred{i} = conv(xpred{i},model.data(xpt).mscorr{i});
        xk = xpred{i}(1:msind(4));
        xpred{i} = h(ctr)*xpred{i}(1:msind(4));
        if simsens
            dxpred{i} = sX{msind(2),xpt}{msind(3)};
            dxpred{i} = conv2(dxpred{i},model.data(xpt).mscorr{i});
            dxpred{i} = dxpred{i}(:,1:msind(4));
            dh = zeros(sum(nh),msind(4));
            dh(ctr,:) = xk;
            dxpred{i} = [h(ctr)*dxpred{i};dh];
        end
        mserr = model.data(xpt).mswt';
    end


    res = [res',(vpred-vmeas)',[cell2mat(xpred)-cell2mat(xmeas)]]';
    std = [std',verr',cell2mat(mserr)]';

    if simsens
         J = [J',dvpred',cell2mat(dxpred)]';
    end
end
W = diag(1./(std.^2));
end


function [y,sy] = emuconv(X,esize,C,y,sX,sy)
% helper function to convolute EMUs of smaller sizes to serve as inputs to
% the EMU network of larger size

if nargin<5
    simsens = false;
else
    simsens = true;
end

%X = X(1:(esize-1));
if simsens
    %sX = sX(1:(esize-1));
    np = length(sX{1}(1,:))/2;
end
for i = 1:(esize-1)
    [nx,nmdv] = size(X{i});
    X{i} = mat2cell(X{i},ones(1,nx),nmdv);
    X{i} = X{i}';
    if simsens
        sX{i} = sX{i}';
        sX{i} = reshape(sX{i},np,nx*nmdv);
        sX{i} = mat2cell(sX{i},np,nmdv*ones(1,nx));
        %sX{i} = cell2mat(sX{i}(:)');
        %sX{i} = mat2cell(sX{i},1,(i+1)*ones(1,nx));
    end
end
%X = cell2mat(X');
X1 = X;
sX1 = sX;
X = cell(1,0);
sX = cell(1,0);
for i = 1:esize-1
    X = [X,X1{i}];
    if simsens
        sX = [sX,sX1{i}];
    end
    %X{i} = [];
end
nmdv = length(X1{esize}(1,:));
%convoluting EMUs
ny = size(C,1);
%y = zeros(ny,length(X1{esize}(1,:)));
%if simsens
%    sy = zeros(ny,np*(esize+1));
%end
inds = 1:ny;
inds = inds';
i1 = inds(any(C,2));
for i = 1:length(i1)
    ix = i1(i);
    in2 = find(C(ix,:));
    ncv = C(ix,in2);
    i2 = zeros(1,sum(ncv));
    cv = 0;
    for j = 1:length(in2)
        i2(cv+1:cv+ncv(j)) = in2(j);
        cv = cv+ncv(j);
    end
    y0 = 1;
    if simsens
        sy0 = zeros(np,nmdv);
    end
    for j = 1:length(i2)
        y0 = conv(y0,X{i2(j)});
        if simsens
            syk = 1;
            for k = 1:length(i2)
                if k == j
                    syk = conv2(syk,sX{i2(k)});
                else
                    syk = conv2(syk,X{i2(k)});
                end
            end
            sy0 = sy0 + syk;
        end
    end
    y{ix,:} = y0;
    if simsens
        sy(ix,:) = sy0(:)';
    end
end
y = cell2mat(y);
end

