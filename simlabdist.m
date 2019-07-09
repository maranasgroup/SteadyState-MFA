function [ r,W,J,X ] = simlabdist( x,model )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if model.options.ss
    [r,W,J,X] = stsim(x,model);
else
    [X,dXdp] = instsim(x,prm,opt);
end

%[r,J] = measextract(X,dXdp.data);

end