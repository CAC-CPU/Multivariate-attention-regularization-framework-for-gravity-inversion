function inversion()

% Main part of the inversion code
% Inversion result is stored in showpar for future re-loading.
% showinv(showpar): display inversion result

addpath lib;
load data.dat;
load topo.dat;
inverpar = makeinverpar();
[A d inverpar2] = GetInverpar(inverpar,data,topo);

switch inverpar.Modelspace
    case 1
        showpar = linearinv(A,d,inverpar2);
    case 2
        showpar = loginv(A,d,inverpar2);
end
save showpar showpar;
 
if inverpar.dispflag == 1
    showinv(showpar);
end
if inverpar.isdrawsclice == 1
    draw
end