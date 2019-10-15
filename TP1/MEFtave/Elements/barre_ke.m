function [Ke,Fe] = barre_ke(iel)
% Calcul de la matrice raideur Ke et de la force generalisee Fe 
% pour un element (iel) d'une structure treillis
%
% appel [Ke,Fe] = barre_ke(iel)
%    ou [Ke,Fe] = feval('barre_ke',iel)
% en entree iel : numero de l'element
% en sortie Ke  : matrice raideur elementaire (2*ndim,2*ndim)
%           Fe  : force generalisee elementaire (2*ndim) 
%
%  H.Oudin 
global ndim
global Coord Connec Nprop Prop 

ES=Prop(Nprop(iel),1);
X  = Coord(Connec(iel,:),:);    % ----- Coordonnees des 2 noeuds de l'element
dX = X(2,:) - X(1,:);           % ----- x2-x1  et y2-y1
if ndim == 1
    L = abs(dX);
    Ke = (ES/L)*[  1  -1
                  -1   1 ];
    f=Prop(Nprop(iel),2);
    Fe = (f*L/2)*[1;1];
elseif ndim == 2
    L = sqrt(dX(1)^2 + dX(2)^2);
    c = dX(1)/L; s = dX(2)/L;   % ----- Cosinus directeurs de l'element
    cc = c*c; ss = s*s; cs = c*s;
    Ke = (ES/L)*[  cc  cs -cc -cs
                   cs  ss -cs -ss
                  -cc -cs  cc  cs
                  -cs -ss  cs  ss];
     fx=Prop(Nprop(iel),2); fy=Prop(Nprop(iel),3);
     Fe = (L/2)*[fx;fy;fx;fy];
     %Fe =(L/2)*[fx*c-fy*s; fx*s+fy*c; fx*c-fy*s; fx*s+fy*c];
 elseif ndim == 3
     disp('================================================ ');
     disp('       element non programme ');
     disp('================================================ ');
 end
return