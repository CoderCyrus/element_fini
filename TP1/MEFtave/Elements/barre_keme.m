function [Ke,Me] = barre_keme(iel)
% Calcul de la matrice raideur Ke et de la matrice masse 
% pour un element (iel) d'une structure treillis
%
% appel [Ke,Me] = barre_keme(iel)
%    ou [Ke,Me] = feval('barre_keme',iel)
% en entree iel : numero de l'element
% en sortie Ke  : matrice raideur elementaire (2*ndim,2*ndim)
%           Me  : matrice masse   elementaire (2*ndim,2*ndim)
%
%  H.Oudin 
global ndim
global Coord Connec Nprop Prop 

ES=Prop(Nprop(iel),1);
roS=Prop(Nprop(iel),2);
X  = Coord(Connec(iel,:),:);    % ----- Coordonnees des 2 noeuds de l'element
dX = X(2,:) - X(1,:);           % ----- x2-x1  et y2-y1
if ndim == 1
    L = abs(dX);
    Ke = (ES/L)*[  1  -1
                  -1   1 ];  
    Me = (roS*L/6)*[2  1
                    1  2 ];
elseif ndim == 2
    L = sqrt(dX(1)^2 + dX(2)^2);
    c = dX(1)/L; s = dX(2)/L;   % ----- Cosinus directeurs de l'element
    cc = c*c; ss = s*s; cs = c*s;
    Ke = (ES/L)*[  cc  cs -cc -cs
                   cs  ss -cs -ss
                  -cc -cs  cc  cs
                  -cs -ss  cs  ss];
    M = (roS*L/6)*[2  1
                   1  2 ];
    P = [  c  s  0  0  
            0  0  c  s ];
    Me =P' * M * P;            
 elseif ndim == 3
     disp('================================================ ');
     disp('       element non programme ');
     disp('================================================ ');
 end
return