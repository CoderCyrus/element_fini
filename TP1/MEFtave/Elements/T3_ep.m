function [Ke,Fe] = T3_ep(iel)
% Calcul de la matrice raideur Ke et de la force generalisee Fe 
% pour un element T3 d'une structure en elasticite plane
% 
%   Nous n'utilisons pas l'integration numerique, le calcul des
%   matrices est presente dans le cours chapitre V-2.4 nous utilisons
%   ici ces resultats.
%
% appel [Ke,Fe] = T3_ke(iel)
%    ou [Ke,Fe] = feval('T3_ke',iel)
% en entree iel : numero de l'element
% en sortie Ke  : matrice raideur elementaire (6,6)
%           Fe  : force generalisee elementaire (6,1)
%
%  H.Oudin 
global Coord Connec Nprop Prop 

ndle = 6;                %----- initialisations
Ke = zeros(ndle); Fe = zeros(ndle,1); 

E=Prop(Nprop(iel),1);    %----- matrice d'elasticite D
nu=Prop(Nprop(iel),2);
ep=Prop(Nprop(iel),3);
if ep > 0  a = 0 ; else a = 1 ; ep = 1;  end
coef = ep * E * (1-a*nu)/((1+nu)*(1-nu-a*nu));
D = coef * [     1       nu/(1-a*nu)            0             ;...
            nu/(1-a*nu)      1                  0             ;...
                 0           0       .5*(1-nu-a*nu)/(1-a*nu)];

X  = Coord(Connec(iel,[1:3]),:);    % coordonnees des noeuds de l'element

                                    %-----  determinant detj = 2A
detj=(X(2,1)-X(1,1))*(X(3,2)-X(1,2))-(X(3,1)-X(1,1))*(X(2,2)-X(1,2));
        
B=zeros(3,6);                       %-----  matrice B(3,6)  Matrice B des déformations
B(1,[1 3 5])=(X([2 3 1],2)-X([3 1 2],2) )' /detj;
B(2,[2 4 6])=(X([3 1 2],1)-X([2 3 1],1) )' /detj ;
B(3,[1 3 5 2 4 6]) =[ B(2,[2 4 6]), B(1,[1 3 5]) ];
        
Ke = (B'*D*B)*.5*detj;              %----- matrice Ke(6,6)
                                    %----- vecteur Fe(6,1)
fx=Prop(Nprop(iel),4); fy=Prop(Nprop(iel),5);
Fe([1 3 5],1) = (ep*fx*detj/6)*[1 1 1]';
Fe([2 4 6],1) = (ep*fy*detj/6)*[1 1 1]';

%disp(Ke),disp(Fe)
return