function [Ne,Te,Mfe] = poutre_stress(iel,Ue)
% Calcul des contraintes dans un element poutre d'une structure portique
%
% appel poutre_stress(iel,Ue)
%    ou feval('poutre_stress',iel,Ue)
% en entree iel : numero de l'element
%           Ue  : vecteur des depalements nodaux de l'element
%
%  H.Oudin 
global ndim
global Coord Connec Nprop Prop 

X  = Coord(Connec(iel,:),:);
dX = X(2,:) - X(1,:);

if ndim == 1
    EI=Prop(Nprop(iel),1);
    L = abs(dX); Ne = 0;
    Te = -(EI/L^3)*(12*Ue(1)+6*L*Ue(2)-12*Ue(3)+6*L*Ue(4));
    Mfe= (EI/L^2)*[-6*Ue(1)-4*L*Ue(2)+6*Ue(3)-2*L*Ue(4) ; 6*Ue(1)+2*L*Ue(2)-6*Ue(3)+4*L*Ue(4)];
    fprintf('Dans l''element %3i\n',iel)
    fprintf('Ne = %8.3e  Te = %8.3e Mfe au noeud i  %8.3e  en noeud j %8.3e \n',Ne,Te,Mfe(1),Mfe(2))
elseif ndim == 2
    ES=Prop(Nprop(iel),1); EI=Prop(Nprop(iel),2);
    L = sqrt(dX(1)^2 + dX(2)^2);
    c = dX(1)/L; s = dX(2)/L; 
    P = [  c  s  0  0  0  0
          -s  c  0  0  0  0
           0  0  1  0  0  0
           0  0  0  c  s  0
           0  0  0 -s  c  0
           0  0  0  0  0  1
         ];
    Ue = P*Ue;
    Ne = (ES/L)*(Ue(4)-Ue(1));
    Te = -(EI/L^3)*(12*Ue(2)+6*L*Ue(3)-12*Ue(5)+6*L*Ue(6));
    Mfe= (EI/L^2)*[-6*Ue(2)-4*L*Ue(3)+6*Ue(5)-2*L*Ue(6) ; 6*Ue(2)+2*L*Ue(3)-6*Ue(5)+4*L*Ue(6)];
    fprintf('Dans l''element %3i\n',iel)
    fprintf('Ne = %8.3e  Te = %8.3e Mfe au noeud i  %8.3e  en noeud j %8.3e \n',Ne,Te,Mfe(1),Mfe(2))
 elseif ndim == 3
     disp('================================================ ');
     disp('       element non programme ');
     disp('================================================ ');
end

return