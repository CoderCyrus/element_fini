function Ne = barre_stress(iel,Ue)
% Calcul de l'effort normal dans un element d'une structure treillis
%
% appel Ne=feval('barre_stress',iel,Ue)
% en entree iel : numero de l'element
%           Ue  : vecteur des depalements nodaux de l'element
% en sortie l'effort normal dans la barre
% H.Oudin 
global ndim
global Coord Connec Nprop Prop 

ES=Prop(Nprop(iel),1);
X  = Coord(Connec(iel,:),:);
dX = X(2,:) - X(1,:);

if ndim == 1
    L = abs(dX);
    Ne = (ES/L)*(Ue(2)-Ue(1));
elseif ndim == 2
    L = sqrt(dX(1)^2 + dX(2)^2);
    c = dX(1)/L; s = dX(2)/L; 
    Ne=(ES/L)*(c*(Ue(3)-Ue(1))+s*(Ue(4)-Ue(2)));
 elseif ndim == 3
     disp('================================================ ');
     disp('       element non programme ');
     disp('================================================ ');
end
return