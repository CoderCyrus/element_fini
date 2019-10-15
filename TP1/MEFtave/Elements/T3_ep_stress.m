function sig = T3_ep_stress(iel,Ue)
% Calcul du champ de contrainte sur un element T3
% en elasticite plane  (c'est une constante)
% 
% appel [sig] = T3_ep_stress(iel,Ue)
%    ou [sig] = feval('T3_stress',iel,Ue)
% en entree iel : numero de l'element
%           Ue  : vecteur des depalements nodaux de l'element
% en sortie Tenseur des contraintes moyennes sur l'element
%
%  H.Oudin 
global Coord Connec Nprop Prop 

nnode=3; ndle = 6;                %----- initialisations

E=Prop(Nprop(iel),1);    %----- matrice d'elasticite D
nu=Prop(Nprop(iel),2);
ep=Prop(Nprop(iel),3);
if ep > 0  a = 0 ; else a = 1 ; ep = 1;  end
coef = E * (1-a*nu)/((1+nu)*(1-nu-a*nu));
D = coef * [     1       nu/(1-a*nu)            0             ;...
            nu/(1-a*nu)      1                  0             ;...
                 0           0       .5*(1-nu-a*nu)/(1-a*nu)];

X  = Coord(Connec(iel,[1:nnode]),:);    % coordonnees des noeuds de l'element

                                    %-----  determinant detj = 2A
detj=(X(2,1)-X(1,1))*(X(3,2)-X(1,2))-(X(3,1)-X(1,1))*(X(2,2)-X(1,2));
        
B=zeros(nnode,ndle);                       %-----  matrice B(3,6) 
B(1,[1 3 5])=(X([2 3 1],2)-X([3 1 2],2) )' /detj;
B(2,[2 4 6])=(X([3 1 2],1)-X([2 3 1],1) )' /detj ;
B(3,[1 3 5 2 4 6]) =[ B(2,[2 4 6]), B(1,[1 3 5]) ];
        
sig = (D*B*Ue)';              %----- vecteur des contraintes
return