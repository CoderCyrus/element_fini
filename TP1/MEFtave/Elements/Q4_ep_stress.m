function sig = Q4_ep_stress(iel,Ue)
% Calcul du tenseur "contrainte moyenne" sur un element Q4
% en elasticite plane  
% 
% appel [sig] = Q4_ep_stress(iel,Ue)
%    ou [sig] = feval('Q4_stress',iel,Ue)
% en entree iel : numero de l'element
%           Ue  : vecteur des depalements nodaux de l'element
% en sortie Tenseur des contraintes moyennes sur l'element
%
%  H.Oudin 
global Coord Connec Nprop Prop 

                  %-----  position des points de calcul de la contrainte
npg = 4; c = 1/sqrt(3); posg = [ -c -c ; c -c ; c  c ; -c  c ]; %points de gauss
%npg = 4;  c = 1;         posg = [ -c -c ; c -c ; c  c ; -c  c ]; %noeuds de l'element
%npg = 1; posg = [ 0 0 ]; %au centre


E=Prop(Nprop(iel),1);    %----- matrice d'elasticite D
nu=Prop(Nprop(iel),2);
ep=Prop(Nprop(iel),3);
if ep > 0  a = 0 ; else a = 1 ; ep = 1;  end
coef = E * (1-a*nu)/((1+nu)*(1-nu-a*nu));
D = coef * [     1       nu/(1-a*nu)            0             ;...
            nu/(1-a*nu)      1                  0             ;...
                 0           0       .5*(1-nu-a*nu)/(1-a*nu)];

ndle = 8;            %----- initialisations
sigpg = zeros(4,3);sig = zeros(1,3);

%fprintf('Valeur des contraintes en des points de l''element %3i\n',iel)
for ipg=1:npg       %----- boucle sur les points 
   s = posg(ipg,1); t = posg(ipg,2);
    %----- vecteur <N(s,t)>
   N = .25*[(1-s)*(1-t)  (1+s)*(1-t)  (1+s)*(1+t)  (1-s)*(1+t)]; 
    %----- matrice [dN/ds ;dN/dt]
   dN = .25*[-(1-t)  (1-t) (1+t)  -(1+t)
   		     -(1-s) -(1+s) (1+s)   (1-s)]; 
	%----- matrice jacobienne
   J = dN*Coord(Connec(iel,[1:4]),:);
   detj = J(1,1)*J(2,2)-J(1,2)*J(2,1);
   J_1 = [J(2,2) -J(1,2); -J(2,1) J(1,1)]/detj ;
    %----- matrice [dN/dx ;dN/dy]
   dNx = J_1*dN;
	%----- matrice B(3x8)
   B=zeros(3,8);
   B(1,[1 3 5 7])=dNx(1,:);
   B(2,[2 4 6 8])=dNx(2,:);			
   B(3,[1 3 5 7,2 4 6 8])=[dNx(2,:),dNx(1,:)];
	%----- tableau des contraintes aux points de Gauss   
   sigpg(ipg,:) = (D*B*Ue)';
   sig(1,:) = sig(1,:) + sigpg(ipg,:);
  % fprintf('point (s,t) %5.2f %5.2f ',s,t)
  % fprintf('sigxx = %8.3e  sigyy = %8.3e   sigxy = %8.3e  \n',sigpg(ipg,1),sigpg(ipg,2),sigpg(ipg,3))
end
%----- Valeur moyenne du tenseur des contraintes
sig = sig/npg;
return