function [Ke,Fe] = Q4_th(iel)
% Calcul de la matrice comportement Ke et du vecteur flux Fe 
% pour un element Q4 de thermique lineaire en regime permanent
% 
% appel [Ke,Fe] = Q4_ke(iel)
%    ou [Ke,Fe] = feval('Q4_ke',iel)
% en entree iel : numero de l'element
% en sortie Ke  : matrice raideur elementaire (4,4)
%           Fe  : force generalisee elementaire (4,1)
%
%  H.Oudin  
global Coord Connec Nprop Prop 

npg = 4;          %----- integration a 4 points de Gauss
wg = [1,1,1,1];   %----- poids et position
c = 1/sqrt(3); posg = [ -c -c ; c -c ; c  c ; -c  c ];

D=Prop(Nprop(iel),1);    %----- carateristiques thermiques
r=Prop(Nprop(iel),2);

ndle = 4;            %----- initialisations
Ke = zeros(ndle); Fe = zeros(ndle,1); %  aire=0
for ipg=1:npg       %----- boucle d'integration
   s = posg(ipg,1); t = posg(ipg,2); poids = wg(ipg);
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
	%----- matrice B(2x4)
   B=zeros(2,4);
   B(1,[1 2 3 4])=dNx(1,:);
   B(2,[1 2 3 4])=dNx(2,:);			
	%----- matrice Ke(4x4)   %aire=aire+detj*poids;  
   Ke=Ke+(B'*D*B)*detj*poids;
	%----- vecteur Fe(4,1) 
   Fe([1 2 3 4],1) = Fe([1 2 3 4],1)+ r*detj*poids*N';
end
%disp(Ke),disp(Fe)
return
