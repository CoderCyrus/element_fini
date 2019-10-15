function [Ke,Fe] = Q4_ep(iel)
% Calcul de la matrice raideur Ke et de la force generalisee Fe 
% pour un element Q4 d'une structure en elasticite plane
% 
% appel [Ke,Fe] = Q4_ep(iel)
%    ou [Ke,Fe] = feval('Q4_ke',iel)
% en entree iel : numero de l'element
% en sortie Ke  : matrice raideur elementaire (8,8)
%           Fe  : force generalisee elementaire (8,1)
%
%  H.Oudin  
global Coord Connec Nprop Prop 

%X  = Coord(Connec(iel,[1:4]),:);

npg = 4;          %----- integration a 4 points de Gauss
wg = [1,1,1,1];   %----- poids et position
c = 1/sqrt(3); posg = [ -c -c ; c -c ; c  c ; -c  c ];
%npg = 3; wg = [4/3, 4/3, 4/3]; posg = [sqrt(2/3) 0 ; -sqrt(1/6) sqrt(1/2) ; -sqrt(1/6) -sqrt(1/2)];
%npg = 7; wg = [8/7, 20/63, 20/63, 20/36, 20/36, 20/36, 20/36]; 
%c = sqrt(3/5); d = sqrt(14/5); posg = [0 0 ; 0 d ; 0 -d ; -c -c ; c -c ; c c ; -c  c ];

E=Prop(Nprop(iel),1);    %----- matrice d'elasticite D
nu=Prop(Nprop(iel),2);
ep=Prop(Nprop(iel),3);
if ep > 0  a = 0 ; else a = 1 ; ep = 1;  end
coef = ep * E * (1-a*nu)/((1+nu)*(1-nu-a*nu));
D = coef * [     1       nu/(1-a*nu)            0             ;...
            nu/(1-a*nu)      1                  0             ;...
                 0           0       .5*(1-nu-a*nu)/(1-a*nu)];

ndle = 8;            %----- initialisations
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
	%----- matrice B(3x8)
   B=zeros(3,8);
   B(1,[1 3 5 7])=dNx(1,:);
   B(2,[2 4 6 8])=dNx(2,:);			
   B(3,[1 3 5 7,2 4 6 8])=[dNx(2,:),dNx(1,:)];
	%----- matrice Ke(8x8)   %aire=aire+detj*poids;  
   Ke=Ke+(B'*D*B)*detj*poids;
	%----- vecteur Fe(8,1) 
   fx=Prop(Nprop(iel),4); fy=Prop(Nprop(iel),5);
   Fe([1 3 5 7],1) = Fe([1 3 5 7],1)+ ep*fx*detj*poids*N';
   Fe([2 4 6 8],1) = Fe([2 4 6 8],1)+ ep*fy*detj*poids*N';
end
%disp(Ke),disp(Fe)
return