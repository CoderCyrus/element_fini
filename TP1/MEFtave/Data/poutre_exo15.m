close all ; clc;
%  Script pour la poutre de l'exercice de cours exo15
%         poutre sur 2 appuis, modelisee par nelt elements
% 
%  Fonctions utilisees
%   statiqueUR      : calcul de la reponse statique [U,R]
%   plotstr         : trace du maillage avec numero noeuds et elements
%   plotdef         : trace de la deformee
%   resultante      : calcul de la resultante en x, y, z d'un vecteur nodal
%   poutre_stress   : calcul de la contrainte dans un element poutre
%   poutre_compar   : comparaison avec la solution analytique
%   
global nddln nnod nddlt nelt nnode ndim ncld
global Coord Connec Typel Nprop Prop Ncl Vcl F 
disp(' ');
disp('structure etudiee : poutre de l''exercice de cours exo15');
disp('==================');
% definition du maillage
nelt = input('donner le nombre d''elements ne ? [2]: ');
if isempty(nelt) nelt=2; end 
L = 1;
Coord=[]; 
for j=0:nelt Coord=[Coord; j*L/nelt]; end 
[nnod,ndim]=size(Coord);
nddln=2;  nddlt=nddln*nnod;
Connec=[]; nnode = 2;
for j=1:nelt Connec=[Connec;[j  j+1]]; end       
% definition du modele EF : type des elements
Typel = 'poutre_ke';       
for i=1:nelt Typel = char('poutre_ke',Typel); end
% definition des caracteristiques mecaniques elementaires (EI f)  (en 1D)
Nprop = ones(nelt);   % pour chaque element numero de la propriete
Prop=[ 1 -384 ];         % tableau des differentes valeurs de EI fy    
% definition des CL en deplacement
CL=[ 1 , 1 , 0 ; ...     % numero du noeud, (1 ddl impose ,0 ddl libre)
    nelt+1 , 1 , 0 ];
Ncl=zeros(1,nddlt);ncld=0;
Vcl=zeros(1,nddlt);         % Valeurs imposees nulles
for i=1:size(CL,1)
    for j=1:nddln 
    if CL(i,1+j)==1 Ncl(1,(CL(i,1)-1)*nddln+j)=1;ncld=ncld+1; end
    end
end
% definition des charges nodales
F=zeros(nddlt,1);	   
[Fx,Fy,Fz] = feval('resultante',F); 

plotstr  % trace du maillage pour validation des donnees
% ----- resolution du probleme
U = zeros(nddlt,1);
R = zeros(nddlt,1);
[U(:,1),R(:,1)] = statiqueUR;   
plotdef(U)
%-----	post-traitement
form =' %8.3e   %8.3e   %8.3e  '; format = [form(1:8*nddln),' \n']; 
disp(' ');disp('------- deplacements nodaux sur (x,y,z) ----------');
fprintf(format,U)                                     
disp(' ');disp('------- Efforts aux appuis  ----------');
fprintf(format,R(:,1));
[Rx,Ry,Rz] = feval('resultante',R);     %----- resultantes et reactions
disp(' ');
fprintf('La resultante des charges nodales    en (x,y,z) est : %8.3e   %8.3e   %8.3e \n',Fx,Fy,Fz);                    
fprintf('La resultante des charges reparties  en (x,y,z) est : %8.3e   %8.3e   %8.3e \n',-Rx-Fx,-Ry-Fy,-Rz-Fz);
fprintf('La resultante des efforts aux appuis en (x,y,z) est : %8.3e   %8.3e   %8.3e \n',Rx,Ry,Rz);
disp(' ');disp('------- Contraintes sur les elements ----------');
for iel=1:nelt          %----- boucle sur les elements
  loce=[]; for i=1:nnode loce=[loce,(Connec(iel,i)-1)*nddln+[1:nddln]];end
  Ue=U(loce);
  feval('poutre_stress',iel,Ue);
  end
reponse = input('Voulez-vous comparer avec la solution analytique? O/N [O]: ','s');
if isempty(reponse) | reponse =='O'
 feval('poutre_compar',U); %----- comparaison avec la solution analytique
 end                                               
clear all
return