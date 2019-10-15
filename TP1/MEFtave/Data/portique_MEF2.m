close all ; clc;
%  script du portique de l'exercice MEF2
%                 
%  Fonctions utilisees
%   statiqueUR      : calcul de la reponse statique [U,R]
%   plotstr         : trace du maillage avec numero noeuds et elements
%   plotdef         : trace de la deformee
%   resultante      : calcul de la resultante en x, y, z d'un vecteur nodal
%   poutre_stress   : calcul de la contrainte dans un element barre
%  
global nddln nnod nddlt nelt nnode ndim ncld
global Coord Connec Typel Nprop Prop Ncl Vcl F 
disp(' ');
disp('structure etudiee : portique de l''exercice MEF2');
disp('==================');
% definition du maillage
h = 1; L = 1;
Coord=[ 0 , 0 ; ...         % definition des coordonnees des noeuds X , Y
        0 , h ; ...
        L , h ; ...
        L , 0];
[nnod,ndim]=size(Coord);
nddln=3;  nddlt=nddln*nnod;     
Connec=[ 1 , 2 ; ...        % definition de la matrice de connectivite i , j
         2 , 3 ; ...
         3 , 4];
[nelt,nnode]=size(Connec);

% definition du modele EF : type des elements
Typel = 'poutre_ke';     
for i=1:nelt Typel = char('poutre_ke',Typel); end

% definition des caracteristiques mecaniques elementaires (ES EI fx fy)
Nprop=[1;2;1];          
Prop=[ 1000 1 0 0 ; 1000 1 0 0 ]; 
     
% definition des CL en deplacement
CL=[ 1 , 1 , 1 , 1; ...        % 1 ddl impose ,0 ddl libre
     4 , 1 , 1 , 1];
Ncl=zeros(1,nddlt);ncld=0;
Vcl=zeros(1,nddlt);         % Valeurs imposees nulles
for i=1:size(CL,1)
   for j=1:nddln 
       if CL(i,1+j)==1 Ncl(1,(CL(i,1)-1)*nddln+j)=1;ncld=ncld+1; end
   end
end

% definition des charges nodales
Charg=[ 2  84  0  0 ];      %  numero du noeud , Fx,Fy,Mz      
F=zeros(nddlt,1);	          %  vecteur sollicitation
for iclf=1:size(Charg,1)           
	noeud=Charg(iclf,1);  
	for i=1:nddln
       F((noeud-1)*nddln+i)=F((noeud-1)*nddln+i) + Charg(iclf,i+1);
    end
 end
[Fx,Fy,Fz] = feval('resultante',F);      %----- resultante des charges nodales

plotstr  % trace du maillage pour validation des donnees 

U = zeros(nddlt,1);
R = zeros(nddlt,1);
[U(:,1),R(:,1)] = statiqueUR;   % ----- resolution du probleme
%-----	post-traitement
plotdef(U)
%----- format d'impression des vecteurs
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
for iel=1:nelt   
    loce=[]; for i=1:nnode loce=[loce,(Connec(iel,i)-1)*nddln+[1:nddln]];end
    Ue=U(loce);
    feval('poutre_stress',iel,Ue);
end                     
return