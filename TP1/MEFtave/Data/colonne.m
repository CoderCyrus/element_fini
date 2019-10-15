close all ; clc;
%  Donnees de la colonne etudiee dans les exercices de cours
%  Colonne
%  avec h=100 , E = 1000 , S=1 , mg = 6, maillage : 3 elements de degre 1
% 
%  Fonctions utilisees
%   statiqueUR      : calcul de la reponse statique [U,R]
%   plotstr         : trace du maillage avec ne noeuds et elements
%   plotdef         : trace de la deformee
%   resultante      : calcul de la resultante en x, y, z d'un vecteur nodal
%   barre_stress    : calcul de la contrainte dans un element barre
%   barre_compar     : comparaison avec la solution analytique
%
% Initialisation des variables globales
%   nddln : nb de ddl par noeud
% 	nnod  : nb de noeuds 
%  	nddlt : nb de ddl total(=ndln*nnod)   
%	  nelt  : nb d'elements  
% 	nnode : nb de noeuds par element (2)
% 	ndim  : dimension du probleme (1D,2D ou 3D)
%   ncld  : nb de conditions de champ impose (dirichlet)
%
% 	Coord(nnod,ndim): coordonnees des noeuds
%   Connec(nelt,2)	: connectivites	des elements
%   Typel(nelt)     : Type des elements (barre_ke)
% 	Nprop(nelt)		: Ne de caracteristique pour chaque element
%   Prop(nprop,ncar): Tableau des caracteristiques mecaniques (ES, f)     
%	  Ncl(1,nddlt)	: vaut 1 si le ddl est impose (deplacements imposes)
%	  Vcl(1,nddlt)	: valeur du deplacement impose 
%	  F (nddlt,1)		: vecteur des charges nodales donnees

global nddln nnod nddlt nelt nnode ndim ncld
global Coord Connec Typel Nprop Prop Ncl Vcl F
disp(' ');
disp('Structure etudiee : Colonne traitee en exercice de cours N 8');
disp('==================');
% definition du maillage
h = 100;                    % definition des coordonnees des noeuds
Coord=[0 ; h ; 3*h ; 6*h ];
[nnod,ndim]=size(Coord);
nddln=1;  nddlt=nddln*nnod;     
 
Connec=[ 1 , 2 ; ...        % definition de la matrice de connectivite
         2 , 3 ; ...
         3 , 4 ];
[nelt,nnode]=size(Connec);

% definition du modele EF : type des elements
Typel = 'barre_ke';               % definition du type des elements
for i=1:nelt
    Typel = char('barre_ke',Typel);
end

% definition des caracteristiques elementaires (ES fx)
Nprop=[1;1;1];              % pour chaque element Ne de la propriete
Prop =[1000 -0.01];         % tableau des valeurs de ES fx     

% definition des CL en deplacement
CL=[ 1 , 1 ]; ...        % Ne du noeud, type   (1 ddl impose, 0 ddl libre)
Ncl=zeros(1,nddlt);ncld=0;
Vcl=zeros(1,nddlt);      % Valeurs des deplacements imposes
for i=1:size(CL,1)
   for j=1:nddln 
       if CL(i,1+j)==1 Ncl(1,(CL(i,1)-1)*nddln+j)=1; ncld=ncld+1; end
   end
end

% definition des charges nodales
Charg=[ 4  0 ];          %  Ne du noeud , Fx      
F=zeros(nddlt,1);	     %  vecteur sollicitation
for iclf=1:size(Charg,1)           
	noeud=Charg(iclf,1);  
	for i=1:nddln
       F((noeud-1)*nddln+i)=F((noeud-1)*nddln+i) + Charg(iclf,i+1);
    end
end
[Fx,Fy,Fz] = feval('resultante',F);      %----- resultante des charges nodales

% trace du maillage pour validation des donnees 
plotstr   
reponse = '';                  
%%reponse = input('Voulez-vous continuer? O/N [O]: ','s');
if isempty(reponse) | reponse =='O'
    %clc; echo on 
    % echo off
    U = zeros(nddlt,1);
    R = zeros(nddlt,1);
    [U(:,1),R(:,1)] = statiqueUR;   % ----- resolution du probleme

                                %----- format d'impression des vecteurs
    form =' %8.3e   %8.3e   %8.3e  '; format = [form(1:8*nddln),' \n']; 
    disp(' ');disp('------- deplacements nodaux sur (x,y,z) ----------');
    fprintf(format,U)
    %plotdef(U)   
                                            %-----	post-traitement
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
        Ne=feval('barre_stress',iel,Ue);
        fprintf('Dans l''element %3i l''effort normal est %8.3e\n',iel,Ne)
    end

   reponse = input('Voulez-vous comparer avec la solution analytique? O/N [O]: ','s');
   if isempty(reponse) | reponse =='O'
    feval('barre_compar',U);   %----- comparaison avec la solution analytique
   end                          
return
else
    disp(' ');disp('---------------- arret du calcul----------------');
    clear all
end