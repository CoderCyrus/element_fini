close all ; clc;
%  Script de calcul statique du treillis traite en exemple dans le cours
%  avec a=100 , ES = 100*sqrt(2) , F = 40.
% 
% Fonctions utilisees
%   statiqueUR      : calcul de la reponse statique [U,R]
%   plotstr         : trace du maillage avec ne noeuds et elements
%   plotdef         : trace de la deformee
%   resultante      : calcul de la resultante en x, y, z d'un vecteur nodal
%   barre_stress    : calcul de la contrainte dans un element barre
%

% Initialisation des variables globales pour un treillis
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
% 	Nprop(nelt)		  : Ne de caracteristique pour chaque element
%   Prop(nprop,ncar): Tableau des caracteristiques mecaniques (ES, f)     
%	  Ncl(1,nddlt)	  : vaut 1 si le ddl est impose (deplacements imposes)
%	  Vcl(1,nddlt)	  : valeur du deplacement impose 
%	  F (nddlt,1)		  : vecteur des charges nodales donnees
%
% L'objectif de scripts de donnees est d'initialiser ses variables globales,
% avant de lancer les calculs et d'exploiter les resultats (post-traitement).
% 
global nddln nnod nddlt nelt nnode ndim ncld
global Coord Connec Typel Nprop Prop Ncl Vcl F tab_Lc
disp(' ');
disp('structure etudiee : treillis traite en exemple dans le cours');
disp('==================');
% definition du maillage

Coord=[ 0 , 0 ; ...         % definition des coordonnees des noeuds X , Y
        3 , -1.5 ; ...
        5 , 0 ];
[nnod,ndim]=size(Coord);
nddln=2;  nddlt=nddln*nnod;     
 
Connec=[ 1 , 2 ; ...        % definition de la matrice de connectivite i , j
         1 , 3 ; ...
         2 , 3 ];

[nelt,nnode]=size(Connec);

%matrice pour N
global Nbis;
Nbis=[];
%definition du modele EF : type des elements
Typel = 'barre_ke';               % definition du type des elements
for i=1:nelt
    Typel = char('barre_ke',Typel);
end
% definition des caracteristiques mecaniques elementaires (ES fx fy)
Nprop=[1;2;1];              % pour chaque element Ne de la propriete
Prop=[ 210*10^9*49*10^-4 0 0; ...     % tableau des differentes valeurs de ES fx fy
      210*10^9*25*10^-4 0 0];
      
% definition des CL en deplacement
CL=[ 1 , 0 , 1 ; ...  % Ne du noeud, type sur u et v (1 ddl impose ,0 ddl libre)
     3 , 1 , 1 ];
Ncl=zeros(1,nddlt); ncld=0;
Vcl=zeros(1,nddlt);  % Valeurs des deplacements imposes
%Vcl(2)=1;           % à utiliser pour imposer une valeur non nulle sur un ddl i
for i=1:size(CL,1)
   for j=1:nddln 
       if CL(i,1+j)==1 
           Ncl(1,(CL(i,1)-1)*nddln+j)=1;
           ncld=ncld+1; 
       end
   end
end

% definition des charges nodales
Charg=[ 2  0  -950000                 %  Ne du noeud , Fx , Fy
      ];
F=zeros(nddlt,1);	    %----- vecteur sollicitation
for iclf=1:size(Charg,1)           
	noeud=Charg(iclf,1);  
	for i=1:nddln
       F((noeud-1)*nddln+i)=F((noeud-1)*nddln+i) + Charg(iclf,i+1);
    end
end
[Fx,Fy,Fz] = feval('resultante',F);      %----- resultante des charges nodales

disp('Les variables globales sont initialisees');
disp('Fin de lecture des donnees');
% trace du maillage pour validation des donnees 
plotstr                     
reponse = input('Voulez-vous continuer? O/N [O]: ','s');
if isempty(reponse) | reponse =='O'
    U = zeros(nddlt,1);
    R = zeros(nddlt,1);
    [U(:,1),R(:,1)] = statiqueUR;   % ----- resolution du probleme
%----- format d'impression des vecteurs
    form =' %8.3e   %8.3e   %8.3e  '; format = [form(1:8*nddln),' \n']; 
    disp(' ');disp('------- deplacements nodaux sur (x,y,z) ----------');
    fprintf(format,U)
    plotdef(U)
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
        Ne = feval('barre_stress',iel,Ue);
        Nbis=[Nbis,Ne];
        fprintf('Dans l''element %3i l''effort normal est %8.3e\n',iel,Ne)
    end 
    E=210*10^9*ones(nelt,1);
    Re=350*10^6*ones(1,nelt);
  
    % function g�nerale avec cylindre
  %  sigmamax =  verification(Nbis,Nprop,E,Prop,Re)
[critere1 , critere2] =  verification(Nbis,Nprop,E,Prop,Re,Connec,Coord)
%critere1=plastification
    for i=1:nelt  
        %----- boucle sur les elements
        if critere1(1,i)
            disp(i)
              fprintf('L''element'  ,i,' ne plastifie pas car sigma <\n', Re(1,i), 'Pa')
        else
            fprintf('L''element ',iel,' plastifie car sigma >', Re(1,iel), 'Pa')
        end
        if critere2(1,iel)
            fprintf('L''element ',iel,' ne flambe pas')
        else
            fprintf('L''element',iel,' va flamber')
        end
    end
    clear all
    return
    else
disp(' ');disp('---------------- arret du calcul----------------');
clear all
    end