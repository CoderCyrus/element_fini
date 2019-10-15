close all ; clc;
%  Exercice 1 du chapitre de cours NUM - elements T3
%  Maillage en 2 element d'une plaque rectangulaire
%
%  Fonctions utilisees
%   statique        : calcul de la reponse statique [U,R]
%   plotstr         : trace du maillage avec numero noeuds et elements
%   plotdef         : trace de la deformee
%   resultante      : calcul de la resultante en x, y, z d'un vecteur nodal
%   T3_ep_stress    : calcul de la contrainte dans un element T3
%   plot_sig        : visualisation de l'etat de contrainte
%
global nddln nnod nddlt nelt nnode ndim ncld
global Coord Connec Typel Nprop Prop Ncl Vcl F
disp(' ');
disp('Maillage en deux T3 de plaque de l''exercice 1 du chapitre NUM');
disp('===================');
% definition du maillage
L = 2 ; h = 1;
Coord=[ 0  0; L  0; L  h;0  h ];   % definition des coordonnees des noeuds
[nnod,ndim]=size(Coord);    
nddln=2;  nddlt=nddln*nnod;  
Connec=[ 1 2 4; ...            % Tableau de connectivite i , j
         2 3 4 ];
[nelt,nnode]=size(Connec);
Typel = 'T3_ep';               % definition du type des elements
for i=1:nelt   Typel = char('T3_ep',Typel); end

% definition des caracteristiques mecaniques
Nprop=[1;1;1;1];           
Prop=[ 210000 0.3 2 0 0 ];     % valeurs  (E nu e fx fy)

% CL en deplacement                                     
CL=[ 1 , 1 , 1; ...            
     2 , 0 , 1; ...
     4 , 1 , 1;];
Ncl=zeros(1,nddlt);ncld=0;
Vcl=zeros(1,nddlt);                        
for i=1:size(CL,1)
   for j=1:nddln 
       if CL(i,1+j)==1 Ncl(1,(CL(i,1)-1)*nddln+j)=1;ncld=ncld+1; end
   end
end
% charges nodales: numero du noeud , Fx,Fy
Charg=[ 3  5  0 ];              
F=zeros(nddlt,1);	   
for iclf=1:size(Charg,1)           
	noeud=Charg(iclf,1);  
	for i=1:nddln
       F((noeud-1)*nddln+i)=F((noeud-1)*nddln+i) + Charg(iclf,i+1);
    end
 end
[Fx,Fy,Fz] = feval('resultante',F);  %----- resultante des charges nodales

plotstr  % trace du maillage 
U = zeros(nddlt,1);
R = zeros(nddlt,1);
[U(:,1),R(:,1)] = statique;   % ----- resolution du probleme (pivot = 1)
%-----	post-traitement
plotdef(U)
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
    loce=[]; 
    for i=1:nnode 
		if Connec(iel,i) > 0 loce=[loce,(Connec(iel,i)-1)*nddln+[1:nddln]]; end
	end;  
	Ue=U(loce);
    sig = feval([deblank(Typel(iel,:)),'_stress'],iel,Ue);
    % contrainte moyenne de Von Mises sur l'element
    sVM = sqrt(sig(1)^2+sig(2)^2-sig(1)*sig(2)+3*sig(3)^2);   
    VM(iel)=sVM ;
    fprintf('Valeur Moyenne des contraintes dans l''element %3i\n',iel)
    fprintf('sigxx = %8.3e  sigyy = %8.3e   sigxy = %8.3e  sigVM = %8.3e \n',sig(1),sig(2),sig(3),sVM)
end
plot_sig(VM)
return