close all ; clc;
%  Treillis de l'exercice MEF6 en tenant compte de la symetrie
%  H.Oudin 
global nddln nnod nddlt nelt nnode ndim ncld
global Coord Connec Typel Nprop Prop Ncl Vcl F
disp(' ');
disp('Structure etudiee  : Exercice MEF6 avec symetrie');
disp('=================');
% definition du maillage
h = 100*sqrt(2)/2;  
Coord=[ 0   0; h  -h ; 2*h  0 ; 2*h  -h ];
[nnod,ndim]=size(Coord);
nddln=2;  nddlt=nddln*nnod;     
Connec=[ 1 2 ; 1 3 ; 2 3 ; 2 4 ];
[nelt,nnode]=size(Connec); 
Typel = 'barre_ke';              
for i=1:nelt
    Typel = char('barre_ke',Typel);
end
Nprop=ones(nelt);
Prop=[100*sqrt(2) 0 0 ];   
CL=[ 1 , 0 , 1 ; ...     
     3 , 1 , 0 ; ...
     4 , 1 , 1 ; ...
     ]; 
Ncl=zeros(1,nddlt); ncld=0;
Vcl=zeros(1,nddlt);     
for i=1:size(CL,1)
   for j=1:nddln 
       if CL(i,1+j)==1 
           Ncl(1,(CL(i,1)-1)*nddln+j)=1;
           ncld=ncld+1; 
       end
   end
end
Charg=[ 2    0  -1.];
F=zeros(nddlt,1);	    
for iclf=1:size(Charg,1)           
	noeud=Charg(iclf,1);  
	for i=1:nddln
       F((noeud-1)*nddln+i)=F((noeud-1)*nddln+i) + Charg(iclf,i+1);
    end
end
[Fx,Fy,Fz] = feval('resultante',F);     

% trace du maillage 
plotstr                     
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
   Ne=feval('barre_stress',iel,Ue);
   fprintf('Dans l''element %3i l''effort normal est %8.3e\n',iel,Ne)
   end  
disp(' ');   
disp('Expliquer pourquoi les resultats sont justes,');
disp('alors que le noeud 4 est bloque !!');
return


