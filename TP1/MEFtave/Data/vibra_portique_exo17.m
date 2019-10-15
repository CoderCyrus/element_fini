close all ; clc;
%  Script de calcul des modes de vibration du portique de l'exo 17
%  portique 2D  avec 2*ne elements
% 
%  fonction utilisee
%    plotstr  : trace du maillage avec numero noeuds et elements
%    vibrations : resolution du probleme
%    plotmodes : trace des modes de vibrations d'un portique 2D
%
%  H.Oudin 

global nddln nnod nddlt nelt nnode ndim
global Coord Connec Nprop Prop Ncl Vcl F
disp(' ');
disp('modes de vibration du portique de l''exercice de cours 17');
disp('==================');
disp('les poutres sont maillees en ne elements');
ne = input('donner le nombre d''elements ne ? [2]: ');
if isempty(ne) ne=2; end 
% definition du maillage
h = 1; nelt=2*ne;
Coord=[]; 
for j=0:ne Coord=[Coord;[0  j*h/ne]]; end 
for j=1:ne Coord=[Coord;[j*h/ne h ]]; end 
[nnod,ndim]=size(Coord);
nddln=3;  nddlt=nddln*nnod;  
Connec=[]; nnode = 2;
for j=1:nelt Connec=[Connec;[j  j+1]]; end    
% tableau des proprietees ES EI rhoS
Nprop = ones(nelt);
Prop=[ 21000000 210 0.78 ];      
% definition des CL en deplacement
CL=[ 1 , 1 , 1 , 1; ...    % numero du noeud, type (1 ddl impose ,0 ddl libre)
     nnod , 1 , 1 , 1];
Ncl=zeros(1,nddlt);ncld=0;
Vcl=zeros(1,nddlt);         % deplacements imposes nuls
for i=1:size(CL,1)
   for j=1:nddln 
       if CL(i,1+j)==1 Ncl(1,(CL(i,1)-1)*nddln+j)=1; end
   end
end
F=zeros(nddlt,1);	    %----- vecteur sollicitation
%plotstr  % trace du maillage                    
nmode = input('combien de frequences faut-il calculer ? [3]: ');
if isempty(nmode) nmode=3; end   
[f,U] = feval('vibrations',nmode);
for j= 1: nmode 
   figure('Name','modes de vibration du portique',...
        'Position',[taille(3)/2.01 taille(4)/2.4 taille(3)/2 taille(4)/2])
   plotmodes(U(:,nmode+1-j),j,f(nmode+1-j))
end
return