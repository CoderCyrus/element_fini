close all ; clear; clc;
%  Script de calcul des modes de vibration de la barre encastree - libre
%  exemple du cours modelisee par nelt elements
% 
%  fonction utilisee
%    plotstr  : trace du maillage avec numero noeuds et elements
%    vibrations : resolution du probleme
%
%  H.Oudin 
disp(' ');
disp('modes de vibration d''une barre encastree - libre');
disp('==================');
global nddln nnod nddlt nelt nnode ndim
global Coord Connec Typel Nprop Prop Ncl Vcl F

% definition du maillage
nelt = input('donner le nombre d''elements nelt ? [2]: ');
if isempty(nelt) nelt=2; end 
L = 1;
Coord=[]; 
for j=0:nelt Coord=[Coord; j*L/nelt]; end 
[nnod,ndim]=size(Coord);
nddln=1;  nddlt=nddln*nnod;
Connec=[]; nnode = 2;
for j=1:nelt Connec=[Connec;[j  j+1]]; end 

% definition du modele EF : type des elements
Typel = 'barre_keme';         
for i=1:nelt Typel = char('barre_keme',Typel); end

% definition des caracteristiques mecaniques  (en 1D)
Nprop = ones(nelt);   % pour chaque element numero de la propriete
Prop=[ 1 1 ]; % tableau de ES RhoS   

% definition des CL en deplacement
CL=[ 1 , 1 ];
Ncl=zeros(1,nddlt);
Vcl=zeros(1,nddlt);         % deplacements imposes nuls
for i=1:size(CL,1)
   for j=1:nddln 
       if CL(i,1+j)==1 Ncl(1,(CL(i,1)-1)*nddln+j)=1; end
   end
end
F=zeros(nddlt,1);
%plotstr  % trace du maillage                    
nmode = input('combien de frequences faut-il calculer ? [2]: ');
if isempty(nmode) nmode=2; end   
f=[];  U=[]; 
[f,U] = feval('vibrations',nmode);
disp('frequences obtenues par le modele EF');
f
if f(1)>f(2) 
  ft=fliplr(f');
  U=fliplr(U);
  f=ft';
end

ES=Prop(1,1);RS=Prop(1,2);
norme = RS*L/2; 
taille = get(0,'ScreenSize'); 
figure('Name','Comparaison analytique / numerique des modes de vibration de la poutre',...
        'Position',[taille(3)/2.01 taille(4)/2.4 taille(3)/2 taille(4)/2])

for j= 1: nmode 
   x=[0:0.01*L:L];  y = sin((2*j-1)*pi*x/(2*L));%----- sol analytique appuyeee - appuyee
   fanal = (2*j-1)*sqrt(ES/(RS*(L^2)))/4;err=100*(f(j)-fanal)/fanal;
   subplot(nmode,1,j),title([int2str(j),' mode de vibration a la frequence ',...
   num2str(f(j),'%7.2f'),'  / sol analytique ', num2str(fanal,'%7.2f'), ...
   '  soit une erreur de ', num2str(err,'%4.2f'), ' %']),hold on
   plot(x,y,'g'); plot(x,-y,'g'); axis([0  L  -1.5  1.5])
   s = sqrt(norme);%---- facteur d'echelle sur la deformee EF a partir de la norme analytique
   for iel = 1:nelt
    Pos  = Coord(Connec(iel,:),:); Le=Pos(2,1)-Pos(1,1);
    je=(Connec(iel,1)-1)*nddln+1;
    x = [0:0.1:1];
    y = s*( (1-x)*U(je,j) + x.*U(je+1,j));
    x = Pos(1,1)+ x*Le;
    plot(x,y,'b'); plot(x,-y,'b')
   end       
end
return