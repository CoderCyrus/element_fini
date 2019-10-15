close all ; clc;
%  Script de calcul des modes de vibration de la poutre de l'exo 15
%  exemple du cours modelisee par nelt elements
% 
%  fonction utilisee
%    plotstr  : trace du maillage avec numero noeuds et elements
%    vibrations : resolution du probleme
%
%  H.Oudin 
disp(' ');
disp('modes de vibration d''une poutre appuyee - appuyee');
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
nddln=2;  nddlt=nddln*nnod;
Connec=[]; nnode = 2;
for j=1:nelt Connec=[Connec;[j  j+1]]; end 

% definition du modele EF : type des elements
Typel = 'poutre_keme';         
for i=1:nelt Typel = char('poutre_keme',Typel); end

% definition des caracteristiques mecaniques elementaires (EI f)  (en 1D)
Nprop = ones(nelt);   % pour chaque element numero de la propriete
Prop=[ 21000000 210 0.78 ]; % tableau de ES EI RhoS   

% definition des CL en deplacement
CL=[ 1 , 1 , 0 ; ...     % numero du noeud, (1 ddl impose ,0 ddl libre)
    nelt+1 , 1 , 0 ];
Ncl=zeros(1,nddlt);
Vcl=zeros(1,nddlt);         % deplacements imposes nuls
for i=1:size(CL,1)
   for j=1:nddln 
       if CL(i,1+j)==1 Ncl(1,(CL(i,1)-1)*nddln+j)=1; end
   end
end
F=zeros(nddlt,1);
%plotstr  % trace du maillage                    
nmode = input('combien de frequences faut-il calculer ? [3]: ');
if isempty(nmode) nmode=3; end
f=[];  U=[]; 
[f,U] = feval('vibrations',nmode);
disp('frequences obtenues par le modele EF');
f
if f(1)>f(2) 
  ft=fliplr(f');
  U=fliplr(U);
  f=ft';
  end

EI=Prop(1,2);RS=Prop(1,3);
norme = RS*L/2;
disp('en bleu les modes EF / en vert la solution analytique');
taille = get(0,'ScreenSize'); 
figure('Name','Comparaison analytique / numerique des modes de vibration de la poutre',...
        'Position',[taille(3)/2.01 taille(4)/2.4 taille(3)/2 taille(4)/2])

for j= 1: nmode 
   x=[0:0.01*L:L];  y = sin(j*pi*x/L);%----- sol analytique appuyeee - appuyee
   fanal = (j^2)*sqrt(EI/(RS*(L^4)))*pi/2;
   subplot(nmode,1,j),title([int2str(j),' mode de vibration a la frequence ',...
   num2str(f(j),'%7.2f'),'  / sol analytique ', ...
   num2str(fanal,'%7.2f')]),hold on
   plot(x,y,'g'); plot(x,-y,'g'); axis([0  L  -1.2  1.2])
   s = sqrt(norme);%---- facteur d'echelle sur la deformee EF a partir de la norme analytique
   for iel = 1:nelt
    Pos  = Coord(Connec(iel,:),:); Le=Pos(2,1)-Pos(1,1);
    je=(Connec(iel,1)-1)*nddln+1;
    x = [0:0.1:1]; x2=x.*x ; x3=x.*x2;
    y = s*( (1-3*x2+2*x3)*U(je,j) + Le*(x-2*x2+x3)*U(je+1,j) + ...
    (3*x2-2*x3)*U(je+2,j) + Le*(-x2+x3)*U(je+3,j) );
    x = Pos(1,1)+ x*Le;
    plot(x,y,'b'); plot(x,-y,'b')
   end       
end
return