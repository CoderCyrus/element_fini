close all ; clc;
%  NUM_exo2
%  Poutre console avec un mailllage automatique en Q4
%  H.Oudin 

global nddln nnod nddlt nelt nnode ndim ncld
global Coord Connec Typel Nprop Prop Ncl Vcl F
disp(' ');
disp('Elasticite plane : Poutre console avec un mailllage automatique en Q4');
disp('=================');

% donnees unites N , mm
L = 3000 ; h =200 ; e=10;    
E=210.e+3; nu=0.3;
poid=-200;

% ------ maillage
nex = input('donner le nombre d''element en x (longueur) ? [5]: ');
if isempty(nex) nex=5; end  
ney = input('donner le nombre d''element en y (hauteur) ? [3]: ');
if isempty(ney) ney=3; end           
Coord=[];                                                
for i=0:nex
    for j=0:ney Coord=[Coord;[i*L/nex  j*h/ney]]; end 
 end
Connec=[]; pas=ney+1;   
for i=1:nex
    for j=1:ney 
    i1=(i-1)*pas+j;
    Connec=[Connec;[i1  i1+pas i1+pas+1 i1+1]];
    end 
 end
[nnod,ndim]=size(Coord);
[nelt,nnode]=size(Connec);
nddln=2;  nddlt=nddln*nnod;               
for i=1:nelt  Typel = char('Q4_ep',Typel); end

% ------ caracteristiques mecaniques (E nu e fx fy)
Nprop=ones(nex*ney);   
Prop=[ E nu e 0 0 ];

% ------ definition des CL en deplacement
CL=[];          
    for i=0:ney  CL=[CL;[i+1,  1 , 1 ]];  end
Ncl=zeros(1,nddlt);ncld=0;
Vcl=zeros(1,nddlt);        
for i=1:size(CL,1)
   for j=1:nddln 
       if CL(i,1+j)==1 Ncl(1,(CL(i,1)-1)*nddln+j)=1;ncld=ncld+1; end
   end
end

% ------ definition des charges nodales
Charg=[(ney+1)*(nex+1) , 0, poid ];     
F=zeros(nddlt,1);	 
for iclf=1:size(Charg,1)           
	noeud=Charg(iclf,1);  
	for i=1:nddln
       F((noeud-1)*nddln+i)=F((noeud-1)*nddln+i) + Charg(iclf,i+1);
    end
 end
% ----- trace du maillage
%plotstr  
% ----- resolution du probleme
U = zeros(nddlt,1);
R = zeros(nddlt,1);
[U(:,1),R(:,1)] = statiqueUR;   
plotdef(U)
VM = [];
% -----	post-traitement

sige=[];   % tableau des contraintes sigxx
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
    %fprintf('Valeur Moyenne des contraintes dans l''element %3i\n',iel)
    %fprintf('sigxx = %8.3e  sigyy = %8.3e   sigxy = %8.3e  sigVM = %8.3e \n',sig(1),sig(2),sig(3),sVM)
    sige=[sige;[sig(1)]];
end
disp(' ');disp('------- Contraintes moyenne de Von Mises sur les elements ----------');
for iel=1:nelt         
 fprintf('element %3i sigma VM : %8.3e \n',iel,VM(iel));
end
taille = get(0,'ScreenSize'); 
figure('Name','Contrainte moyenne de Von Mises sur les elements',...
        'Position',[taille(3)/2.01 taille(4)/2.8 taille(3)/2 taille(4)/2])  
plot_sig(VM) 

% comparaison avec la solution poutre
I=e*h^3/12; EI=E*I ;        % solution poutre
v=poid*(L^3)/(3*EI); sigxx=poid*L*h/(2*I);
disp(' ');
fprintf('Solution analytique poutre\n')
fprintf('la fleche en bout est = %8.3e mm et la contrainte max sigxx = %8.3e  MPa \n',v,sigxx)
disp(' ');
T=linspace(0,L,1000);
for i=1:size(T,2)
    Sol(i)=poid*(T(i)^2*(3*L-T(i)))/(6*EI);
end
taille = get(0,'ScreenSize');
figure('Name','comparaison des fleches 2D-Q4 (rouge) / poutre (bleu)',...
       'Position',[taille(3)/2.02 taille(4)/2.6 taille(3)/2 taille(4)/2] )  
plot(T,Sol);
hold on
Uv=U([2:2*(ney+1):2*(nex+1)*(ney+1)]);
Tv=linspace(0,L,size(Uv,1));
plot(Tv,Uv,'r');
figure('Name','Contraintes max sigxx 2D-Q4 (rouge) / poutre (bleu)',...
       'Position',[taille(3)/2.02 taille(4)/2.6 taille(3)/2 taille(4)/2] )   
hold on,
% La contrainte sigxx poutre est calculée au centre des éléments de la rangée du bas
x=0:L; y=poid*h*(1-1/ney)*(L-x)/(2*I); plot(x,y,'b'),
for i=1:nex      %contrainte sigxx sur les éléments de la rangée du bas
    x1=(i-1)*L/nex ; x2 = i*L/nex ;  
    sign=sige(1+ney*(i-1));
    x=x1:(x2-x1)/20: x2;
    y=sign*ones(1,length(x));
    plot(x,y,'r')
end
grid 
return