clc;clear all;close all;
% Mapping d'un element Q4 etude de la transformation geometrique
% Exercice NUM4 du chapitre NUM du site
%  H.Oudin 
taille = get(0,'ScreenSize');           
%------ mappage de l'element de reference
[S,T] = meshgrid(-1:0.2:1);
[n,n] = size(S);
%------ donnee de la geometrie de l'element reel (position des 4 noeuds)
disp('---------------------------------------');
disp('Mapping d''un element Q4 etude de la transformation geometrique');
disp('cliquez les 4 noeuds de l''element reel dans cette figure');
disp('---------------------------------------');
figure('Name','donnees de l''element reel',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])    
axis([0 1 0 1])
grid on
[Coord(:,1),Coord(:,2),button] = ginput(4); 
close all;                                      
%Coord=[0 0;8 0;8 2;4 7 ];% peut etre utilise comme donnee
%------ mapping de l'element   
for i=1:n
    for j=1:n
        s=S(j,i);t=T(j,i);
        N = .25*[(1-s)*(1-t)  (1+s)*(1-t)  (1+s)*(1+t)  (1-s)*(1+t)];
        %------ calcul du mapping de l'element reel 
        Xr(j,i)=N*Coord(:,1);Yr(j,i)=N*Coord(:,2);
     %------ calcul de detJ sur le mapping
     dN = .25*[-(1-t)  (1-t) (1+t)  -(1+t)
   		       -(1-s) -(1+s) (1+s)   (1-s)]; 
     J = dN*Coord;
     det(j,i) = J(1,1)*J(2,2)-J(1,2)*J(2,1);
    end
end
C = max(det); maxi=max(C); D = min(det);mini=min(D);
disp('Figures donnant l''evolution de la valeur de detJ');
disp('Le mapping represente les lignes a s ou t constant');
disp('Attention l''element reel est plan, ne vous laissez pas influencer par le trace');
disp('-------------------------------------------------------------------------------');
figure('Name','Valeur de detJ sur le mapping de l''element de reference',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])    
subplot(1,2,1), surf(S,T,det),title('element de reference'),axis equal ...
    , colormap jet, colorbar, view([0 0 1])
subplot(1,2,2), surf(Xr,Yr,det),title('element reel'), axis equal, ...
      view([0 0 1])
% calcul de la surface de l'element par integration numerique
npg = 4;          %----- integration a 4 points de Gauss
wg = [1,1,1,1];   %----- poids et position
c = 1/sqrt(3); posg = [ -c -c ; c -c ; c  c ; -c  c ];
aire=0;
for ipg=1:npg       %----- boucle d'integration
   s = posg(ipg,1); t = posg(ipg,2); poids = wg(ipg);
   dN = .25*[-(1-t)  (1-t) (1+t)  -(1+t)
   		     -(1-s) -(1+s) (1+s)   (1-s)]; 
   J = dN*Coord;
   detj = J(1,1)*J(2,2)-J(1,2)*J(2,1);
   aire=aire+detj*poids;  
end
Coord 
if sign(maxi) == sign(mini) disp('---------------------------------');
    disp('detJ ne s''annule pas sur l''element reel');
    disp('la transformation geometrique est bijective ');
    fprintf('Par integration numerique la surface de l''element reel est : %6.4f\n',aire);
else disp('---------------------------------');
    disp('detJ change de signe sur l''element');
    disp('la transformation geometrique n''est pas bijective ');
    fprintf('La valeur mini du determinant est : %6.4f  et le max : %6.4f \n',mini,maxi);
    fprintf('L''integration numerique donne une surface calculee de %6.4f\n',aire);
end
