function  barre_compar(U)
% comparaison de la solution numerique et analytique pour 
% une colonne soumise à son poids propre
%  H.Oudin 
global nddln nelt nnode
global Coord Connec Nprop Prop

fa1=@(x) -(6*x/1000).*(1-x/1200);  % solution analytique
fa2=@(x) -6*(1-x/600);

taille = get(0,'ScreenSize'); 
figure('Name','comparaison avec la solution analytique',...
        'Position',[taille(3)/2.01 taille(4)/2.4 taille(3)/2 taille(4)/2])          
           %----- solution analytique pour l'exemple du chapitre II_3 du cours
subplot(2,1,1), title('deplacement u'), hold on,
fplot(fa1,[0 600],'r')
for iel=1:nelt          %----- boucle sur les elements solution numerique
    loce=[]; for i=1:nnode loce=[loce,(Connec(iel,i)-1)*nddln+[1:nddln]]; end
	Ue=U(loce);
    X = Coord(Connec(iel,:)); x2 = X(2); x1=X(1); L = abs(x2-x1);  
    x=x1:(x2-x1)/5: x2; y=Ue(1)+(Ue(2)-Ue(1))*(x-x1)/L; 
    plot(x,y,'b')
end
grid
subplot(2,1,2), title('effort normal N'), hold on 
fplot(fa2,[0 600],'r') , 
for iel=1:nelt          %----- boucle sur les elements solution numerique
    loce=[]; for i=1:nnode loce=[loce,(Connec(iel,i)-1)*nddln+[1:nddln]];end
	Ue=U(loce);
    X = Coord(Connec(iel,:)); x2 = X(2); x1=X(1); L = abs(x2-x1);  
    ES=Prop(Nprop(iel),1); Ne=(ES/L)*(Ue(2)-Ue(1));
    x=[x1:(x2-x1)/5:x2];
    y=Ne*ones(1,length(x));
    plot(x,y,'b')
end
return