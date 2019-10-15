function  poutre_compar(U)
% comparaison de la solution numerique et analytique pour 
% une poutre uniformement chargee sur deux appuis
%  H.Oudin 
global nddln nelt nnode
global Coord Connec Nprop Prop

fa1=@(x) Prop(1,2)*x.*(x-1)/2 ;  % solution analytique
fa2=@(x) -Prop(1,2)*(2*x-1)/2;
taille = get(0,'ScreenSize'); 
figure('Name','comparaison avec la solution analytique',...
      'Position',[taille(3)/2.02 taille(4)/2.6 taille(3)/2 taille(4)/2]) 

subplot(2,1,1),title('moment de flexion Mf'), hold on,     
fplot(fa1,[0 1],'r')   %----- solution analytique
for iel=1:nelt          %----- boucle sur les elements solution numerique
  loce=[]; for i=1:nnode loce=[loce,(Connec(iel,i)-1)*nddln+[1:nddln]]; end
	Ue=U(loce);
  EI=Prop(Nprop(iel),1);
  X = Coord(Connec(iel,:)); x2 = X(2); x1=X(1); L = abs(x2-x1); 
  Mfe= (EI/L^2)*[-6*Ue(1)-4*L*Ue(2)+6*Ue(3)-2*L*Ue(4),... 
                  6*Ue(1)+2*L*Ue(2)-6*Ue(3)+4*L*Ue(4)];
  x=x1:(x2-x1)/5:x2; y=Mfe(1)+(Mfe(2)-Mfe(1))*(x-x1)/L;
  plot(x,y,'b')
end
grid


subplot(2,1,2)
hold on,       
fplot(fa2,[0 1],'r')
title('effort tranchant'),
for iel=1:nelt          %----- boucle sur les elements solution numerique
    loce=[]; for i=1:nnode loce=[loce,(Connec(iel,i)-1)*nddln+[1:nddln]];end
	Ue=U(loce);
    EI=Prop(Nprop(iel),1);
    X = Coord(Connec(iel,:)); x2 = X(2); x1=X(1); L = abs(x2-x1); 
    Te = -(EI/L^3)*(12*Ue(1)+6*L*Ue(2)-12*Ue(3)+6*L*Ue(4));
    x=[x1:(x2-x1)/5:x2];
    y=Te*ones(1,length(x));
    plot(x,y,'b'),
end
grid
return