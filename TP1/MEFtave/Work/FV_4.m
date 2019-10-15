clc;clear all;close all;
% comparaison des solutions numeriques et analytique de l'exercice FV4 
% Barre bi-encastree
%  H.Oudin 
%----------------------------------------------------------
% donnees du probleme
ES= 1 ; roS=1; L= 1;
%----------------------------------------------------------
%----- Methodes d'approximation sur la formulation forte
%----- avec la fonction de forme x(L-x)
omeg1=pi*sqrt(ES/roS/L/L);    % solution analytique
taille = get(0,'ScreenSize'); 
figure('Name','comparaison fonction de forme - premier mode',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])            
fplot(@(x)[sin(pi*x/L) , 4*x.*(L-x)],[0 L]);
disp('---------------------------------------');
disp('comparaison de la fonction de forme et du premier mode');
disp('---------------------------------------');
pause(2);
close all;
%Collocation
f=@(x) x.*(L-x); % D2f=@(x) -2*x./x;
disp('Methode de Collocation');
xp=0.5;
while xp > 0
k = 2*ES ; m = roS*f(xp) ; sol = sqrt(k/m);
disp('---------------------------------------');
fprintf('Pour un point de collocation situe en %3.1f L \n',xp);
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
disp('pour quitter la methode de Collocation ne pas donner de nouvelle valeur');
xp = input('donner une nouvelle position du point de collocation sur [0,1]: ');
if isempty(xp) xp=0; end 
end
disp('---------------------------------------');
disp('Methode de la valeur moyenne');
k = 2*ES*L; m = roS*quad(f,0,L); sol = sqrt(k/m);
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
disp('---------------------------------------');
disp('Methode de galerkin');
f2 = @(x) (x.*(L-x)).^2;
k = 2*ES*quad(f,0,L); m = roS*quad(f2,0,L); sol = sqrt(k/m);
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
close all;
%----------------------------------------------------------
%------ Methode des elements finis
disp('---------------------------------------');
disp('Methode des elements finis');
ne=2;
while ne > 1
    K=zeros(ne-1);	M=zeros(ne-1);	 
figure('Name',['comparaison EF - analytique pour un mailage de '...
    ,num2str(ne),' elements'],...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])            
Le=L/ne;
if ne <= 2 K(1,1)=2*ES/Le; M(1,1)=4*roS*Le/6;end
for i=1:ne-2         
    K(i,i)=2*ES/Le;  M(i,i)=4*roS*Le/6;
    K(i,i+1)=-ES/Le; M(i,i+1)=roS*Le/6;
    K(i+1,i)=-ES/Le; M(i+1,i)=roS*Le/6;
    K(i+1,i+1)=2*ES/Le;  M(i+1,i+1)=4*roS*Le/6;
end
if ne >= 4 nmod=3; else nmod =ne-1; end
[modes,omega] = eigs(K,M,ne-1,'sm'); 
omeg = diag(sqrt(omega));
disp('---------------------------------------');
fprintf('Pour un maillage de %2d elements finis \n',ne);
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',omeg(1),omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(omeg(1)-omeg1)/omeg1);

for imod = 1:nmod
    U= modes(:,imod)  ; maxi = max(abs(U));
    if U(1)<0  U = -U ;end 
    U = [ 0 ;U ;0]/maxi;
% trace des modes 
    for i=1:ne
        x1=(i-1)*Le; x2=i*Le;
        x=x1:Le/10:x2; y=U(i)+(U(i+1)-U(i))*(x-x1)/Le; 
        subplot(nmod,1,imod),hold on, plot(x,y,'b');
    end
sol = omeg(imod); anal=imod*pi*sqrt(ES/roS/L/L);
err = 100*abs(sol-anal)/anal;
subplot(nmod,1,imod), fplot(@(x) sin(imod*pi*x/L),[0 L],'r'),...
title([num2str(imod),'ieme mode : pulsation approchee ',num2str(sol),...
  ' / analytique ',num2str(anal),' soit une erreur de ',num2str(err),...
  ' % ' ]), legend('sol. analytique en rouge','Sol. EF en bleu'), grid
end
disp('pour quitter l''application ne pas donner de nouvelle valeur');
ne = input('donner le nombre d''elements ne ? [>2]: ');
if isempty(ne) ne=0; end 
end
close all;