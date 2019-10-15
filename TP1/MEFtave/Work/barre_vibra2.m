clc;clear all;close all;
% Methode d'approximation appliquee aux vibrations des barres
% comparaison des solutions numeriques et analytique de l'exercice vibra2 
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
%figure('Name','comparaison fonction de forme - premier mode',...
%        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])            
%fplot(@(x)[sin(pi*x/L) , 4*x.*(L-x)],[0 L]);
disp('Methodes d''approximation sur la formulation forte pour une barre bi-encastree');
disp('------------------------------------------------------------------------------');
disp('Methode de la valeur moyenne');
f=@(x) x.*(L-x);
k = 2*ES*L; m = roS*quad(f,0,L); sol = sqrt(k/m);
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
disp('---------------------------------------');
disp('Methode des residus ponderes avec p(x)=x');
f=@(x) x.^2*(L-x);
k = 2*ES*quad(@(x) x,0,L); m = roS*quad(f,0,L); sol = sqrt(k/m);
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
disp('---------------------------------------');
disp('Methode des residus ponderes avec p(x)=x^2');
f=@(x) x.^3*(L-x);
k = 2*ES*quad(@(x) x^2,0,L); m = roS*quad(f,0,L); sol = sqrt(k/m);
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
disp('---------------------------------------');
disp('Methode de Collocation');
f=@(x) x.*(L-x); xp=0.5;
while xp > 0
  k = 2*ES ; m = roS*f(xp) ; sol = sqrt(k/m);
  fprintf('Pour un point de collocation situe en %3.1f L \n',xp);
  fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
  fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
  disp('---------------------------------------');
  disp('pour quitter la methode de Collocation ne pas donner de nouvelle valeur');
  xp = input('donner une nouvelle position du point de collocation sur [0,1]: ');
  if isempty(xp) xp=0; end 
end
disp('---------------------------------------');
disp('Methode de galerkin');
f=@(x) x.*(L-x); f2 = @(x) (x.*(L-x)).^2;
k = 2*ES*quad(f,0,L); m = roS*quad(f2,0,L); sol = sqrt(k/m);
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
close all;
%----------------------------------------------------------
close all;