clc;clear all;close all;
%----- Methodes d'approximation sur la formulation forte
% comparaison des solutions numeriques et analytique de l'exercice de cours 03
% Vibration d'une barre encastree - ressort
%  H.Oudin 
%----------------------------------------------------------
% donnees du probleme
ES= 1 ; roS=1; L= 1;
alpa = input('donner la valeur du rapport ES/kl ? [10]: ');
if isempty(alpa) alpa=10; end
a= pi/2+0.01;                    % solution analytique
b=3*pi/2-0.01;
lx = fzero(@(x) tan(x)+alpa*x,[a b]);
omeg1=lx*sqrt(ES/roS)/L
%----------------------------------------------------------
disp('Methodes d''approximation sur la formulation forte');
disp('avec une fonction de comparaison polynomiale de degre 2');
disp('---------------------------------------');
beta = (alpa+1)/(2*alpa+1)
f=@(x) x.*(1-beta*x/L);
disp('Methode de Collocation');
xp=0.5;
while xp > 0
  k = 2*beta*ES/L ; m = roS*f(xp*L) ; sol = sqrt(k/m);
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
k = 2*beta*ES ; m = roS*quad(f,0,L); sol = sqrt(k/m);
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
disp('---------------------------------------');
disp('Methode de galerkin');
f2 = @(x) (x.*(1-beta*x/L)).^2;
k = 2*beta*ES*quad(f,0,L)/L ; m = roS*quad(f2,0,L); sol = sqrt(k/m);
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
%----------------------------------------------------------
disp('------------------------------------------');
disp('amelioration de la fonction de comparaison de forme sinusoidale ');
disp('------------------------------------------');
beta=-(1+alpa*pi/2)
f=@(x) 1 - cos(pi*x/2/L)+beta*sin(pi*x/2/L);
D2f=@(x) (pi/2/L)^2*(cos(pi*x/2/L)-beta*sin(pi*x/2/L));
disp('Methode de Collocation');
xp=0.5;
while xp > 0
  k = -ES*D2f(xp*L) ; m = roS*f(xp*L) ; sol = sqrt(k/m);
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
k = -ES*quad(D2f,0,L) ; m = roS*quad(f,0,L); sol = sqrt(k/m);
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
disp('---------------------------------------');
disp('Methode de galerkin');
f1 = @(x) (1 - cos(pi*x/2/L)+beta*sin(pi*x/2/L)).*(pi/2/L)^2*(cos(pi*x/2/L)-beta*sin(pi*x/2/L));
f2 = @(x) (1 - cos(pi*x/2/L)+beta*sin(pi*x/2/L)).*(1 - cos(pi*x/2/L)+beta*sin(pi*x/2/L));
k = -ES*quad(f1,0,L) ; m = roS*quad(f2,0,L); sol = sqrt(k/m);
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
