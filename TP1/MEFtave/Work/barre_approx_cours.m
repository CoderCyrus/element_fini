clc;clear all;close all;
% comparaison des solutions numeriques et analytique de l'exercice de cours 
% Vibration d'une barre encastree - libre
%  H.Oudin 
%----------------------------------------------------------
% donnees du probleme
ES= 1 ; roS=1; L= 1;
%----------------------------------------------------------
%----- Methodes d'approximation sur la formulation forte
%----- avec la fonction de forme x(2x-L)
disp('---------------------------------------');
disp('Vibration d''une barre encastree - libre');
disp('---------------------------------------');
disp('Methodes d''approximation avec la fonction de forme x(2x-L)');
disp('Utilise la formulation forte du probleme');
omeg1=pi*sqrt(ES/roS/L/L)/2;    % solution analytique
taille = get(0,'ScreenSize'); 
figure('Name','comparaison de la fonction de forme et du premier mode',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])            
fplot(@(x)[sin(pi*x/2/L) , x.*(2*L-x)],[0 L]);
%Collocation
f=@(x) x.*(2*L-x); % D2f=@(x) -2*x./x;
disp('---------------------------------------');
disp('Methode de Collocation');
xp=0.5;
while xp > 0
  k = 2*ES ; m = roS*f(xp) ; sol = sqrt(k/m);
  disp('---------------------------------------');
  fprintf('Pour un point de collocation situe en %3.1f L \n',xp);
  fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
  fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
  disp('---------------------------------------');
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
f2 = @(x) (x.*(2*L-x)).^2;
k = 2*ES*quad(f,0,L); m = roS*quad(f2,0,L); sol = sqrt(k/m);
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
close all;
%----------------------------------------------------------
%----- Méthode de Galerkin sur la formulation faible pour un polynome de degre n
disp('---------------------------------------------------------------------------');
disp('---------------------------------------------------------------------------');
n=2;
fprintf('Methode de Galerkin sur la formulation faible pour un polynome de degre %3i\n',n);
while n > 1
K=zeros(n);	M=zeros(n);	
for i=1:n
    for j=1:n
     F1 = @(x) i*j*x.^(i+j-2); K(i,j)= ES*quad(F1,0,L);
     F2 = @(x) x.^(i+j); M(i,j)= roS*quad(F2,0,L);
    end
end
%disp('systeme matriciel');K M
if n >= 4 nmod=3; else nmod =n; end
[modes,omega] = eigs(K,M,n,'sm'); 
omeg = diag(sqrt(omega));
%modes
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',omeg(1),omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(omeg(1)-omeg1)/omeg1);

figure('Name','comparaison Galerkin - analytique',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2]) 

for imod = 1:nmod
  U= modes(:,imod)  ; %
  %if U(1)<0  U = -U ;end 
  y =0;   x=0:L/100:L;
  for i=1:n  y=y+U(i)*x.^i;  end
  maxi = max(abs(y)); y=y/maxi;
  subplot(nmod,1,imod),hold on, plot(x,y,'b'),plot(x,-y,'b')
  sol = omeg(imod); anal=(2*imod-1)*pi*sqrt(ES/roS/L/L)/2;
  err = 100*abs(sol-anal)/anal;
  subplot(nmod,1,imod), fplot(@(x) -sin((2*imod-1)*pi*x/L/2),[0 L],'r'),
  subplot(nmod,1,imod), fplot(@(x) sin((2*imod-1)*pi*x/L/2),[0 L],'r'),...
  title([num2str(imod),'ieme mode : pulsation approchee ',num2str(sol),...
  ' / analytique ',num2str(anal),' soit une erreur de ',num2str(err),...
  ' % ' ]), legend('sol. analytique en rouge','Sol. approchee en bleu'), grid
end
disp('---------------------------------------');
disp('pour quitter l''application ne pas donner de nouvelle valeur');
n = input('donner le degre n de l''approximation ? [2]: ');
if isempty(n) n=0; end 
end
close all;