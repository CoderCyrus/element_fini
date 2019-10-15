clc;clear all;close all;
% comparaison des solutions numeriques et analytique de l'exercice de cours exo9
% Vibration d'une poutre encastree - appuyee
%  H.Oudin 
%----------------------------------------------------------
% donnees du probleme
EI= 1 ; roS=1; L= 1;

% solution analytique pour le premier mode
lx1 = fzero(@(x) tan(x)-tanh(x),[pi/2+.1 3*pi/2-.1]);
omeg1=(lx1^2)*sqrt(EI/roS/(L^4));
aa = (cos(lx1)+cosh(lx1))/(sin(lx1)+sinh(lx1));
mode1 = @(x) -cos(lx1*x/L)+cosh(lx1*x/L)+aa*(sin(lx1*x/L)-sinh(lx1*x/L));
disp('---------------------------------------');
disp('Methode de Galerkin a un parametre');
%  Methode de Galerkin a un parametre (formulation forte)
%  avec les fonctions de forme w1 = x^2(L-x)(3L-2x) et w2 = x^3(L-x)(4L-3x)
w1 = @(x) ((x.^2).*(L-x)).*(3*L-2*x); 
k(1,1) = 48*EI*quad(w1,0,L);
w1w1 = @(x) (((x.^2).*(L-x)).*(3*L-2*x)).*(((x.^2).*(L-x)).*(3*L-2*x));
m11 = roS*quad(w1w1,0,L);
w1 = @(x) (((x.^2).*(L-x)).*(3*L-2*x))/sqrt(m11);
m(1,1)= m11;
sol = sqrt(k(1,1)/m11);
disp('---------------------------------------');
fprintf('Pour W1 l''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
%sol*sol
w2 = @(x) ((x.^3).*(L-x)).*(4*L-3*x);
k(2,2) = EI*quad(@(x) (((x.^3).*(L-x)).*(4*L-3*x)).*(360*x-168*L),0,L);
w2w2 = @(x) (((x.^3).*(L-x)).*(4*L-3*x)).*(((x.^3).*(L-x)).*(4*L-3*x));
m22 = roS*quad(w2w2,0,L); m(2,2) = m22;
w2 = @(x) (((x.^3).*(L-x)).*(4*L-3*x))/sqrt(m22);
sol = sqrt(k(2,2)/m22);
disp('---------------------------------------');
fprintf('Pour W2 l''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',sol,omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(sol-omeg1)/omeg1);
%sol*sol
taille = get(0,'ScreenSize'); 
figure('Name','comparaison des fonctions de forme et du premier mode',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2]) 
hold on
fplot(mode1,[0 L],'r'); fplot(w1,[0 L],'b');fplot(w2,[0 L],'g'),...
,legend('sol. analytique en rouge , fonction W1 en bleu ,et W2 en vert '), grid
pause(1);
disp('---------------------------------------');
disp('---------------------------------------');
disp('Methode de Galerkin a deux parametres');
disp('---------------------------------------');
% La solution analytique pour le second mode
lx2 = fzero(@(x) tan(x)-tanh(x),[3*pi/2+.1 5*pi/2-.1]);
omeg2=(lx2^2)*sqrt(EI/roS/(L^4));
aa = (cos(lx2)+cosh(lx2))/(sin(lx2)+sinh(lx2));
mode2 = @(x) -cos(lx2*x/L)+cosh(lx2*x/L)+aa*(sin(lx2*x/L)-sinh(lx2*x/L));
%----- Méthode de Galerkin a 2 parametres
k(1,2) = EI*quad(@(x) (((x.^2).*(L-x)).*(3*L-2*x)).*(360*x-168*L),0,L);
k(2,1) = 48*EI*quad(@(x) ((x.^3).*(L-x)).*(4*L-3*x),0,L);
%disp('matrice K');k
w1w2 = @(x) (((x.^2).*(L-x)).*(3*L-2*x)).*(((x.^3).*(L-x)).*(4*L-3*x));
m(1,2) = roS*quad(w1w2,0,L); m(2,1)= m(1,2);
%disp('matrice M'); m
[modes,omega] = eigs(k,m,2,'sm'); 
omeg = diag(sqrt(omega));
fprintf('L''approximation de la premiere pulsation propre est %6.3f / analytique %6.3f\n',omeg(1),omeg1);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(omeg(1)-omeg1)/omeg1);
fprintf('L''approximation de la seconde pulsation propre est %6.3f / analytique %6.3f\n',omeg(2),omeg2);
fprintf('Soit une erreur de %2.1f %%\n',100*abs(omeg(2)-omeg2)/omeg2);
Z1 = @(x) -5.58*modes(1,1)*(((x.^2).*(L-x)).*(3*L-2*x))-...
5.58*modes(2,1)*((x.^3).*(L-x)).*(4*L-3*x);
Z2 = @(x) -33.5*modes(1,2)*(((x.^2).*(L-x)).*(3*L-2*x))-...
33.5*modes(2,2)*((x.^3).*(L-x)).*(4*L-3*x);
figure('Name','comparaison Analytique / Galerkin pour les deux modes ',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2]) 
hold on
fplot(mode1,[0 L],'r');fplot(mode2,[0 L],'r'); fplot(Z1,[0 L],'b');fplot(Z2,[0 L],'b'),...
,legend('sol. analytique en rouge , Approximation en bleu '), grid