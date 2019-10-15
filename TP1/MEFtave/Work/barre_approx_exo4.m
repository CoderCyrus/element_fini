clc;clear all;close all;
%----- Methodes d'approximation sur la formulation faible
% comparaison des solutions numeriques et analytique de l'exercice de cours 04
% Vibration d'une barre encastree - ressort
%  H.Oudin 
%----------------------------------------------------------
% donnees du probleme
ES= 1 ; roS=1; L= 1;
alpa = input('donner la valeur du rapport ES/kl ? [10]: ');
if isempty(alpa) alpa=10; end 
%taille = get(0,'ScreenSize'); 
%figure('Name','solution analytique : 3 premiers modes de vibration de la barre',...
%        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])                
for imod = 1:3    % solution analytique
  a= (2*imod-1)*pi/2+0.01;
  b=(2*imod+1)*pi/2-0.01;
  lx = fzero(@(x) tan(x)+alpa*x,[a b]);
  anal(imod)=lx*sqrt(ES/roS/L/L);
  %subplot(3,1,imod), fplot(@(x) sin(lx*x/L),[0 L],'r'),...
  %title(['barre encastree - ressort : ',num2str(imod),...
  %' ieme mode de pulsation ',num2str(anal(imod)),' Hz ' ]), grid
end
n=2;
while n > 1
K=zeros(n);	M=zeros(n);	            
for i=1:n 
  for j=1:n  
    F1 = @(x) i*j*x.^(i+j-2);F2 = @(x) x.^(i+j);
    K(i,j)= ES*(quad(F1,0,L)+(L^(i+j-1))/alpa); 
    M(i,j)= roS*quad(F2,0,L); 
  end    
end
if n >= 4 nmod=3; else nmod =n; end
[modes,omega] = eigs(K,M,nmod,'sm'); 
omeg = diag(sqrt(omega));
disp('--------------------------------------------------');
fprintf('Methode de galerkin pour un polynome de degre %2d \n',n);
for imod = 1:nmod
  fprintf('%2d pulsation propre : Approximation : %6.3f / analytique %6.3f\n',imod,omeg(imod),anal(imod));
  fprintf('Soit une erreur de %2.1f %%\n',100*abs(omeg(1)-anal(1))/anal(1));
end
disp('---------------------------------------');
disp('pour quitter l''application ne pas donner de nouvelle valeur');
n = input('donner le degre du polynome ? [>2]: ');
if isempty(n) n=0; end 
end
close all;