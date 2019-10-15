clc;clear all;close all;
% Solution analytique frequence et mode de vibration d'une barre 
% Barre bi-encastree ;libre -libre ; encastree - ressort
%  H.Oudin 
%----------------------------------------------------------
% Caracteristiques de la barre
ES= 1 ; roS=1; L= 1;
%----------------------------------------------------------
k = 0;
taille = get(0,'ScreenSize'); 

while k < 5
k = menu('Choix des conditions aux limites',...
    'bi-encastree','encastree - libre','libre -libre',...
    'encastree - ressort', 'sortir') ;
tr = [1,2,3,4,5];

nmod = 3;
 switch tr(k)
  case 1, disp('Barre bi-encastree'); 
  figure('Name','solution analytique : 3 premiers modes de vibration d''une barre',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2]) 
    for imod = 1:nmod
     anal=imod*pi*sqrt(ES/roS/L/L);
     subplot(nmod,1,imod), fplot(@(x) sin(imod*pi*x/L),[0 L],'r'),...
     title(['barre bi-encastree : ',num2str(imod),'ieme mode de pulsation ',num2str(anal),...
     ' Hz ' ]), grid
    end     
  case 2, disp('encastree -libre');
  figure('Name','solution analytique : 3 premiers modes de vibration d''une barre',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2]) 
    for imod = 1:nmod
     aa=(2*imod-1)*pi/2;
     anal=aa*sqrt(ES/roS/L/L);
     subplot(nmod,1,imod), fplot(@(x) sin(aa*x/L),[0 L],'r'),...
     title(['barre encastree -libre : ',num2str(imod),'ieme mode de pulsation ',num2str(anal),...
     ' Hz ' ]), grid
    end     
  case 3, disp('libre - libre');
  figure('Name','solution analytique : 3 premiers modes de vibration d''une barre',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2]) 
    for imod = 1:nmod
     anal=imod*pi*sqrt(ES/roS/L/L);
     subplot(nmod,1,imod), fplot(@(x) cos(imod*pi*x/L),[0 L],'r'),...
     title(['barre libre -libre : ',num2str(imod),'ieme mode de pulsation ',num2str(anal),...
     ' Hz ' ]), grid
    end
  case 4, disp('encastree - ressort');
  alpa = input('donner la valeur du rapport ES/kl ? [1/2]: ');
  if isempty(alpa) alpa=0.5; end 
  figure('Name',['equation implicite tan(x)+alpa*x = 0  pour alpa =',...
      num2str(alpa)],'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2]) 
  subplot(nmod+1,1,1),hold on, fplot(@(x) tan(x),[0,3*pi,-6,6],'r'),...
           fplot(@(x) -alpa*x,[0,3*pi,-6,6],'g'),
       title('representation de l''equation implicite'), grid
 %figure('Name','solution analytique : 3 premiers modes de vibration d''une barre',...
%        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])      
  for imod = 1:nmod
     a= (2*imod-1)*pi/2+0.01;
     b=(2*imod+1)*pi/2-0.01;
     lx = fzero(@(x) tan(x)+alpa*x,[a b]);
     anal=lx*sqrt(ES/roS/L/L);
     subplot(nmod+1,1,imod+1), fplot(@(x) sin(lx*x/L),[0 L],'r'),...
     title(['barre encastree - ressort : ',num2str(imod),'ieme mode de pulsation ',num2str(anal),...
     ' Hz ' ]), grid
    end
 end
end