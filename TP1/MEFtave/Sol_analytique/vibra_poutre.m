clc;clear all;close all;
% Solution analytique frequence et mode de vibration d'une poutre 
%  H.Oudin 
%----------------------------------------------------------
% Caracteristiques de la poutre
EI= 1 ; roS=1; L= 1;
%----------------------------------------------------------
k = 0;
taille = get(0,'ScreenSize'); 

while k < 6
k = menu('Choix des conditions aux limites',...
    'libre -libre','bi-appuyee','bi-encastree',...
    'appuyee - libre','encastree - appuyee', 'sortir') ;
tr = [1,2,3,4,5,6];
nmod = 3;
switch tr(k)
case 1, disp('poutre libre - libre');
  disp('c''est l''exemple traite dans le cours');
  disp('il existe deux modes rigide une translation une rotation');
  figure('Name','3 premiers modes de vibration d''une poutre libre - libre',...
  'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2]) 
  subplot(nmod+1,1,1),hold on, fplot('[cos(x) , 1/cosh(x)]',[0,5*pi,-1,1],'b'),...
  title('representation de l''equation implicite cos(x)cosh(x)-1 = 0'), grid   
  for imod = 1:nmod
     a= imod*pi;
     b=(imod+1)*pi;
     lx = fzero(@(x) cos(x)*cosh(x)-1,[a b])
     anal=(lx^2)*sqrt(EI/roS/(L^4));
     aa = -(cos(lx)-cosh(lx))/(sin(lx)-sinh(lx));
     mode = @(x) cos(lx*x/L)+cosh(lx*x/L)+aa*(sin(lx*x/L)+sinh(lx*x/L));
     %mode2= @(x) (cos(lx*x/L)+cosh(lx*x/L)+aa*(sin(lx*x/L)+sinh(lx*x/L))).*...
     %(cos(lx*x/L)+cosh(lx*x/L)+aa*(sin(lx*x/L)+sinh(lx*x/L)));
     %m = roS*quad(mode2,0,L);  % Les modes sont M norme  m=1 
     subplot(nmod+1,1,imod+1), fplot(mode,[0 L],'r'),...
     title(['poutre libre - libre ',num2str(imod),'ieme mode de pulsation ',num2str(anal),...
     ' Hz ' ]), grid
   end  
case 2, disp('poutre bi-appuyee les solutions sont lx = i*pi'); 
    figure('Name','3 premiers modes de vibration d''une poutre bi-appuyee',...
       'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2]) 
    for imod = 1:nmod
     lx = imod*pi
     anal=(lx^2)*sqrt(EI/roS/(L^4));
     subplot(nmod,1,imod), fplot(@(x) sin(imod*pi*x/L),[0 L],'r'),...
     title(['poutre bi-appuyee ',num2str(imod),'ieme mode de pulsation ',num2str(anal),...
     ' Hz ' ]), grid
    end     
case 3, disp('poutre bi-encastree');
  figure('Name','3 premiers modes de vibration d''une poutre bi-encastree',...
  'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2]) 
  subplot(nmod+1,1,1),hold on, fplot('[cos(x) , 1/cosh(x)]',[0,5*pi,-1,1],'b'),...
  title('representation de l''equation implicite cos(x)cosh(x)-1 = 0'), grid    
  for imod = 1:nmod
     a= imod*pi;
     b=(imod+1)*pi;
     lx = fzero(@(x) cos(x)*cosh(x)-1,[a b])
     anal=(lx^2)*sqrt(EI/roS/(L^4));
     aa = -(cos(lx)-cosh(lx))/(sin(lx)-sinh(lx));
     mode = @(x) cos(lx*x/L)-cosh(lx*x/L)+aa*(sin(lx*x/L)-sinh(lx*x/L));
     subplot(nmod+1,1,imod+1), fplot(mode,[0 L],'r'),...
     title(['poutre bi - encastree ',num2str(imod),'ieme mode de pulsation ',num2str(anal),...
     ' Hz ' ]), grid
  end  
case 4, disp('poutre appuyee - libre');
  disp('il existe un mode rigide de rotation');
  figure('Name','3 premiers modes de vibration d''une poutre appuyee - libre',...
  'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2]) 
  subplot(nmod+1,1,1),hold on, fplot('[tan(x) , tanh(x)]',[0,4*pi,-1.5,1.5],'b'),...
  title('representation de l''equation implicite tan(x) = th(x)'), grid   
  for imod = 1:nmod
     a= imod*pi/2+.1;
     b=(2*imod+1)*pi/2-.1;
     lx = fzero(@(x) tan(x)-tanh(x),[a b])
     anal=(lx^2)*sqrt(EI/roS/(L^4));
     aa = sin(lx)/sinh(lx);
     mode = @(x) sin(lx*x/L)+aa*sinh(lx*x/L);
     subplot(nmod+1,1,imod+1), fplot(mode,[0 L],'r'),...
     title(['poutre appuyee - libre ',num2str(imod),'ieme mode de pulsation ',num2str(anal),...
     ' Hz ' ]), grid
   end
 case 5, disp('poutre encastree - appuyee');
  figure('Name','3 premiers modes de vibration d''une poutre encastree - appuyee',...
  'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2]) 
  subplot(nmod+1,1,1),hold on, fplot('[tan(x) , tanh(x)]',[0,4*pi,-1.5,1.5],'b'),...
  title('representation de l''equation implicite tan(x) = th(x)'), grid   
  for imod = 1:nmod
     a= imod*pi/2+.1;
     b=(2*imod+1)*pi/2-.1;
     lx = fzero(@(x) tan(x)-tanh(x),[a b])
     anal=(lx^2)*sqrt(EI/roS/(L^4));
     aa = (cos(lx)+cosh(lx))/(sin(lx)+sinh(lx));
     mode = @(x) cos(lx*x/L)-cosh(lx*x/L)-aa*(sin(lx*x/L)-sinh(lx*x/L));
     subplot(nmod+1,1,imod+1), fplot(mode,[0 L],'r'),...
     title(['poutre encastree - appuyee ',num2str(imod),'ieme mode de pulsation ',num2str(anal),...
     ' Hz ' ]), grid
   end    
end
end
disp('les autres cas seront simples a programmer en utilisant le cours');