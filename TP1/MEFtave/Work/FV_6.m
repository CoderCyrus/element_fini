clc;clear all;close all;
% comparaison de la solution numerique et analytique de l'exercice FV6 
% Probleme de convection-diffusion (stationaire)avec stabilisation SUPG
% conditions de dirichlet nulles aux deux extremites
%  H.Oudin 
L=1;r=1;V=1; 
Pe = input('donner le nombre de peclet du probleme ? [10]: ');
if isempty(Pe) Pe=10; end   
Dk=V*L/Pe;   % diffusivite du milieu en fonction de Pe
  
nelt = input('donner le nombre d''element en x (longueur) ? [3]: ');
if isempty(nelt) nelt=3; end  
n=nelt-1;
Le=L/nelt;
K=zeros(n); F=r*Le*ones(n,1);
taille = get(0,'ScreenSize'); 
figure('Name','comparaison sol EF - analytique',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])            
hold on,
%----- solution analytique il faut utiliser la solution ue dans fa1
%ue = dsolve('-k*D2u + V*Du = r','u(0)=0','u(L)=0');
%ue = simple(ue);
fa1=@(t) r*(-exp(V/Dk*t)*L-t+t*exp(V/Dk*L)+L)/V/(-1+exp(V/Dk*L));
fplot(fa1,[0 L],'r');
%fa=@(x) x/L-(exp(Pe*x/L)-1)/(exp(Pe)-1);
%fplot(fa,[0 L],'r'), 
%----- solution numerique sans stabilisation
Peh=V*Le/(2*Dk);
for i=1:n-1          
    K(i,i)=2/Peh;
    K(i,i+1)=1-1/Peh;
    K(i+1,i)=-1-1/Peh;
end
K(n,n)=2/Peh;    %disp(K); disp(F);
K=V*K/2.;
U = K\F ;
x1=0;x2=Le;
x=x1:Le/10:x2; y=U(1)*x/Le; plot(x,y,'b');
for i=2:n
    x1=(i-1)*Le; x2=i*Le;
    x=x1:Le/10:x2; y=U(i-1)+(U(i)-U(i-1))*(x-x1)/Le; plot(x,y,'b');
end
x1=n*Le; x2=L;
x=x1:Le/10:x2; y=U(n)-U(n)*(x-x1)/Le; plot(x,y,'b');
legend('sol. analytique en rouge','Sol. EF en bleu',2);
grid
%----- solution numerique avec stabilisation SUPG


beta = input('donner la valeur du coefficient de stabilisation choisie? [nominal]: ');
if isempty(beta) beta=Le*(coth(Peh)-1/Peh)/(2*V)% valeur optimale de stabilisation 
end   
figure('Name','comparaison sol EF - analytique avec stabilisation',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])            
hold on,
fa1=@(t) r*(-exp(V/Dk*t)*L-t+t*exp(V/Dk*L)+L)/V/(-1+exp(V/Dk*L));
fplot(fa1,[0 L],'r');
x1=0;x2=Le;
x=x1:Le/10:x2; y=U(1)*x/Le; plot(x,y,'b');
for i=2:n
    x1=(i-1)*Le; x2=i*Le;
    x=x1:Le/10:x2; y=U(i-1)+(U(i)-U(i-1))*(x-x1)/Le; plot(x,y,'b');
end
x1=n*Le; x2=L;
x=x1:Le/10:x2; y=U(n)-U(n)*(x-x1)/Le; plot(x,y,'b');

Ksup=zeros(n);
for i=1:n-1          
    Ksup(i,i)=2/Le;
    Ksup(i,i+1)=-1/Le;
    Ksup(i+1,i)=-1/Le;
end
Ksup(n,n)=2/Le;    %disp(K); disp(F);
K = K + V*V*beta*Ksup;
U1 = K\F ;
x1=0;x2=Le;
x=x1:Le/10:x2; y=U1(1)*x/Le; plot(x,y,'g');
for i=2:n
    x1=(i-1)*Le; x2=i*Le;
    x=x1:Le/10:x2; y=U1(i-1)+(U1(i)-U1(i-1))*(x-x1)/Le; plot(x,y,'g');
end
x1=n*Le; x2=L;
x=x1:Le/10:x2; y=U1(n)-U1(n)*(x-x1)/Le; plot(x,y,'g');
legend('sol. analytique en rouge','Sol. EF en bleu','avec stabilisation en vert',2);
grid