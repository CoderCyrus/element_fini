clc;clear all;close all;
% comparaison des solutions numeriques et analytique de l'exercice FV5 
% Mur climatise
%  H.Oudin 
%----------------------------------------------------------
% donnees du probleme
T0= 50 ; phie= 10; ro=-10*pi; lamda=1; e=1;
%----- solution analytique
%ue = dsolve('D2u = -ro*sin(pi*t)/lamda','u(0)=T0','Du(e)=-phie/lamda');
%ue = simple(ue)
fa=@(x) T0+x*(-phie+ro/pi)+ro*sin(pi*x)/pi/pi ;
%----------------------------------------------------------
%----- Methode de Galerkin
disp('---------------------------------------');
disp('Methode de Galerkin');
n = input('donner le degre n de l''approximation ? [2]: ');
if isempty(n) n=2; end 
disp('---------------------------------------');
taille = get(0,'ScreenSize'); 
figure('Name','comparaison Galerkin - analytique',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])            
hold on,
fplot(fa,[0 e],'r');
K=[];F=[];
for i=1:n
    for j=1:n
    F1 = @(x) i*j*x.^(i+j-2);
    K(i,j)= quad(F1,0,e);
    end
    F2 = @(x) sin(x*pi/e).*x.^i;
    F(i)= ro*quad(F2,0,e)-phie*e^i;
end
disp('---------------------------------------');
fprintf('Methode de galerkin pour un polynome de degre %3i\n',n);
disp('systeme matriciel');
disp(K); disp(F);
U = K\F' ;
y =T0;
x=0:e/100:e;
for i=1:n   y=y+U(i)*x.^i;  end
plot(x,y,'b');
legend('sol. analytique en rouge','approx par Galerkin en bleu',2);
grid
%----------------------------------------------------------
%------ Methode des elements finis
disp('---------------------------------------');
disp('Methode des elements finis');
ne = input('donner le nombre d''elements ne ? [2]: ');
if isempty(ne) ne=2; end 
figure('Name','comparaison EF - analytique',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])            
hold on, fplot(fa,[0 e],'r');
Le=e/ne;
K=[];F=[];
for i=1:ne-1         
    K(i,i)=2*lamda/Le;
    K(i,i+1)=-lamda/Le;
    K(i+1,i)=-lamda/Le;
    F1 = @(x) (sin(pi*(i/ne+x/e)) - 2*sin(pi/(2*ne))*...
    cos(pi*((i-0.5)/ne+x/e)).*x./Le);
    F(i)= ro*quad(F1,0,Le);
end
F(1)=F(1)+T0*lamda/Le;
K(ne,ne)=lamda/Le;
F2 = @(x) sin(pi*((ne-1)/ne+x/e)).*x./Le;
F(ne)= ro*quad(F2,0,Le)-phie;
disp('---------------------------------------');
fprintf('systeme matriciel pour %2d elements finis\n',ne);
    disp(K); disp(F);
U = K\F' ;
x1=0;x2=Le;
x=x1:Le/10:x2; y=T0+(U(1)-T0)*(x-x1)/Le; plot(x,y,'b');
for i=2:ne
    x1=(i-1)*Le; x2=i*Le;
    x=x1:Le/10:x2; y=U(i-1)+(U(i)-U(i-1))*(x-x1)/Le; plot(x,y,'b');
end
legend('sol. analytique en rouge','Sol. EF en bleu',2);
grid
%----------------------------------------------------------
%----- Methode de Collocation
disp('---------------------------------------');
disp('Methode de Collocation');
n = input('donner le degre n de l''approximation ? [2]: ');
if isempty(n) n=2; end 
figure('Name','comparaison Collocation - analytique',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])            
hold on,
fplot(fa,[0 e],'r');
K=[];F=[];
F(1)=-phie ;
for j=1:n  K(1,j)= j ; end
for i=2:n
    xi=e*(i-1)/n;
    for j=1:n   K(i,j)= j*(j-1)*xi^(j-2); end
    F(i)= -ro*sin(pi*xi/e)/lamda;
end
disp('---------------------------------------');
fprintf('Methode de collocation pour un polynome de degre %3i\n',n);
disp('systeme matriciel');
disp(K); disp(F);
U = K\F' ;
y =T0;
x=0:e/100:e;
for i=1:n
    y=y+U(i)*x.^i;
end
plot(x,y,'b');
legend('sol. analytique en rouge','approx par Collocation en bleu',2);
grid
%----------------------------------------------------------
%----- Methode de la valeur moyenne par sous domaine
disp('---------------------------------------');
disp('Methode de la valeur moyenne');
n = input('donner le degre n de l''approximation ? [2]: ');
if isempty(n) n=2; end 
figure('Name','comparaison Valeur moyenne - analytique',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])            
hold on,
fplot(fa,[0 e],'r');
nsd=n-1; Ld=e/nsd;
K=[];F=[];
F(1)=-phie ;
for j=1:n  K(1,j)= j; end
for i=2:n
    e1=(i-2)*Ld; e2=e1+Ld;
    for j=1:n   K(i,j)= j*(e2^(j-1)-e1^(j-1)); end
    F2 = @(x) sin(x*pi/e);
    F(i)= -ro*quad(F2,e1,e2)/lamda;
end
disp('---------------------------------------');
fprintf('Methode de la valeur moyenne pour un polynome de degre %3i\n',n);
disp('systeme matriciel');
disp(K); disp(F);
U = K\F' ;
y =T0;
x=0:e/100:e;
for i=1:n
    y=y+U(i)*x.^i;
end
plot(x,y,'b');
legend('sol. analytique en rouge','approx par valeurs moyennes en bleu',2);
grid