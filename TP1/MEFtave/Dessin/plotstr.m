function plotstr
% trace du maillage d'une structure avec les numeros des noeuds et des elements 
% et les conditions aux limites (Charges nodales, deplacements imposes) 
%
% appel plotstr
%  H.Oudin
%==========================================================================
global nddln nnod nddlt nelt nnode ndim
global Coord Connec Ncl F

taille = get(0,'ScreenSize'); 
figure('Name','maillage de la structure',...
        'Position',[taille(3)/2.01 taille(4)/2.4 taille(3)/2 taille(4)/2])
hold on
title('maillage et CL de la structure etudiee')
axis equal
    grid on
    xlabel('X');
    ylabel('Y');
    zlabel('Z');    
XY = Coord;
Cu = Ncl(1:nddln:nddlt); Cv = Ncl(2:nddln:nddlt); Cw = Ncl(3:nddln:nddlt);
Fx = F(1:nddln:nddlt);   Fy = F(2:nddln:nddlt); Fz = F(3:nddln:nddlt);
if ndim == 1 
    XY(:,2) = zeros(nnod,1);XY(:,3) = zeros(nnod,1); 
    Cv = zeros(nnod,1); Fy = zeros(nnod,1);Cw = zeros(nnod,1); Fz = zeros(nnod,1);
end
if ndim == 2 
    XY(:,3) = zeros(nnod,1); 
    Cw = zeros(nnod,1); Fz = zeros(nnod,1);
end 


for inod=1:nnod     % ----- Noeuds , conditions aux limites et charges nodales
    text(XY(inod,1),XY(inod,2),XY(inod,3),[' ',int2str(inod)],'color','m','FontSize',14)
    if Cu(inod)==1 plot3(XY(inod,1),XY(inod,2),XY(inod,3),'>','color','k','MarkerSize',10);end
    if Cv(inod)==1 plot3(XY(inod,1),XY(inod,2),XY(inod,3),'^','color','k','MarkerSize',10);end
    if Cw(inod)==1 plot3(XY(inod,1),XY(inod,2),XY(inod,3),'*','color','k','MarkerSize',10);end
    if Fx(inod)~= 0  
        text(XY(inod,1),XY(inod,2),XY(inod,3),['  \rightarrow',num2str(Fx(inod))],...
            'HorizontalAlignment','left','color','g','FontSize',14);
    end
    if Fy(inod)~= 0  
        %text(XY(inod,1),XY(inod,2),XY(inod,3),['\uparrow',num2str(Fy(inod))],...
         %   'VerticalAlignment','bottom','color','g','FontSize',14);
            text(XY(inod,1),XY(inod,2),XY(inod,3),['\uparrow',num2str(Fy(inod))],...
            'VerticalAlignment','bottom','color','g','FontSize',14);
    end   
    if Fz(inod)~= 0  
        text(XY(inod,1),XY(inod,2),XY(inod,3),['\uparrow',num2str(Fz(inod))],...
            'VerticalAlignment','bottom','color','g','FontSize',14);
    end   

end
plot3(XY(:,1),XY(:,2),XY(:,3),'.','color','m','MarkerSize',12);
for iel = 1:nelt        %----- visualisation du maillage
    loce=[];            %----- table de localisation pour l'element
    for i=1:nnode 
		if Connec(iel,i) > 0 loce=[loce,Connec(iel,i)]; end
	end;  
    Pos  = XY(loce,:);
    X = [[Pos(:,1)]; Pos(1,1)];Y = [[Pos(:,2)]; Pos(1,2)];Z = [[Pos(:,3)]; Pos(1,3)];
    line(X,Y,Z,'color','b')
    text(mean(Pos(:,1)),mean(Pos(:,2)),mean(Pos(:,3)),[' ',int2str(iel)],'color','b','FontSize',14)
end 
return