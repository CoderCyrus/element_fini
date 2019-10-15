function plot2D
% trace du maillage d'une structure plane 
% avec les numero des noeuds et des elements
% les conditions aux limites (deplacements imposes) 
% et les Charges nodales
%
% appel plot2D
%
%  H.Oudin 
global nddln nnod nddlt nelt nnode ndim
global Coord Connec Ncl F

hold on
title('maillage de la structure')
axis equal

XY = Coord;
Cu = Ncl(1:nddln:nddlt); Cv = Ncl(2:nddln:nddlt);
Fx = F(1:nddln:nddlt);   Fy = F(2:nddln:nddlt);
if ndim == 1 
    XY(:,2) = zeros(nnod,1); 
    if (nddln == 1) Cv = zeros(nnod,1); Fy = zeros(nnod,1); end
end
  
for inod=1:nnod     % ----- Noeuds , conditions aux limites et charges nodales
    text(XY(inod,1),XY(inod,2),[' ',int2str(inod)],'color','m','FontSize',14)
    if Cu(inod)==1 plot(XY(inod,1),XY(inod,2),'>','color','k','MarkerSize',10);end
    if Cv(inod)==1 plot(XY(inod,1),XY(inod,2),'^','color','k','MarkerSize',10);end
    if Fx(inod)~= 0  
        text(XY(inod,1),XY(inod,2),['  \rightarrow',num2str(Fx(inod))],...
            'HorizontalAlignment','left','color','g','FontSize',14);
    end
    if Fy(inod)~= 0  
        text(XY(inod,1),XY(inod,2),['\uparrow',num2str(Fy(inod))],...
            'VerticalAlignment','bottom','color','g','FontSize',14);
    end   
end
plot(XY(:,1),XY(:,2),'.','color','m','MarkerSize',12);

for iel = 1:nelt        %----- visualisation du maillage
    loce=[];            %----- table de localisation pour l'element
    for i=1:nnode 
		if Connec(iel,i) > 0 loce=[loce,Connec(iel,i)]; end
	end;  
    Pos  = XY(loce,:);
    X = [[Pos(:,1)]; Pos(1,1)];Y = [[Pos(:,2)]; Pos(1,2)];
    line(X,Y,'color','b')
    text(mean(Pos(:,1)),mean(Pos(:,2)),[' ',int2str(iel)],'color','b','FontSize',14)
end 
return