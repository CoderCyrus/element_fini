function  plotdef(Sol)
% dessin de la deformee d'une structure / a la position initiale 
%
% appel plotdef(Sol)
%       en entree Sol : vecteur des deplacements nodaux dimension (nddlt)
%
%  H.Oudin
%==========================================================================
global nddln  nnod nddlt nelt nnode ndim
global Coord Connec
taille = get(0,'ScreenSize'); 
figure('Name','deformee de la structure',...
        'Position',[taille(3)/2.02 taille(4)/2.6 taille(3)/2 taille(4)/2])
axis equal
    grid on
    xlabel('X');
    ylabel('Y');
    zlabel('Z');    

XY = Coord;
if ndim == 1 
    XY(:,2) = zeros(nnod,1);XY(:,3) = zeros(nnod,1); 
    if nddln == 2 
        u(:,1) = zeros(nnod,1);u(:,2) = Sol(1:nddln:nddlt);u(:,3) = zeros(nnod,1);
    else
        u(:,1) = Sol(1:nddln:nddlt) ;u(:,2) = zeros(nnod,1); u(:,3) = zeros(nnod,1); 
    end
end
if ndim == 2 
    XY(:,3) = zeros(nnod,1); 
    u(:,1) = Sol(1:nddln:nddlt);u(:,2) = Sol(2:nddln:nddlt);u(:,3) = zeros(nnod,1);
end
if ndim == 3  
    u(:,1) = Sol(1:nddln:nddlt);u(:,2) = Sol(2:nddln:nddlt);u(:,3) = Sol(3:nddln:nddlt);
end
                 %---- calcul du facteur d'echelle et de la position deformee                    
dX = max(XY(:,1)) - min(XY(:,1)); dY = max(XY(:,2)) - min(XY(:,2));dZ = max(XY(:,3)) - min(XY(:,3));
                 %---- s : facteur d'echelle
s = max([dX dY dZ])/(20.*max(abs([u(:,1) ; u(:,2) ; u(:,3)])));
Def = XY + s * [u(:,1),u(:,2),u(:,3)];  %----- position deformee

title(['deformee de la structure avec un facteur d''echelle de ',...
        num2str(s,'%8.4f')],'Color','b')
for iel = 1:nelt
    loce=[];            %----- table de localisation pour l'element
    for i=1:nnode 
		if Connec(iel,i) > 0 loce=[loce,Connec(iel,i)]; end
	end;  
    Pos  = XY(loce,:);
    X = [[Pos(:,1)]; Pos(1,1)]; Y = [[Pos(:,2)]; Pos(1,2)];Z = [[Pos(:,3)]; Pos(1,3)];
    line(X,Y,Z,'color','g','LineStyle','--')
    Pos = Def(loce,:);    
    X = [[Pos(:,1)]; Pos(1,1)]; Y = [[Pos(:,2)]; Pos(1,2)];Z = [[Pos(:,3)]; Pos(1,3)];
    line(X,Y,Z,'color','b')  
end 
return