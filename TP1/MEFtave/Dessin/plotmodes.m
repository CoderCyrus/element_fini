function  plotmodes(Sol,k,f)
% Trace des modes de vibrations d'un portique 2D  
%
% appel plotmodes(Sol,k,f)
% en entree Sol : vecteur des deplacements nodaux dimension (nddlt)
%           k   : numero du mode
%           f   : valeur de la frequence
%
%  H.Oudin
%==========================================================================
global nddln nddlt nelt nnode
global Coord Connec

axis equal
                 %---- calcul du facteur d'echelle et de la position deformee                    
dX = max(Coord(:,1)) - min(Coord(:,1)); dY = max(Coord(:,2)) - min(Coord(:,2));
%s = max([dX dY])/(10.*max(abs(Sol))) ;     %---- facteur d'echelle 
s = max([dX dY])/(20.*max(abs([Sol(1:nddln:nddlt) ; Sol(2:nddln:nddlt)])));
Def = Coord + s * [Sol(1:nddln:nddlt),Sol(2:nddln:nddlt)];  %----- postion deformee

title([int2str(k),' mode de vibration a la frequence ',...
        num2str(f,'%7.2f')],'Color','b')

for iel = 1:nelt
    loce=[];            %----- table de localisation pour l'element
    for i=1:nnode 
		if Connec(iel,i) > 0 loce=[loce,Connec(iel,i)]; end
	end;  
    Pos  = Coord(loce,:);
    X = [[Pos(:,1)]; Pos(1,1)]; Y = [[Pos(:,2)]; Pos(1,2)];
    line(X,Y,'color','g','LineStyle','--')
    Pos = Def(loce,:);    
    X = [[Pos(:,1)]; Pos(1,1)]; Y = [[Pos(:,2)]; Pos(1,2)];
    line(X,Y,'color','b')   
end 
return