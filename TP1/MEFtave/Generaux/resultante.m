function [Fx,Fy,Fz] = resultante(F)
% Calcul des resultantes d'un vecteur F 
% en fonction de la dimension 1D 2D ou 3D du probleme
%
% appel [Fx,Fy,Fz] = resultante(F)
%    ou [Fx,Fy,Fz] = feval('resultante',F)
% en entree F   : vecteur de dimension nddlt
% en sortie les resltantes Fx, 0, O en dim 1
%                          Fx,Fy, 0 en dim 2
%                          Fx,Fy,Fz en dim 3
% 
%  H.Oudin 
global nddln nddlt ndim

if ndim == 3     
    Fz=sum(F(3:nddln:nddlt)); Fy=sum(F(2:nddln:nddlt)); Fx=sum(F(1:nddln:nddlt));
elseif ndim == 2 
    Fz = 0 ; Fy=sum(F(2:nddln:nddlt)); Fx=sum(F(1:nddln:nddlt));
else             
    Fz = 0 ; Fy = 0 ; Fx=sum(F(1:nddln:nddlt));
end

return