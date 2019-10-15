function [Ke,Me] = poutre_keme(iel)
% Calcul des matrices elementaire Ke et Me de l'element (iel) d'un portique
% 
% appel [Ke,Me] = poutre_keme(iel)
%    ou [Ke,Me] = feval('poutre_keme',iel)
% en entree iel : numero de l'element
% en sortie Ke  : matrice raideur elementaire (4,4)en 1D; (6,6) en 2D 
%           Me  : matrice masse elementaire 
%
%  H.Oudin 
global ndim
global Coord Connec Nprop Prop 

X  = Coord(Connec(iel,:),:);
dX = X(2,:) - X(1,:);  

if ndim == 1  
    EI=Prop(Nprop(iel),2);RS=Prop(Nprop(iel),3);
    L = abs(dX);    %----- matrice raideur elementaire sur (v,teta)
    Ke = (EI/L^3)*[  12   6*L   -12   6*L
                     6*L  4*L^2 -6*L  2*L^2
                    -12  -6*L    12  -6*L
                     6*L  2*L^2 -6*L  4*L^2
                   ]; 
                   
     Me = (RS*L)*[  13/35     11*L/210    9/70   -13*L/420
                    11*L/210   L*L/105  13*L/420  -L*L/140 
                    9/70      13*L/420    13/35   -11*L/210
                  -13*L/420   -L*L/140  -11*L/210  L*L/105
               ];
elseif ndim == 2 
 ES=Prop(Nprop(iel),1); EI=Prop(Nprop(iel),2);RS=Prop(Nprop(iel),3);
 L = sqrt(dX(1)^2 + dX(2)^2);
 c = dX(1)/L; s = dX(2)/L;  
 K = (ES/L)*[  1  0  0 -1  0  0
               0  0  0  0  0  0
               0  0  0  0  0  0
              -1  0  0  1  0  0
               0  0  0  0  0  0
               0  0  0  0  0  0
            ];
 K = K + (EI/L^3)*[  0  0    0     0   0    0
                     0  12   6*L   0  -12   6*L
                     0  6*L  4*L^2 0  -6*L  2*L^2
                     0  0    0     0   0    0
                     0 -12  -6*L   0   12  -6*L
                     0  6*L  2*L^2 0  -6*L  4*L^2
                   ]; 
  P = [  c  s  0  0  0  0
        -s  c  0  0  0  0
         0  0  1  0  0  0
         0  0  0  c  s  0
         0  0  0 -s  c  0
         0  0  0  0  0  1
            ];
  Ke =P' * K * P;    %----- matrice raideur elementaire sur (u,v,teta)

  M = (RS*L)*[  1/3      0            0      1/6     0         0
                 0     13/35     11*L/210     0     9/70    -13*L/420
                 0    11*L/210      L*L/105   0    13*L/420  -L*L/140 
                1/6      0             0     1/3     0          0
                 0      9/70      13*L/420    0     13/35    -11*L/210
                 0   -13*L/420     -L*L/140   0    -11*L/210   L*L/105
               ];
  Me =P' * M * P;    %----- matrice massse elementaire sur (u,v,teta)
  
elseif ndim == 3
     disp('================================================ ');
     disp('       element non programme ');
     disp('================================================ ');
end
return