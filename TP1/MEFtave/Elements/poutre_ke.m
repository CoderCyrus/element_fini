function [Ke,Fe] = poutre_ke(iel)
% Calcul de la matrice raideur Ke et de la force generalisee Fe 
% pour un element (iel) d'une structure portique
% 
% appel [Ke,Fe] = poutre_ke(iel)
%    ou [Ke,Fe] = feval('poutre_ke',iel)
% en entree iel : numero de l'element
% en sortie Ke  : matrice raideur elementaire 
%           Fe  : force generalisee elementaire  
%
%  H.Oudin & D.Motte (EF 3D)
global ndim
global Coord Connec Nprop Prop 

X  = Coord(Connec(iel,:),:);
dX = X(2,:) - X(1,:);  

if ndim == 1  
    EI=Prop(Nprop(iel),1);
    L = abs(dX);    %----- matrice raideur elementaire sur (v,teta)
    Ke = (EI/L^3)*[  12   6*L   -12   6*L
                     6*L  4*L^2 -6*L  2*L^2
                    -12  -6*L    12  -6*L
                     6*L  2*L^2 -6*L  4*L^2
                   ]; 
                   
    f=Prop(Nprop(iel),2);
    Fe = (f*L/2)*[1; L/6; 1; -L/6];
elseif ndim == 2 
 ES=Prop(Nprop(iel),1); EI=Prop(Nprop(iel),2);
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

  f  = [ c*Prop(Nprop(iel),3)+s*Prop(Nprop(iel),4)
        -s*Prop(Nprop(iel),3)+c*Prop(Nprop(iel),4)];
  fe = (f(1)*L/2)*[1; 0; 0; 1; 0; 0];
  fe = fe + (f(2)*L/2)*[0; 1; L/6; 0; 1; -L/6];
  Fe = P' * fe;     %----- force generalisee elementaire sur (u,v,teta)
  
elseif ndim == 3
%-------------------------------------------------------------------------------------------------
ES = Prop(Nprop(iel),1);				% Matrice de raideur locale
GJ = Prop(Nprop(iel),2);
EIy = Prop(Nprop(iel),3);
EIz = Prop(Nprop(iel),4);
L = sqrt(dX(1)^2+dX(2)^2+dX(3)^2);

K = [ ES/L	0			0			0		0			0			-ES/L	0			0			0		0			0
      0		12*EIz/L^3	0			0		0			6*EIz/L^2	0		-12*EIz/L^3	0			0		0			6*EIz/L^2
	  0		0			12*EIy/L^3	0		-6*EIy/L^2	0			0		0			-12*EIy/L^3	0		-6*EIy/L^2	0
	  0		0			0			GJ/L	0			0			0		0			0			-GJ/L	0			0
	  0		0			-6*EIy/L^2	0		4*EIy/L		0			0		0			6*EIy/L^2	0		2*EIy/L		0
      0		6*EIz/L^2	0			0		0			4*EIz/L		0		-6*EIz/L^2	0			0		0			2*EIz/L
	  -ES/L	0			0			0		0			0			ES/L	0			0			0		0			0
      0		-12*EIz/L^3	0			0		0			-6*EIz/L^2	0		12*EIz/L^3	0			0		0			-6*EIz/L^2
	  0		0			-12*EIy/L^3	0		6*EIy/L^2	0			0		0			12*EIy/L^3	0		6*EIy/L^2	0
	  0		0			0			-GJ/L	0			0			0		0			0			GJ/L	0			0
	  0		0			-6*EIy/L^2	0		2*EIy/L		0			0		0			6*EIy/L^2	0		4*EIy/L		0
      0		6*EIz/L^2	0			0		0			2*EIz/L		0		-6*EIz/L^2	0			0		0			4*EIz/L
	];

tx=dX(1)/L;								% parametres du repere local / repere global
ty=dX(2)/L;
tz=dX(3)/L;
t=[tx ; ty ; tz];
prodzt=cross([0,0,1],t);
tn=prodzt/sqrt(prodzt(1)^2+prodzt(2)^2+prodzt(3)^2);
tm=cross(t,tn);
Q=[ t(1)	,	t(2)	,	t(3)		% Matrice de passage 3x3
	tn(1)	,	tn(2)	,	tn(3)
	tm(1)	,	tm(2)	,	tm(3)
  ];
nul=zeros(3);
P=[	Q	,	nul	,	nul	,	nul			% Matrice de passage 12x12
	nul	,	Q	,	nul	,	nul
	nul	,	nul	,	Q	,	nul
	nul	,	nul	,	nul	,	Q
  ];

Ke = P' * K * P;		%----- matrice raideur elementaire sur (u,v,w,teta1,teta2,teta3)

Fe = zeros(12,1);		%----- pas de force generalisee elementaire
%-------------------------------------------------------------------------------------------------
end
return