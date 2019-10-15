function [critere1 , critere2] =  verification(Nbis,Nprop,E,Prop,Re,Connec,Coord)
%critere1
s=zeros(1,length(Nbis));
length(Nbis)
for i= 1:length(Nbis)
    s(i)=Prop(Nprop(i),1) ./ E(i,1);
end
sigma=Nbis./s;
critere1 = sigma < Re;

%critere2
%longueur de Lc
  tab_Lc=zeros(length(Connec),1)
  for i = 1:length(Connec)  
  Lccarre =(Coord(Connec(i,1),1)-Coord(Connec(i,2),1))^2 +(Coord(Connec(i,1),2)-Coord(Connec(i,2),2))^2;
  tab_Lc(i,1) = sqrt (Lccarre)
  end
  tab_Lc
 
  % Force critique d'euler

  I = s.^2 /4*pi
  I
  E
  E'
  pi*pi*E'.*I 
  tab_Lc'.^2

 F=pi*pi*E'.*I ./ (tab_Lc'.^2)
critere2 = abs(Nbis) < F;
F
Nbis

end