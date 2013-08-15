procedure SUNMOON(T: real; var RS,RM: real; var S,M: MS3);
const AE=1.4959787e5; RE=6.37814; EP0=23.439291111;
      EP1=-46.815/3600.0; EP2=-0.00059/3600.0; EP3=0.001813/3600.0;
var DJ,TJ,L,B,C,D,Y,Z,EP,SE,CE,GRAD: real;
    I,JD: longint;
{$I SUN200}
{$I MOON}
begin
  GRAD:=pi/180;
  JD:=JD0-2451545; DJ:=(T0+T)/86.4-0.5;
  I:=trunc(DJ); JD:=JD+I; DJ:=DJ-I;
  TJ:=(JD+DJ)/36525.0;
  SUN200(TJ,L,B,D);
  L:=L*GRAD; B:=B*GRAD; RS:=D*AE;
  C:=cos(B); S[1]:=C*cos(L); Y:=C*sin(L); Z:=sin(B);
  EP:=GRAD*(EP0+TJ*(EP1+TJ*(EP2+TJ*EP3)));
  CE:=cos(EP); SE:=sin(EP);
  S[2]:=Y*CE-Z*SE; S[3]:=Y*SE+Z*CE;
  MOON(TJ,L,B,D);
  L:=L*GRAD; B:=B*GRAD; RM:=D*RE;
  C:=RM*cos(B);
  M[1]:=C*cos(L); Y:=C*sin(L); Z:=RM*sin(B);
  M[2]:=Y*CE-Z*SE; M[3]:=Y*SE+Z*CE
end;
