procedure NUTAC1(T: real);
var DJ,D,S,Z,U,CS,CZ,CU,SS,SZ,SU,SEC,DP,DE,DS: real;
    I,J,K: integer;
    JD: longint;
    NUT,TEM: MS33;
{$I NUTEQU}
begin
  SEC:=pi/180.0/3600.0;
  JD:=JD0-2451545; DJ:=(T0+T)/86.4-0.5;
  D:=(JD+DJ)/36525.0;
  Z:=84381.448+D*(-46.815+D*(-0.00059+D*0.001813)); Z:=Z*SEC;
  CZ:=cos(Z); SZ:=sin(Z);
  NUTEQU(D,DE,DP);
  CU:=cos(DP); SU:=sin(DP);
  S:=Z+DE; CS:=cos(S); SS:=sin(S);
  {writeln('dpsi,de ',DP,' ',DE);}
  NUT[1,1]:=CU;    NUT[1,2]:=-SU*CZ;         NUT[1,3]:=-SU*SZ;
  NUT[2,1]:=CS*SU; NUT[2,2]:=CS*CU*CZ+SS*SZ; NUT[2,3]:=CS*CU*SZ-SS*CZ;
  NUT[3,1]:=SS*SU; NUT[3,2]:=SS*CU*CZ-CS*SZ; NUT[3,3]:=SS*CU*SZ+CS*CZ;
  Z:=DP*CZ;
  CZ:=cos(Z); SZ:=sin(Z);
  TEM[1,1]:= CZ; TEM[1,2]:=SZ; TEM[1,3]:=0;
  TEM[2,1]:=-SZ; TEM[2,2]:=CZ; TEM[2,3]:=0;
  TEM[3,1]:= 0;  TEM[3,2]:=0;  TEM[3,3]:=1;
  for I:=1 to 3 do for J:=1 to 3 do
  begin
    MATR[I,J]:=0;
    for K:=1 to 3 do MATR[I,J]:=MATR[I,J]+TEM[I,K]*NUT[K,J]
  end
end;
