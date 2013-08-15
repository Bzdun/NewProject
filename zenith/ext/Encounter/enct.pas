program ENCT;
uses graph, GRAPHICS;
const SC=3.14159265359/1.8e2/60;
type real=double;
     MS3=array[1..3] of real;
     MS5=array[1..5] of real;
     MS6=array[1..6] of real;
     MASN=array[1..84] of real;
     MS33=array[1..3,1..3] of real;
     MS66=array[1..6,1..6] of real;
     MASN9=array[1..9] of MASN;
     INTPR=procedure(T: real; var Y,G: MASN);
     INTOU= procedure(T,H: real; var Z: MASN; var P: MS5; var F: MASN9);

var VARI: boolean;
var CM0,CM1,A,B: MS66;
    Y,Y0: MASN;
    FQ: MASN9;
    P: MS5;
    T0,RSOL,RMOO,SZ0,W,BK,BK0,BK1,S,DL,D,XN,XV,
    SMIN,HMIN,U1,U2,H,S0: real;
    JD0: longint;
    I,J,K: integer;
    N,NN: longint;
    SOL,MOO: MS3;
    MATR: MS33;
    Z: MS6;
    LV,EV,EN: MS3;
    FL1,FL2,FL_KA1, FL_KA2: text;
    inputKA1,inputKA2,outputKA1,outputKA2: string;

{$I HARMON36}
{$I ROATMST}
{$I NUTAC1}
{$I SUNMOON}
{���� �ࠢ�� ��� ����. ��-��� �������� ��
Y[1..3] - ࠤ. �����, Y[4..6] - ᪮����, dY/dT=G(T,Y)}
procedure PRAV0(T: real; var Y,G: MASN);
const GC=398.60044; R0=6.378136; J2EARTH=1082625.7E-9;
      GCS=332946.0*GC; GCM=GC/81.30068; NGAM=16;
var A,B,C,DA,DC,DH,DR,H,R,R2,X1,X2,X3,V1,V2,V3,Y12,Y3,RV,
    SZ,SA,CA,VX,VY,DRM,RS,RM,SR,JDATE: real;
    X,F,LN,SL: MS3;
    I,J: integer;
begin
  G[1]:=Y[4]; G[2]:=Y[5]; G[3]:=Y[6];
  Y12:=sqr(Y[1])+sqr(Y[2]); Y3:=sqr(Y[3]);
  R2:=Y12+Y3; R:=sqrt(R2); H:=-GC/(R*R2);
  {Earth gravitational field}
  SZ:=SZ0+W*T; SA:=sin(SZ); CA:=cos(SZ);
  X[1]:=Y[1]*CA+Y[2]*SA; X[2]:=-Y[1]*SA+Y[2]*CA; X[3]:=Y[3];
  HARMO(NGAM,GC,R0,A,X,F);
  X[1]:=F[1]*CA-F[2]*SA; X[2]:=F[1]*SA+F[2]*CA; X[3]:=F[3];
  G[4]:=H*Y[1]+X[1];
  G[5]:=H*Y[2]+X[2];
  G[6]:=H*Y[3]+X[3];
  {atmosphere}
  RV:=1e6*ROATMST(Y[1],Y[2],Y[3]);
  {
  JDATE:=(T0+T)/86.4-0.5+JD0;
  RV:=1e6*ROA01(JDATE,F10_7,F81,KP,Y[1],Y[2],Y[3],false);
  }
  VX:=Y[4]+W*Y[2]; VY:=Y[5]-W*Y[1];
  RV:=RV*sqrt(sqr(VX)+sqr(VY)+sqr(Y[6]));
  A:=BK*RV;
  G[4]:=G[4]-A*VX; G[5]:=G[5]-A*VY; G[6]:=G[6]-A*Y[6];
  {Sun and Moon}
  SL:=SOL; LN:=MOO; RS:=RSOL; RM:=RMOO;
  DRM:=sqr(LN[1]-Y[1])+sqr(LN[2]-Y[2])+sqr(LN[3]-Y[3]);
  SR:=3*(SL[1]*Y[1]+SL[2]*Y[2]+SL[3]*Y[3]);
  DRM:=GCM/(DRM*sqrt(DRM));
  RM:=GCM/(RM*sqr(RM));
  RS:=GCS/(RS*sqr(RS));
  G[4]:=G[4]+DRM*(LN[1]-Y[1])-RM*LN[1]+RS*(SR*SL[1]-Y[1]);
  G[5]:=G[5]+DRM*(LN[2]-Y[2])-RM*LN[2]+RS*(SR*SL[2]-Y[2]);
  G[6]:=G[6]+DRM*(LN[3]-Y[3])-RM*LN[3]+RS*(SR*SL[3]-Y[3]);
  if not VARI then exit;
  {Variation equations, only c_20 is taken into account}
  B:=1.5*J2EARTH*sqr(R0/R2);
  A:=H*(1+B*(Y12-4*Y3)); C:=H*(1+B*(3*Y12-2*Y3));
  for J:=1 to 6 do
  begin
    I:=6*J;
    X1:=Y[I+1]; X2:=Y[I+2]; X3:=Y[I+3];
    V1:=Y[I+4]; V2:=Y[I+5]; V3:=Y[I+6];
    DR:=(Y[1]*X1+Y[2]*X2+Y[3]*X3)/R2;
    DH:=-3*H*DR;
    DA:=DH-7*(A-H)*DR+2*H*B*(Y[1]*X1+Y[2]*X2-4*Y[3]*X3);
    DC:=DH-7*(C-H)*DR+2*H*B*(3*(Y[1]*X1+Y[2]*X2)-2*Y[3]*X3);
    G[I+1]:=V1; G[I+2]:=V2; G[I+3]:=V3;
    G[I+4]:=A*X1+Y[1]*DA;
    G[I+5]:=A*X2+Y[2]*DA;
    G[I+6]:=C*X3+Y[3]*DC
  end
end;

{���� �ࠢ�� ��� ����. ��-��� �������� ���� ��;
Y[1..6] - 䠧��� ����� ��-����⥫�,
Y[7..12] - 䠧��� ����� ��-楫�}
procedure PRAVC1(T: real; var Y,G: MASN);
var Z,V: MASN;
    I: integer;
begin
  VARI:=false;
  BK:=BK0;
  for I:=1 to 6 do Z[I]:=Y[I];
  PRAV0(T,Z,V);
  for I:=1 to 6 do begin G[I]:=V[I]; Z[I]:=Y[I+6] end;
  BK:=BK1;
  PRAV0(T,Z,V);
  for I:=1 to 6 do G[I+6]:=V[I]
end;

{���� �ࠢ�� ��� ����. ��-��� �������� ���� �� � ᮮ⢥�������
���室��� �����; Y[1..42] - 䠧�� ����� � ���室��� 6*6-����� ��
�⮫�栬 ��-����⥫�, Y[43..84] - 䠧�� ����� � ����� ��-楫�}
procedure PRAVC2(T: real; var Y,G: MASN);
var Z,V: MASN;
    I: integer;
begin
  VARI:=true;
  BK:=BK0;
  for I:=1 to 42 do Z[I]:=Y[I];
  PRAV0(T,Z,V);
  for I:=1 to 42 do begin G[I]:=V[I]; Z[I]:=Y[I+42] end;
  BK:=BK1;
  PRAV0(T,Z,V);
  for I:=1 to 42 do G[I+42]:=V[I]
end;

{OUTP, OUT1 - �ணࠬ�� ��ࠡ�⪨ १���⮢ ��⥣�஢����,
��뢠���� �� ��� 㧫�� �⪨ ��⥣�஢����}
procedure OUTP(T,H: real; var Z: MASN; var P: MS5; var F: MASN9);
var I,J: integer;
    S,M: MS3;
begin
  {����� ������ ���樨, ��������� ��� � �����}
  NUTAC1(T);
  SUNMOON(T,RSOL,RMOO,S,M);
  for I:=1 to 3 do
  begin
    SOL[I]:=0; MOO[I]:=0;
    for J:=1 to 3 do
    begin
      SOL[I]:=SOL[I]+MATR[I,J]*S[J];
      MOO[I]:=MOO[I]+MATR[I,J]*M[J]
    end
  end
end;

{�ணࠬ�� ���௮��樨}
procedure INTER(Z: real; var D0,D1,D2,D3,D4,D5: real);
begin
  D0:=4*sqr((Z-1)*(Z-0.5))*(6*Z+1);
  D1:=4*Z*sqr((Z-1)*(Z-0.5));
  D2:=4*sqr(Z*(Z-0.5))*(7-6*Z);
  D3:=4*(Z-1)*sqr(Z*(Z-0.5));
  D4:=16*sqr(Z*(Z-1));
  D5:=16*sqr(Z*(Z-1))*(Z-0.5)
end;

function INTERP(Z: real; var Y: MASN; var F: MASN9): real;
var D,D0,D1,D2,D3,D4,D5,S: real;
    I: integer;
begin
  INTER(Z,D0,D1,D2,D3,D4,D5);
  for I:=1 to 6 do
    Y[I]:=D0*F[9,I]+D2*F[2,I]+D4*F[4,I]+D1*F[1,I]+D3*F[3,I]+D5*F[5,I];
  S:=0;
  for I:=1 to 3 do S:=S+Y[I+3]*Y[I];
  INTERP:=S
end;

{$I utils}
procedure OUT1(T,H: real; var Y: MASN; var P: MS5; var F: MASN9);
const B0=121/1536; B1=317125/861696; B2=397/6912;
      B3=17/704; B4=125/12096; B5=-817/14784;
      C30=0.866; C10=0.985;
var I,J,VID: integer;
    Y05,G05,G1: MASN;
    D,ZL,ZR,ZM,FM: real;
    S,M: MS3;
    X1,X2,X3,TZ: real;
begin
  {��ॢ�� �� ��᮫�⭮� ��⥬� ���न��� � �ਭ������}
  {TZ:=T+(ZM-1)*H;}
  TZ:=T;
  AbsToGsk(Y[1],Y[2],Y[3], W,T,SZ0, X1,X2,X3);
  writeln(FL_KA1, TZ,' ',X1,' ',X2,' ',X3);
  AbsToGsk(Y[1+6],Y[2+6],Y[3+6], W,T,SZ0, X1,X2,X3);
  writeln(FL_KA2, TZ,' ',X1,' ',X2,' ',X3);
  {����� ������ ���樨, ��������� ��� � �����}
  NUTAC1(T);
  SUNMOON(T,RSOL,RMOO,S,M);
  for I:=1 to 3 do
  begin
    SOL[I]:=0; MOO[I]:=0;
    for J:=1 to 3 do
    begin
      SOL[I]:=SOL[I]+MATR[I,J]*S[J];
      MOO[I]:=MOO[I]+MATR[I,J]*M[J]
    end
  end;
  {������ ����ﭨ� ����� ��, ��������� ����ﭨ� �
  �६� � ���ᨢ�� ��� �ᮢ���� ��䨪�}
  D:=0; for I:=1 to 3 do D:=D+sqr(Y[I+6]-Y[I]); D:=sqrt(D);
  ZL:=sqrt(sqr(Y[1])+sqr(Y[2])+sqr(Y[3]));
  ZR:=sqrt(sqr(Y[7])+sqr(Y[8])+sqr(Y[9]));
  writeln(FL1,T+(ZM-1)*H:15,' ',D:15,' ',ZL:13,' ',ZR:13); N:=N+1;
  {��� ��६��� ����� � - �� + �ந������� �� �६���
  �����ﭨ� ����� ��}
  if H=0 then
  begin
    U1:=0;
    for I:=1 to 3 do U1:=U1+(Y[I+6]-Y[I])*(Y[I+9]-Y[I+3]);
    exit
  end
  else
  begin
    U2:=0;
    for I:=1 to 3 do U2:=U2+(Y[I+6]-Y[I])*(Y[I+9]-Y[I+3]);
    {�᫨ �㦭�� ��६��� ����� ���, ��宦� �� ��楤���}
    if (U1>0) or (U2<0) then begin U1:=U2; exit end;
    {�㦭�� ��६��� ����� ����, ���� ���௮��樮��� �������,
    ���பᨬ����騩 ࠧ����� 䠢����� ����஢ �� ����� �������
    蠣� ��⥣�஢����}
    PRAVC1(T,Y,G1);
    for I:=1 to 12 do
      Y05[I]:=F[9,I]+H*(B0*F[1,I]+B1*F[3,I]+B2*F[4,I]-F[5,I]/68+
                        B3*F[6,I]+B4*F[7,I]+B5*F[8,I]+G1[I]/32);
    PRAVC1(T-H*0.5,Y05,G05);
    for I:=1 to 6 do
    begin
      F[1,I]:=H*(F[1,I+6]-F[1,I]); F[2,I]:=Y[I+6]-Y[I];
      F[3,I]:=H*(G1[I+6]-G1[I]); F[4,I]:=Y05[I+6]-Y05[I];
      F[5,I]:=H*(G05[I+6]-G05[I]); F[9,I]:=F[9,I+6]-F[9,I]
    end;
    {��⮤�� ������� ��१�� ������� (���⮬��) ��� ���
    �ந������� �� �६��� �����ﭨ� ����� ��}
    ZL:=0; ZR:=1;
    repeat
      ZM:=(ZL+ZR)/2; FM:=INTERP(ZM,Y05,F);
      if FM*(U1-U2)>0 then ZL:=ZM else ZR:=ZM
    until ZR-ZL<1e-8;
    D:=0; for I:=1 to 3 do D:=D+sqr(Y05[I]); D:=sqrt(D);
    {��襫 ������� ������ ����ﭨ�, �஢���� ��� ��
    ������쭮��� � �� �᫮��� 㤠������� ���ࠢ����� ��0-��1
    �� ���ࠢ����� �� ��� � �����, �����뢠� � 䠩�� ���
    �ᮢ���� ��䨪��}
    writeln(FL1,T+(ZM-1)*H:15,' ',D:15); N:=N+1;
    {�஢�ઠ �᫮��� �� ������ � �㭥}
    ZL:=(SOL[1]*Y05[1]+SOL[2]*Y05[2]+SOL[3]*Y05[3])/D;
    ZR:=(MOO[1]*Y05[1]+MOO[2]*Y05[2]+MOO[3]*Y05[3])/D/RMOO;
    if (ZL<C30) and (ZR<C10) then VID:=1 else VID:=0;
    writeln(FL2,T+(ZM-1)*H:15,' ',D:15,' ',VID); NN:=NN+1;
    {�஢�ઠ ����ﭨ� �� ������}
    if (D<SMIN) and (VID=1) then
    begin SMIN:=D; HMIN:=T+(ZM-1)*H end;
    {����஫쭠� �����}
    writeln(NN,' ',U1:12,' ',U2:12,' ',FM:13,' ',D:13);
    U1:=U2
  end
end;

{���᫥��� ��᪮� ����, I-���, J-�����, K-�᫮}
function JDATE(I,J,K: longint): longint;
begin
  JDATE:=K-32075+1461*(I+4800+(J-14) div 12 ) div 4 +
  367*(J-2-((J-14) div 12)*12) div 12 -
  3*((I+4900+(J-14) div 12) div 100) div 4
end;

{���寧� 䠮���� ����� � ����ਠ樮���� ������}
{$I RKV56}
procedure CONT(JD: longint; T,DL: real; var X: MS6; var CM: MS66);
var A,B: MS66;
    I,J,K: integer;
    D: real;
begin
  JD0:=JD; T0:=T;
  {���᫨� �������� �६� � 㣫���� ᪮���� �����}
  D:=(JD0-2451545.5)/36525.0;
  W:=pi/43.2*(1.002737909350795+D*(5.9006e-11-D*5.9e-15));
  SZ0:=24110.54841+D*(8640184.812866+D*(0.093104-D*6.2e-6));
  SZ0:=SZ0/86400.0; SZ0:=2*pi*(SZ0-trunc(SZ0))+W*T0;
  {����� ��砫�� �᫮���}
  P[1]:=0; P[2]:=DL; P[3]:=0.3; P[4]:=1e-9; P[5]:=0;
  VARI:=true; for I:=7 to 42 do Y[I]:=0;
  for I:=1 to 6 do begin Y[I]:=X[I]; Y[7*I]:=1 end;
  RKV56(42,Y,FQ,P,PRAV0,OUTP);
  for I:=1 to 6 do X[I]:=Y[I];
  {����� ����� ����ਠ樮���� ������}
  for I:=1 to 6 do for J:=1 to 6 do A[I,J]:=Y[I+6*J];
  for I:=1 to 6 do for J:=1 to 6 do
  begin
    B[I,J]:=0; for K:=1 to 6 do B[I,J]:=B[I,J]+A[I,K]*CM[K,J]
  end;
  for I:=1 to 6 do for J:=1 to 6 do
  begin
    CM[I,J]:=0; for K:=1 to 6 do CM[I,J]:=CM[I,J]+B[I,K]*A[J,K]
  end
end;

{���뢠��� ��室��� ������}
procedure PREP0;
var D,SEC,T1,T2: real;
    FL: text;
    JD1,JD2: longint;
    YY,MM,DD,HR,MIN,I: integer;
    X0,X1: MS6;
begin
  assign(FL,inputKA1);
  reset(FL);
  {��⠫ ��砫쭮� �६� ��0}
  readln(FL,YY,MM,DD,HR,MIN,SEC);
  {�����᪠� ���}
  JD1:=JDATE(YY,MM,DD);
  {�६� �� �ਭ���᪮� ���㭮�, 1000 �}
  T1:=0.001*SEC+0.06*MIN+3.6*HR;
  {��⠫ ��������᪨� �����樥�� ��0}
  readln(FL,BK0);
  {��⠫ ��砫�� �᫮��� ��0}
  for I:=1 to 6 do readln(FL,X0[I]);
  {��⠫ ����ਠ樮���� ������ ��0 (�� ��ப��)}
  for I:=1 to 6 do for J:=1 to 6 do readln(FL,CM0[I,J]);
  close(FL);
  writeln('ka0: JD, T, bk ',JD1,' ',T1:15,' ',BK0:11);
  for I:=1 to 6 do writeln(I,' ',X0[I]:15);

  assign(FL,inputKA2);
  reset(FL);
  {��⠫ ��砫쭮� �६� ��1}
  readln(FL,YY,MM,DD,HR,MIN,SEC);
  {�����᪠� ���}
  JD2:=JDATE(YY,MM,DD);
  {�६� �� �ਭ���᪮� ���㭮�, 1000 �}
  T2:=0.001*SEC+0.06*MIN+3.6*HR;
  {��⠫ ��������᪨� �����樥�� ��1}
  readln(FL,BK1);
  {��⠫ ��砫�� �᫮��� ��1}
  for I:=1 to 6 do readln(FL,X1[I]);
  {��⠫ ����ਠ樮���� ������ ��1 (�� ��ப��))}
  for I:=1 to 6 do for J:=1 to 6 do readln(FL,CM1[I,J]);
  close(FL);
  writeln('ka1: JD, T, bk ',JD2,' ',T2:15,' ',BK1:10);
  for I:=1 to 6 do writeln(I,' ',X1[I]:15);

  VARI:=true;
  if JD2+T2/86.4>JD1+T1/86.4 then
  begin
    CONT(JD1,T1,86.4*(JD2-JD1)+T2-T1,X0,CM0);
    JD0:=JD2; T0:=T2
  end
  else
  begin
    CONT(JD2,T2,86.4*(JD1-JD2)+T1-T2,X1,CM1);
    JD0:=JD1; T0:=T1
  end;
  {���᫨� �������� �६� � 㣫���� ᪮���� �����}
  D:=(JD0-2451545-0.5)/36525.0;
  W:=pi/43.2*(1.002737909350795+D*(5.9006e-11-D*5.9e-15));
  SZ0:=24110.54841+D*(8640184.812866+D*(0.093104-D*6.2e-6));
  SZ0:=SZ0/86400.0; SZ0:=2*pi*(SZ0-trunc(SZ0))+W*T0;
  {����� ��砫�� �᫮��� ᮢ���⭮�� ��⥣�஢����}
  for I:=1 to 6 do begin Y0[I]:=X0[I]; Y0[I+6]:=X1[I] end;
  {��⠫ ��ࠬ���� ��⥣�஢����}
  assign(FL,'parint.pas');
  reset(FL);
  readln(FL,P[2],P[3],P[4],S0);
  writeln('P[2..4],S0 ',P[2]:11,' ',P[3]:11,' ',P[4]:11,' ',S0:11);
  {readln;}
  P[2]:=P[2]*86.4; P[5]:=0;
  close(FL)
end;

begin
  if(ParamCount <> 3) then
  begin
    writeln('Syntax: enct.exe <ka> <km> <output ka> <output km>');
    halt(1)
  end
  else begin
    inputKA1 := ParamStr(1);
    inputKA2 := ParamStr(2);
    outputKA1 := ParamStr(3);
    outputKA2 := ParamStr(4);
  end;

  PREP0;
  {� 䠩� FL1 �����뢠���� ����� ��� �ᮢ���� ����ﭨ�}
  assign(FL1,'dist.pas');
  rewrite(FL1);
  {� 䠩� FL2 �����뢠���� ����� ��� �ᮢ���� ������饩
  �����㬮� ����ﭨ�}
  assign(FL2,'dist_m.pas');
  rewrite(FL2);
  {� 䠩� FL_KA1 �����뢠���� ����� �ࠥ��ਨ ��1}
  assign(FL_KA1,outputKA1);
  rewrite(FL_KA1);
  {� 䠩� FL_KA2 �����뢠���� ����� �ࠥ��ਨ ��2}
  assign(FL_KA2,outputKA2);
  rewrite(FL_KA2);
  {�ண�� �������� ��� ���᪠ �������쭮�� ����ﭨ�}
  P[1]:=0; Y:=Y0; N:=0; NN:=0; SMIN:=1e6;
  RKV56(12,Y,FQ,P,PRAVC1,OUT1);
  writeln(N,' ',NN);
  writeln('glob.min(km), time(d)',1e3*SMIN:14,' ',HMIN/86.4:14);
  close(FL1);
  close(FL2);
  close(FL_KA1); close(FL_KA2);
  reset(FL2);
  {��� � �뢮�� �� �࠭ �� ���稬� ������� (<150 ��) ����ﭨ�}
  while not eof(FL2) do
  begin
    readln(FL2,H,S,I);
    if S<0.15 then writeln(1e3*S:14,' ',H/86.4:14,' ',I)
  end;
  close(FL2);
  {writeln('Press Enter to continue!');
  readln;}
  assign(FL1,'results.pas');
  rewrite(FL1);
  writeln(FL1,'glob.min(km) time(d) ',1e3*SMIN:14,' ',HMIN/86.4:14);

  {����� ᡫ������, �����뢠� �襭�� � ���室�� ������}
  P[1]:=0; P[2]:=HMIN;
  for I:=1 to 6 do begin Y[I]:=Y0[I]; Y[I+42]:=Y0[I+6] end;
  for I:=7 to 42 do begin Y[I]:=0; Y[I+42]:=0 end;
  for I:=1 to 6 do begin Y[7*I]:=1; Y[7*I+42]:=1 end;
  RKV56(84,Y,FQ,P,PRAVC2,OUTP);
  {����� ����ਠ樮���� �����}
  for I:=1 to 6 do for J:=1 to 6 do A[I,J]:=Y[I+6*J];
  for I:=1 to 6 do for J:=1 to 6 do
  begin
    B[I,J]:=0; for K:=1 to 6 do B[I,J]:=B[I,J]+A[I,K]*CM0[K,J]
  end;
  for I:=1 to 6 do for J:=1 to 6 do
  begin
    CM0[I,J]:=0; for K:=1 to 6 do CM0[I,J]:=CM0[I,J]+B[I,K]*A[J,K]
  end;
  for I:=1 to 6 do for J:=1 to 6 do A[I,J]:=Y[I+6*J+42];
  for I:=1 to 6 do for J:=1 to 6 do
  begin
    B[I,J]:=0; for K:=1 to 6 do B[I,J]:=B[I,J]+A[I,K]*CM1[K,J]
  end;
  for I:=1 to 6 do for J:=1 to 6 do
  begin
    CM1[I,J]:=0; for K:=1 to 6 do CM1[I,J]:=CM1[I,J]+B[I,K]*A[J,K]
  end;
  {���᫥��� ����ਠ樮���� ������ �_0+�_1, ����� �� ENCOUNTER}
  for J:=1 to 6 do
    for K:=1 to 6 do A[J,K]:=CM0[J,K]+CM1[J,K];
  {���᫥��� ����஢ R � V}
  for I:=1 to 6 do Z[I]:=Y[I+42]-Y[I];
  DL:=sqrt(sqr(Z[1])+sqr(Z[2])+sqr(Z[3]));
  writeln('chek distance(km) ',1e3*DL:14);
  writeln(FL1,'chek distance(km) ',1e3*DL:14);
  {���᫥��� ��� s}
  for I:=1 to 3 do LV[I]:=Z[I]/DL;
  S:=LV[1]*Z[4]+LV[2]*Z[5]+LV[3]*Z[6];
  D:=sqrt(sqr(Z[4])+sqr(Z[5])+sqr(Z[6])-sqr(S));
  {���᫥��� ��� e_1=EV}
  for I:=1 to 3 do EV[I]:=(Z[I+3]-S*LV[I])/D;
  {���᫥��� ��� e_2=EN}
  EN[1]:=LV[2]*EV[3]-LV[3]*EV[2];
  EN[2]:=LV[3]*EV[1]-LV[1]*EV[3];
  EN[3]:=LV[1]*EV[2]-LV[2]*EV[1];
  D:=0;
  for J:=1 to 3 do
    for K:=1 to 3 do D:=D+A[J,K]*EN[J]*EN[K];
  {���᫥��� sigma ��� �訡�� � 䠧���� ������ � 㣫. ���.}
  D:=sqrt(D)/DL/SC;
  writeln('sigma p.v.(arcmin) ',D:12);
  writeln(FL1,'sigma p.v.(arcmin) ',D:12);

  Y:=Y0; P[1]:=0; S:=HMIN+S0; P[2]:=HMIN-S0;
  while P[2]<=S do
  begin
    RKV56(12,Y,FQ,P,PRAVC1,OUTP);
    for I:=1 to 3 do Z[I]:=Y[I+6]-Y[I];
    D:=sqrt(sqr(Z[1])+sqr(Z[2])+sqr(Z[3]));
    for I:=1 to 3 do Z[I]:=Z[I]/D;
    XN:=0; XV:=0;
    for I:=1 to 3 do
    begin
      XV:=XV+EV[I]*Z[I]; XN:=XN+EN[I]*Z[I]
    end;
    writeln(P[2]:20,' ',XV/SC:13,' ',XN/SC:13);
    writeln(FL1,P[2]:20,' ',XV/SC:13,' ',XN/SC:13);
    P[1]:=P[2]; P[2]:=P[2]+S0/10
  end;
  close(FL1)
end.
 {
�ணࠬ�� ��⥣�஢���� ��⥬� ���
procedure RKV56(N: integer; var Y: MASN; var F: MASN9; var P: MS5;
                PRAVCH: INTPR; OUTP: INTOU);

PRAVCH - ��楤�� ���� �ࠢ�� ��⥩ ��� (�. ���)
OUTP - ��楤�� ��ࠡ�⪨ १���⮢ ������� 蠣� ��⥣�஢����
       (�. ���)
N - ࠧ��୮��� ��⥬�
Y - �� �室� - ��砫쭮� ���祭�� 䠧����� ����� ��⥬�,
    �� ��室� - ���祭�� 䠧����� ����� � ����筮� �窥 ��⥣�஢����
P - �㦥��� ���ᨢ
    P[1] - ��砫쭠� �窠 ��⥣�஢����
    P[2] - ����筠� �窠 ��⥣�஢����
    P[3] - ��砫�� 蠣 ��⥣�஢����, ��� 䠪��᪮� ���祭��
           �롨ࠥ��� ��⮬���᪨ �� �����쭮� �筮��
    P[4] - �����쭠� �筮���
    P[5] - �ࠢ���騩 ��ࠬ��� ���砫� P[5]=0, � OUTP ��� �����
           ��������, ⮣�� ��⥣�஢���� ��ࢥ��� (���ਬ��, ���ᨢ�
           ��� ��䨪�� ���������)
}
