program VIEW_E;
uses graph, GRAPHICS;
const GC=398.60044; NMAX=17000;
type real=double;
     MS3=array[1..3] of real;
     MS5=array[1..5] of real;
     MS6=array[1..6] of real;
     MASN=array[1..12] of real;
     MS33=array[1..3,1..3] of real;
     MSGR=array[1..NMAX] of single;
     GLDARRAY=MSGR;
     INTPR=procedure(T: real; var Y,G: MASN);
     INTOU= procedure(T,H: real; var Z: MASN; var P: MS5);

var X,X1: MS6;
    Y,Y0: MASN;
    P: MS5;
    T0,RSOL,RMOO,SZ0,W,BK,BK0,BK1,S,H,F,F1,STEPF: real;
    I,J,N,NN: integer;
    JD0: longint;
    Z0,Z1,Z2: MSGR;
    U: GLDARRAY;
    SOL,MOO: MS3;
    MATR: MS33;

function AMAX1(X,Y: real): real;
begin
  if X>Y then AMAX1:=X else AMAX1:=Y
end;

function AMIN1(X,Y: real): real;
begin
  if X>Y then AMIN1:=Y else AMIN1:=X
end;

procedure GRAFIK1(N: integer; var X,Y: MSGR; X1,Y1,X2,Y2,
                 IX,JX,KX,LX,MX,IY,JY,KY,LY,MY: integer);
var GR: array [1..NMAX] of VECT2;
    Pmin,Pmax: VECT2;
    I: integer;
begin
  Pmin.X:=1e20; Pmax.X:=-1e20;
  Pmin.Y:=1e20; Pmax.Y:=-1e20;
  for I:=1 to N do
  begin
    Pmin.X:=AMIN1(Pmin.X,X[I]); Pmax.X:=AMAX1(Pmax.X,X[I]);
    Pmin.Y:=AMIN1(Pmin.Y,Y[I]); Pmax.Y:=AMAX1(Pmax.Y,Y[I]);
    GR[I].X:=X[I]; GR[I].Y:=Y[I]
  end;
  PlotAxis (Pmin, Pmax, X1, Y1, X2, Y2,
            14, 10,
            IX, JX, KX, LX, MX, IY, JY, KY, LY, MY );
  PlotLine (GR, N, 1, 11)
end;


{$I HARMON36}
{$I ROATMST}
{$I NUTAC1}
{$I SUNMOON}
procedure PRAV0(T: real; var Y,G: MS6);
const R0=6.378136; J2EARTH=1082625.7E-9; NGAM=8;
      GCS=332946.0*GC; GCM=GC/81.30068;
var A,H,R2,RV,SZ,SA,CA,VX,VY,DRM,RS,RM,SR: real;
    X,F,LN,SL: MS3;
    I: integer;
begin
  G[1]:=Y[4]; G[2]:=Y[5]; G[3]:=Y[6];
  R2:=sqr(Y[1])+sqr(Y[2])+sqr(Y[3]);
  H:=-GC/R2/sqrt(R2);
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
  G[6]:=G[6]+DRM*(LN[3]-Y[3])-RM*LN[3]+RS*(SR*SL[3]-Y[3])
end;

procedure PRAVCH(T: real; var Y,G: MASN);
var Z,V: MS6;
    I: integer;
begin
  BK:=BK0;
  for I:=1 to 6 do Z[I]:=Y[I];
  PRAV0(T,Z,V);
  for I:=1 to 6 do begin G[I]:=V[I]; Z[I]:=Y[I+6] end;
  BK:=BK1;
  PRAV0(T,Z,V);
  for I:=1 to 6 do G[I+6]:=V[I]
end;

procedure OUTP(T,H: real; var Z: MASN; var P: MS5);
var I,J: integer;
    S,M: MS3;
begin
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

{$I RKDP78}
{$I ARCTAN2}
procedure KEPLER(T: real; var Y,X: MS6);
var CX,CY,CZ,C,FX,FY,FZ,F,PX,PY,PZ,QX,QY,QZ,RX,RY,RZ,R,R2,A,
    E,E1,N0,SE,CE,INC,UZL,ARG,M0,Q: real;
    I: integer;
begin
  CX:=Y[2]*Y[6]-Y[3]*Y[5];
  CY:=Y[3]*Y[4]-Y[1]*Y[6];
  CZ:=Y[1]*Y[5]-Y[2]*Y[4];
  C:=sqrt(sqr(CX)+sqr(CY)+sqr(CZ));
  RX:=CX/C; RY:=CY/C; RZ:=CZ/C;
  R2:=sqr(Y[1])+sqr(Y[2])+sqr(Y[3]);
  R:=sqrt(R2);
  FX:=Y[5]*CZ-Y[6]*CY-GC*Y[1]/R;
  FY:=Y[6]*CX-Y[4]*CZ-GC*Y[2]/R;
  FZ:=Y[4]*CY-Y[5]*CX-GC*Y[3]/R;
  F:=sqrt(sqr(FX)+sqr(FY)+sqr(FZ));
  PX:=FX/F; PY:=FY/F; PZ:=FZ/F;
  QX:=RY*PZ-RZ*PY; QY:=RZ*PX-RX*PZ; QZ:=RX*PY-RY*PX;
  A:=sqr(Y[4])+sqr(Y[5])+sqr(Y[6])-2*GC/R; A:=-GC/A;
  E:=F/GC; E1:=sqrt(1-sqr(E));
  N0:=sqrt(GC/A)/A;
  CE:=(PX*Y[1]+PY*Y[2]+PZ*Y[3])/A+E;
  SE:=(QX*Y[1]+QY*Y[2]+QZ*Y[3])/A/E1;
  M0:=ATAN2(SE,CE)-E*SE-N0*T;
  while M0<0 do M0:=M0+2*pi;
  while M0>2*pi do M0:=M0-2*pi;
  Q:=sqrt(sqr(RX)+sqr(RY));
  INC:=ATAN2(Q,RZ); UZL:=ATAN2(RX,-RY); ARG:=ATAN2(PZ,QZ);
  X[1]:=A; X[2]:=E; X[3]:=M0; X[4]:=UZL; X[5]:=INC; X[6]:=ARG
end;

function ANG(var Y: MASN): real;
var C,D,CX,CY,CZ,DX,DY,DZ: real;
begin
  CX:=Y[2]*Y[6]-Y[3]*Y[5];
  CY:=Y[3]*Y[4]-Y[1]*Y[6];
  CZ:=Y[1]*Y[5]-Y[2]*Y[4];
  C:=sqrt(sqr(CX)+sqr(CY)+sqr(CZ));
  DX:=Y[8]*Y[12]-Y[9]*Y[11];
  DY:=Y[9]*Y[10]-Y[7]*Y[12];
  DZ:=Y[7]*Y[11]-Y[8]*Y[10];
  D:=sqrt(sqr(DX)+sqr(DY)+sqr(DZ));
  C:=(CX*DX+CY*DY+CZ*DZ)/C/D;
  D:=sqrt(1-sqr(C));
  ANG:=ATAN2(D,C)*180/pi
end;

procedure PLANES;
{рисует 4 графика; по осям абсцисс - время в сут.;
по осям ординат:
1 - угол между нормалями к плоскостям орбит, направление нормалей
    согласовано с орб. движением;
2 - долгота восх. узла орбиты 0;
3 - долгота восх. узла орбиты 1;
4 - разность этих долгот (1-0) }
var X,Z: MS6;
    N,I: integer;
begin
  Y:=Y0; Z0[1]:=0; U[1]:=ANG(Y);
  for J:=1 to 6 do X[J]:=Y[J];
  KEPLER(0,X,Z); Z1[1]:=Z[4]*180/pi;
  for J:=1 to 6 do X[J]:=Y[J+6];
  KEPLER(0,X,Z); Z2[1]:=Z[4]*180/pi;
  H:=43.2; N:=200;
  for I:=1 to N do
  begin
    P[1]:=H*(I-1); P[2]:=H*I;
    RKDP78(12,Y,P,PRAVCH,OUTP);
    Z0[I+1]:=P[2]/86.4; U[I+1]:=ANG(Y);
    for J:=1 to 6 do X[J]:=Y[J];
    KEPLER(P[2],X,Z); Z1[I+1]:=Z[4]*180/pi;
    for J:=1 to 6 do X[J]:=Y[J+6];
    KEPLER(P[2],X,Z); Z2[I+1]:=Z[4]*180/pi;
    if I mod 20 =1 then writeln(I,' ',P[2]:10:3)
  end;
  N:=N+1;
  for I:=2 to N do
  begin
    while U[I]<U[I-1]-100 do U[I]:=U[I]+360;
    while U[I]>U[I-1]+100 do U[I]:=U[I]-360;
    while Z1[I]<Z1[I-1]-100 do Z1[I]:=Z1[I]+360;
    while Z1[I]>Z1[I-1]+100 do Z1[I]:=Z1[I]-360;
    while Z2[I]<Z2[I-1]-100 do Z2[I]:=Z2[I]+360;
    while Z2[I]>Z2[I-1]+100 do Z2[I]:=Z2[I]-360;
  end;
  GrInit;
  GRAFIK1(N,Z0,U,80,15,600,130,10,0,2,3,1,10,0,2,6,1);
  GRAFIK1(N,Z0,Z1,80,175,600,290,10,0,2,3,1,10,0,2,6,1);
  GRAFIK1(N,Z0,Z2,80,335,600,450,10,0,2,3,1,10,0,2,6,1);
  for I:=1 to N do Z1[I]:=Z2[I]-Z1[I];
  GRAFIK1(N,Z0,Z1,80,495,600,610,10,0,2,3,1,10,0,2,6,1);
  writeln('Press Enter to continue!');
  readln;
  CloseGraph
end;

function JDATE(I,J,K: longint): longint;
begin
  JDATE:=K-32075+1461*(I+4800+(J-14) div 12 ) div 4 +
  367*(J-2-((J-14) div 12)*12) div 12 -
  3*((I+4900+(J-14) div 12) div 100) div 4
end;

{$I RV_TO_EL}
{$I EL_TO_RV}
procedure PREP;
var D,W: real;
    FL: text;
begin
  {Задал начальное время - 04:00:00.0 UTC 05.02.2011,
  одинаковое для обоих КА}
  JD0:=JDATE(2011,2,5);
  T0:=4*3.6;
  {Вычислил угловую скорость вращения Земли и звездное время}
  D:=(JD0-2451545-0.5)/36525.0;
  W:=pi/43.2*(1.002737909350795+D*(5.9006e-11-D*5.9e-15));
  SZ0:=24110.54841+D*(8640184.812866+D*(0.093104-D*6.2e-6));
  SZ0:=SZ0/86400.0; SZ0:=2*pi*(SZ0-trunc(SZ0))+W*T0;
  {Задал нач. условия орбиты 0}
  Y[1]:= 1.593010645;
  Y[2]:=-0.8282701672;
  Y[3]:= 6.960280794;
  Y[4]:= 2.911963361-W*Y[2];
  Y[5]:=-6.78221909 +W*Y[1];
  Y[6]:=-1.472114314;
  for I:=1 to 6 do X1[I]:=Y[I]; BK0:=0.0016;
  {Вычислил элементы орбиты 0}
  RV_EL(X1,X);
  {Вычислил кепл.период и циклич. орбит. частоту}
  H:=2*pi*X[1]*sqrt(X[1]/GC); F:=1/H;
  writeln('0: P,f ',H:14,' ',F:14);
  {Формирую кепл.элементы орбиты 1}
  X[1]:=X[1]+0.06; X[5]:=pi*51.6/180; X[4]:=X[4]+pi/2;
  {Вычислил кепл.период и циклич. орбит. частоту}
  H:=2*pi*X[1]*sqrt(X[1]/GC); F1:=1/H;
  {Задал нач. условия орбиты 1}
  EL_RV(X,X1);
  for I:=1 to 6 do Y[I+6]:=X1[I]; BK1:=0.002;
  writeln('1: P,f ',H:14,' ',F1:14);
  writeln('freq. ',abs(F1-F):14,' ',F+F1:14);
  P[3]:=0.3; P[4]:=1e-9; P[5]:=0; Y0:=Y;
  {
  assign(FL,'list_n.pas');
  rewrite(FL);
  for I:=1 to 12 do writeln(FL,Y[I]);
  close(FL)
  }
end;

procedure PREP0;
var D,W,SEC: real;
    FL: text;
    YY,MM,DD,HR,MIN,I: integer;
begin
  assign(FL,'list_n.pas');
  reset(FL);
  {Считал начальное время, одинаковое для обоих КА}
  readln(FL,YY,MM,DD,HR,MIN,SEC);
  {Юлианская дата}
  JD0:=JDATE(YY,MM,DD);
  {Время от гринвичской полуночи, 1000 с}
  T0:=0.001*SEC+0.06*MIN+3.6*HR;
  readln(FL,BK0,BK1);
  {Считал начальные условия для обоих КА}
  for I:=1 to 12 do readln(FL,Y0[I]);
  close(FL);
  for I:=1 to 6 do X1[I]:=Y0[I];
  {Вычислил элементы орбиты 0}
  RV_EL(X1,X);
  {Вычислил кепл.период и циклич. орбит. частоту}
  H:=2*pi*X[1]*sqrt(X[1]/GC); F:=1/H;
  writeln('0: Period, freq. ',H:14,' ',F:14);
  for I:=1 to 6 do X1[I]:=Y0[I+6];
  {Вычислил элементы орбиты 1}
  RV_EL(X1,X);
  {Вычислил кепл.период и циклич. орбит. частоту}
  H:=2*pi*X[1]*sqrt(X[1]/GC); F1:=1/H;
  writeln('1: Period, freq. ',H:14,' ',F1:14);
  {Вычислил комбинационные частоты}
  writeln('comb. freq. ',abs(F1-F):14,' ',F+F1:14);
  writeln('Press Enter to continue!');
  readln;
  {Задал параметры интегрирования}
  P[3]:=0.3; P[4]:=1e-9; P[5]:=0
end;

{$I FOUR2}
{$I REALFFT}
begin
  {PREP;}
  PREP0;
  {Строю график расстояния между КА (1е3 км) в функции времени (сут)}
  P[2]:=0; Y:=Y0; N:=0;
  repeat
    P[1]:=P[2]; P[2]:=P[1]+H;
    RKDP78(12,Y,P,PRAVCH,OUTP);
    S:=0; for I:=1 to 3 do S:=S+sqr(Y[I]-Y[I+6]);
    S:=sqrt(S); N:=N+1; Z0[N]:=P[2]/86.4; Z1[N]:=S;
    H:=0.2;
    if S<2.0 then H:=0.1;
    if S<1.0 then H:=0.01;
    if S<0.3 then H:=0.0001;
    if N mod 200=1 then writeln(N,' ',P[2]:10:2,' ',S:14)
  until N>=NMAX;
  GrInit;
  GRAFIK1(N,Z0,Z1,80,15,600,130,10,0,2,4,1,10,0,2,9,3);
  {Ищу значимые минимумы (<150 км) расстояния}
  for J:=2 to N-1 do
    if (Z1[J]<Z1[J-1]) and (Z1[J]<Z1[J+1]) and (Z1[J]<0.15) then
    begin
      writeln(Z0[J]:14,' ',Z1[J]:14)
    end;
  writeln('Press Enter to continue!');
  readln;
  {Графики спектров, график 1 - общая картина, график 2 - окр. нуля}
  Y:=Y0;
  NN:=4096*2; {должно быть NN<NMAX}
  N:=NN div 2;
  H:=0.47; STEPF:=1/H/NN; {шаги по времени и частоте}
  for I:=1 to NN do
  begin
    P[1]:=H*(I-1); P[2]:=H*I;
    RKDP78(12,Y,P,PRAVCH,OUTP);
    S:=0; for J:=1 to 3 do S:=S+sqr(Y[J]-Y[J+6]);
    U[I]:=sqrt(S);
    if I mod 1000=1 then writeln(I,' ',P[2]:10:2,' ',S:14)
  end;
  REALFFT(U,N,-1);
  for I:=1 to NN do U[I]:=U[I]/NN;
  Z0[1]:=0; Z1[1]:=0; Z0[N+1]:=STEPF*N; Z1[N+1]:=0;
  for J:=1 to N-1 do
  begin
    Z0[J+1]:=STEPF*J;
    Z1[J+1]:=2*sqrt(sqr(U[2*J+1])+sqr(U[2*J+2]));
  end;
  {Оценки значимых max: номер в массиве, частота и амплитуда}
  for J:=2 to N-2 do
    if (Z1[J]>Z1[J-1]) and (Z1[J]>Z1[J+1]) and (Z1[J]>1) then
    begin
      F:=J-1-0.5*(Z1[J+1]-Z1[J-1])/(Z1[J+1]-2*Z1[J]+Z1[J-1]);
      F1:=0.5*(F-J+1)*((Z1[J+1]-2*Z1[J]+Z1[J-1])*(F-J+1)+
          Z1[J+1]-Z1[J-1])+Z1[J];
      writeln(J,' ',F*STEPF:14,' ',F1:14);
    end;
  GRAFIK1(N+1,Z0,Z1,80,175,600,290,10,0,2,5,4,10,0,2,9,3);
  N:=50;
  GRAFIK1(N+1,Z0,Z1,80,335,600,450,10,0,2,5,4,10,0,2,9,3);
  writeln('Press Enter to continue!');
  readln;
  CloseGraph;
  PLANES
end.
