program ENCOUNTER3;
uses graph, GRAPHICS;
const NMAX=17000; S0=1e-5; SC=3.14159265359/1.8e2/60;
type real=double;
     MS3=array[1..3] of real;
     MS5=array[1..5] of real;
     MS6=array[1..6] of real;
     MASN=array[1..84] of real;
     MS33=array[1..3,1..3] of real;
     MS66=array[1..6,1..6] of real;
     MASN9=array[1..9] of MASN;
     MSGR=array[1..NMAX] of single;
     INTPR=procedure(T: real; var Y,G: MASN);
     INTOU= procedure(T,H: real; var Z: MASN; var P: MS5; var F: MASN9);

var VARI: boolean;
var CM0,CM1,A,B: MS66;
    Y,Y0: MASN;
    FQ: MASN9;
    P: MS5;
    VID: array[1..NMAX] of boolean;
    Z0,Z1,ZZ0,ZZ1: MSGR;
    T0,RSOL,RMOO,SZ0,W,BK,BK0,BK1,S,DL,D,XN,XV,
    SMIN,HMIN,U1,U2: real;
    JD0: longint;
    I,J,K,N,NN: integer;
    SOL,MOO: MS3;
    MATR: MS33;
    Z: MS6;
    LV,EV,EN: MS3;


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
{расчет правой части дифф. ур-ния движения КА
Y[1..3] - рад. вектор, Y[4..6] - скорость, dY/dT=G(T,Y)}
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

{расчет правой части дифф. ур-ния движения двух КА;
Y[1..6] - фазовый вектор КА-наблюдателя,
Y[7..12] - фазовый вектор КА-цели}
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

{расчет правой части дифф. ур-ния движения двух КА и соответствующих
переходных матриц; Y[1..42] - фазвый вектор и переходная 6*6-матрица по
столбцам КА-наблюдателя, Y[43..84] - фазвый вектор и матрица КА-цели}
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

{OUTP, OUT1 - программы обработки результатов интегрирования,
вызываются во всех узлах сетки интегрировагия}
procedure OUTP(T,H: real; var Z: MASN; var P: MS5; var F: MASN9);
var I,J: integer;
    S,M: MS3;
begin
  {Расчет матрицы нутации, положений Луны и Солнца}
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

{Программы интерполяции}
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

procedure OUT1(T,H: real; var Y: MASN; var P: MS5; var F: MASN9);
const B0=121/1536; B1=317125/861696; B2=397/6912;
      B3=17/704; B4=125/12096; B5=-817/14784;
      C30=0.866; C10=0.985;
var I,J: integer;
    Y05,G05,G1: MASN;
    D,ZL,ZR,ZM,FM: real;
    S,M: MS3;
begin
  {Расчет матрицы нутации, положений Луны и Солнца}
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
  {вычисляю расстояние между КА, запоминаю расстояние и
  время в массивах для рисования графика}
  D:=0; for I:=1 to 3 do D:=D+sqr(Y[I+6]-Y[I]);
  D:=sqrt(D);
  N:=N+1; Z0[N]:=T/86.4; Z1[N]:=D;
  if N>=NMAX then begin P[5]:=1; exit end;
  {Ищу перемену знака с - на + производной по времени
  расстряния между КА}
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
    {если нужной перемены знака нет, выхожу из процедуры}
    if (U1>0) or (U2<0) then begin U1:=U2; exit end;
    {нужная перемена знака есть, строю интерполяционный полином,
    аппроксимирующий разность фавзовых векторов КА внутри данного
    шага интегрирования}
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
    {Методом деления отрезка пополам (дихотомия) ищу нуль
    производной по времени расстряния между КА}
    ZL:=0; ZR:=1;
    repeat
      ZM:=(ZL+ZR)/2; FM:=INTERP(ZM,Y05,F);
      if FM*(U1-U2)>0 then ZL:=ZM else ZR:=ZM
    until ZR-ZL<1e-8;
    D:=0;
    for I:=1 to 3 do D:=D+sqr(Y05[I]);
    D:=sqrt(D);
    {Нашел локальный минимум расстояния, проверяю его на
    глобальность и на условия удаленности направления КА0-КА1
    от направлений на Луну и Солнце, записываю в массивы для
    рисования графиков}
    NN:=NN+1; ZZ0[NN]:=(T+(ZM-1)*H)/86.4; ZZ1[NN]:=D;
    N:=N+1; Z0[N]:=ZZ0[NN]; Z1[N]:=D;
    ZL:=(SOL[1]*Y05[1]+SOL[2]*Y05[2]+SOL[3]*Y05[3])/D;
    ZR:=(MOO[1]*Y05[1]+MOO[2]*Y05[2]+MOO[3]*Y05[3])/D/RMOO;
    {Проверка условий по Солнцу и Луне}
    VID[NN]:=(ZL<C30) and (ZR<C10);
    {Проверка расстояния на минимум}
    if D<SMIN then begin SMIN:=D; HMIN:=T+(ZM-1)*H end;
    {контрольная печать}
    writeln(NN,' ',U1:12,' ',U2:12,' ',FM:13,' ',D:13);
    U1:=U2
  end
end;

{Вычисление юлианской даты, I-год, J-месяц, K-число}
function JDATE(I,J,K: longint): longint;
begin
  JDATE:=K-32075+1461*(I+4800+(J-14) div 12 ) div 4 +
  367*(J-2-((J-14) div 12)*12) div 12 -
  3*((I+4900+(J-14) div 12) div 100) div 4
end;

{Считывание исходных данных}
procedure PREP0;
var D,SEC: real;
    FL: text;
    YY,MM,DD,HR,MIN,I: integer;
begin
  assign(FL,'list.pas');
  reset(FL);
  {Считал начальное время, одинаковое для обоих КА}
  readln(FL,YY,MM,DD,HR,MIN,SEC);
  {Юлианская дата}
  JD0:=JDATE(YY,MM,DD);
  {Время от гринвичской полуночи, 1000 с}
  T0:=0.001*SEC+0.06*MIN+3.6*HR;
  {Вычислил угловую скорость вращения Земли и звездное время}
  D:=(JD0-2451545-0.5)/36525.0;
  W:=pi/43.2*(1.002737909350795+D*(5.9006e-11-D*5.9e-15));
  SZ0:=24110.54841+D*(8640184.812866+D*(0.093104-D*6.2e-6));
  SZ0:=SZ0/86400.0; SZ0:=2*pi*(SZ0-trunc(SZ0))+W*T0;
  readln(FL,BK0,BK1);
  {Считал начальные условия для обоих КА}
  for I:=1 to 12 do readln(FL,Y0[I]);
  {Надо еще считать ковариационные матрицы, но пока они
  задаются явно и примитивно}
  for I:=1 to 6 do for J:=1 to 6 do
  begin CM0[I,J]:=0; CM1[I,J]:=0 end;
  for I:=1 to 6 do
  begin CM0[I,I]:= CM1[I,I]:=1e-10 end;

  close(FL);
  {Задал параметры интегрирования}
  P[3]:=0.3; P[4]:=1e-9; P[5]:=0
end;

{$I RKV56}
begin
  PREP0;
  {Прогон движения для поиска минимального расстояния}
  P[1]:=0; P[2]:=1e5; Y:=Y0; N:=0; NN:=0; SMIN:=1e6;
  RKV56(12,Y,FQ,P,PRAVC1,OUT1);
  writeln(N,' ',NN);
  {Рисую график взаимного расстояния на отрезке времени, длина которого
  определяется длиной соответствующих массивов, и график "огибающей"
  локальных минимумов этого расстояния}
  GrInit;
  GRAFIK1(N,Z0,Z1,80,15,600,130,10,0,2,4,1,10,0,2,9,3);
  GRAFIK1(NN,ZZ0,ZZ1,80,175,600,290,10,0,2,5,4,10,0,2,9,3);
  {Ищу и вывожу на экран все значимые минимумы (<150 км) расстояния}
  for J:=1 to NN do if ZZ1[J]<0.15 then writeln(ZZ1[J]:14,' ',ZZ0[J]:14);
  writeln('glob. min.',SMIN:14,' ',HMIN/86.4:14);
  SMIN:=0.15;
  {Ищу и вывожу н экран все минимумы, допустимые по Солнце и Луне}
  for J:=1 to NN do if (ZZ1[J]<SMIN) and VID[NN] then
  begin
    HMIN:=ZZ0[J]; SMIN:=ZZ1[J]; writeln(SMIN:14,' ',HMIN:14)
  end;
  writeln('glob. min. with additional conditions ',SMIN:14,' ',HMIN:14);
  HMIN:=HMIN*86.4;
  writeln('Press Enter to continue!');
  readln;

  {Расчет сближения, рассчитываю решение и переходные матрицы}
  P[1]:=0; P[2]:=HMIN; P[5]:=0;
  for I:=1 to 6 do begin Y[I]:=Y0[I]; Y[I+42]:=Y0[I+6] end;
  for I:=7 to 42 do begin Y[I]:=0; Y[I+42]:=0 end;
  for I:=1 to 6 do begin Y[7*I]:=1; Y[7*I+42]:=1 end;
  RKV56(84,Y,FQ,P,PRAVC2,OUTP);
  {Расчет ковариационных матриц}
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
  {Вычисление ковариационной матрицы С_0+С_1, далее по ENCOUNTER}
  for J:=1 to 6 do
    for K:=1 to 6 do A[J,K]:=CM0[J,K]+CM1[J,K];
  {Вычисление векторов R и V}
  for I:=1 to 6 do Z[I]:=Y[I+42]-Y[I];
  DL:=sqrt(sqr(Z[1])+sqr(Z[2])+sqr(Z[3]));
  writeln('distance ',DL:14);
  {Вычисление орта s}
  for I:=1 to 3 do LV[I]:=Z[I]/DL;
  S:=LV[1]*Z[4]+LV[2]*Z[5]+LV[3]*Z[6];
  D:=sqrt(sqr(Z[4])+sqr(Z[5])+sqr(Z[6])-sqr(S));
  {Вычисление орта e_1=EV}
  for I:=1 to 3 do EV[I]:=(Z[I+3]-S*LV[I])/D;
  {Вычисление орта e_2=EN}
  EN[1]:=LV[2]*EV[3]-LV[3]*EV[2];
  EN[2]:=LV[3]*EV[1]-LV[1]*EV[3];
  EN[3]:=LV[1]*EV[2]-LV[2]*EV[1];
  D:=0;
  for J:=1 to 3 do
    for K:=1 to 3 do D:=D+A[J,K]*EN[J]*EN[K];
  {Вычисление sigma для ошибок в фазовых векторах в угл. мин.}
  D:=sqrt(D)/DL/SC;
  writeln('sigma p.v,',D:12);


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
    P[1]:=P[2]; P[2]:=P[2]+S0/10
  end
end.

Программа интегрирования системы ОДУ
procedure RKV56(N: integer; var Y: MASN; var F: MASN9; var P: MS5;
                PRAVCH: INTPR; OUTP: INTOU);

PRAVCH - процедура расчета правых частей ОДУ (см. выше)
OUTP - процедура обработки результатов каждого шага интегрирования
       (см. выше)
N - размерность системы
Y - на входе - начальное значение фазового вектора системы,
    на выходе - значение фазового вектора в конечной точке интегрирования
P - служебный массив
    P[1] - начальная точка интегрирования
    P[2] - конечная точка интегрирования
    P[3] - начальный шаг интегрирования, его фактическое значение
           выбирается автоматически по локальной точности
    P[4] - локальная точность
    P[5] - управляющий параметр вначале P[5]=0, в OUTP его можно
           изменить, тогда интегрирование прервется (например, массивы
           для графиков заполнены)
