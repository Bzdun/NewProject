{Считывание исходных данных в требуемом формате,
но ков. матр. - для абс. с.к.}
procedure PREP0;
var D,SEC,X,Y,T1,T2,SA,CA: real;
    FL: text;
    JD1,JD2: longint;
    YY,MM,DD,HR,MIN,I: integer;
    X0,X1: MS6;
begin
  assign(FL,'ka0.pas');
  reset(FL);
  {Считал начальное время КА0}
  readln(FL,YY,MM,DD,HR,MIN,SEC);
  {Юлианская дата}
  JD1:=JDATE(YY,MM,DD);
  {Время от гринвичской полуночи, 1000 с}
  T1:=0.001*SEC+0.06*MIN+3.6*HR;
  {Считал баллистический коэффициент КА0}
  readln(FL,BK0);
  {Считал начальные условия КА0 (гск, м, м/с)}
  for I:=1 to 6 do readln(FL,X0[I]);
  for I:=1 to 3 do
  begin X0[I]:=1e-6*X0[I]; X0[I+3]:=1e-3*X0[I+3] end;
  {Вычислил звездное время}
  D:=(JD1-2451545.5)/36525.0;
  W:=pi/43.2*(1.002737909350795+D*(5.9006e-11-D*5.9e-15));
  SZ0:=24110.54841+D*(8640184.812866+D*(0.093104-D*6.2e-6));
  SZ0:=SZ0/86400.0; SZ0:=2*pi*(SZ0-trunc(SZ0))+W*T1;
  CA:=cos(SZ0); SA:=sin(SZ0);
  {Перешел к абс. сист.коор.}
  X:=X0[1]*CA-X0[2]*SA; Y:=X0[1]*SA+X0[2]*CA;
  X0[1]:=X; X0[2]:=Y;
  X:=X0[4]*CA-X0[5]*SA-W*X0[2]; Y:=X0[4]*SA+X0[5]*CA+W*X0[1];
  X0[4]:=X; X0[5]:=Y;
  {Считал ковариационную матрицу КА0 (по строкам)}
  for I:=1 to 6 do for J:=1 to 6 do readln(FL,CM0[I,J]);
  close(FL);
  writeln('ka0: JD, T, bc ',JD1,' ',T1:15,' ',BK0:11);
  for I:=1 to 6 do writeln(I,' ',X0[I]:15);

  assign(FL,'ka1.pas');
  reset(FL);
  {Считал начальное время КА1}
  readln(FL,YY,MM,DD,HR,MIN,SEC);
  {Юлианская дата}
  JD2:=JDATE(YY,MM,DD);
  {Время от гринвичской полуночи, 1000 с}
  T2:=0.001*SEC+0.06*MIN+3.6*HR;
  {Считал баллистический коэффициент КА1}
  readln(FL,BK1);
  {Считал начальные условия КА1 (гск, м, м/с)}
  for I:=1 to 6 do readln(FL,X1[I]);
  for I:=1 to 3 do
  begin X1[I]:=1e-6*X1[I]; X1[I+3]:=1e-3*X1[I+3] end;
  {Вычислил звездное время}
  D:=(JD2-2451545.5)/36525.0;
  W:=pi/43.2*(1.002737909350795+D*(5.9006e-11-D*5.9e-15));
  SZ0:=24110.54841+D*(8640184.812866+D*(0.093104-D*6.2e-6));
  SZ0:=SZ0/86400.0; SZ0:=2*pi*(SZ0-trunc(SZ0))+W*T2;
  CA:=cos(SZ0); SA:=sin(SZ0);
  {Перешел к абс. сист.коор.}
  X:=X1[1]*CA-X1[2]*SA; Y:=X1[1]*SA+X1[2]*CA;
  X0[1]:=X; X0[2]:=Y;
  X:=X1[4]*CA-X1[5]*SA-W*X1[2]; Y:=X1[4]*SA+X1[5]*CA+W*X1[1];
  X1[4]:=X; X1[5]:=Y;
  {Считал ковариационную матрицу КА1 (по строкам))}
  for I:=1 to 6 do for J:=1 to 6 do readln(FL,CM1[I,J]);
  close(FL);
  writeln('ka1: JD, T, bc ',JD2,' ',T2:15,' ',BK1:10);
  for I:=1 to 6 do writeln(I,' ',X1[I]:15);

  VARI:=true;
  if JD2+T2>JD1+T1 then
  begin
    CONT(JD1,T1,86.4*(JD2-JD1)+T2-T1,X0,CM0);
    JD0:=JD2; T0:=T2
  end
  else
  begin
    CONT(JD2,T2,86.4*(JD1-JD2)+T1-T2,X1,CM1);
    JD0:=JD1; T0:=T1
  end;
  {Вычислил звездное время и угловую скорость Земли}
  D:=(JD0-2451545-0.5)/36525.0;
  W:=pi/43.2*(1.002737909350795+D*(5.9006e-11-D*5.9e-15));
  SZ0:=24110.54841+D*(8640184.812866+D*(0.093104-D*6.2e-6));
  SZ0:=SZ0/86400.0; SZ0:=2*pi*(SZ0-trunc(SZ0))+W*T0;
  {Задал начальные условия совместного интегрирования}
  for I:=1 to 6 do begin Y0[I]:=X0[I]; Y0[I+6]:=X1[I] end;
  {Считал параметры интегрирования}
  assign(FL,'parint.pas');
  reset(FL);
  readln(FL,P[2],P[3],P[4],S0);
  writeln('P[2..4],S0 ',P[2]:11,' ',P[3]:11,' ',P[4]:11,' ',S0:11);
  readln;
  P[2]:=P[2]*86.4; P[5]:=0;
  close(FL)
end;
