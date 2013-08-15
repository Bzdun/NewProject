{
Переход от фазового вектора r_0,v_0 к
кеплеровым элементам a,e,E_0,\Omega,i,\omega
}
procedure RV_EL(var Y,X: MS6);
var CX,CY,CZ,C,FX,FY,FZ,F,R,CE,SE,Z,
    PX,PY,PZ,QX,QY,QZ,RX,RY,RZ,A,E,E0: real;
begin
  CX:=Y[2]*Y[6]-Y[3]*Y[5];
  CY:=Y[3]*Y[4]-Y[1]*Y[6];
  CZ:=Y[1]*Y[5]-Y[2]*Y[4];
  C:=sqrt(sqr(CX)+sqr(CY)+sqr(CZ));
  RX:=CX/C; RY:=CY/C; RZ:=CZ/C;
  R:=sqrt(sqr(Y[1])+sqr(Y[2])+sqr(Y[3]));
  Z:=GC/R;
  FX:=Y[5]*CZ-Y[6]*CY-Z*Y[1];
  FY:=Y[6]*CX-Y[4]*CZ-Z*Y[2];
  FZ:=Y[4]*CY-Y[5]*CX-Z*Y[3];
  F:=sqrt(sqr(FX)+sqr(FY)+sqr(FZ));
  PX:=FX/F; PY:=FY/F; PZ:=FZ/F;
  QX:=RY*PZ-RZ*PY; QY:=RZ*PX-RX*PZ; QZ:=RX*PY-RY*PX;
  A:=sqr(Y[4])+sqr(Y[5])+sqr(Y[6])-2*Z;
  A:=-GC/A; E:=F/GC;
  CE:=(PX*Y[1]+PY*Y[2]+PZ*Y[3])/A+E;
  SE:=(QX*Y[1]+QY*Y[2]+QZ*Y[3])/A/sqrt(1-sqr(E));
  X[1]:=A; X[2]:=E;
  X[3]:=ATAN2(SE,CE);
  X[4]:=ATAN2(RX,-RY);
  X[5]:=ATAN2(sqrt(sqr(RX)+sqr(RY)),RZ);
  X[6]:=ATAN2(PZ,QZ)
end;
