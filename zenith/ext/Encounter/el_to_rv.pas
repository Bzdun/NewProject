{
Переход от кеплеровых элементов a,e,E_0,\Omega,i,\omega
к фазовому вектору r_0,v_0
}
procedure EL_RV(var X,Y: MS6);
var CA,SA,CU,SU,CI,SI,CE,SE,PX,PY,PZ,QX,QY,QZ,E1,
    RP,RQ,V0,VP,VQ: real;
begin
  CE:=cos(X[3]); SE:=sin(X[3]);
  CU:=cos(X[4]); SU:=sin(X[4]);
  CI:=cos(X[5]); SI:=sin(X[5]);
  CA:=cos(X[6]); SA:=sin(X[6]);
  PX:=CA*CU-SA*CI*SU; QX:=-SA*CU-CA*CI*SU;
  PY:=CA*SU+SA*CI*CU; QY:=-SA*SU+CA*CI*CU;
  PZ:=SA*SI;          QZ:= CA*SI;
  RP:=X[1]*(CE-X[2]); E1:=sqrt(1-sqr(X[2]));
  RQ:=X[1]*E1*SE; V0:=sqrt(GC/X[1]);
  VP:=-V0*SE/(1-X[2]*CE); VQ:=V0*E1*CE/(1-X[2]*CE);
  Y[1]:=RP*PX+RQ*QX; Y[4]:= VP*PX+VQ*QX;
  Y[2]:=RP*PY+RQ*QY; Y[5]:= VP*PY+VQ*QY;
  Y[3]:=RP*PZ+RQ*QZ; Y[6]:= VP*PZ+VQ*QZ
end;
