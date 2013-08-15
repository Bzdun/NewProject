function ROATMST(X,Y,Z: real): real;
const R0=6378140; ALFA=0.00335289187;
var H,H0,A,K1,K2,R: real;

begin
  R:=sqrt(sqr(X)+sqr(Y)+sqr(Z));
  H:=1E6*R-R0*(1-ALFA*sqr(Z/R));
  if H<3E5 then
  begin H0:=1.5E5; A:=0.2173E-9; K1:=0.8004E-10; K2:=0.3734E-4 end
  else
  begin
    if H<6E5 then
    begin H0:=3E5; A:=0.4861E-11; K1:=0.7111E-11; K2:=0.1547E-4 end
    else
    begin
      if H<9E5 then
      begin H0:=6E5; A:=0.8904E-13; K1:=0.1831E-11; K2:=0.9275E-5 end
      else begin H0:=9E5; A:=0.6497E-14; K1:=0; K2:=0.954E-5 end
    end
  end;
  H:=H-H0;
  ROATMST:=9.80665*A*exp((K1*H-K2)*H)
end;
