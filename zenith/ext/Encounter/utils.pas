procedure AbsToGsk(aX,aY,aZ: real; W,T,SZ0: real; var gX, gY, gZ: real);
var
  Y12,Y3,SZ,CA,SA:real;
begin
  {Перевод из абсолютной системы координат в гринвичскую}
  SZ:=SZ0+W*T;
  SA:=sin(SZ);
  CA:=cos(SZ);

  gX:=aX*CA+aY*SA;
  gY:=-aX*SA+aY*CA;
  gZ:=aZ;
end;

