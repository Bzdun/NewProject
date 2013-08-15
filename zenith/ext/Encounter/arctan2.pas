function ATAN2(S,C: real): real;
var X: real;
begin
  if C=0 then
  begin
    if S<0 then X:=-0.5*pi else X:=0.5*pi
  end
  else
  begin
    X:=arctan(S/C);
    if C<0 then
    begin
      if S>0 then X:=X+pi else X:=X-pi
    end
  end;
  ATAN2:=X
end;
