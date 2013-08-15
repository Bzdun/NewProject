procedure REALFFT(var DATA: GLDARRAY; N,ISIGN: integer);
{ N - a power of 2 = a half of number of points,
  ISIGN=-1 for direct transform,
  ISIGN=1 for inverse transform,
  direct*inverse = N }
var Q,R0,R1,R2,S0,S1,S2,THETA: real;
    A,B,C1,C2,H1R,H1I,H2R,H2I: single;
    I,I1,I2,I3,I4: integer;
begin
  THETA:=-pi/N; C1:=0.5;
  if ISIGN=-1 then begin C2:=-0.5; FOUR2(DATA,N,-1) end
    else begin C2:=0.5; THETA:=-THETA end;
  R1:=cos(THETA); S1:=sin(THETA);
  R0:=1; S0:=0; Q:=2*R1;
  for I:=2 to (N div 2)+1 do
  begin
    I1:=I+I-1; I2:=I1+1; I3:=N+N+3-I2; I4:=I3+1;
    H1R:= C1*(DATA[I1]+DATA[I3]); H1I:=C1*(DATA[I2]-DATA[I4]);
    H2R:=-C2*(DATA[I2]+DATA[I4]); H2I:=C2*(DATA[I1]-DATA[I3]);
    A:=R1*H2R-S1*H2I; B:=R1*H2I+S1*H2R;
    DATA[I1]:=H1R+A; DATA[I2]:=B+H1I;
    DATA[I3]:=H1R-A; DATA[I4]:=B-H1I;
    R2:=Q*R1-R0; S2:=Q*S1-S0; R0:=R1; R1:=R2; S0:=S1; S1:=S2
  end;
  if ISIGN=-1 then
  begin
    H1R:=DATA[1]; DATA[1]:=H1R+DATA[2]; DATA[2]:=H1R-DATA[2]
  end
  else
  begin
    H1R:=DATA[1]; DATA[1]:=C1*(H1R+DATA[2]);
    DATA[2]:=C1*(H1R-DATA[2]); FOUR2(DATA,N,1)
  end
 end;
