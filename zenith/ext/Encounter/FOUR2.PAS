procedure FOUR2(var DATA: GLDARRAY; NN,ISIGN: integer);
{ NN - a power of 2 = number of points,
  ISIGN=-1 for direct transform,
  ISIGN=1 for inverse transform,
  direct*inverse = NN }
var II,JJ,N,NMAX,M,J,ISTEP,I: integer;
    C0,C1,C2,S0,S1,S2,Q,THETA: real;
    TEMPR,TEMPI: single;
begin
  N:=NN+NN; J:=1;
  for II:=1 to NN do
  begin
    I:=II+II-1;
    if J>I then
    begin
      TEMPR:=DATA[J]; TEMPI:=DATA[J+1]; DATA[J]:=DATA[I];
      DATA[J+1]:=DATA[I+1]; DATA[I]:=TEMPR; DATA[I+1]:=TEMPI
    end;
    M:=N div 2;
    while (M>=2) and (J>M) do
    begin
      J:=J-M; M:=M div 2
    end;
    J:=J+M
  end;
  NMAX:=2;
  while N>NMAX do
  begin
    ISTEP:=NMAX+NMAX; THETA:=2*pi/(ISIGN*NMAX);
    C1:=cos(THETA); S1:=sin(THETA);
    C0:=1; S0:=0; Q:=2*C1;
    for II:=1 to (NMAX div 2) do
    begin
      M:=II+II-1; I:=M-ISTEP;
      for JJ:=0 to ((N-M) div ISTEP) do
      begin
        I:=I+ISTEP; J:=I+NMAX;
        TEMPR:=C0*DATA[J]-S0*DATA[J+1];
        TEMPI:=C0*DATA[J+1]+S0*DATA[J];
        DATA[J]:=DATA[I]-TEMPR; DATA[J+1]:=DATA[I+1]-TEMPI;
        DATA[I]:=DATA[I]+TEMPR; DATA[I+1]:=DATA[I+1]+TEMPI
      end;
      C2:=Q*C1-C0; S2:=Q*S1-S0; C0:=C1; C1:=C2; S0:=S1; S1:=S2;
    end;
    NMAX:=ISTEP
  end
end;
