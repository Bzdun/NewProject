procedure RKV56(N: integer; var Y: MASN; var F: MASN9; var P: MS5;
                PRAVCH: INTPR; OUTP: INTOU);
type MAS8=array[1..8] of real;
     MAS77=array[2..8,1..7] of real;
const B: MAS8=(3/40, 0, 875/2244, 23/72, 264/1955, 0, 125/11592, 43/616);
      C: MAS8=(0, 1/6, 4/15, 2/3, 5/6, 1, 1/15, 1);
      D: MAS8=(-1/160, 0, -125/17952, 1/144, -12/1955, -3/44, 125/11592, 43/616);
      A: MAS77=((1/6, 0, 0, 0, 0, 0, 0), (4/75, 16/75, 0, 0, 0, 0, 0),
                (5/6, -8/3, 5/2, 0, 0, 0, 0),
                (-165/64, 55/6, -425/64, 85/96, 0, 0, 0),
                (12/5, -8, 4015/612, -11/36, 88/255, 0, 0),
                (-8263/15000, 124/75, -643/680, -81/250, 2484/10625, 0, 0),
                (3501/1720, -300/43, 297275/52632, -319/2322, 24068/84065, 0, 3850/26703));
      TOL=0.117649;
var I,J,K: integer;
    S,T,H,H1,TE,TE1: real;
    BOL: boolean;
begin
  if (P[1]=P[2]) or (P[3]<=0) or (P[4]<=0) or (P[5]<>0) or (N<=0) then
  begin
    writeln('bad initial information');
    for I:=1 to 5 do write('  ',P[I]:10);
    writeln(' ');
    exit
  end;
  T:=P[1];
  if P[2]<T then H:=-P[3] else H:=P[3]; H1:=0;
  while abs(P[2]-T)>=1E-11 do
  begin
    if ((P[2]>P[1]) and (H>P[2]-T)) or ((P[2]<P[1]) and (H<P[2]-T)) then
      H:=P[2]-T;
    OUTP(T,H1,Y,P,F);
    if P[5]<>0 then begin P[2]:=T; exit end;
    PRAVCH(T,Y,F[1]);
    repeat
      for I:=2 to 8 do
      begin
        for K:=1 to N do
        begin
          S:=0;
          for J:=1 to I-1 do S:=S+A[I,J]*F[J,K];
          F[9,K]:=Y[K]+H*S;
        end;
        PRAVCH(T+H*C[I],F[9],F[I]);
      end;
      TE:=0;
      for K:=1 to N do
      begin
        TE1:=0; for J:=1 to 8 do TE1:=TE1+D[J]*F[J,K];
        TE1:=abs(H*TE1); if TE1>TE then TE:=TE1
      end;
      if TE>P[4] then begin BOL:=false; H:=0.7*H end else BOL:=true
    until BOL;
    for K:=1 to N do
    begin
      S:=0;
      for I:=1 to 8 do S:=S+B[I]*F[I,K];
      F[9,K]:=Y[K]; Y[K]:=Y[K]+S*H
    end;
    T:=T+H; H1:=H;
    if TE<(TOL*P[4]) then H:=H/0.7
  end;
  OUTP(T,H1,Y,P,F)
end;

