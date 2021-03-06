procedure RKDP78(N: integer; var Y: MASN; var P: MS5;
                PRAVCH: INTPR; OUTP: INTOU);
var BOL: boolean;
    I,J,K: integer;
    S,T,H,TE,TE1: real;
    B,D: array[1..13] of extended;
    C: array[2..13] of real;
    A: array[2..13,1..12] of extended;
    F: array[1..14] of MASN;
  begin
  if (P[1]=P[2]) or (P[3]<=0) or (P[4]<=0) or (P[5]<>0) or (N<=0) then
  begin
    writeln('bad initial information');
    for I:=1 to 5 do write('  ',P[I]:10);
    writeln(' ');
    exit
  end;

  for I:=2 to 13 do
  begin B[I]:=0; D[I]:=0; for J:=1 to 12 do A[I,J]:=0 end;

  B[1]:=1400.5451/33548.0064; B[6]:=-5923.8493/106827.7825;
  B[7]:=18160.6767/75886.7731; B[8]:=56129.2985/79784.5732;
  B[9]:=-104189.143/137134.3529; B[10]:=76041.7239/115116.5299;
  B[11]:=11882.0643/75113.8087; B[12]:=-52874.7749/222060.717; B[13]:=1/4;

  D[1]:=1345.1932/45517.6623; D[6]:=-80871.9846/97600.0145;
  D[7]:=175700.4468/564515.9321; D[8]:=65604.5339/26589.1186;
  D[9]:=-386757.4721/151851.7206; D[10]:=46588.5868/32273.6535;
  D[11]:=5301.1238/66751.6719; D[12]:=2/45;

  C[2]:=1/18; C[3]:=1/12; C[4]:=1/8; C[5]:=5/16; C[6]:=3/8;
  C[7]:=59/400; C[8]:=93/200; C[9]:=549002.3248/971916.9821;
  C[10]:=13/20; C[11]:=120114.6811/129901.9798; C[12]:=1; C[13]:=1;

  A[2,1]:=1.0/18.0;
  A[3,1]:=1.0/48.0; A[3,2]:=1.0/16.0;
  A[4,1]:=1.0/32.0; A[4,3]:=3.0/32.0;
  A[5,1]:=5.0/16.0; A[5,3]:=-75.0/64.0; A[5,4]:=-A[5,3];
  A[6,1]:=3.0/80.0; A[6,4]:=3.0/16.0; A[6,5]:=3.0/20.0;
  A[7,1]:=2944.3841/61456.3906; A[7,4]:=7773.6538/69253.8347;
  A[7,5]:=-2869.3883/112500.0; A[7,6]:=2312.4283/180000.0;
  A[8,1]:=1601.6141/94669.2911; A[8,4]:=6156.418/15873.2637;
  A[8,5]:=2278.9713/63344.5777; A[8,6]:=54581.5736/277105.7229;
  A[8,7]:=-18019.3667/104330.7555;
  A[9,1]:=3963.2708/57359.1083; A[9,4]:=-43363.6366/68370.1615;
  A[9,5]:=-42173.9975/261629.2301; A[9,6]:=10030.2831/72342.3059;
  A[9,7]:=79020.4164/83981.3087; A[9,8]:=80063.531/378307.1287;
  A[10,1]:=24612.1993/134084.7787; A[10,4]:=-3769504.2795/1526876.6246;
  A[10,5]:=-30912.1744/106122.7803; A[10,6]:=-1299.2083/49076.6935;
  A[10,7]:=600594.3493/210894.7869; A[10,8]:=39300.6217/139667.3457;
  A[10,9]:=12387.2331/100102.9789;
  A[11,1]:=-102846.8189/84618.0014; A[11,4]:=847823.5783/50851.2852;
  A[11,5]:=131172.9495/143242.2823; A[11,6]:=-1030412.9995/170130.4382;
  A[11,7]:=-4877792.5059/304793.956; A[11,8]:=1533672.6248/103282.4649;
  A[11,9]:=-4544286.8181/339846.7696; A[11,10]:=306599.3473/59717.2653;
  A[12,1]:=18589.2177/71811.6043; A[12,4]:=-318509.4517/66710.7341;
  A[12,5]:=-47775.5414/109805.3517; A[12,6]:=-70363.5378/23073.9211;
  A[12,7]:=573156.6787/102754.5527; A[12,8]:=523286.6602/85006.6563;
  A[12,9]:=-409366.4535/80868.8257; A[12,10]:=396213.7247/180595.7418;
  A[12,11]:=6568.6358/48791.0083;
  A[13,1]:=40386.3854/49106.3109; A[13,4]:=-506849.2393/43474.0067;
  A[13,5]:=-41142.1997/54304.3805; A[13,6]:=65278.3627/91429.6604;
  A[13,7]:=1117396.2825/92532.0556; A[13,8]:=-1315899.0841/618472.7034;
  A[13,9]:=393664.7629/197804.968; A[13,10]:=-16052.8059/68517.8525;
  A[13,11]:=24863.8103/141353.1060; A[13,12]:=0;

(*-------------------- Cheking --------------------------------------------
  for I:=2 to 13 do
  begin
    S:=-C[I]; for J:=1 to I-1 do S:=S+A[I,J];
    writeln(S:20)
  end;
  readln;
  for I:=3 to 13 do
  begin
    S:=-sqr(C[I])/2;
    for J:=2 to I-1 do S:=S+A[I,J]*C[J];
    writeln(S:20)
  end;
  readln;
  for I:=3 to 13 do
  begin
    S:=-sqr(C[I])*C[I]/3;
    for J:=2 to I-1 do S:=S+A[I,J]*sqr(C[J]);
    writeln(S:20)
  end;
  readln;
  for I:=6 to 13 do
  begin
    S:=-sqr(sqr(C[I]))/4;
    for J:=2 to I-1 do S:=S+A[I,J]*sqr(C[J])*C[J];
    writeln(S:20)
  end;
  readln;
  for I:=1 to 8 do
  begin
    S:=-1/(I+1);
    for J:=2 to 13 do S:=S+B[J]*exp(I*ln(C[J]));
    writeln(S:20)
  end;
  readln;
  for I:=1 to 8 do
  begin
    S:=-1/(I+1);
    for J:=2 to 13 do S:=S+D[J]*exp(I*ln(C[J]));
    writeln(S:20)
  end;
  readln;

  S:=-1; for I:=1 to 13 do S:=S+B[I];
  writeln(S:20);
  S:=-1; for I:=1 to 13 do S:=S+D[I];
  writeln(S:20);
  readln;
  exit;
------------------------------------------------------------------------ *)
  for I:=1 to 13 do D[I]:=B[I]-D[I];

  T:=P[1];
  if P[2]<T then H:=-P[3] else H:=P[3];
  while abs(P[2]-T)>=1E-11 do
  begin
    if ((P[2]>P[1]) and (H>P[2]-T)) or ((P[2]<P[1]) and (H<P[2]-T)) then
      H:=P[2]-T;
    OUTP(T,H,Y,P);
    if P[5]<>0 then begin P[2]:=T; exit end;
    PRAVCH(T,Y,F[1]);
    repeat
      for I:=2 to 13 do
      begin
        for K:=1 to N do
        begin
          S:=0;
          for J:=1 to I-1 do S:=S+A[I,J]*F[J,K];
          F[14,K]:=Y[K]+H*S;
        end;
        PRAVCH(T+H*C[I],F[14],F[I]);
      end;
      TE:=0;
      for K:=1 to N do
      begin
        TE1:=0; for J:=1 to 13 do TE1:=TE1+D[J]*F[J,K];
        TE1:=abs(TE1*H); if TE1>TE then TE:=TE1
      end;
      if TE>P[4] then begin BOL:=false; H:=0.7*H end
        else BOL:=true
    until BOL;
    for K:=1 to N do
    begin
      S:=0;
      for I:=1 to 13 do S:=S+B[I]*F[I,K];
      Y[K]:=Y[K]+S*H
    end;
    T:=T+H; H:=0.9*H*exp(0.125*ln(P[4]/TE))
  end;
  OUTP(T,H,Y,P)
end;
