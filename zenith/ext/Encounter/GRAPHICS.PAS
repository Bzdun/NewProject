{
浜様様様様様様様様様様様様様様様様様様様様様様様様様様様様�
�    23.02.93                                             �
�                  M.M. Komarov                           �
�                                                         �
藩様様様様様様様様様様様様様様様様様様様様様様様様様様様様�
}
unit Graphics;
interface
uses    Graph;
type
     real=single;
     Vect2 = record X, Y: real end;

procedure GrInit;
procedure SetMinMax (Xn, Xx, Yn, Yx: real);
procedure MinMaxXY (var A; N: word; var Pmin, Pmax: vect2);
procedure SetMarker (M, C, D: integer);
procedure PlotGraph (var XY; N: integer; AxColor, GrColor: word);
procedure PlotAxis (var Pmin, Pmax: vect2;
                        XA, YA, XB, YB: word;
                        AxColor, NbColor: word;
                        Stx, Lfx, Lsx, Fx1, Fx2,
                        Sty, Lfy, Lsy, Fy1, Fy2: word);
procedure PlotLine (var XY; N: integer; Flag: integer; GrColor: integer);

implementation
type Arp = array [1..50000] of vect2;

var
  MinMaxFlag: boolean;
  MarkerType, MarkerColor, MarkerD: integer;
  SCR_REG: record
              SXMIN, SXMAX, SYMIN, SYMAX: integer; {SCREAN;}
           end;
  XYRATIO, XMIN, XMAX, YMIN, YMAX: real;
  Xasp, Yasp: word;

{様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様�}

procedure Set_Min_Max_Flag (Flag: boolean);
begin
MinMaxFlag:=Flag;
end;
{様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様�}
procedure SetMinMax (Xn, Xx, Yn, Yx: real);
begin
   XMIN:=Xn;
   XMAX:=Xx;
   YMIN:=Yn;
   YMAX:=Yx;
end;
{様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様�}

function ScreanX (X: real): integer;
var  FX: real;
begin
with SCR_REG do begin
FX:=(XMAX-XMIN)/(SXMAX-SXMIN)*XYRATIO;
if X < XMIN then
   begin  ScreanX:=SXMIN end
   else if X > XMAX
           then
           begin ScreanX:=SXMAX end
           else  ScreanX:=round((X-XMIN)*XYRATIO/FX)+SXMIN;
end;
end;
{様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様�}

function ScreanY (Y: real): integer;
var  FY: real;
begin
with SCR_REG do begin
FY:=(YMAX-YMIN)/(SYMAX-SYMIN);
if Y < YMIN
   then
   begin  ScreanY:=SYMAX end
   else if Y > YMAX
           then
           begin  ScreanY:=SYMIN end
           else  ScreanY:=SYMAX-round((Y-YMIN)/FY);
end;
end;
{様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様�}

procedure ScreanXY (X, Y: real; var SX: integer; SY: integer);
var  FX, FY: real;
begin
with SCR_REG do begin
FX:=(XMAX-XMIN)/(SXMAX-SXMIN)*XYRATIO;
FY:=(YMAX-YMIN)/(SYMAX-SYMIN);
if X < XMIN
   then SX:=SXMIN
   else if X > XMAX
           then  SX:=SXMAX
           else  SX:=round((X-XMIN)*XYRATIO/FX)+SXMIN;
if Y < YMIN
   then SY:=SYMAX-SYMIN
   else if Y > YMAX
           then  SY:=SYMIN
           else  SY:=SYMAX-round((Y-YMIN)/FY);
end;
end;
{様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様�}
procedure SetMarker (M, C, D: integer);
begin
  MarkerType:= M;
  MarkerColor:= C;
  MarkerD := D;
end;
{様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様�}
procedure DrawMarker(var XY; N: integer);
var J, Sx, Sy: integer;
begin
SetColor(MarkerColor);
   for J:=1 to N do
   begin
       Sx := ScreanX(Arp(XY)[J].X);
       Sy := ScreanY(Arp(XY)[J].Y);
       MoveTo (Sx, Sy);
       case MarkerType of
       1: begin
             SetFillStyle(1, MarkerColor);
             Bar(Sx-MarkerD, Sy-MarkerD, Sx+MarkerD, Sy+MarkerD);
          end;
       2: begin
             SetFillStyle(1, 0);
             Bar(Sx-MarkerD, Sy-MarkerD, Sx+MarkerD, Sy+MarkerD);
             Rectangle(Sx-MarkerD, Sy-MarkerD, Sx+MarkerD, Sy+MarkerD);
          end;
       3: begin
             Circle(Sx, Sy, MarkerD);
             SetFillStyle(1, MarkerColor);
             FloodFill(Sx, Sy, MarkerColor);
          end;
       4: begin
             SetFillStyle(1, 0);
             Circle(Sx, Sy, MarkerD);
             FloodFill(Sx, Sy, MarkerColor);
          end;
       5: begin
             Line(Sx-MarkerD, Sy, Sx+MarkerD, Sy);
          end;
       6: begin
             Line(Sx, Sy-MarkerD, Sx, Sy+MarkerD);
          end;
       7: begin
             Line(Sx-MarkerD, Sy, Sx+MarkerD, Sy);
             Line(Sx, Sy-MarkerD, Sx, Sy+MarkerD);
          end;
       8: begin
             Line(Sx-MarkerD-1, Sy, Sx+MarkerD+1, Sy);
             Line(Sx-MarkerD+1, Sy-MarkerD, Sx+MarkerD-1, Sy+MarkerD);
             Line(Sx+MarkerD-1, Sy-MarkerD, Sx-MarkerD+1, Sy+MarkerD);
          end;
       9: begin
             Line(Sx-MarkerD-1, Sy, Sx+MarkerD+1, Sy);
             Line(Sx-MarkerD, Sy-MarkerD, Sx+MarkerD, Sy+MarkerD);
             Line(Sx+MarkerD, Sy-MarkerD, Sx-MarkerD, Sy+MarkerD);
          end;
      end; {case}
   end;
end;
{様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様�}
procedure MinMaxXY (var A; N: word; var Pmin, Pmax: vect2);
var
    J: integer;
begin
Pmin.X:=Arp(A)[1].X;  Pmax.X:=Arp(A)[1].X;
Pmin.Y:=Arp(A)[1].Y;  Pmax.Y:=Arp(A)[1].Y;
for J:=2 to N do
   begin
   if Arp(A)[J].X < Pmin.X then Pmin.X:=Arp(A)[J].X
                  else if Arp(A)[J].X > Pmax.X then Pmax.X:=Arp(A)[J].X;
   if Arp(A)[J].Y < Pmin.Y then Pmin.Y:=Arp(A)[J].Y
                  else if Arp(A)[J].Y > Pmax.Y then Pmax.Y:=Arp(A)[J].Y;
   end;
end;
{様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様�}

procedure GrInit;
{ Initialize graphics and report any errors that may occur }
var
  GraphDriver : integer;  { The Graphics device driver }
  GraphMode   : integer;  { The Graphics mode value }
  ErrorCode   : integer;  { Reports any graphics errors }
begin
  {GraphMode:=VGAHi;}

  GraphDriver := Detect;                 { use autodetection }
  {GraphDriver := VGA;}                                        {--!!!--}
  InitGraph(GraphDriver, GraphMode,'c:\pp\BGI'); { activate graphics }
  ErrorCode := GraphResult;               { error? }
  if ErrorCode <> grOk then
  begin
    Writeln('!!!  GrInit / Graphics error: ', GraphErrorMsg(ErrorCode));
    Halt(1);
  end;
  GetAspectRatio (Xasp, Yasp);
  XYRatio:=Yasp/Xasp;
  with SCR_REG do begin
     SXMIN:=0;  SXMAX:=GetMaxX;
     SYMIN:=0;  SYMAX:=GetMaxY;
     end;
 MinMaxFlag:=true;
end; { Initialize }
{様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様�}

procedure PlotGraph (var XY; N: integer; AxColor, GrColor: word);
var J: integer;
    Pmin, Pmax: vect2;
    S: string;
begin
   MinMaxXY(XY, N, Pmin, Pmax);
   XMIN:=PMin.X;
   XMAX:=PMax.X;
   YMIN:=PMin.Y;
   YMAX:=PMax.Y;
S:='Xmin:=';
OutTextXY(SCR_REG.SXmin+10, SCR_REG.SYmax+10, S);
Str (Xmin:10:4, S);
OutTextXY(SCR_REG.SXmin+50, SCR_REG.SYmax+10, S);
S:='Xmax:=';
OutTextXY(SCR_REG.SXmin+150, SCR_REG.SYmax+10, S);
Str (Xmax:10:4, S);
OutTextXY(SCR_REG.SXmin+200, SCR_REG.SYmax+10, S);
S:='Ymin:=';
OutTextXY(SCR_REG.SXmin+10, SCR_REG.SYmax+20, S);
Str (Ymin:10:4, S);
OutTextXY(SCR_REG.SXmin+50, SCR_REG.SYmax+20, S);
S:='Ymax:=';
OutTextXY(SCR_REG.SXmin+150, SCR_REG.SYmax+20, S);
Str (Ymax:10:4, S);
OutTextXY(SCR_REG.SXmin+200, SCR_REG.SYmax+20, S);
SetColor(AxColor);
with SCR_REG do
   begin
   Line(SXmin, SYmax, SXmax, SYmax);
   Line(SXmin, SYmin, SXmin, SYmax);
   Line(SXmin, SYmin, SXmin-2, SYmin+4);
   Line(SXmin, SYmin, SXmin+2, SYmin+4);
   Line(SXmax, SYmax, SXmax-4, SYmax-2);
   Line(SXmax, SYmax, SXmax-4, SYmax+2);
   for J:=1 to 9 do
   begin
   Line(SXmin+round((SXmax-SXmin)*J/10), SYmax+2,
        SXmin+round((SXmax-SXmin)*J/10), SYmax-2);
   end;
   for J:=1 to 9 do
   begin
   Line(SXmin-2, SYmin+round((SYmax-SYmin)*J/10),
        SXmin+2, SYmin+round((SYmax-SYmin)*J/10));
   end;
   end;
SetColor(GrColor);
MoveTo (ScreanX(Arp(XY)[1].X), ScreanY(Arp(XY)[1].Y));
for J:=2 to N do
begin
    LineTo (ScreanX(Arp(XY)[J].X), ScreanY(Arp(XY)[J].Y));
end;
end;
{様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様�}
procedure PlotAxis (var Pmin, Pmax: vect2;
                        XA, YA, XB, YB: word;
                        AxColor, NbColor: word;
                        Stx, Lfx, Lsx, Fx1, Fx2,
                        Sty, Lfy, Lsy, Fy1, Fy2: word);
{
   Pmin -    �┃━��讚襯 Л�腑��� �� X � Y
   Pmax -    ���瓱���讚襯 Л�腑��� �� X � Y
   XA, YA -  ���爐┃�硅 ▲爿ィ� �ア��� 磽�� ���筥�� (!) �� 蹣����
   XB, YB -  ���爐┃�硅 �└�ィ� ���〓�� 磽�� ���筥�� (!) �� 蹣����
   AxColor - 罌モ 爬甌���� �瓮�
   NbColor - 罌モ 爬甌����� Л�腑�┤ ぅ�キ┤ �� �瘴�
   Stx -     腮甄� ぅ�キ┤ �� �瓱 X
   Lfx -     ���ム �ム〓�� ぅ�キ�� � ��皰牀�� ����瓱矚� Л�腑�┘ ぅ�キ�� �� �瓱 X
   Lsx -     ��� ���メキ�� Л�腑�┤ �� �瓱 X
   Fx1 -     腮甄� ��Ж腮� ��� Л�腑�┘ ぅ�キ�� �� �瓱 X
   Fx2 -     腮甄� 罔籥 ��甄� ���閧�� � Л�腑�┬ ぅ�キ�� �� �瓱 X
   Sty, Lfy, Lsy, Fy1, Fy2 - 皰Ε �� �瓱 Y
}
var J, K: integer;
    S: string;
    R: real;
begin
   XMIN:=PMin.X;
   XMAX:=PMax.X;
   YMIN:=PMin.Y;
   YMAX:=PMax.Y;
SetColor(AxColor);
with SCR_REG do
   begin
   SXMIN:=XA;
   SXMAX:=XB;
   SYMIN:=YA;
   SYMAX:=YB;
   Line(SXmin, SYmax, SXmax, SYmax);
   Line(SXmin, SYmin, SXmin, SYmax);
   Line(SXmin, SYmin, SXmin-2, SYmin+4);
   Line(SXmin, SYmin, SXmin+2, SYmin+4);
   Line(SXmax, SYmax, SXmax-4, SYmax-2);
   Line(SXmax, SYmax, SXmax-4, SYmax+2);
   K:=Lfx;
   SetTextStyle(DefaultFont, HorizDir, 1);
   for J:=0 to Stx do
   begin
   if (J <> 0) and (J <> Stx) then
        Line(SXmin+round((SXmax-SXmin)*J/Stx), SYmax+2,
             SXmin+round((SXmax-SXmin)*J/Stx), SYmax-2);
   if J=K then
      begin
         R:=(Xmax-Xmin)/Stx*J+Xmin;
         Str (R:Fx1:Fx2, S);
         SetColor(NbColor);
         OutTextXY(SXmin-Fx1*4+((SXmax-SXmin)*J) div Stx, SCR_REG.SYmax+10, S);
         SetColor(AxColor);
         K:=K+Lsx;
      end;
   end;
   K:=Lfy;
   for J:=0 to Sty do
   begin
   if (J <> 0) and (J <> Sty) then
      Line(SXmin-2, SYmin+round((SYmax-SYmin)*J/Sty),
           SXmin+2, SYmin+round((SYmax-SYmin)*J/Sty));
   if J=K then
      begin
         R:=(Ymax-Ymin)/Sty*J+Ymin;
         Str (R:Fy1:Fy2, S);
         SetColor(NbColor);
         OutTextXY(SXmin-Fy1*8-4,SYmax-4-((SYmax-SYmin)*J) div Sty, S);
         SetColor(AxColor);
         K:=K+Lsy;
      end;
   end;
   end;
end;
{様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様�}
procedure PlotLine (var XY; N: integer; Flag: integer; GrColor: integer);
var J, K, M: integer;
    Y, YL, YG, YG2, X, DX, A: real;
begin
SetColor(GrColor);
case Flag of
1: begin
   MoveTo (ScreanX(Arp(XY)[1].X), ScreanY(Arp(XY)[1].Y));
   for J:=2 to N do
       LineTo (ScreanX(Arp(XY)[J].X), ScreanY(Arp(XY)[J].Y));
   end; { 1 }
2: begin
   M:=100;
   DX:=(Arp(XY)[N].X-Arp(XY)[1].X)/M;
   X:=Arp(XY)[1].X;
   K:=1;
   MoveTo (ScreanX(Arp(XY)[1].X), ScreanY(Arp(XY)[1].Y));
   for J:=2 to M-1 do
      begin
      X:=X+DX;
      while not( (X >= Arp(XY)[K].X) and (X < Arp(XY)[K+1].X)) do K:=K+1;
      if K > 1 then YL:=Arp(XY)[K-1].Y else YL:=2*Arp(XY)[1].Y-Arp(XY)[2].Y;
      YG:=Arp(XY)[K+1].Y;
      if K <> N-1 then YG2:=Arp(XY)[K+2].Y else YG2:=2*Arp(XY)[N].Y-Arp(XY)[N-1].Y;
      A:=(X-Arp(XY)[K].X)/(Arp(XY)[K+1].X-Arp(XY)[K].X);
      Y:=(1-A)*Arp(XY)[K].Y+A*YG+(1-A)*sqr(1-A)*
         (YG-2*Arp(XY)[K].Y+YL)/6+
         A*sqr(A)*(Arp(XY)[K].Y-2*YG+YG2)/6;
      LineTo (ScreanX(X), ScreanY(Y));
      end;
   LineTo (ScreanX(Arp(XY)[N].X), ScreanY(Arp(XY)[N].Y));
end; {2}
-9..-1: begin
   SetMarker(Abs(Flag), GrColor, 2);
   DrawMarker(XY, N);
   SetMarker(0, GrColor, 2);
   end; { -9..-1 }
end; {case }
if MarkerType > 0 then DrawMarker(XY, N);
end;
{様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様様�}
end.
