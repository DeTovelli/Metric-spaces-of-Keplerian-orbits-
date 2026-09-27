Program D_criteria;

{$APPTYPE CONSOLE}

{for verhin, metod delenia+}

uses
  Math,
  SysUtils;

Label MET1,AGAIN;

Const N=85641074;  {���������� ������������ ������}
    {IAU 1994}
    KV=1/29.7846918325927450;  {��� �������� �������� � ��/� ��� mu=1}
    rad=Pi/180.0;
    au_km=1.49597870691000015E+8; {1 au in km; DE406}
    day_sec = 86400;  {����� � �������}
    KMAU=1/au_km;     {�� � ��}
    GM=1.32712440018E+11;  {km3*s-2}
    MU=2.9591220828559110225E-4; {= k**2; k = 0.01720209895 [(AU^3 d^-2])^1/2}

var

c,c1,kn,MU_red,n_,JD,JDi:Extended;
a,e,i,om,w,omega,M,cv,sv,v,EE,r,r1,r2:Extended;
mass,rpart,dens,RC, beta: Extended;
AI,EI,II,OMI,VI,WI,Mi,EEi,se,ce: Extended;
RE,XE,YE,ZE,VE,VXE,VYE,VZE,VV,VXC,VYC,VZC: Extended;
TJD0,TJD,x,y,z,xx,yy,zz: Extended;
Cx,Cy,Cz,A_m: Extended;
JD1,t1,a1,e1,i1,om1,w1,M1,q1,p1,JD2,t2,a2,e2,i2,om2,w2,M2,q2,p2: Extended;
alpha,F,CoT,SiT,T,coT1,KSI1,KSI2,LXc,LYc,LZc,LX,LY,LZ,L: Extended;
CosI,CosP,ro1,ro2,ro5,PI_big,D_SH,D_D,Pheta:Extended;
j,jj:integer;
MOID, CurrentDist: Extended;
v1_step, v2_step: Extended;
x1_m, y1_m, z1_m, v_dummy: Extended; // Coordinates of the point on the 1st orbit
x2_m, y2_m, z2_m: Extended;         // Coordinates of the point on the second orbit
// Variables for two-stage MOID optimization
best_v1, best_v2: Extended;
grid_step, opt_step: Extended;
test_v1, test_v2: Extended;
improved: Boolean;

infile, infile1, infile2, outfile, outdebug: text;

{$I kepl1.fnc}
{$I coor.prc}
{$I velor3.prc}
{$I Convert_.prc}
{$I power.fnc}
{$I r4.fnc}
{$I arctg.fnc}

begin {main}

{Assign(infile,'ref.orb'); Reset(infile);
{��������� ����� ���������, ������� � �������� � JD}
{Read(infile, a,e,i,om,w,M,JD);
{��������� ������� �������� ������� � �������}
{i:=i*rad; om:=om*rad; w:=w*rad; M:=M*rad;
JDi:=JD;

Close(infile);}

{Assign(outdebug,'example.out');
Rewrite(outdebug);  }



Assign(outfile,'file.TXT');
Rewrite(outfile);

Writeln(outfile, 'JD t ro1 ro2 ro5 D_SH D_D MOID');
RandSeed:=1;
{Randomize;}

{RC:=1.9;     {km, ������ ����������}
{dens:=3.0;  {g cm-3, ��������� ����������}
{mass:=1.0000; {g, ����� ����������}
{rpart:=power(0.75*mass/(Pi*dens),1.0/3.0);  {cm, ������ ����������}
{beta:=5.8E-5/(rpart*dens); {Fr/Fgr -- ������������ �����������}
{mu_red:=MU*(1-beta); {�������� ����� ������ �� �������� ��������}


{kn:=sqrt(GM/(au_km*au_km*au_km))*day_sec;
n_:=kn/(a*sqrt(a));    {���/���}


{�� ������� �������� ������� ��������}
{EE:=kepl1(M,e);
r:=a*(1-e*cos(EE));
sv:=a*sqrt(1-e*e)*sin(EE)/r;
cv:=a*(cos(EE)-e)/r;
v:=ArcTg(sv,cv);

{������� � ���� ��������� ������}
{j:=0;
Writeln(outfile, a:15:10,e:15:10,i/rad:15:9,om/rad:15:9,w/rad:15:9,M/rad:15:9,
           ' ',Jdi:14:5,'  ',v/rad:5:1,'  ',j:5);


{��� � ���������� � ������}
{alpha:=180*rad;  {����������� ������, 180 ��� }
//TJD0:= 2459586.50;
Assign(infile1,'elem_point1_2003EH1_-7430.dat');
Reset(infile1);
{��������� JD, ����� ���������, ������� � �������� }
Readln(infile1);
Readln(infile1);
Readln(infile1, JD1,t1,a1,e1,i1,om1,w1,M1,q1);
i1:=i1*rad; om1:=om1*rad; w1:=w1*rad; M1:=M1*rad; {��������� ������� �������� ������� � �������}
Writeln(JD1:6:2,t1:6:2,a1:6:2,e1:6:2,i1:6:2,om1:6:2,w1:6:2,M1:6:2,q1:6:2);
//Readln;

Assign(infile2,'elem_point2_2019SG18_-7430.dat');
Reset(infile2);
{��������� ����� ���������, ������� � �������� � JD}
Readln(infile2);
Readln(infile2);
Readln(infile2, JD2,t2,a2,e2,i2,om2,w2,M2,q2);
i2:=i2*rad; om2:=om2*rad; w2:=w2*rad; M2:=M2*rad; {��������� ������� �������� ������� � �������}

Writeln(JD2:6:2,t2:6:2,a2,e2,i2,om2,w2,M2,q2);
//Readln;
L:=1;

while not Eof(infile1) do
begin
Readln(infile1, JD1,t1,a1,e1,i1,om1,w1,M1,q1);
// Converting angles from degrees to radians
i1:=i1*rad; om1:=om1*rad; w1:=w1*rad; M1:=M1*rad; {��������� ������� �������� ������� � �������}
p1:=a1*(1-sqr(e1));

Writeln('TJD1 = ',JD1:6:2,' t1 = ',t1:6:2,t1:6:2);
//Readln;

for j:=1 to N do
 begin

  Readln(infile2, JD2,t2,a2,e2,i2,om2,w2,M2,q2);
  if JD1=JD2 then
   begin
   // Converting angles from degrees to radians
    i2:=i2*rad; om2:=om2*rad; w2:=w2*rad; M2:=M2*rad; {��������� ������� �������� ������� � �������}
    p2:=a2*(1-sqr(e2));
    // Protection against floating-point rounding errors
    CosI:=(Cos(i1)*Cos(i2))+(Sin(i1)*Sin(i2)*Cos(om1-om2));
    CosI:=EnsureRange(CosI, -1.0, 1.0);

    CosP:=Sin(i1)*Sin(i2)*Sin(w1)*Sin(w2)+(Cos(w1)*Cos(w2)+Cos(i1)*Cos(i2)*Sin(w1)*Sin(w2))*Cos(om1-om2)+(Cos(i2)*Cos(w1)*Sin(w2)-Cos(i1)*Sin(w1)*Cos(w1))*Sin(om1-om2);
    CosP:=EnsureRange(CosP, -1.0, 1.0);

    I:=ArcCos(CosI);

    if Abs(om1-om2) > Pi then
      begin
       PI_big:=om2-om1-2*ArcSin(Cos(i2+i1)*Sin((om2-om1)/2)*Sec(I/2));
      end
    else
      begin
       PI_big:=om2-om1+2*ArcSin(Cos(i2+i1)*Sin((om2-om1)/2)*Sec(I/2));
      end;

    //Pheta:=ArcCos(Sin(i1)*Sin(i2)*Sin(w1)+(Cos(w1)*Cos(w2)+Cos(i1)*Cos(i2)*Sin(w1)*Sin(w2))*Cos(om1-om2)+(Cos(i1)*Cos(w1)*Cos(w2)-Cos(i1)*Sin(w1)*Sin(w2))*Sin(om1-om2));
    Pheta:=ArcCos(EnsureRange(Sin(i1)*Sin(i2)*Sin(w1)+(Cos(w1)*Cos(w2)+Cos(i1)*Cos(i2)*Sin(w1)*Sin(w2))*Cos(om1-om2)+(Cos(i1)*Cos(w1)*Cos(w2)-Cos(i1)*Sin(w1)*Sin(w2))*Sin(om1-om2), -1.0, 1.0));

      ro1:=sqrt((1/L)*(p1+p2-2*sqrt(p1*p2)*CosI)+(sqr(e1)+sqr(e2)-2*e1*e2*CosP));
      ro2:=sqrt((1+sqr(e1))*p1+(1+sqr(e2))*p2-2*sqrt(p1*p2)*(CosI+e1*e2*CosP));
      ro5:=sqrt((1+sqr(e1))*p1+(1+sqr(e2))*p2-2*sqrt(p1*p2)*(e1*e2+Cos(i1-i2)));
      D_SH:=sqrt(sqr(q1-q2)+sqr(e1-e2)+4*sqr(Sin(I/2))+sqr(e1+e2)*sqr(Sin(PI_big)));
      D_D:=sqrt(sqr((e1-e2)/(e1+e2))+sqr((q1-q2)/(q1+q2))+sqr(I/Pi)+sqr((e1+e2)/2)*sqr(Pheta/Pi));

      // ======= START OF OPTIMIZED TWO-STAGE MOID CALCULATION BLOCK =======
      MOID := 1e10;
      best_v1 := 0;
      best_v2 := 0;

      // STAGE 1: Coarse Scan
      grid_step := 15.0 * rad;
      v1_step := 0;
      while v1_step < 2 * Pi do
      begin
        COOR(q1, e1, i1, om1, w1, v1_step, MU, re, x1_m, y1_m, z1_m, vxe, vye, vze, v_dummy);
        v2_step := 0;
        while v2_step < 2 * Pi do
        begin
          COOR(q2, e2, i2, om2, w2, v2_step, MU, re, x2_m, y2_m, z2_m, vxe, vye, vze, v_dummy);
          CurrentDist := Sqrt(sqr(x1_m - x2_m) + sqr(y1_m - y2_m) + sqr(z1_m - z2_m));

          if CurrentDist < MOID then
          begin
            MOID := CurrentDist;
            best_v1 := v1_step;
            best_v2 := v2_step;
          end;
          v2_step := v2_step + grid_step;
        end;
        v1_step := v1_step + grid_step;
      end;

      // STAGE 2: Gradient Descent with Angle Normalization
      opt_step := 5.0 * rad;
      while opt_step > 1e-6 do 
      begin
        improved := False;
        
        // ???? ??????????? +v1
        test_v1 := best_v1 + opt_step; if test_v1 >= 2*Pi then test_v1 := test_v1 - 2*Pi;
        test_v2 := best_v2;
        COOR(q1, e1, i1, om1, w1, test_v1, MU, re, x1_m, y1_m, z1_m, vxe, vye, vze, v_dummy);
        COOR(q2, e2, i2, om2, w2, test_v2, MU, re, x2_m, y2_m, z2_m, vxe, vye, vze, v_dummy);
        CurrentDist := Sqrt(sqr(x1_m - x2_m) + sqr(y1_m - y2_m) + sqr(z1_m - z2_m));
        if CurrentDist < MOID then begin MOID := CurrentDist; best_v1 := test_v1; improved := True; end;

        // ???? ??????????? -v1
        test_v1 := best_v1 - opt_step; if test_v1 < 0 then test_v1 := test_v1 + 2*Pi;
        test_v2 := best_v2;
        COOR(q1, e1, i1, om1, w1, test_v1, MU, re, x1_m, y1_m, z1_m, vxe, vye, vze, v_dummy);
        COOR(q2, e2, i2, om2, w2, test_v2, MU, re, x2_m, y2_m, z2_m, vxe, vye, vze, v_dummy);
        CurrentDist := Sqrt(sqr(x1_m - x2_m) + sqr(y1_m - y2_m) + sqr(z1_m - z2_m));
        if CurrentDist < MOID then begin MOID := CurrentDist; best_v1 := test_v1; improved := True; end;

        // ???? ??????????? +v2
        test_v1 := best_v1;
        test_v2 := best_v2 + opt_step; if test_v2 >= 2*Pi then test_v2 := test_v2 - 2*Pi;
        COOR(q1, e1, i1, om1, w1, test_v1, MU, re, x1_m, y1_m, z1_m, vxe, vye, vze, v_dummy);
        COOR(q2, e2, i2, om2, w2, test_v2, MU, re, x2_m, y2_m, z2_m, vxe, vye, vze, v_dummy);
        CurrentDist := Sqrt(sqr(x1_m - x2_m) + sqr(y1_m - y2_m) + sqr(z1_m - z2_m));
        if CurrentDist < MOID then begin MOID := CurrentDist; best_v2 := test_v2; improved := True; end;

        // ???? ??????????? -v2
        test_v1 := best_v1;
        test_v2 := best_v2 - opt_step; if test_v2 < 0 then test_v2 := test_v2 + 2*Pi;
        COOR(q1, e1, i1, om1, w1, test_v1, MU, re, x1_m, y1_m, z1_m, vxe, vye, vze, v_dummy);
        COOR(q2, e2, i2, om2, w2, test_v2, MU, re, x2_m, y2_m, z2_m, vxe, vye, vze, v_dummy);
        CurrentDist := Sqrt(sqr(x1_m - x2_m) + sqr(y1_m - y2_m) + sqr(z1_m - z2_m));
        if CurrentDist < MOID then begin MOID := CurrentDist; best_v2 := test_v2; improved := True; end;

        if not improved then
          opt_step := opt_step * 0.5;
      end;
      // ======= END OF OPTIMIZED MOID CALCULATION BLOCK =======

      Writeln(JD1:6:2,' ',JD2:6:2,' ',t1:6:2,' ',t2:6:2);
      Writeln('ro2 = ',ro2:6:6);
      Writeln('D_SH = ',D_SH:6:6);
      Writeln('MOID = ',MOID:6:6);

      Writeln(outfile, JD1:7:5,' ',t1:4:5,' ', ro1,' ', ro2,' ', ro5,' ', D_SH,' ', D_D,' ', MOID);
      Break; 
    end;
  end;
end;

Close(infile1);
Close(infile2);
Close(outfile);
Writeln('Calculation finished successfully!');
end.
??????????? ??? ? ?????????????.????????????? ?????? ?????????? ????. ??? ??????????? ?????, ??????? ??????? ?? ???? ? ????????? ?????????? ? ?????? ???????? ????????????? ??? ?? ???? t=4238.31389, ????????? ?????????????.???? ? ??? ????????? ????????? ? ???????????, ????????:????? ?????? ??????????? ?? ??????????? (FPC, Delphi, Lazarus)?????? ?? ????? ?????????? ??? ?????? ??????? ?????????? (????????, Range Check Error)?? ?????? ???????????? ?????????? ?????????????????? ??????!

    ro1:=sqrt((1/L)*(p1+p2-2*sqrt(p1*p2)*CosI)+(sqr(e1)+sqr(e2)-2*e1*e2*CosP));
    ro2:=sqrt((1+sqr(e1))*p1+(1+sqr(e2))*p2-2*sqrt(p1*p2)*(CosI+e1*e2*CosP));
    ro5:=sqrt((1+sqr(e1))*p1+(1+sqr(e2))*p2-2*sqrt(p1*p2)*(e1*e2+Cos(i1-i2)));
    D_SH:=sqrt(sqr(q1-q2)+sqr(e1-e2)+4*sqr(Sin(I/2))+sqr(e1+e2)*sqr(Sin(PI_big)));
    D_D:=sqrt(sqr((e1-e2)/(e1+e2))+sqr((q1-q2)/(q1+q2))+sqr(I/Pi)+sqr((e1+e2)/2)*sqr(Pheta/Pi));

    // ======= START OF OPTIMIZED TWO-STAGE MOID CALCULATION BLOCK =======
    MOID := 1e10;
    best_v1 := 0;
    best_v2 := 0;

    // STAGE 1: Coarse Scan (??????????: ?????? ?????????? q1 ? q2)
    grid_step := 15.0 * rad;
    v1_step := 0;
    while v1_step < 2 * Pi do
    begin
      COOR(q1, e1, i1, om1, w1, v1_step, MU, re, x1_m, y1_m, z1_m, vxe, vye, vze, v_dummy);
      v2_step := 0;
      while v2_step < 2 * Pi do
      begin
        COOR(q2, e2, i2, om2, w2, v2_step, MU, re, x2_m, y2_m, z2_m, vxe, vye, vze, v_dummy);
        CurrentDist := Sqrt(sqr(x1_m - x2_m) + sqr(y1_m - y2_m) + sqr(z1_m - z2_m));

        if CurrentDist < MOID then
        begin
          MOID := CurrentDist;
          best_v1 := v1_step;
          best_v2 := v2_step;
        end;
        v2_step := v2_step + grid_step;
      end;
      v1_step := v1_step + grid_step;
    end;

    // STAGE 2: Gradient Descent with Angle Normalization
    opt_step := 5.0 * rad;
    while opt_step > 1e-6 do 
    begin
      improved := False;
      
      // ???? +v1
      test_v1 := best_v1 + opt_step; if test_v1 >= 2*Pi then test_v1 := test_v1 - 2*Pi;
      test_v2 := best_v2;
      COOR(q1, e1, i1, om1, w1, test_v1, MU, re, x1_m, y1_m, z1_m, vxe, vye, vze, v_dummy);
      COOR(q2, e2, i2, om2, w2, test_v2, MU, re, x2_m, y2_m, z2_m, vxe, vye, vze, v_dummy);
      CurrentDist := Sqrt(sqr(x1_m - x2_m) + sqr(y1_m - y2_m) + sqr(z1_m - z2_m));
      if CurrentDist < MOID then begin MOID := CurrentDist; best_v1 := test_v1; improved := True; end;

      // ???? -v1
      test_v1 := best_v1 - opt_step; if test_v1 < 0 then test_v1 := test_v1 + 2*Pi;
      test_v2 := best_v2;
      COOR(q1, e1, i1, om1, w1, test_v1, MU, re, x1_m, y1_m, z1_m, vxe, vye, vze, v_dummy);
      COOR(q2, e2, i2, om2, w2, test_v2, MU, re, x2_m, y2_m, z2_m, vxe, vye, vze, v_dummy);
      CurrentDist := Sqrt(sqr(x1_m - x2_m) + sqr(y1_m - y2_m) + sqr(z1_m - z2_m));
      if CurrentDist < MOID then begin MOID := CurrentDist; best_v1 := test_v1; improved := True; end;

      // ???? +v2
      test_v1 := best_v1; 
      test_v2 := best_v2 + opt_step; if test_v2 >= 2*Pi then test_v2 := test_v2 - 2*Pi;
      COOR(q1, e1, i1, om1, w1, test_v1, MU, re, x1_m, y1_m, z1_m, vxe, vye, vze, v_dummy);
      COOR(q2, e2, i2, om2, w2, test_v2, MU, re, x2_m, y2_m, z2_m, vxe, vye, vze, v_dummy);
      CurrentDist := Sqrt(sqr(x1_m - x2_m) + sqr(y1_m - y2_m) + sqr(z1_m - z2_m));
      if CurrentDist < MOID then begin MOID := CurrentDist; best_v2 := test_v2; improved := True; end;

      // ???? -v2
      test_v1 := best_v1; 
      test_v2 := best_v2 - opt_step; if test_v2 < 0 then test_v2 := test_v2 + 2*Pi;
      COOR(q1, e1, i1, om1, w1, test_v1, MU, re, x1_m, y1_m, z1_m, vxe, vye, vze, v_dummy);
      COOR(q2, e2, i2, om2, w2, test_v2, MU, re, x2_m, y2_m, z2_m, vxe, vye, vze, v_dummy);
      CurrentDist := Sqrt(sqr(x1_m - x2_m) + sqr(y1_m - y2_m) + sqr(z1_m - z2_m));
      if CurrentDist < MOID then begin MOID := CurrentDist; best_v2 := test_v2; improved := True; end;

      if not improved then
        opt_step := opt_step * 0.5;
    end;
    // ======= END OF OPTIMIZED MOID CALCULATION BLOCK =======

    Writeln(JD1:6:2,' ',JD2:6:2,' ',t1:6:2,' ',t2:6:2);
    Writeln('CosI = ',CosI:6:2,' ',p1:6:2,' ',p2:6:2);
    Writeln('CosP = ',CosP:6:6);
    Writeln('ro2 = ',ro2:6:6);
    Writeln('D_SH = ',D_SH:6:6);
    Writeln('D_D = ',D_D:6:6);
    Writeln('MOID = ',MOID:6:6);

    //Readln;


    //XE:=x;
    //YE:=y;
    //ZE:=z;
    //VXE:=xx;
    //VYE:=yy;
    //VZE:=zz;

    Writeln(outfile, JD1:7:5,' ',t1:4:5,' ', ro1,' ', ro2,' ', ro5,' ', D_SH,' ', D_D,' ', MOID);
    begin
        Break; // ????????? ????
    end;
   end;
 end;

 { ve:=v; {�???�? ?� ?�???}
{  if ve<0 then ve:=ve+2*pi;
  {�?�?�??? �???��}
{   COOR(A,E,I,OM,W,VE,mu,RE,XE,YE,ZE,VXE,VYE,VZE,VV);
   C:=SQRT((1/(rpart*dens*EXP(LN(RE)*2.25))-0.013*RC) *RC)*656;  {cm/sec}
{   C:=C*1E-5*(day_sec/au_km); {cm/sec => au/day}


{  KSI1:=RANDOM;KSI2:=RANDOM;
  F:=2*Pi*KSI1;CoT:=1-(1-Cos(alpha))*KSI2;
  SiT:=Sqrt(1-CoT*CoT);
  T:=ArcTg(SiT,CoT);


  LXc:=CoT;              {�?�? ? ?��??��� ?????�?}
 { LYc:=SiT*COS(F);
  LZc:=SiT*SIN(F);

  {?? ��???�?? ?��??��??�?�?? �?? ��? ?����?� ��?�?� �?????� !}
 { Convert_(-xe,-ye,-ze,vze*ye-vye*ze,vxe*ze-vze*xe,vye*xe-vxe*ye,
  LXc,LYc,LZc,LX,LY,LZ);  {� ??�?????????? ?� ?��??���}

{  CX:=C*LX;CY:=C*LY;CZ:=C*LZ;
  VXC:=VXE+CX;VYC:=VYE+CY;VZC:=VZE+CZ;
  VELOR3(XE,YE,ZE,VXC,VYC,VZC,AI,EI,II,OMI,WI,VI,MU);
  if ai<0 then continue;  {�?????�?� ??????��?}
{  if omi<0 then omi:=omi+2*pi;

  sE:=sin(vi)*re/(ai*sqrt(1-sqr(ei)));
  cE:=cos(vi)*re/ai+ei;
  EEi:=ArcTg(sE,cE);
  MI:=EEi-ei*sE;
  if Mi<0 then Mi:=Mi+2*PI;   }

  (*  ??�? ???�? �?�?�, ??�? ????�? �? ? ?�???

	   {?� ve ��?��?� ��?? ????�?�}
	  sE:=sin(ve)*re/(a*sqrt(1-sqr(e)));
	  cE:=cos(ve)*re/a+e;
	  EE:=ArcTg(sE,cE);
	  M:=EE-e*sE;

	  {�???�?�, ??� �?�?��? �????� �?��????? ? ���?�?? ??�?���?�?? ??????�??}
	  if M<Pi then
	    JDi:=JD+M/n_   {M in rad}
	  else JDi:=JD-(2*Pi-M)/n_;

	  if omi<0 then omi:=omi+2*Pi;

  *)

 { Writeln(outfile, ai:15:10,ei:15:10,ii/rad:15:9,omi/rad:15:9,wi/rad:15:9,Mi/rad:15:9,
           ' ',Jdi:14:5,'  ',ve/rad:5:1,'  ',j:5);
(*     '  ',T/k:5:1,'  ',F/k:5:1,'  ',c/(kv*1e-3):5:1,' ',ve/k:5:1,' ',*)
                                      (* m/s  *)

  {$IFDEF DEBUG}
 (* Writeln(outdebug,'ve = ', ve);
  Writeln(outdebug,'T  = ', T/k);
  coT1:=(-xe*Lx-ye*Ly-ze*Lz)/re;
  Writeln(outdebug,'coT1 = ', coT1);
  Writeln(outdebug);
  Writeln(outdebug,'xe = ', xe);
  Writeln(outdebug,'ye = ', ye);
  Writeln(outdebug,'ze = ', ze);
  Writeln(outdebug,'Vxe = ', Vxe);
  Writeln(outdebug,'Vye = ', Vye);
  Writeln(outdebug,'Vze = ', Vze);
  Writeln(outdebug,'cx = ', cx);
  Writeln(outdebug,'cy = ', cy);
  Writeln(outdebug,'cz = ', cz);
  Close(outdebug);
  Halt;
  {$ENDIF}
          *)
end; {j}

Close(infile1);
Close(infile2);
Close(outfile);
end.
