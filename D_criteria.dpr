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
  // Mathematical range error protection variables (NaN)
  ArgArcSin, ArgPheta: Extended;

infile, infile1, infile2, outfile, outdebug: text;

{$I kepl1.fnc}
{$I coor.prc}
{$I velor3.prc}
{$I Convert_.prc}
{$I power.fnc}
{$I r4.fnc}
{$I arctg.fnc}

// Auxiliary function for protection against exceeding the range [-1, 1]
function Clamp(Value: Extended): Extended;
begin
  if Value > 1.0 then Result := 1.0
  else if Value < -1.0 then Result := -1.0
  else Result := Value;
end;

begin {main}
  // Check the physical presence of files on the disk before launching
  if not FileExists('elem_point1_2003EH1.dat') then
  begin
    Writeln('Error: Input file "elem_point1_2003EH1.dat" not found!');
    Exit;
  end;

  if not FileExists('elem_point2_2009GS18.dat') then
  begin
    Writeln('Error: Input file "elem_point2_2009GS18.dat" not found!');
    Exit;
  end;

Assign(outfile, 'file.TXT');
  Rewrite(outfile);
  Writeln(outfile, 'JD t ro1 ro2 ro5 D_SH D_D');

  Assign(infile1, 'elem_point1_2003EH1.dat');
  Reset(infile1);
  
  // Skipping the first three header lines of the first file
  Readln(infile1);
  Readln(infile1);
  Readln(infile1, JD1, t1, a1, e1, i1, om1, w1, M1, q1);
  
  i1 := i1 * rad; om1 := om1 * rad; w1 := w1 * rad; M1 := M1 * rad;
  Writeln('First object initial entry loaded.');

  Assign(infile2, 'elem_point2_2009GS18.dat');
  
  L := 1;

while not Eof(infile1) do
begin
Readln(infile1, JD1,t1,a1,e1,i1,om1,w1,M1,q1);
// Converting angles from degrees to radians
i1:=i1*rad; om1:=om1*rad; w1:=w1*rad; M1:=M1*rad; {��������� ������� �������� ������� � �������}
p1:=a1*(1-sqr(e1));

Writeln('TJD1 = ',JD1:6:2,' t1 = ',t1:6:2);
//Readln;

// Reset infile2 to the beginning before searching
Reset(infile2);
Readln(infile2);
Readln(infile2);
// Read infile2 to the end or until we find N records
while not Eof(infile2) do
//for j:=1 to N do
 begin

  Readln(infile2, JD2,t2,a2,e2,i2,om2,w2,M2,q2);
  //if JD1=JD2 then
  if Abs(JD1 - JD2) < 1e-6 then// Comparing Float numbers for equality using delta
   begin
   // Converting angles from degrees to radians
    i2:=i2*rad; om2:=om2*rad; w2:=w2*rad; M2:=M2*rad; {��������� ������� �������� ������� � �������}
    p2:=a2*(1-sqr(e2));

    // 1. Protection for the orbital inclination function
    CosI:=(Cos(i1)*Cos(i2))+(Sin(i1)*Sin(i2)*Cos(om1-om2));
    CosI := Clamp(CosI); // We guarantee that CosI lies in [-1; 1]
    CosP:=Sin(i1)*Sin(i2)*Sin(w1)*Sin(w2)+(Cos(w1)*Cos(w2)+Cos(i1)*Cos(i2)*Sin(w1)*Sin(w2))*Cos(om1-om2)+(Cos(i2)*Cos(w1)*Sin(w2)-Cos(i1)*Sin(w1)*Cos(w1))*Sin(om1-om2);

    I:=ArcCos(CosI);

    // 2. Protection for calculating the ArcSin argument (for PI_big)
    // Check that the secant does not cause division by zero if Cos(I/2) is near zero
        if Abs(Cos(I / 2)) > 1e-9 then
        begin
          ArgArcSin := Cos(i2 + i1) * Sin((om2 - om1) / 2) * (1.0 / Cos(I / 2));
          ArgArcSin := Clamp(ArgArcSin);  // Protection against going beyond [-1; 1]

          if Abs(om1 - om2) > Pi then
            PI_big := om2 - om1 - 2 * ArcSin(ArgArcSin)
          else
            PI_big := om2 - om1 + 2 * ArcSin(ArgArcSin);
        end
        else
        begin
          PI_big := 0;  // Prevent division by 0 in collinear orbits
        end;

    // 3. Protection for Pheta
    ArgPheta := Sin(i1)*Sin(i2)*Sin(w1)+(Cos(w1)*Cos(w2)+Cos(i1)*Cos(i2)*Sin(w1)*Sin(w2))*Cos(om1-om2)+(Cos(i1)*Cos(w1)*Cos(w2)-Cos(i1)*Sin(w1)*Sin(w2))*Sin(om1-om2);
    ArgPheta := Clamp(ArgPheta); // Protection against going beyond [-1; 1]
    Pheta := ArcCos(ArgPheta);

    // Calculation of the Kholshevnikov, Southworth-Hawkins, and Drummond metrics
    ro1:=sqrt((1/L)*(p1+p2-2*sqrt(p1*p2)*CosI)+(sqr(e1)+sqr(e2)-2*e1*e2*CosP));
    ro2:=sqrt((1+sqr(e1))*p1+(1+sqr(e2))*p2-2*sqrt(p1*p2)*(CosI+e1*e2*CosP));
    ro5:=sqrt((1+sqr(e1))*p1+(1+sqr(e2))*p2-2*sqrt(p1*p2)*(e1*e2+Cos(i1-i2)));
    D_SH:=sqrt(sqr(q1-q2)+sqr(e1-e2)+4*sqr(Sin(I/2))+sqr(e1+e2)*sqr(Sin(PI_big)));
    D_D:=sqrt(sqr((e1-e2)/(e1+e2))+sqr((q1-q2)/q1+q2)+sqr(I/Pi)+sqr((e1+e2)/2)*sqr(Pheta/Pi));

    Writeln(JD1:6:2,' ',JD2:6:2,' ',t1:6:2,' ',t2:6:2);
    Writeln('CosI = ',CosI:6:2,' ',p1:6:2,' ',p2:6:2);
    Writeln('CosP = ',CosP:6:6);
    Writeln('ro2 = ',ro2:6:6);
    Writeln('D_SH = ',D_SH:6:6);
    Writeln('D_D = ',D_D:6:6);

    // // Write the results to the output file
    Writeln(outfile, JD1:7:5,' ',t1:4:5,' ', ro1,' ', ro2,' ', ro5,' ', D_SH,' ', D_D,' '{, VXE,' ', VYE,' ', VZE });

    Break;  // Found a match for the current JD1 - exit infile2 and take the next line infile1

    end;
   end;
 end;

Close(infile1);
Close(infile2);
Close(outfile);
end.
