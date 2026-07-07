Program D_criteria;

{$APPTYPE CONSOLE}

{for verhin, metod delenia+}

uses
  Math,
  SysUtils;

Label MET1,AGAIN;

Const N=85641074;  { оличество моделируемых частиц}
    {IAU 1994}
    KV=1/29.7846918325927450;  {ƒл€ перевода скорости в км/с при mu=1}
    rad=Pi/180.0;
    au_km=1.49597870691000015E+8; {1 au in km; DE406}
    day_sec = 86400;  {сутки в секенду}
    KMAU=1/au_km;     {км в ј≈}
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
CosI,CosP,ro1,ro2,ro5:Extended;
j,jj:integer;

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
{—читываем шесть элементов, угловые в градусах и JD}
{Read(infile, a,e,i,om,w,M,JD);
{переводим угловые элементы опорной в радианы}
{i:=i*rad; om:=om*rad; w:=w*rad; M:=M*rad;
JDi:=JD;

Close(infile);}

{Assign(outdebug,'example.out');
Rewrite(outdebug);  }



Assign(outfile,'file.TXT');
Rewrite(outfile);

Writeln(outfile, 'JD t ro1 ro2 ro5');
RandSeed:=1;
{Randomize;}

{RC:=1.9;     {km, радиус €дракометы}
{dens:=3.0;  {g cm-3, плотность метеороида}
{mass:=1.0000; {g, масса метеороида}
{rpart:=power(0.75*mass/(Pi*dens),1.0/3.0);  {cm, радиус метеороида}
{beta:=5.8E-5/(rpart*dens); {Fr/Fgr -- безразмерный коэффициент}
{mu_red:=MU*(1-beta); {–едукци€ массы —олнца за —ветовое давление}


{kn:=sqrt(GM/(au_km*au_km*au_km))*day_sec;
n_:=kn/(a*sqrt(a));    {рад/сут}


{по средней аномалии находим истинную}
{EE:=kepl1(M,e);
r:=a*(1-e*cos(EE));
sv:=a*sqrt(1-e*e)*sin(EE)/r;
cv:=a*(cos(EE)-e)/r;
v:=ArcTg(sv,cv);

{¬ыводим в файл начальную орбиту}
{j:=0;
Writeln(outfile, a:15:10,e:15:10,i/rad:15:9,om/rad:15:9,w/rad:15:9,M/rad:15:9,
           ' ',Jdi:14:5,'  ',v/rad:5:1,'  ',j:5);


{ќсь ’ направлена к —олнцу}
{alpha:=180*rad;  {ѕолураствор конуса, 180 дл€ }
//TJD0:= 2459586.50;
Assign(infile1,'elem_point1_2003EH1.dat');
Reset(infile1);
{—читываем JD, шесть элементов, угловые в градусах }
Readln(infile1);
Readln(infile1);
Readln(infile1, JD1,t1,a1,e1,i1,om1,w1,M1,q1);
i1:=i1*rad; om1:=om1*rad; w1:=w1*rad; M1:=M1*rad; {переводим угловые элементы опорной в радианы}
Writeln(JD1:6:2,t1:6:2,a1:6:2,e1:6:2,i1:6:2,om1:6:2,w1:6:2,M1:6:2,q1:6:2);
//Readln;

Assign(infile2,'elem_point2_2009GS18.dat');
Reset(infile2);
{—читываем шесть элементов, угловые в градусах и JD}
Readln(infile2);
Readln(infile2);
Readln(infile2, JD2,t2,a2,e2,i2,om2,w2,M2,q2);
i2:=i2*rad; om2:=om2*rad; w2:=w2*rad; M2:=M2*rad; {переводим угловые элементы опорной в радианы}

Writeln(JD2:6:2,t2:6:2,a2,e2,i2,om2,w2,M2,q2);
//Readln;
L:=1;

while not Eof(infile1) do
begin
Readln(infile1, JD1,t1,a1,e1,i1,om1,w1,M1,q1);
i1:=i1*rad; om1:=om1*rad; w1:=w1*rad; M1:=M1*rad; {переводим угловые элементы опорной в радианы}
p1:=a1*(1-sqr(e2));

Writeln('TJD1 = ',JD1:6:2,' t1 = ',t1:6:2,t1:6:2);
//Readln;

for j:=1 to N do
 begin

  Readln(infile2, JD2,t2,a2,e2,i2,om2,w2,M2,q2);
  if JD1=JD2 then
   begin
    i2:=i2*rad; om2:=om2*rad; w2:=w2*rad; M2:=M2*rad; {переводим угловые элементы опорной в радианы}
    p2:=a2*(1-sqr(e2));
    CosI:=(Cos(i1)*Cos(i2))+(Sin(i1)*Sin(i2)*Cos(om1-om2));
    CosP:=Sin(i1)*Sin(i2)*Sin(w1)*Sin(w2)+(Cos(w1)*Cos(w2)+Cos(i1)*Cos(i2)*Sin(w1)*Sin(w2))*Cos(om1-om2)+(Cos(i2)*Cos(w1)*Sin(w2)-Cos(i1)*Sin(w1)*Cos(w1))*Sin(om1-om2);
    ro1:=sqrt((1/L)*(p1+p2-2*sqrt(p1*p2)*CosI)+(sqr(e1)+sqr(e2)-2*e1*e2*CosP));
    ro2:=sqrt((1+sqr(e1))*p1+(1+sqr(e2))*p2-2*sqrt(p1*p2)*(CosI+e1*e2*CosP));
    ro5:=sqrt((1+sqr(e1))*p1+(1+sqr(e2))*p2-2*sqrt(p1*p2)*(e1*e2+Cos(i1-i2)));
    Writeln(JD1:6:2,' ',JD2:6:2,' ',t1:6:2,' ',t2:6:2);
    Writeln('CosI = ',CosI:6:2,' ',p1:6:2,' ',p2:6:2);
    Writeln('CosP = ',CosP:6:6);
    Writeln('ro2 = ',ro2:6:6);
    //Readln;


    //XE:=x;
    //YE:=y;
    //ZE:=z;
    //VXE:=xx;
    //VYE:=yy;
    //VZE:=zz;

    Writeln(outfile, JD1:7:5,' ',t1:4:5,' ', ro1,' ', ro2,' ', ro5,' '{, VXE,' ', VYE,' ', VZE });
    begin
        Break; // ????????? ????
    end;
   end;
 end;

 { ve:=v; {В???Ѓ? ?І ?Ѓ???}
{  if ve<0 then ve:=ve+2*pi;
  {С?Ѓ?Ѓ??? У???Ђ†}
{   COOR(A,E,I,OM,W,VE,mu,RE,XE,YE,ZE,VXE,VYE,VZE,VV);
   C:=SQRT((1/(rpart*dens*EXP(LN(RE)*2.25))-0.013*RC) *RC)*656;  {cm/sec}
{   C:=C*1E-5*(day_sec/au_km); {cm/sec => au/day}


{  KSI1:=RANDOM;KSI2:=RANDOM;
  F:=2*Pi*KSI1;CoT:=1-(1-Cos(alpha))*KSI2;
  SiT:=Sqrt(1-CoT*CoT);
  T:=ArcTg(SiT,CoT);


  LXc:=CoT;              {У?Ђ? ? ?Ѓђ??≠Ѓ© ?????ђ?}
 { LYc:=SiT*COS(F);
  LZc:=SiT*SIN(F);

  {?? І†???†?? ?Ѓ≠??ЃЂ??Ѓ?†?? Ѓ?? §Ђ? ?†¶§Ѓ?Ѓ ≠Ѓ?Ѓ?Ѓ Ѓ?????† !}
 { Convert_(-xe,-ye,-ze,vze*ye-vye*ze,vxe*ze-vze*xe,vye*xe-vxe*ye,
  LXc,LYc,LZc,LX,LY,LZ);  {В ??Ђ?????????? ?І ?Ѓђ??≠Ѓ©}

{  CX:=C*LX;CY:=C*LY;CZ:=C*LZ;
  VXC:=VXE+CX;VYC:=VYE+CY;VZC:=VZE+CZ;
  VELOR3(XE,YE,ZE,VXC,VYC,VZC,AI,EI,II,OMI,WI,VI,MU);
  if ai<0 then continue;  {Ѓ?????†?ђ ??????ЃЂ?}
{  if omi<0 then omi:=omi+2*pi;

  sE:=sin(vi)*re/(ai*sqrt(1-sqr(ei)));
  cE:=cos(vi)*re/ai+ei;
  EEi:=ArcTg(sE,cE);
  MI:=EEi-ei*sE;
  if Mi<0 then Mi:=Mi+2*PI;   }

  (*  ??Ѓ? ???Ѓ? ≠?¶?≠, ??Ђ? ????Ѓ? ≠? ? ?Ѓ???

	   {?Ѓ ve ≠†?Ѓ§?ђ §†?? ????Ѓ?†}
	  sE:=sin(ve)*re/(a*sqrt(1-sqr(e)));
	  cE:=cos(ve)*re/a+e;
	  EE:=ArcTg(sE,cE);
	  M:=EE-e*sE;

	  {С???†?ђ, ??Ѓ Ѓ?Ѓ?≠†? Ѓ????† Ѓ?≠Ѓ????? ? ђЃђ?≠?? ??Ѓ?Ѓ¶§?≠?? ??????Ђ??}
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
