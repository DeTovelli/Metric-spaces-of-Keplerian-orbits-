function kepl1(M,e: Extended): Extended;

{M - средняя аномалия, е - эксцентриситет }

{ Решает уравнение Кеплера  M=E - e sin E
стартер типа Лежандра, итератор Галлея.
Используется как стандартная в Royal Aircraft
Establishment (England) под именем EKEPL1.
Взято из статьи Odell A.W., Gooding R.H.
Cel.Mech. 1986,38,4,307-336 }

var c,eta,f,fd,fdd,psi,s,xi: Extended;

begin
c:=e*cos(M);
s:=e*sin(M);
psi:=s/sqrt(1-c-c+e*e);

repeat
  xi:=cos(psi);
  eta:=sin(psi);
  fd:=(1-c*xi)+s*eta;
  fdd:=c*eta+s*xi;
  f:=psi-fdd;
  psi:=psi-f*fd/(fd*fd-0.5*f*fdd);
until f*f <= 1e-12;

kepl1:=M+psi;

end; {kepl1}