function kepl1(M,e: Extended): Extended;

{M - средн¤¤ аномали¤, е - эксцентриситет }

{ –ешает уравнение  еплера  M=E - e sin E
стартер типа Ћежандра, итератор √алле¤.
»спользуетс¤ как стандартна¤ в Royal Aircraft
Establishment (England) под именем EKEPL1.
¬з¤то из статьи Odell A.W., Gooding R.H.
Cel.Mech. 1986,38,4,307-336 }

var 
  c,eta,f,fd,fdd,psi,s,xi: Extended;
  iter_count: Integer; // »—ѕ–ј¬Ћ≈Ќќ: добавили счетчик итераций

begin
c:=e*cos(M);
s:=e*sin(M);
psi:=s/sqrt(1-c-c+e*e);

iter_count := 0; // »нициализируем счетчик

repeat
  xi:=cos(psi);
  eta:=sin(psi);
  fd:=(1-c*xi)+s*eta;
  fdd:=c*eta+s*xi;
  f:=psi-fdd;
  
  // «јў»“ј: провер¤ем знаменатель, чтобы избежать делени¤ на ноль при сбо¤х метода √алле¤
  if Abs(fd*fd - 0.5*f*fdd) < 1e-15 then 
    break;
    
  psi:=psi-f*fd/(fd*fd-0.5*f*fdd);
  
  Inc(iter_count);
  // «јў»“ј: если за 100 шагов метод не сошелс¤ (дл¤ e -> 1), 
  // принудительно выходим, предотвраща¤ вечное зависание программы
  if iter_count > 100 then 
    break;

until f*f <= 1e-12;

kepl1:=M+psi;

end; {kepl1}
