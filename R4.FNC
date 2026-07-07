function r4(a,e:Extended): Extended;
{Генерирует v по закону r^(-4);
  v in [-Pi,+Pi];
  v для r из двух возможных выбирается случайно}

var alpha,Q3,q3_,p,r:Extended;
    beta: Integer;
begin
Q3:=power(a*(1+e),-3);
q3_:=power(a*(1-e),-3);
p:=a*(1-sqr(e));
alpha:=Random;
beta:=Random(2);
r:=power(1/((Q3-q3_)*alpha+q3_),1.0/3.0);
if beta=0 then r4:=ArcCos(abs(p-r)/(e*r))
else r4:=-ArcCos(abs(p-r)/(e*r));
end; {r4}