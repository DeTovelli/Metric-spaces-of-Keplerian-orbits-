function power(a:Extended; b:Extended): Extended;
{ Возведение числа А в степень В }

begin
  if a<>0 then  power:=exp(b*ln(a))
  else power:=0;
end; (* power *)