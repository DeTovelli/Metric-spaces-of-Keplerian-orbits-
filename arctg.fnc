function ArcTg(x,y: Extended):Extended;
{Arctang круговой угла alpha
 в диапазоне [0,2*Pi];
 x = c*sin(alpha)
 y = c*cos(alpha),
 с - положительная неизвестная постоянная }

 var a:Extended;
 begin
  if abs(y)<1e-18 then ArcTg:=sign(x)*0.5*Pi
  else
    begin
      a:=ArcTan(x/y);
      if y<0 then a:=a+Pi;
      if a<0 then a:=a+2*Pi;
      ArcTg:=a;
    end;
 end; {ArcTg}
