procedure Convert_(x1,x2,x3,z1,z2,z3,xc,yc,zc:EXTENDED;var xg,yg,zg:EXTENDED);
{Переводит из гелиоцентрической оси которой определяются векторами
 (x1,x2,x3) и (z1,z2,z3) в гелиоцентрическую эклиптическую }

Var    t11,t12,t13,t21,t22,t23,t31,t32,t33 : EXTENDED;
       lz1,ly1,lx1: EXTENDED;
begin
    lx1:=sqrt(x1*x1+x2*x2+x3*x3);
    lz1:=sqrt(z1*z1+z2*z2+z3*z3);
    ly1:=lz1*lx1;
    t11:=x1/lx1;
    t21:=x2/lx1;
    t31:=x3/lx1;
    t13:=z1/lz1;
    t23:=z2/lz1;
    t33:=z3/lz1;
    t12:=(x3*z2-x2*z3)/ly1;
    t22:=(x1*z3-x3*z1)/ly1;
    t32:=(x2*z1-x1*z2)/ly1;

    xg:=t11*xc+t12*yc+t13*zc;
    yg:=t21*xc+t22*yc+t23*zc;
    zg:=t31*xc+t32*yc+t33*zc;

end; {Convert}