PROCEDURE  COOR(Q,E,I,OM,W,V,MU: Extended;
                Var R,X,Y,Z,VX,VY,VZ,VV: EXTENDED);
VAR    CI,COM,CU,SI,SOM,SU,P,P1,R1,VR,VN,VRR,U: EXTENDED;
BEGIN
U:=V+W;

// БЕЗОПАСНО: Вычисляем фокальный параметр через q, а не через a
P:=Q*(1+E); 

// Защита от деления на ноль в корне, если q или e повреждены
if P <= 1e-14 then P := 1e-14; 
P1:=SQRT(MU/P);

R1:=1+E*COS(V);
// Защита от деления на ноль в апоцентре при e -> 1
if Abs(R1) < 1e-14 then 
  R1 := 1e-14; 

R:=P/R1;
VR:=P1*E*SIN(V);
VN:=P1*R1;
VRR:=VR/R;

COM:=COS(OM);SOM:=SIN(OM);
CU:=COS(U);SU:=SIN(U);
CI:=COS(I);SI:=SIN(I);

X:=R*(CU*COM-SU*SOM*CI);
Y:=R*(CU*SOM+SU*COM*CI);
Z:=R*SU*SI;

VX:=X*VRR+(-SU*COM-CU*SOM*CI)*VN;
VY:=Y*VRR+(-SU*SOM+CU*COM*CI)*VN;
VZ:=Z*VRR+CU*SI*VN;
VV:=SQRT(VX*VX+VY*VY+VZ*VZ);
END; (*COOR*)
