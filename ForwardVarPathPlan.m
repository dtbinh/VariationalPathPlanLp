function [ z ] = ForwardVarPathPlan(t,x)
global sigma tau m J obs obsOrder obsScale;
z = zeros(3*4+1,1); % 12 dimension [th, x, y ,lambda]
TH = x(1);
THd = x(2);
THdd = x(3);
THddd = x(4);

XX = x(5);
XXd = x(6);
XXdd = x(7);
XXddd = x(8);

YY = x(9);
YYd = x(10);
YYdd = x(11);
YYddd = x(12);

MU = x(13);

z(1) = THd;
z(2) = THdd;
z(3) = THddd; %-- th

z(5) = XXd;
z(6) = XXdd;
z(7) = XXddd; % - x

z(9) = YYd;
z(10) = YYdd;
z(11) = YYddd; % - y
Vx=tau*obsOrder/2*(XX-obs(1))^(obsOrder-1)/obsScale(1)/(((XX-obs(1))/obsScale(1))^obsOrder+((YY-obs(2))/obsScale(2))^obsOrder-1^obsOrder)^2;
Vy=tau*obsOrder/2*(YY-obs(2))^(obsOrder-1)/obsScale(2)/(((XX-obs(1))/obsScale(1))^obsOrder+((YY-obs(2))/obsScale(2))^obsOrder-1^obsOrder)^2;
MUd = m*(-3*(XXddd*cos(TH)+YYddd*sin(TH))+(sigma-2)*(XXd*cos(TH)+YYd*sin(TH))-Vx*sin(TH)+Vy*cos(TH));
z(4) = sigma*THdd-MU/J*(XXd*cos(TH)+YYd*sin(TH));
z(8) = sigma*XXdd+Vx+1/m*(MUd*sin(TH)+MU*THd*cos(TH));
z(12) = sigma*YYdd+Vy+1/m*(-MUd*cos(TH)+MU*THd*sin(TH));
z(13) = MUd;

end

