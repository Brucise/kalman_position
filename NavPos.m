function  [x,y,z,delt_sat]=NavPos(t,m,mts)

global  NAV  c OMEGAe GM
%t为估计的观测卫星发送的时间
%m为最佳的观测卫星在NAV中的位置

% t=T(i);%1.727999301119927e+05;
% m=kk(i);%13;
% mts=MTS(i)%0.0699;
%计算角平速度
n0=sqrt(GM/NAV(m).A^3);

%计算时间差
tk=t-NAV(m).toe;
if  tk >302400
    tk=tk-604800;
end
if  tk <-302400
    tk=tk+604800;
end;

%计算平近点角
n=n0+NAV(m).dn;   
Mk=NAV(m).M0+n*tk;  

%迭代计算偏近点角
Ek=0;
while 1
    Ek_=(Mk+NAV(m).e*sin(Ek));
    if (abs(Ek_-Ek)<10^(-12))
        Ek=Ek_;
        break;
    else
        Ek=Ek_;
    end;
end;  

%计算真近点角
vk=2*atan(sqrt((1.0+NAV(m).e)/(1-NAV(m).e))*tan(Ek/2));
%vk=atan(sin(Ek)*sqrt(1-NAV(m).e^2)/(cos(Ek)-NAV(m).e));

%计算信号发送时刻交升点焦距
phik=vk+NAV(m).w;    

%2.7 计算周期改正项
duk=NAV(m).Cus*sin(2*phik)+NAV(m).Cuc*cos(2*phik); %17行为Cus，15行为Cuc；
drk=NAV(m).Crs*sin(2*phik)+NAV(m).Crc*cos(2*phik); %12行为Crs，24行为Crc；
dik=NAV(m).Cis*sin(2*phik)+NAV(m).Cic*cos(2*phik); %22行为Cis，20行为Cic；

%计算改正后的纬度参数、向经，倾角
uk=phik+duk;
rk=(1-NAV(m).e*cos(Ek))*NAV(m).A+drk;
ik=NAV(m).i0+dik+NAV(m).idot*tk;    

%计算卫星在轨道平面内的坐标
yk_=rk*sin(uk);
xk_=rk*cos(uk);

%改正升交点的经度
omegak=NAV(m).OMEGA0+(NAV(m).OMEGAdot-OMEGAe)*tk-OMEGAe*NAV(m).toe;  

%计算卫星在地固系中的坐标
x_dot=xk_*cos(omegak)-yk_*cos(ik)*sin(omegak);
y_dot=xk_*sin(omegak)+yk_*cos(ik)*cos(omegak);
z_dot=yk_*sin(ik); %计算卫星在WGS84地固系坐标

%进行地球自转改正
aph=OMEGAe*mts;
temp_xyz=[cos(aph),sin(aph),0;-sin(aph),cos(aph),0;0,0,1]*[x_dot,y_dot,z_dot]';
x=temp_xyz(1);
y=temp_xyz(2);
z=temp_xyz(3);

% %3计算相对论效应和卫星钟差
rela=-2*sqrt(GM)*NAV(m).e*sqrt(NAV(m).A)*sin(Ek)/c^2;
clk_sat=NAV(m).clkt(1)+NAV(m).clkt(2)*tk+NAV(m).clkt(3)*tk^2;
delt_sat=rela+clk_sat;
end