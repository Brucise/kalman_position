function  [x,y,z,delt_sat]=NavPos(t,m,mts)

global  NAV  c OMEGAe GM
%tΪ���ƵĹ۲����Ƿ��͵�ʱ��
%mΪ��ѵĹ۲�������NAV�е�λ��

% t=T(i);%1.727999301119927e+05;
% m=kk(i);%13;
% mts=MTS(i)%0.0699;
%�����ƽ�ٶ�
n0=sqrt(GM/NAV(m).A^3);

%����ʱ���
tk=t-NAV(m).toe;
if  tk >302400
    tk=tk-604800;
end
if  tk <-302400
    tk=tk+604800;
end;

%����ƽ�����
n=n0+NAV(m).dn;   
Mk=NAV(m).M0+n*tk;  

%��������ƫ�����
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

%����������
vk=2*atan(sqrt((1.0+NAV(m).e)/(1-NAV(m).e))*tan(Ek/2));
%vk=atan(sin(Ek)*sqrt(1-NAV(m).e^2)/(cos(Ek)-NAV(m).e));

%�����źŷ���ʱ�̽����㽹��
phik=vk+NAV(m).w;    

%2.7 �������ڸ�����
duk=NAV(m).Cus*sin(2*phik)+NAV(m).Cuc*cos(2*phik); %17��ΪCus��15��ΪCuc��
drk=NAV(m).Crs*sin(2*phik)+NAV(m).Crc*cos(2*phik); %12��ΪCrs��24��ΪCrc��
dik=NAV(m).Cis*sin(2*phik)+NAV(m).Cic*cos(2*phik); %22��ΪCis��20��ΪCic��

%����������γ�Ȳ������򾭣����
uk=phik+duk;
rk=(1-NAV(m).e*cos(Ek))*NAV(m).A+drk;
ik=NAV(m).i0+dik+NAV(m).idot*tk;    

%���������ڹ��ƽ���ڵ�����
yk_=rk*sin(uk);
xk_=rk*cos(uk);

%����������ľ���
omegak=NAV(m).OMEGA0+(NAV(m).OMEGAdot-OMEGAe)*tk-OMEGAe*NAV(m).toe;  

%���������ڵع�ϵ�е�����
x_dot=xk_*cos(omegak)-yk_*cos(ik)*sin(omegak);
y_dot=xk_*sin(omegak)+yk_*cos(ik)*cos(omegak);
z_dot=yk_*sin(ik); %����������WGS84�ع�ϵ����

%���е�����ת����
aph=OMEGAe*mts;
temp_xyz=[cos(aph),sin(aph),0;-sin(aph),cos(aph),0;0,0,1]*[x_dot,y_dot,z_dot]';
x=temp_xyz(1);
y=temp_xyz(2);
z=temp_xyz(3);

% %3���������ЧӦ�������Ӳ�
rela=-2*sqrt(GM)*NAV(m).e*sqrt(NAV(m).A)*sin(Ek)/c^2;
clk_sat=NAV(m).clkt(1)+NAV(m).clkt(2)*tk+NAV(m).clkt(3)*tk^2;
delt_sat=rela+clk_sat;
end