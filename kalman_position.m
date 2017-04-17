clc;
clear all
global  NAV c GM OMEGAe dtr %定义全局变量，用来计算卫星位置
c=2.997925*10^8;  %光速
GM=3.986004418e+14;        %给出万有引力常数
OMEGAe=7.292115*10^(-5);   %给出地球自转角速度（弧度制）
dtr = pi/180;              %角度转弧度
%==========读导航文件，版本是302============================================
if exist('NAV.mat')~=2
    [NAV,Header]= ReadNavfile_302;
else
    load('NAV.mat');
end
%==========读观测值文件，版本是302=====================================================
if exist('OBS.mat')~=2
    [OBS,XYZ]= Read_17O_files();
else
    load('OBS.mat')
end
%===========开始定位=======================================================
POS=zeros(length(OBS),8);                  %给初始矩阵装结果
dt=0.1;                                    %历元间间隔
Phic=[1 dt;0 1];                           %钟差转移矩阵
Phi=[eye(3) dt*eye(3) zeros(3,2);...       
    zeros(3,3) eye(3) zeros(3,2);...
    zeros(2,6) Phic];                                     %状态转移矩阵
Sv=[10^-8 10^-8 10^-8];St=0.05;Sf=0.05;     %速度噪声,钟差噪声和钟漂噪声功率谱密度
Sv=diag(Sv);
Qxyz=[Sv*dt^3 Sv*dt^2/2;...
      Sv*dt^2/2 Sv*dt];                        %坐标分量对应的协方差矩阵
Qc=[St*dt+Sf*dt^3/3 Sf*dt^2/2;Sf*dt^2/2 Sf*dt];%接收机时钟的过程噪声协方差矩阵
Qk=[Qxyz zeros(6,2);...                          %噪声方差
    zeros(2,6) Qc];
for epoch=1:length(OBS)
%================剔除当前历元下，.N中不存在但是在.O文件中存在的卫星===========
    sats_O=OBS(epoch).prn;
    for i=1:length(NAV)
    sats_N(i)=NAV(i).prn;
    end
    Sats = intersect(sats_O,sats_N);
    obs=OBS(epoch).C1;
    lia = ismember(sats_O,Sats);
    sats_O(lia==0) = [];
    sats=sats_O;
    obs(lia==0) = [];  
%===============剔除完毕====================================================
    sow=OBS(epoch).sow;                 %second of week
    ion=Header.ion_gps;                 %计算电离层延迟参数
    tepo=OBS(epoch).tepo;               %jd和sod，也用来计算k8电离层延迟
    tt(1)=OBS(epoch).sow;               %用来计算对流层延迟
    tt(2)=OBS(epoch).www;
    m = size(obs,1);                    % 当前历元下卫星数目
    A=[];                               %最小二乘里的G矩阵，G*[delta_x]=b
    omc=[];                             %最小二乘里的b矩阵
    w=[];                               %权重阵
    R=[];                               %观测值方差阵
    if epoch==1                %最小二乘法确认初始状态方差阵
             [Q,coor]=recpo_ls_klman(obs,sats,sow,ion,tepo);
             sigma0=900/(Q(1,1)+Q(2,2)+Q(3,3));
             P=[sigma0*Q(1,1) sigma0*Q(2,2) sigma0*Q(3,3)...
                2*sigma0*Q(1,1)/dt^2 2*sigma0*Q(2,2)/dt^2 2*sigma0*Q(3,3)/dt^2 ...
                sigma0*Q(4,4) 2*sigma0*Q(4,4)/dt^2];%P矩阵定义出自论文--卡尔曼滤波算法在GPS非差相位精密单点定位中的应用研究_祁芳.caj
             P=diag(P); 
             X=[coor(1) coor(2) coor(3)  0 0 0 coor(4) 0]';%给初始坐标，速度和钟漂都给零。
     end
    for i=1:m
%================在NAV中寻找最适宜的导航卫星=================================
        j=1;comp_mat=[];
        for k=1:size(NAV,2)
                if(NAV(k).prn==sats(i))%在NAV中找到该卫星。
                    dtime=abs(NAV(k).toe-OBS(epoch).sow);%把该卫星在NAV中出现的时间记下来。
                    comp_mat(:,j)=[k;dtime];%comp_mat为2*n的矩阵
                    j=j+1;
                end
        end
                    if isempty(comp_mat);continue;end
                    [wu,ind]=min(comp_mat(2,:));%寻找最小的时间差值所在位置index
                    index=comp_mat(1,ind);  %得到最佳的观测卫星在NAV中的位置；
            if(OBS(epoch).C1(i)<1e-3)%把某个历元中某颗卫星没有的伪距C1去除
                continue
            end
%=============寻找完毕======================================================
%==============计算卫星位置=================================================
                rho2=obs(i)^2;
                t0 = sow - sqrt(rho2)/c;%信号发送时刻初值
                while 1%迭代求卫星位置
                    [x,y,z,delt_sat]=NavPos(t0,index,sqrt(rho2)/c);
                    rho2 = (x-X(1))^2+(y-X(2))^2+(z-X(3))^2;
                    tx_GPS=sow-sqrt(rho2)/c-delt_sat-X(7)/c;%更新卫星信号发送时刻
                    diff=abs(tx_GPS-t0);
                    if diff<1e-15
                        break;
                    end
                    t0=tx_GPS;
                end
%===============卫星位置计算完毕============================================                
                Rot_X(1,1)=x;
                Rot_X(2,1)=y;
                Rot_X(3,1)=z;
                [az,el,dist] = topocent(X(1:3,:),Rot_X-X(1:3,:));%计算卫星高度角   
%============进行对流层改正=================================================
            site.tropmodel='GPT';
            site.tropmap='GMF';
            site.ztdmodel='SAAS';
            [BLH]=transC_xyz2blh(X(1),X(2),X(3));
   %         dtrop=trop_delay(epoch,tt,BLH,el,site);
          %  dtrop=2.47/sind(el)+0.0121;         %简化的hopefiled模型
          dtrop = tropo(sin(el*dtr),0.0,1013.0,293.0,50.0,...
                0.0,0.0,0.0);
                w(i,i)=sind(el);                %观测值权重
%============计算k8电离层改正===============================================
            azel=[az*pi/180,el*pi/180];
            dion=ion_k8_gps(tepo,ion,BLH,azel);
%=================电离层改正完毕============================================
%===========构建关系矩阵====================================================
            omc = [omc; obs(i)-sqrt(rho2)-X(7)+c*delt_sat-dtrop-dion];
            A = [A; (-(Rot_X(1)-X(1)))/sqrt(rho2)...                                            
                (-(Rot_X(2)-X(2)))/sqrt(rho2)...
                (-(Rot_X(3)-X(3)))/sqrt(rho2) 0 0 0 1 0];
%==========关系矩阵构建完毕=================================================
            R(i,i)=4/sind(el);        % 观测噪声协方差阵
    end
%=========卡尔曼滤波方程组==================================================
            A=w*A;
            omc=w*omc;
            X=Phi*X;
            P = Phi*P*Phi'+Qk;
            K = P*A'*inv(A*P*A'+R);
            X = X + K*omc;
            P = (eye(8)-K*A)*P;
%==========滤波结束========================================================           
            POS(epoch,:)=X;
            BLH=transC_xyz2blh(X(1),X(2),X(3));
            B(epoch)=BLH.B;L(epoch)=BLH.L;H(epoch)=BLH.H;
            [b(epoch),l(epoch),h(epoch)]=XYZ2BLH(X(1),X(2),X(3));
            disp(epoch);
end
save('NAV.mat','NAV','Header');
save('OBS.mat','OBS','XYZ');
save('klman_position','b','l','h','POS');
figure;
plot(l,b,'.');
%axis([2.1196812 2.1196840 0.5416814 0.541683]);
%axis([max(l)-0.000008 max(l) max(b)-0.000008 max(b)]);

figure;
plot(1:epoch,(POS(:,1)-X(1)*ones(epoch,1))','-',...
    1:epoch,(POS(:,2)-X(2)*ones(epoch,1))','-',...
    1:epoch,(POS(:,3)-X(3)*ones(epoch,1))','-','linewidth',.25)
title('kinematic Positions over time','fontsize',10)
axis([1 epoch  -8 8]);
legend('X','Y','Z')
xlabel('Epochs [1 s interval]','fontsize',10)
ylabel('Changes in {\itX}, {\itY}, {\itZ} since the first epoch [m]','fontsize',10)
set(gca,'fontsize',10)
set(gca,'XTick',0:500:epoch);
set(gca,'YTick',-8:1:8);
legend










