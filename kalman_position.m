clc;
clear all
global  NAV c GM OMEGAe dtr %����ȫ�ֱ�����������������λ��
c=2.997925*10^8;  %����
GM=3.986004418e+14;        %����������������
OMEGAe=7.292115*10^(-5);   %����������ת���ٶȣ������ƣ�
dtr = pi/180;              %�Ƕ�ת����
%==========�������ļ����汾��302============================================
if exist('NAV.mat')~=2
    [NAV,Header]= ReadNavfile_302;
else
    load('NAV.mat');
end
%==========���۲�ֵ�ļ����汾��302=====================================================
if exist('OBS.mat')~=2
    [OBS,XYZ]= Read_17O_files();
else
    load('OBS.mat')
end
%===========��ʼ��λ=======================================================
POS=zeros(length(OBS),8);                  %����ʼ����װ���
dt=0.1;                                    %��Ԫ����
Phic=[1 dt;0 1];                           %�Ӳ�ת�ƾ���
Phi=[eye(3) dt*eye(3) zeros(3,2);...       
    zeros(3,3) eye(3) zeros(3,2);...
    zeros(2,6) Phic];                                     %״̬ת�ƾ���
Sv=[10^-8 10^-8 10^-8];St=0.05;Sf=0.05;     %�ٶ�����,�Ӳ���������Ư�����������ܶ�
Sv=diag(Sv);
Qxyz=[Sv*dt^3 Sv*dt^2/2;...
      Sv*dt^2/2 Sv*dt];                        %���������Ӧ��Э�������
Qc=[St*dt+Sf*dt^3/3 Sf*dt^2/2;Sf*dt^2/2 Sf*dt];%���ջ�ʱ�ӵĹ�������Э�������
Qk=[Qxyz zeros(6,2);...                          %��������
    zeros(2,6) Qc];
for epoch=1:length(OBS)
%================�޳���ǰ��Ԫ�£�.N�в����ڵ�����.O�ļ��д��ڵ�����===========
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
%===============�޳����====================================================
    sow=OBS(epoch).sow;                 %second of week
    ion=Header.ion_gps;                 %���������ӳٲ���
    tepo=OBS(epoch).tepo;               %jd��sod��Ҳ��������k8������ӳ�
    tt(1)=OBS(epoch).sow;               %��������������ӳ�
    tt(2)=OBS(epoch).www;
    m = size(obs,1);                    % ��ǰ��Ԫ��������Ŀ
    A=[];                               %��С�������G����G*[delta_x]=b
    omc=[];                             %��С�������b����
    w=[];                               %Ȩ����
    R=[];                               %�۲�ֵ������
    if epoch==1                %��С���˷�ȷ�ϳ�ʼ״̬������
             [Q,coor]=recpo_ls_klman(obs,sats,sow,ion,tepo);
             sigma0=900/(Q(1,1)+Q(2,2)+Q(3,3));
             P=[sigma0*Q(1,1) sigma0*Q(2,2) sigma0*Q(3,3)...
                2*sigma0*Q(1,1)/dt^2 2*sigma0*Q(2,2)/dt^2 2*sigma0*Q(3,3)/dt^2 ...
                sigma0*Q(4,4) 2*sigma0*Q(4,4)/dt^2];%P�������������--�������˲��㷨��GPS�ǲ���λ���ܵ��㶨λ�е�Ӧ���о�_�.caj
             P=diag(P); 
             X=[coor(1) coor(2) coor(3)  0 0 0 coor(4) 0]';%����ʼ���꣬�ٶȺ���Ư�����㡣
     end
    for i=1:m
%================��NAV��Ѱ�������˵ĵ�������=================================
        j=1;comp_mat=[];
        for k=1:size(NAV,2)
                if(NAV(k).prn==sats(i))%��NAV���ҵ������ǡ�
                    dtime=abs(NAV(k).toe-OBS(epoch).sow);%�Ѹ�������NAV�г��ֵ�ʱ���������
                    comp_mat(:,j)=[k;dtime];%comp_matΪ2*n�ľ���
                    j=j+1;
                end
        end
                    if isempty(comp_mat);continue;end
                    [wu,ind]=min(comp_mat(2,:));%Ѱ����С��ʱ���ֵ����λ��index
                    index=comp_mat(1,ind);  %�õ���ѵĹ۲�������NAV�е�λ�ã�
            if(OBS(epoch).C1(i)<1e-3)%��ĳ����Ԫ��ĳ������û�е�α��C1ȥ��
                continue
            end
%=============Ѱ�����======================================================
%==============��������λ��=================================================
                rho2=obs(i)^2;
                t0 = sow - sqrt(rho2)/c;%�źŷ���ʱ�̳�ֵ
                while 1%����������λ��
                    [x,y,z,delt_sat]=NavPos(t0,index,sqrt(rho2)/c);
                    rho2 = (x-X(1))^2+(y-X(2))^2+(z-X(3))^2;
                    tx_GPS=sow-sqrt(rho2)/c-delt_sat-X(7)/c;%���������źŷ���ʱ��
                    diff=abs(tx_GPS-t0);
                    if diff<1e-15
                        break;
                    end
                    t0=tx_GPS;
                end
%===============����λ�ü������============================================                
                Rot_X(1,1)=x;
                Rot_X(2,1)=y;
                Rot_X(3,1)=z;
                [az,el,dist] = topocent(X(1:3,:),Rot_X-X(1:3,:));%�������Ǹ߶Ƚ�   
%============���ж��������=================================================
            site.tropmodel='GPT';
            site.tropmap='GMF';
            site.ztdmodel='SAAS';
            [BLH]=transC_xyz2blh(X(1),X(2),X(3));
   %         dtrop=trop_delay(epoch,tt,BLH,el,site);
          %  dtrop=2.47/sind(el)+0.0121;         %�򻯵�hopefiledģ��
          dtrop = tropo(sin(el*dtr),0.0,1013.0,293.0,50.0,...
                0.0,0.0,0.0);
                w(i,i)=sind(el);                %�۲�ֵȨ��
%============����k8��������===============================================
            azel=[az*pi/180,el*pi/180];
            dion=ion_k8_gps(tepo,ion,BLH,azel);
%=================�����������============================================
%===========������ϵ����====================================================
            omc = [omc; obs(i)-sqrt(rho2)-X(7)+c*delt_sat-dtrop-dion];
            A = [A; (-(Rot_X(1)-X(1)))/sqrt(rho2)...                                            
                (-(Rot_X(2)-X(2)))/sqrt(rho2)...
                (-(Rot_X(3)-X(3)))/sqrt(rho2) 0 0 0 1 0];
%==========��ϵ���󹹽����=================================================
            R(i,i)=4/sind(el);        % �۲�����Э������
    end
%=========�������˲�������==================================================
            A=w*A;
            omc=w*omc;
            X=Phi*X;
            P = Phi*P*Phi'+Qk;
            K = P*A'*inv(A*P*A'+R);
            X = X + K*omc;
            P = (eye(8)-K*A)*P;
%==========�˲�����========================================================           
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










