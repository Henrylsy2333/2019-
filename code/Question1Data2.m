% This is the program of the first question about the second data file.
% You can run it successfully by setting the correct load path as below.
% At the end of this programme, you can get the plot of trajactory.
%%
clc;clear;
%position_raw2 = xlsread('����2�����ݼ�2-�ո�.xlsx','B3:E329');
%save position_raw2
dist2=inline('sum((a-b).^2,2).^(1/2)','a','b'); % ������֮��ľ��� ��������
%% Data processing
load Raw_Pos2.mat%����ʱ��Ҫ�ı�·��
position_raw=[position_raw [1:size(position_raw,1)]'];
x=position_raw(:,1);y=position_raw(:,2);z=position_raw(:,3);mark=position_raw(:,4);
vtcpos=position_raw(find(mark==1),:); % ��ֱ���У����
hrzpos=position_raw(find(mark==0),:); % ˮƽ���У����
% ��ֱУ׼�㼯
vx=vtcpos(:,1);
vy=vtcpos(:,2);
vz=vtcpos(:,3);
% ˮƽУ׼�㼯
hx=hrzpos(:,1);
hy=hrzpos(:,2);
hz=hrzpos(:,3);
% �߶ȱ任����
% dx=max(x)-min(x);
% dy=max(y)-min(y);
% dz=max(z)-min(z);
dx=1;
dy=1;
dz=1;
%%  ����һ�����������趨
a1=20;
a2=10;
b1=15;
b2=20;
theta=20;
delta=0.001;
rref = 2.5*10^4; %�����뾶�����ܶ��й�
%��Ⱥ�㷨����
N=80; %��N���켣
Q=100000; %ȫֵ
G=20; %��������
tau=1.*(ones(327,327)-eye(327));%�洢��Ϣ�صľ���,�ȸ�һ����ֵ���мɣ����ɶ���0��
%%
for j=1:G
    PointIndex = zeros(50,N); %ÿһ��������ţ���n�����Ͼ�����·��
    L=[]; %L��ʾ���ص��������ݣ���������Ϣ���ã���һ�б�ʾ���ȣ��ڶ��б�ʾת�۸���
    n=1;
    while n<=N  %������Ϸֱ��ܣ�һ������ͳһ������Ϣ��
        %% ������ʼ��
        
        eh=0; % horizen error
        ev=0; % vertical error
        e=zeros(1,2);%��eh,ev��
        crtpoint=zeros(1,3);% current point
        tgtpoint=zeros(1,3);%target point
        dde=0;%Ҫǰ���ľ���
        dflag=0;%Ŀ�꼯 ����ֵ��ˮƽУ������1����ֱУ������2
        x0=[0 50000 5000]'; %��ǰ��������
        xend = [100000	74860.5	5499.61]';
        %���ڹ켣�͵�ı���
        ip=2; %index of point
        PointRecord =zeros(50,6);%����ÿ���켣�����50����
        PointRecord(1,:)=[x0;0;1;0]';
        %% ��Ⱥ�㷨����
        alpha=5; %��Ϣ�ص���Ҫ�Գ̶�
        beta=2.5;  %����ʽ��������Ҫ�Գ̶�   ע�����Բ��ñ������˼�룬һ��ʼ������ʽ����Ϊ����Ȼ������Ϣ��Ϊ��
        current_i=1;
        target_j=1;
        %%
        
        while ~((dist2(x0',xend')*delta+eh)<=theta & (dist2(x0',xend')*delta+ev)<=theta)
            %% ���������ǰ������ֵ�Լ������ Flag
            %��ʱ����������Ȼ�ǣ�eh,0���ͣ�0,ev)�е�һ�����
            if eh==0 & ev==0
                dde = 15;
                dflag = 1; 
            elseif eh>0 && ev==0
                deh=b2-eh; %ˮƽ�����ϵ��������
                if deh<0  %������������
                    disp('Error��deС��0��ˮƽ������Χ');
                    pause;
                elseif deh>b1 %���������ֵ
                    deh = b1;
                end
                dde = deh;
                dflag = 1;   %Ŀ���ǽ���ˮƽУ��
            elseif ev>0 && eh==0
                dev=a1-ev; %��ֱ�����ϵ��������
                if dev<0  %������������
                    disp('Error��deС��0����ֱ������Χ');
                    pause;
                elseif dev>a2 %���������ֵ
                    dev = a2;
                end
                dde = dev;
                dflag = 2;   %Ŀ���ǽ��д�ֱУ��
            end
            % ����Ŀ���
            r = dde*(1/delta);
            R = sum((xend-x0).^2)^(1/2);
            xe = xend-x0;
            xref = x0 + r/R * eye(3) * xe;
            %% ������Χ
            if dflag == 2 % ��ֱУ׼������
                ResearchIndex = vtcpos(:,1)>=(xref(1)-rref*dx) & vtcpos(:,1)<=(xref(1)+rref*dx) ...
                    & vtcpos(:,2)>=(xref(2)-rref*dy) & vtcpos(:,2)<=(xref(2)+rref*dy) ...
                    & vtcpos(:,3)>=(xref(3)-rref*dz) & vtcpos(:,3)<=(xref(3)+rref*dz);    %�����˲���
                
                ResPoint = vtcpos(ResearchIndex,:); %����Ҫ��Ĵ�ֱ������
            elseif dflag == 1 % ˮƽУ׼������
                ResearchIndex = hrzpos(:,1)>=(xref(1)-rref*dx) & hrzpos(:,1)<=(xref(1)+rref*dx) ...
                    & hrzpos(:,2)>=(xref(2)-rref*dy) & hrzpos(:,2)<=(xref(2)+rref*dy) ...
                    & hrzpos(:,3)>=(xref(3)-rref*dz) & hrzpos(:,3)<=(xref(3)+rref*dz);
                ResPoint = hrzpos(ResearchIndex,:); %����Ҫ���ˮƽ������
            end
            Filtercir = ResPoint(dist2(ResPoint(:,[1:3]),xref')< rref,:);       %Բ���˲���
            Filterhalf = Filtercir(dist2(Filtercir(:,[1:3]),x0')<= dde*(1/delta),:);  %�����˲���
            %% �������������ÿ�����Ӧ�ĸ��ʣ�����ѡ������㣨��Ҫ����������޸ģ�
            % �������ָ��Ĵ�С
            RandomPoint=Filterhalf;
            % �����ڻ��ĸ�������
%              dr2o=dist2(RandomPoint(:,[1:3]),x0');
%              dd3=(RandomPoint(:,[1:3])-x0')*(xref-x0);%���ߵ��ڻ���ֵԽ��Խ�ӽ�
%              rdvct=dd3; %random vector
            % ���ڼнǴ�С�ĸ�������
%             dcos_theta = dd3./(dr2o*dist2(xref',x0'));
%             dcos = [dcos_theta RandomPoint(:,[4,5])];
%             rdvct=dcos; %random vector
            %%ʹ��ָ���һ�����ɸ���
            % ���ڵ�Ŀ������ĸ�������
            dd1=dist2(RandomPoint(:,[1:3]),xref');
            rdvct=1000./dd1; %random vector
            rdvct(find(rdvct(:,1)<=0),1)=0; %������cos_thetaС��0
            if sum(rdvct(:,1))==0
                break
            end
            randd=rdvct(:,1)./sum(rdvct(:,1)); %����ʽ����
            %��������Ժ�����
            %%  ��Ⱥ�㷨������ʣ�������Ϣ�غ������������ں����ݣ�
            % ����ʽ����
            xx=[randd RandomPoint(:,5)];
            % ��Ϣ
            P_cpt=[tau(current_i,xx(:,2)')',xx];% ��һ��ʽ��Ϣ�أ��ڶ���ʽ��������ֵ���������ǵ������
            P = ((P_cpt(:,1).^alpha).*(P_cpt(:,2).^beta))./((P_cpt(:,1).^alpha)'*(P_cpt(:,2).^beta));
            Index = min(find((cumsum(P)-rand)>0)); %���ѡȡ��������Ӧ����
            PointRecord(ip,:)=[RandomPoint(Index,:),dflag]; %��¼���ѡȡ�ĵ����Ϣ
            tgtpoint=RandomPoint(Index,[1:3])';
            current_i = RandomPoint(Index,5);
            %% ���£�eh��ev)
            eh1=0;
            ev1=0;
            do2t=dist2(x0',tgtpoint');%
            if dflag == 1  %ˮƽУ��
                eh1 = eh+delta*do2t;
                ev1 = ev+delta*do2t;
                eh = 0;
                ev = ev+delta*do2t;
            elseif dflag == 2   %��ֱУ��
                ev1 = ev + delta*do2t;
                eh1 = eh+delta*do2t;
                ev = 0;
                eh = eh + delta * do2t;
            end
            et(ip,:)=[eh1,ev1,eh,ev,dflag];
            x0 = tgtpoint;
            ip=ip+1;
        end
        %% ���㵽�յ�ľ���
        if (dist2(x0',xend')*delta+eh)<=theta & (dist2(x0',xend')*delta+ev)<=theta
            PointRecord(ip,:)=[xend;0;327;0]';
            %disp('�����յ�');
        else
            continue
        end
        %% ����·�����Ⱥ͸ı����
        PointRecordA = PointRecord;
        PointRecordA(all(PointRecordA==0,2),:)=[];
        D = pdist(PointRecordA(:,[1:3]));
        T = squareform(D);
        A=diag(ones(1,length(T)-1),1);
        L(n)=sum(sum(A.*T,2)); %�����켣�ĳ���
        PointIndex(:,n)=PointRecord(:,5);      
        Lt = L;  
        n = n+1;
    end
    %% ������Ϣ��
    p=0.02;
    tau = (1-p)*tau;
    PointIndex2 = [PointIndex([2:end],:);zeros(1,N)];
    for n=1:N
        temp=PointIndex(find(PointIndex(:,n)~=0),n);
        for i=1:length(temp)-1
            tau(temp(i),temp(i+1))= tau(temp(i),temp(i+1))+Q./Lt(n);
        end
    end
    LJ(j)=mean(Lt);
    disp(['���µ��˵�',num2str(j),'��']);
    %disp(['���µ��˵�',num2str(j),'��']);
    %%
     % tau(a,b)=Q./L
end
disp('��������');

%% ���ӻ�����
figure
plot3([x(1),x(end)],[y(1),y(end)],[z(1),z(end)],'.m','MarkerSize',20);
hold on
plot3(vx,vy,vz,'.r');
plot3(hx,hy,hz,'.b');
grid on
xlabel('x����(m)');
ylabel('y����(m)');
zlabel('z����(m)');
hold on
plot3(PointRecordA(:,1),PointRecordA(:,2),PointRecordA(:,3),'r-');
plot3(PointRecordA(:,1),PointRecordA(:,2),PointRecordA(:,3),'*m','MarkerSize',5);
grid on
%% ��Ⱥ�㷨����ͼ 
figure
plot([1:length(LJ)],LJ,'b');
xlabel('��������');
ylabel('����(m)');