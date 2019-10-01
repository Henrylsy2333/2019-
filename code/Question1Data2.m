% This is the program of the first question about the second data file.
% You can run it successfully by setting the correct load path as below.
% At the end of this programme, you can get the plot of trajactory.
%%
clc;clear;
%position_raw2 = xlsread('附件2：数据集2-终稿.xlsx','B3:E329');
%save position_raw2
dist2=inline('sum((a-b).^2,2).^(1/2)','a','b'); % 求两点之间的距离 内联函数
%% Data processing
load Raw_Pos2.mat%运行时需要改变路径
position_raw=[position_raw [1:size(position_raw,1)]'];
x=position_raw(:,1);y=position_raw(:,2);z=position_raw(:,3);mark=position_raw(:,4);
vtcpos=position_raw(find(mark==1),:); % 垂直误差校正点
hrzpos=position_raw(find(mark==0),:); % 水平误差校正点
% 垂直校准点集
vx=vtcpos(:,1);
vy=vtcpos(:,2);
vz=vtcpos(:,3);
% 水平校准点集
hx=hrzpos(:,1);
hy=hrzpos(:,2);
hz=hrzpos(:,3);
% 尺度变换因子
% dx=max(x)-min(x);
% dy=max(y)-min(y);
% dz=max(z)-min(z);
dx=1;
dy=1;
dz=1;
%%  问题一附件二参数设定
a1=20;
a2=10;
b1=15;
b2=20;
theta=20;
delta=0.001;
rref = 2.5*10^4; %搜索半径，和密度有关
%蚁群算法参数
N=80; %跑N条轨迹
Q=100000; %全值
G=20; %迭代次数
tau=1.*(ones(327,327)-eye(327));%存储信息素的矩阵,先给一个初值（切忌！不可都给0）
%%
for j=1:G
    PointIndex = zeros(50,N); %每一列用来存放，第n个蚂蚁经过的路径
    L=[]; %L表示返回的数据内容，作修正信息素用！第一列表示长度；第二列表示转折个数
    n=1;
    while n<=N  %多个蚂蚁分别跑，一遍跑完统一更新信息素
        %% 变量初始化
        
        eh=0; % horizen error
        ev=0; % vertical error
        e=zeros(1,2);%（eh,ev）
        crtpoint=zeros(1,3);% current point
        tgtpoint=zeros(1,3);%target point
        dde=0;%要前进的距离
        dflag=0;%目标集 变量值：水平校正集：1，垂直校正集：2
        x0=[0 50000 5000]'; %当前向量坐标
        xend = [100000	74860.5	5499.61]';
        %关于轨迹和点的变量
        ip=2; %index of point
        PointRecord =zeros(50,6);%假设每条轨迹最多有50个点
        PointRecord(1,:)=[x0;0;1;0]';
        %% 蚁群算法参数
        alpha=5; %信息素的重要性程度
        beta=2.5;  %启发式函数的重要性程度   注：可以采用变参数的思想，一开始以启发式函数为主，然后以信息素为主
        current_i=1;
        target_j=1;
        %%
        
        while ~((dist2(x0',xend')*delta+eh)<=theta & (dist2(x0',xend')*delta+ev)<=theta)
            %% 计算理想的前进距离值以及归零的 Flag
            %此时的误差坐标必然是（eh,0）和（0,ev)中的一种情况
            if eh==0 & ev==0
                dde = 15;
                dflag = 1; 
            elseif eh>0 && ev==0
                deh=b2-eh; %水平方向上的误差冗余
                if deh<0  %误差出现了问题
                    disp('Error：de小于0，水平误差超出范围');
                    pause;
                elseif deh>b1 %误差超出了最大值
                    deh = b1;
                end
                dde = deh;
                dflag = 1;   %目标是进行水平校正
            elseif ev>0 && eh==0
                dev=a1-ev; %垂直方向上的误差冗余
                if dev<0  %误差出现了问题
                    disp('Error：de小于0，垂直误差超出范围');
                    pause;
                elseif dev>a2 %误差超出了最大值
                    dev = a2;
                end
                dde = dev;
                dflag = 2;   %目标是进行垂直校正
            end
            % 计算目标点
            r = dde*(1/delta);
            R = sum((xend-x0).^2)^(1/2);
            xe = xend-x0;
            xref = x0 + r/R * eye(3) * xe;
            %% 搜索范围
            if dflag == 2 % 垂直校准点搜索
                ResearchIndex = vtcpos(:,1)>=(xref(1)-rref*dx) & vtcpos(:,1)<=(xref(1)+rref*dx) ...
                    & vtcpos(:,2)>=(xref(2)-rref*dy) & vtcpos(:,2)<=(xref(2)+rref*dy) ...
                    & vtcpos(:,3)>=(xref(3)-rref*dz) & vtcpos(:,3)<=(xref(3)+rref*dz);    %方块滤波器
                
                ResPoint = vtcpos(ResearchIndex,:); %符合要求的垂直搜索点
            elseif dflag == 1 % 水平校准点搜索
                ResearchIndex = hrzpos(:,1)>=(xref(1)-rref*dx) & hrzpos(:,1)<=(xref(1)+rref*dx) ...
                    & hrzpos(:,2)>=(xref(2)-rref*dy) & hrzpos(:,2)<=(xref(2)+rref*dy) ...
                    & hrzpos(:,3)>=(xref(3)-rref*dz) & hrzpos(:,3)<=(xref(3)+rref*dz);
                ResPoint = hrzpos(ResearchIndex,:); %符合要求的水平搜索点
            end
            Filtercir = ResPoint(dist2(ResPoint(:,[1:3]),xref')< rref,:);       %圆球滤波器
            Filterhalf = Filtercir(dist2(Filtercir(:,[1:3]),x0')<= dde*(1/delta),:);  %半球滤波器
            %% 求解搜索区域内每个点对应的概率，并且选择随机点（主要在这里进行修改）
            % 先求各个指标的大小
            RandomPoint=Filterhalf;
            % 基于内积的概率生成
%              dr2o=dist2(RandomPoint(:,[1:3]),x0');
%              dd3=(RandomPoint(:,[1:3])-x0')*(xref-x0);%两者的内积，值越大越接近
%              rdvct=dd3; %random vector
            % 基于夹角大小的概率生成
%             dcos_theta = dd3./(dr2o*dist2(xref',x0'));
%             dcos = [dcos_theta RandomPoint(:,[4,5])];
%             rdvct=dcos; %random vector
            %%使用指标归一化生成概率
            % 基于到目标点距离的概率生成
            dd1=dist2(RandomPoint(:,[1:3]),xref');
            rdvct=1000./dd1; %random vector
            rdvct(find(rdvct(:,1)<=0),1)=0; %让所有cos_theta小于0
            if sum(rdvct(:,1))==0
                break
            end
            randd=rdvct(:,1)./sum(rdvct(:,1)); %启发式函数
            %构造非线性函数法
            %%  蚁群算法计算概率（考虑信息素和启发函数的融合数据）
            % 启发式函数
            xx=[randd RandomPoint(:,5)];
            % 信息
            P_cpt=[tau(current_i,xx(:,2)')',xx];% 第一列式信息素，第二列式启发函数值，第三列是点的序列
            P = ((P_cpt(:,1).^alpha).*(P_cpt(:,2).^beta))./((P_cpt(:,1).^alpha)'*(P_cpt(:,2).^beta));
            Index = min(find((cumsum(P)-rand)>0)); %随机选取，产生对应索引
            PointRecord(ip,:)=[RandomPoint(Index,:),dflag]; %记录随机选取的点的信息
            tgtpoint=RandomPoint(Index,[1:3])';
            current_i = RandomPoint(Index,5);
            %% 更新（eh，ev)
            eh1=0;
            ev1=0;
            do2t=dist2(x0',tgtpoint');%
            if dflag == 1  %水平校正
                eh1 = eh+delta*do2t;
                ev1 = ev+delta*do2t;
                eh = 0;
                ev = ev+delta*do2t;
            elseif dflag == 2   %垂直校正
                ev1 = ev + delta*do2t;
                eh1 = eh+delta*do2t;
                ev = 0;
                eh = eh + delta * do2t;
            end
            et(ip,:)=[eh1,ev1,eh,ev,dflag];
            x0 = tgtpoint;
            ip=ip+1;
        end
        %% 计算到终点的距离
        if (dist2(x0',xend')*delta+eh)<=theta & (dist2(x0',xend')*delta+ev)<=theta
            PointRecord(ip,:)=[xend;0;327;0]';
            %disp('到达终点');
        else
            continue
        end
        %% 计算路径长度和改变次数
        PointRecordA = PointRecord;
        PointRecordA(all(PointRecordA==0,2),:)=[];
        D = pdist(PointRecordA(:,[1:3]));
        T = squareform(D);
        A=diag(ones(1,length(T)-1),1);
        L(n)=sum(sum(A.*T,2)); %该条轨迹的长度
        PointIndex(:,n)=PointRecord(:,5);      
        Lt = L;  
        n = n+1;
    end
    %% 更新信息素
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
    disp(['更新到了第',num2str(j),'次']);
    %disp(['更新到了第',num2str(j),'次']);
    %%
     % tau(a,b)=Q./L
end
disp('结束运行');

%% 可视化程序
figure
plot3([x(1),x(end)],[y(1),y(end)],[z(1),z(end)],'.m','MarkerSize',20);
hold on
plot3(vx,vy,vz,'.r');
plot3(hx,hy,hz,'.b');
grid on
xlabel('x方向(m)');
ylabel('y方向(m)');
zlabel('z方向(m)');
hold on
plot3(PointRecordA(:,1),PointRecordA(:,2),PointRecordA(:,3),'r-');
plot3(PointRecordA(:,1),PointRecordA(:,2),PointRecordA(:,3),'*m','MarkerSize',5);
grid on
%% 蚁群算法迭代图 
figure
plot([1:length(LJ)],LJ,'b');
xlabel('迭代次数');
ylabel('长度(m)');