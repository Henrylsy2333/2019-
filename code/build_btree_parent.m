% Input=[1  ;522;65 ;81 ;171;279;370;215;398;613]
% op=   [NaN;0  ;1  ;0  ;1  ;0  ;1  ;0  ;1  ;NaN];   1:垂直校正 0:水平校正
% uc   =[0;1;0;1;0;1;0;1;NaN]; %轨迹点的不确定信息  1:不正常点 0：正常点
% dde  =[]
%np也取1，2，3仅表示Input的索引下标
%函数传递的参数有：节点对应的概率p
%                 父节点的eh和ev，为更新前误差

function tree_out=build_btree_parent(parent,p,eh,ev,np)
%np表示上一个节点
global node_number
global tree_cell size_input
global op de uc
theta=30;
a1=25;
a2=15;
b1=20;
b2=25;
% a1=20;
% a2=10;
% b1=15;
% b2=20;
% theta=20;
%delta=0.001;
    %% 更新数据    
    node_number=node_number+1;
%     assigned_nn=node_number;% 最新的节点
%     parent=assigned_nn;
    tree_cell(node_number).p=p;        %上一个确定点或不确定点决定的概率，用来进行误差更新
    tree_cell(node_number).ev=ev;      %上一个轨迹点遗留下来的未更新垂直误差
    tree_cell(node_number).eh=eh;      %上一个轨迹点遗留下来的未更新水平误差
    tree_cell(node_number).np=np;      %上一个轨迹点
    tree_cell(node_number).parent=parent;
%end
assigned_nn=node_number;
    %% Step1：按照轨点类型，按概率更新轨点的误差
    if op(np)==0 %|| np==2%水平校正，要从第np-1个点开始，才是第一个轨点
        %if np>2
            if p==1 || p==0.8%父节点的更新误差计算，按照分配的概率
                eht=0;
            elseif p==0.2
                eht=min(5,eh);
            end
         np=np+1; %到达了新的点
        eh=eht+de(np); %新的点更新前误差
        ev=ev+de(np);
        %% 垂直校正
    elseif op(np)==1 
        if p==1 || p==0.8%父节点的更新误差计算，按照分配的概率
            evt=0;
        elseif p==0.2
            evt=min(5,eh);
        end
        np=np+1; %到达了新的点
        ev=evt+de(np); %新的点更新前误差
        eh=eh+de(np);
    end
    tree_cell(node_number).ev=ev;      %上一个轨迹点遗留下来的未更新垂直误差
    tree_cell(node_number).eh=eh;      %上一个轨迹点遗留下来的未更新水平误差
%% 判断是否到达终点
        if np==size_input %
            if eh<=theta & ev<=theta
                tree_out=node_number;
                %disp(['成功到达终点'])
                tree_cell(node_number).record = 1;
                return;       %成功到达终点，这里是迭代返回值，怎么用暂时不清楚
            else             %到达终点失败
                %disp(['到达终点失败'])
                tree_out=5555;
                tree_cell(node_number).record = 3;
                return;
            end
        elseif np<size_input
            if op(np)==0  %该点为水平校正   判断新的点的垂直，平行类型
                if ~(ev<=b1 & eh<=b2) % 如果不满足水平约束
                    tree_out=node_number;
                    tree_cell(node_number).record = 4;
                    return;
                end
            elseif op(np)==1 %该点为垂直校正
                if ~(ev<=a1 & eh<=a2) % 如果不满足垂直约束
                    tree_out=node_number;
                    tree_cell(node_number).record = 2;
                    return;%这里就到头了，凉凉
                end
            end
            tree_cell(node_number).record = 5;
            %按照概率大小创建分支
            if uc(np)==1 %不正常点
                tree_cell(assigned_nn).left=build_btree_parent(assigned_nn,0.8,eh,ev,np);
                tree_cell(assigned_nn).right=build_btree_parent(assigned_nn,0.2,eh,ev,np);
            elseif uc(np)==0  %正常点
                tree_cell(assigned_nn).left=build_btree_parent(assigned_nn,1,eh,ev,np);
            end
        end
    tree_out=assigned_nn;
    %tree_cell(node_number).record = 5;


%%


