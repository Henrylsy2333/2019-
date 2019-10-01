%此代码用了二叉树的结构，用来生成指定轨迹的概率

function [Sum]=ComputeSize(Input)
%global  position_raw
global op de uc size_input
 load Raw_Pos13.mat
 position_raw=[position_raw [1:size(position_raw,1)]'];
% x=position_raw(:,1);y=position_raw(:,2);z=position_raw(:,3);mark=position_raw(:,4);
size_input=length(Input);
op=[0;position_raw(Input(2:end),4)];
uc= position_raw(Input,5);
D=pdist(position_raw(Input,[1:3]));
T = squareform(D);
A=diag(ones(1,length(T)-1),1);
dder=sum(A.*T,2);
de=[0;dder(1:end-1)*0.001];
M=[Input op uc de];
%%
% op=   [NaN;0  ;1  ;0  ;1  ;0  ;1  ;0  ;1  ;NaN];   1:垂直校正 0:水平校正
% uc   =[0;1;0;1;0;1;0;1;NaN]; %轨迹点的不确定信息  1:不正常点 0：正常点
% dde  =[]
%np也取1，2，3仅表示Input的索引下标
%函数传递的参数有：节点对应的概率p
%                 父节点的eh和ev，为更新前误差
%%
global node_number
global tree_cell
%节点1的参数
node_number = 1;

tree_cell(node_number).p=0;
tree_cell(node_number).ev=0;
tree_cell(node_number).eh=0;
tree_cell(node_number).np=0;
tree_cell(node_number).parent=0;
tree_cell(node_number).record=5;
%%
tree_cell(1).left=build_btree_parent(node_number,1,de(1),de(1),1); %节点是2(parent,p,eh,ev,np)
%节点2的父亲是节点1，从1中继承得到1的概率，这里假设起点是平行校准的，无伤大雅，因为eh,ev起点始终为0，np为父节点np
%% 计算概率
P = cat(1,tree_cell.p);
Record= cat(1,tree_cell.record);
Parent= cat(1,tree_cell.parent);
clear global 
Sum=0;
%
while ~isempty(find(Record==1))
    PT=find(Record==1);
    PP=PT(1);
    Record(PP)=0;
    Rt=[];
    while PP~=1
        Rt=[Rt PP];
        PP=Parent(PP);
    end
    K=P(Rt);
    PK=prod(K);
    Sum=Sum+PK;
end
