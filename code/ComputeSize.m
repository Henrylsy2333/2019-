%�˴������˶������Ľṹ����������ָ���켣�ĸ���

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
% op=   [NaN;0  ;1  ;0  ;1  ;0  ;1  ;0  ;1  ;NaN];   1:��ֱУ�� 0:ˮƽУ��
% uc   =[0;1;0;1;0;1;0;1;NaN]; %�켣��Ĳ�ȷ����Ϣ  1:�������� 0��������
% dde  =[]
%npҲȡ1��2��3����ʾInput�������±�
%�������ݵĲ����У��ڵ��Ӧ�ĸ���p
%                 ���ڵ��eh��ev��Ϊ����ǰ���
%%
global node_number
global tree_cell
%�ڵ�1�Ĳ���
node_number = 1;

tree_cell(node_number).p=0;
tree_cell(node_number).ev=0;
tree_cell(node_number).eh=0;
tree_cell(node_number).np=0;
tree_cell(node_number).parent=0;
tree_cell(node_number).record=5;
%%
tree_cell(1).left=build_btree_parent(node_number,1,de(1),de(1),1); %�ڵ���2(parent,p,eh,ev,np)
%�ڵ�2�ĸ����ǽڵ�1����1�м̳еõ�1�ĸ��ʣ�������������ƽ��У׼�ģ����˴��ţ���Ϊeh,ev���ʼ��Ϊ0��npΪ���ڵ�np
%% �������
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
