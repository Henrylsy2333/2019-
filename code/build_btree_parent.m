% Input=[1  ;522;65 ;81 ;171;279;370;215;398;613]
% op=   [NaN;0  ;1  ;0  ;1  ;0  ;1  ;0  ;1  ;NaN];   1:��ֱУ�� 0:ˮƽУ��
% uc   =[0;1;0;1;0;1;0;1;NaN]; %�켣��Ĳ�ȷ����Ϣ  1:�������� 0��������
% dde  =[]
%npҲȡ1��2��3����ʾInput�������±�
%�������ݵĲ����У��ڵ��Ӧ�ĸ���p
%                 ���ڵ��eh��ev��Ϊ����ǰ���

function tree_out=build_btree_parent(parent,p,eh,ev,np)
%np��ʾ��һ���ڵ�
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
    %% ��������    
    node_number=node_number+1;
%     assigned_nn=node_number;% ���µĽڵ�
%     parent=assigned_nn;
    tree_cell(node_number).p=p;        %��һ��ȷ�����ȷ��������ĸ��ʣ���������������
    tree_cell(node_number).ev=ev;      %��һ���켣������������δ���´�ֱ���
    tree_cell(node_number).eh=eh;      %��һ���켣������������δ����ˮƽ���
    tree_cell(node_number).np=np;      %��һ���켣��
    tree_cell(node_number).parent=parent;
%end
assigned_nn=node_number;
    %% Step1�����չ�����ͣ������ʸ��¹������
    if op(np)==0 %|| np==2%ˮƽУ����Ҫ�ӵ�np-1���㿪ʼ�����ǵ�һ�����
        %if np>2
            if p==1 || p==0.8%���ڵ�ĸ��������㣬���շ���ĸ���
                eht=0;
            elseif p==0.2
                eht=min(5,eh);
            end
         np=np+1; %�������µĵ�
        eh=eht+de(np); %�µĵ����ǰ���
        ev=ev+de(np);
        %% ��ֱУ��
    elseif op(np)==1 
        if p==1 || p==0.8%���ڵ�ĸ��������㣬���շ���ĸ���
            evt=0;
        elseif p==0.2
            evt=min(5,eh);
        end
        np=np+1; %�������µĵ�
        ev=evt+de(np); %�µĵ����ǰ���
        eh=eh+de(np);
    end
    tree_cell(node_number).ev=ev;      %��һ���켣������������δ���´�ֱ���
    tree_cell(node_number).eh=eh;      %��һ���켣������������δ����ˮƽ���
%% �ж��Ƿ񵽴��յ�
        if np==size_input %
            if eh<=theta & ev<=theta
                tree_out=node_number;
                %disp(['�ɹ������յ�'])
                tree_cell(node_number).record = 1;
                return;       %�ɹ������յ㣬�����ǵ�������ֵ����ô����ʱ�����
            else             %�����յ�ʧ��
                %disp(['�����յ�ʧ��'])
                tree_out=5555;
                tree_cell(node_number).record = 3;
                return;
            end
        elseif np<size_input
            if op(np)==0  %�õ�ΪˮƽУ��   �ж��µĵ�Ĵ�ֱ��ƽ������
                if ~(ev<=b1 & eh<=b2) % ���������ˮƽԼ��
                    tree_out=node_number;
                    tree_cell(node_number).record = 4;
                    return;
                end
            elseif op(np)==1 %�õ�Ϊ��ֱУ��
                if ~(ev<=a1 & eh<=a2) % ��������㴹ֱԼ��
                    tree_out=node_number;
                    tree_cell(node_number).record = 2;
                    return;%����͵�ͷ�ˣ�����
                end
            end
            tree_cell(node_number).record = 5;
            %���ո��ʴ�С������֧
            if uc(np)==1 %��������
                tree_cell(assigned_nn).left=build_btree_parent(assigned_nn,0.8,eh,ev,np);
                tree_cell(assigned_nn).right=build_btree_parent(assigned_nn,0.2,eh,ev,np);
            elseif uc(np)==0  %������
                tree_cell(assigned_nn).left=build_btree_parent(assigned_nn,1,eh,ev,np);
            end
        end
    tree_out=assigned_nn;
    %tree_cell(node_number).record = 5;


%%


