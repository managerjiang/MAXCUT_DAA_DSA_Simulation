function V=local_com_dis(inum,V2,V,W,J,Rset,Sparset)

m=6;
num=139;
local_n=9600;
theta=1;
for i=1:local_n
    if (1/norm(W(i,:),1)<theta)
        theta=1/norm(W(i,:),1);
    end
end
%     V=cell(m,1);
%     for i=1:m
%         V{i,1}=V2{i,1};%%%个体i的第1次迭代
%     end
    
    if Rset(inum,1)~=0
        for j=1:sum(Rset(inum,:)~=0) %%%对于余集部分向量
            sum_inum=zeros(num,1);
            for chi=1:inum  %%%%子集和自身
                if ismember(Rset(inum,j),J(chi,:))%%%判断是否在子集chi中
                    J_unique=unique(J(chi,1:sum(J(chi,:)~=0)));
                     for chij=1:length(J_unique) 
                         sum_inum=sum_inum+W(Rset(inum,j),J_unique(1,chij))*V2{chi,1}(:,J_unique(1,chij));
                     end
                end
            end
            interV=V2{inum,1}(:,Rset(inum,j))-theta*sum_inum;
            if (norm(interV)==0)
                V{inum,1}(:,Rset(inum,j))=zeros(num,1);
            else
                V{inum,1}(:,Rset(inum,j))=interV/norm(interV);
            end
        end
    end
    if Sparset(inum,1)~=0
        for j=1:sum(Sparset(inum,:)~=0)  %%%对于交集部分向量
            for pari=inum:m
                if ismember(Sparset(inum,j),J(pari,:))%%%判断向量是否在父辈pari中               
                    V{inum,1}(:,Sparset(inum,j))=V2{pari,1}(:,Sparset(inum,j));
                end               
            end
        end 
    end
end