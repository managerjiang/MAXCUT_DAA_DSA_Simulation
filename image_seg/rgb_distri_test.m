clc;
close all;
clear;
% global Rset;
% global Sparset;
% global J;
% global m;
% global num;
% global W;
% global theta;
%����ԭͼ��
image = imread('chur_b.jpg');%%chur.jpg  320*480
s = size(image);
%s=[98,97,3],��ʾimage��3��98*97�ľ���,�ֱ�������ͼ��ÿ�����ص��R��G��Bֵ
%�ֱ��ȡRGB
% image_r = image(:,:,1);
% image_g = image(:,:,2);
% image_b = image(:,:,3);
%����RGB���
% subplot(2,2,1),imshow(image_r),title('Red component');  
% subplot(2,2,2),imshow(image_g),title('green component');  
% subplot(2,2,3),imshow(image_g),title('blue component');  

imshow(image,'InitialMagnification','fit');
%%title('original image');
%%%���ؾ���
image=im2double(image);
pixl=cell(s(1),s(2));

for i=1:s(1)
    for j=1:s(2)        
        for d=1:s(3)
        pixl{i,j}(1,d)=image(i,j,d);
      
        end
    end
end
pixl_vec=cell(1,s(1)*s(2));
for i=1:s(1)
    for j=1:s(2)
        pixl_vec{1,(i-1)*s(2)+j}=pixl{i,j};
    end
end
iblock=16;
 num=s(1)/iblock*s(2);%%%ͼ��ֿ� ά��Ϊ9600
cs=cell(1,iblock);
for i=1:iblock
    cs{1,i}=zeros(num,1);
end
parfor_progress(iblock);

parfor block=1:iblock
    block_before =(block-1)*(s(1)/iblock*s(2));
    block_start=(block-1)*(s(1)/iblock*s(2))+1;
    block_end=block*(s(1)/iblock*s(2));
    
       
    t=0.1;%%%%%%0.1  red_brid1=0.15
    W=zeros(num,num); %%%%%%%%%%%%
    %'LineWidth',G.Edges.Weight*10);%%%%%������������ͼ   *********��ֵt
     %%ew(vi; vj) = (2 �� [||rgb(vi) - rgb(vj)|| > t] - 1) �� ||rgb(vi) - rgb(vj)||:
      for i=block_start:block_end%%%%%%W
             if i<=block_end-s(2)
                if rem(i,s(2))==0
                    W(i-block_before,i-block_before+s(2))=(2*(norm(pixl_vec{1,i}-pixl_vec{1,i+s(2)})>t)-1)*norm(pixl_vec{1,i}-pixl_vec{1,i+s(2)});
                    W(i-block_before+s(2),i-block_before)=W(i-block_before,i-block_before+s(2));
                else
                   W(i-block_before,i-block_before+1)=(2*(norm(pixl_vec{1,i}-pixl_vec{1,i+1})>t)-1)*norm(pixl_vec{1,i}-pixl_vec{1,i+1});
                   W(i-block_before,i-block_before+s(2))=(2*(norm(pixl_vec{1,i}-pixl_vec{1,i+s(2)})>t)-1)*norm(pixl_vec{1,i}-pixl_vec{1,i+s(2)});
                   W(i-block_before+1,i-block_before)=W(i-block_before,i-block_before+1);
                   W(i-block_before+s(2),i-block_before)=W(i-block_before,i-block_before+s(2));
                end
            else
                if i<block_end
                    W(i-block_before,i-block_before+1)=(2*(norm(pixl_vec{1,i}-pixl_vec{1,i+1})>t)-1)*norm(pixl_vec{1,i}-pixl_vec{1,i+1});
                    W(i-block_before+1,i-block_before)=W(i-block_before,i-block_before+1);
                end
            end
       end
        figure(2);
%         names=(1:num);
%         namestr=cell(1,num);
%         for i=1:num
%         namestr{i}=num2str(names(i));
%         end
        G=graph(W);
        G.Edges;
        G.Nodes;
        plot(G);
         title('��������ͼ');
        %%%%%%�����طָ�������� *********��������

        m=6;%%%%agent number  6
        local_n=num/m;%
        J=zeros(m,num);
        indi=zeros(m,num);
        for i=1:num
            if(rem(i,local_n)==1)
                local_i=1;%%%%J����Ԫ������
            else
                local_i=local_i+1;
            end
            c=0;
            if (rem(i,local_n)==0)%%%������
                c=fix(i/local_n);
            else
                c=fix(i/local_n)+1;
            end
            J(c,local_i)=i;%%%%�ٶ�M����û��ȫ���У���J���ϵ�һ��Ԫ��������ص�c_i
            for j=1:num
                if j==1
                    indi(c,i-local_n*(c-1))=local_i;
                end
                if W(i,j)~=0&&(j~=i)
                    local_i=local_i+1;
                    J(c,local_i)=j;  %%%%��һ������c����ڼ������壬�ڶ�����������������Ԫ��
                end
            end
            indi(c,local_n+1)=local_i;
        end
        Sset=cell(m,1);%%%agent֮���Ԫ�ؽ���
        for i=1:m
            Sset{i,1}=zeros(m,m);
        end
        for i=1:m
            for j=1:m
                if (j~=i)
                    insec=intersect(J(i,:),J(j,:));
                    Sset{i,1}(j,1:size(insec,2))=insec;
                end
            end
        end

        Rset=zeros(m,num);
        Sparset=zeros(m,num);
        for i=1:m  %%%���˸���m��
            %if i<m       %%%��m��û�и���
                uset=Sset{i,1}(i,:); %%%0Ԫ�ع��ɵ�����
                for j=i+1:m  %%%%ֻ���������븸���޽���
                    uset=union(uset,Sset{i,1}(j,:));
                end
                uedi=setdiff(J(i,:),uset);
                Rset(i,1:length(uedi))=uedi;%%%%�븸�����༯
            %end

        %     if i>1 %%%��һ������û���ӱ�
        %         cset=Sset{i,1}(i,:); %%%0Ԫ�ع��ɵ�����
        %         for j=1:i-1
        %             cset=union(cset,Sset{i,1}(j,:));
        %         end
        %         cedi=setdiff(Rset(i,:),cset);
        %         Schset(i,1:length(cedi))=cedi;%%%%%%%%%%%%���ӱ��Ľ���
        %     end
        %     
        %     ll=setdiff(Rset(i,:),Schset(i,:));%%%%������Ԫ��=�븸�����༯��ȥ�ӱ�����
        %     Lset(i,1:length(ll))=ll;

            spar=setdiff(J(i,:),Rset(i,:));
            Sparset(i,1:length(spar))=spar;   %%%%�븸���Ľ�����
        end
        iterations=60;%%%60
        theta=1;
        for i=1:num
            if (1/norm(W(i,:),1)<theta)
                theta=1/norm(W(i,:),1);
            end
        end
        V=cell(m,1);
        V_2=cell(m,1);
        for i=1:m
            V{i,1}=rand(139,num);%%%����i�ĵ�1�ε���
            V_2{i,1}=V{i,1};
        end
        %%%%%%%%%%%%%%%%%%%%���п�ʼ
        for k=1:iterations
            inum=randperm(m,1);%%ag(mod(k,m)+1);
            V=local_com_dis(inum,V_2,V,W,J,Rset,Sparset);
            for i=1:m
                V_2{i,1}=V{i,1};
            end
        end
        V_variable=zeros(139,num);

        for i=1:m
           if Rset(i,1)~=0
               for j=1:sum(Rset(i,:)~=0)
                   V_variable(:,Rset(i,j))=V{i,1}(:,Rset(i,j));
               end
           end
        end
        rr=rand(139,1);
        rrv=rr/norm(rr);
        %cs=zeros(num,1);
        for i=1:num
            cs{1,block}(i,1)=dot(rrv, V_variable(:,i))>0;
        end
        parsave(sprintf('output%d.mat', block), cs{1,block});
        parfor_progress;
end
parfor_progress(0);
%%%ff=trace(W'*(V_variable'*V_variable));         
for block=1:iblock
    block_start=(block-1)*(s(1)/iblock*s(2))+1;
    block_end=block*(s(1)/iblock*s(2));
    
    Aa=reshape(cs{1,block},s(2),s(1)/iblock)';   %%%%%%%%%%%%************reshape(cs,s(2),s(1))'
    for i=1:s(1)/iblock
        for j=1:s(2) 
            if Aa(i,j)==0
                for d=1:s(3)
                pixl{(block-1)*(s(1)/iblock)+i,j}(1,d)=0;
                end
            else
                for d=1:s(3)
                pixl{(block-1)*(s(1)/iblock)+i,j}(1,d)=255;%image(i,j,d);
                end
            end
        end
    end
    for i=(block-1)*(s(1)/iblock)+1:(s(1)/iblock)*block
        for j=1:s(2)
            for d=1:s(3)
                image(i,j,d)=pixl{i,j}(1,d);
            end
        end
    end
end
figure(3);
imshow(im2uint8(image),'InitialMagnification','fit');




