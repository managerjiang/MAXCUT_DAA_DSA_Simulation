clc;
clear;
close all;
global A;
global bi;
global yv;
global sigma;
global C;
global n;
global C1;
 global C2;
 global C3;
 global C4;
 global C5;
 global C6;
 global m;
 global m1;
 global m2;
 global m3;
 global m4;
 global m5;
 global m2f1;
 global m2f2;
 global iteration;
 
  num=8;
 v1=rand(num,num);
 v1_update=cell(1,iteration);
 v1_update{1,1}=v1;
 
v2=v1;
 v2_update=cell(1,iteration);
 v2_update{1,1}=v2;
 
 v3=v1;
 v3_update=cell(1,iteration);
 v3_update{1,1}=v3;

  v4=v1;
 v4_update=cell(1,iteration);
 v4_update{1,1}=v4;
 
  v5=v1;
 v5_update=cell(1,iteration);
 v5_update{1,1}=v5;
  
  v6=v1;
 v6_update=cell(1,iteration);
 v6_update{1,1}=v6;
 
 
 
C1=zeros(8,8);
C1(1,3)=1/3;
C1(3,1)=1/3;


 C2=[  0   1/2   0   1/4  0 0 0 0
       1/2   0   0  1/3   0 0 0 0
       0     0   0    0   0 0 0 0
       1/4  1/3 0    0    0 0 0 0
       0     0   0    0  0 0 0 0
       0     0   0    0  0 0 0 0
       0     0   0    0  0 0 0 0
       0     0   0    0  0 0 0 0
     ];
 C3=[  0     0  0   0  0 0 0 0
       0     0  0   0  0 0 0 0
       0     0  0   0  0 0 0 0
       0     0  0   0  1/2 0 0 0
       0     0  0   1/2  0 0 0 0
       0     0   0    0  0 0 0 0
       0     0   0    0  0 0 0 0
       0     0   0    0  0 0 0 0
     ];
 C4=[  0     0  0    0  0 0 0 0
       0     0  0    0  0 0 0 0
       0     0  0    1/4  0 0 0 0
       0     0  1/4  0   0 0 0 0
       0     0  0    0  0 0 0 0
       0     0   0   0  0 0 0 0
       0     0   0   0  0 0 0 0
       0     0   0   0  0 0 0 0
     ];
  C5=[  0     0  0    0   0  0  0 0
       0     0  0    0    0  0  0 0
       0     0  0    0    0 1/7 1/2 0
       0     0  0    0    0  0  0   0
       0     0  0    0    0  0  0   0
       0     0  1/7  0    0  0  1/4 0
       0     0  1/2  0    0  1/4 0  0
       0     0   0   0    0  0   0  0
     ];
 C6=[  0     0  0    0  0 0 0 0
       0     0  0    0  0 0 0 0
       0     0  0    0   0 0 0 1/6
       0     0  0    0   0 0 0 0
       0     0  0    0  0 0 0 0
       0     0   0   0  0 0 0 0
       0     0   0    0  0 0 0 0
       0     0   1/6  0  0 0 0 0
     ];
 C=C1+C2+C3+C4+C5+C6;

 m=8;
 m4=[3 6 7];
 m5=[3 8];
 m2=[1 3 4];
 m2f1=[1 3];
 m2f2=[3 4];
 m3=[4 5];
 m1=[1 2 4];
  %%%%%%%%%%%%%%%%%%%
iteration=140;
 V_update=cell(1,iteration);
 vvupdate=cell(6,1);
 vvupdate{1,1}=v1_update{1,1};
 vvupdate{2,1}=v2_update{1,1};
 vvupdate{3,1}=v3_update{1,1};
 vvupdate{4,1}=v4_update{1,1};
 vvupdate{5,1}=v5_update{1,1};
 vvupdate{6,1}=v6_update{1,1};
 V_update{1,1}=vvupdate;
 
 %%DSA算法
  t_seq_syn=zeros(1,iteration)+1;
  tic;   %%%%开始计时
  k=1;
  ffseq=zeros(1,iteration+1);
   tic;
  
   t_seq_syn(1,1)=toc;   
   ffseq(1,1)=trace(C'*(v1'*v1));
 for k=1:iteration
     v1_update{1,k+1}=v1_update{1,k};
     v2_update{1,k+1}=v2_update{1,k};
     v3_update{1,k+1}=v3_update{1,k};
     v4_update{1,k+1}=v4_update{1,k};
     v5_update{1,k+1}=v5_update{1,k};
     v6_update{1,k+1}=v6_update{1,k};
     
     m42=setdiff(m4,m2);
     for i=1:length(m42)       %v4
         g4=zeros(m,1);
         for j=1:length(m4)
             g4=g4+C5(m42(i),m4(j))*v5_update{1,k}(:,m4(j));
         end
             g4=v5_update{1,k}(:,m42(i))-(1/(norm(C(m42(i),:),1))-0.001)*g4;   
         if(norm(g4)==0)
              v5_update{1,k+1}(:,m42(i))=zeros(num,1);
         else
             v5_update{1,k+1}(:,m42(i))= g4/norm(g4);
         end
     end 

     m52=setdiff(m5,m2);   %v5
      for i=1:length(m52)
         g5=zeros(m,1);
         for j=1:length(m5)
             g5=g5+C6(m52(i),m5(j))*v6_update{1,k}(:,m5(j));
         end
             g5=v6_update{1,k}(:,m52(i))-(1/(norm(C(m52(i),:),1))-0.001)*g5;
         if norm(g5)==0
              v6_update{1,k+1}(:,m52(i))=zeros(num,1);
         else
              v6_update{1,k+1}(:,m52(i))=g5/norm(g5);
         end
      end

      iim24=intersect(m2,m4);  %v2  共有部分为3，而且3同时在F1，F4里面，2和1的不同部分恰好也为3
      m21=setdiff(m2,m1);
      for i=1:length(m21)     
         g2=zeros(m,1);
         for j=1:length(m2f1)
             g2=g2+C1(m21(i),m2f1(j))*v1_update{1,k}(:,m2f1(j));
         end
         for j=1:length(m2f2)
             g2=g2+C4(m21(i),m2f2(j))*v4_update{1,k}(:,m2f2(j));
         end
         
          if ismember( m21(i),intersect(m2,m4))
                 for p=1:length(m4)
                     g2=g2+C5(m21(i),m4(p))*v5_update{1,k}(:,m4(p));   
                 end
          end
          if ismember( m21(i),intersect(m2,m5))
                 for p=1:length(m5)
                     g2=g2+C6(m21(i),m5(p))*v6_update{1,k}(:,m5(p));   
                 end
          end
          g2=v1_update{1,k}(:,m21(i))-(1/(norm(C(m21(i),:),1))-0.001)*g2;  
          if norm(g2)==0
              v1_update{1,k+1}(:,m21(i))=zeros(num,1);
              v4_update{1,k+1}(:,m21(i))=zeros(num,1);
          else
              v1_update{1,k+1}(:,m21(i))=g2/norm(g2);
              v4_update{1,k+1}(:,m21(i))=g2/norm(g2);
          end

      end
      for i=1:length(m2)
           if ismember(m2(i),intersect(m4,m2f1))
               v5_update{1,k+1}(:,m2(i))=v1_update{1,k+1}(:,m2(i));
           end
           if ismember(m2(i),intersect(m5,m2f1))
               v6_update{1,k+1}(:,m2(i))=v1_update{1,k+1}(:,m2(i));
           end
       end
   
         m31=setdiff(m3,m1);
         for i=1:length(m31)       %v3
             g3=zeros(m,1);
             for j=1:length(m3)
                 g3=g3+C3(m31(i),m3(j))*v3_update{1,k}(:,m3(j));
             end
                 g3=v3_update{1,k}(:,m31(i))-(1/(norm(C(m31(i),:),1))-0.001)*g3;  
             if(norm(g3)==0)
                  v3_update{1,k+1}(:,m31(i))=zeros(num,1);
             else
                 v3_update{1,k+1}(:,m31(i))= g3/norm(g3);
             end
         end 
      
       for i=1:length(m1)      %v1
         g1=zeros(m,1);
         for j=1:length(m1)
             g1=g1+C2(m1(i),m1(j))*v2_update{1,k}(:,m1(j));
         end
          if ismember( m1(i),intersect(m2f1,m1))
                 for p=1:length(m2f1)
                     g1=g1+C1(m1(i),m2f1(p))*v1_update{1,k}(:,m2f1(p));   
                 end
          end
          if ismember( m1(i),intersect(m2f2,m1))
                 for p=1:length(m2f2)
                     g1=g1+C4(m1(i),m2f2(p))*v4_update{1,k}(:,m2f2(p));  
                 end
          end
          if ismember( m1(i), intersect(m1,m3))
              for p=1:length(m3)
                   g1=g1+C3(m1(i),m3(p))*v3_update{1,k}(:,m3(p));  
              end
          end
          
          g1=v2_update{1,k}(:,m1(i))-(1/(norm(C(m1(i),:),1))-0.001)*g1;
          if norm(g1)==0
              v2_update{1,k+1}(:,m1(i))=zeros(num,1);
          else
              v2_update{1,k+1}(:,m1(i))=g1/norm(g1);
          end

      end
      
       for i=1:length(m1)
           if ismember(m1(i),intersect(m1,m3))
               v3_update{1,k+1}(:,m1(i))= v2_update{1,k+1}(:,m1(i));
           end
           if ismember(m1(i), intersect(m1,m2f1))
               v1_update{1,k+1}(:,m1(i))=v2_update{1,k+1}(:,m1(i));
           end
           if ismember(m1(i), intersect(m1,m2f2))
               v4_update{1,k+1}(:,m1(i))=v2_update{1,k+1}(:,m1(i));
           end
       end 
      vvupdate=cell(6,1);
     vvupdate{1,1}=v1_update{1,k+1};
     vvupdate{2,1}=v2_update{1,k+1};
     vvupdate{3,1}=v3_update{1,k+1};
     vvupdate{4,1}=v4_update{1,k+1};
     vvupdate{5,1}=v5_update{1,k+1};
     vvupdate{6,1}=v6_update{1,k+1};
     V_update{1,k+1}=vvupdate;
 t_seq_syn(1,k+1)=toc;           
 end
V1=cell(1,iteration);
V2=cell(1,iteration);
V3=cell(1,iteration);
V4=cell(1,iteration);
V5=cell(1,iteration);
V6=cell(1,iteration);
 for k=1:iteration            %%%%V_update{1,k+1}
       for q=1:length(m1)
              V2{1,k}(:,m1(q))=V_update{1,k}{2,1}(:,m1(q));
       end
          for q=1:length(m2f1)
              V1{1,k}(:,m2f1(q))=V_update{1,k}{1,1}(:,m2f1(q));
          end

          for q=1:length(m2f2)
              V4{1,k}(:,m2f2(q))=V_update{1,k}{4,1}(:,m2f2(q));
          end
          for q=1:length(m3)
              V3{1,k}(:,m3(q))=V_update{1,k}{3,1}(:,m3(q));
          end

          for q=1:length(m4)
              V5{1,k}(:,m4(q))=V_update{1,k}{5,1}(:,m4(q));
          end

          for q=1:length(m5)
              V6{1,k}(:,m5(q))=V_update{1,k}{6,1}(:,m5(q));
          end
 end
 Vseq=cell(1,iteration);
for k=1:iteration
  for i=1:length(m1)
          Vseq{1,k}(:,m1(i))= V2{1,k}(:,m1(i));
  end

  for i=1:length(m2f1)
          Vseq{1,k}(:,m2f1(i))= V1{1,k}(:,m2(i));
  end

  for i=1:length(m2f2)
          Vseq{1,k}(:,m2f2(i))= V4{1,k}(:,m2(i));
  end
  
  for i=1:length(m3)
          Vseq{1,k}(:,m3(i))= V3{1,k}(:,m3(i));
  end

 for i=1:length(m4)
          Vseq{1,k}(:,m4(i))= V5{1,k}(:,m4(i));
 end
  
 for i=1:length(m5)
          Vseq{1,k}(:,m5(i))= V6{1,k}(:,m5(i));
 end
end
XXseq=cell(1,iteration);
%%ffseq=zeros(1,iteration);
nablafseq=zeros(1,iteration);
xgap=zeros(1,iteration+1);
flim=-4.880952;
fgi=zeros(1,iteration+1);
for k=1:iteration
    XXseq{1,k}=Vseq{1,k}'*Vseq{1,k};
    ffseq(1,k+1)=trace(C'*XXseq{1,k});%%-flim;
    if k>2
        xgap(1,k-1)=norm(XXseq{1,k}-XXseq{1,k-1},'fro');
    end
   %%%gradf
    for i=1:num
        for j=1:num
            g_i=C(i,j)*Vseq{1,k}(:,j);
        end
        inng=(Vseq{1,k}(:,i)'*g_i);
        fgi(1,k)=fgi(1,k)+norm(g_i)^2-inng^2;
    end
end
figure(1);
plot(t_seq_syn,ffseq,'g-o','linewidth',1.5);
hold on;
figure(2);
plot(t_seq_syn,xgap,'g-o','linewidth',1.5);
hold on;
figure(3);
plot(t_seq_syn,fgi,'g-o','linewidth',1.5);
hold on;

%%%DAA算法
 iteration=800; 
 V_update=cell(1,iteration);
 vvupdate=cell(6,1);
 vvupdate{1,1}=v1_update{1,1};
 vvupdate{2,1}=v2_update{1,1};
 vvupdate{3,1}=v3_update{1,1};
 vvupdate{4,1}=v4_update{1,1};
 vvupdate{5,1}=v5_update{1,1};
 vvupdate{6,1}=v6_update{1,1};
 V_update{1,1}=vvupdate;
 a=[1,5,4,3,5,4,2];%%a=[2,1,3]  【45231】
 %a=[1,4,5,2,3];
  ffseq=zeros(1,iteration+1);
   tic;
   t_seq=zeros(iteration+1,1);
   t_seq(1,1)=toc;   
   ffseq(1,1)=trace(C'*(v1'*v1));
 for k=1:iteration
   knum=a(mod(k,7)+1); %%%knum为更新哪个个体，k为变量迭代次数  (k,5)
   V_update{1,k+1}=randombigupdate(knum,k,V_update);
 t_seq(k+1,1)=toc;        
 end
V1=cell(1,iteration);
V2=cell(1,iteration);
V3=cell(1,iteration);
V4=cell(1,iteration);
V5=cell(1,iteration);
V6=cell(1,iteration);
 for k=1:iteration            %%%%V_update{1,k+1}
       for q=1:length(m1)
              V2{1,k}(:,m1(q))=V_update{1,k}{2,1}(:,m1(q));
       end

          for q=1:length(m2f1)
              V1{1,k}(:,m2f1(q))=V_update{1,k}{1,1}(:,m2f1(q));
          end

          for q=1:length(m2f2)
              V4{1,k}(:,m2f2(q))=V_update{1,k}{4,1}(:,m2f2(q));
          end


          for q=1:length(m3)
              V3{1,k}(:,m3(q))=V_update{1,k}{3,1}(:,m3(q));
          end

          for q=1:length(m4)
              V5{1,k}(:,m4(q))=V_update{1,k}{5,1}(:,m4(q));
          end

          for q=1:length(m5)
              V6{1,k}(:,m5(q))=V_update{1,k}{6,1}(:,m5(q));
          end
 end
 Vseq=cell(1,iteration);
%  V=zeros(num,num);
for k=1:iteration
  for i=1:length(m1)
          Vseq{1,k}(:,m1(i))= V2{1,k}(:,m1(i));
  end

  for i=1:length(m2f1)
          Vseq{1,k}(:,m2f1(i))= V1{1,k}(:,m2(i));
  end

  for i=1:length(m2f2)
          Vseq{1,k}(:,m2f2(i))= V4{1,k}(:,m2(i));
  end
  
  for i=1:length(m3)
          Vseq{1,k}(:,m3(i))= V3{1,k}(:,m3(i));
  end

 for i=1:length(m4)
          Vseq{1,k}(:,m4(i))= V5{1,k}(:,m4(i));
 end
  
 for i=1:length(m5)
          Vseq{1,k}(:,m5(i))= V6{1,k}(:,m5(i));
 end
end
XXseq=cell(1,iteration);
nablafseq=zeros(1,iteration);
xgap=zeros(1,iteration+1);
flim=-4.880952;
fgi=zeros(1,iteration+1);
for k=1:iteration
    XXseq{1,k}=Vseq{1,k}'*Vseq{1,k};
    ffseq(1,k+1)=trace(C'*XXseq{1,k});
    if k>2
        xgap(1,k-1)=norm(XXseq{1,k}-XXseq{1,k-1},'fro');
    end
   %%%gradf
    for i=1:num
        for j=1:num
            g_i=C(i,j)*Vseq{1,k}(:,j);
        end
        inng=(Vseq{1,k}(:,i)'*g_i);
        fgi(1,k)=fgi(1,k)+norm(g_i)^2-inng^2;
    end
end
figure(1);
plot(t_seq,ffseq,'m-.','linewidth',2);
hold on;
figure(2);
plot(t_seq,xgap,'m-.','linewidth',2);
hold on; 
figure(3);
plot(t_seq,fgi,'m-.','linewidth',2);
hold on;

%%SDPLR算法
n=num;
bi=1;
Em=eye(n);
A=cell(n,1);
for i=1:n
    A{i,1}=Em(:,i)*Em(:,i)';
end
yv=rand(n,1);
sigma=0.3;%%%可调节
yita=8;%%>1
gamma=0.2;%%<1
iterations=15;
x0=v1;%rand(n,r);
tic;
tol_k=0; 
for j=1:n
    tol_k=tol_k+(trace(A{j,1}'*(x0*x0'))-bi)^2;
end
tseq=zeros(iterations+1,1);
tseq(1,1)=toc;
ffseq2=zeros(iterations+1,1);
ffseq2(1,1)=trace(C'*(v1*v1'));
xxseq=cell(1,iterations);
for i=1:iterations

    options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
    fun = @(x)f(x);
    [x,fval,exitflag,output,grad] = fminunc(fun,x0,options);
       tol=0; 
    for j=1:n
        tol=tol+(trace(A{j,1}'*(x*x'))-bi)^2;
    end
    if tol<gamma*tol_k
        for iny=1:n
            yv(iny,1)=yv(iny,1)-sigma*(trace(A{iny,1}'*(x*x'))-bi);
        end
        tol_k=tol;
    else
        sigma=sigma*yita;
        tol_k=tol;
    end
    xxseq{1,i}=x;
    x0=x;
    tseq(i+1,1)=toc;
end
t=toc
xgap2=zeros(1,iterations+1);
flim=-4.880952;
fgi2=zeros(1,iterations+1);

for i=1:iterations
    xxseq2{1,i}=xxseq{1,i}*xxseq{1,i}';
    ffseq2(i+1,1)=trace(C'*xxseq2{1,i});
	if i>2
        xgap2(1,i-1)=norm(xxseq2{1,i}-xxseq2{1,i-1},'fro');
    end
   %%%gradf
    for ii=1:num
        for j=1:num
            g_i2=C(ii,j)*Vseq{1,i}(:,j);
        end
        inng=(Vseq{1,i}(:,ii)'*g_i2);
        fgi2(1,i)=fgi2(1,i)+norm(g_i2)^2-inng^2;
    end
end
figure(1);
semilogy(tseq,ffseq2,'b-x','linewidth',1);
hold on;
flimseq=ones(iterations+1,1)*flim;
plot(tseq,flimseq,'color',[1 0.5 0],'linewidth',1,'linestyle','--');
hold on;
h=legend('{\rm DSA}','{\rm DAA}','{\rm SDPLR}');
set(h,'Interpreter','latex');
set(h,'FontName','Times New Roman','FontSize',12,'FontWeight','normal')
xlabel('Time (s)');
ylabel('$$log(f(k))$$','Interpreter','latex');
XX=x*x';
ff=trace(C'*(xxseq{1,iterations}*xxseq{1,iterations}'))
figure(2);
%%plot(t_seq,xgap,'m--','linewidth',2);
%hold on;
semilogy(tseq,xgap2,'b-x','linewidth',1);
hold on;
h=legend('{\rm DSA}','{\rm DAA}','{\rm SDPLR}');
set(h,'Interpreter','latex');
set(h,'FontName','Times New Roman','FontSize',12,'FontWeight','normal')
xlabel('Time (s)')
ylabel('$$log(\left\|X(k+1)-X(k)\right\|_F)$$','Interpreter','latex');
figure(3);
%plot(t_seq,fgi,'m-','linewidth',2);
%hold on;
semilogy(tseq,fgi2,'b-x','linewidth',1);
hold on
zeronabf=ones(iterations+1,1)*0;
semilogy(tseq,zeronabf,'color',[1 0.5 0],'linewidth',1,'linestyle','--');
h=legend('{\rm DSA}','{\rm DAA}','{\rm SDPLR}');
set(h,'Interpreter','latex');
set(h,'FontName','Times New Roman','FontSize',12,'FontWeight','normal')
 xlabel('Time (s)')
ylabel('$$log(\left\|{\rm grad}f(k)\right\|_F)$$','Interpreter','latex');

