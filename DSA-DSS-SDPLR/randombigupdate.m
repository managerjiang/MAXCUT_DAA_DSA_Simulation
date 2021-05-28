 function vvupdate=randombigupdate(knum,k,v0_update) 
 global C1;
 global C2;
 global C3;
 global C4;
 global C5;
 global C6;
 global C;
 global m;
 global m1;
 global m2;
 global m3;
 global m4;
 global m5;
 global m2f1;
 global m2f2;
 global theta;
 global iteration;
v1_update=cell(1,iteration);
v2_update=cell(1,iteration);
v3_update=cell(1,iteration);
v4_update=cell(1,iteration);
v5_update=cell(1,iteration);
v6_update=cell(1,iteration);
num=8;
  step=zeros(8,1);
  vvupdate=cell(6,1);
  
     v1_update{1,k}=v0_update{1,k}{1,1};
     v2_update{1,k}=v0_update{1,k}{2,1};
     v3_update{1,k}=v0_update{1,k}{3,1};
     v4_update{1,k}=v0_update{1,k}{4,1};
     v5_update{1,k}=v0_update{1,k}{5,1};
     v6_update{1,k}=v0_update{1,k}{6,1};

     v1_update{1,k+1}=v1_update{1,k};
     v2_update{1,k+1}=v2_update{1,k};
     v3_update{1,k+1}=v3_update{1,k};
     v4_update{1,k+1}=v4_update{1,k};
     v5_update{1,k+1}=v5_update{1,k};
     v6_update{1,k+1}=v6_update{1,k};
     if knum==4
         step=zeros(8,1);
         for i=1:length(m4)
              step(i)=norm(C5(m4(i),:));
         end
         ss=max(step);
         theta=1/ss;
         m42=setdiff(m4,m2);
         for i=1:length(m42)       %v4
             g4=zeros(m,1);
             for j=1:length(m4)
                 g4=g4+C5(m42(i),m4(j))*v5_update{1,k}(:,m4(j));
             end
          
                 g4=v5_update{1,k}(:,m42(i))-(1/(norm(C(m42(i),:),1))-0.001)*g4;   %%%%%%%%%%%%%%%
             if(norm(g4)==0)
                  v5_update{1,k+1}(:,m42(i))=zeros(num,1);
             else
                 v5_update{1,k+1}(:,m42(i))= g4/norm(g4);
             end
         end 
     end
     if knum==5
         step=zeros(8,1);
          for i=1:length(m5)
              step(i)=norm(C6(m5(i),:));
         end
         ss=max(step);
         theta=1/ss;
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
     end

     if knum==2
         step=zeros(8,1);
          for i=1:length(m2f1)
              step(i)=norm(C1(m2f1(i),:),1);
          end
         for i=1:length(m2f2)
             step(length(m2f1)+i)=norm(C4(m2f2(i),:));
         end
         ss=max(step);
         theta=1/ss;
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
                         g2=g2+C5(m21(i),m4(p))*v5_update{1,k}(:,m4(p));   %%%%%
                     end
              end
              if ismember( m21(i),intersect(m2,m5))
                     for p=1:length(m5)
                         g2=g2+C6(m21(i),m5(p))*v6_update{1,k}(:,m5(p));   %%%%
                     end
              end
              g2=v1_update{1,k}(:,m21(i))-(1/(norm(C(m21(i),:),1))-0.001)*g2;  %%%只有第一次时候，v1和v4的初始值有可能不同
              if norm(g2)==0
                  v1_update{1,k+1}(:,m21(i))=zeros(num,1);
                  v4_update{1,k+1}(:,m21(i))=zeros(num,1);   %%%%%这里面v1_update和v4_update只是为了写程序方便
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
     end
     
     if knum==3
         step=zeros(8,1);
         for i=1:length(m3)
              step(i)=norm(C3(m3(i),:));
         end
         ss=max(step);
         theta=1/ss;
         m31=setdiff(m3,m1);
         for i=1:length(m31)       %v3
             g3=zeros(m,1);
             for j=1:length(m3)
                 g3=g3+C3(m31(i),m3(j))*v3_update{1,k}(:,m3(j));
             end
                 g3=v3_update{1,k}(:,m31(i))-(1/(norm(C(m31(i),:),1))-0.001)*g3;   %%%%%%%%%%%%%%%
             if(norm(g3)==0)
                  v3_update{1,k+1}(:,m31(i))=zeros(num,1);
             else
                 v3_update{1,k+1}(:,m31(i))= g3/norm(g3);
             end
         end 
     end
     if knum==1
         step=zeros(8,1);
         for i=1:length(m3)
              step(i)=norm(C3(m3(i),:));
         end
         ss=max(step);
         theta=1/ss;
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
     end
     vvupdate{1,1}=v1_update{1,k+1};
     vvupdate{2,1}=v2_update{1,k+1};
     vvupdate{3,1}=v3_update{1,k+1};
     vvupdate{4,1}=v4_update{1,k+1};
     vvupdate{5,1}=v5_update{1,k+1};
     vvupdate{6,1}=v6_update{1,k+1};
 end