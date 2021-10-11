% This is the source code of article:
% DANIEL: A Fast and Robust Consensus Maximization Method for Point Cloud Registration with High Outlier Ratios
% Copyright by Lei Sun (leisunjames@126.com) only for academic use.

clc;
clear all;
close all;

%% Set outlier ratio:
outlier_ratio=0.99; %0-0.99 in percentage, maximum: 0.99 (99% outliers)


%% Initiate Environments

n_ele=1000;
noise=0.01;
show_figure=0;
[pts_3d,pts_3d_,R_gt,t_gt]=Build_Scenario_(n_ele,noise,outlier_ratio,1,show_figure);


%% Daniel

tic;

itr_num=0;max_num=459;
R_th=0.133;t_th=0.07;
c_max_final=0;c_max_set_final=[];
min_size=0.01*n_ele;
while 1
    
    itr_num=itr_num+1;
    
    ele_this=randperm(n_ele,1);
    
    comp_set=zeros(1,1);count=0;
    for j=1:n_ele
        res=abs(norm(pts_3d(ele_this,:)-pts_3d(j,:))-norm(pts_3d_(ele_this,:)-pts_3d_(j,:)));
        if res<=5*noise && res>0
            count=count+1;
            comp_set(count)=j;
        end
    end
    
c_max=0;c_max_set=[];in_samp=0;in_samp_=0;

% inner sampling:
if count>=min_size
    
in_num=2*round(log(1-0.995)/log(1-(min_size/count)^2));%min([1055,2*round(log(1-0.995)/log(1-(min_size/count)^2))]);

store1=zeros(9,in_num);store2=zeros(3,in_num);store_num=zeros(in_num,3);
% for in_samp_=1:in_num
while 1
in_samp_=in_samp_+1;
% for in_samp1=1:count-1
% for in_samp2=in_samp1+1:count
% ele_in_this=comp_set([in_samp1,in_samp2]);

   while 1 
   ele_in_this=comp_set(randperm(count,2)); 
   if ele_in_this(1)~=ele_in_this(2) && ele_in_this(1)~=ele_this && ele_in_this(2)~=ele_this
      break 
   end
   end
   
if abs(norm(pts_3d(ele_in_this(1),:)-pts_3d(ele_in_this(2),:))-norm(pts_3d_(ele_in_this(1),:)-pts_3d_(ele_in_this(2),:)))<=5*noise
in_samp=in_samp+1;

   ele_in_this(3)=ele_this;
   [R_raw,t_raw]=Horn_minimal(pts_3d(ele_in_this,:),pts_3d_(ele_in_this,:));
   
   store1(:,in_samp)=R_raw(:);
    
   store2(:,in_samp)=t_raw;

   store_num(in_samp,:)=ele_in_this;
   
cc=0;ok_set=zeros(1,3);
if in_samp>=2
       for i_=1:in_samp-1
           
            t1=abs(store2(1,in_samp)-store2(1,i_));
            if t1<=t_th
            t2=abs(store2(2,in_samp)-store2(2,i_));
            if t2<=t_th
            t3=abs(store2(3,in_samp)-store2(3,i_));
            if t3<=t_th          
            
            r1=abs(store1(1,in_samp)-store1(1,i_));
            if r1<=R_th
            r2=abs(store1(2,in_samp)-store1(2,i_));
            if r2<=R_th
            r3=abs(store1(3,in_samp)-store1(3,i_));
            if r3<=R_th
            r4=abs(store1(4,in_samp)-store1(4,i_));
            if r4<=R_th
            r5=abs(store1(5,in_samp)-store1(5,i_));
            if r5<=R_th
            r6=abs(store1(6,in_samp)-store1(6,i_));
            if r6<=R_th
            r7=abs(store1(7,in_samp)-store1(7,i_));
            if r7<=R_th
            r8=abs(store1(8,in_samp)-store1(8,i_));
            if r8<=R_th
            r9=abs(store1(9,in_samp)-store1(9,i_));
            if r9<=R_th

            if length(unique([store_num(i_,:),ele_in_this]))>=5
            
            cc=cc+1;
            ok_set(cc,:)=store_num(i_,:);
            
if cc>=1
       
%        Ok_set=unique([ok_set(:);ele_in_this']');
%        c_this=length(Ok_set);
           
R_1=reshape(store1(:,i_),3,[]);R_2=reshape(store1(:,in_samp),3,[]);
R_raw=R_1*expm(0.5*logm(R_1'*R_2));
            
t_raw=(store2(:,in_samp)+store2(:,i_))/2;

con_set=zeros(1,1);

res=zeros(1,count);coun=0;
    
    for i=1:count
        res(i)=norm(1*R_raw*((pts_3d(comp_set(i),:)))'+t_raw-(pts_3d_(comp_set(i),:))');
        if res(i)<=6*noise 
            coun=coun+1;
            con_set(coun)=comp_set(i);
        end
    end
    
       if coun>c_max
           c_max=coun;
           c_max_set=con_set;
           
           in_num=log(1-0.99)/log(1-(c_max/count)^2);
       end
end
                
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end
       end
            
   
end
       
       if in_samp_>=in_num
           break
       end
end

% end
end

end

if c_max>c_max_final && c_max>3
    
opt_set=c_max_set;

q_=zeros(3,1);
p_=zeros(3,1);

len_o=length(opt_set);

for i=1:len_o

q_(1)=q_(1)+pts_3d_(opt_set(i),1);
q_(2)=q_(2)+pts_3d_(opt_set(i),2);
q_(3)=q_(3)+pts_3d_(opt_set(i),3);

p_(1)=p_(1)+pts_3d(opt_set(i),1);
p_(2)=p_(2)+pts_3d(opt_set(i),2);
p_(3)=p_(3)+pts_3d(opt_set(i),3);

end

p_=p_/len_o;
q_=q_/len_o;

s_best=1;

H=zeros(3,3);
for i=1:len_o
    H=H+(pts_3d(opt_set(i),:)'-p_)*(pts_3d_(opt_set(i),:)'-q_)';
end

[U,~,V]=svd(H);

R_opt=V*U';

t_opt=q_-s_best*R_opt*p_;


con_set_=zeros(1,1);

res=zeros(1,n_ele);coun=0;
    
    for i=1:n_ele
        res(i)=(s_best*R_opt*((pts_3d(i,:)))'+t_opt-(pts_3d_(i,:))')'*(s_best*R_opt*((pts_3d(i,:)))'+t_opt-(pts_3d_(i,:))');
        if sqrt(res(i))<=6*noise 
            coun=coun+1;
            con_set_(coun)=i;
        end
    end
    
    c_max_final=coun;
    c_max_set_final=con_set_;
    max_num=log(1-0.99)/log(1-(c_max_final/n_ele));
    
end
        
    if itr_num>=max_num 
        break
    end
    
end

best_set=unique(c_max_set_final(:)');

opt_set=best_set;

q_=zeros(3,1);
p_=zeros(3,1);

opt_set=unique(opt_set);

len_o=length(opt_set);

for i=1:len_o

q_(1)=q_(1)+pts_3d_(opt_set(i),1);
q_(2)=q_(2)+pts_3d_(opt_set(i),2);
q_(3)=q_(3)+pts_3d_(opt_set(i),3);

p_(1)=p_(1)+pts_3d(opt_set(i),1);
p_(2)=p_(2)+pts_3d(opt_set(i),2);
p_(3)=p_(3)+pts_3d(opt_set(i),3);

end

p_=p_/len_o;
q_=q_/len_o;

s_best=1;

H=zeros(3,3);
for i=1:len_o
    H=H+(pts_3d(opt_set(i),:)'-p_)*(pts_3d_(opt_set(i),:)'-q_)';
end

[U,~,V]=svd(H);

R_opt=V*U';

t_opt=q_-s_best*R_opt*p_;


best_set=zeros(1,1);

res=zeros(1,n_ele);s_set=zeros(1,1);coun=0;
    
    for i=1:n_ele
        res(i)=(s_best*R_opt*((pts_3d(i,:)))'+t_opt-(pts_3d_(i,:))')'*(s_best*R_opt*((pts_3d(i,:)))'+t_opt-(pts_3d_(i,:))');
        if sqrt(res(i))<=6*noise 
            coun=coun+1;
            best_set(coun)=i;
        end
    end


q_=zeros(3,1);
p_=zeros(3,1);

best_size=length(best_set);

for i=1:best_size

q_(1)=q_(1)+pts_3d_(best_set(i),1);
q_(2)=q_(2)+pts_3d_(best_set(i),2);
q_(3)=q_(3)+pts_3d_(best_set(i),3);

p_(1)=p_(1)+pts_3d(best_set(i),1);
p_(2)=p_(2)+pts_3d(best_set(i),2);
p_(3)=p_(3)+pts_3d(best_set(i),3);

end

p_=p_/best_size;
q_=q_/best_size;

s_best=1;

H=zeros(3,3);
for i=1:best_size
    H=H+(pts_3d(best_set(i),:)'-p_)*(pts_3d_(best_set(i),:)'-q_)';
end

[U,~,V]=svd(H);

R_opt=V*U';

t_opt=q_-s_best*R_opt*p_;

time=toc();

R_error=getAngularError(R_gt,R_opt)*180/pi;

t_error=norm(t_opt - t_gt');


%% Display Results: 

disp(['Rotation Error (in degrees): ', num2str(R_error)]);

disp(['Translation Error : ', num2str(t_error)]);

disp(['Runtime of RANSIC (in seconds): ', num2str(time)]);



%% Function Horn's :

function [R,t]=Horn_minimal(pts_3d,pts_3d_)

    v12=pts_3d(2,:)-pts_3d(1,:);
    X_axis=v12'/norm(v12);
    v13=pts_3d(3,:)-pts_3d(1,:);
    v23=cross(v12,v13);
    Y_axis=v23'/norm(v23);
    Z_axis=cross(X_axis,Y_axis);
    
    v12=pts_3d_(2,:)-pts_3d_(1,:);
    X_axis_=v12'/norm(v12);
    v13=pts_3d_(3,:)-pts_3d_(1,:);
    v23=cross(v12,v13);
    Y_axis_=v23'/norm(v23);
    Z_axis_=cross(X_axis_,Y_axis_);

    R=[X_axis_,Y_axis_,Z_axis_]*[X_axis,Y_axis,Z_axis]';
    
q_=zeros(3,1);
p_=zeros(3,1);
    
for i=1:3

q_(1)=q_(1)+pts_3d_(i,1);
q_(2)=q_(2)+pts_3d_(i,2);
q_(3)=q_(3)+pts_3d_(i,3);

p_(1)=p_(1)+pts_3d(i,1);
p_(2)=p_(2)+pts_3d(i,2);
p_(3)=p_(3)+pts_3d(i,3);

end

p_=p_/3;
q_=q_/3;

t=q_-R*p_;

end
