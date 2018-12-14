%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%作者：鲁尚宗 时间：2018年12月14日 
%本程序实现了通过无人机位置信息计算无人机的瞬时速度
%并且比较了不同求导的步长、不同半径、不同插值方法对结果的影响
%还实现了对数据的预处理，即删除不可用数据的功能
%不同的情况通过不同的注释组合来实现
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %清空内存

%步长相当于时间间隔
step=1;
%step=0.5;
%step=0.1;
%step=0.01;

%对原始数据进行插值，作为结果的一个参考值
method='spline';  %样条插值
%method='linear';  %线性插值


%地球半径，用来计算距离
%r=6371e3;  %平均值
r=6378e3;   %赤道半径

%读取数据，将每一种数据读入到数组中
filename='POS2_data.txt';
a=importdata(filename);
time=a.data(:,1); %时间
latitude=a.data(:,2); %纬度
longitude=a.data(:,3); %经度
height=a.data(:,4); %高度
v_east=a.data(:,6); %东速度
v_north=a.data(:,7); %北速度
v_up=a.data(:,8);  %上升速度

%检查数据是否齐全，如果有某个数据为空，则记录下这一行的编号，最后所有数据删除此行
num=size(time);
k=zeros(100,1);  %预计空值的行不会大于100个
j=1;
%如果遇到空值就将这一行的编号存下来
for i=1:num
    if isnan(time(i))||isnan(latitude(i))||isnan(longitude(i))||isnan(height(i))...
            ||isnan(v_east(i))||isnan(v_north(i))||isnan(v_up(i))
        k(j)=i;
        j=j+1;
    end
end
%将所有数据的空值行删去
time(k(1:j-1),:)=[];
latitude(k(1:j-1),:)=[];
longitude(k(1:j-1),:)=[];
height(k(1:j-1),:)=[];
v_east(k(1:j-1),:)=[];
v_north(k(1:j-1),:)=[];
v_up(k(1:j-1),:)=[];

%再次获取数据大小
num=size(time);

%对时间进行扩充，为等差序列
time_inter=(time(1):step:time(num(1)))';

%和速度大小为三个分量的平方和再开根
v=(v_east.^2+v_north.^2+v_up.^2).^0.5;

%对原始数据进行插值，一方面解决数据间隔不相等的问题，另一方面对计算结果做个比较
v_inter=interp1(time,v,time_inter,method);
lat_inter=interp1(time,latitude,time_inter,method);
lon_inter=interp1(time,longitude,time_inter,method);
height_inter=interp1(time,height,time_inter,method);

%获取插值后的数据大小
number=size(lat_inter);
num=number(1);

%创建一个地球半径数组，用于计算路程
%每一个值为地球半径加上对应高度
r_list=height_inter;
for i=1:num
    r_list(i)=r+height_inter(i);
end

%计算路程
distance=r_list;
distance(1)=0; %第一个路程为0


for i=2:num
    %分别计算经度、维度、高度的距离，实际路程为三者的平方和再开根
    %需要进行角度和弧度的转换
    dis_lat=(lat_inter(i)-lat_inter(i-1))/360*2*pi*r_list(i);
    dis_lon=(lon_inter(i)-lon_inter(i-1))/360*2*pi*r_list(i)*cos(lat_inter(i)*pi/180);
    dis_height=height_inter(i)-height_inter(i-1);
    distance(i)=distance(i-1)+(dis_lat^2+dis_lon^2+dis_height^2)^0.5;
end

%采用两点公式、三点公式、五点公式、样条求导的方法计算速度
v_two=two_point(distance,step);
v_three=three_point(distance,step);
v_five=five_point(distance,step);
v_spline=spline(distance,step,v_inter(1),v_inter(num));

%用插值后的速度减去计算出来的结果，进行比较
error_two=v_inter-v_two;
error_three=v_inter-v_three;
error_five=v_inter-v_five;
error_spline=v_inter-v_spline;

%绘制图形
%第一张图，绘制计算结果
figure(1); 
subplot(1,2,1);
plot(time_inter,v_inter);
hold on;
plot(time_inter,v_two);
hold on;
plot(time_inter,v_three);
hold on;
plot(time_inter,v_five);
hold on;
plot(time_inter,v_spline);
hold on;
scatter(time,v);
title('计算结果')
xlabel('时间');
ylabel('速度');
legend('插值结果','两点公式','三点公式','五点公式','样条求导');

%第二张图，绘制误差
subplot(1,2,2);
plot(time_inter,error_two);
hold on;
plot(time_inter,error_three);
hold on;
plot(time_inter,error_five);
hold on;
plot(time_inter,error_spline);
title('误差比较')
xlabel('时间');
ylabel('误差');
legend('两点公式','三点公式','五点公式','样条求导');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%两点公式，计算速度
%传入路程的数组，返回相同大小的速度数组
%路程的数组应该是递增的数组
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result=two_point(list,step)
%获取数组大小
num=size(list);
number=num(1);

result=list;

%除了最后一个速度，其它的速度都是后面的减去前面的路程再除以时间间隔
for i=1:number-1
    result(i)= (list(i+1)-list(i))/step;
end
%最后一个速度等于前面的速度，这是两点法的缺陷
result(number)=result(number-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%三点公式，计算速度
%传入路程的数组，返回相同大小的速度数组
%路程的数组应该是递增的数组
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result=three_point(list,step)
%获取数组大小
num=size(list);
number=num(1);

result=list;

%除了最后一个和第一个其它的都要用相邻两个点的数据来算
for i=2:number-1
    result(i)=(list(i+1)-list(i-1))/(2*step);
end

%第一个点和最后一个点要分开计算
result(1)=(list(1)*(-3)+list(2)*4-list(3))/(2*step);
result(number)=(list(number-2)-4*list(number-1)+3*list(number))/(2*step);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%五点公式，计算速度
%传入路程的数组，返回相同大小的速度数组
%路程的数组应该是递增的数组
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result=five_point(list,step)

%获取数组大小
num=size(list);
number=num(1);

result=list;

%除了最后两个和开始两个其它的都要用相邻四个个点的数据来算
for i=3:number-2
    result(i)=(list(i-2)-8*list(i-1)+8*list(i+1)-list(i+2))/(12*step);
end

%最后两个和开始两个单独计算
result(1)=(-25*list(1)+48*list(2)-36*list(3)+16*list(4)-3*list(5))/(12*step);
result(2)=(-3*list(1)-10*list(2)+18*list(3)-6*list(4)+list(5))/(12*step);
result(number-1)=(-list(number-4)+6*list(number-3)-...
    18*list(number-2)+10*list(number-1)+3*list(number))/(12*step);
result(number)=(3*list(number-4)-16*list(number-3)+...
    36*list(number-2)-48*list(number-1)+25*list(number))/(12*step);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%样条求导，计算速度
%传入路程的数组和边界条件，返回速度数组
%路程的数组应该是递增的数组
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result=spline(list,step,v0,vn)

%获取数组大小
num=size(list);
number=num(1);

%系数矩阵
A=zeros(number-2,number-2);

%创建系数矩阵，第一行和最后一行单独赋值，其余可用循环来做
A(1,1)=4;
A(1,2)=1;
A(number-2,number-3)=1;
A(number-2,number-2)=4;
for i=2:number-3
    A(i,i-1)=1;
    A(i,i)=4;
    A(i,i+1)=1;
end

%创建方程组中最右边的值
g=zeros(number-2,1);
%首尾两个值单独赋值，其余可用循环
g(1)=3*(list(3)-list(1))/step-v0;
g(number-2)=3*(list(number)-list(number-2))/step-vn;
for i=2:number-3
    g(i)=3*(list(i+1)-list(i-1))/step;
end

%解线性方程组，采用对称矩阵的方法
opts.SYM = true;
tem=linsolve(A,g,opts);
result=list;
result(1)=v0;
result(number)=vn;
result(2:number-1)=tem;
end

