%6对称光照方向，方向不再平行，非平面-斜面
clear; close all; clc;
addpath(genpath(pwd));  % 添加当前路径下的所有子目录

g_pic_num = 6;
shadowThresh = 0.01;

%% 光照方向
L = zeros(3, g_pic_num);
Slant = 45;
Slant_sin = sin(Slant/180*pi);
Slant_cos = cos(Slant/180*pi);
for i = 1:g_pic_num
    Tilt = (i-1) * 360/g_pic_num;
    %%本实验的坐标系
    L(1,i) = Slant_sin * cos(Tilt/180*pi);   
    L(2,i) = Slant_sin * sin(Tilt/180*pi);
    L(3,i) = Slant_cos;
end

%% 图像仿真,仿真一个平面，法向（0，0，1）
w1 = 600; h1 = 600; %相机分辨率
[X,Y]= meshgrid(1:w1, 1:h1);

scale = 3; %6-50,3-100,2-150,1.5-200 %观测的真实尺度大小
x = (X - 300)/scale; %+-150 length(长)
y = (Y - 300)/scale; %+-150 width(宽)

% 斜面
z_a = 0; z_b = -0.3;
z = - z_a*x - z_b*y;
figure;mesh(z);title('z');

% 标准值，真实尺度
z_gt = z;
z_gt_pixel = z * scale; % 像素尺度，两个尺度不混淆
N_gt = [z_a, z_b, 1];
N_temp = sqrt(z_a.^2 + z_b.^2 + 1);
N_gt = N_gt/N_temp;

r = 400; h = 400;
% 光源到物体的距离变化 % 真实尺度
d1 = sqrt( (x - r).^2 + y.^2 + (h-z).^2 );
d2 = sqrt( (x - r/2).^2 + (y-sqrt(3)*r/2).^2 + (h-z).^2 );
d3 = sqrt( (x + r/2).^2 + (y-sqrt(3)*r/2).^2 + (h-z).^2 );
d4 = sqrt( (x + r).^2 + y.^2 + (h-z).^2 );
d5 = sqrt( (x + r/2).^2 + (y+sqrt(3)*r/2).^2 + (h-z).^2 );
d6 = sqrt( (x - r/2).^2 + (y+sqrt(3)*r/2).^2 + (h-z).^2 );

I0 = 10000; para_attenuation = 3;
% 斜面的仿真光强 N*L/|L|,由于|N|相同未计算 % 真实尺度
I1 = I0*(z_a .*(r-x) + z_b.*(-y) + (h-z))./(d1.^para_attenuation);
I2 = I0*(z_a .*(r/2-x) + z_b.*(sqrt(3)*r/2-y) + (h-z))./(d2.^para_attenuation);
I3 = I0*(z_a .*(-r/2-x) + z_b.*(sqrt(3)*r/2-y) + (h-z))./(d3.^para_attenuation);
I4 = I0*(z_a .*(-r-x) + z_b.*(-y) + (h-z))./(d4.^para_attenuation);
I5 = I0*(z_a .*(-r/2-x) + z_b.*(-sqrt(3)*r/2-y) + (h-z))./(d5.^para_attenuation);
I6 = I0*(z_a .*(r/2-x) + z_b.*(-sqrt(3)*r/2-y) + (h-z))./(d6.^para_attenuation);

% 显示，调整到[0 255]
max_I = max([prctile(I1(:), 99),prctile(I2(:), 99),prctile(I3(:), 99),prctile(I4(:), 99),prctile(I5(:), 99),prctile(I6(:), 99)]);
I1 = I1/max_I * 255;
I2 = I2/max_I * 255;
I3 = I3/max_I * 255;
I4 = I4/max_I * 255;
I5 = I5/max_I * 255;
I6 = I6/max_I * 255;
figure; subplot(2,3,1); imshow(uint8(I1)); subplot(2,3,2); imshow(uint8(I2)); subplot(2,3,3); imshow(uint8(I3));
subplot(2,3,4); imshow(uint8(I4)); subplot(2,3,5); imshow(uint8(I5)); subplot(2,3,6); imshow(uint8(I6));

%% 读入图片
[M, N, C] = size(I1);
I = ones(M, N, g_pic_num);
I(:,:,1) = I1;
I(:,:,2) = I2;
I(:,:,3) = I3;
I(:,:,4) = I4;
I(:,:,5) = I5;
I(:,:,6) = I6;

%% 图像的尺度信息
[g_rows, g_cols] = size(I1);
g_length = g_rows * g_cols;

% Create a shadow mask.
shadow_mask = (I > shadowThresh);
se = strel('disk', 2);
for i = 1:g_pic_num
  % Erode the shadow map to handle de-mosaiking artifact near shadow boundary.
  shadow_mask(:,:,i) = imerode(shadow_mask(:,:,i), se);
end

[rho, n] = PhotometricStereo(I, shadow_mask, L); %法向没有尺度

%% Visualize the normal map. axis xy;
N_RGB(:,:,1)  = (n(:,:,1) + 1) / 2;
N_RGB(:,:,2)  = (n(:,:,2) + 1) / 2;
N_RGB(:,:,3)  = n(:,:,3);
figure; imshow(N_RGB); 
% figure; imshow(rho);

%% Estimate depth map from the normal vectors.
fprintf('Estimating depth map from normal vectors...\n');
p = -n(:,:,1) ./ n(:,:,3);
q = -n(:,:,2) ./ n(:,:,3);
p(isnan(p)) = 0; %判断数组的元素是否是NaN。 NaN 即 Not a Number 的缩写。
q(isnan(q)) = 0;

figure; subplot(1,2,1); mesh(p); subplot(1,2,2); mesh(q);

%% integration
Height_poisson =poisson_solver_function_neumann(p, q);  % 像素尺度
% Height_poisson = flipud(Height_poisson);
% figure;mesh(Height_poisson); title('积分高度(像素)')

%% The propsed method to fit the bias. (拟合 y - namda * H  = a x.^2 + b y.^2 + c x + d y + f;)
index = floor(randi(length(Height_poisson(:)),1000,1));
Height_diff = Height_poisson(index) - z_gt_pixel(index); % 像素尺度
[y, x] = ind2sub(size(Height_poisson), index); %ind2sub得到的是行列，注意 行列 与 xy 相反
AA = [x.^2, y.^2, x, y, ones(length(x),1)];

% 方法1：矩阵除法
para = AA \ Height_diff; % 像素尺度

% Compute the corrected height (根据参数计算真实结果)
[m, n] = size(Height_poisson); %size 得到的是行列
[xx, yy] = meshgrid(1:n, 1:m);          %注意行、列，和x，y的对应关系 % 注意像素尺度，不能混淆！！
fitted_height = [xx(:).^2, yy(:).^2, xx(:), yy(:), ones(m*n,1)] * para(:);
Height_correct = Height_poisson(:) - fitted_height;
Height_correct = reshape(Height_correct, m, n);
fitted_height = reshape(fitted_height, m, n);


figure; subplot(1,3,1); mesh(Height_poisson); title('Initial Height'); xlabel('x(pixel)'); ylabel('y(pixel)'); zlabel('Initial Height(pixel)'); 
subplot(1,3,2);  mesh(fitted_height); title('Fitted Bias Height'); xlabel('x(pixel)'); ylabel('y(pixel)'); zlabel('Fitted Bias Height(pixel)');
subplot(1,3,3); mesh(Height_correct); title('Correct Height'); xlabel('x(pixel)'); ylabel('y(pixel)'); zlabel('Correct Height(pixel)');
set(gcf,'unit','centimeters','position',[15 15 19 6]);

% Height Error Analysis
Height_d = Height_correct - z_gt_pixel; % 像素尺度
disp(['平均高度误差：' num2str(mean(mean(abs(Height_d)))) '像素']);

%figure; mesh(Height_d/scale); title('Height error(mm)');colorbar;
disp(['平均高度误差：' num2str(mean(mean(abs(Height_d/scale)))) '毫米']);

% Angle error Analysis (偏差角度的误差分析)
% 法向没有尺度，但是要确认光度立体高度积分的尺度都是像素，或都是真实尺度；混淆后计算的法向是错误的！！！
[p_correct,q_correct] = gradient(Height_correct);
bb = sqrt(p_correct.^2 + q_correct.^2 + 1);
N_co(:,:,1) = - p_correct./bb;
N_co(:,:,2) = - q_correct./bb;
N_co(:,:,3) = 1./bb;
N_angle = N_co(:,:,1).*N_gt(1) + N_co(:,:,2).*N_gt(2) + N_co(:,:,3).*N_gt(3);
Height_n_error = acos(N_angle)/pi * 180;

figure; subplot(1,2,1); mesh(Height_d); title('Height error'); xlabel('x(pixel)'); ylabel('y(pixel)'); zlabel('Height error(pixel)');
subplot(1,2,2);mesh(Height_n_error); title('Angle error');  xlabel('x(pixel)'); ylabel('y(pixel)'); zlabel('Angle error(degree)');
set(gcf,'unit','centimeters','position',[10 5 13 6]);

disp(['平均角度误差：' num2str(mean(mean(abs(Height_n_error)))) '度']);