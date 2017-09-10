% Check if certain points are on the same plain
function [points] = generate_points(num_points,epsilon)


%the point on the plane
x = zeros(7,1);

%choosing a random pont x(s,t) on the plane
s = rand(1)/4;
t = rand(1)/4;

x(1) = s;
x(2) = 1/4 - t;
x(3) = t;
x(4) = 1/4 - s;
x(5) = s;
x(6) = 1/4 - t;
x(7) = t;

x_s_plus_1 = x;
x_s_plus_1(1) = x_s_plus_1(1) + 1;
x_s_plus_1(4) = x_s_plus_1(4) - 1;
x_s_plus_1(5) = x_s_plus_1(5) + 1;
x_t_plus_1 = x;
x_t_plus_1(2) = x_t_plus_1(2) - 1;
x_t_plus_1(3) = x_t_plus_1(3) + 1;
x_t_plus_1(6) = x_t_plus_1(6) - 1;
x_t_plus_1(7) = x_t_plus_1(7) + 1;

%the vectors with direction of the plane. These vectors dtart from (0,...,0)
x_1 = x_s_plus_1 - x;
x_2 = x_t_plus_1 - x;

A = zeros(5,7); %the equetions that define the normal n
b = zeros(7,1); %the solutions of the equations

for i=1:7
    A(1,i) = x_1(i);
    A(2,i) = x_2(i);
end
A(3,2) = 1;
A(3,3) = 1;
A(4,4) = 1;
A(4,5) = 1;
A(5,6) = 1;
A(5,7) = 1;
%we choose n(1) and n(2) to be 1.
A(6,1) = 1; 
A(7,2) = 1;

b(3) = 1/4;
b(4) = 1/4;
b(5) = 1/4;
b(6) = 1;
b(7) = 1;

%the normal
n = A\b;

%y is the point that its distance from x on the plane is epsilon, and the
%direction of the distance is the normal n. y is on the line x+n*t. We look
%for the right t s.t. ||y - x|| = epsilon.
y = zeros(7,1);
%epsilon = the distance between x and y. 

n2 = n.^2;

sum_of_n2 = 0;
for i=1:7
    sum_of_n2 = sum_of_n2 + n2(i);
end

sum_of_n2 = sqrt(sum_of_n2);

t = epsilon/sum_of_n2;

y = x + n*t;


y_vec = zeros(7,num_points);

% Randomize points
for i = 1:num_points 
    y_vec(i,:) = x+n*t*i/points;
end

points = y_vec;