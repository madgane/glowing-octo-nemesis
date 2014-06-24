
function cardinalityTest()

cvx_precision high
cvx_solver('sedumi')

cvx_begin

variables x(3,1) t

minimize(t)

subject to

Q = [1 -2 -3 ; -2 5 4 ; -3 4 17];

x' * Q * x + x(3,1) <= t;

norm(x,1) <= 2;

x <= 1;
x >= -1;

cvx_end

t = 0.01;
y = randn(3,1);
re_iterate = 1;
epsilon = zeros(3,1);
h_value = zeros(3,1);

while re_iterate
    
    for i = 1:3
        epsilon(i,1) = sub_grad(y(i,1),t);
        h_value(i,1) = h_func(y(i,1),t);
    end    
    
    cvx_begin
    
    variables x(3,1) m
    expression sum_g
    
    minimize(m)
    
    subject to
    
    Q = [1 -2 -3 ; -2 5 4 ; -3 4 17];
    
    x' * Q * x + x(3,1) <= m;
    
    sum_g = 0;
    for i = 1:3
        sum_g = sum_g + h_value(i,1) + epsilon(i,1) * (x(i,1) - y(i,1));
        (1/t) * abs(x(i,1)) - (1/t) * (h_value(i,1) + epsilon(i,1) * (x(i,1) - y(i,1))) <= 1;
    end
    
    (1/t) * (norm(x,1) - sum_g) <= 2;
    
    x <= 1;
    x >= -1;
    
    cvx_end
    
    if strfind(cvx_status,'Solved')
        y = x;
    else
        y = y / 2;
    end
    
end

end

function [out_val] = h_func(u,v)

out_val = max(0,(u - v)) + max(0,(-u - v));

end

function [sg] = sub_grad(u,v)

if u < -v
    sg = -1.0;
    return;
end

if u == -v
    sg = -0.5;
    return;
end

if u > -v
    if u < v
        sg = 0;
        return;
    end
end    

if u == v
    sg = 0.5;
    return;
end

if u > v
    sg = 1.0;
    return;
end

end

