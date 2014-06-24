function [sum_flops,prod_flops] = est_flops(mat_1,mat_2,varargin)

if mat_1(1,2) ~= mat_2(1,1)
    error('unequal matrices');
end

sum_flops = mat_1(1,1) * (mat_1(1,2) - 1) * mat_2(1,2);
prod_flops = mat_1(1,1) * mat_1(1,2) * mat_2(1,2);

if ~isempty(varargin)
    if varargin{1,1}(1,1)
        est_flops([mat_2(1,2),mat_1(1,1)],[mat_1(1,1),mat_2(1,2)]);
    end
end

end