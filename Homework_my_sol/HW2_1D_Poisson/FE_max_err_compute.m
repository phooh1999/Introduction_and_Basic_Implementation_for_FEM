function error = FE_max_err_compute(ref_fun,Pb_trial,solution)
% 

num = size(Pb_trial,2);

err = sparse(1,num);

for i = 2:num-1
    
    ref = ref_fun(Pb_trial(1,i));
    
%     err(1,i) = (ref-solution(i,1))/ref;
    
    err(1,i) = (ref-solution(i,1));

    err(1,i) = abs(err(1,i));
    
end

error = max(err);

error = full(error);

end

