function [ x, niters ] = Method_of_Steepest_Descent_ichol( A,b,x0 )
    % preconditioner from the instruction
    % 1e-2 given in the instruction does not help. So, I reduced it to
    % 1e-10 and it give niters=2
    % 1e-4 gives niters=18
    L = ichol(sparse(A), struct('type','ict','droptol',1e-10,'michol','off'));
    
    % Followed the middle algorithm from Figure 8.2.5.1
    M=L*transpose(L);
    A_t=inv(L)*A*transpose(inv(L));
    b_t=inv(L)*b;
    x_t=transpose(L)*x0;
    r_t=b_t-A_t*x_t;
    niters=0;
    while norm(r_t)>= eps * norm(b_t)
        p_t=r_t;
        q_t=A_t*p_t;
        alpha_t=dot(p_t,r_t)/dot(p_t,q_t);
        x_t=x_t+alpha_t*p_t;
        r_t=r_t-alpha_t*q_t;
        x=transpose(inv(L))*x_t;
        niters=niters+1
    end
    
end


    