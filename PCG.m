function [ x, niters ] = PCG( A,b,x0 )
    % I adopted the algorithm from HW 8.3.6.1. which is directly about PCG.
    % This has an initial r(0) setting different from what Jeff
    % confirmed @936. 
    % 1e-10 and it give niters=2
    % 1e-4 gives niters=11.
    L = ichol(sparse(A), struct('type','ict','droptol',1e-4,'michol','off'));
    x=x0;
    A_t=inv(L)*A*transpose(inv(L));
    r_t=inv(L)*b;
    niters=0;
    while norm(r_t) >= eps * norm(b)
        if niters==0
            p_t=r_t;
        else
            gamma_t=dot(r_t,r_t)/denominator;
            p_t=r_t+gamma_t*p_t;
        end    
        alpha_t=dot(r_t,r_t)/dot(p_t,A_t*p_t);
        x=x+alpha_t*p_t;
        denominator=dot(r_t,r_t);
        r_t=r_t-alpha_t*A_t*p_t;
        niters=niters+1;
    end
end