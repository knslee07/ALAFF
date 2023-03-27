function [ x, niters ] = Method_of_Steepest_Descent( A,b,x0 )
    x=x0;
    niters=0;
    r=b-A*x;
    % I first used the next line.
    % while  any(abs(r)>=1e-24) 
    % computed r may not be real zeros. compute until the r differential is small enough.
    % takes more iterations than the below code.
    
    while any(abs(r)>=1e-24) 
        p=r;
        q=A*p;
        alpha=dot(p,r)/dot(p,q);
        x=x+alpha*p;
        r=r-alpha*q;
        niters=niters+1;
    end
    x

end
