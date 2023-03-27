function [ x, niters ] = CG( A,b,x0 )
    % gives niters=19, which I saw another student have the exactly same
    % number on Piazza.
    x=x0;
    r=b-A*x;
    niters=0;
    while norm(r) >= eps * norm(b)
        if niters==0
            p=r;
        else
            gamma=-dot(p,A*r)/dot(p,A*p);
            p=r+gamma*p;
        end    
        alpha=dot(p,r)/dot(p,A*p);
        x=x+alpha*p;
        r=r-alpha*A*p;
        niters=niters+1;
    end
end
