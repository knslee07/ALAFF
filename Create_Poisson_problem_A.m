function [ A ] = Create_Poisson_problem_A( N )

  % 1. Create the archtypical matrix A for an N x N Poisson problem (5-point stencil)
    A=zeros(N^2,N^2);
  % 2. Set the diagonal
    % we iterate over all the rows
    for i = 1:N^2
        for j=1:N^2
            if i==j
                A(i,j)=4;
            end
        end
    end
  % 3. Set the entries of the first sub and super diagonals
    for i=1:N^2
        R = mod(i, N);
        if R==1
            idx=[i+1];
        elseif R==0
            idx=[i-1];
        else
            idx=[i-1,i+1];
        end
        A(i,idx)=-1;
    end
     
  % 4. Set the other off-diagonal entries
    for i=1:N^2
        offdiag=[i-N;i+N];
        idx=offdiag(offdiag>0 &offdiag<=N^2);
        A(i,idx)=-1;
    end
  
end




