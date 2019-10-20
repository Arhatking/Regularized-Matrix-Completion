function [X] = alm_WMC(M, Omega,Q,maxIter,mu1,mu2,isDisp)

%mu1 = 0.01;mu2 = 0.01;
tol = 1e-3;tolProj = 1e-5;
[m n] = size(M);

Omega_c = ones(m,n) - Omega;

I = eye(n);


% initialize

X = zeros( m, n);
E = zeros( m, n);
Z = zeros( m, n);
Y1 = zeros( m, n);
Y2 = zeros( m, n);
R1 = zeros( m, n);
R2 = zeros( m, n);
invM = inv(mu1*I+mu2*Q*Q');

% dnorm = norm(M, 'fro');
% tolProj = tolProj * dnorm;
total_svd = 0;


iter = 0;
converged = false;
stopCriterion = 1;
sv = 5;
svp = sv;
max_sv = min(m,n);
while ~converged       
    iter = iter + 1;
    
    % solve the primal problem by alternative projection
    primal_converged = false;
    primal_iter = 0;
    sv = min(sv + round(max_sv * 0.1),max_sv);
    while primal_converged == false && primal_iter <10
        
        
        temp_X = (mu1*M+mu1*E-Y1+mu2*Z*Q'+Y2*Q')*invM;
        E = Omega_c.*(temp_X-M+Y1/mu1);
        
        
        
        if choosvd(n, sv) == 1
            [U Sig V] = lansvd(temp_X*Q - (1/mu2)*Y2, sv, 'L');
        else
            [U Sig V] = svd(temp_X*Q - (1/mu2)*Y2, 'econ');
        end
        diagSig = diag(Sig);
        svp = length(find(diagSig > 1/mu2));
        if svp < sv
            sv = min(svp + 1, max_sv);
        else
            sv = min(svp + round(0.05*max_sv), max_sv);
        end
        Z = U(:,1:svp)*diag(diagSig(1:svp)-1/mu2)*V(:,1:svp)';    
        
        
        if norm(X - temp_X, 'fro') < tolProj*norm(X,'fro')
            primal_converged = true;
        end
        X = temp_X;
        primal_iter = primal_iter + 1;
        total_svd = total_svd + 1;
               
    end
     
    temp_R1 = X - M  - E;
    temp_R2 = Z-X*Q;        
    Y1 = Y1 + mu1*temp_R1;
    Y2 = Y2 + mu2*temp_R2;
    
    %% stop Criterion    
    stopCriterion1 = norm(temp_R1 - R1, 'fro');
    stopCriterion2 = norm(temp_R2 - R2, 'fro');
    if stopCriterion1 < tol && stopCriterion2 < tol
        converged = true;
    end    
    R1 = temp_R1;
    R2 = temp_R2;
    
   if isDisp
       disp(['Iteration' num2str(iter) ' #svd ' num2str(total_svd) ' Rank(X) ' num2str(svp)...
            ' stopCriterion1 ' num2str(stopCriterion1)  ' stopCriterion2 ' num2str(stopCriterion2)]);
   end
    
    if ~converged && iter >= maxIter
        %disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end


