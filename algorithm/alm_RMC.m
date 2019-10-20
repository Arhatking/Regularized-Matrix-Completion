function [X] = alm_RMC(M, Omega,Q,lambda,maxIter,mu1,mu2,isDisp)

%mu1 = 0.5;mu2 = 0.01;
tol = 1e-3;tolProj = 1e-3;



% initialize
[m n] = size(M);
Omega_c = ones(m,n) - Omega;
I = eye(n);
invM = inv(lambda*Q*Q'+mu2*I);


X = zeros( m, n);
E = zeros( m, n);
Z = zeros( m, n);
Y1 = zeros( m, n);
Y2 = zeros( m, n);
R1 = zeros( m, n);
R2 = zeros( m, n);

%dnorm = norm(M, 'fro');
%tolProj = tolProj * dnorm;
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
        

       
        B1=M+E-(1/mu1)*Y1;
        B2 = Z+(1/mu2)*Y2;
        tau = 1/(mu1+mu2);
        A = (mu1*B1+mu2*B2)*tau;
        
        if choosvd(n, sv) == 1
            [U Sig V] = lansvd(A, sv, 'L');
        else
            [U Sig V] = svd(A, 'econ');
        end
        diagSig = diag(Sig);
        svp = length(find(diagSig > tau));
        if svp < sv
            sv = min(svp + 1, max_sv);
        else
            sv = min(svp + round(0.05*max_sv), max_sv);
        end
		
        temp_X = U(:,1:svp)*diag(diagSig(1:svp)-tau)*V(:,1:svp)';    
        
        E = Omega_c.*(temp_X-M+(1/mu1)*Y1);
        
        Z = (mu2*temp_X - Y2)*invM;
       
            
        
        if norm(X - temp_X, 'fro') < tolProj*norm(X, 'fro')
            primal_converged = true;
        end
        X = temp_X;
        primal_iter = primal_iter + 1;
        total_svd = total_svd + 1;
               
    end
     
    temp_R1 = X - M  - E;
    temp_R2 = Z-X;        
    Y1 = Y1 + mu1*temp_R1;
    Y2 = Y2 + mu2*temp_R2;
    
    %% stop Criterion    
    stopCriterion1 = norm(temp_R1 - R1, 'fro');
    stopCriterion2 = norm(temp_R2 - R2, 'fro');
    
    if stopCriterion1 < tol && stopCriterion2 < tol
        %disp(['end at the ' num2str(iter) 'th iteration']) ;
        converged = true;
    end   
    
    R1 = temp_R1;
    R2 = temp_R2;
   
    
   if isDisp
       disp(['Iter ' num2str(iter) ' primal_iter ' num2str(primal_iter) ' #svd ' num2str(total_svd) ' Rank(X) ' num2str(svp)...
            ' stopCriterion1 ' num2str(stopCriterion1)  ' stopCriterion2 ' num2str(stopCriterion2)]);
   end
    
    if ~converged && iter >= maxIter
        %disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end


