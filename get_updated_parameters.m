function  [new_params] = get_updated_parameters(a,b,history,params)
    
    function [C1,C2,C] = get_C ()
		% Calculate C1, C2, C [ref.: Eq(15) in paper]
        
         C1 = exp(history(1))/(exp(history(1)) +exp(history(2))); %+ 1/2* (params.sigma(a)+params.sigma(b))*params.emu(a)*params.emu(b)*(params.emu(b)-params.emu(a))/(params.emu(a)+params.emu(b))^3;
        
%         C1 = params.emu(a)/(params.emu(a) +params.emu(b)); %+ 1/2* (params.sigma(a)+params.sigma(b))*params.emu(a)*params.emu(b)*(params.emu(b)-params.emu(a))/(params.emu(a)+params.emu(b))^3;
        C2 = 1-C1;
        % C1 = params.emu(a)/(params.emu(a) +params.emu(b));
        C = (C1*params.alpha(1)+C2*params.beta(1))/(params.alpha(1)+params.beta(1));
    end
    
    
    %one annotator
    %term = params.alpha(1)*params.emu(a)/(params.alpha(1)*params.emu(a)+params.beta(1)*params.emu(b))-params.emu(a)/(params.emu(a)+params.emu(b));
	
    term =1 -params.emu(a)/(params.emu(a)+params.emu(b));

%      term = exp(-(params.mu(a)-params.mu(b)));
    
    %Eq 11: mu_i^(i >_k j) = mu_i + sigma_i * E(z_i)
    new_params.mu_a = params.mu(a) + params.sigma(a)*term;
    
    %Eq 12: mu_j^(i >_k j) = mu_j + sigma_j * E(z_j)
    new_params.mu_b = params.mu(b) - params.sigma(b)*term;
    
    
    
    % if(new_params.mu_a > params.mu_max)
    %     new_params.mu_a = params.mu_max;
    % end
    % if(new_params.mu_b < params.mu_min)
    %     new_params.mu_b = params.mu_min;
    % end
    
    
    %one annotator
    %     term = params.alpha(1)*params.emu(a)*params.beta(1)*params.emu(b)/(params.alpha(1)*params.emu(a)+params.beta(1)*params.emu(b))^2-params.emu(a)*params.emu(b)/(params.emu(a)+params.emu(b))^2;
    term = - params.emu(a)*params.emu(b)/(params.emu(a)+params.emu(b))^2;
    
	%Parameter kappa is a small constant to ensure the positivity of variance
    kappa = 10e-4;
	
    %Eq 13 : [sigma_i^(i >_k j)]^2
    new_params.sigma_a =  params.sigma(a)* max(1+ params.sigma(a)*term,kappa);
    %Eq 14 : [sigma_j^(i >_k j)]^2
    new_params.sigma_b = params.sigma(b)* max(1+ params.sigma(b)*term,kappa);
    
%     new_params.sigma_a =  params.sigma(a);
%     %Eq 15
%     new_params.sigma_b = params.sigma(b);
    
    %Update alpha_k and beta_k
    [C1,C2,C] = get_C ();
	
    new_params.pr = C1;
    
    
    
    
end
