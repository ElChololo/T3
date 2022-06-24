classdef T3_obj
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        base
        Y
        regresores
    end
    
    methods
        function obj = T3_obj(xlsx)

            df =  xlsread(xlsx);
            %obtenemos las variables relevantes
            obj.base = [df(:,2) df(:,27) df(:,28) df(:,5) df(:,4)];
            %revisiom de missing values
            assert(sum(isnan(obj.base),'all') == 0)
            %asignaciom de variable dependiente del modelo
            obj.Y = obj.base(:,1);
            %asignacion de constante y regresores del modelo
            obj.regresores = [ones(size(obj.base,1),1) obj.base(:,2:end)];
            
            
        end
        
        function mu = mu_i(obj,x,beta)

            mu = exp(x*beta);
        end
        function log_likelihood = log_likelihood(obj,x,y,beta)
            log_likelihood = sum(obj.mu_i(x,beta) + y * (x.*beta) - ln(factorial(y)));
            
            
        end
        
        function gradiente = gradiente(obj,x,y,beta)
           gradiente = x'*(y-obj.mu_i(x,beta));
           %gradiente = sum(diag(y-obj.mu_i(x,beta)) * x,1)'; 
        end
        
        function hessiano = hessiano(obj,x,beta)
            hessiano = (-1)* x'*diag(obj.mu_i(x,beta))*x ;
        end
        function varscore = varscore(obj,x,y,beta)
            x=x';
            input_score = @(x,y,beta) (y-exp(x'*beta))*x;
            varscore = zeros(5,5);
            for ii = 1:length(y)
                sum_aux = input_score(x(:,ii),y(ii),beta);
                auxvarscore = sum_aux*sum_aux';
                varscore = varscore + auxvarscore;
            end

            
        end
        function beta_est = MetNum(obj,tipo,iter,guess,x,y)
            beta_t = guess';
            if tipo == "NR"
                for ii=1:iter
                    beta_t1 = beta_t - (obj.hessiano(x,beta_t)\ obj.gradiente(x,y,beta_t));

                    if sum(abs(beta_t1-beta_t)<1e-12)
                        beta_est = beta_t1;
                        fprintf("Convergencia en %d iteraciones \n",ii)
                        fprintf("Resultado Iteracion [");
                        fprintf("%g, ", beta_t1(1:end-1));
                        fprintf("%g]\n",beta_t1(end));
                        return
                    else
                        beta_t = beta_t1;
                    end
               
                end
            elseif tipo == "BHHH"
                for ii=1:iter
                    gradi= obj.gradiente(x,y,beta_t);
                    varscore = obj.varscore(x,y,beta_t);
                    beta_t1 = beta_t + varscore \ gradi;

                    if sum(abs(beta_t1-beta_t)<1e-8)
                        beta_est = beta_t1;
                        fprintf("Convergencia en %d iteraciones \n",ii)
                        fprintf("Resultado Iteracion [");
                        fprintf("%g, ", beta_t1(1:end-1));
                        fprintf("%g]\n",beta_t1(end));
                        return

                    else
                        beta_t = beta_t1;
                       

                    end
               
                end
            else
                error("Not numerical Method provided")
            end

        end
        
        function b_0 = NR_R(obj,iter,y,beta)
           beta_t = beta;
           for ii = 1:iter
              prim_derivada = mean(y)-exp(beta_t);
              seg_derivada = -exp(beta_t);
              beta_t1 = beta_t + prim_derivada/ ((-1)*seg_derivada);
              if abs(beta_t1-beta_t)<1e-8
                b_0 = beta_t1;
                fprintf("Convergencia en %d iteraciones \n",ii)
                fprintf("Resultado Iteracion [");
                fprintf("%g, ", beta_t1(1:end-1));
                fprintf("%g]\n",beta_t1(end));
                return

              else
                 beta_t = beta_t1;


              end
           end
        end
        
        
            
    end
end
