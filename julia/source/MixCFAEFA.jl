module MixCFAEFA


using Distributions, LinearAlgebra, Clustering, Base.Threads, JLD2, Random, CategoricalArrays, Combinatorics

export CFAmixEFAstandardno, Entropy, care_sim, CarelessGen, CFAmixEFA, boots_CFAmEFA_careless, boots_CFAmEFA_faking, classperf, measuresp, inv_lemma_1_fun, inv_lemma_2_fun, det_lemma_1_fun, det_lemma_2_fun, CFA


########################################################################################## Simulation functions #######################################################################################
function Entropy(res,Y,p,q,n,K,ncov)
    Lambda_1 = reshape(res[1:(p*q)],p,q)
    Lambda_2 = reshape(res[(p*q+1):((p*q)+(p*K))],p,K)
    Theta_delta = Diagonal(res[((p*q)+(p*K)+1):((p*q)+(p*K)+p)])
    Psi_delta = Diagonal(res[((p*q)+(p*K)+p+1):((p*q)+(p*K)+p+p)])
    Phi = reshape(res[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)],q,q)
    kappa = res[((p*q)+(p*K)+p+p+q^2+1):((p*q)+(p*K)+p+p+q^2+n)]
    nu = res[((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+K)]


    inv_lemma_1 = inv_lemma_1_fun(q,Lambda_1,Theta_delta,Phi)
    inv_lemma_2 = inv_lemma_2_fun(K,Lambda_2,Psi_delta)      
    det_lemma_1 = det_lemma_1_fun(Lambda_1,Theta_delta,Phi)
    det_lemma_2 = det_lemma_2_fun(K,Lambda_2,Psi_delta)
    b1 = Phi*Lambda_1'*inv(Lambda_1*Phi*Lambda_1'+Theta_delta);
    b2 = Lambda_2'*inv(Lambda_2*Lambda_2' + Psi_delta);
    P1p = zeros(n)
    P1 = zeros(n)
    for i in 1:n
        P1p[i] = ((kappa[i] .* ((exp((-1/2) * Y[:,i]'*inv_lemma_1*Y[:,i]))/(sqrt((2 .* pi)^(p) *det_lemma_1))))/((kappa[i] .* ((exp((-1/2) * Y[:,i]'*inv_lemma_1*Y[:,i]))/(sqrt((2 .* pi)^(p) *det_lemma_1)))).+((exp((-1/2) * (Y[:,i].-Lambda_2*nu)'*inv_lemma_2*(Y[:,i].-Lambda_2*nu)))/(sqrt((2 .* pi)^(p) .*det_lemma_2))).-(kappa[i] .*((exp((-1/2) *(Y[:,i].-Lambda_2*nu)'*inv_lemma_2*(Y[:,i].-Lambda_2*nu)))/(sqrt((2 .* pi)^(p) .*det_lemma_2))))))
    end
    P1p = P1p .+ 1e-7
    P1p[findall(P1p.>1)].=0.9999999999999
    P0 = zeros(n)
    P1p[findall(P1p.>1)].=0.9999999999999
    for j=1:n
        P0[j] = (1-P1p[j])*log(1-P1p[j])
        P1[j] = P1p[j]*log(P1p[j])
    end
    P0[isnan.(P0)].=0.0
    Ek = - (sum(P1) + sum(P0))
    Entropy = 1-(Ek/(n*log(2)))
    return Entropy
end

    
function CFAmixEFAstandardno(Y::Matrix{Float64}, q::Int64,  K::Int64, maxIter::Int64, X::Matrix{Float64}, maxNumRestart::Int64, mu_hat::Int64=1,L_str=nothing)
    local re::Int64, cal::Int64, non_conv::Bool, err2::Bool, sum_diff_tf::Bool, retr::Int64, res::Vector{Float64}, start_vals::Vector{Float64}, DomainError_warn::Int64, p::Int64, n::Int64
   
    
    ncov = size(X,2); p = size(Y,1); n = size(Y,2); iter_err = 1; code_error = 0.0;  err_0_count = 0; err_1_count = 0; err_2_count = 0; err_3_count = 0; err_4_count = 0; err_5_count = 0;
    while code_error != 8.0
        cal = (p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n
        
        if iter_err >= maxNumRestart
            code_error = 7.0 
            return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n),code_error, err_0_count, err_1_count, err_2_count, err_3_count, err_4_count, err_5_count) 
        end
        
        # Starting values
        start_vals = starting_val(n,p,q,K,X,mu_hat,L_str)
         
        # EM 
        res = EM2(Y,X,q,K,p,n,ncov,maxIter,reshape(start_vals[1:(p*q)],p,q), reshape(start_vals[(p*q+1):((p*q)+(p*K))],p,K), Diagonal(start_vals[((p*q)+(p*K)+1):((p*q)+(p*K)+p)]), Diagonal(start_vals[((p*q)+(p*K)+p+1):((p*q)+(p*K)+p+p)]), reshape(start_vals[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)],q,q), vec(start_vals[((p*q)+(p*K)+p+p+q^2+1):((p*q)+(p*K)+p+p+q^2+n)]), vec(start_vals[((p*q)+(p*K)+p+p+q^2+n+1):((p*q)+(p*K)+p+p+q^2+n+ncov)]), reshape(start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q)],p,q),start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q)],start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+K)],mu_hat)
  
        

        # Starting values
        # repli = 10
        # ress = Matrix{Float64}(undef,cal,repli); liks = Vector{Float64}(undef,repli); 
        # for i in 1:repli
            
        #     start_vals = starting_val(n,p,q,K,X,mu_hat,L_str)
        #     ress[:,i] = EM(Y,X,q,K,p,n,ncov,maxIter,reshape(start_vals[1:(p*q)],p,q), reshape(start_vals[(p*q+1):((p*q)+(p*K))],p,K), Diagonal(start_vals[((p*q)+(p*K)+1):((p*q)+(p*K)+p)]), Diagonal(start_vals[((p*q)+(p*K)+p+1):((p*q)+(p*K)+p+p)]), reshape(start_vals[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)],q,q), vec(start_vals[((p*q)+(p*K)+p+p+q^2+1):((p*q)+(p*K)+p+p+q^2+n)]), vec(start_vals[((p*q)+(p*K)+p+p+q^2+n+1):((p*q)+(p*K)+p+p+q^2+n+ncov)]), reshape(start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q)],p,q),start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q)],start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+K)],mu_hat)[1:cal]
        #     liks[i] = ress[(p*q)+(p*K)+p+p+q^2+n+1,i]
        # end
        # println(liks[Bool.(1 .-(isnan.(liks)))])
        # liktr = liks[Bool.(1 .-(isnan.(liks)))]

        # if sum(isnan.(liks)) != repli
        #     rex = ress[:,Bool.(1 .-(isnan.(liks)))][:,findall((liktr).==maximum(liktr))]
        #     # EM 
        #     res = EM(Y,X,q,K,p,n,ncov,maxIter,reshape(rex[1:(p*q)],p,q), reshape(rex[(p*q+1):((p*q)+(p*K))],p,K), Diagonal(rex[((p*q)+(p*K)+1):((p*q)+(p*K)+p)]), Diagonal(rex[((p*q)+(p*K)+p+1):((p*q)+(p*K)+p+p)]), reshape(rex[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)],q,q), vec(rex[((p*q)+(p*K)+p+p+q^2+1):((p*q)+(p*K)+p+p+q^2+n)]), vec(rex[((p*q)+(p*K)+p+p+q^2+n+1):((p*q)+(p*K)+p+p+q^2+n+ncov)]), reshape(rex[((p*q)+(p*K)+p+p+q^2+n+ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q)],p,q),rex[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q)],rex[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+K)],mu_hat)
        # elseif sum(isnan.(liks)) == repli
        #     # EM 
        #     start_vals = starting_val(n,p,q,K,X,mu_hat,L_str)
        #     res = EM(Y,X,q,K,p,n,ncov,maxIter,reshape(start_vals[1:(p*q)],p,q), reshape(start_vals[(p*q+1):((p*q)+(p*K))],p,K), Diagonal(start_vals[((p*q)+(p*K)+1):((p*q)+(p*K)+p)]), Diagonal(start_vals[((p*q)+(p*K)+p+1):((p*q)+(p*K)+p+p)]), reshape(start_vals[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)],q,q), vec(start_vals[((p*q)+(p*K)+p+p+q^2+1):((p*q)+(p*K)+p+p+q^2+n)]), vec(start_vals[((p*q)+(p*K)+p+p+q^2+n+1):((p*q)+(p*K)+p+p+q^2+n+ncov)]), reshape(start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q)],p,q),start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q)],start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+K)],mu_hat)
        # end
        
        # Controls for non covnergence
        cal = (p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n+1
        sum_diff_tf = !(res[(cal+1):(cal+1)][1]==1.0)
        non_conv = res[(cal+1+1):(cal+1+1)][1]>=maxIter
        code_error = res[(cal):(cal)][1]
        if code_error == 0.0  && iter_err < maxNumRestart
            err_0_count += 1 
            iter_err += 1
        elseif code_error == 1.0  && iter_err < maxNumRestart
            err_1_count += 1 
            iter_err += 1
        elseif code_error == 2.0  && iter_err < maxNumRestart
            err_2_count += 1 
            iter_err += 1
        elseif code_error == 3.0 && iter_err < maxNumRestart
            err_3_count += 1 
            iter_err += 1
        elseif non_conv && iter_err < maxNumRestart
            err_4_count += 1 
            iter_err += 1
        elseif sum_diff_tf && iter_err < maxNumRestart
            err_5_count += 1 
            iter_err += 1
        end
    end 
    return res[1:(cal-1)]
end

# Expectation-Maximization algorithm
function EM2(Y::Matrix{Float64},X::Matrix{Float64},q::Int64, K::Int64, p::Int64, n::Int64, ncov::Int64, maxIter::Int64, Lambda_1::Matrix{Float64}, Lambda_2::Matrix{Float64}, Theta_delta::Diagonal{Float64, Vector{Float64}}, Psi_delta::Diagonal{Float64, Vector{Float64}}, Phi::Matrix{Float64}, kappa::Vector{Float64}, betas::Vector, L_str::Matrix, mu::Vector{Float64}, nu::Vector{Float64}, mu_hat::Int64)
    local llik::Vector{Float64}, Ez1::Vector{Float64}, diff_llik::Float64, oldllik::Float64, dtol::Float64, t::Int64, err::Vector{Bool}, err2::Bool, sum_diff_tf::Bool, sum_diff::Vector{Bool}, old_diff::Float64,  ESxixi::Matrix{Float64}, ESetaeta::Matrix{Float64}, ESyxi::Matrix{Float64}, ESyeta::Matrix{Float64}, ESyy0::Matrix{Float64}, ESyy1::Matrix{Float64}, Ek0::Float64, Ek1::Float64, Ez0::Vector{Float64}, B2::Matrix{Float64}, b2::Matrix{Float64}, B1::Matrix{Float64}, b1::Matrix{Float64}, Exi4::Matrix{Float64}, Eeta4::Matrix{Float64}
    
    zero_indices = [i + (j - 1) * p for j in 1:K for i in 1:(j-1)]
            
    t = 1; err2 = false; llik = zeros(Float64, maxIter); oldllik = -Inf; err=zeros(Bool,maxIter); sum_diff = zeros(Bool,maxIter); sum_diff_tf = true;  diff_llik = 1.0; dtol = 0.001; code_error = 8.0;
    while diff_llik>dtol && t<=maxIter && sum_diff_tf
            inv_lemma_1 = inv_lemma_1_fun(q,Lambda_1,Theta_delta,Phi)
            inv_lemma_2 = inv_lemma_2_fun(K,Lambda_2,Psi_delta)      
            det_lemma_1 = det_lemma_1_fun(Lambda_1,Theta_delta,Phi)
            det_lemma_2 = det_lemma_2_fun(K,Lambda_2,Psi_delta)

            if det_lemma_1<=0 || det_lemma_2<=0
                code_error = 0.0
                return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n), code_error, sum_diff_tf, t-1) 
            end
            # println(t)
            #E-step
            # b1 = b1_fun(q,Lambda_1,Phi,Theta_delta,mu)
            # B1 = B1_fun(Lambda_1,Phi,Theta_delta,b1,mu)
            # b2 = b2_fun(K,Lambda_2,Psi_delta,nu)
            # B2 = B2_fun(K,Lambda_2,Psi_delta,b2,nu)
            
            b1 = b1_fun(q,Lambda_1,Phi,Theta_delta)
            B1 = B1_fun(Lambda_1,Phi,Theta_delta,b1)
            b2 = b2_fun(K,Lambda_2,Psi_delta)
            B2 = B2_fun(K,Lambda_2,Psi_delta,b2)
            
            Ez1 = Expectz(n,p,Y,kappa,inv_lemma_1,inv_lemma_2, det_lemma_1,det_lemma_2,Lambda_1,Lambda_2,mu,nu)

            Ez0 = Expectz0(Ez1)
            Ek1 = Ekz1(Ez1)
            Ek0 = Ekz0(Ez0)
            
            ESyy1 = ESyyz1(n,p,Y,Ez1)  
            ESyy0 = ESyyz0(n,p,Y,Ez0)
            
            Ezy1 = ESzy1(n,p,Y,Ez1,mu,Lambda_1)
            Ezy0 = ESzy0(n,p,Y,Ez0,nu,Lambda_2)
            Ey1 = ESyone(n,p,Y,Ez1)
            Ey0 = ESyzero(n,p,Y,Ez0)
            
            Eeta = Eeta1(p,n,Ek1,Ez1,Y,b1, mu, Lambda_1)
            Exi = Exi0(p,n,Ek0,Ez0,Y,b2, nu, Lambda_2)
            
            Eetay1 =  ESetay1(n,p,Y,Ez1,Lambda_1,mu,b1,Ey1)
            Exiy0 =  ESxiy0(n,p,Y,Ez0,Lambda_2,nu,b2,Ey0)

            ESetaeta =  ESetaeta1(n,p,Y,Ez1,Lambda_1,mu,b1,B1)
            ESxixi = ESxixi0(n,p,Y,Ez0,Lambda_2,nu,b2,B2)
         
            if isequal(sum(ESetaeta),NaN) || isequal(sum(Eetay1),NaN) || isnan(ESxixi[1,1]) || isequal(sum(ESxixi),0)
                code_error = 0.0
                return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n), code_error, sum_diff_tf, t-1) 
            else
                #M-step
                Lambda_1 = Lambda_1_mstep(Eetay1,ESetaeta,L_str)
                Lambda_2 = Lambda_2_mstep(Exiy0,ESxixi)
                L2 = vec(Lambda_2)
                L2[zero_indices,:].=0.
                Lambda_2 = reshape(L2,p,K)

                Phi = Phi_mstep(ESetaeta,Eeta,Ek1,mu,n)

                Theta_delta =  Diagonal(1 .-Diagonal((Lambda_1*Phi*Lambda_1')));# Theta_delta_mstep(n,Lambda_1,ESyy1,Eetay1,Ek1)
                Psi_delta = Psi_delta_mstep(n,Lambda_2,ESyy0,Exiy0,Ek0)#Psi_delta_mstep(Y,b2,B2,nu,n,p,K,Lambda_2,ESyy0,Exiy0,Ez0)
                

                if mu_hat == 1
                    mu = ((sum(Eeta,dims=2)/n)/(Ek1/n))[1:q,1]
                elseif mu_hat == 0
                    mu = zeros(q)
                end
                nu = zeros(K)
                
            
            
                #If no covariates are provided
                if ncov == 1
                    kappa_1 = Ek1/n
                    kappa = fill(kappa_1,n)
                    betas = log(kappa_1/(1 - kappa_1))
                else
                    betas = NewtonRaphson(betas_gradient, betas, X, Ez1, Ez0, ncov)
                    kappa = exp.(X*betas)./(exp.(X*betas).+1)
                end         
            
                inv_lemma_1 = inv_lemma_1_fun(q,Lambda_1,Theta_delta,Phi)
                inv_lemma_2 = inv_lemma_2_fun(K,Lambda_2,Psi_delta)      
                det_lemma_1 = det_lemma_1_fun(Lambda_1,Theta_delta,Phi)
                det_lemma_2 = det_lemma_2_fun(K,Lambda_2,Psi_delta)
                #Controls for non-convergence
                if sum(Theta_delta.<0) > 0 || det(Psi_delta) < 0 || isequal(sum(Theta_delta),NaN) || isequal(sum(Ez0),NaN) || sum(diag(inv_lemma_1).<0) > 0 || sum(diag(inv_lemma_2).<0) > 0 || det_lemma_1 < 0 || det_lemma_2<0 
                    code_error = 0.0
                    return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n), code_error, sum_diff_tf, t-1) 
                elseif isnan(betas[1]) 
                    code_error = 1.0
                    return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n), code_error, sum_diff_tf, t-1) 
                elseif det(Phi)<0 
                    code_error = 2.0
                    return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n), code_error, sum_diff_tf, t-1) 
                else
                    old_diff = diff_llik
                    llik[t] = llik_compute(Y,n,p,q,K,Lambda_1, Lambda_2, kappa,mu,nu,inv_lemma_1,inv_lemma_2,det_lemma_1,det_lemma_2)

                    diff_llik = abs(oldllik - llik[t])
                    if !(oldllik<llik[t])
                        code_error = 3.0
                        return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n), code_error, sum_diff_tf, t-1) 
                    end

                    sum_diff[t] = diff_llik>old_diff               
                    oldllik = llik[t]
                    t=t+1;
                    sum_diff_tf = sum(sum_diff) < 150            
                end
            end
        end
        Exitrue=Exi3(p,n,Y,b2, nu, Lambda_2)
        Eetatrue=Eeta3(p,n,Y,b1, mu, Lambda_1)
    return vcat(vec(Lambda_1), vec(Lambda_2), vec(diag(Theta_delta)), vec(diag(Psi_delta)), vec(Phi), kappa, llik[t-1], Ez1, betas, mu, nu, vec(Exitrue), vec(Eetatrue), code_error, sum_diff_tf, t-1);
end     

function care_sim(design,B,maxit, maxNumRestart)
    local res_B
    ics = JLD2.load("/home/niccolo.cao/storage/careless/data/ics.jld2")["ics"] 
    for d in 1:size(design,1)
        @threads for b in 1:B
            println("Replication: n. ",b)
            n = 1000
            p = 9
            q = 3
            c = design[d,1]
            K = design[d,2]
            ncov = 2
            
            A =  Matrix{Int64}(I,q,q)
            str = ones(Int64,Int64(p/q))
            L_str = kron(A, str);
            
            X = JLD2.load(string("/home/niccolo.cao/storage/careless/data/data_",d,"_",b,".jld2"))["data"]#X = JLD2.load(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"))["data"]
            
            Y = Matrix{Float64}(X);

            #some minimal variation
            for i in 1:size(Y,1)
                Y[i,:] = Y[i,:] + rand(Uniform(-0.02,0.02),size(Y,2))
            end
            #standarization of observed variables
            for i in 1:size(Y,2)
                Y[:,i]=(Y[:,i].-mean(Y,dims=1)[i])/std(Y,dims=1)[i]
            end

            Y = Matrix(Y');
            repli = 10
            ress = Matrix{Float64}(undef,(p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n,repli); liks = Vector{Float64}(undef,repli);  
            for i in 1:repli
                ress[:,i] = CFAmixEFA(Y, q, K, maxit, ics, maxNumRestart, 0,L_str)
                liks[i] = ress[((p*q)+(p*K)+p+p+q^2+n+1):((p*q)+(p*K)+p+p+q^2+n+1),i][1]
            end
            liktr = liks[Bool.(1 .-(isnan.(liks)))]
            res_B = ress[:,Bool.(1 .-(isnan.(liks)))][:,findall((liktr).==maximum(liktr))]
            
            Ez1 = res_B[((p*q)+(p*K)+p+p+q^2+n+1+1):((p*q)+(p*K)+p+p+q^2+n+1+n)]
            bets = res_B[((p*q)+(p*K)+p+p+q^2+n+1+n+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov)]
    
            if sum(Ez1)/n<0.5
                res_B[((p*q)+(p*K)+p+p+q^2+n+1+1):((p*q)+(p*K)+p+p+q^2+n+1+n)] .= 1 .- Ez1
                res_B[((p*q)+(p*K)+p+p+q^2+n+1+n+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov)] .= .- bets
                println("label unswitched: ",bets)
            end
            
            jldsave(string("/home/niccolo.cao/storage/careless/res/res_",d,"_",b,".jld2"),res=res_B)#jldsave(string("/home/nc/MEGA/LatentVariables/res/res_",d,"_",b,".jld2"),res=res_B)
        end
    end
end




function CarelessGen(design,B)
    local kappa0, t10, TruePar::Vector{Float64}, Lambda0::Matrix{Float64}, Theta_delta0::Diagonal{Float64, Vector{Float64}}, Phi0::Matrix{Float64}
    
    seed = vcat(322011611885,11955714332337,14889268790486,19321661366842, 4539386885746,17304354949517, 1981879088971,14248930143082,16766904013492, 1434926186704,22695264190645,17517986428637,17678641845771, 1581593402381, 3465694164649, 9106731443143, 9729393064462,12160146297624, 4420744722171,16184395927257,19766499429427,20201406362254,19233052068296, 4730430435062)
    n = 1000
    p = 9
    q = 3

    ## True CFA
    leftthresh4 = vcat(-1.83,-1.67,-1.86,-1.75,-1.8,-1.62,-1.68,-1.72,-1.66,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.83,1.67,1.86,1.75,1.8,1.62,1.68,1.72,1.66)
    # leftthresh5 = vcat(-1.09,-1.13,-1.68,-2.23,-2.21,-1.61,-2.03,-1.55,-1.52,0.65,0.62,0.06,-0.49,-0.47,0.14,-0.28,0.19,0.22,2.40,2.36,1.81,1.26,1.28,1.88,1.46,1.94,1.97)
    # leftthresh6 = vcat(-1.09,-1.13,-1.68,-2.23,-2.21,-1.61,-2.03,-1.55,-1.52,0.65,0.62,0.06,-0.49,-0.47,0.14,-0.28,0.19,0.22,2.40,2.36,1.81,1.26,1.28,1.88,1.46,1.94,1.97)
    leftthresh10 = vcat(-2.17,-1.66,-2.23,-2.26,-2.19,-1.33,-2.09,-1.28,-1.67,-1.73,-1.22,-1.79,-1.82,-1.75,-0.9,-1.65,-0.84,-1.23,-1.29,-0.78,-1.36,-1.38,-1.32,-0.46,-1.22,-0.4,-0.79,-0.86,-0.35,-0.92,-0.94,-0.88,-0.02,-0.78,0.03,-0.35,-0.42,0.09,-0.48,-0.51,-0.44,0.42,-0.34,0.47,0.08,0.02,0.53,-0.04,-0.07,-0.0,0.85,0.1,0.91,0.52,0.46,0.97,0.39,0.37,0.43,1.29,0.53,1.35,0.96,0.89,1.4,0.83,0.81,0.87,1.73,0.97,1.78,1.4,1.33,1.84,1.27,1.24,1.31,2.17,1.41,2.22,1.83)
    # t10 = hcat(fill(-Inf,p),leftthresh,leftthresh.+1,leftthresh.+2,leftthresh.+3.5, fill(Inf,p))
    A = reshape(Diagonal(ones(q)),q,q);
    L_str = kron(A,ones(Int64(p/q)));
    Lambda0 = [0.69, 0.43, 0.7, 0.36, 0.78, 0.80, 0.56, 0.42, 0.67] .* L_str;
    Phi0 = [1 0.12 0.13;0.12 1 -0.21;.13 -0.21 1]; #Diagonal(ones(q,q))
    Theta_delta0 = Diagonal(1 .-Diagonal((Lambda0*Phi0*Lambda0'))); #Diagonal([0.56,1.1,0.53,1.2,0.46,0.4,0.73,1.15,59])
    betas085 = vcat(log(0.85/(1 - .85)),0.8)#,rand(Uniform(-1.5,1.5),ncov-1))
    betas095 = vcat(log(0.95/(1 - .95)),0.8)#,rand(Uniform(-1.5,1.5),ncov-1))
    # npartot = Int(p*(c-1-2) + p + p + q*(q-1)/2 + p) #per avere vettori output con dimensioni uniformi nonostante la differenza tra disegni sperimentali
    TruePar = vcat(leftthresh4,leftthresh10,vec(sum(Lambda0,dims=2)),Phi0[tril!(trues(size(Phi0)), -1)],betas085,betas095)#, vec(zeros(numpar-npartot))) 
    
    jldsave("/home/niccolo.cao/storage/careless/data/TruePar.jld2",TruePars=TruePar) #jldsave("/home/nc/MEGA/LatentVariables/data/TruePar.jld2",TruePars=TruePar) 

    
    ics = hcat(ones(Float64,n),rand(Normal(0,1),n)); 
    jldsave(string("/home/niccolo.cao/storage/careless/data/ics.jld2"),ics=ics) 

    for d in 1:size(design,1)
        c = design[d,1]
        K = design[d,2]
        propo = design[d,3]
        careless = design[d,4]
            

        # Generation
        Random.seed!(seed[d])
        @threads for b in 1:B
            Eta0 = rand(MvNormal(zeros(q),Phi0),n);
            Delta0 = rand(MvNormal(zeros(p),Theta_delta0),n);
            Ystarcfa = Eta0'*Lambda0' + Delta0';
                
            

            if c == 4.
                t10 = hcat(fill(-Inf,p),reshape(leftthresh4,p,c-1),fill(Inf,p))
            elseif c == 10.
                t10 = hcat(fill(-Inf,p),reshape(leftthresh10,p,c-1),fill(Inf,p))
            end

            Ycfa = Matrix{Int64}(undef,n,p);
            for j in 1:p
                Ycfa[:,j] = cut(Ystarcfa[:,j],t10[j,:],labels=1:c) 
            end
            


            ncov = 2
            if propo == 0.85
                kappa0 = exp.(ics*betas085)./(exp.(ics*betas085).+1)
            elseif propo == 0.95
                kappa0 = exp.(ics*betas095)./(exp.(ics*betas095).+1)
            end
            z = rand.(Binomial.(1, kappa0))

            carel = Matrix{Int64}(undef,n,p);
            if careless == 0. #random
                carel = reshape(wsample([1,2,3,4],[0.25,0.25,0.25,0.25],n*p,replace=true),n,p)
            elseif careless == 1. #complete straightlining 
                carel = reshape(repeat(wsample([1,2,3,4],[0.25,0.25,0.25,0.25],n,replace=true),outer=p),n,p)
            end 

            X=zeros(Int64,n,p)
            for i=1:n
                if z[i]==1
                    X[i,:] = vec(Ycfa[i,:]);
                elseif z[i]==0
                    X[i,:] = vec(carel[i,:]);
                end
            end 
            jldsave(string("/home/niccolo.cao/storage/careless/data/data_",d,"_",b,".jld2"),data=X)
            jldsave(string("/home/niccolo.cao/storage/careless/data/z_",d,"_",b,".jld2"),z=z)
            #jldsave(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"),data=X)  
        end
    end
end




function measuresp(design,B)
    local betas0
    D = size(design,1)
    Lambda_1_rmse = Array{Float64}(undef,D,B); Theta_delta_rmse = Array{Float64}(undef,D,B); Phi_rmse = Array{Float64}(undef,D,B); kappa_rmse = Array{Float64}(undef,D,B); 
    Lambda_1_bias = Array{Float64}(undef,D,B); Theta_delta_bias = Array{Float64}(undef,D,B); Phi_bias = Array{Float64}(undef,D,B); kappa_bias = Array{Float64}(undef,D,B); 
    BACC = Array{Float64}(undef,D,B);  SE = Array{Float64}(undef,D,B);  SP = Array{Float64}(undef,D,B);  MCC = Array{Float64}(undef,D,B); ARI = Array{Float64}(undef,D,B);
    bL1 = Array{Float64}(undef,D,B); bPh = Array{Float64}(undef,D,B); bTh = Array{Float64}(undef,D,B); bbe = Array{Float64}(undef,D,B); 
    Lambda_1_rmse = Array{Float64}(undef,D,B); Phi_rmse = Array{Float64}(undef,D,B); Thet_rmse = Array{Float64}(undef,D,B); betas_rmse = Array{Float64}(undef,D,B);  
    seed = vcat(322011611885,11955714332337,14889268790486,19321661366842, 4539386885746,17304354949517, 1981879088971,14248930143082,16766904013492, 1434926186704,22695264190645,17517986428637,17678641845771, 1581593402381, 3465694164649, 9106731443143, 9729393064462,12160146297624, 4420744722171,16184395927257,19766499429427,20201406362254,19233052068296, 4730430435062)
    

    ics = JLD2.load("/home/niccolo.cao/storage/careless/data/ics.jld2")["ics"] 
    ds = 0
    for d in 1#:size(design,1)
        ds += 1
        n = 1000
        p = 9
        q = 3
        c = design[d,1]
        K = design[d,2]
        propo = design[d,3]
        ncov = 2
            
        A = reshape(Diagonal(ones(q)),q,q);
        L_str = kron(A,ones(Int64(p/q)));
        
        Random.seed!(seed[d])
        ok = 0
        for b in 1:B
            res = JLD2.load(string("/home/niccolo.cao/storage/careless/res/res_",d,"_",b,".jld2"))["res"]#X = JLD2.load(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"))["data"]
            Lambda_1 = sum(reshape(res[1:(p*q)],p,q),dims=2)
            Theta_delta = res[((p*q)+(p*K)+1):((p*q)+(p*K)+p)]
            Phi = reshape(res[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)],q,q)
            Phi = Phi[tril(trues(size(Phi)),-1)]
            kappa = res[((p*q)+(p*K)+p+p+q^2+1):((p*q)+(p*K)+p+p+q^2+n)]
            Ez1 = round.(Int64,res[((p*q)+(p*K)+p+p+q^2+n+1+1):((p*q)+(p*K)+p+p+q^2+n+1+n)])
            betas = res[((p*q)+(p*K)+p+p+q^2+n+1+n+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov)]

            Lambda0 = [0.69, 0.43, 0.7, 0.36, 0.78, 0.80, 0.56, 0.42, 0.67];
            Phi0vec = [0.12, 0.13, -0.21]; #Diagonal(ones(q,q))
            Phi0 = [1 0.12 0.13;0.12 1 -0.21;.13 -0.21 1]
            Theta_delta0 = diag(Diagonal(1 .-Diagonal(((Lambda0.*L_str)*Phi0*(Lambda0.*L_str)'))));
            if propo == 0.85
                betas0 = vcat(log(0.85/(1 - .85)),0.8)
            elseif propo == 0.95
                betas0 = vcat(log(0.95/(1 - .95)),0.8)
            end
            kappa0 = exp.(ics*betas0)./(exp.(ics*betas0).+1)
            z = JLD2.load(string("/home/niccolo.cao/storage/careless/data/z_",d,"_",b,".jld2"))["z"]

            
            if sum(Ez1)/n<0.5
                Ez1 = 1 .- Ez1
                betas = .- betas
                println("label unswitched: ",betas)
            end
            
            if Ez1[1] !== NaN
                ok += 1
                cls = classperf(z,Ez1)
                println(cls)
                BACC[ds,b] = Float64(cls[2,1])
                MCC[ds,b] = Float64(cls[2,2])
                SE[ds,b] = Float64(cls[2,3])
                SP[ds,b] = Float64(cls[2,4])
                ARI[ds,b] = Float64(cls[2,5])
                bL1[ds,b] = mean(Lambda_1.-Lambda0)
                bPh[ds,b] = mean(Phi.-Phi0vec)
                bTh[ds,b] = mean(Theta_delta.-Theta_delta0)
                bbe[ds,b] = mean(betas.-betas0)
                Lambda_1_rmse[ds,b] = sqrt(mean((Lambda_1 .- Lambda0)*(Lambda_1.-Lambda0)'))
                Phi_rmse[ds,b] = sqrt(mean((Phi .- Phi0vec)*(Phi.-Phi0vec)'))
                Thet_rmse[ds,b] = sqrt(mean((Theta_delta .- Theta_delta0)*(Theta_delta .- Theta_delta0)'))
                betas_rmse[ds,b] = sqrt(mean((betas .- betas0)*(betas.-betas0)'))
                # println(Lambda_1,Lambda0," ",Phi,Phi0vec," ",Theta_delta,Theta_delta0," ",betas,betas0)
            else
                BACC[ds,b] = NaN
                MCC[ds,b] = NaN
                SE[ds,b] = NaN
                SP[ds,b] = NaN
                ARI[ds,b] = NaN
                bL1[ds,b] = NaN
                bPh[ds,b] = NaN
                bTh[ds,b] = NaN
                bbe[ds,b] = NaN
                Lambda_1_rmse[ds,b] = NaN
                Phi_rmse[ds,b] = NaN
                Thet_rmse[ds,b] = NaN
                betas_rmse[ds,b] = NaN
            end
        end
        # println("Design ",d,"; ",ok," converged")
    end
    BACC = replace(BACC, NaN => missing); BACC_mean = mean(reshape(collect(skipmissing(BACC)),D,B),dims=2)[:,1]; BACC_sd = std(reshape(collect(skipmissing(BACC)),D,B),corrected=true,dims=2)[:,1]
    MCC = replace(MCC, NaN => missing); MCC_mean = mean(reshape(collect(skipmissing(MCC)),D,B),dims=2)[:,1]; MCC_sd = std(reshape(collect(skipmissing(MCC)),D,B),corrected=true,dims=2)[:,1]
    SE = replace(SE, NaN => missing); SE_mean = mean(reshape(collect(skipmissing(SE)),D,B),dims=2)[:,1]; SE_sd = std(reshape(collect(skipmissing(SE)),D,B),corrected=true,dims=2)[:,1]
    SP = replace(SP, NaN => missing); SP_mean = mean(reshape(collect(skipmissing(SP)),D,B),dims=2)[:,1]; SP_sd = std(reshape(collect(skipmissing(SP)),D,B),corrected=true,dims=2)[:,1]
    ARI = replace(ARI, NaN => missing); ARI_mean = mean(reshape(collect(skipmissing(ARI)),D,B),dims=2)[:,1]; ARI_sd = std(reshape(collect(skipmissing(ARI)),D,B),corrected=true,dims=2)[:,1]
    # bL1 = replace(bL1, NaN => missing); bL1_mean = mean(reshape(collect(skipmissing(bL1)),D,B),dims=2)[:,1]; bL1_sd = std(reshape(collect(skipmissing(bL1)),D,B),corrected=true,dims=2)[:,1]
    # bPh = replace(bPh, NaN => missing); bPh_mean = mean(reshape(collect(skipmissing(bPh)),D,B),dims=2)[:,1]; bPh_sd = std(reshape(collect(skipmissing(bPh)),D,B),corrected=true,dims=2)[:,1]
    # bTh = replace(bTh, NaN => missing); bTh_mean = mean(reshape(collect(skipmissing(bTh)),D,B),dims=2)[:,1]; bTh_sd = std(reshape(collect(skipmissing(bTh)),D,B),corrected=true,dims=2)[:,1]
    # bbe = replace(bbe, NaN => missing); bbe_mean = mean(reshape(collect(skipmissing(bbe)),D,B),dims=2)[:,1]; bbe_sd = std(reshape(collect(skipmissing(bbe)),D,B),corrected=true,dims=2)[:,1]
    # Lambda_1_rmse = replace(Lambda_1_rmse, NaN => missing); Lambda_1_rmse_mean = mean(reshape(collect(skipmissing(Lambda_1_rmse)),D,B),dims=2)[:,1]; Lambda_1_rmse_sd = std(reshape(collect(skipmissing(Lambda_1_rmse)),D,B),corrected=true,dims=2)[:,1]
    # Phi_rmse = replace(Phi_rmse, NaN => missing); Phi_rmse_mean = mean(reshape(collect(skipmissing(Phi_rmse)),D,B),dims=2)[:,1]; Phi_rmse_sd = std(reshape(collect(skipmissing(Phi_rmse)),D,B),corrected=true,dims=2)[:,1]
    # Thet_rmse = replace(Thet_rmse, NaN => missing); Thet_rmse_mean = mean(reshape(collect(skipmissing(Thet_rmse)),D,B),dims=2)[:,1]; Thet_rmse_sd = std(reshape(collect(skipmissing(Thet_rmse)),D,B),corrected=true,dims=2)[:,1]
    # betas_rmse = replace(betas_rmse, NaN => missing); betas_rmse_mean = mean(reshape(collect(skipmissing(betas_rmse)),D,B),dims=2)[:,1]; betas_rmse_sd = std(reshape(collect(skipmissing(betas_rmse)),D,B),corrected=true,dims=2)[:,1]

    jldsave("/home/niccolo.cao/storage/careless/resclass.jld2",res=hcat(BACC_mean,BACC_sd,MCC_mean,MCC_sd,SE_mean,SE_sd,SP_mean,SP_sd,ARI_mean,ARI_sd))
    # jldsave("/home/niccolo.cao/storage/careless/resbias.jld2",res=hcat(bL1_mean,bL1_sd,bPh_mean,bPh_sd,bTh_mean,bTh_sd,bbe_mean,bbe_sd))
    # jldsave("/home/niccolo.cao/storage/careless/resrmse.jld2",res=hcat(Lambda_1_rmse_mean,Lambda_1_rmse_sd,Phi_rmse_mean,Phi_rmse_sd,Thet_rmse_mean,Thet_rmse_sd,betas_rmse_mean,betas_rmse_sd))
end




function measuresp2(design,B)
   
    D = size(design,1)
    Lambda_1_rmse = Array{Float64}(undef,D,B); Theta_delta_rmse = Array{Float64}(undef,D,B); Phi_rmse = Array{Float64}(undef,D,B); kappa_rmse = Array{Float64}(undef,D,B); 
    Lambda_1_bias = Array{Float64}(undef,D,B); Theta_delta_bias = Array{Float64}(undef,D,B); Phi_bias = Array{Float64}(undef,D,B); kappa_bias = Array{Float64}(undef,D,B); 
    BACC = Array{Float64}(undef,D,B);  SE = Array{Float64}(undef,D,B);  SP = Array{Float64}(undef,D,B);  MCC = Array{Float64}(undef,D,B); ARI = Array{Float64}(undef,D,B);
    bL1 = Array{Float64}(undef,D,B); bPh = Array{Float64}(undef,D,B); bTh = Array{Float64}(undef,D,B); bbe = Array{Float64}(undef,D,B); 
    Lambda_1_rmse = Array{Float64}(undef,D,B); Phi_rmse = Array{Float64}(undef,D,B); Thet_rmse = Array{Float64}(undef,D,B); betas_rmse = Array{Float64}(undef,D,B);  
    seed = vcat(322011611885,11955714332337,14889268790486,19321661366842, 4539386885746,17304354949517, 1981879088971,14248930143082,16766904013492, 1434926186704,22695264190645,17517986428637,17678641845771, 1581593402381, 3465694164649, 9106731443143, 9729393064462,12160146297624, 4420744722171,16184395927257,19766499429427,20201406362254,19233052068296, 4730430435062)
    

    ds = 0
    for d in 1:size(design,1)
        ds += 1
        n = 1000
        p = 9
        q = 3
        c = design[d,1]
        K = design[d,2]
        propo = design[d,3]
        ncov = 2
            
        A = reshape(Diagonal(ones(q)),q,q);
        L_str = kron(A,ones(Int64(p/q)));
        
        Random.seed!(seed[d])
        ok = 0
        for b in 1:B
            res = JLD2.load(string("/home/nc/Documents/careless_res_sim/res/res_",d,"_",b,".jld2"))["res"]#X = JLD2.load(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"))["data"]
            Lambda_1 = sum(reshape(res[1:(p*q)],p,q),dims=2)
            Theta_delta = res[((p*q)+(p*K)+1):((p*q)+(p*K)+p)]
            Phi = reshape(res[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)],q,q)
            Phi = Phi[tril(trues(size(Phi)),-1)]
            kappa = res[((p*q)+(p*K)+p+p+q^2+1):((p*q)+(p*K)+p+p+q^2+n)]
            Ez1 = round.(Int64,res[((p*q)+(p*K)+p+p+q^2+n+1+1):((p*q)+(p*K)+p+p+q^2+n+1+n)])
            betas = res[((p*q)+(p*K)+p+p+q^2+n+1+n+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov)]

            Lambda0 = [0.69, 0.43, 0.7, 0.36, 0.78, 0.80, 0.56, 0.42, 0.67];
            Phi0vec = [0.12, 0.13, -0.21]; #Diagonal(ones(q,q))
            Phi0 = [1 0.12 0.13;0.12 1 -0.21;.13 -0.21 1]
            Theta_delta0 = diag(Diagonal(1 .-Diagonal(((Lambda0.*L_str)*Phi0*(Lambda0.*L_str)'))));
            # if propo == 0.85
            #     betas0 = vcat(log(0.85/(1 - .85)),0.8)
            # elseif propo == 0.95
            #     betas0 = vcat(log(0.95/(1 - .95)),0.8)
            # end
            z = JLD2.load(string("/home/nc/Documents/careless_res_sim/data/z_",d,"_",b,".jld2"))["z"]

            
            if sum(Ez1)/n<0.5
                Ez1 = 1 .- Ez1
                # betas = .- betas
                # println("label unswitched: ",betas)
            end
            
            if Ez1[1] !== NaN
                ok += 1
                cls = classperf(z,Ez1)
                # println(cls)
                BACC[ds,b] = Float64(cls[2,1])
                MCC[ds,b] = Float64(cls[2,2])
                SE[ds,b] = Float64(cls[2,3])
                SP[ds,b] = Float64(cls[2,4])
                ARI[ds,b] = Float64(cls[2,5])
                bL1[ds,b] = mean(Lambda_1.-Lambda0)
                bPh[ds,b] = mean(Phi.-Phi0vec)
                bTh[ds,b] = mean(Theta_delta.-Theta_delta0)
                bbe[ds,b] = mean(betas.-betas0)
                Lambda_1_rmse[ds,b] = sqrt(mean((Lambda_1 .- Lambda0)*(Lambda_1.-Lambda0)'))
                Phi_rmse[ds,b] = sqrt(mean((Phi .- Phi0vec)*(Phi.-Phi0vec)'))
                Thet_rmse[ds,b] = sqrt(mean((Theta_delta .- Theta_delta0)*(Theta_delta .- Theta_delta0)'))
                # betas_rmse[ds,b] = sqrt(mean((betas .- betas0)*(betas.-betas0)'))
                # println(Lambda_1,Lambda0," ",Phi,Phi0vec," ",Theta_delta,Theta_delta0," ",betas,betas0)
            else
                println("uno")
                BACC[ds,b] = NaN
                MCC[ds,b] = NaN
                SE[ds,b] = NaN
                SP[ds,b] = NaN
                ARI[ds,b] = NaN
                bL1[ds,b] = NaN
                bPh[ds,b] = NaN
                bTh[ds,b] = NaN
                bbe[ds,b] = NaN
                Lambda_1_rmse[ds,b] = NaN
                Phi_rmse[ds,b] = NaN
                Thet_rmse[ds,b] = NaN
                # betas_rmse[ds,b] = NaN
            end
        
        end

    end
    writedlm(string("/home/nc/MEGA/CFA+EFA/res_simulazioni/distrclassCareless/BACC.csv"),BACC)
    writedlm(string("/home/nc/MEGA/CFA+EFA/res_simulazioni/distrclassCareless/SE.csv"),SE)
    writedlm(string("/home/nc/MEGA/CFA+EFA/res_simulazioni/distrclassCareless/SP.csv"),SP)
    writedlm(string("/home/nc/MEGA/CFA+EFA/res_simulazioni/distrclassCareless/MCC.csv"),MCC)
    writedlm(string("/home/nc/MEGA/CFA+EFA/res_simulazioni/distrclassCareless/ARI.csv"),ARI)


end







################################################################################################ Model ################################################################################################

#Y: data matrix in items x subjects form (p x n)
#q: CFA latent factors
#K: EFA latent factors
#maxIter: maximum number of rerun of the EM algorithm in case of non-convergence
#mu_hat: optional argument with 1 for estimation of CFA latent means, 0 for setting to zero the CFA means
#L_str: "Lambda structure" matrix of the same size of CFA factor loadings matrix (p x q), which indicates the items loading to a CFA latent factors by specifying ones in case of item-factor loading and zeros in absence of loading.
#       The matrix is needed if the items in Y are not partitioned in subsets loading to a single factor or if the number of items loading to different factors is heterogeneous. It automatically provides a Lambda structure matrix with no cross-loadings.


function CFAmixEFA(Y::Matrix{Float64}, q::Int64,  K::Int64, maxIter::Int64, X::Matrix{Float64}, maxNumRestart::Int64, mu_hat::Int64=1,L_str=nothing)
    local re::Int64, cal::Int64, non_conv::Bool, err2::Bool, sum_diff_tf::Bool, retr::Int64, res::Vector{Float64}, start_vals::Vector{Float64}, DomainError_warn::Int64, p::Int64, n::Int64
   
    
    ncov = size(X,2); p = size(Y,1); n = size(Y,2); iter_err = 1; code_error = 0.0;  err_0_count = 0; err_1_count = 0; err_2_count = 0; err_3_count = 0; err_4_count = 0; err_5_count = 0;
    while code_error != 8.0
        cal = (p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n
        
        if iter_err >= maxNumRestart
            code_error = 7.0 
            return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n),code_error, err_0_count, err_1_count, err_2_count, err_3_count, err_4_count, err_5_count) 
        end
        
        # Starting values
        start_vals = starting_val(n,p,q,K,X,mu_hat,L_str)
         
        # EM 
        res = EM(Y,X,q,K,p,n,ncov,maxIter,reshape(start_vals[1:(p*q)],p,q), reshape(start_vals[(p*q+1):((p*q)+(p*K))],p,K), Diagonal(start_vals[((p*q)+(p*K)+1):((p*q)+(p*K)+p)]), Diagonal(start_vals[((p*q)+(p*K)+p+1):((p*q)+(p*K)+p+p)]), reshape(start_vals[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)],q,q), vec(start_vals[((p*q)+(p*K)+p+p+q^2+1):((p*q)+(p*K)+p+p+q^2+n)]), vec(start_vals[((p*q)+(p*K)+p+p+q^2+n+1):((p*q)+(p*K)+p+p+q^2+n+ncov)]), reshape(start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q)],p,q),start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q)],start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+K)],mu_hat)
  
        

        # Starting values
        # repli = 10
        # ress = Matrix{Float64}(undef,cal,repli); liks = Vector{Float64}(undef,repli); 
        # for i in 1:repli
            
        #     start_vals = starting_val(n,p,q,K,X,mu_hat,L_str)
        #     ress[:,i] = EM(Y,X,q,K,p,n,ncov,maxIter,reshape(start_vals[1:(p*q)],p,q), reshape(start_vals[(p*q+1):((p*q)+(p*K))],p,K), Diagonal(start_vals[((p*q)+(p*K)+1):((p*q)+(p*K)+p)]), Diagonal(start_vals[((p*q)+(p*K)+p+1):((p*q)+(p*K)+p+p)]), reshape(start_vals[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)],q,q), vec(start_vals[((p*q)+(p*K)+p+p+q^2+1):((p*q)+(p*K)+p+p+q^2+n)]), vec(start_vals[((p*q)+(p*K)+p+p+q^2+n+1):((p*q)+(p*K)+p+p+q^2+n+ncov)]), reshape(start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q)],p,q),start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q)],start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+K)],mu_hat)[1:cal]
        #     liks[i] = ress[(p*q)+(p*K)+p+p+q^2+n+1,i]
        # end
        # println(liks[Bool.(1 .-(isnan.(liks)))])
        # liktr = liks[Bool.(1 .-(isnan.(liks)))]

        # if sum(isnan.(liks)) != repli
        #     rex = ress[:,Bool.(1 .-(isnan.(liks)))][:,findall((liktr).==maximum(liktr))]
        #     # EM 
        #     res = EM(Y,X,q,K,p,n,ncov,maxIter,reshape(rex[1:(p*q)],p,q), reshape(rex[(p*q+1):((p*q)+(p*K))],p,K), Diagonal(rex[((p*q)+(p*K)+1):((p*q)+(p*K)+p)]), Diagonal(rex[((p*q)+(p*K)+p+1):((p*q)+(p*K)+p+p)]), reshape(rex[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)],q,q), vec(rex[((p*q)+(p*K)+p+p+q^2+1):((p*q)+(p*K)+p+p+q^2+n)]), vec(rex[((p*q)+(p*K)+p+p+q^2+n+1):((p*q)+(p*K)+p+p+q^2+n+ncov)]), reshape(rex[((p*q)+(p*K)+p+p+q^2+n+ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q)],p,q),rex[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q)],rex[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+K)],mu_hat)
        # elseif sum(isnan.(liks)) == repli
        #     # EM 
        #     start_vals = starting_val(n,p,q,K,X,mu_hat,L_str)
        #     res = EM(Y,X,q,K,p,n,ncov,maxIter,reshape(start_vals[1:(p*q)],p,q), reshape(start_vals[(p*q+1):((p*q)+(p*K))],p,K), Diagonal(start_vals[((p*q)+(p*K)+1):((p*q)+(p*K)+p)]), Diagonal(start_vals[((p*q)+(p*K)+p+1):((p*q)+(p*K)+p+p)]), reshape(start_vals[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)],q,q), vec(start_vals[((p*q)+(p*K)+p+p+q^2+1):((p*q)+(p*K)+p+p+q^2+n)]), vec(start_vals[((p*q)+(p*K)+p+p+q^2+n+1):((p*q)+(p*K)+p+p+q^2+n+ncov)]), reshape(start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q)],p,q),start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q)],start_vals[((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+1):((p*q)+(p*K)+p+p+q^2+n+ncov+p*q+n*ncov+q+K)],mu_hat)
        # end
        
        # Controls for non covnergence
        cal = (p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n+1
        sum_diff_tf = !(res[(cal+1):(cal+1)][1]==1.0)
        non_conv = res[(cal+1+1):(cal+1+1)][1]>=maxIter
        code_error = res[(cal):(cal)][1]
        if code_error == 0.0  && iter_err < maxNumRestart
            err_0_count += 1 
            iter_err += 1
        elseif code_error == 1.0  && iter_err < maxNumRestart
            err_1_count += 1 
            iter_err += 1
        elseif code_error == 2.0  && iter_err < maxNumRestart
            err_2_count += 1 
            iter_err += 1
        elseif code_error == 3.0 && iter_err < maxNumRestart
            err_3_count += 1 
            iter_err += 1
        elseif non_conv && iter_err < maxNumRestart
            err_4_count += 1 
            iter_err += 1
        elseif sum_diff_tf && iter_err < maxNumRestart
            err_5_count += 1 
            iter_err += 1
        end
    end 
    return res[1:(cal-1)]
end

# Expectation-Maximization algorithm
function EM(Y::Matrix{Float64},X::Matrix{Float64},q::Int64, K::Int64, p::Int64, n::Int64, ncov::Int64, maxIter::Int64, Lambda_1::Matrix{Float64}, Lambda_2::Matrix{Float64}, Theta_delta::Diagonal{Float64, Vector{Float64}}, Psi_delta::Diagonal{Float64, Vector{Float64}}, Phi::Matrix{Float64}, kappa::Vector{Float64}, betas::Vector, L_str::Matrix, mu::Vector{Float64}, nu::Vector{Float64}, mu_hat::Int64)
    local llik::Vector{Float64}, Ez1::Vector{Float64}, diff_llik::Float64, oldllik::Float64, dtol::Float64, t::Int64, err::Vector{Bool}, err2::Bool, sum_diff_tf::Bool, sum_diff::Vector{Bool}, old_diff::Float64,  ESxixi::Matrix{Float64}, ESetaeta::Matrix{Float64}, ESyxi::Matrix{Float64}, ESyeta::Matrix{Float64}, ESyy0::Matrix{Float64}, ESyy1::Matrix{Float64}, Ek0::Float64, Ek1::Float64, Ez0::Vector{Float64}, B2::Matrix{Float64}, b2::Matrix{Float64}, B1::Matrix{Float64}, b1::Matrix{Float64}, Exi4::Matrix{Float64}, Eeta4::Matrix{Float64}
    
    zero_indices = [i + (j - 1) * p for j in 1:K for i in 1:(j-1)]
            
    t = 1; err2 = false; llik = zeros(Float64, maxIter); oldllik = -Inf; err=zeros(Bool,maxIter); sum_diff = zeros(Bool,maxIter); sum_diff_tf = true;  diff_llik = 1.0; dtol = 0.001; code_error = 8.0;
    while diff_llik>dtol && t<=maxIter && sum_diff_tf
            inv_lemma_1 = inv_lemma_1_fun(q,Lambda_1,Theta_delta,Phi)
            inv_lemma_2 = inv_lemma_2_fun(K,Lambda_2,Psi_delta)      
            det_lemma_1 = det_lemma_1_fun(Lambda_1,Theta_delta,Phi)
            det_lemma_2 = det_lemma_2_fun(K,Lambda_2,Psi_delta)

            if det_lemma_1<=0 || det_lemma_2<=0
                code_error = 0.0
                return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n), code_error, sum_diff_tf, t-1) 
            end
            # println(t)
            #E-step
            # b1 = b1_fun(q,Lambda_1,Phi,Theta_delta,mu)
            # B1 = B1_fun(Lambda_1,Phi,Theta_delta,b1,mu)
            # b2 = b2_fun(K,Lambda_2,Psi_delta,nu)
            # B2 = B2_fun(K,Lambda_2,Psi_delta,b2,nu)
            
            b1 = b1_fun(q,Lambda_1,Phi,Theta_delta)
            B1 = B1_fun(Lambda_1,Phi,Theta_delta,b1)
            b2 = b2_fun(K,Lambda_2,Psi_delta)
            B2 = B2_fun(K,Lambda_2,Psi_delta,b2)
            
            Ez1 = Expectz(n,p,Y,kappa,inv_lemma_1,inv_lemma_2, det_lemma_1,det_lemma_2,Lambda_1,Lambda_2,mu,nu)

            Ez0 = Expectz0(Ez1)
            Ek1 = Ekz1(Ez1)
            Ek0 = Ekz0(Ez0)
            
            ESyy1 = ESyyz1(n,p,Y,Ez1)  
            ESyy0 = ESyyz0(n,p,Y,Ez0)
            
            Ezy1 = ESzy1(n,p,Y,Ez1,mu,Lambda_1)
            Ezy0 = ESzy0(n,p,Y,Ez0,nu,Lambda_2)
            Ey1 = ESyone(n,p,Y,Ez1)
            Ey0 = ESyzero(n,p,Y,Ez0)
            
            Eeta = Eeta1(p,n,Ek1,Ez1,Y,b1, mu, Lambda_1)
            Exi = Exi0(p,n,Ek0,Ez0,Y,b2, nu, Lambda_2)
            
            Eetay1 =  ESetay1(n,p,Y,Ez1,Lambda_1,mu,b1,Ey1)
            Exiy0 =  ESxiy0(n,p,Y,Ez0,Lambda_2,nu,b2,Ey0)

            ESetaeta =  ESetaeta1(n,p,Y,Ez1,Lambda_1,mu,b1,B1)
            ESxixi = ESxixi0(n,p,Y,Ez0,Lambda_2,nu,b2,B2)
         
            if isequal(sum(ESetaeta),NaN) || isequal(sum(Eetay1),NaN) || isnan(ESxixi[1,1]) || isequal(sum(ESxixi),0)
                code_error = 0.0
                return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n), code_error, sum_diff_tf, t-1) 
            else
                #M-step
                Lambda_1 = Lambda_1_mstep(Eetay1,ESetaeta,L_str)
                Lambda_2 = Lambda_2_mstep(Exiy0,ESxixi)
                L2 = vec(Lambda_2)
                L2[zero_indices,:].=0.
                Lambda_2 = reshape(L2,p,K)

                Phi = Phi_mstep(ESetaeta,Eeta,Ek1,mu,n)

                Theta_delta =  Diagonal(1 .-Diagonal((Lambda_1*Phi*Lambda_1')));# Theta_delta_mstep(n,Lambda_1,ESyy1,Eetay1,Ek1)
                Psi_delta = Psi_delta_mstep(n,Lambda_2,ESyy0,Exiy0,Ek0)#Psi_delta_mstep(Y,b2,B2,nu,n,p,K,Lambda_2,ESyy0,Exiy0,Ez0)
                

                if mu_hat == 1
                    mu = ((sum(Eeta,dims=2)/n)/(Ek1/n))[1:q,1]
                elseif mu_hat == 0
                    mu = zeros(q)
                end
                nu = ((sum(Exi,dims=2)/n)/(Ek0/n))[1:K,1]
                
            
            
                #If no covariates are provided
                if ncov == 1
                    kappa_1 = Ek1/n
                    kappa = fill(kappa_1,n)
                    betas = log(kappa_1/(1 - kappa_1))
                else
                    betas = NewtonRaphson(betas_gradient, betas, X, Ez1, Ez0, ncov)
                    kappa = exp.(X*betas)./(exp.(X*betas).+1)
                end         
            
                inv_lemma_1 = inv_lemma_1_fun(q,Lambda_1,Theta_delta,Phi)
                inv_lemma_2 = inv_lemma_2_fun(K,Lambda_2,Psi_delta)      
                det_lemma_1 = det_lemma_1_fun(Lambda_1,Theta_delta,Phi)
                det_lemma_2 = det_lemma_2_fun(K,Lambda_2,Psi_delta)
                #Controls for non-convergence
                if sum(Theta_delta.<0) > 0 || det(Psi_delta) < 0 || isequal(sum(Theta_delta),NaN) || isequal(sum(Ez0),NaN) || sum(diag(inv_lemma_1).<0) > 0 || sum(diag(inv_lemma_2).<0) > 0 || det_lemma_1 < 0 || det_lemma_2<0 
                    code_error = 0.0
                    return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n), code_error, sum_diff_tf, t-1) 
                elseif isnan(betas[1]) 
                    code_error = 1.0
                    return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n), code_error, sum_diff_tf, t-1) 
                elseif det(Phi)<0 
                    code_error = 2.0
                    return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n), code_error, sum_diff_tf, t-1) 
                else
                    old_diff = diff_llik
                    llik[t] = llik_compute(Y,n,p,q,K,Lambda_1, Lambda_2, kappa,mu,nu,inv_lemma_1,inv_lemma_2,det_lemma_1,det_lemma_2)

                    diff_llik = abs(oldllik - llik[t])
                    if !(oldllik<llik[t])
                        code_error = 3.0
                        return vcat(fill(NaN,p*q+p*K+p+p+q^2+n+1+n+ncov+q+K+K*n+q*n), code_error, sum_diff_tf, t-1) 
                    end

                    sum_diff[t] = diff_llik>old_diff               
                    oldllik = llik[t]
                    t=t+1;
                    sum_diff_tf = sum(sum_diff) < 150            
                end
            end
        end
        Exitrue=Exi3(p,n,Y,b2, nu, Lambda_2)
        Eetatrue=Eeta3(p,n,Y,b1, mu, Lambda_1)
    return vcat(vec(Lambda_1), vec(Lambda_2), vec(diag(Theta_delta)), vec(diag(Psi_delta)), vec(Phi), kappa, llik[t-1], Ez1, betas, mu, nu, vec(Exitrue), vec(Eetatrue), code_error, sum_diff_tf, t-1);
end        

#Beta gradient
function betas_gradient(betas::Vector,X::Matrix{Float64},Ez1::Vector{Float64},Ez0::Vector{Float64})
    fval = - X'*((-Ez1.+Ez0.*exp.(X*betas))./(exp.(X*betas).+1))#( sum((-Ez1.+Ez0.*exp.(X*betas)).*X./(exp.(X*betas).+1), dims=1)) #(.- sum(Ez1.*X./(exp.(X*betas).+1), dims=1) .+ sum(Ez0.*X.*exp.(X*betas)./(exp.(X*betas).+1),dims=1))
    return vec(fval)
end




function NewtonRaphson(betas_gradient, betas::Vector, X::Matrix{Float64}, Ez1::Vector{Float64}, Ez0::Vector{Float64}, ncov::Int64)
    local betas_old::Vector{Float64}, iter::Int64, tol::Float64, maxIter::Int64, grad::Vector{Float64}, hess::Matrix{Float64}
    betas_old = repeat([-Inf],size(X,2)); iter = 0;  tol=1e-8; maxIter = 1000
    while sum(abs.(betas.-betas_old)) > tol && iter < maxIter
        grad = betas_gradient(betas,X,Ez1,Ez0)
        hess = -(X.*exp.(X*betas)./(exp.(X*betas).+1).^2)'*X
        # println(hess)
        if det(hess) == 0.0 || isnan(sum(hess)) || isinf(sum(hess))
            return vec(fill(NaN,ncov))
        else
        betas_old = betas
        betas = betas - hess\grad
        iter += 1
        end
        if iter >= maxIter
            betas = vec(fill(NaN,ncov))
        end
    end
    return betas
end

# Starting values
function starting_val(n::Int64, p::Int64, q::Int64, K::Int64, X::Matrix, mu_hat::Int64,L_str=nothing)
    local Lambda_1::Matrix{Float64}, Lambda_2::Matrix{Float64}, Theta_delta::Diagonal{Float64, Vector{Float64}}, Psi_delta::Diagonal{Float64, Vector{Float64}}, Phi::Matrix{Float64}, kappa::Vector{Float64}

    #CFA
    if L_str === nothing
        L_str = L_str_fun(p,q)
    end
    Lambda_1 = rand(Uniform(0,1),p,q) .* L_str
    Phi = rand(LKJ(q,1.0))
    Theta_delta = Diagonal(1 .-Diagonal((Lambda_1*Phi*Lambda_1')));
    if mu_hat == 1
        mu = rand(Uniform(-2,2),q) 
    elseif mu_hat == 0
        mu = zeros(q)
    end
        
    #EFA
    Lambda_2 = rand(Uniform(0,1),p,K)
    Psi_delta = Diagonal(rand(Uniform(0.05,2),p))
    nu = rand(Uniform(-2,4),K) 

    #Mixture    
    betas = rand(Uniform(-1.5,1.5), size(X,2))
    kappa = Vector{Float64}(undef,n)
    if size(X,2) == 1
        kappa = vec(exp.(X*betas)./(exp.(X*betas).+1))
    else
        kappa = exp.(X*betas)./(exp.(X*betas).+1)
    end

    return vcat(vec(Lambda_1), vec(Lambda_2), vec(diag(Theta_delta)), vec(diag(Psi_delta)), vec(Phi), kappa, betas, vec(L_str), vec(X), mu, nu)
end

# Construction of Lambda structure matrix
function L_str_fun(p::Int64,q::Int64)
    local A::Matrix{Int64}, str::Vector{Int64}, L_str::Matrix{Float64}
    A =  Matrix{Int64}(I,q,q)
    str = ones(Int64,Int64(p/q))
    L_str = kron(A, str);
    return L_str
end

# More convenient matrix calculus
function inv_lemma_1_fun(q::Int64,Lambda_1::Matrix{Float64},Theta_delta::Diagonal{Float64, Vector{Float64}},Phi::Matrix{Float64})
    local inv_lemma::Matrix{Float64}
    inv_lemma = (inv(Theta_delta)-(Theta_delta\Lambda_1)*((Matrix{Float64}(I,q,q)+Phi*Lambda_1'*(Theta_delta\Lambda_1))\Phi)*Lambda_1'*inv(Theta_delta))
    return inv_lemma
end

function inv_lemma_2_fun(K::Int64,Lambda_2::Matrix{Float64},Psi_delta::Diagonal{Float64, Vector{Float64}})
    local inv_lemma::Matrix{Float64}
    inv_lemma = (inv(Psi_delta)-(Psi_delta\Lambda_2)*((Matrix{Float64}(I,K,K)+Lambda_2'*(Psi_delta\Lambda_2))\Lambda_2')*inv(Psi_delta))
    return inv_lemma
end

function det_lemma_1_fun(Lambda_1::Matrix{Float64},Theta_delta::Diagonal{Float64, Vector{Float64}},Phi::Matrix{Float64})
    local det_lemma::Float64
    det_lemma = det(inv(Phi)+Lambda_1'*(Theta_delta\Lambda_1))*det(Phi)*det(Theta_delta)
    return det_lemma
end

function det_lemma_2_fun(K::Int64,Lambda_2::Matrix{Float64},Psi_delta::Diagonal{Float64, Vector{Float64}})
    local det_lemma::Float64
    det_lemma = det(Matrix{Float64}(I,K,K)+Lambda_2'*(Psi_delta\Lambda_2))*det(Psi_delta)
    return det_lemma
end



# Functions for E-step

function b1_fun(q::Int64,Lambda_1::Matrix{Float64},Phi::Matrix{Float64},Theta_delta::Diagonal{Float64, Vector{Float64}})
    local b1::Matrix{Float64}
    b1 = Phi*Lambda_1'*inv(Theta_delta)-Phi*Lambda_1'*(Theta_delta\Lambda_1)*((Matrix{Float64}(I,q,q)+Phi*Lambda_1'*(Theta_delta\Lambda_1))\Phi)*Lambda_1'*inv(Theta_delta)
    return b1
end

function B1_fun(Lambda_1::Matrix{Float64},Phi::Matrix{Float64},Theta_delta::Diagonal{Float64, Vector{Float64}}, b1::Matrix{Float64})
    local B1::Matrix{Float64}
    B1 = Phi-b1*Lambda_1*Phi
    return B1
end

function b2_fun(K::Int64,Lambda_2::Matrix{Float64},Psi_delta::Diagonal{Float64, Vector{Float64}})
    local b2::Matrix{Float64}
    b2 =  Lambda_2'*inv(Psi_delta)-Lambda_2'*(Psi_delta\Lambda_2)*((Matrix{Float64}(I,K,K)+Lambda_2'*(Psi_delta\Lambda_2))\Lambda_2')*inv(Psi_delta)
    return b2
end

function B2_fun(K::Int64,Lambda_2::Matrix{Float64},Psi_delta::Diagonal{Float64, Vector{Float64}}, b2::Matrix{Float64})
    local B2::Matrix{Float64}
    B2 = Matrix{Float64}(I,K,K)-b2*Lambda_2
    return B2
end


# function b1_fun(q::Int64,Lambda_1::Matrix{Float64},Phi::Matrix{Float64},Theta_delta::Diagonal{Float64, Vector{Float64}},mu::Vector{Float64})
#     local b1::Matrix{Float64}
#     b1 = (Phi*Lambda_1'-mu*mu'*Lambda_1')*inv(Lambda_1*Phi*Lambda_1'+Theta_delta-Lambda_1*mu*mu'*Lambda_1')
#     return b1
# end

# function B1_fun(Lambda_1::Matrix{Float64},Phi::Matrix{Float64},Theta_delta::Diagonal{Float64, Vector{Float64}}, b1::Matrix{Float64},mu::Vector{Float64})
#     local B1::Matrix{Float64}
#     B1 = Phi-mu*mu'-b1*(Lambda_1*Phi-Lambda_1*mu*mu')
#     return B1
# end

# function b2_fun(K::Int64,Lambda_2::Matrix{Float64},Psi_delta::Diagonal{Float64, Vector{Float64}},nu::Vector{Float64})
#     local b2::Matrix{Float64}
#     b2 = (Lambda_2'-nu*nu'*Lambda_2')*inv(Lambda_2*Lambda_2'+Psi_delta-Lambda_2*nu*nu'*Lambda_2')
#     return b2
# end

# function B2_fun(K::Int64,Lambda_2::Matrix{Float64},Psi_delta::Diagonal{Float64, Vector{Float64}}, b2::Matrix{Float64},nu::Vector{Float64})
#     local B2::Matrix{Float64}
#     B2 = Matrix{Float64}(I,K,K)-nu*nu'-b2*(Lambda_2-Lambda_2*nu*nu')
#     return B2
# end


function Expectz(n::Int64, p::Int64, Y::Matrix{Float64},kappa::Vector{Float64},inv_lemma_1::Matrix{Float64},inv_lemma_2::Matrix{Float64}, det_lemma_1::Float64,det_lemma_2::Float64, Lambda_1::Matrix{Float64}, Lambda_2::Matrix{Float64}, mu::Vector{Float64}, nu::Vector{Float64})
    local Ez1::Vector{Float64}
    
    Ez1 = zeros(Float64,n)
    @inbounds  for i in 1:n
        Ez1[i] = ((kappa[i] * ((exp((-1/2) * (Y[:,i]-Lambda_1*mu)'*inv_lemma_1*(Y[:,i]-Lambda_1*mu))/(sqrt((2 * pi)^(p) *det_lemma_1)))))/
        ((kappa[i] * ((exp((-1/2) * (Y[:,i]-Lambda_1*mu)'*inv_lemma_1*(Y[:,i]-Lambda_1*mu))/(sqrt((2 * pi)^(p) *det_lemma_1)))) 
        +((exp((-1/2) * (Y[:,i]-Lambda_2*nu)'*inv_lemma_2*(Y[:,i]-Lambda_2*nu)))/(sqrt((2 * pi)^(p) *det_lemma_2))) 
        -(kappa[i] *((exp((-1/2) * (Y[:,i]-Lambda_2*nu)'*inv_lemma_2*(Y[:,i]-Lambda_2*nu)))/(sqrt((2 * pi)^(p) *det_lemma_2))))))) 
    end
    return Ez1
end

function Expectz0(Ez1::Vector{Float64})
    local Ez0::Vector{Float64}
    Ez0 = @. 1 - Ez1
    return Ez0
end

function Ekz1(Ez1::Vector{Float64})
    local Ek1::Float64
    Ek1 = sum(Ez1)
    return Ek1
end

function Ekz0(Ez0::Vector{Float64})
    local Ek0::Float64
    Ek0 = sum(Ez0)
    return Ek0
end

function ESyyz1(n::Int64,p::Int64,Y::Matrix{Float64},Ez1::Vector{Float64})
    local ESyy1::Matrix{Float64}

    ESyy1 = (Y.*Ez1')*Y'/n   
end


function ESyyz0(n::Int64,p::Int64,Y::Matrix{Float64},Ez0::Vector{Float64})
    local ESyy0::Matrix{Float64}

    ESyy0 = (Y.*Ez0')*Y'/n  
end

function ESzy1(n::Int64,p::Int64,Y::Matrix{Float64},Ez1::Vector{Float64},mu::Vector{Float64},Lambda_1::Matrix{Float64})
    local Ezy1::Matrix{Float64}

    Ezy1 = sum((Y.-Lambda_1*mu).*Ez1',dims=2)/n     
end

function ESzy0(n::Int64,p::Int64,Y::Matrix{Float64},Ez0::Vector{Float64},nu::Vector{Float64},Lambda_2::Matrix{Float64})
    local Ezy0::Matrix{Float64}

    Ezy0 = sum((Y.-Lambda_2*nu).*Ez0',dims=2)/n
end

function ESyone(n::Int64,p::Int64,Y::Matrix{Float64},Ez1::Vector{Float64})
    local Ey1::Matrix{Float64}

    Ey1 = sum(Y.*Ez1',dims=2)/n     
end

function ESyzero(n::Int64,p::Int64,Y::Matrix{Float64},Ez0::Vector{Float64})
    local Ey0::Matrix{Float64}

    Ey0 = sum(Y.*Ez0',dims=2)/n
end

function Eeta1(p::Int64,n::Int64,Ek1::Float64,Ez1::Vector{Float64},Y::Matrix{Float64},b1::Matrix{Float64}, mu::Vector{Float64}, Lambda_1::Matrix{Float64})
    local Eeta::Matrix{Float64}
   
    Eeta = b1*(Ez1'.*(Y.-Lambda_1*mu)) + mu*Ez1'   
    return Eeta
end

function Eeta3(p::Int64,n::Int64,Y::Matrix{Float64},b1::Matrix{Float64}, mu::Vector{Float64}, Lambda_1::Matrix{Float64})
    local Eeta4::Matrix{Float64}

    Eeta4 =  b1*(Y.-Lambda_1*mu) .+ mu
    return Eeta4
end

function Exi0(p::Int64,n::Int64,Ek0::Float64,Ez0::Vector{Float64},Y::Matrix{Float64},b2::Matrix{Float64}, nu::Vector{Float64}, Lambda_2::Matrix{Float64})
    local Exi::Matrix{Float64}

    Exi =  b2*(Ez0'.*(Y.-Lambda_2*nu)) + nu*Ez0' 
    return Exi
end

function Exi3(p::Int64,n::Int64,Y::Matrix{Float64},b2::Matrix{Float64}, nu::Vector{Float64}, Lambda_2::Matrix{Float64})
    local Exi4::Matrix{Float64}

    Exi4 =  b2*(Y.-Lambda_2*nu) .+ nu
    return Exi4
end

function ESetay1(n::Int64,p::Int64,Y::Matrix{Float64},Ez1::Vector{Float64},Lambda_1::Matrix{Float64},mu::Vector{Float64},b1::Matrix{Float64},Ey1::Matrix{Float64})
    local Eetay1::Matrix{Float64}

    Eetay1 = (b1*((Ez1'.*(Y.-Lambda_1*mu)) *Y'))'/n .+ Ey1*mu' 
end

function ESxiy0(n::Int64,p::Int64,Y::Matrix{Float64},Ez0::Vector{Float64},Lambda_2::Matrix{Float64},nu::Vector{Float64},b2::Matrix{Float64},Ey0::Matrix{Float64})
    local Exiy0::Matrix{Float64}

    Exiy0 = (b2*((Ez0'.*(Y.-Lambda_2*nu)) *Y'))'/n .+ Ey0*nu' 
end

function ESetaeta1(n::Int64,p::Int64,Y::Matrix{Float64},Ez1::Vector{Float64},Lambda_1::Matrix{Float64},mu::Vector{Float64},b1::Matrix{Float64},B1::Matrix{Float64})
    local ESetaeta::Matrix{Float64}

    ESetaeta =  B1*(sum(Ez1)/n) + (b1*(((Y.-Lambda_1*mu).*Ez1')*(Y.-Lambda_1*mu)')/n)*b1' + mu*mu'*(sum(Ez1)/n) + mu*sum((Ez1.*(Y.-Lambda_1*mu)'),dims=1)/n*b1'+b1*(sum((Ez1.*(Y.-Lambda_1*mu)'),dims=1)/n)'*mu'
end

function ESxixi0(n::Int64,p::Int64,Y::Matrix{Float64},Ez0::Vector{Float64},Lambda_2::Matrix{Float64},nu::Vector{Float64},b2::Matrix{Float64},B2::Matrix{Float64})
    local ESxixi::Matrix{Float64}

    ESxixi =  B2*(sum(Ez0)/n) + (b2*(((Y.-Lambda_2*nu).*Ez0')*(Y.-Lambda_2*nu)')/n)*b2' + nu*nu'*(sum(Ez0)/n) + nu*sum((Ez0.*(Y.-Lambda_2*nu)'),dims=1)/n*b2'+b2*(sum((Ez0.*(Y.-Lambda_2*nu)'),dims=1)/n)'*nu'
end

# (Matrix{Float64}(I,K,K)-nu*nu'+(Lambda_2'-nu*nu'*Lambda_2')*inv(Lambda_2*Lambda_2'+Psi_delta-Lambda_2*nu*nu'*Lambda_2')*(Lambda_2-Lambda_2*nu*nu')*(Ek0/n)+(nu.+(Lambda_2'-nu*nu'Lambda_2')*inv(Lambda_2*Lambda_2'+Psi_delta-Lambda_2*nu*nu'*Lambda_2')*(Y-Lambda_2*nu))*(nu.+(Lambda_2'-nu*nu'Lambda_2')*inv(Lambda_2*Lambda_2'+Psi_delta-Lambda_2*nu*nu'*Lambda_2')*(Y-Lambda_2*nu))')

# Functions for M-step

function Lambda_1_mstep(Eetay1::Matrix{Float64},ESetaeta::Matrix{Float64},L_str::Matrix)
    local Lambda_1::Matrix{Float64}
    Lambda_1 = Eetay1*inv(ESetaeta) .* L_str
    return Lambda_1
end

function Lambda_2_mstep(Exiy0::Matrix{Float64},ESxixi::Matrix{Float64})
    local Lambda_2::Matrix{Float64}
    Lambda_2 = Exiy0*inv(ESxixi) 
    return Lambda_2
end

function Theta_delta_mstep(n::Int64,Lambda_1::Matrix{Float64},ESyy1::Matrix{Float64},Eetay1::Matrix{Float64},Ek1::Float64)
    local Theta_delta::Diagonal{Float64, Vector{Float64}}
    Theta_delta = Diagonal((ESyy1 - Lambda_1*Eetay1')/(Ek1/n))
    return Theta_delta
end

function Psi_delta_mstep(n::Int64,Lambda_2::Matrix{Float64},ESyy0::Matrix{Float64},Exiy0::Matrix{Float64},Ek0::Float64)
    local Psi_delta::Diagonal{Float64, Vector{Float64}}
    Psi_delta = Diagonal((ESyy0 - Lambda_2*Exiy0')/(Ek0/n))
    return Psi_delta
end

# function Psi_delta_mstep(Y::Matrix{Float64},b2::Matrix{Float64},B2::Matrix{Float64},nu::Vector{Float64},n::Int64,p::Int64,K::Int64,Lambda_2::Matrix{Float64},ESyy0::Matrix{Float64},Exiy0::Matrix{Float64},Ez0::Vector{Float64})
#     local Psi_delta::Diagonal{Float64, Vector{Float64}}, sigma::Float64
#     sigma = 0.0
#     for i in 1:n
#         sigma = sigma + (Ez0[i].*norm((Y[:,i]))^2-2*(b2*(Ez0[i].*(Y[:,i].-Lambda_2*nu)))'*Lambda_2'*(Y[:,i]) +tr((B2*((Ez0[i])) + (b2*(((Y[:,i].-Lambda_2*nu).*Ez0[i]')*(Y[:,i].-Lambda_2*nu)'))*b2' + nu*nu'*((Ez0[i])) + nu*((Ez0[i].*(Y[:,i].-Lambda_2*nu)'))*b2'+b2*(((Ez0[i].*(Y[:,i].-Lambda_2*nu)')))'*nu')*Lambda_2'*Lambda_2))/(n*p)
#     end
#     Psi_delta = Diagonal(fill(sigma,p))
#     return Psi_delta
# end



function Phi_mstep(ESetaeta::Matrix{Float64},Eeta::Matrix{Float64},Ek1::Float64,mu::Vector{Float64},n::Int64)
    local Phi::Matrix{Float64}
    Phi = (ESetaeta-sum(Eeta,dims=2)*mu'/n)/(Ek1/n)
    Phi = sqrt(inv(Diagonal(Phi)))*Phi*sqrt(inv(Diagonal(Phi)))
    return Phi
end



# Log-likelihood computation

# function llik_compute(Y::Matrix{Float64},n::Int64,p::Int64,q::Int64,K::Int64,Lambda_1::Matrix{Float64}, Lambda_2::Matrix{Float64}, Theta_delta::Diagonal{Float64, Vector{Float64}}, Psi_delta::Diagonal{Float64, Vector{Float64}}, Phi::Matrix{Float64}, kappa::Vector{Float64}, Ez1::Vector{Float64}, Ez0::Vector{Float64},Ek0::Float64,Ek1::Float64,mu::Vector{Float64},nu::Vector{Float64},Eeta::Matrix{Float64},Exi::Matrix{Float64},Eetay1::Matrix{Float64},Exiy0::Matrix{Float64},ESyy1::Matrix{Float64},ESyy0::Matrix{Float64},ESetaeta::Matrix{Float64},ESxixi::Matrix{Float64},)
#     local llik::Float64, transpY::Adjoint{Float64, Matrix{Float64}}, tr1::Float64, tr2::Float64, sum1::Float64, tr3::Float64, sum2::Float64, tr4::Float64, tr5::Float64, tr6::Float64, sum3::Float64, tr7::Float64, sum4::Float64, tr8::Float64, B1_intra::Matrix{Float64}, B2_intra::Matrix{Float64}, b1_intra::Matrix{Float64}, b2_intra::Matrix{Float64}
#     transpY = Adjoint{Float64, Matrix{Float64}}(Y)
    
#     tr1 = tr(Theta_delta\ESyy1)*n
#     tr2 = 2*tr(Lambda_1'*(Theta_delta\Eetay1))*n
#     tr3 = tr(Lambda_1'*(Theta_delta\Lambda_1)*ESetaeta)*n
#     tr4 = tr(Phi\ESetaeta)*n
#     tr5 = 2*tr(inv(Phi)*sum(Eeta,dims=2)*mu')
#     tr6 = tr((Phi)\(Ek1*mu*mu'))
    
#     tr7 = tr(Psi_delta\ESyy0)*n
#     tr8 = 2*tr(Lambda_2'*(Psi_delta\Exiy0))*n
#     tr9 = tr(Lambda_2'*(Psi_delta\Lambda_2)*ESxixi)*n
#     tr10 = tr(ESxixi)*n
#     tr11 = 2*tr(sum(Exi,dims=2)*nu')
#     tr12 = tr((Ek0*nu*nu'))

#     lk1 = sum(Ez1.*log.(kappa))
#     lk0 = sum(Ez0.*log.(1 .-kappa))
#     lT = sum(Ez1)*(1/2)*log(det(Theta_delta))
#     lPh = sum(Ez1)*(1/2)*log(det(Phi))
#     lPs = sum(Ez0)*(1/2)*log(det(Psi_delta))

#     llik = (lk1-sum(Ez1)*(p/2)*log(2*pi)-lT-((1/2)*(tr1-tr2+tr3))-sum(Ez1)*(q/2)*log(2*pi)-lPh-(1/2)*(tr4-tr5+tr6))+
#     (lk0-sum(Ez0)*(p/2)*log(2*pi)-lPs-((1/2)*(tr7-tr8+tr9))-sum(Ez0)*(K/2)*log(2*pi)-(1/2)*(tr10-tr11+tr12))
#     return llik
# end

function llik_compute(Y::Matrix{Float64},n::Int64,p::Int64,q::Int64,K::Int64,Lambda_1::Matrix{Float64}, Lambda_2::Matrix{Float64}, kappa::Vector{Float64}, mu::Vector{Float64},nu::Vector{Float64},inv_lemma_1,inv_lemma_2,det_lemma_1,det_lemma_2)
    local llik::Float64
    
    llik = sum(log.((kappa .* (exp.((.-1/2) .* diag((Y.-Lambda_1*mu)'*inv_lemma_1*(Y.-Lambda_1*mu)))./(sqrt((2 * pi)^(p) *det_lemma_1)))).+((1 .-kappa) .*((exp.((.-1/2) .* diag((Y.-Lambda_2*nu)'*inv_lemma_2*(Y.-Lambda_2*nu))))./(sqrt((2 * pi)^(p) *det_lemma_2))))))
    return llik
end


function classperf(z,Ez)
    z = Int.(z); Ez = Int.(Ez);
    TrueNegatives = sum(Ez.== 1 .&& z.==  1)
    FalseNegatives = sum(Ez.!= 1 .&& z .== 1)
    FalsePositives = sum(Ez.== 1 .&& z .!= 1)
    TruePositives = sum(Ez.!= 1 .&& z .!= 1)


    SP = TrueNegatives/(TrueNegatives+FalsePositives)
    SE = TruePositives/(TruePositives+FalseNegatives)
    MCC = (TruePositives*TrueNegatives-FalsePositives*FalseNegatives)/sqrt((TruePositives+FalsePositives)*(TruePositives+FalseNegatives)*(TrueNegatives+FalsePositives)*(TrueNegatives+FalseNegatives))
    BACC = (TruePositives/(TruePositives+FalseNegatives)+(TrueNegatives/(TrueNegatives+FalsePositives))) / 2
    ARI = Clustering.randindex(z,Int.(Ez))[1]

    return vcat(hcat("BACC","MCC","SE","SP","ARI"),hcat(BACC,MCC,SE,SP,ARI))
end


function boots_CFAmEFA_careless(Y,X,n,p,q,K,B,std_predictors,L_str,mu,str=missing)
    std_preds = vcat(1,1,std_predictors);numnan = 0
    ncov = size(X,2); Lambda_1_boots = Array{Float64}(undef,p,B); Theta_delta_boots = Array{Float64}(undef,p,B); Phi_boots = Array{Float64}(undef,q*q,B); Lambda_2_boots = Array{Float64}(undef,p*K,B); Psi_delta_boots = Array{Float64}(undef,p,B); beta_boots = Array{Float64}(undef,ncov,B); mu_boots = Array{Float64}(undef,q,B); nu_boots = Array{Float64}(undef,K,B); Ez1_boots = Array{Float64}(undef,n,B);
    for b in 1:B
        println(b)
        boots = sample(1:n,n,replace=true)
        
        
        Y = Y[:,boots]
        X = X[boots,:]
        
        res =  CFAmixEFA(Y, q, K, 500, X, 2000, mu, L_str)
      

        Lambda_1_boots[:,b] = sum(reshape(res[1:(p*q)],p,q),dims=2)
        Lambda_2_boots[:,b] = res[(p*q+1):((p*q)+(p*K))]
        Theta_delta_boots[:,b] = res[((p*q)+(p*K)+1):((p*q)+(p*K)+p)]
        Psi_delta_boots[:,b] = res[((p*q)+(p*K)+p+1):((p*q)+(p*K)+p+p)]
        Phi_boots[:,b] = res[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)]
        Ez1_boots[:,b] = res[((p*q)+(p*K)+p+p+q^2+n+1+1):((p*q)+(p*K)+p+p+q^2+n+1+n)]
        beta_boots[:,b] = res[((p*q)+(p*K)+p+p+q^2+n+1+n+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov)]./std_preds
        mu_boots[:,b] = res[((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q)]
        nu_boots[:,b] = res[((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+K)]

        boo = vcat(Lambda_1_boots[:,b], Theta_delta_boots[:,b], Phi_boots[:,b], Lambda_2_boots[:,b], Psi_delta_boots[:,b], beta_boots[:,b], mu_boots[:,b], nu_boots[:,b])
        if str !== missing
            jldsave(string(str,"/boots_",b,".jld2"),res=boo)
        end
        numnan = numnan + Int(isnan(boo[1,1]))
    end

    
    
    
    Theta_boots = vcat(Lambda_1_boots, Theta_delta_boots, Phi_boots, Lambda_2_boots, Psi_delta_boots, beta_boots, mu_boots, nu_boots)
    Theta_boots = Theta_boots[:,isnan.(Theta_boots[1,:])]
    if size(Theta_boots,2) == 0
        println("All NaN, possibly increase the number of restarts")
        return Theta_boots
    end
    n_theta = size(Theta_boots,1)
    Theta_boots_mean = sum(Theta_boots,dims=2)/B
    cov = zeros(n_theta,n_theta)
    for b in 1:B
        cov = cov + (Theta_boots[:,b] - Theta_boots_mean)*(Theta_boots[:,b] - Theta_boots_mean)'/(B-1)
    end

    se = sqrt.(diag(cov))

    return vcat(se,numnan)
end



function boots_CFAmEFA_faking(Y,X,n,p,q,restarts,K,B,std_predictors,L_str,mu,str=missing)
    std_preds = vcat(1,std_predictors); numnan = 0
    ncov = size(X,2); Lambda_1_boots = Array{Float64}(undef,p,B); Theta_delta_boots = Array{Float64}(undef,p,B); Phi_boots = Array{Float64}(undef,q*q,B); Lambda_2_boots = Array{Float64}(undef,p*K,B); Psi_delta_boots = Array{Float64}(undef,p,B); beta_boots = Array{Float64}(undef,ncov,B); mu_boots = Array{Float64}(undef,q,B); nu_boots = Array{Float64}(undef,K,B); Ez1_boots = Array{Float64}(undef,n,B);
    for b in 1:B
        println(string("boot: ",b))
        boots0 = sample(1:220,220,replace=true)
        boots1 = sample(221:735,515,replace=true)

        Y = Y[:,vcat(boots0,boots1)]
        X = X[vcat(boots0,boots1),:]
        
        
        res =  CFAmixEFA(Y, q, K, restarts, X, 2000, mu, L_str)
      

        Lambda_1_boots[:,b] = sum(reshape(res[1:(p*q)],p,q),dims=2)
        Lambda_2_boots[:,b] = res[(p*q+1):((p*q)+(p*K))]
        Theta_delta_boots[:,b] = res[((p*q)+(p*K)+1):((p*q)+(p*K)+p)]
        Psi_delta_boots[:,b] = res[((p*q)+(p*K)+p+1):((p*q)+(p*K)+p+p)]
        Phi_boots[:,b] = res[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)]
        Ez1_boots[:,b] = res[((p*q)+(p*K)+p+p+q^2+n+1+1):((p*q)+(p*K)+p+p+q^2+n+1+n)]
        beta_boots[:,b] = res[((p*q)+(p*K)+p+p+q^2+n+1+n+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov)]./std_preds
        mu_boots[:,b] = res[((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q)]
        nu_boots[:,b] = res[((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+K)]

        boo = vcat(Lambda_1_boots[:,b], Theta_delta_boots[:,b], Phi_boots[:,b], Lambda_2_boots[:,b], Psi_delta_boots[:,b], beta_boots[:,b], mu_boots[:,b], nu_boots[:,b])
        if str !== missing
            jldsave(string(str,"/boots_",b,".jld2"),res=boo)
        end
        numnan = numnan + Int(isnan(boo[1,1]))
    end

    
    
    
    Theta_boots = vcat(Lambda_1_boots, Theta_delta_boots, Phi_boots, Lambda_2_boots, Psi_delta_boots, beta_boots, mu_boots, nu_boots)
    Theta_boots = Theta_boots[:,isnan.(Theta_boots[1,:])]
    if size(Theta_boots,2) == 0
        println("All NaN, possibly increase the number of restarts")
        return Theta_boots
    end
    n_theta = size(Theta_boots,1)
    Theta_boots_mean = sum(Theta_boots,dims=2)/B
    cov = zeros(n_theta,n_theta)
    for b in 1:B
        cov = cov + (Theta_boots[:,b] - Theta_boots_mean)*(Theta_boots[:,b] - Theta_boots_mean)'/(B-1)
    end

    se = sqrt.(diag(cov))

    return vcat(se,numnan)
end






function CFA(Y::Matrix{Float64}, q::Int64, maxIter::Int64=4000, L_str=nothing)
    local re::Int64, non_conv::Bool, err2::Bool, sum_diff_tf::Bool, retr::Int64, res::Vector{Float64}, start_vals::Vector{Float64}, DomainError_warn::Int64, p::Int64, n::Int64
    n = size(Y,2); p = size(Y,1)
    iter_err = 1; code_error = 0.0;
    while iter_err <= 800 && code_error != 8.0
        #Starting values
        start_vals = starting_val_cfa(n,p,q,L_str)
        
        # EM
        res = EM_cfa(Y,q,p,n,maxIter,reshape(start_vals[1:(p*q)],p,q), Diagonal(start_vals[((p*q)+1):((p*q)+p)]), reshape(start_vals[((p*q)+p+1):((p*q)+p+q^2)],q,q), reshape(start_vals[((p*q)+p+q^2+1):((p*q)+p+q^2+p*q)],p,q))
        
        

        # Controls
        sum_diff_tf = !(res[((p*q)+p+q^2+1+1+1):((p*q)+p+q^2+1+1+1)][1]==1.0)
        non_conv = res[((p*q)+p+q^2+1+1+1+1):((p*q)+p+q^2+1+1+1+1)][1]>=maxIter
        code_error = res[((p*q)+p+q^2+1+1):((p*q)+p+q^2+1+1)][1]
        
        if (code_error == 3.0 || code_error == 0.0) && iter_err < 800
            iter_err += 1
            println(iter_err)
        elseif code_error == 3.0 && iter_err > 800
            return res[1:(p*q)+p+q^2+1+1]
        elseif code_error >= 0.0 && code_error <= 2.0 
            return res[1:(p*q)+p+q^2+1+1]
        elseif non_conv 
            return vcat(res[1:(p*q)+p+q^2+1],4.0)
        elseif sum_diff_tf
            return vcat(res[1:(p*q)+p+q^2+1],5.0)
        end
    end 
    return res[1:(p*q)+p+q^2+1+1]
end





function starting_val_cfa(n::Int64, p::Int64, q::Int64,L_str=nothing)
    local Lambda_1::Matrix{Float64}, Theta_delta::Diagonal{Float64, Vector{Float64}}, Phi::Matrix{Float64}
    if L_str === nothing
        L_str = L_str_fun(p,q)
    end
    Lambda_1 = rand(Uniform(0,1),p,q) .* L_str
    Phi = rand(LKJ(q,1.0))
    
    Theta_delta = Diagonal(vec(1 .- sum(Lambda_1,dims=2).^2));

    return vcat(vec(Lambda_1),vec(diag(Theta_delta)), vec(Phi), vec(L_str))
end

function EM_cfa(Y::Matrix{Float64},q::Int64, p::Int64, n::Int64, maxIter::Int64, Lambda_1::Matrix{Float64}, Theta_delta::Diagonal{Float64, Vector{Float64}}, Phi::Matrix{Float64}, L_str::Matrix)
    local llik::Vector{Float64}, C_yy::Matrix{Float64}, diff_llik::Float64, oldllik::Float64, dtol::Float64, t::Int64, err::Vector{Bool}, err2::Bool, sum_diff_tf::Bool, sum_diff::Vector{Bool}, old_diff::Float64, ESetaeta::Matrix{Float64}, ESyeta::Matrix{Float64}, B1::Matrix{Float64}, b1::Matrix{Float64}  
    t = 1; err2 = false; llik = zeros(Float64, maxIter); oldllik = -Inf; err=zeros(Bool,maxIter); sum_diff = zeros(Bool,maxIter); sum_diff_tf = true;  diff_llik = 1.0; dtol = 0.05; code_error = 8.0;
    while diff_llik>dtol && t<=maxIter && sum_diff_tf
            b1 = Phi*Lambda_1'*inv(Theta_delta)-Phi*Lambda_1'*(Theta_delta\Lambda_1)*((Matrix{Float64}(I,q,q)+Phi*Lambda_1'*(Theta_delta\Lambda_1))\Phi)*Lambda_1'*inv(Theta_delta)
            B1 = Phi-b1*Lambda_1*Phi
            
                
            C_yy = Y*Y'/n 
            ESyeta =   C_yy*b1'
            ESetaeta = B1+b1*ESyeta
        
            # m step
            Lambda_1 = ESyeta*inv(ESetaeta) .* L_str
            
            Theta_delta = Diagonal(C_yy-ESyeta*inv(ESetaeta)*(ESyeta)')
            Phi = sqrt(inv(Diagonal(ESetaeta)))*ESetaeta*sqrt(inv(Diagonal(ESetaeta)))


            
            if sum(Theta_delta.<0) > 0 
                #println("Error Theta_delta")
                code_error = 0.0
                return vcat(vec(fill(NaN,p,q)), vec(diag(fill(NaN,p,p))), vec(fill(NaN,q,q)), NaN, code_error, sum_diff_tf, t-1) 
            elseif det(Phi)<0
                #println("Error Phi")
                code_error = 2.0
                return vcat(vec(fill(NaN,p,q)), vec(diag(fill(NaN,p,p))), vec(fill(NaN,q,q)), NaN, code_error, sum_diff_tf, t-1) 
            else
                old_diff = diff_llik
                llik[t] = llik_compute_cfa(Y,n,p,q,Lambda_1,Theta_delta,Phi)
                    
                diff_llik = abs(oldllik - llik[t])
                if !(oldllik<llik[t])
                    code_error = 3.0
                    return vcat(vec(fill(NaN,p,q)), vec(diag(fill(NaN,p,p))), vec(fill(NaN,q,q)), NaN, code_error, sum_diff_tf, t-1) 
                end
                # println([t llik[t] diff_llik])  

                sum_diff[t] = diff_llik>old_diff               
                oldllik = llik[t]
                t=t+1;
                sum_diff_tf = sum(sum_diff) < 150            
            end
        end
    return vcat(vec(Lambda_1), vec(diag(Theta_delta)), vec(Phi), llik[t-1], code_error, sum_diff_tf, t-1);
end        



function llik_compute_cfa(Y,n,p,q,Lambda_1, Theta_delta,  Phi)
    local llik::Float64, transpY::Adjoint{Float64, Matrix{Float64}}, tr1::Float64, tr2::Float64, sum1::Float64, tr3::Float64, sum2::Float64, tr4::Float64, tr5::Float64, tr6::Float64, sum3::Float64, tr7::Float64, sum4::Float64, tr8::Float64, B1_intra::Matrix{Float64}, B2_intra::Matrix{Float64}, b1_intra::Matrix{Float64}, b2_intra::Matrix{Float64}
    

    b1= Phi*Lambda_1'*inv(Theta_delta)-Phi*Lambda_1'*(Theta_delta\Lambda_1)*((Matrix{Float64}(I,q,q)+Phi*Lambda_1'*(Theta_delta\Lambda_1))\Phi)*Lambda_1'*inv(Theta_delta)
   
    B1 = Phi-b1*Lambda_1*Phi
    llik = 0
    for i=1:n 
        llik = llik + (-(p/2)*log(2*pi)-(1/2)*log(det(Theta_delta))-(1/2)*(Y[:,i]'*inv(Theta_delta)*Y[:,i]-2*Y[:,i]'*inv(Theta_delta)*Lambda_1*(b1*Y[:,i])+tr(Lambda_1'*inv(Theta_delta)*Lambda_1*((Phi - b1*Lambda_1*Phi + b1*Y[:,i]*Y[:,i]'*b1'))))-(q/2)*log(2*pi)-(1/2)*log(det(Phi))-(1/2)*tr(inv(Phi)*((Phi - b1*Lambda_1*Phi + b1*Y[:,i]*Y[:,i]'*b1'))))
    end
    return llik
end




end