module PMLCFAEFA

const PML_VERBOSE = Ref(false)
const PML_FAST = Ref(false)

function pml_em_log(branch, iteration, logcl, delta_stop, delta_aitken, delta_abs, tol, increasing)
    if PML_VERBOSE[]
        println("EM branch=", branch, " iter=", iteration,
                " logCL=", logcl,
                " delta_stop=", delta_stop,
                " delta_abs=", delta_abs,
                " delta_aitken=", delta_aitken,
                " tol=", tol, " increasing=", increasing)
        flush(stdout)
    end
    return nothing
end

#slurm98.


using MixCFAEFA, Sgr, LinearAlgebra, JLD2, Distributions, CategoricalArrays, SpecialFunctions, Optim, Combinatorics, StatsBase, Random, Base.Threads, Clustering, ForwardDiff


export bvnP, aberr_sim, AberrancyGen, CarelessGen, PMLcfaefa, care_sim, classifications, reconstruct_correlation_matrix, reconstruct_L2_matrix, count_and_reduce_rows, pmlcfa, pmlefa, EMparsMix, pmllast, contingency_table, pmlefalm, find_ones_columns, obslik
 

function aberr_sim(design,B,maxit, maxNumRestart)
    local res_B
    tutte=[string("/home/niccolo.cao/pmlmix/res/res_",d,"_",b,".jld2") for d in 1:size(design,1) for b in 1:B]

    strs = readdir("/home/niccolo.cao/pmlmix/res"; join=true)

    results = map(s -> match(r"res_.*", s), strs)

    matches = filter(s -> occursin(r"res_.*", s), strs)
    extracted = map(s -> match(r"res_(.*)", s).captures[1], matches)

    giafatte = string.("/home/niccolo.cao/pmlmix/data/data_",extracted)

    left = setdiff(tutte,strs)


    pattern = r"res_(\d+)_(\d+)"
    numeri = [match(pattern, p).captures for p in left]

    designsn  = [parse(Int, x[1]) for x in numeri]
    replicationsn = [parse(Int, x[2]) for x in numeri]

    desxrep = hcat(designsn,replicationsn)

    designs = unique(designsn)

    ics = JLD2.load("/home/niccolo.cao/pmlmix/data/ics.jld2")["ics"]
    for d in 10:11#d in designs#1:size(design,1)
        repss = desxrep[desxrep[:,1].==d,2]
        @threads for b in repss #1:B
            println("Replication: n. ",b)
            n = 400
            p = 6
            q = 2
            c = 4
            K = Int64(design[d,1])
            propo = design[d,2]
            aberrancy = design[d,3]
            ncov = 2

            # if in(string("/home/niccolo.cao/pmlmix/data/data_",d,"_",b,".jld2"),giafatte)
            #     continue
            # end

            X = JLD2.load(string("/home/niccolo.cao/pmlmix/data/data_",d,"_",b,".jld2"))["data"]#X = JLD2.load(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"))["data"]
            # println("Not done yet: Design n. ",d,", replication n. ",b)


            Y = Matrix{Float64}(X);

            jldsave(string("/home/niccolo.cao/pmlmix/pr/print_res_",d,"_",b,".jld2"),res=repss)


            res_B = PMLcfaefa(X,p,q,K,c,ics,1e-2,1000, 200,1,0.,false,1.)

            jldsave(string("/home/niccolo.cao/pmlmix/res/res_",d,"_",b,".jld2"),res=res_B)#jldsave(string("/home/nc/MEGA/LatentVariables/res/res_",d,"_",b,".jld2"),res=res_B)
        end
    end
end


# function aberr_sim(design,B,maxit, maxNumRestart)
#     local res_B
#     tutte=[string("/home/niccolo.cao/pmlmix/res/res_",d,"_",b,".jld2") for d in 1:size(design,1) for b in 1:B]
    
#     strs = readdir("/home/niccolo.cao/pmlmix/res"; join=true)

#     results = map(s -> match(r"res_.*", s), strs)

#     matches = filter(s -> occursin(r"res_.*", s), strs)
#     extracted = map(s -> match(r"res_(.*)", s).captures[1], matches)

#     giafatte = string.("/home/niccolo.cao/pmlmix/data/data_",extracted) 

#     left = setdiff(tutte,strs)

    
#     pattern = r"res_(\d+)_(\d+)"
#     numeri = [match(pattern, p).captures for p in left]

#     designsn  = [parse(Int, x[1]) for x in numeri]
#     replicationsn = [parse(Int, x[2]) for x in numeri]

#     desxrep = hcat(designsn,replicationsn)

#     designs = unique(designsn)

#     ics = JLD2.load("/home/niccolo.cao/pmlmix/data/ics.jld2")["ics"] 
#     for d in designs#1:size(design,1)
#         repss = desxrep[desxrep[:,1].==d,2]
#         @threads for b in repss #1:B
#             println("Replication: n. ",b)
#             n = 400
#             p = 6
#             q = 2
#             c = 4
#             K = Int64(design[d,1])
#             propo = design[d,2]
#             aberrancy = design[d,3]
#             ncov = 2
          
#             # if in(string("/home/niccolo.cao/pmlmix/data/data_",d,"_",b,".jld2"),giafatte)
#             #     continue
#             # end
          
#             X = JLD2.load(string("/home/niccolo.cao/pmlmix/data/data_",d,"_",b,".jld2"))["data"]#X = JLD2.load(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"))["data"]
#             # println("Not done yet: Design n. ",d,", replication n. ",b)

            
#             Y = Matrix{Float64}(X);

#             jldsave(string("/home/niccolo.cao/pmlmix/pr/print_res_",d,"_",b,".jld2"),res=repss)


#             res_B = PMLcfaefa(X,p,q,K,c,ics,1e-2,1000, 200,1,0.,false,1.)
            
#             jldsave(string("/home/niccolo.cao/pmlmix/res/res_",d,"_",b,".jld2"),res=res_B)#jldsave(string("/home/nc/MEGA/LatentVariables/res/res_",d,"_",b,".jld2"),res=res_B)
#         end
#     end
# end




function AberrancyGen(design,B)
    local kappa0, t10, TruePar::Vector{Float64}, Lambda0::Matrix{Float64}, Theta_delta0::Diagonal{Float64, Vector{Float64}}, Phi0::Matrix{Float64}
    
    seed = vcat(322011611885,11955714332337,14889268790486,19321661366842, 4539386885746,17304354949517, 1981879088971,14248930143082,16766904013492, 1434926186704,22695264190645,17517986428637,17678641845771, 1581593402381, 3465694164649, 9106731443143)
    n = 400
    p = 6
    q = 2
    c = 4 

    ## True CFA
    leftthresh4 =  vcat(-1.9,-1.9,-1.8,-1.75,-1.7,-2.05,0.05,0.05,-0.1,-0.1,0.0,0.1,2.0,1.8,1.85,1.85,1.8,1.95);
    
    A = reshape(Diagonal(ones(q)),q,q);
    L_str = kron(A,ones(Int64(p/q)));
    Lambda0 = [1.1,1.1 ,1.2, 1.05,  1.1, 1.2] .* L_str;
    Phi0 = [1 0.15; 0.15 1]; #Diagonal(ones(q,q))
    Deltas = Diagonal(Diagonal((Lambda0*Phi0*Lambda0')) .+ 1)

    betas08 = vcat(log(0.8/(1 - .8)),0.8)#,rand(Uniform(-1.5,1.5),ncov-1))
    betas09 = vcat(log(0.9/(1 - .9)),0.8)#,rand(Uniform(-1.5,1.5),ncov-1))
    
    # npartot = Int(p*(c-1-2) + p + p + q*(q-1)/2 + p) #per avere vettori output con dimensioni uniformi nonostante la differenza tra disegni sperimentali
    TruePar = vcat(leftthresh4,vec(sum(Lambda0,dims=2)),Phi0[tril!(trues(size(Phi0)), -1)],betas08,betas09)#, vec(zeros(numpar-npartot))) 
    
    jldsave("/home/niccolo.cao/pmlmix/data/TruePar.jld2",TruePars=TruePar) #jldsave("/home/nc/MEGA/LatentVariables/data/TruePar.jld2",TruePars=TruePar) 
   
    ics = hcat(ones(Float64,n),rand(Normal(0,1),n))#vcat(ones(Float64,round(Int,propo*n)),zeros(Float64,round(Int,(1-propo)*n))); 
    jldsave(string("/home/niccolo.cao/pmlmix/data/ics.jld2"),ics=ics) 
     
    # jldsave(string("/home/niccolo.cao/pmlmix/data/ics.jld2"),ics=ics) 

    for d in 1:size(design,1)
        K = Int64(design[d,1])
        propo = design[d,2]
        aberrancy = design[d,3]
            
      

        # Generation
        Random.seed!(seed[d])
        for b in 1:B
            
            Ystarcfa = rand(MvNormal(zeros(p),Hermitian((Deltas^(-1/2))*(Lambda0*Phi0*Lambda0' + diagm(ones(p)))*(Deltas^(-1/2)))),n)'

            

            t10 = hcat(fill(-Inf,p),reshape(leftthresh4,p,c-1),fill(Inf,p))
          
            Ycfa = Matrix{Int64}(undef,n,p);
            for j in 1:p
                Ycfa[:,j] = cut(Ystarcfa[:,j],t10[j,:],labels=1:c) 
            end
            


            ncov = 2
            kappa0 = Vector{Float64}(undef,n);
            if propo == 0.8
                kappa0 = exp.(ics*betas08)./(exp.(ics*betas08).+1)
            elseif propo == 0.9
                kappa0 = exp.(ics*betas09)./(exp.(ics*betas09).+1)
            end
            z = rand.(Binomial.(1, kappa0))

            aberr = Matrix{Int64}(undef,n,p);
            X=zeros(Int64,n,p);
            if aberrancy == 0 #faking slight
                fakparams = Dict("p"=>[1 0],"gam"=>[1.5 0],"del"=>[4 0])  
                PM = repeat(reshape(1 .-z,n,1),outer=(1,p)) # matrix in which 0 means "do not replace" and 1 "replace", to replace a subject's ratings a vector of ones is needed           
                fak = partition_replacement(Ycfa,PM,fakparams);
                X = fak
            elseif aberrancy == 1 #faking extreme
                fakparams = Dict("p"=>[1 0],"gam"=>[4 0],"del"=>[1.5 0]) 
                PM = repeat(reshape(1 .-z,n,1),outer=(1,p)) # matrix in which 0 means "do not replace" and 1 "replace", to replace a subject's ratings a vector of ones is needed           
                fak = partition_replacement(Ycfa,PM,fakparams);
                X = fak
            elseif aberrancy == 2 #careless random
                aberr = reshape(wsample(1:c,fill(0.25,c),n*p,replace=true),n,p) # reshape(repeat(wsample([1,2,3,4],[0.25,0.25,0.25,0.25],n,replace=true),outer=p),n,p)
                for i=1:n
                    if z[i]==1
                        X[i,:] = vec(Ycfa[i,:]);
                    elseif z[i]==0
                        X[i,:] = vec(aberr[i,:]);
                    end
                end
            elseif aberrancy == 3 #careless straightlining
                aberr = reshape(repeat(wsample(1:c,fill(0.25,c),n,replace=true),outer=p),n,p) 
                for i=1:n
                    if z[i]==1
                        X[i,:] = vec(Ycfa[i,:]);
                    elseif z[i]==0
                        X[i,:] = vec(aberr[i,:]);
                    end
                end
            end 
            jldsave(string("/home/niccolo.cao/pmlmix/data/data_",d,"_",b,".jld2"),data=X)
            jldsave(string("/home/niccolo.cao/pmlmix/data/z_",d,"_",b,".jld2"),z=z)
            #jldsave(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"),data=X)  
        end
    end
end





function aberr_simloc(design,B,maxit, maxNumRestart)
    local res_B
    ics = JLD2.load("/home/niccolo.cao/pmlmix/data/ics.jld2")["ics"] 
    for d in 1:size(design,1)
        @threads for b in 1:B
            println("Replication: n. ",b)
            n = design[d,1]
            p = 9
            q = 3
            c = 4
            K = Int64(design[d,2])
            aberrancy = design[d,4]
            ncov = 1
          
            X = JLD2.load(string("/home/niccolo.cao/pmlmix/data/data_",d,"_",b,".jld2"))["data"]#X = JLD2.load(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"))["data"]
            
            Y = Matrix{Float64}(X);




            # ET = reshape(res[(Int(((c-1)*p+p+(q*(q-1))/2+p*K-K*(K-1)/2)+p)+K+ncov+1):(Int(((c-1)*p+p+(q*(q-1))/2+p*K-K*(K-1)/2)+p)+K+ncov+n*size(pairsindex,2)*2)],n,2,size(pairsindex,2))
            # post=sum(ET,dims=3)[:,:,1]
            # post=repeat(1 ./sum(post,dims=2),1,2).*post
            # Ez1 = round.(Int64,post[:,1])

            # if sum(Ez1)/n<0.5
            #     res_B[((p*q)+(p*K)+p+p+q^2+n+1+1):((p*q)+(p*K)+p+p+q^2+n+1+n)] .= 1 .- Ez1
            #     res_B[((p*q)+(p*K)+p+p+q^2+n+1+n+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov)] .= .- bets
            #     println("label unswitched: ",bets)
            # end

            res_B = PMLcfaefa(X,p,q,K,c,ics,1e-6,1000, 20,20,3.,false,1.)
            
            jldsave(string("/home/niccolo.cao/pmlmix/res/res_",d,"_",b,".jld2"),res=res_B)#jldsave(string("/home/nc/MEGA/LatentVariables/res/res_",d,"_",b,".jld2"),res=res_B)
        end
    end
end




function measures(design,B,model)
    local betas0
    D = size(design,1)
    Thresh_rmse = Array{Float64}(undef,D,B); Lambda_1_rmse = Array{Float64}(undef,D,B); Phi_rmse = Array{Float64}(undef,D,B); betas_rmse = Array{Float64}(undef,D,B); 
    Thresh_bias = Array{Float64}(undef,D,B); Lambda_1_bias = Array{Float64}(undef,D,B);  Phi_bias = Array{Float64}(undef,D,B); betas_bias = Array{Float64}(undef,D,B); 
    BACC = Array{Float64}(undef,D,B);  SE = Array{Float64}(undef,D,B);  SP = Array{Float64}(undef,D,B);  MCC = Array{Float64}(undef,D,B); AUC = Array{Float64}(undef,D,B); ARI = Array{Float64}(undef,D,B);
    bTh = Array{Float64}(undef,D,B); bL1 = Array{Float64}(undef,D,B); bPh = Array{Float64}(undef,D,B); bbe = Array{Float64}(undef,D,B); 
    


    ds = 0
    for d in 1:size(design,1)
        ds += 1
        n = design[d,1]
        p = 9
        q = 3
        c = 4
        K = Int64(design[d,2])
        propo = design[d,3]
        aberrancy = design[d,4]
        ncov = 1            
        pairsindex=mapreduce(permutedims,hcat,collect(combinations(1:p,2))')

            
        ics = hcat(ones(Float64,n))

        A = reshape(Diagonal(ones(q)),q,q);
        L_str = kron(A,ones(Int64(p/q)));

        Thresh0 =  vcat(-1.9,-1.9,-1.8,-1.75,-1.7,-2.05,-1.75,-1.75,-2.0,0.05,0.05,-0.1,-0.1,0.0,0.1,-0.1,0.05,-0.05,2.0,1.8,1.85,1.85,1.8,1.95,1.75,1.8,2.1);
        Lambda0 = [1.1,1.1 ,1.2, 1.05,  1.1, 1.2, 1.15,  1.05, 1.3] .* L_str;
        Phi0 = [1 .3 0.35;.3 1 .3; .35 0.3 1]; 
        if propo == 0.8
            betas0 = vcat(log(0.8/(1 - .8)))
        elseif propo == 0.9
            betas0 = vcat(log(0.9/(1 - .9)))
        end

        ok = 0
        for b in 1:B
            z = JLD2.load(string("/home/nc/Downloads/data (1)/data/z_",d,"_",b,".jld2"))["z"]
            res = JLD2.load(string("/home/nc/Downloads/res (1)/res/res_",d,"_",b,".jld2"))["res"]#X = JLD2.load(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"))["data"]
            Thresh = res[1:((c-1)*p)]
            Lambda_1 = res[((c-1)*p+1):((c-1)*p+(p))] 
            Phi = res[Int((c-1)*p+p+1):Int((c-1)*p+p+(q*(q-1))/2)]
            Ez1 = round.(Int64,res[((p*q)+(p*K)+p+p+q^2+n+1+1):((p*q)+(p*K)+p+p+q^2+n+1+n)])

            if model == "inv"
                betas = res[(Int(((c-1)*p+p+(q*(q-1))/2+p*K-K*(K-1)/2)+p)+K+1):(Int(((c-1)*p+p+(q*(q-1))/2+p*K-K*(K-1)/2)+p)+K+ncov)]
                ET = reshape(res[(Int(((c-1)*p+p+(q*(q-1))/2+p*K-K*(K-1)/2)+p)+K+ncov+1):(Int(((c-1)*p+p+(q*(q-1))/2+p*K-K*(K-1)/2)+p)+K+ncov+n*size(pairsindex,2)*2)],n,2,size(pairsindex,2))
            elseif model == "Millsap"
                betas = res[(Int(((c-1)*p+p+(q*(q-1))/2+((c-1-2)*p)+p*K-K*(K-1)/2)+p)+K+1):(Int(((c-1)*p+p+(q*(q-1))/2+((c-1-2)*p)+p*K-K*(K-1)/2)+p)+K+ncov)]
                ET = reshape(res[(Int(((c-1)*p+p+(q*(q-1))/2+((c-1-2)*p)+p*K-K*(K-1)/2)+p)+K+ncov+1):(Int(((c-1)*p+p+(q*(q-1))/2+((c-1-2)*p)+p*K-K*(K-1)/2)+p)+K+ncov+n*size(pairsindex,2)*2)],n,2,size(pairsindex,2))
            end

            post=sum(ET,dims=3)[:,:,1]
            post=repeat(1 ./sum(post,dims=2),1,2).*post
            Ez1 = round.(Int64,post[:,1])
            
            if sum(Ez1)/n<0.5
                Ez1 = 1 .- Ez1
                betas = .- betas
                println("label unswitched: ",betas)
            end
            
            if Ez1[1] !== NaN
                ok += 1
                cls = classifications(z,Ez1)
                println(d," ",b," ",cls)
                BACC[ds,b] = Float64(cls[2,1])
                MCC[ds,b] = Float64(cls[2,2])
                SE[ds,b] = Float64(cls[2,3])
                SP[ds,b] = Float64(cls[2,4])
                
                bTh[ds,b] = mean(Thresh.-Thresh0)
                bL1[ds,b] = mean(Lambda_1.-Lambda0)
                bPh[ds,b] = mean(Phi.-Phi0)
                bbe[ds,b] = mean(betas.-betas0)
                Thresh_rmse[ds,b] = sqrt(mean((Thresh .- Thresh0)*(Thresh.-Thresh0)'))
                Lambda_1_rmse[ds,b] = sqrt(mean((Lambda_1 .- Lambda0)*(Lambda_1.-Lambda0)'))
                Phi_rmse[ds,b] = sqrt(mean((Phi .- Phi0)*(Phi.-Phi0)'))
                betas_rmse[ds,b] = sqrt(mean((betas .- betas0)*(betas.-betas0)'))
                # println(Lambda_1,Lambda0," ",Phi,Phi0vec," ",Theta_delta,Theta_delta0," ",betas,betas0)
            else
                BACC[ds,b] = NaN
                MCC[ds,b] = NaN
                SE[ds,b] = NaN
                SP[ds,b] = NaN
                
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
    BACC = replace(BACC, NaN => 0.0); BACC_mean = mean(reshape(collect(skipmissing(BACC)),D,B),dims=2)[:,1]; BACC_sd = std(reshape(collect(skipmissing(BACC)),D,B),corrected=true,dims=2)[:,1]
    MCC = replace(MCC, NaN => 0.0); MCC_mean = mean(reshape(collect(skipmissing(MCC)),D,B),dims=2)[:,1]; MCC_sd = std(reshape(collect(skipmissing(MCC)),D,B),corrected=true,dims=2)[:,1]
    SE = replace(SE, NaN => 0.0); SE_mean = mean(reshape(collect(skipmissing(SE)),D,B),dims=2)[:,1]; SE_sd = std(reshape(collect(skipmissing(SE)),D,B),corrected=true,dims=2)[:,1]
    SP = replace(SP, NaN => 0.0); SP_mean = mean(reshape(collect(skipmissing(SP)),D,B),dims=2)[:,1]; SP_sd = std(reshape(collect(skipmissing(SP)),D,B),corrected=true,dims=2)[:,1]
    bTh = replace(bTh, NaN => 0.0); bTh_mean = mean(reshape(collect(skipmissing(bTh)),D,B),dims=2)[:,1]; bTh_sd = std(reshape(collect(skipmissing(bTh)),D,B),corrected=true,dims=2)[:,1]
    bL1 = replace(bL1, NaN => 0.0); bL1_mean = mean(reshape(collect(skipmissing(bL1)),D,B),dims=2)[:,1]; bL1_sd = std(reshape(collect(skipmissing(bL1)),D,B),corrected=true,dims=2)[:,1]
    bPh = replace(bPh, NaN => 0.0); bPh_mean = mean(reshape(collect(skipmissing(bPh)),D,B),dims=2)[:,1]; bPh_sd = std(reshape(collect(skipmissing(bPh)),D,B),corrected=true,dims=2)[:,1]
    bbe = replace(bbe, NaN => 0.0); bbe_mean = mean(reshape(collect(skipmissing(bbe)),D,B),dims=2)[:,1]; bbe_sd = std(reshape(collect(skipmissing(bbe)),D,B),corrected=true,dims=2)[:,1]
    Thresh_rmse = replace(Thresh_rmse, NaN => 0.0); Thresh_rmse_mean = mean(reshape(collect(skipmissing(Thresh_rmse)),D,B),dims=2)[:,1]; Thresh_rmse_sd = std(reshape(collect(skipmissing(Thresh_rmse)),D,B),corrected=true,dims=2)[:,1]
    bLambda_1_rmse = replace(Lambda_1_rmse, NaN => 0.0); Lambda_1_rmse_mean = mean(reshape(collect(skipmissing(Lambda_1_rmse)),D,B),dims=2)[:,1]; Lambda_1_rmse_sd = std(reshape(collect(skipmissing(Lambda_1_rmse)),D,B),corrected=true,dims=2)[:,1]
    Phi_rmse = replace(Phi_rmse, NaN => 0.0); Phi_rmse_mean = mean(reshape(collect(skipmissing(Phi_rmse)),D,B),dims=2)[:,1]; Phi_rmse_sd = std(reshape(collect(skipmissing(Phi_rmse)),D,B),corrected=true,dims=2)[:,1]
    etas_rmse = replace(betas_rmse, NaN => 0.0); betas_rmse_mean = mean(reshape(collect(skipmissing(betas_rmse)),D,B),dims=2)[:,1]; betas_rmse_sd = std(reshape(collect(skipmissing(betas_rmse)),D,B),corrected=true,dims=2)[:,1]

    jldsave("/home/MEGA/CFA+EFA/PmlmixSim/resclass.jld2",res=hcat(BACC_mean,BACC_sd,MCC_mean,MCC_sd,SE_mean,SE_sd,SP_mean,SP_sd))
    jldsave("/home/MEGA/CFA+EFA/PmlmixSim/resbias.jld2",res=hcat(bTh_mean,bTh_sd,bL1_mean,bL1_sd,bPh_mean,bPh_sd,bbe_mean,bbe_sd))
    jldsave("/home/MEGA/CFA+EFA/PmlmixSim/resrmse.jld2",res=hcat(Thresh_rmse_mean,Thresh_rmse_sd,Lambda_1_rmse_mean,Lambda_1_rmse_sd,Phi_rmse_mean,Phi_rmse_sd,betas_rmse_mean,betas_rmse_sd))
end






function AberrancyGenloc(design,B)
    local kappa0, t10, TruePar::Vector{Float64}, Lambda0::Matrix{Float64}, Theta_delta0::Diagonal{Float64, Vector{Float64}}, Phi0::Matrix{Float64}
    
    seed = vcat(322011611885,11955714332337,14889268790486,19321661366842, 4539386885746,17304354949517, 1981879088971,14248930143082,16766904013492, 1434926186704,22695264190645,17517986428637,17678641845771, 1581593402381, 3465694164649, 9106731443143)
    p = 9
    q = 3
    c = 4 

    ## True CFA
    leftthresh4 =  vcat(-1.9,-1.9,-1.8,-1.75,-1.7,-2.05,-1.75,-1.75,-2.0,0.05,0.05,-0.1,-0.1,0.0,0.1,-0.1,0.05,-0.05,2.0,1.8,1.85,1.85,1.8,1.95,1.75,1.8,2.1);

    # leftthresh5 = vcat(-1.09,-1.13,-1.68,-2.23,-2.21,-1.61,-2.03,-1.55,-1.52,0.65,0.62,0.06,-0.49,-0.47,0.14,-0.28,0.19,0.22,2.40,2.36,1.81,1.26,1.28,1.88,1.46,1.94,1.97)
    # leftthresh6 = vcat(-1.09,-1.13,-1.68,-2.23,-2.21,-1.61,-2.03,-1.55,-1.52,0.65,0.62,0.06,-0.49,-0.47,0.14,-0.28,0.19,0.22,2.40,2.36,1.81,1.26,1.28,1.88,1.46,1.94,1.97)
    
    # leftthresh10 = vcat(-2.17,-1.66,-2.23,-2.26,-2.19,-1.33,-2.09,-1.28,-1.67,-1.73,-1.22,-1.79,-1.82,-1.75,-0.9,-1.65,-0.84,-1.23,-1.29,-0.78,-1.36,-1.38,-1.32,-0.46,-1.22,-0.4,-0.79,-0.86,-0.35,-0.92,-0.94,-0.88,-0.02,-0.78,0.03,-0.35,-0.42,0.09,-0.48,-0.51,-0.44,0.42,-0.34,0.47,0.08,0.02,0.53,-0.04,-0.07,-0.0,0.85,0.1,0.91,0.52,0.46,0.97,0.39,0.37,0.43,1.29,0.53,1.35,0.96,0.89,1.4,0.83,0.81,0.87,1.73,0.97,1.78,1.4,1.33,1.84,1.27,1.24,1.31,2.17,1.41,2.22,1.83)
   
    # t10 = hcat(fill(-Inf,p),leftthresh,leftthresh.+1,leftthresh.+2,leftthresh.+3.5, fill(Inf,p))
    A = reshape(Diagonal(ones(q)),q,q);
    L_str = kron(A,ones(Int64(p/q)));
    Lambda0 = [1.1,1.1 ,1.2, 1.05,  1.1, 1.2, 1.15,  1.05, 1.3] .* L_str;
    Phi0 = [1 .3 0.35;.3 1 .3; .35 0.3 1]; #Diagonal(ones(q,q))
    # Theta_delta0 = Diagonal(1 .-Diagonal((Lambda0*Phi0*Lambda0'))); #Diagonal([0.56,1.1,0.53,1.2,0.46,0.4,0.73,1.15,59])
    Deltas = Diagonal(Diagonal((Lambda0*Phi0*Lambda0')) .+ 1)

    betas08 = vcat(log(0.8/(1 - .8)))#,rand(Uniform(-1.5,1.5),ncov-1))
    betas09 = vcat(log(0.9/(1 - .9)))#,rand(Uniform(-1.5,1.5),ncov-1))
    
    # npartot = Int(p*(c-1-2) + p + p + q*(q-1)/2 + p) #per avere vettori output con dimensioni uniformi nonostante la differenza tra disegni sperimentali
    TruePar = vcat(leftthresh4,vec(sum(Lambda0,dims=2)),Phi0[tril!(trues(size(Phi0)), -1)],betas08,betas09)#, vec(zeros(numpar-npartot))) 
    
    jldsave("/home/niccolo.cao/pmlmix/data/TruePar.jld2",TruePars=TruePar) #jldsave("/home/nc/MEGA/LatentVariables/data/TruePar.jld2",TruePars=TruePar) 

    
     
    # jldsave(string("/home/niccolo.cao/pmlmix/data/ics.jld2"),ics=ics) 

    for d in 1:size(design,1)
        n = design[d,1]
        K = Int64(design[d,2])
        propo = design[d,3]
        aberrancy = design[d,4]
            
        ics = hcat(ones(Float64,n)); 

        # Generation
        Random.seed!(seed[d])
        @threads for b in 1:B
            
            Ystarcfa = rand(MvNormal(zeros(p),Hermitian((Deltas^(-1/2))*(Lambda0*Phi0*Lambda0' + diagm(ones(p)))*(Deltas^(-1/2)))),n)'

            

            t10 = hcat(fill(-Inf,p),reshape(leftthresh4,p,c-1),fill(Inf,p))
          
            Ycfa = Matrix{Int64}(undef,n,p);
            for j in 1:p
                Ycfa[:,j] = cut(Ystarcfa[:,j],t10[j,:],labels=1:c) 
            end
            


            ncov = 1
            if propo == 0.8
                kappa0 = exp.(ics*betas08)./(exp.(ics*betas08).+1)
            elseif propo == 0.9
                kappa0 = exp.(ics*betas09)./(exp.(ics*betas09).+1)
            end
            z = rand.(Binomial.(1, kappa0))

            aberr = Matrix{Int64}(undef,n,p);
            X=zeros(Int64,n,p);
            if aberrancy == 0. #faking
                fakparams = Dict("p"=>[1 0],"gam"=>[4 0],"del"=>[1.5 0]) #extreme
                PM = repeat(reshape(1 .-z,n,1),outer=(1,p)) # matrix in which 0 means "do not replace" and 1 "replace", to replace a subject's ratings a vector of ones is needed           
                fak = partition_replacement(Ycfa,PM,fakparams);
                X = fak
            elseif aberrancy == 1. #random
                aberr = reshape(wsample([1,2,3,4],[0.25,0.25,0.25,0.25],n*p,replace=true),n,p)# reshape(repeat(wsample([1,2,3,4],[0.25,0.25,0.25,0.25],n,replace=true),outer=p),n,p)
                for i=1:n
                    if z[i]==1
                        X[i,:] = vec(Ycfa[i,:]);
                    elseif z[i]==0
                        X[i,:] = vec(aberr[i,:]);
                    end
                end
            end 
            jldsave(string("/home/niccolo.cao/pmlmix/data/data_",d,"_",b,".jld2"),data=X)
            jldsave(string("/home/niccolo.cao/pmlmix/data/z_",d,"_",b,".jld2"),z=z)
            #jldsave(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"),data=X)  
        end
    end
end




function care_sim(design,B,maxit, maxNumRestart)
    local res_B
    ics = JLD2.load("/home/niccolo.cao/careless/data/ics.jld2")["ics"] 
    for d in 1:size(design,1)
        @threads for b in 1:B
            println("Replication: n. ",b)
            n = 1000
            p = 9
            q = 3
            c = design[d,1]
            K = design[d,2]
            ncov = 2
          
            X = JLD2.load(string("/home/niccolo.cao/careless/data/data_",d,"_",b,".jld2"))["data"]#X = JLD2.load(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"))["data"]
            
            Y = Matrix{Float64}(X);


            
            jldsave(string("/home/niccolo.cao/careless/res/res_",d,"_",b,".jld2"),res=res_B)#jldsave(string("/home/nc/MEGA/LatentVariables/res/res_",d,"_",b,".jld2"),res=res_B)
        end
    end
end




function CarelessGen(design,B)
    local kappa0, t10, TruePar::Vector{Float64}, Lambda0::Matrix{Float64}, Theta_delta0::Diagonal{Float64, Vector{Float64}}, Phi0::Matrix{Float64}
    
    seed = vcat(322011611885,11955714332337,14889268790486,19321661366842, 4539386885746,17304354949517, 1981879088971,14248930143082,16766904013492, 1434926186704,22695264190645,17517986428637,17678641845771, 1581593402381, 3465694164649, 9106731443143, 9729393064462,12160146297624, 4420744722171,16184395927257,19766499429427,20201406362254,19233052068296, 4730430435062)
    n = 250
    p = 9
    q = 3

    ## True CFA
    leftthresh4 = 
    # leftthresh5 = vcat(-1.09,-1.13,-1.68,-2.23,-2.21,-1.61,-2.03,-1.55,-1.52,0.65,0.62,0.06,-0.49,-0.47,0.14,-0.28,0.19,0.22,2.40,2.36,1.81,1.26,1.28,1.88,1.46,1.94,1.97)
    # leftthresh6 = vcat(-1.09,-1.13,-1.68,-2.23,-2.21,-1.61,-2.03,-1.55,-1.52,0.65,0.62,0.06,-0.49,-0.47,0.14,-0.28,0.19,0.22,2.40,2.36,1.81,1.26,1.28,1.88,1.46,1.94,1.97)
    leftthresh10 = vcat(-2.17,-1.66,-2.23,-2.26,-2.19,-1.33,-2.09,-1.28,-1.67,-1.73,-1.22,-1.79,-1.82,-1.75,-0.9,-1.65,-0.84,-1.23,-1.29,-0.78,-1.36,-1.38,-1.32,-0.46,-1.22,-0.4,-0.79,-0.86,-0.35,-0.92,-0.94,-0.88,-0.02,-0.78,0.03,-0.35,-0.42,0.09,-0.48,-0.51,-0.44,0.42,-0.34,0.47,0.08,0.02,0.53,-0.04,-0.07,-0.0,0.85,0.1,0.91,0.52,0.46,0.97,0.39,0.37,0.43,1.29,0.53,1.35,0.96,0.89,1.4,0.83,0.81,0.87,1.73,0.97,1.78,1.4,1.33,1.84,1.27,1.24,1.31,2.17,1.41,2.22,1.83)
    # t10 = hcat(fill(-Inf,p),leftthresh,leftthresh.+1,leftthresh.+2,leftthresh.+3.5, fill(Inf,p))
    A = reshape(Diagonal(ones(q)),q,q);
    L_str = kron(A,ones(Int64(p/q)));
    Lambda0 = [0.69, 0.43, 0.7, 0.36, 0.78, 0.80, 0.56, 0.42, 0.67] .* L_str;
    Phi0 = [1 0.12 0.13;0.12 1 -0.21;.13 -0.21 1]; #Diagonal(ones(q,q))
    Theta_delta0 = Diagonal(1 .-Diagonal((Lambda0*Phi0*Lambda0'))); #Diagonal([0.56,1.1,0.53,1.2,0.46,0.4,0.73,1.15,59])
    betas08 = vcat(log(0.85/(1 - .85)),0.8)#,rand(Uniform(-1.5,1.5),ncov-1))
    betas09 = vcat(log(0.95/(1 - .95)),0.8)#,rand(Uniform(-1.5,1.5),ncov-1))
    # npartot = Int(p*(c-1-2) + p + p + q*(q-1)/2 + p) #per avere vettori output con dimensioni uniformi nonostante la differenza tra disegni sperimentali
    TruePar = vcat(leftthresh4,leftthresh10,vec(sum(Lambda0,dims=2)),Phi0[tril!(trues(size(Phi0)), -1)],betas08,betas09)#, vec(zeros(numpar-npartot))) 
    
    jldsave("/home/niccolo.cao/careless/data/TruePar.jld2",TruePars=TruePar) #jldsave("/home/nc/MEGA/LatentVariables/data/TruePar.jld2",TruePars=TruePar) 

    
    ics = hcat(ones(Float64,n),rand(Normal(0,1),n)); 
    jldsave(string("/home/niccolo.cao/careless/data/ics.jld2"),ics=ics) 

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
                kappa0 = exp.(ics*betas09)./(exp.(ics*betas09).+1)
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
            jldsave(string("/home/niccolo.cao/careless/data/data_",d,"_",b,".jld2"),data=X)
            jldsave(string("/home/niccolo.cao/careless/data/z_",d,"_",b,".jld2"),z=z)
            #jldsave(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"),data=X)  
        end
    end
end







function CarelessGen(design,B,CFAparametrization,carenofake)
    local TruePar::Vector{Float64}, Lambda0::Matrix{Float64}, Theta_delta0::Diagonal{Float64, Vector{Float64}}, Phi0::Matrix{Float64}
    
    seed = vcat(4352746652536,8401675438322,4933401147650,3759925580001,4726536670900,5703165301676,7833702708054,7792120499866, 648698602044,6859789824488, 814784558462,1263536034076)
    n = 1000
    p = 9
    q = 3
    c = 5

    ## True CFA 
    # Intercepts yes, thresholds 0,1
    # if CFAparametrization == 0. 
    #     leftthresh = vcat(-1.6,-1.4,-1.65,-2,-2.35,-1.45,-1.9,-2.5,-1.55)
    #     t10 = hcat(fill(-Inf,p),leftthresh,leftthresh.+1,leftthresh.+2,leftthresh.+3.5, fill(Inf,p))
    #     A = reshape(Diagonal(ones(q)),q,q);
    #     L_str = kron(A,ones(Int64(p/q)));
    #     Lambda0 = [0.69, 0.43, 0.7, 0.36, 0.78, 0.80, 0.56, 0.42, 0.67] .* L_str#[0.69, 0.43, 0.87, 0.36, 0.78, 0.89, 0.56, 0.42, 0.67] .* L_str;
    #     Phi0 = [1 0.12 0.13;0.12 1 -0.21;.13 -0.21 1]#[1 0.12 0.3;0.12 1 -0.21;.3 -0.21 1]; #Diagonal(ones(q,q))
    #     Theta_delta0 =  Diagonal([1.2,1.3,0.7,0.8,1.1,0.5,0.8,.9,0.92])#Diagonal(1 .-Diagonal((Lambda0*Phi0*Lambda0')));
    #     interc10 = [1.3,1.9,2.3,2.6,2.9,3.2,3.5,4.3,5]
    #     # npartot = Int(p*(c-1-2) + p + p + q*(q-1)/2 + p) #per avere vettori output con dimensioni uniformi nonostante la differenza tra disegni sperimentali
    #     TruePar = vcat(t10[:,4],vec(sum(Lambda0,dims=2)),diag(Theta_delta0),Phi0[tril!(trues(size(Phi0)), -1)],interc10)#, vec(zeros(numpar-npartot))) 
    
    # # no intercepts, thresholds 0,1
    # elseif CFAparametrization == 1. 
    #     leftthresh = vcat(-1.6,-1.4,-1.65,-2,-2.35,-1.45,-1.9,-2.5,-1.55)
    #     t10 = hcat(fill(-Inf,p),leftthresh,leftthresh.+1,leftthresh.+2,leftthresh.+3.5, fill(Inf,p))
    #     A = reshape(Diagonal(ones(q)),q,q);
    #     L_str = kron(A,ones(Int64(p/q)));
    #     Lambda0 = [0.69, 0.43, 0.7, 0.36, 0.78, 0.80, 0.56, 0.42, 0.67] .* L_str;
    #     Phi0 = [1 0.12 0.13;0.12 1 -0.21;.13 -0.21 1]; #Diagonal(ones(q,q))
    #     Theta_delta0 =  Diagonal([1.2,1.3,0.7,0.8,1.1,0.5,0.8,.9,0.92])#Diagonal(1 .-Diagonal((Lambda0*Phi0*Lambda0')));
    #     # npartot = Int(p*(c-1-2) + p + p + q*(q-1)/2 + p) #per avere vettori output con dimensioni uniformi nonostante la differenza tra disegni sperimentali
    #     TruePar = vcat(t10[:,4],vec(sum(Lambda0,dims=2)),diag(Theta_delta0),Phi0[tril!(trues(size(Phi0)), -1)])#, vec(zeros(numpar-npartot))) 
    
    # # threshold parametrization
    # elseif CFAparametrization == 2.
    #     leftthresh = vcat(-1.6,-1.4,-1.65,-2,-2.35,-1.45,-1.9,-2.5,-1.55)
    #     t10 = hcat(fill(-Inf,p),leftthresh,leftthresh.+1,leftthresh.+2,leftthresh.+3.5, fill(Inf,p))
    #     A = reshape(Diagonal(ones(q)),q,q);
    #     L_str = kron(A,ones(Int64(p/q)));
    #     Lambda0 = [0.69, 0.43, 0.7, 0.36, 0.78, 0.80, 0.56, 0.42, 0.67] .* L_str;
    #     Phi0 = [1 0.12 0.13;0.12 1 -0.21;.13 -0.21 1]; #Diagonal(ones(q,q))
    #     Theta_delta0 = Diagonal(1 .-Diagonal((Lambda0*Phi0*Lambda0')));
    #     # npartot = Int(p*(c-1-2) + p + p + q*(q-1)/2 + p) #per avere vettori output con dimensioni uniformi nonostante la differenza tra disegni sperimentali
    #     TruePar = vcat(vec(t10[:,2:(c)]),vec(sum(Lambda0,dims=2)),Phi0[tril!(trues(size(Phi0)), -1)])#, vec(zeros(numpar-npartot))) 
    # end
    leftthresh5  = vcat(-1.09,-1.13,-1.68,-2.23,-2.21,-1.61,-2.03,-1.55,-1.52,0.65,0.62,0.06,-0.49,-0.47,0.14,-0.28,0.19,0.22,2.40,2.36,1.81,1.26,1.28,1.88,1.46,1.94,1.97)
    leftthresh10 = vcat(-2.17,-1.66,-2.23,-2.26,-2.19,-1.33,-2.09,-1.28,-1.67,-1.73,-1.22,-1.79,-1.82,-1.75,-0.9,-1.65,-0.84,-1.23,-1.29,-0.78,-1.36,-1.38,-1.32,-0.46,-1.22,-0.4,-0.79,-0.86,-0.35,-0.92,-0.94,-0.88,-0.02,-0.78,0.03,-0.35,-0.42,0.09,-0.48,-0.51,-0.44,0.42,-0.34,0.47,0.08,0.02,0.53,-0.04,-0.07,-0.0,0.85,0.1,0.91,0.52,0.46,0.97,0.39,0.37,0.43,1.29,0.53,1.35,0.96,0.89,1.4,0.83,0.81,0.87,1.73,0.97,1.78,1.4,1.33,1.84,1.27,1.24,1.31,2.17,1.41,2.22,1.83)
    # t10 = hcat(fill(-Inf,p),leftthresh,leftthresh.+1,leftthresh.+2,leftthresh.+3.5, fill(Inf,p))
    A = reshape(Diagonal(ones(q)),q,q);
    L_str = kron(A,ones(Int64(p/q)));
    Lambda0 = [0.69, 0.43, 0.7, 0.36, 0.78, 0.80, 0.56, 0.42, 0.67] .* L_str;
    Phi0 = [1 0.12 0.13;0.12 1 -0.21;.13 -0.21 1]; #Diagonal(ones(q,q))
    Theta_delta0 = Diagonal(1 .-Diagonal((Lambda0*Phi0*Lambda0')));
    betas085 = vcat(log(0.85/(1 - .85)),0.8)#,rand(Uniform(-1.5,1.5),ncov-1))
    betas09 = vcat(log(0.95/(1 - .95)),0.8)#,rand(Uniform(-1.5,1.5),ncov-1))
    # npartot = Int(p*(c-1-2) + p + p + q*(q-1)/2 + p) #per avere vettori output con dimensioni uniformi nonostante la differenza tra disegni sperimentali
    TruePar = vcat(leftthresh4,leftthresh10,vec(sum(Lambda0,dims=2)),Phi0[tril!(trues(size(Phi0)), -1)],betas085,betas09)#, vec(zeros(numpar-npartot))) 
    

    jldsave("/home/niccolo.cao/careless/data/TruePar.jld2",TruePars=TruePar) #jldsave("/home/nc/MEGA/LatentVariables/data/TruePar.jld2",TruePars=TruePar) 

    
    if carenofake == true
        for d in 1:size(design,1)
            K = design[d,1]
            propo = design[d,2]
            careless = design[d,3]
            

            # Generation
            Random.seed!(seed[d])
            for b in 1:B
                Eta0 = rand(MvNormal(zeros(q),Phi0),n);
                Delta0 = rand(MvNormal(zeros(p),Theta_delta0),n);

                if CFAparametrization == 0. 
                    Ystarcfa = interc10' .+  (Eta0'*Lambda0' + Delta0');

                elseif CFAparametrization == 1. 
                    Ystarcfa = Eta0'*Lambda0' + Delta0';

                elseif CFAparametrization == 2. 
                    Ystarcfa = Eta0'*Lambda0' + Delta0';
                end

                Ycfa = Matrix{Int64}(undef,n,p);
                for j in 1:p
                    Ycfa[:,j] = cut(Ystarcfa[:,j],t10[j,:],labels=1:c) 
                end
                


                ncov = 2
                ics = hcat(ones(Float64,n),rand(Normal(0,1),n)); 
                betas0 = vcat(log(propo/(1 - propo)),0.8)#,rand(Uniform(-1.5,1.5),ncov-1))
                kappa0 = exp.(ics*betas0)./(exp.(ics*betas0).+1)
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
                jldsave(string("/home/niccolo.cao/careless/data/data_",d,"_",b,".jld2"),data=X)
                jldsave(string("/home/niccolo.cao/careless/data/z_",d,"_",b,".jld2"),z=z)  #jldsave(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"),data=X)  
            end
        end
        
    elseif carenofake == false
        for d in 1:size(design,1)
            K = design[d,1]
            propo = design[d,2]
            faking = design[d,3]
            
    
            # Generation
            Random.seed!(seed[d])
            for b in 1:B
                Eta0 = rand(MvNormal(zeros(q),Phi0),n);
                Delta0 = rand(MvNormal(zeros(p),Theta_delta0),n);
        
                if CFAparametrization == 0. 
                    Ystarcfa = interc10' .+  (Eta0'*Lambda0' + Delta0');
        
                elseif CFAparametrization == 1. 
                    Ystarcfa = Eta0'*Lambda0' + Delta0';
        
                elseif CFAparametrization == 2. 
                    Ystarcfa = Eta0'*Lambda0' + Delta0';
                end
        
                Ycfa = Matrix{Int64}(undef,n,p);
                for j in 1:p
                    Ycfa[:,j] = cut(Ystarcfa[:,j],t10[j,:],labels=1:c) 
                end
                    
                fakparams = Matrix{Int64}(undef,n,p);
                if faking == 0. #slight
                    fakparams =  Dict("p"=>[1 0],"gam"=>[1.5 0],"del"=>[4 0])
                elseif faking == 1. #extreme
                    fakparams = Dict("p"=>[1 0],"gam"=>[4 0],"del"=>[1.5 0])
                end 
            
                PM = vcat(zeros(Int(propo*n),p),ones(Int(round((1-propo)*n)),p)) # matrix in which 0 means "do not replace" and 1 "replace", to replace a subject's ratings a vector of ones is needed           
                fak = partition_replacement(Ycfa,PM,fakparams);
                
                jldsave(string("/home/niccolo.cao/careless/data/data_",d,"_",b,".jld2"),data=fak)#jldsave(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"),data=fak)  
            end
        end
    end
end



###################################################################################### Functions for Thresholds CFA parametrization #######################################################################################


function PMLcfaefa(X,p,q,K,c,ics,tol,maxit, maxNumRestart::Int64,replix=10,fixthrehsolds=0.,L_str=false,latentmeans=1.,parametr=false,verbose=false,fast=false)#param_initialii=false,
    local res::Vector{Float64}, Ez1, bets, trcfainit, trefainit,param_initialii
    PML_VERBOSE[] = Bool(verbose)
    PML_FAST[] = Bool(fast)
    if fast && !(latentmeans == 1.0 && fixthrehsolds == 0.0)
        throw(ArgumentError("fast supports latentmeans=1 and fixthrehsolds=0 only"))
    end
    n = size(X,1); ncov = size(ics,2)
    
    if L_str == false A = reshape(Diagonal(ones(q)),q,q);L_str = kron(A,ones(Int64(p/q))) end
    pairsindex=mapreduce(permutedims,hcat,collect(combinations(1:p,2))');
    iter_err = 1; code_error = 0.0;  err_0_count = 0; A2 = diagm(ones(K));  indices = [i + (j - 1) * p for j in 1:K for i in j:p]   
    while code_error == 0.0
        if PML_VERBOSE[]
            println("PML restart=", iter_err)
            flush(stdout)
        end
        # println("Restart: n. ", iter_err)

        # replix = 10
        # liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+p)+ncov,replix);
        # for i in 1:replix
        #     ## Starting values: 
        #     # 1. CFA+EFA contiunous
        #     Y = Matrix{Float64}(X);

        #     #some minimal variation
        #     for i in 1:size(Y,1)
        #         Y[i,:] = Y[i,:] + rand(Uniform(-0.01,0.01),size(Y,2))
        #     end
        #     #standarization of observed variables
        #     for i in 1:size(Y,2)
        #         Y[:,i]=(Y[:,i].-mean(Y,dims=1)[i])/std(Y,dims=1)[i]
        #     end

        #     Y = Matrix(Y');
        #     # ics = ones(n,1)
        #     # ncov = 1
        #     repli = 20
        #     Ez1s = Matrix{Float64}(undef,n,repli); liks = Vector{Float64}(undef,repli); betss = Matrix{Float64}(undef,ncov,repli); Entropy = Vector{Float64}(undef,repli); 
        #     for i in 1:repli
        #         res = CFAmixEFA(Y, q, K, 500, ics, 20000, 0,L_str)

        #         Lambda_1 = reshape(res[1:(p*q)],p,q)
        #         Lambda_2 = reshape(res[(p*q+1):((p*q)+(p*K))],p,K)
        #         Theta_delta = Diagonal(res[((p*q)+(p*K)+1):((p*q)+(p*K)+p)])
        #         Psi_delta = Diagonal(res[((p*q)+(p*K)+p+1):((p*q)+(p*K)+p+p)])
        #         Phi = reshape(res[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)],q,q)
        #         kappa = res[((p*q)+(p*K)+p+p+q^2+1):((p*q)+(p*K)+p+p+q^2+n)]
        #         nu = res[((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+K)]


        #         inv_lemma_1 = inv_lemma_1_fun(q,Lambda_1,Theta_delta,Phi)
        #         inv_lemma_2 = inv_lemma_2_fun(K,Lambda_2,Psi_delta)      
        #         det_lemma_1 = det_lemma_1_fun(Lambda_1,Theta_delta,Phi)
        #         det_lemma_2 = det_lemma_2_fun(K,Lambda_2,Psi_delta)
        #         b1 = Phi*Lambda_1'*inv(Lambda_1*Phi*Lambda_1'+Theta_delta);
        #         b2 = Lambda_2'*inv(Lambda_2*Lambda_2' + Psi_delta);
        #         P1p = zeros(n)
        #         P1 = zeros(n)
        #         for i in 1:n
        #             P1p[i] = ((kappa[i] .* ((exp((-1/2) * Y[:,i]'*inv_lemma_1*Y[:,i]))/(sqrt((2 .* pi)^(p) *det_lemma_1))))/((kappa[i] .* ((exp((-1/2) * Y[:,i]'*inv_lemma_1*Y[:,i]))/(sqrt((2 .* pi)^(p) *det_lemma_1)))).+((exp((-1/2) * (Y[:,i].-Lambda_2*nu)'*inv_lemma_2*(Y[:,i].-Lambda_2*nu)))/(sqrt((2 .* pi)^(p) .*det_lemma_2))).-(kappa[i] .*((exp((-1/2) *(Y[:,i].-Lambda_2*nu)'*inv_lemma_2*(Y[:,i].-Lambda_2*nu)))/(sqrt((2 .* pi)^(p) .*det_lemma_2))))))
        #         end
        #         P1p = P1p .+ 1e-7
        #         P1p[findall(P1p.>1)].=0.9999999999999
        #         P0 = zeros(n)
        #         P1p[findall(P1p.>1)].=0.9999999999999
        #         for j=1:n
        #             P0[j] = (1-P1p[j])*log(1-P1p[j])
        #             P1[j] = P1p[j]*log(P1p[j])
        #         end
        #         P0[isnan.(P0)].=0.0
        #         Ek = - (sum(P1) + sum(P0))
        #         Entropy[i] = 1-(Ek/(n*log(2)))

        #         Ez1s[:,i] = res[((p*q)+(p*K)+p+p+q^2+n+1+1):((p*q)+(p*K)+p+p+q^2+n+1+n)]
        #         betss[:,i] = res[((p*q)+(p*K)+p+p+q^2+n+1+n+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov)]
        #         liks[i] = res[((p*q)+(p*K)+p+p+q^2+n+1):((p*q)+(p*K)+p+p+q^2+n+1)][1]
        #     end

        #     Ez1 = round.(Int64,Ez1s[:,findall((liks).==maximum(liks))])[:,1]
        #     bets = betss[:,findall((liks).==maximum(liks))]
            
        #     #### Entropy
        #     # if sum((Entropy).==maximum(Entropy))==1
        #     #     Ez1 = round.(Int64,Ez1s[:,findall((Entropy).==maximum(Entropy))])[:,1]
        #     #     bets = betss[:,findall((Entropy).==maximum(Entropy))]
        #     # elseif sum((Entropy).==maximum(Entropy))>1
        #     #     Ezone = Ez1s[:,findall((Entropy).==maximum(Entropy))]
        #     #     betsone = betss[:,findall((Entropy).==maximum(Entropy))]
        #     #     likone = liks[findall((Entropy).==maximum(Entropy))]
                
        #     #     Ez1 = round.(Int64,Ezone[:,findall(likone.==maximum(likone))])[:,1]
        #     #     bets = betss[:,findall(likone.==maximum(likone))]
        #     # end



        #     if sum(Ez1)/n<0.5
        #         Ez1 = 1 .- Ez1
        #         bets = .- bets
        #         println("label unswitched: ",bets)
        #     end

        #     # 2. CFA ordinal PML
        #     if c == 4
        #         trcfainit = rand(Uniform(-1.8,-0.8),1)
        #         trcfainit = hcat(trcfainit,zeros(1),.-trcfainit)
        #     elseif c > 4
        #         trcfainit = rand(Uniform(-1.6,-1.5),1).+reshape(repeat(collect(range(0,3.5,c-1)),inner=1),1,(c-1))
        #     end
            
        #     paramcfa = vcat(vec(trcfainit),rand(Uniform(0.1,1),p),rand(Uniform(1e-3,1),Int(q*(q-1)/2)))
        #     lb=vec(vcat(fill(-5,(c-1)),fill(0.001,p),fill(-1.,Int((q*(q-1))/2))))#,fill(0.0001,p)));fill(-5,(c-1)).+repeat(collect(0:(c-2)),inner=p)
        #     ub=vec(vcat(fill(5,(c-1)),fill(1.,p),fill(1.,Int((q*(q-1))/2))))#,fill(4,p))); vcat(fill(5,(c-1)).-repeat(collect((c-2):-1:0),inner=p)
        #     parcfa = optimize(param ->pmlcfath(param,p,q,c,X[Bool.(round.(Ez1)),:],L_str),lb,ub,paramcfa,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)
    

        #     tau1=parcfa.minimizer[(1):((c-1))]
        #     tau1=vec(sort(reshape(tau1,1,c-1),dims=2))

        #     # 3. EFA ordinal PML
        #     paramefa = vcat(rand(Uniform(0.1,1.),Int(p*K-(K*(K-1)/2))),rand(Uniform(0.1,3),p),rand(Uniform(0.5,2),p))# sort(rand(Uniform(1.2,4),p*(c-1-2)).+repeat(collect(0:(c-2-2)),inner=p))
        #     lb=vec(vcat(fill(1e-3,Int((p*K-(K*(K-1)/2)))),fill(0.001,p),fill(-10,p)))#,fill(0.0001,p))); .+repeat(collect(0:(c-2-2)),inner=p)
        #     ub=vec(vcat(fill(1.,Int((p*K-(K*(K-1)/2)))),fill(3.,p),fill(10.,p)))#,fill(4,p)));  .-repeat(collect((c-2-2):-1:0),inner=p)
        #     parefa = optimize(param ->pmlefathfixed(param,p,K,c,X[Bool.(1 .-Ez1),:],indices,tau1),lb,ub,paramefa,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

        #     # in case of very strong straightlining
        #     # if parefa.minimizer[((c-1-2)*p+1):Int((c+1+2)*p+(p*K-(K*(K-1)/2)))] .== ones(Float64,Int(p*K-(K*(K-1)/2)))
        #     #     parefa.minimizer[((c-1-2)*p+1):Int((c+1+2)*p+(p*K-(K*(K-1)/2)))]  .= rand(Uniform(0.8,0.98),Int(p*K-(K*(K-1)/2)))
        #     # end
        #     # if parefa.minimizer[Int((c+1+2)*p+(p*K-(K*(K-1)/2))+1):Int((c+1+2)*p+(p*K-(K*(K-1)/2))+p)] .== fill(0.01,p)
        #     #     parefa.minimizer[Int((c+1+2)*p+(p*K-(K*(K-1)/2))+1):Int((c+1+2)*p+(p*K-(K*(K-1)/2))+p)]  .= rand(Uniform(0.02,0.2),p)
        #     # end
        
            
           
        #     lambda1=parcfa.minimizer[((c-1)+1):((c-1)+p)]
        #     phi=parcfa.minimizer[((c-1)+p+1):((c-1)+p+Int(q*(q-1)/2))]
        #     lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
        #     lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
        #     phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
        #     lambda2=parefa.minimizer[(1):(Int((p*K-(K*(K-1)/2))))]
        #     psi=parefa.minimizer[(Int((p*K-(K*(K-1)/2)))+1):(Int((p*K-(K*(K-1)/2)))+p)]
        #     interc=parefa.minimizer[(Int((p*K-(K*(K-1)/2)))+p+1):(Int((p*K-(K*(K-1)/2)))+p+p)]
        #     lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
        #     lambda2[lambda2.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.99))
        #     psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
            
        #     param_inits = vcat(tau1,lambda1,phi,lambda2,psi,interc,bets)
        
        #     # # 2. CFA ordinal PML
        #     # paramcfa = vcat(sort(rand(Uniform(-4.9,-2.3),p*(c-1))).+repeat(collect(0:(c-2)),inner=p),rand(Uniform(0.11,0.9),p),rand(Uniform(-0.9,0.9),Int(q*(q-1)/2)))
        #     # lb=vec(vcat(fill(-5,p*(c-1)).+repeat(collect(0:(c-2)),inner=p),fill(0.1001,p),fill(-0.9,Int((q*(q-1))/2))))#,fill(0.0001,p)));
        #     # ub=vec(vcat(fill(3,p*(c-1)).-repeat(collect((c-2):-1:0),inner=p),fill(0.9,p),fill(0.9,Int((q*(q-1))/2))))#,fill(4,p)));
        #     # parcfa = optimize(param ->pmlcfa(param,p,q,c,X[Bool.(Ez1),:],L_str),lb,ub,paramcfa,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

        #     # # 3. EFA ordinal PML
        #     # paramefa = vcat(sort(rand(Uniform(1.2,4),p*(c-1-2)).+repeat(collect(0:(c-2-2)),inner=p)),rand(Uniform(0.11,0.9),Int(p*K-(K*(K-1)/2))),rand(Uniform(0.1,2.9),p),rand(Uniform(-9.90,9.90),p)) 
        #     # lb=vec(vcat(fill(1.1,p*(c-1-2)).+repeat(collect(0:(c-2-2)),inner=p),fill(0.10001,Int((p*K-(K*(K-1)/2)))),fill(0.101,p),fill(-9.90,p)))#,fill(0.0001,p)));
        #     # ub=vec(vcat(fill(6,p*(c-1-2)).-repeat(collect((c-2-2):-1:0),inner=p),fill(0.9,Int((p*K-(K*(K-1)/2)))),fill(2.9,p),fill(9.90,p)))#,fill(4,p)));
        #     # parefa = optimize(param ->pmlefa(param,p,K,c,X[Bool.(1 .-Ez1),:]),lb,ub,paramefa,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)


        #     # param_init = vcat(parcfa.minimizer.+rand(Uniform(-0.1,0.1),size(parcfa.minimizer,1)),parefa.minimizer.+rand(Uniform(-0.1,0.1),size(parefa.minimizer,1)),sum(Ez1)/n)
        #     println(round.(param_inits,digits=3),sum(Ez1s[:,findall((liks).==maximum(liks))])/n," ",sum(Ez1)/n," ", maximum(liks)," ",maximum(liks))
            

        #     # replis = 10
        #     # paramss = zeros(Float64,size(param_inits,1),replis); likss = zeros(Float64,replis) 
        #     # for i in 1:replis
        #     #     println(i)
        #     #     plimix = E_stepMix(param_inits,q,K,p,n,c,X,pairsindex,L_str,true);
        #     #     T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        #     #     t = time()
        #     #     paramss[:,i] .= M_parsMix(T,param_inits,X,n,p,q,K,c,pairsindex,L_str,indices,1,t,plimix,true); #solo primo argomento di output
        #     #     likss[i] = E_stepMix(paramss[:,i],q,K,p,n,c,X,pairsindex,L_str,false,T);
        #     # end 
        #     # println(likss,likss[findall((likss).==maximum(likss))],likss[findall((likss).==minimum(likss))])
        #     # paramsss = paramss[:,findall((likss).==maximum(likss))]
        #     # println("started")

        #     # first iteration
        #     plimix = E_stepMixonlyth(param_inits,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        #     T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        #     t = time()
        #     params[:,i] = M_parsMixonlyth(T,param_inits,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        #     dt = time() - t
        #     liksmix[i] = E_stepMixonlyth(params[:,i],q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
        # end
        # param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
        # lik=liksmix[findall(liksmix.==maximum(liksmix))]
    
        # tau1=param_initials0[(1):((c-1))]
        # tau1=vec(sort(reshape(tau1,1,c-1),dims=2))
        # lambda1=param_initials0[((c-1)+1):((c-1)+p)]
        # phi=param_initials0[((c-1)+p+1):((c-1)+p+Int(q*(q-1)/2))]
        # lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
        # lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
        # phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
        # lambda2=param_initials0[((c-1)+p+Int(q*(q-1)/2)+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2))))]
        # psi=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+p)]
        # interc=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+p+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+p+p)]
        # lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
        # lambda2[lambda2.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.99))
        # psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))

        # param_initialii = vcat(tau1,lambda1,phi,lambda2,psi,interc,bets)

        # res =  EMparsMixonlyth(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
        
        
        # # if param_initialii != false && iter_err == 1
            
        #     if fixthrehsolds == 0.
        #         res =  EMparsMix(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
        #     elseif fixthrehsolds == 1.
        #         res =  EMparsMixthefa(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
        #     elseif fixthrehsolds == 2.
        #         res =  EMparsMixth(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
        #     end
        #     if res[size(res,1)] == 0.0
        #         return "change statring values"
        #     end
        # else

            ### Random start + 1 iteration
            # EM
            if latentmeans == 0.
                if fixthrehsolds == 0.
                    param_initialii = startvals(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMix(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                elseif fixthrehsolds == 1.
                    param_initialii = startvalsthefa(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMixthefa(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                elseif fixthrehsolds == 2.
                    # thresholds unici per items ma non invarianti
                    param_initialii = startvalsth(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMixth(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                elseif fixthrehsolds == 3.
                    ## NON ANCORA FATTO
                    param_initialii = startvalsth(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMixth(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                elseif fixthrehsolds == 4.
                    # thresholds invariati e unici 
                    param_initialii = startvalsonlyth(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMixonlyth(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                elseif fixthrehsolds == 5.
                    # thresholds invariati e unici con 0 e 1... poco senso con una cfa senza intercetta
                    param_initialii = startvalsonlyth01(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMixonlyth01(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                elseif fixthrehsolds == 6.
                    # thresholds invariati e unici con 0
                    param_initialii = startvalsonlyth0intercCFA(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMixonlyth0intercCFA(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                elseif fixthrehsolds == 7.
                    # thresholds invariati e unici + intercetta CFA, problemi in identificabilità della mistura
                    param_initialii = startvalsonlyth01intercCFA(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMixonlyth01intercCFA(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                end
            elseif latentmeans == 1.
                if fixthrehsolds == 0.
                    # println("millsap")
                    # Millsap parametrization: 2 first EFA thresholds fixed to the CFA ones
                    param_initialii = startvalslm(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMixlm(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                elseif fixthrehsolds == -1.
                    PML_VERBOSE[] && println("millsap binary")
                    # Millsap parametrization: invariant EFA thresholds fixed to the CFA ones
                    param_initialii = startvalslmbin(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMixlmbin(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
            elseif fixthrehsolds == 1.
                    ## NON ANCORA FATTO
                    param_initialii = startvalsthefa(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMixthefa(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                elseif fixthrehsolds == 2.
                   
                    param_initialii = startvalsth(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMixth(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                elseif fixthrehsolds == 3.
                    # Thresholds invarianti + latent means EFA + Psi
                    PML_VERBOSE[] && println("t invariant but one set of per item")
                    if parametr != false
                        param_initialii = parametr
                    else
                        param_initialii = startvalsinv(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    end
                    res =  EMparsMixinv(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                elseif fixthrehsolds == 4.
                    # SOLO PER K=1: Threhsolds invarianti + latent means EFA + no Psi
                    param_initialii = startvalslminv(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMixlminv(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                elseif fixthrehsolds == 5.
                    PML_VERBOSE[] && println("t invariant and equal between-items")
                    # Threhsolds invarianti e unici + latent means EFA + Psi
                    param_initialii = startvalsthinv(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                    res =  EMparsMixthinv(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
                end
            elseif latentmeans == 2.
                # SOLO PER K=1: Threhsolds invarianti + no latent means EFA + no Psi  
                param_initialii = startvalsl2inv(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
                res =  EMparsMixl2inv(param_initialii,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
            end
        # end
        # Controls
        code_error = res[size(res,1)]
        
        if code_error == 0.0  && iter_err < maxNumRestart
            err_0_count += 1.
            iter_err += 1.
        end
    end 
    return vcat(res,err_0_count) #numero parametri + codice errore + numero restart
end



## first version of rational starting values

# replix = 10
#         liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)*p+p+(q*(q-1))/2+(c-1-2)*p+(p*K-(K*(K-1)/2))+p+p)+ncov,replix);
#         for i in 1:replix
#             ## Starting values: 
#             # 1. CFA+EFA contiunous
#             Y = Matrix{Float64}(X);

#             #some minimal variation
#             for i in 1:size(Y,1)
#                 Y[i,:] = Y[i,:] + rand(Uniform(-0.35,0.35),size(Y,2))
#             end
#             #standarization of observed variables
#             for i in 1:size(Y,2)
#                 Y[:,i]=(Y[:,i].-mean(Y,dims=1)[i])/std(Y,dims=1)[i]
#             end

#             Y = Matrix(Y');
#             # ics = ones(n,1)
#             # ncov = 1
#             repli = 40
#             Ez1s = Matrix{Float64}(undef,n,repli); liks = Vector{Float64}(undef,repli); betss = Matrix{Float64}(undef,ncov,repli); Entropy = Vector{Float64}(undef,repli); 
#             for i in 1:repli
#                 res = CFAmixEFA(Y, q, K, 500, ics, 20000, 0,L_str)

#                 Lambda_1 = reshape(res[1:(p*q)],p,q)
#                 Lambda_2 = reshape(res[(p*q+1):((p*q)+(p*K))],p,K)
#                 Theta_delta = Diagonal(res[((p*q)+(p*K)+1):((p*q)+(p*K)+p)])
#                 Psi_delta = Diagonal(res[((p*q)+(p*K)+p+1):((p*q)+(p*K)+p+p)])
#                 Phi = reshape(res[((p*q)+(p*K)+p+p+1):((p*q)+(p*K)+p+p+q^2)],q,q)
#                 kappa = res[((p*q)+(p*K)+p+p+q^2+1):((p*q)+(p*K)+p+p+q^2+n)]
#                 nu = res[((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov+q+K)]


#                 inv_lemma_1 = inv_lemma_1_fun(q,Lambda_1,Theta_delta,Phi)
#                 inv_lemma_2 = inv_lemma_2_fun(K,Lambda_2,Psi_delta)      
#                 det_lemma_1 = det_lemma_1_fun(Lambda_1,Theta_delta,Phi)
#                 det_lemma_2 = det_lemma_2_fun(K,Lambda_2,Psi_delta)
#                 b1 = Phi*Lambda_1'*inv(Lambda_1*Phi*Lambda_1'+Theta_delta);
#                 b2 = Lambda_2'*inv(Lambda_2*Lambda_2' + Psi_delta);
#                 P1p = zeros(n)
#                 P1 = zeros(n)
#                 for i in 1:n
#                     P1p[i] = ((kappa[i] .* ((exp((-1/2) * Y[:,i]'*inv_lemma_1*Y[:,i]))/(sqrt((2 .* pi)^(p) *det_lemma_1))))/((kappa[i] .* ((exp((-1/2) * Y[:,i]'*inv_lemma_1*Y[:,i]))/(sqrt((2 .* pi)^(p) *det_lemma_1)))).+((exp((-1/2) * (Y[:,i].-Lambda_2*nu)'*inv_lemma_2*(Y[:,i].-Lambda_2*nu)))/(sqrt((2 .* pi)^(p) .*det_lemma_2))).-(kappa[i] .*((exp((-1/2) *(Y[:,i].-Lambda_2*nu)'*inv_lemma_2*(Y[:,i].-Lambda_2*nu)))/(sqrt((2 .* pi)^(p) .*det_lemma_2))))))
#                 end
#                 P1p = P1p .+ 1e-7
#                 P1p[findall(P1p.>1)].=0.9999999999999
#                 P0 = zeros(n)
#                 P1p[findall(P1p.>1)].=0.9999999999999
#                 for j=1:n
#                     P0[j] = (1-P1p[j])*log(1-P1p[j])
#                     P1[j] = P1p[j]*log(P1p[j])
#                 end
#                 P0[isnan.(P0)].=0.0
#                 Ek = - (sum(P1) + sum(P0))
#                 Entropy[i] = 1-(Ek/(n*log(2)))

#                 Ez1s[:,i] = res[((p*q)+(p*K)+p+p+q^2+n+1+1):((p*q)+(p*K)+p+p+q^2+n+1+n)]
#                 betss[:,i] = res[((p*q)+(p*K)+p+p+q^2+n+1+n+1):((p*q)+(p*K)+p+p+q^2+n+1+n+ncov)]
#                 liks[i] = res[((p*q)+(p*K)+p+p+q^2+n+1):((p*q)+(p*K)+p+p+q^2+n+1)][1]
#             end

#             # Ez1 = round.(Int64,Ez1s[:,findall((liks).==maximum(liks))])[:,1]
#             # bets = betss[:,findall((liks).==maximum(liks))]
#             if sum((Entropy).==maximum(Entropy))==1
#                 Ez1 = round.(Int64,Ez1s[:,findall((Entropy).==maximum(Entropy))])[:,1]
#                 bets = betss[:,findall((Entropy).==maximum(Entropy))]
#             elseif sum((Entropy).==maximum(Entropy))>1
#                 Ezone = Ez1s[:,findall((Entropy).==maximum(Entropy))]
#                 betsone = betss[:,findall((Entropy).==maximum(Entropy))]
#                 likone = liks[findall((Entropy).==maximum(Entropy))]
                
#                 Ez1 = round.(Int64,Ezone[:,findall(likone.==maximum(likone))])[:,1]
#                 bets = betss[:,findall(likone.==maximum(likone))]
#             end


#             if sum(Ez1)/n<0.5
#                 Ez1 = 1 .- Ez1
#                 bets = .- bets
#                 println("label unswitched: ",bets)
#             end

#             # 2. CFA ordinal PML
#             if c == 4
#                 trcfainit = rand(Uniform(-1.8,-0.8),p)
#                 trcfainit = hcat(trcfainit,zeros(p),.-trcfainit)
#             elseif c > 4
#                 trcfainit = rand(Uniform(-1.9,-1.6),p).+reshape(repeat(collect(range(0,3.5,c-1)),inner=p),p,(c-1))
#             end
            
#             paramcfa = vcat(vec(trcfainit),rand(Uniform(0.1,1),p),rand(Uniform(1e-3,1),Int(q*(q-1)/2)))
#             lb=vec(vcat(fill(-5,p*(c-1)),fill(0.001,p),fill(1e-3,Int((q*(q-1))/2))))#,fill(0.0001,p)));fill(-5,p*(c-1)).+repeat(collect(0:(c-2)),inner=p)
#             ub=vec(vcat(fill(5,p*(c-1)),fill(1.,p),fill(1.,Int((q*(q-1))/2))))#,fill(4,p))); vcat(fill(5,p*(c-1)).-repeat(collect((c-2):-1:0),inner=p)
#             parcfa = optimize(param ->pmlcfa(param,p,q,c,X[Bool.(Ez1),:],L_str),lb,ub,paramcfa,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)
    
#             # 3. EFA ordinal PML
#             if c > 4
#                 trefainit = rand(Uniform(1.2,2),p).+reshape(repeat(collect(range(0,2,c-1-2)),inner=p),p,(c-1-2))
#             elseif c == 4
#                 trefainit = rand(Uniform(1.4,2.5),p)
#             end
            
#             paramefa = vcat(vec(trefainit),rand(Uniform(0.1,1.),Int(p*K-(K*(K-1)/2))),rand(Uniform(0.1,3),p),rand(Uniform(0.5,2),p))# sort(rand(Uniform(1.2,4),p*(c-1-2)).+repeat(collect(0:(c-2-2)),inner=p))
#             lb=vec(vcat(fill(1,p*(c-1-2)),fill(1e-3,Int((p*K-(K*(K-1)/2)))),fill(0.001,p),fill(-10,p)))#,fill(0.0001,p))); .+repeat(collect(0:(c-2-2)),inner=p)
#             ub=vec(vcat(fill(7,p*(c-1-2)),fill(1.,Int((p*K-(K*(K-1)/2)))),fill(3.,p),fill(10.,p)))#,fill(4,p)));  .-repeat(collect((c-2-2):-1:0),inner=p)
#             parefa = optimize(param ->pmlefa(param,p,K,c,X[Bool.(1 .-Ez1),:],indices),lb,ub,paramefa,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

#             # in case of very strong straightlining
#             # if parefa.minimizer[((c-1-2)*p+1):Int((c+1+2)*p+(p*K-(K*(K-1)/2)))] .== ones(Float64,Int(p*K-(K*(K-1)/2)))
#             #     parefa.minimizer[((c-1-2)*p+1):Int((c+1+2)*p+(p*K-(K*(K-1)/2)))]  .= rand(Uniform(0.8,0.98),Int(p*K-(K*(K-1)/2)))
#             # end
#             # if parefa.minimizer[Int((c+1+2)*p+(p*K-(K*(K-1)/2))+1):Int((c+1+2)*p+(p*K-(K*(K-1)/2))+p)] .== fill(0.01,p)
#             #     parefa.minimizer[Int((c+1+2)*p+(p*K-(K*(K-1)/2))+1):Int((c+1+2)*p+(p*K-(K*(K-1)/2))+p)]  .= rand(Uniform(0.02,0.2),p)
#             # end
        
            
#             tau1=parcfa.minimizer[(1):((c-1)*p)]
#             tau1=vec(sort(reshape(tau1,p,c-1),dims=2))
#             lambda1=parcfa.minimizer[((c-1)*p+1):((c-1)*p+p)]
#             phi=parcfa.minimizer[((c-1)*p+p+1):((c-1)*p+p+Int(q*(q-1)/2))]
#             lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
#             lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
#             phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
#             tau2=parefa.minimizer[(1):((c-1-2)*p)]
#             tau2=vec(sort(reshape(tau2,p,c-1-2),dims=2))
#             tau2[tau2.<=1.01].=rand(Uniform(1.05,3),sum(tau2.<=1.01))
#             lambda2=parefa.minimizer[((c-1-2)*p+1):((c-1-2)*p+Int((p*K-(K*(K-1)/2))))]
#             psi=parefa.minimizer[((c-1-2)*p+Int((p*K-(K*(K-1)/2)))+1):((c-1-2)*p+Int((p*K-(K*(K-1)/2)))+p)]
#             interc=parefa.minimizer[((c-1-2)*p+Int((p*K-(K*(K-1)/2)))+p+1):((c-1-2)*p+Int((p*K-(K*(K-1)/2)))+p+p)]
#             lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
#             lambda2[lambda2.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.99))
#             psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
            
#             param_inits = vcat(tau1,lambda1,phi,tau2,lambda2,psi,interc,bets)
        
#             # # 2. CFA ordinal PML
#             # paramcfa = vcat(sort(rand(Uniform(-4.9,-2.3),p*(c-1))).+repeat(collect(0:(c-2)),inner=p),rand(Uniform(0.11,0.9),p),rand(Uniform(-0.9,0.9),Int(q*(q-1)/2)))
#             # lb=vec(vcat(fill(-5,p*(c-1)).+repeat(collect(0:(c-2)),inner=p),fill(0.1001,p),fill(-0.9,Int((q*(q-1))/2))))#,fill(0.0001,p)));
#             # ub=vec(vcat(fill(3,p*(c-1)).-repeat(collect((c-2):-1:0),inner=p),fill(0.9,p),fill(0.9,Int((q*(q-1))/2))))#,fill(4,p)));
#             # parcfa = optimize(param ->pmlcfa(param,p,q,c,X[Bool.(Ez1),:],L_str),lb,ub,paramcfa,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

#             # # 3. EFA ordinal PML
#             # paramefa = vcat(sort(rand(Uniform(1.2,4),p*(c-1-2)).+repeat(collect(0:(c-2-2)),inner=p)),rand(Uniform(0.11,0.9),Int(p*K-(K*(K-1)/2))),rand(Uniform(0.1,2.9),p),rand(Uniform(-9.90,9.90),p)) 
#             # lb=vec(vcat(fill(1.1,p*(c-1-2)).+repeat(collect(0:(c-2-2)),inner=p),fill(0.10001,Int((p*K-(K*(K-1)/2)))),fill(0.101,p),fill(-9.90,p)))#,fill(0.0001,p)));
#             # ub=vec(vcat(fill(6,p*(c-1-2)).-repeat(collect((c-2-2):-1:0),inner=p),fill(0.9,Int((p*K-(K*(K-1)/2)))),fill(2.9,p),fill(9.90,p)))#,fill(4,p)));
#             # parefa = optimize(param ->pmlefa(param,p,K,c,X[Bool.(1 .-Ez1),:]),lb,ub,paramefa,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)


#             # param_init = vcat(parcfa.minimizer.+rand(Uniform(-0.1,0.1),size(parcfa.minimizer,1)),parefa.minimizer.+rand(Uniform(-0.1,0.1),size(parefa.minimizer,1)),sum(Ez1)/n)
#             println(round.(param_inits,digits=3),sum(Ez1s[:,findall((liks).==minimum(liks))])/n,sum(Ez1)/n)
            

#             # replis = 10
#             # paramss = zeros(Float64,size(param_inits,1),replis); likss = zeros(Float64,replis) 
#             # for i in 1:replis
#             #     println(i)
#             #     plimix = E_stepMix(param_inits,q,K,p,n,c,X,pairsindex,L_str,true);
#             #     T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
#             #     t = time()
#             #     paramss[:,i] .= M_parsMix(T,param_inits,X,n,p,q,K,c,pairsindex,L_str,indices,1,t,plimix,true); #solo primo argomento di output
#             #     likss[i] = E_stepMix(paramss[:,i],q,K,p,n,c,X,pairsindex,L_str,false,T);
#             # end 
#             # println(likss,likss[findall((likss).==maximum(likss))],likss[findall((likss).==minimum(likss))])
#             # paramsss = paramss[:,findall((likss).==maximum(likss))]
#             # println("started")

#             # first iteration
#             plimix = E_stepMix(param_inits,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
#             T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
#             t = time()
#             params[:,i] = M_parsMix(T,param_inits,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
#             dt = time() - t
#             liksmix[i] = E_stepMix(params[:,i],q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
#         end
#         param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
#         lik=liksmix[findall(liksmix.==maximum(liksmix))]
    
#         tau1=param_initials0[(1):((c-1)*p)]
#         tau1=vec(sort(reshape(tau1,p,c-1),dims=2))
#         lambda1=param_initials0[((c-1)*p+1):((c-1)*p+p)]
#         phi=param_initials0[((c-1)*p+p+1):((c-1)*p+p+Int(q*(q-1)/2))]
#         lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
#         lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
#         phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
#         tau2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p)]
#         tau2=vec(sort(reshape(tau2,p,c-1-2),dims=2))
#         tau2[tau2.<=1.01].=rand(Uniform(1.05,3),sum(tau2.<=1.01))
#         lambda2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2))))]
#         psi=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+p)]
#         interc=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+p+p)]
#         lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
#         lambda2[lambda2.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.99))
#         psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))

function startvals(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
    liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)*p+p+(q*(q-1))/2+(c-1-2)*p+(p*K-(K*(K-1)/2))+(p)+p)+ncov,replix);#p-1
    for i in 1:replix
        # println(i)

        tau1= vec(hcat(rand(Uniform(-2.3,-1.),p).+reshape(repeat(collect(range(0,3.5,c-1)),inner=p),p,(c-1))))
        lambda1 = rand(Uniform(0.5,.90),p);
        phi=rand(LKJ(q,1))
        phi=phi[tril!(trues(size(phi)), -1)]
        tau2=vec(reshape(rand(Uniform(1.2,4),p*(c-1-2)),p,(c-1-2)).+reshape(repeat(collect(0:(c-2-2)),inner=p),p,(c-1-2)))
        tau2=vec(sort(reshape(tau2,p,c-1-2),dims=2))
        lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(0.15,.90),Int((p*K-(K*(K-1)/2))))))
        lambda2=lamb2[lamb2.!=0.]
        psi=rand(Uniform(0.5,4.5),p)
        interc=rand(Uniform(0.8,3),p)



        prop = rand(Uniform(0.65,0.95))
        bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.1,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


        param_initials = vcat(tau1,lambda1,phi,tau2,lambda2,psi,interc,bets)

        # first iteration
        plimix = E_stepMix(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        params[:,i] = M_parsMix(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        dt = time() - t
        prrot=params[:,i]
        liksmix[i] = E_stepMix(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
    end
    param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    lik=liksmix[findall(liksmix.==maximum(liksmix))]

    tau1=param_initials0[(1):((c-1)*p)]
    tau1=vec(sort(reshape(tau1,p,c-1),dims=2))
    lambda1=param_initials0[((c-1)*p+1):((c-1)*p+p)]
    phi=param_initials0[((c-1)*p+p+1):((c-1)*p+p+Int(q*(q-1)/2))]
    lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    tau2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p)]
    tau2=vec(sort(reshape(tau2,p,c-1-2),dims=2))
    tau2[tau2.<=1.01].=rand(Uniform(1.05,3),sum(tau2.<=1.01))
    lambda2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2))))]
    psi=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p))]# psi=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1))]
    interc=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p)+p)]#interc=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p)+p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p)+p+ncov)]#    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


    param_initialii = vcat(tau1,lambda1,phi,tau2,lambda2,psi,interc,bets)

    # println(param_initialii,liksmix,"first iteration: done")
    return param_initialii
end


function EMparsMix(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
    local plimix::Array, T,oldlik
    n = size(X,1)
    er=0;
    # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
    it=0;
    dev=Inf; dev1=Inf; dev2=Inf;

    
    oldlik=-1000000.;
    oldlik2=1000000.;
    oldlikA=-1000000.;
    lik=0.
    

    noerror = true; devnostop = true;  sum_devs = 0;
    while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
        it = it+1;
        param=param01;

        T00 = 1
        plimix = E_stepMix(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        param01 = M_parsMix(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
        dt = time() - t
        if param01 == false
            return return(vcat(param01,vec(plimix),vec(T),0.,0.))
        end
        lik = E_stepMix(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);

        #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
        con = (lik-oldlik)/(oldlik-oldlik2);
        likA = oldlik + (lik-oldlik)/(1-con)
        old_dev = dev1
        dev1 = abs(likA-oldlikA)
        dev2 = abs(lik-oldlik)
        # println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
        # println(param01)
        dev = abs(lik-oldlik)/abs(oldlik)

        noerror = oldlik<lik

        oldlikA=likA
        oldlik2=oldlik
        oldlik=lik;
        pml_em_log("EMparsMix", it, lik, dev, dev1, abs(lik - oldlik2), tol, noerror)

        # sum_devs += dev > old_dev
        # devstop = sum_devs <= 3
    end
    
    # classificazioni dopo ultima massimizzazione
    plimix = E_stepMix(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

    if noerror == true && it<=maxit 
        return vcat(param01,vec(plimix),vec(T),lik,1.)
    else
        return vcat(param01,vec(plimix),vec(T),0.,0.)
    end
end

function reconstruct_correlation_matrix(corr_vector::AbstractVector{T}) where T <: Real
    n = Int((1 + sqrt(1 + 8 * length(corr_vector))) ÷ 2)
    corr_matrix = Matrix{T}(I, n, n)
    k = 1
    for i in 1:n-1
        for j in i+1:n
            value = corr_vector[k]
            corr_matrix[i, j] = value
            corr_matrix[j, i] = value
            k += 1
        end
    end
    return corr_matrix
end

function reconstruct_L2_matrix(p::Int, K::Int, indices::Vector{Int}, elements::AbstractVector{T}) where T <: Real
    # Initialize the matrix with zeros of type T
    matrix = zeros(T, p, K)
    
    # Place the elements in the matrix according to the indices
    for (i, idx) in enumerate(indices)
        row = (idx - 1) % p + 1
        col = div(idx - 1, p) + 1
        matrix[row, col] = elements[i]
    end
    return matrix
end




function E_stepMix(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
    if sum(param .== 0.) == size(param,1)
        return -Inf
    end
    trs = (c-1-2)
    trsi = (c-1)
    t1s = param[1:trsi*p]
    L1 = param[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    
    t2s = param[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs*p)]
    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi*p+p+(q*(q-1))/2+trs*p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))] )

    
    Psi = diagm(param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    interc2 = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+p)]   # interc2 = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p)]
    betas = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+p+ncov)]  #  betas = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)
    
    Sigms=ones(p,p,2);
    Sigms[:,:,1]=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sigms[:,:,2]=L2*L2'+Psi 
    
    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
    br2 = hcat(zeros(Float64,p),ones(Float64,p),reshape(t2s,p,trs))
    t2=hcat(fill(-Inf,p),br2,fill(Inf,p))

    tcfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tcfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
    tefa1=[t2[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tefa=[t2[j,X[i,j]] for i in 1:n, j in 1:p]

    tlow = cat(tcfa,tefa,dims=3)
    tup = cat(tcfa1,tefa1,dims=3)
    
    if estep == true
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));
    
        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tup[:,pairsindex[:,idx],g])
                    low=(tlow[:,pairsindex[:,idx],g])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tup[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                    low=(tlow[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3.*pg[:,g];
            end
        end
    # % probabilities in each row must sum to 1 (application of Bayes law)
        return plimix
    elseif estep == false
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));

        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tup[:,pairsindex[:,idx],g])
                    low=(tlow[:,pairsindex[:,idx],g])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tup[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                    low=(tlow[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3
            end
        end
        lik1=(log.(plimix));
        lik1=sum(lik1,dims=3);
        U=exp.(lik1)[:,:,1].*pg; 
        cllik=sum(log.(sum(U,dims=2)))
        return cllik
    end
end
    

    


function M_parsMix(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)

    trsi = (c-1)*p #attenzione p è già qui
    trs = (c-1-2)*p

    # PML - CFA & EFA
    theta_init = theta[1:(size(theta,1)-ncov)]#0.000001
    # # -1 loadings
    # lb=vec(vcat(fill(-5.1,trsi),fill(0.00001,p),fill(-1.,Int((q*(q-1))/2)),fill(1.,trs),fill(-1.,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-10,p)));
    # # 0. loadings
    lb=vec(vcat(fill(-5.1,trsi),fill(-1,p),fill(-1.,Int((q*(q-1))/2)),fill(1.,trs),fill(-Inf,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-5,p)));#1.0e-5
    ub=vec(vcat(fill(5.1,trsi),fill(1.,p),fill(1.,Int((q*(q-1))/2)),fill(7.,trs),fill(Inf,Int(p*K-K*(K-1)/2)),fill(Inf,p),fill(10.,p)));

    par = zeros(Float64,size(lb))
    if it == 1
        par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min ,Optim.Options(iterations=1000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0))
    elseif it > 1
        par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min
    end 


    # par = optimize(x->pmllast(x,X,T,p,q,K,c,n,pairsindex,L_str,indices),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

    # par = zeros(Float64,size(lb))
    # if it == 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 300.0)) #20 min
    # elseif it > 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 1200.0)) #20 min
    # end 
    
    # mixture parameter
    betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+p+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
    Ez1 = T[:,1,:]
    Ez0 = T[:,2,:]
    
    pairs = size(pairsindex,2)
    betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)

    
    # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
   
    # println(par.minimizer,pg[1])
    if !(Optim.converged(par)) && isnan(betas[1]) 
        println(betas[1])
        if start == true
            return zeros(Float64,size(theta,1))
        elseif start == false
            return Optim.converged(par)
        end
    end
    return vcat(par.minimizer,betasnew)#,pg[1])
end

function betas_gradient(betas::Vector,ics::Matrix{Float64},Ez1s::Vector{Float64},Ez0s::Vector{Float64},pairs)
    fval = - ics'*((-Ez1s/pairs .+Ez0s/pairs .*exp.(ics*betas))./(exp.(ics*betas).+1))#( sum((-Ez1.+Ez0.*exp.(X*betas)).*X./(exp.(X*betas).+1), dims=1)) #(.- sum(Ez1.*X./(exp.(X*betas).+1), dims=1) .+ sum(Ez0.*X.*exp.(X*betas)./(exp.(X*betas).+1),dims=1))
    return vec(fval)
end


function NewtonRaphson(betas_gradient, betas::Vector, ics::Matrix{Float64}, Ez1::Matrix{Float64}, Ez0::Matrix{Float64}, ncov::Int64,pairs)
    local betas_old::Vector{Float64}, iter::Int64, tol::Float64, maxIter::Int64, grad::Vector{Float64}, hess::Matrix{Float64}
    betas_old = fill(-Inf,ncov); iter = 0;  tol=1e-8; maxIter = 1000; Ez1s = sum(Ez1,dims=2)[:,1];  Ez0s = sum(Ez0,dims=2)[:,1];
    while sum(abs.(betas.-betas_old)) > tol && iter < maxIter
        grad = betas_gradient(betas,ics,Ez1s,Ez0s,pairs)
        hess = -(ics.*exp.(ics*betas)./(exp.(ics*betas).+1).^2)'*ics
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

function betas_gradient(betas::Vector,ics::Matrix{Float64},Ez1::Matrix{Float64},Ez0::Matrix{Float64},pairs)
    fval = - ics'*((-sum(Ez1,dims=2)/pairs .+sum(Ez0,dims=2)/pairs .*exp.(ics*betas))./(exp.(ics*betas).+1))#( sum((-Ez1.+Ez0.*exp.(X*betas)).*X./(exp.(X*betas).+1), dims=1)) #(.- sum(Ez1.*X./(exp.(X*betas).+1), dims=1) .+ sum(Ez0.*X.*exp.(X*betas)./(exp.(X*betas).+1),dims=1))
    return vec(fval)
end


# function pmlmix(param,X,T,p,q,K,c,n,pairsindex,L_str,indices)
#     # local trs::Int64, trsi::Int64, plimix1::Vector{Float64}, plimix2::Vector{Float64}, Tt1::Vector{Float64}, Tt2::Vector{Float64}, B::Matrix{Float64}
#     trs = (c-1-2)
#     trsi = (c-1)
#     t1 = param[1:trsi*p]
#     L1 = param[(trsi*p+1):(trsi*p+p)].*L_str 
#     Phi = reconstruct_correlation_matrix(param[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
   
    
#     t2 = param[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs*p)]
#     L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi*p+p+(q*(q-1))/2+trs*p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))] )
#     Psi = diagm(param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+p)])#K
#     interc2 = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p)]
   
#     br1 = reshape(t1s,p,trsi)
#     t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
#     br2 = hcat(zeros(Float64,p),ones(Float64,p),reshape(t2s,p,trs))
#     t2=hcat(fill(-Inf,p),br2,fill(Inf,p))

#     Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
#     Sig2=L2*L2'+Psi

#     plimix1=zeros(Float64,n);plimix2=zeros(Float64,n);
#     tcfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     tcfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
#     tefa1=[t2[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     tefa=[t2[j,X[i,j]] for i in 1:n, j in 1:p]
   
   
#     low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#     for idx=1:size(pairsindex,2) #parallelizzabile
#         for g=1:2
#             if g==1
#                 s2c=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
#                 # s2c=s2c+0.0001*Matrix(I,2,2)
#                 r2c=cov2cor(s2c);
#                 B=hcat(tcfa[:,pairsindex[:,idx]],tcfa1[:,pairsindex[:,idx]]);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 Tt1 = T[:,1,idx]
#                 plimix1 += log.(t3).*Tt1;
#             else 
#                 s2c=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
#                 # s2c=s2c+0.0001*Matrix(I,2,2)
#                 r2c=cov2cor(s2c);
#                 B=hcat((tefa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(s2c)'),(tefa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(s2c)'));
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 Tt2 = T[:,2,idx]
#                 plimix2 += log.(t3).*Tt2;
#             end
            
#         end
#     end
#     f=-sum(plimix1+plimix2)
#     # @info "Current objective: $f"
#     return f
# end


function pmllast(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    trs = (c-1-2)
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])


    t2s = x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs*p)]
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+trs*p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))] )

    Psi = diagm(x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    interc2 = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+p)]    #interc2 = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p)]
    # prop = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]

    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
    br2 = hcat(zeros(Float64,p),ones(Float64,p),reshape(t2s,p,trs))
    t2=hcat(fill(-Inf,p),br2,fill(Inf,p))


    Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi



    plimix=zeros(Float64,n);
    tcfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tcfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
    tefa1=[t2[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tefa=[t2[j,X[i,j]] for i in 1:n, j in 1:p]
    
    for idx=1:size(pairsindex,2) #parallelizzabile
        # CFA
        Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
        Bcfa=hcat(tcfa[:,pairsindex[:,idx]],tcfa1[:,pairsindex[:,idx]]);
        tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
        Tt1 = T[:,1,idx]
        
        #EFA
        Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
        Refa=cov2cor(Sefa);
        Befa=hcat((tefa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tefa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
        tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
        Tt2 = T[:,2,idx]

        plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
    end
    f=-sum(plimix)
return f
end

function fg!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trs = (c-1-2)
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
   

    t2s = x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs*p)]
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+trs*p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))] )

    Psi = diagm(x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    interc2 = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+p)]    #interc2 = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p)]
    # prop = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]

    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
    br2 = hcat(zeros(Float64,p),ones(Float64,p),reshape(t2s,p,trs))
    t2=hcat(fill(-Inf,p),br2,fill(Inf,p))


    Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi


    if G !== nothing
        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradintercefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)
    
        gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)
    
        gradlambci = zeros(Float64,2,q,size(pairsindex,2))
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
        
        idx = 0
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1
                rho1 = Sig1[i,j]
                cormat1 = ones(Float64,2,2); cormat1[1,2]=rho1; cormat1[2,1]=rho1
                sigs = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sigs)
                rho2 = cormat2[1,2]
                muij = interc2[vcat(i,j)]
                sigsij = diag(sigs)
                sigij = sqrt.(sigsij)
                sigij13 = sigij[1]^(3)
                sigij23 = sigij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]
                @inbounds for N in 1:n
                    xni = X[N,i]
                    xnj = X[N,j]
                    t1i = t1[i,xni]
                    t1i1 = t1[i,xni+1]
                    t1j = t1[j,xnj]
                    t1j1 = t1[j,xnj+1]
                    T1 = T[N,1,idx][1]
                    T2 = T[N,2,idx][1]
                    t2i = t2[i,xni]
                    t2i1 = t2[i,xni+1]
                    t2j = t2[j,xnj]
                    t2j1 = t2[j,xnj+1]
                    t2ist = (t2i-muij[1])/sigij[1]
                    t2i1st = (t2i1-muij[1])/sigij[1]
                    t2jst = (t2j-muij[2])/sigij[2]
                    t2j1st = (t2j1-muij[2])/sigij[2]
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                   
                    pijicfa = bvnP(t1i, t1i1, t1j, t1j1, rho1, false)
                    gradtaucfa[xni+1,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i1)*(cdf(Normal(0,1),(t1j1-rho1*t1i1)/sqrtrij1)-cdf(Normal(0,1),(t1j-rho1*t1i1)/sqrtrij1))) 
                    gradtaucfa[xni,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i)*(cdf(Normal(0,1),(t1j-rho1*t1i)/sqrtrij1)-cdf(Normal(0,1),(t1j1-rho1*t1i)/sqrtrij1))) 
                    gradtaucfa[xnj+1,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j1)*(cdf(Normal(0,1),(t1i1-rho1*t1j1)/sqrtrij1)-cdf(Normal(0,1),(t1i-rho1*t1j1)/sqrtrij1))) 
                    gradtaucfa[xnj,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j)*(cdf(Normal(0,1),(t1i-rho1*t1j)/sqrtrij1)-cdf(Normal(0,1),(t1i1-rho1*t1j)/sqrtrij1))) 
                    
    
                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))
    
    
                    
                    if !(t1i == -Inf && t1j == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1 == Inf && t1j1 == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j1))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end
        
                    if !(t1i1 == Inf && t1j == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end
        
                    if !(t1i == -Inf && t1j1 == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j1))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    gradlambci[1:2,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)
                    gradphici[1:q,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*(L1[i,:]*L1[j,:]')
                    
                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sigij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sigij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sigij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sigij[2])
    
                    gradintercefa[xni+1,idx,1] += (derivNormi1/(-sigij[1]))
                    gradintercefa[xnj+1,idx,2] += (derivNormj1/(-sigij[2]))
                    gradintercefa[xni,idx,1] += (derivNormi/(-sigij[1]))
                    gradintercefa[xnj,idx,2] += (derivNormj/(-sigij[2]))
    
                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)
    
                    
                    gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sigij13)*(sigij[2])))
                    gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sigij23)*(sigij[1])))
    
                    gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij[1]-t2i1))/(2*(sigij13))
                    gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij[2]-t2j1))/(2*(sigij23))
                    gradpsiefa[xni,idx,1] += (derivNormi*(muij[1]-t2i))/(2*(sigij13))
                    gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij[2]-t2j))/(2*(sigij23))
                    
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sigij13*sigij[2])+L2jk/((sigij[1]*sigij[2])))#((L2j*Psi[i,i])/(sigij13*sigij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sigij23*sigij[1])+L2ik/((sigij[2]*sigij[1])))#((L2[i,:].*Psi[j,j])./(sigij23*sigij[1]))
    
                        gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(((muij[1].-t2i1).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(((muij[2].-t2j1).*L2jk)./((sigij23)))
                        gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(((muij[1].-t2i).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(((muij[2].-t2j).*L2jk)./((sigij23)))
                    end 
                
                end
            end
        end
    
    
        vecpair = vec(pairsindex)
        liktau = zeros(p*(c-1));idpcfa = 0   
        @inbounds for j in 2:c
            @inbounds for i in 1:p
                idpcfa += 1
                liktau[idpcfa] = -sum(vec(vcat(gradtaucfa[j,:,1]',gradtaucfa[j,:,2]'))[vecpair.==i])
            end
        end
    
        liktau2 = zeros(p*(c-1-2));idpefa = 0   
        @inbounds for j in (1+1+2):c
            @inbounds for i in 1:p
                idpefa += 1
                liktau2[idpefa] = -sum(vec(vcat(gradtauefa[j,:,1]',gradtauefa[j,:,2]'))[vecpair.==i])
            end
        end
        
    
        llamb1 = permutedims(gradlambci,(2,3,1)); llam1 = zeros(p); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
        likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
        likinterc2 = zeros(p); sumgradinter = sum(gradintercefa[2:c,:,:],dims=1)[1,:,:]
        likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        @inbounds for i in 1:p
            llam1[i] = -sum(vec(llamb1[ids[i],:,:]')[vecpair.==i])
            likinterc2[i] = -sum(vec(sumgradinter')[vecpair.==i])
            likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
            likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end
    
        likpsitrue = likpsi.+likpsiidx
        liklstrue = likls.+liklsidx
    
        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]
    
        G .= vec(vcat(liktau,llam1,likphida,liktau2,liklstrue[indices],likpsitrue,likinterc2))#likpsitrue[2:end]
      
    end
    if F !== nothing 
  
    
        plimix=zeros(Float64,n);
        tcfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
        tcfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
        tefa1=[t2[j,X[i,j]+1] for i in 1:n, j in 1:p]
        tefa=[t2[j,X[i,j]] for i in 1:n, j in 1:p]
       
       
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
            Bcfa=hcat(tcfa[:,pairsindex[:,idx]],tcfa1[:,pairsindex[:,idx]]);
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
            Tt1 = T[:,1,idx]
            
            #EFA
            Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tefa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tefa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            Tt2 = T[:,2,idx]
            plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2;
            
        end
        f=-sum(plimix)
        # @info "Current objective: $f"
        return f
    end
end


########################################################################### Latent means EFA - Millsap's parametrization ##################################################################################################

function startvalslm(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
    liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)*p+p+(q*(q-1))/2+(c-1-2)*p+(p*K-(K*(K-1)/2))+(p)+K)+ncov,replix);#p-1
    for i in 1:replix
        # println(i)

        tau1 = vec(hcat(rand(Uniform(-1.9,-1.),p).+reshape(repeat(collect(range(0,3.1,c-1)),inner=p),p,(c-1))))
        lambda1 = rand(Uniform(0.5,1.3),p);
        phi=rand(LKJ(q,1))
        phi=phi[tril!(trues(size(phi)), -1)]
        tau2=vec(reshape(rand(Uniform(0.8,3),p*(c-1-2)),p,(c-1-2)).+reshape(repeat(collect(0:(c-2-2)),inner=p),p,(c-1-2)))
        tau2=vec(sort(reshape(tau2,p,c-1-2),dims=2))
        lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(-0.90,.90),Int((p*K-(K*(K-1)/2))))))
        lambda2=lamb2[lamb2.!=0.]
        psi=rand(Uniform(0.2,2.5),p)
        kappa=rand(Uniform(-1.2,1.2),K)



        prop = rand(Uniform(0.65,0.97))
        bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.1,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


        param_initials = vcat(tau1,lambda1,phi,tau2,lambda2,psi,kappa,bets)

        # first iteration
        plimix = E_stepMixlm(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        params[:,i] = M_parsMixlm(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        dt = time() - t
        prrot=params[:,i]
        liksmix[i] = E_stepMixlm(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false);
    end
    param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    lik=liksmix[findall(liksmix.==maximum(liksmix))]

    tau1=param_initials0[(1):((c-1)*p)]
    tau1=vec(sort(reshape(tau1,p,c-1),dims=2))
    lambda1=param_initials0[((c-1)*p+1):((c-1)*p+p)]
    phi=param_initials0[((c-1)*p+p+1):((c-1)*p+p+Int(q*(q-1)/2))]
    # lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    # lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    tau2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p)]
    tau2=vec(sort(reshape(tau2,p,c-1-2),dims=2))
    tau2[tau2.<=1.01].=rand(Uniform(1.05,3),sum(tau2.<=1.01))
    lambda2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2))))]
    psi=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p))]# psi=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1))]
    kappa=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p)+K)]#interc=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    # lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    # lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p)+K+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p)+K+ncov)]#    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]



    param_initialii = vcat(tau1,lambda1,phi,tau2,lambda2,psi,kappa,bets)

    # println(param_initialii,liksmix,"first iteration: done")
    return param_initialii
end


    # liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)*p+p+(q*(q-1))/2+(c-1-2)*p+(p*K-(K*(K-1)/2))+(p)+K)+ncov,replix);#p-1
    # for i in 1:replix
    #     println(i)

    #     tau1 = vec(hcat(rand(Uniform(-1.9,-1.),p).+reshape(repeat(collect(range(0,3.1,c-1)),inner=p),p,(c-1))))
    #     lambda1 = rand(Uniform(0.5,.90),p);
    #     phi=rand(LKJ(q,1))
    #     phi=phi[tril!(trues(size(phi)), -1)]
    #     tau2=vec(reshape(rand(Uniform(0.8,3),p*(c-1-2)),p,(c-1-2)).+reshape(repeat(collect(0:(c-2-2)),inner=p),p,(c-1-2)))
    #     tau2=vec(sort(reshape(tau2,p,c-1-2),dims=2))
    #     lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(0.15,.90),Int((p*K-(K*(K-1)/2))))))
    #     lambda2=lamb2[lamb2.!=0.]
    #     psi=rand(Uniform(0.2,1.5),p)
    #     kappa=rand(Uniform(-1.2,1.2),K)



    #     prop = rand(Uniform(0.65,0.95))
    #     bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.1,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


    #     param_initials = vcat(tau1,lambda1,phi,tau2,lambda2,psi,kappa,bets)

    #     # first iteration
    #     plimix = E_stepMixlm(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    #     T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
    #     t = time()
    #     params[:,i] = M_parsMixlm(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
    #     dt = time() - t
    #     prrot=params[:,i]
    #     liksmix[i] = E_stepMixlm(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false);
    # end
    # param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    # lik=liksmix[findall(liksmix.==maximum(liksmix))]

    # tau1=param_initials0[(1):((c-1)*p)]
    # tau1=vec(sort(reshape(tau1,p,c-1),dims=2))
    # lambda1=param_initials0[((c-1)*p+1):((c-1)*p+p)]
    # phi=param_initials0[((c-1)*p+p+1):((c-1)*p+p+Int(q*(q-1)/2))]
    # lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    # lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    # phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    # tau2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p)]
    # tau2=vec(sort(reshape(tau2,p,c-1-2),dims=2))
    # tau2[tau2.<=1.01].=rand(Uniform(1.05,3),sum(tau2.<=1.01))
    # lambda2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2))))]
    # psi=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p))]# psi=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1))]
    # kappa=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p)+K)]#interc=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    # lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    # lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    # psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
    # bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p)+K+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p)+K+ncov)]#    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


function EMparsMixlm(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
    local plimix::Array, T,oldlik
    n = size(X,1)
    er=0;
    # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
    it=0;
    dev=Inf; dev1=Inf; dev2=Inf;

    
    oldlik=-1000000.;
    oldlik2=1000000.;
    oldlikA=-1000000.;
    lik=0.
    

    noerror = true; devnostop = true;  sum_devs = 0;
    while (dev>tol) && it <= maxit && noerror #&& devnostop
        
        it = it+1;
        param=param01;

        T00 = 1
        plimix = E_stepMixlm(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        param01 = M_parsMixlm(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
        dt = time() - t
        if param01 == false
            return return(vcat(param01,vec(plimix),vec(T),0.,0.))
        end
        lik = E_stepMixlm(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false);

        #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
        con = (lik-oldlik)/(oldlik-oldlik2);
        likA = oldlik + (lik-oldlik)/(1-con)
        old_dev = dev1
        dev1 = abs(likA-oldlikA)
        # dev2 = abs(lik-oldlik)
        
        # println(param01)
        dev = abs(lik-oldlik)#abs(lik-oldlik)/abs(oldlik)
        # println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)

        noerror = oldlik<lik

        oldlikA=likA
        oldlik2=oldlik
        oldlik=lik;
        pml_em_log("EMparsMixlm", it, lik, dev, dev1, abs(lik - oldlik2), tol, noerror)

        # sum_devs += dev > old_dev
        # devstop = sum_devs <= 3
    end
    
    # classificazioni dopo ultima massimizzazione
    plimix = E_stepMixlm(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

    if noerror == true && it<=maxit 
        return vcat(param01,vec(plimix),vec(T),lik,1.)
    else
        return vcat(param01,vec(plimix),vec(T),0.,0.)
    end
end





# function E_stepMixlm(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
#     if sum(param .== 0.) == size(param,1)
#         return -Inf
#     end
#     trs = (c-1-2)
#     trsi = (c-1)
#     t1s = param[1:trsi*p]
#     L1 = param[(trsi*p+1):(trsi*p+p)].*L_str 
#     Phi = reconstruct_correlation_matrix(param[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    
#     t2s = param[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs*p)]
#     L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi*p+p+(q*(q-1))/2+trs*p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))] )

    
#     Psi = diagm(param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p-1))]))#K
#     kappa = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+K)]   # interc2 = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p)]
#     betas = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+K+ncov)]  #  betas = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

#     propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
#     pg = hcat(propo,1. .-propo)

#     mu = L2*kappa
    
#     Sigms=ones(p,p,2);
#     Sigms[:,:,1]=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
#     Sigms[:,:,2]=L2*L2'+Psi 
    
#     br1 = reshape(t1s,p,trsi)
#     t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
#     br2 = hcat(br1[:,1:2],reshape(t2s,p,trs))
#     t2=hcat(fill(-Inf,p),br2,fill(Inf,p))


#     tcfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     tcfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
#     tefa1=[t2[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     tefa=[t2[j,X[i,j]] for i in 1:n, j in 1:p]

#     tlow = cat(tcfa,tefa,dims=3)
#     tup = cat(tcfa1,tefa1,dims=3)
    
#     if estep == true
#         low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#         plimix=ones(n,2,size(pairsindex,2));
    
#         for idx=1:size(pairsindex,2) #parallelizzabile
#             for g=1:2
#                 if g==1
#                     s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                     r2c=cov2cor(s2c);
#                     up=(tup[:,pairsindex[:,idx],g])
#                     low=(tlow[:,pairsindex[:,idx],g])
#                 elseif g==2
#                     s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                     r2c=cov2cor(s2c);
#                     up=(tup[:,pairsindex[:,idx],g].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                     low=(tlow[:,pairsindex[:,idx],g].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 end
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix[:,g,idx]= t3.*pg[:,g];
#             end
#         end
#     # % probabilities in each row must sum to 1 (application of Bayes law)
#         return plimix
#     elseif estep == false
#         low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#         plimix=ones(n,2,size(pairsindex,2));

#         for idx=1:size(pairsindex,2) #parallelizzabile
#             for g=1:2
#                 if g==1
#                     s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                     r2c=cov2cor(s2c);
#                     up=(tup[:,pairsindex[:,idx],g])
#                     low=(tlow[:,pairsindex[:,idx],g])
#                 elseif g==2
#                     s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                     r2c=cov2cor(s2c);
#                     up=(tup[:,pairsindex[:,idx],g].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                     low=(tlow[:,pairsindex[:,idx],g].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 end
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix[:,g,idx]= t3
#             end
#         end
#         lik1=(log.(plimix));
#         lik1=sum(lik1,dims=3);
#         U=exp.(lik1)[:,:,1].*pg; 
#         cllik=sum(log.(sum(U,dims=2)))
#         return cllik
#     end
# end
    


function E_stepMixlm(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep)#::AbstractVector{T},q::Int64,K::Int64,p::Int64,n::Int64,c::Int64,X::Matrix{Int64},ics::Matrix,ncov::Int64,pairsindex::Matrix{Int64},L_str::Matrix,indices,estep::Bool) where T
    # local trs::Int64,trsi::Int64,t1s::Vector{Real}, L1::Matrix{Real}, Phi::Matrix{Real},  t2s::Vector{Real}, L2::Matrix{Real}, Psi::Matrix{Real}, kappa::Vector{Real}, betas::Vector{Real}, propo::Vector{Real}, pg::Matrix{Real}, mu::Vector{Real}, Sigms1::Matrix{Real}, Sigms2::Matrix{Real}, br1::Matrix{Real}, br2::Matrix{Real}, t1::Matrix{Real}, t2::Matrix{Real}, tcfa1::Matrix{Real}, tcfa::Matrix{Real}, tefa1::Matrix{Real}, tefa::Matrix{Real}, plimix::Array{Real}, plimixll::Vector{Real}, Scfa::Matrix{Real}, Rcfa::Matrix{Real}, Bcfa::Matrix{Real}, tcfa3::Vector{Real}, Sefa::Matrix{Real}, Refa::Matrix{Real}, Befa::Matrix{Real}, tefa3::Vector{Real}, cllik::Real
    if sum(param .== 0.) == size(param,1)
        return -Inf
    end
    trs = (c-1-2)
    trsi = (c-1)
    t1s = param[1:trsi*p]
    L1 = param[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    
    t2s = param[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs*p)]
    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi*p+p+(q*(q-1))/2+trs*p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))] )

    
    Psi = diagm(param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    kappa = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+K)]   # interc2 = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p)]
    betas = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+K+ncov)]  #  betas = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)

    mu = L2*kappa
    
    # Sigms=zeros(p,p,2);
    # Deltas = Diagonal(Diagonal((L1*Phi*L1')) .+ 1);
    # Sigms[:,:,1]=(Deltas^(-1/2))*(L1*Phi*L1' + diagm(ones(p)))*(Deltas^(-1/2))
    Sigms1=L1*Phi*L1' + diagm(ones(p))
    Sigms2=L2*L2'+Psi 
    
    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
    br2 = hcat(br1[:,1:2],reshape(t2s,p,trs))
    t2=hcat(fill(-Inf,p),br2,fill(Inf,p))


    tcfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tcfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
    tefa1=[t2[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tefa=[t2[j,X[i,j]] for i in 1:n, j in 1:p]


    # tlow = cat(tcfa,tefa,dims=3)
    # tup = cat(tcfa1,tefa1,dims=3)
    
    # if estep == true
    #     low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
    #     plimix=ones(n,2,size(pairsindex,2));
    
    #     for idx=1:size(pairsindex,2) #parallelizzabile
    #         for g=1:2
    #             if g==1
    #                 s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
    #                 r2c=cov2cor(s2c);
    #                 up=(tcfa1[:,pairsindex[:,idx]]./sqrt.(diag(s2c)'))
    #                 low=(tcfa[:,pairsindex[:,idx]]./sqrt.(diag(s2c)'))
    #             elseif g==2
    #                 s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
    #                 r2c=cov2cor(s2c);
    #                 up=(tefa1[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
    #                 low=(tefa[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
    #             end
    #             B=hcat(low,up);
    #             t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
    #             plimix[:,g,idx]= t3.*pg[:,g];
    #         end
    #     end
    # # % probabilities in each row must sum to 1 (application of Bayes law)
    #     return plimix
    # elseif estep == false
        # low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        # plimix=ones(n,2,size(pairsindex,2));

        # for idx=1:size(pairsindex,2) #parallelizzabile
        #     for g=1:2
        #         if g==1
        #             s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
        #             r2c=cov2cor(s2c);
        #             up=(tcfa1[:,pairsindex[:,idx]]./sqrt.(diag(s2c)'))
        #             low=(tcfa[:,pairsindex[:,idx]]./sqrt.(diag(s2c)'))
        #         elseif g==2
        #             s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
        #             r2c=cov2cor(s2c);
        #             up=(tefa1[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
        #             low=(tefa[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
        #         end
        #         B=hcat(low,up);
        #         t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
        #         plimix[:,g,idx] = t3
        #     end
        # end
        # lik1=(log.(plimix));
        # lik1=sum(lik1,dims=3);
        # U=exp.(lik1)[:,:,1].*pg; 
        # cllik=sum(log.(sum(U,dims=2)))
    #     return cllik
    # end
    
    plimix=zeros(n,2,size(pairsindex,2));
    for idx=1:size(pairsindex,2) #parallelizzabile
        # CFA
        Scfa=Sigms1[pairsindex[:,idx],pairsindex[:,idx]]
        Rcfa=cov2cor(Scfa);
        @views Bcfa=hcat(tcfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tcfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
        @views tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);


        #EFA
        Sefa=Sigms2[pairsindex[:,idx],pairsindex[:,idx]]
        Refa=cov2cor(Sefa);
        @views Befa=hcat((tefa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tefa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
        @views tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);

    
        @views plimix[:,1,idx] = tcfa3 .* pg[:,1]
        @views plimix[:,2,idx] = tefa3 .* pg[:,2]
    end
    # % probabilities in each row must sum to 1 (application of Bayes law)
        
    if estep == true
        return plimix
    elseif estep == false
        # plimixll=zeros(Float64,n);
        # # plimix=zeros(n,2,size(pairsindex,2));
        # for idx=1:size(pairsindex,2) #parallelizzabile
        #     # CFA
        #     Scfa=Sigms1[pairsindex[:,idx],pairsindex[:,idx]]
        #     Rcfa=cov2cor(Scfa);
        #     sdScfa=sqrt.(diag(Scfa)')
        #     tcfaidx = tcfa[:,pairsindex[:,idx]]
        #     tcfa1idx=tcfa1[:,pairsindex[:,idx]]
        #     Bcfa=hcat(tcfaidx./sdScfa,tcfa1idx./sdScfa);
        #     tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            
            
        #     #EFA
        #     Sefa=Sigms2[pairsindex[:,idx],pairsindex[:,idx]]
        #     Refa=cov2cor(Sefa);
        #     sdSefa=sqrt.(diag(Sefa)')
        #     muij=mu[pairsindex[:,idx]]'
        #     tefaidx = tefa[:,pairsindex[:,idx]]
        #     tefa1idx = tefa1[:,pairsindex[:,idx]]
        #     Befa=hcat((tefaidx.-muij)./sdSefa,(tefa1idx.-muij)./sdSefa);
        #     tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            
        #     # plimix[:,1,idx] = tcfa3 .* pg[:,1]
        #     # plimix[:,2,idx] = tefa3 .* pg[:,2]
        #     plimixll += log.(tcfa3.*pg[:,1].+tefa3.*pg[:,2])
        # end
            
       
    # plimix: n × 2 × nPairs  (già include i pesi κ_i e 1-κ_i)
    S = sum(plimix, dims=2)       # n × 1 × nPairs  = κ_i fC + (1-κ_i) fE per ogni coppia
    S2 = @view S[:, 1, :]         # n × nPairs      (view, nessuna copia)
    S2 = clamp.(S2, 1e-300, Inf)  # stabilità numerica
    cllik = sum(log.(S2))         # somma su i e su tutte le coppie
    return cllik
    end
end
    


    

# function E_stepMixlm(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep)#::AbstractVector{T},q::Int64,K::Int64,p::Int64,n::Int64,c::Int64,X::Matrix{Int64},ics::Matrix,ncov::Int64,pairsindex::Matrix{Int64},L_str::Matrix,indices,estep::Bool) where T
#     # local trs::Int64,trsi::Int64,t1s::Vector{Real}, L1::Matrix{Real}, Phi::Matrix{Real},  t2s::Vector{Real}, L2::Matrix{Real}, Psi::Matrix{Real}, kappa::Vector{Real}, betas::Vector{Real}, propo::Vector{Real}, pg::Matrix{Real}, mu::Vector{Real}, Sigms1::Matrix{Real}, Sigms2::Matrix{Real}, br1::Matrix{Real}, br2::Matrix{Real}, t1::Matrix{Real}, t2::Matrix{Real}, tcfa1::Matrix{Real}, tcfa::Matrix{Real}, tefa1::Matrix{Real}, tefa::Matrix{Real}, plimix::Array{Real}, plimixll::Vector{Real}, Scfa::Matrix{Real}, Rcfa::Matrix{Real}, Bcfa::Matrix{Real}, tcfa3::Vector{Real}, Sefa::Matrix{Real}, Refa::Matrix{Real}, Befa::Matrix{Real}, tefa3::Vector{Real}, cllik::Real
#     if sum(param .== 0.) == size(param,1)
#         return -Inf
#     end
#     trs = (c-1-2)
#     trsi = (c-1)
#     t1s = param[1:trsi*p]
#     L1 = param[(trsi*p+1):(trsi*p+p)].*L_str 
#     Phi = reconstruct_correlation_matrix(param[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    
#     t2s = param[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs*p)]
#     L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi*p+p+(q*(q-1))/2+trs*p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))] )

    
#     Psi = diagm(param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p-1))]))#K
#     kappa = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+K)]   # interc2 = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p)]
#     betas = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+K+ncov)]  #  betas = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

#     propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
#     pg = hcat(propo,1. .-propo)

#     mu = L2*kappa
    
#     # Sigms=zeros(p,p,2);
#     # Deltas = Diagonal(Diagonal((L1*Phi*L1')) .+ 1);
#     # Sigms[:,:,1]=(Deltas^(-1/2))*(L1*Phi*L1' + diagm(ones(p)))*(Deltas^(-1/2))
#     Sigms1=L1*Phi*L1' + diagm(ones(p))
#     Sigms2=L2*L2'+Psi 
    
#     br1 = reshape(t1s,p,trsi)
#     t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
#     br2 = hcat(br1[:,1:2],reshape(t2s,p,trs))
#     t2=hcat(fill(-Inf,p),br2,fill(Inf,p))


#     tcfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     tcfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
#     tefa1=[t2[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     tefa=[t2[j,X[i,j]] for i in 1:n, j in 1:p]


#     # tlow = cat(tcfa,tefa,dims=3)
#     # tup = cat(tcfa1,tefa1,dims=3)
    
#     # if estep == true
#     #     low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#     #     plimix=ones(n,2,size(pairsindex,2));
    
#     #     for idx=1:size(pairsindex,2) #parallelizzabile
#     #         for g=1:2
#     #             if g==1
#     #                 s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#     #                 r2c=cov2cor(s2c);
#     #                 up=(tcfa1[:,pairsindex[:,idx]]./sqrt.(diag(s2c)'))
#     #                 low=(tcfa[:,pairsindex[:,idx]]./sqrt.(diag(s2c)'))
#     #             elseif g==2
#     #                 s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#     #                 r2c=cov2cor(s2c);
#     #                 up=(tefa1[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#     #                 low=(tefa[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#     #             end
#     #             B=hcat(low,up);
#     #             t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#     #             plimix[:,g,idx]= t3.*pg[:,g];
#     #         end
#     #     end
#     # # % probabilities in each row must sum to 1 (application of Bayes law)
#     #     return plimix
#     # elseif estep == false
#         # low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#         # plimix=ones(n,2,size(pairsindex,2));

#         # for idx=1:size(pairsindex,2) #parallelizzabile
#         #     for g=1:2
#         #         if g==1
#         #             s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#         #             r2c=cov2cor(s2c);
#         #             up=(tcfa1[:,pairsindex[:,idx]]./sqrt.(diag(s2c)'))
#         #             low=(tcfa[:,pairsindex[:,idx]]./sqrt.(diag(s2c)'))
#         #         elseif g==2
#         #             s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#         #             r2c=cov2cor(s2c);
#         #             up=(tefa1[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#         #             low=(tefa[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#         #         end
#         #         B=hcat(low,up);
#         #         t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#         #         plimix[:,g,idx] = t3
#         #     end
#         # end
#         # lik1=(log.(plimix));
#         # lik1=sum(lik1,dims=3);
#         # U=exp.(lik1)[:,:,1].*pg; 
#         # cllik=sum(log.(sum(U,dims=2)))
#     #     return cllik
#     # end
#     if estep == true
#         plimix=zeros(n,2,size(pairsindex,2));
#         for idx=1:size(pairsindex,2) #parallelizzabile
#             # CFA
#             Scfa=Sigms1[pairsindex[:,idx],pairsindex[:,idx]]
#             Rcfa=cov2cor(Scfa);
#             Bcfa=hcat(tcfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tcfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
#             tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            
            
#             #EFA
#             Sefa=Sigms2[pairsindex[:,idx],pairsindex[:,idx]]
#             Refa=cov2cor(Sefa);
#             Befa=hcat((tefa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tefa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
#             tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            
    
#             plimix[:,1,idx] = tcfa3 .* pg[:,1]
#             plimix[:,2,idx] = tefa3 .* pg[:,2]
#         end
#     # % probabilities in each row must sum to 1 (application of Bayes law)
#         return plimix
#     elseif estep == false
#         plimixll=zeros(Float64,n);
#         # plimix=zeros(n,2,size(pairsindex,2));
#         for idx=1:size(pairsindex,2) #parallelizzabile
#             # CFA
#             Scfa=Sigms1[pairsindex[:,idx],pairsindex[:,idx]]
#             Rcfa=cov2cor(Scfa);
#             sdScfa=sqrt.(diag(Scfa)')
#             tcfaidx = tcfa[:,pairsindex[:,idx]]
#             tcfa1idx=tcfa1[:,pairsindex[:,idx]]
#             Bcfa=hcat(tcfaidx./sdScfa,tcfa1idx./sdScfa);
#             tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            
            
#             #EFA
#             Sefa=Sigms2[pairsindex[:,idx],pairsindex[:,idx]]
#             Refa=cov2cor(Sefa);
#             sdSefa=sqrt.(diag(Sefa)')
#             muij=mu[pairsindex[:,idx]]'
#             tefaidx = tefa[:,pairsindex[:,idx]]
#             tefa1idx = tefa1[:,pairsindex[:,idx]]
#             Befa=hcat((tefaidx.-muij)./sdSefa,(tefa1idx.-muij)./sdSefa);
#             tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            
#             # plimix[:,1,idx] = tcfa3 .* pg[:,1]
#             # plimix[:,2,idx] = tefa3 .* pg[:,2]
#             plimixll += log.(tcfa3.*pg[:,1].+tefa3.*pg[:,2])
#         end
#         # cllik=sum(log.(sum(plimix,dims=2))) 
#         cllik=sum(plimixll) 
#         return cllik
#     end
# end
    





function obslik(param::AbstractVector{T},q::Int64,K::Int64,p::Int64,n::Int64,c::Int64,X::Matrix{Int64},ics::Matrix,ncov::Int64,pairsindex::Matrix{Int64},L_str::Matrix,indices::Vector,estep::Bool) where T
    local trs::Int64,trsi::Int64,t1s::Vector{Real}, L1::Matrix{Real}, Phi::Matrix{Real},  t2s::Vector{Real}, L2::Matrix{Real}, Psi::Matrix{Real}, kappa::Vector{Real}, betas::Vector{Real}, propo::Vector{Real}, pg::Matrix{Real}, mu::Vector{Real}, Sigms1::Matrix{Real}, Sigms2::Matrix{Real}, br1::Matrix{Real}, br2::Matrix{Real}, t1::Matrix{Real}, t2::Matrix{Real}, tcfa1::Matrix{Real}, tcfa::Matrix{Real}, tefa1::Matrix{Real}, tefa::Matrix{Real}, plimixll::Vector{Real}, Scfa::Matrix{Real}, Rcfa::Matrix{Real}, Bcfa::Matrix{Real}, tcfa3::Vector{Real}, Sefa::Matrix{Real}, Refa::Matrix{Real}, Befa::Matrix{Real}, tefa3::Vector{Real}, cllik::Real
    # if sum(param .== 0.) == size(param,1)
    #     return -Inf
    # end
    trs = (c-1-2)*p
    trsi = (c-1)*p

    t1s = param[1:trsi]
    L1 = param[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
    
    t2s = param[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+trs)]
    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi+p+(q*(q-1))/2+trs+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))] )

    
    Psi = diagm(param[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,param[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    kappa = param[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+K)]   # interc2 = param[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p)]
    betas = param[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+K+ncov)]  #  betas = param[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)

    mu = L2*kappa
    
    # Sigms=zeros(p,p,2);
    # Deltas = Diagonal(Diagonal((L1*Phi*L1')) .+ 1);
    # Sigms[:,:,1]=(Deltas^(-1/2))*(L1*Phi*L1' + diagm(ones(p)))*(Deltas^(-1/2))
    Sigms1=L1*Phi*L1' + diagm(ones(p))
    Sigms2=L2*L2'+Psi 
    
    br1 = reshape(t1s,p,c-1)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
    br2 = hcat(br1[:,1:2],reshape(t2s,p,c-1-2))
    t2=hcat(fill(-Inf,p),br2,fill(Inf,p))


    tcfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tcfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
    tefa1=[t2[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tefa=[t2[j,X[i,j]] for i in 1:n, j in 1:p]

        plimixll=zeros(Float64,n);
        # plimix=zeros(n,2,size(pairsindex,2));
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sigms1[pairsindex[:,idx],pairsindex[:,idx]]
            Rcfa=cov2cor(Scfa);
            sdScfa=sqrt.(diag(Scfa)')
            tcfaidx = tcfa[:,pairsindex[:,idx]]
            tcfa1idx=tcfa1[:,pairsindex[:,idx]]
            Bcfa=hcat(tcfaidx./sdScfa,tcfa1idx./sdScfa);
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            
            
            #EFA
            Sefa=Sigms2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            sdSefa=sqrt.(diag(Sefa)')
            muij=mu[pairsindex[:,idx]]'
            tefaidx = tefa[:,pairsindex[:,idx]]
            tefa1idx = tefa1[:,pairsindex[:,idx]]
            Befa=hcat((tefaidx.-muij)./sdSefa,(tefa1idx.-muij)./sdSefa);
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            
            # plimix[:,1,idx] = tcfa3 .* pg[:,1]
            # plimix[:,2,idx] = tefa3 .* pg[:,2]
            plimixll += log.(tcfa3.*pg[:,1].+tefa3.*pg[:,2])
        end
        # cllik=sum(log.(sum(plimix,dims=2))) 
        cllik=sum(plimixll) 
        return cllik
end

    


function M_parsMixlm(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)

    trsi = (c-1)*p #attenzione p è già qui
    trs = (c-1-2)*p

    # PML - CFA & EFA
    theta_init = theta[1:(size(theta,1)-ncov)]#0.000001
    # # -1 loadings
    # lb=vec(vcat(fill(-5.1,trsi),fill(0.00001,p),fill(-1.,Int((q*(q-1))/2)),fill(1.,trs),fill(-1.,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-10,p)));
    # # 0. loadings
    lb=vec(vcat(fill(-Inf,trsi),fill(-Inf,p),fill(-1.,Int((q*(q-1))/2)),fill(-Inf,trs),fill(-Inf,Int(p*K-K*(K-1)/2)),fill(1.e-6,p),fill(-Inf,K)));#1.0e-5
    ub=vec(vcat(fill(Inf,trsi),fill(Inf,p),fill(1.,Int((q*(q-1))/2)),fill(Inf,trs),fill(Inf,Int(p*K-K*(K-1)/2)),fill(Inf,p),fill(Inf,K)));

    par = zeros(Float64,size(lb))
    pair_weights = PML_FAST[] ? pml_pair_category_weights(X,T,pairsindex,c) : nothing
    if it == 1
        par = optimize(Optim.only_fg!((F,G,x)->fglm!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices,pair_weights)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min ,Optim.Options(iterations=1000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0))
    elseif it > 1
        par = optimize(Optim.only_fg!((F,G,x)->fglm!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices,pair_weights)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min
    end 

    # func = TwiceDifferentiable((x)->pmllastlm(x,X,T,p,q,K,c,n,pairsindex,L_str,indices),(G,x)->fun_grad!(G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices),theta_init; autodiff=:forward);
    # dfc = TwiceDifferentiableConstraints(lb, ub)
    # par = optimize(func,dfc,theta_init,IPNewton(),Optim.Options(iterations=10000,f_reltol=1.0000e-06))


    # par = optimize(x->pmllast(x,X,T,p,q,K,c,n,pairsindex,L_str,indices),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

    # par = zeros(Float64,size(lb))
    # if it == 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 300.0)) #20 min
    # elseif it > 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 1200.0)) #20 min
    # end 
    
    # mixture parameter
    betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+K+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
    Ez1 = T[:,1,:]
    Ez0 = T[:,2,:]
    
    pairs = size(pairsindex,2)
    betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)

    
    # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[1,1];
    # betasnew = vcat(log(pg/(1 - pg)))
   
    # println(par.minimizer,pg[1])
    if !(Optim.converged(par)) && isnan(betas[1]) 
        println(betas[1])
        if start == true
            return zeros(Float64,size(theta,1))
        elseif start == false
            return Optim.converged(par)
        end
    end
    return vcat(par.minimizer,betasnew)#,pg[1])
end



function pmllastlm(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    trs = (c-1-2)
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])


    t2s = x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs*p)]
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+trs*p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))] )

    Psi = diagm(x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    kappa = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+K)]    #interc2 = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p)]
    # prop = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]

    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
    br2 = hcat(br1[:,1:2],reshape(t2s,p,trs))
    t2=hcat(fill(-Inf,p),br2,fill(Inf,p))

    mu = L2*kappa
    
    # Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    
    # Deltas = Diagonal(Diagonal((L1*Phi*L1')) .+ 1);
    # Sig1=(Deltas^(-1/2))*(L1*Phi*L1' + diagm(ones(p)))*(Deltas^(-1/2))
    
    Sig1=L1*Phi*L1' + diagm(ones(p))
    Sig2=L2*L2'+Psi



    plimix=zeros(Float64,n);
    tcfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tcfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
    tefa1=[t2[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tefa=[t2[j,X[i,j]] for i in 1:n, j in 1:p]
    
    for idx=1:size(pairsindex,2) #parallelizzabile
        # CFA
        Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
        Rcfa=cov2cor(Scfa);
        Bcfa=hcat(tcfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tcfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
        tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
        Tt1 = T[:,1,idx]
        
        #EFA
        Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
        Refa=cov2cor(Sefa);
        Befa=hcat((tefa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tefa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
        tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
        Tt2 = T[:,2,idx]

        plimix += log.(tcfa3).*Tt1.+log.(tefa3).*Tt2
    end
    f=-sum(plimix)
return f
end

function pml_pair_category_weights(X,T,pairsindex,c)
    n = size(X,1)
    npairs = size(pairsindex,2)
    W = zeros(Float64,c*c,2,npairs)
    @inbounds for pair in 1:npairs
        i = pairsindex[1,pair]
        j = pairsindex[2,pair]
        for r in 1:n
            category = X[r,i] + (X[r,j]-1)*c
            W[category,1,pair] += T[r,1,pair]
            W[category,2,pair] += T[r,2,pair]
        end
    end
    return W
end

function fglm!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices, W=nothing)
    # do common computations here
    # ...
    trs = (c-1-2)
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
   

    t2s = x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs*p)]
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+trs*p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))] )

    Psi = diagm(x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    kappa = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+K)]    #interc2 = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p)]
    # prop = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]

    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
    br2 = hcat(br1[:,1:2],reshape(t2s,p,trs))
    t2=hcat(fill(-Inf,p),br2,fill(Inf,p))

    mu = L2*kappa

    # Deltas = Diagonal(Diagonal((L1*Phi*L1')) .+ 1);
    

    Sig1=(L1*Phi*L1' + diagm(ones(p))) #(Deltas^(-1/2))*(L1*Phi*L1' + diagm(ones(p)))*(Deltas^(-1/2))
    Sig2=L2*L2'+Psi


    if G !== nothing
      
        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradkappacefa = zeros(Float64,K,size(pairsindex,2))
        gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)

        gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)


        gradlscfaidx = zeros(Float64,q,size(pairsindex,2),2)
        gradlambci = zeros(Float64,q,c+1,size(pairsindex,2),2)
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
            
        idx = 0
        nwork = W === nothing ? n : c*c
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1

                sig1s = hcat(vcat(Sig1[i,i],Sig1[i,j]),vcat(Sig1[j,i],Sig1[j,j]))
                cormat1 = cov2cor(sig1s)
                rho1 = cormat1[1,2]
                sig1sij = diag(sig1s)   
                sig1ij = sqrt.(sig1sij)
                sig1ij13 = sig1ij[1]^(3)
                sig1ij23 = sig1ij[2]^(3)
                L1i=L1[i,:]
                L1j=L1[j,:]

                sig2s = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sig2s)
                rho2 = cormat2[1,2]
                muij = mu[vcat(i,j)]
                sig2sij = diag(sig2s)
                sig2ij = sqrt.(sig2sij)
                sig2ij13 = sig2ij[1]^(3)
                sig2ij23 = sig2ij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]
                Psii=Psi[i,i]
                Psij=Psi[j,j]

                @inbounds for N in 1:nwork
                    if W === nothing
                        xni = X[N,i]
                        xnj = X[N,j]
                    else
                        xni = (N-1) % c + 1
                        xnj = div(N-1,c) + 1
                        if W[N,1,idx] == 0.0 && W[N,2,idx] == 0.0
                            continue
                        end
                    end
                    t1i = t1[i,xni]
                    t1i1 = t1[i,xni+1]
                    t1j = t1[j,xnj]
                    t1j1 = t1[j,xnj+1]
                    t1ist = (t1i)/sig1ij[1]
                    t1i1st = (t1i1)/sig1ij[1]
                    t1jst = (t1j)/sig1ij[2]
                    t1j1st = (t1j1)/sig1ij[2]
                    T1 = W === nothing ? T[N,1,idx][1] : W[N,1,idx]
                    T2 = W === nothing ? T[N,2,idx][1] : W[N,2,idx]
                    t2i = t2[i,xni]
                    t2i1 = t2[i,xni+1]
                    t2j = t2[j,xnj]
                    t2j1 = t2[j,xnj+1]
                    t2ist = (t2i-muij[1])/sig2ij[1]
                    t2i1st = (t2i1-muij[1])/sig2ij[1]
                    t2jst = (t2j-muij[2])/sig2ij[2]
                    t2j1st = (t2j1-muij[2])/sig2ij[2]
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                
                    pijicfa = bvnP(t1ist, t1i1st, t1jst, t1j1st, rho1, false)
                    derivNorm1i1 = (T1/pijicfa)*(pdf(Normal(0,1),t1i1st)*(cdf(Normal(0,1),(t1j1st-rho1*t1i1st)/sqrtrij1)-cdf(Normal(0,1),(t1jst-rho1*t1i1st)/sqrtrij1)))
                    derivNorm1i = (T1/pijicfa)*(pdf(Normal(0,1),t1ist)*(cdf(Normal(0,1),(t1jst-rho1*t1ist)/sqrtrij1)-cdf(Normal(0,1),(t1j1st-rho1*t1ist)/sqrtrij1)))
                    derivNorm1j1 = (T1/pijicfa)*(pdf(Normal(0,1),t1j1st)*(cdf(Normal(0,1),(t1i1st-rho1*t1j1st)/sqrtrij1)-cdf(Normal(0,1),(t1ist-rho1*t1j1st)/sqrtrij1)))
                    derivNorm1j = (T1/pijicfa)*(pdf(Normal(0,1),t1jst)*(cdf(Normal(0,1),(t1ist-rho1*t1jst)/sqrtrij1)-cdf(Normal(0,1),(t1i1st-rho1*t1jst)/sqrtrij1)))
                
                    gradtaucfa[xni+1,idx,1] +=  derivNorm1i1/(sig1ij[1])
                    gradtaucfa[xni,idx,1] += derivNorm1i/(sig1ij[1])
                    gradtaucfa[xnj+1,idx,2] += derivNorm1j1/(sig1ij[2])
                    gradtaucfa[xnj,idx,2] += derivNorm1j/(sig1ij[2])


                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))

                    if isnan(derivNormi1)
                        derivNormi1 = 0.
                    end
                    if isnan(derivNormi)
                        derivNormi = 0.
                    end
                    if isnan(derivNormj1)
                        derivNormj1 = 0.
                    end
                    if isnan(derivNormj)
                        derivNormj = 0.
                    end

                    
                    if !(t1ist == -Inf && t1jst == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1jst))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else                  
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1st == Inf && t1j1st == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1j1st))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end

                    if !(t1i1st == Inf && t1jst == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1jst))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end

                    if !(t1ist == -Inf && t1j1st == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1j1st))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    dpidrij1 = (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)#*vcat((-((L1j'*Phi.*t1j1st)./(2*sig1sij[2]))),(-((L1i'*Phi.*t1i1st)./(2*sig1sij[1]))))#vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)#
                

                    gradlscfaidx[1:q,idx,1] .+= dpidrij1*(Phi*L1j)/(sig1ij13*sig1ij[2])#(-Phi*L1i*(L1i'*L1j)/(sig1ij13*sig1ij[2])+Phi*L1j/((sig1ij[1]*sig1ij[2])))
                    gradlscfaidx[1:q,idx,2] .+= dpidrij1*(Phi*L1i)/(sig1ij23*sig1ij[1])#(-Phi*L1j*(L1j'*L1i)/(sig1ij23*sig1ij[1])+Phi*L1i/((sig1ij[2]*sig1ij[1])))
                    
                    gradlambci[1:q,xni+1,idx,1] .+=  derivNorm1i1 *(-((L1i*t1i1st)./(sig1sij[1])))#(-((L1i*t1i1).+(kappak*Psii))./((sig1ij13)))
                    gradlambci[1:q,xnj+1,idx,2] .+=  derivNorm1j1 *(-((L1j*t1j1st)./(sig1sij[2])))#(-((L1j*t1j1).+(kappak*Psij))./((sig1ij23)))
                    gradlambci[1:q,xni,idx,1]   .+=  derivNorm1i  *(-((L1i*t1ist)./(sig1sij[1])))#(-((L1i*t1i).+(kappak*Psii))./((sig1ij13)))
                    gradlambci[1:q,xnj,idx,2]   .+=  derivNorm1j  *(-((L1j*t1jst)./(sig1sij[2])))#(-((L1j*t1j).+(kappak*Psij))./((sig1ij23)))               
                    
                    gradphici[1:q,1:q,idx] .+= dpidrij1*(L1i*L1j')/(sig1ij[1]*sig1ij[2])

                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sig2ij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sig2ij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sig2ij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sig2ij[2])


                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)

                    # gradkappacefa[xni+1,1:K,idx,1] += (derivNormi1*(-L2[i,:]/sig2ij[1]))
                    # gradkappacefa[xnj+1,1:K,idx,2] += (derivNormj1*(-L2[j,:]/sig2ij[2]))
                    # gradkappacefa[xni,1:K,idx,1] += (derivNormi*(-L2[i,:]/sig2ij[1]))
                    # gradkappacefa[xnj,1:K,idx,2] += (derivNormj*(-L2[j,:]/sig2ij[2]))
                    gradkappacefa[1:K,idx] .+= derivNormi1*(-L2[i,:]/sig2ij[1]) + derivNormi*(-L2[i,:]/sig2ij[1]) + derivNormj1*(-L2[j,:]/sig2ij[2]) + derivNormj*(-L2[j,:]/sig2ij[2])
                    
                    gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sig2ij13)*(sig2ij[2])))
                    gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sig2ij23)*(sig2ij[1])))

                    gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij[1]-t2i1))/(2*(sig2ij13))
                    gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij[2]-t2j1))/(2*(sig2ij23))
                    gradpsiefa[xni,idx,1] += (derivNormi*(muij[1]-t2i))/(2*(sig2ij13))
                    gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij[2]-t2j))/(2*(sig2ij23))
                    
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        kappak=kappa[k]
                        gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sig2ij13*sig2ij[2])+L2jk/((sig2ij[1]*sig2ij[2])))#((L2j*Psi[i,i])/(sig2ij13*sig2ij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sig2ij23*sig2ij[1])+L2ik/((sig2ij[2]*sig2ij[1])))#((L2[i,:].*Psi[j,j])./(sig2ij23*sig2ij[1]))

                        gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(-((L2ik*t2i1st)./(sig2sij[1]).+(kappak)./(sig2ij[1])))#(-((L2ik*t2i1st).+(kappak*Psii))./((sig2ij13)))
                        gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(-((L2jk*t2j1st)./(sig2sij[2]).+(kappak)./(sig2ij[2])))#(-((L2jk*t2j1st).+(kappak*Psij))./((sig2ij23)))
                        gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(-((L2ik*t2ist)./(sig2sij[1]).+(kappak)./(sig2ij[1])))#(-((L2ik*t2ist).+(kappak*Psii))./((sig2ij13)))
                        gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(-((L2jk*t2jst)./(sig2sij[2]).+(kappak)./(sig2ij[2])))#(-((L2jk*t2jst).+(kappak*Psij))./((sig2ij23)))
                    end 
                
                end
            end
        end

        vecpair = vec(pairsindex)
        liktaushared = zeros(p*(2));idpcfa = 0   
        @inbounds for j in 2:3
            @inbounds for i in 1:p
                idpcfa += 1
                liktaushared[idpcfa] = -sum(vec(vcat(gradtaucfa[j,:,1]',gradtaucfa[j,:,2]'))[vecpair.==i])-sum(vec(vcat(gradtauefa[j,:,1]',gradtauefa[j,:,2]'))[vecpair.==i])
            end
        end

        liktau = zeros(p*(c-1-2));idpcfa = 0   
        @inbounds for j in 4:c
            @inbounds for i in 1:p
                idpcfa += 1
                liktau[idpcfa] = -sum(vec(vcat(gradtaucfa[j,:,1]',gradtaucfa[j,:,2]'))[vecpair.==i])
            end
        end

        liktau2 = zeros(p*(c-1-2));idpefa = 0   
        @inbounds for j in (4):c
            @inbounds for i in 1:p
                idpefa += 1
                liktau2[idpefa] = -sum(vec(vcat(gradtauefa[j,:,1]',gradtauefa[j,:,2]'))[vecpair.==i])
            end
        end



        likls1 = zeros(p,q); llamb1 = sum(gradlambci[:,2:c,:,:],dims=2)[:,1,:,:]; llam1 = zeros(p,q); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
        likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
        llam1idx  = zeros(p,q); likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        @inbounds for i in 1:p
            @inbounds for qs in 1:q
                llam1[i,qs] = -sum(vec(llamb1[qs,:,:]')[vecpair.==i])
                llam1idx[i,qs] = -sum(vec(gradlscfaidx[qs,:,:]')[vecpair.==i])
            end
            likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
            likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end

        likpsitrue = likpsi.+likpsiidx
        liklstrue = likls.+liklsidx
        liklamtrue = llam1.+llam1idx

        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]

        likappa = -sum(gradkappacefa,dims=2)

        G .= vec(vcat(liktaushared,liktau,sum(liklamtrue.*L_str,dims=2),likphida,liktau2,liklstrue[indices],likpsitrue,likappa))#likpsitrue[2:end]

    end
    if F !== nothing 
  
        if W !== nothing
            fast_value = 0.0
            @inbounds for idx in 1:size(pairsindex,2)
                i = pairsindex[1,idx]
                j = pairsindex[2,idx]
                Scfa = Sig1[pairsindex[:,idx],pairsindex[:,idx]]
                Sefa = Sig2[pairsindex[:,idx],pairsindex[:,idx]]
                rho1 = cov2cor(Scfa)[1,2]
                rho2 = cov2cor(Sefa)[1,2]
                sd11 = sqrt(Scfa[1,1])
                sd12 = sqrt(Scfa[2,2])
                sd21 = sqrt(Sefa[1,1])
                sd22 = sqrt(Sefa[2,2])
                for category in 1:(c*c)
                    weight1 = W[category,1,idx]
                    weight2 = W[category,2,idx]
                    if weight1 == 0.0 && weight2 == 0.0
                        continue
                    end
                    a = (category-1) % c + 1
                    b = div(category-1,c) + 1
                    if weight1 != 0.0
                        prob1 = bvnP(t1[i,a]/sd11,t1[i,a+1]/sd11,
                                     t1[j,b]/sd12,t1[j,b+1]/sd12,rho1,false)
                        fast_value -= weight1*log(prob1)
                    end
                    if weight2 != 0.0
                        prob2 = bvnP((t2[i,a]-mu[i])/sd21,(t2[i,a+1]-mu[i])/sd21,
                                     (t2[j,b]-mu[j])/sd22,(t2[j,b+1]-mu[j])/sd22,rho2,false)
                        fast_value -= weight2*log(prob2)
                    end
                end
            end
            return fast_value
        end
    
        plimix=zeros(Float64,n);
        tcfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
        tcfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
        tefa1=[t2[j,X[i,j]+1] for i in 1:n, j in 1:p]
        tefa=[t2[j,X[i,j]] for i in 1:n, j in 1:p]
              
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
            Rcfa=cov2cor(Scfa);
            Bcfa=hcat(tcfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tcfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            Tt1 = T[:,1,idx]
            
            #EFA
            Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tefa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tefa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            Tt2 = T[:,2,idx]

            plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
        end
        f=-sum(plimix)
        # @info "Current objective: $f"
        return f
    end
end


using LinearAlgebra, Distributions

# densità bivariata standard con correlazione rho (senza allocazioni)
@inline function φ2_std(u::Float64, v::Float64, ρ::Float64)
    ρ2 = ρ*ρ
    inv_s = inv(sqrt(1.0 - ρ2))
    z = (u*u - 2.0*ρ*u*v + v*v) * (0.5 * inv_s^2)
    return (0.15915494309189535) * inv_s * exp(-z)  # 1/(2π)=0.159154...
end

# function fglm2!(F, G, x, X, T, p, q, K, c, n, pairsindex, L_str, indices)
#     trs  = c - 1 - 2
#     trsi = c - 1

#     ---- slice parametri (solo viste, no copie) ----
#     t1s   = @view x[1:trsi*p]
#     λ1raw = @view x[trsi*p+1 : trsi*p+p]
#     ϕvec  = @view x[trsi*p+p+1 : trsi*p+p + (q*(q-1))÷2]
#     t2s   = @view x[trsi*p+p + (q*(q-1))÷2 + 1 : trsi*p+p + (q*(q-1))÷2 + trs*p]
#     l2vec = @view x[trsi*p+p + (q*(q-1))÷2 + trs*p + 1:trsi*p+p + (q*(q-1))÷2 + trs*p + (p*K - (K*(K-1))÷2)]
#     ψvec  = @view x[trsi*p+p + (q*(q-1))÷2 + trs*p + (p*K - (K*(K-1))÷2) + 1 : trsi*p+p + (q*(q-1))÷2 + trs*p + (p*K - (K*(K-1))÷2) + p]
#     κ     = @view x[trsi*p+p + (q*(q-1))÷2 + trs*p + (p*K - (K*(K-1))÷2) + p + 1 :  trsi*p+p + (q*(q-1))÷2 + trs*p + (p*K - (K*(K-1))÷2) + p + K]

#     ---- ricostruzioni parametriche ----
#     t1 = hcat(fill(-Inf, p), reshape(t1s, p, trsi), fill(Inf, p))
#     t2 = hcat(fill(-Inf, p), hcat(reshape(t1s, p, trsi)[:,1:2], reshape(t2s, p, trs)), fill(Inf, p))

#     L1 = λ1raw .* L_str                         # p×q mascherato
#     Φ  = reconstruct_correlation_matrix(ϕvec)   # q×q
#     L2 = reconstruct_L2_matrix(p, K, indices, l2vec)  # p×K
#     μ  = L2 * κ                                 # p

#     Sig1 = L1*Φ*L1' + I   (non materializziamo Sig2)
#     Sig1 = L1 * Φ * L1' .+ I

#     n01 = Normal()  # riutilizzato

#     =================== GRADIENTE ===================
#     if G !== nothing
#         prealloc grad storage
#         gradtaucfa   = zeros(Float64, c+1, size(pairsindex,2), 2)
#         gradtauefa   = zeros(Float64, c+1, size(pairsindex,2), 2)
#         gradkappacefa = zeros(Float64, K, size(pairsindex,2))
#         gradpsiefa   = zeros(Float64, c+1, size(pairsindex,2), 2)

#         gradpsiefaidx = zeros(Float64, size(pairsindex,2), 2)
#         gradlsefaidx  = zeros(Float64, K, size(pairsindex,2), 2)

#         gradlscfaidx  = zeros(Float64, q, size(pairsindex,2), 2)
#         gradlambci    = zeros(Float64, q, c+1, size(pairsindex,2), 2)
#         gradphici     = zeros(Float64, q, q, size(pairsindex,2))
#         gradlamb2ci   = zeros(Float64, K, c+1, size(pairsindex,2), 2)

#         idx = 0
#         @inbounds for i in 1:(p-1)
#             L1i = @view L1[i, :]
#             L2i = @view L2[i, :]
#             ψi  = ψvec[i]
#             s11_1 = dot(L1i, Φ * L1i) + 1.0
#             s11_2 = dot(L2i, L2i) + ψi
#             sd11  = sqrt(s11_1); sd21 = sqrt(s11_2)
#             sd11³ = sd11*sd11*sd11; sd21³ = sd21*sd21*sd21

#             @inbounds for j in (i+1):p
#                 idx += 1
#                 L1j = @view L1[j, :]
#                 L2j = @view L2[j, :]
#                 ψj  = ψvec[j]

#                 CFA: s12, s22, rho1
#                 s12_1 = dot(L1i, Φ * L1j)
#                 s22_1 = dot(L1j, Φ * L1j) + 1.0
#                 sd12  = sqrt(s22_1)
#                 ρ1    = s12_1 / (sd11 * sd12)
#                 one_m_ρ1² = 1.0 - ρ1*ρ1
#                 inv_sqrt_1mρ1² = inv(sqrt(one_m_ρ1²))

#                 EFA: s12, s22, rho2 (senza materializzare Sig2)
#                 s12_2 = dot(L2i, L2j)
#                 s22_2 = dot(L2j, L2j) + ψj
#                 sd22  = sqrt(s22_2)
#                 ρ2    = s12_2 / (sd21 * sd22)
#                 one_m_ρ2² = 1.0 - ρ2*ρ2
#                 inv_sqrt_1mρ2² = inv(sqrt(one_m_ρ2²))

#                 μi = μ[i]; μj = μ[j]

#                 T1v = @view T[:, 1, idx]  # vettori n
#                 T2v = @view T[:, 2, idx]

#                 @inbounds for N in 1:n
#                     xni = X[N,i];  xnj = X[N,j]
#                     thresholds
#                     t1i   = t1[i, xni]
#                     t1i1  = t1[i, xni+1]
#                     t1jv  = t1[j, xnj]
#                     t1j1  = t1[j, xnj+1]

#                     t2i   = t2[i, xni]
#                     t2i1  = t2[i, xni+1]
#                     t2jv  = t2[j, xnj]
#                     t2j1  = t2[j, xnj+1]

#                     standardizzazioni
#                     t1ist  = t1i  / sd11
#                     t1i1st = t1i1 / sd11
#                     t1jst  = t1jv / sd12
#                     t1j1st = t1j1 / sd12

#                     t2ist  = (t2i  - μi) / sd21
#                     t2i1st = (t2i1 - μi) / sd21
#                     t2jst  = (t2jv - μj) / sd22
#                     t2j1st = (t2j1 - μj) / sd22

#                     sq1 = sqrt(one_m_ρ1²)
#                     sq2 = sqrt(one_m_ρ2²)

#                     Prob quadratino bivar (CFA/EFA)
#                     pijicfa = bvnP(t1ist, t1i1st, t1jst, t1j1st, ρ1, false)
#                     pijiefa = bvnP(t2ist, t2i1st, t2jst, t2j1st, ρ2, false)

#                     derivate per soglie (CFA)
#                     T1 = T1v[N]; T2 = T2v[N]

#                     cdf_j1_i1  = cdf(n01, (t1j1st - ρ1*t1i1st)/sq1)
#                     cdf_j_i1   = cdf(n01, (t1jst  - ρ1*t1i1st)/sq1)
#                     cdf_i1_j1  = cdf(n01, (t1i1st - ρ1*t1j1st)/sq1)
#                     cdf_i_j1   = cdf(n01, (t1ist  - ρ1*t1j1st)/sq1)
#                     cdf_j_i    = cdf(n01, (t1jst  - ρ1*t1ist)/sq1)
#                     cdf_j1_i   = cdf(n01, (t1j1st - ρ1*t1ist)/sq1)

#                     ϕ_i1 = pdf(n01, t1i1st);  ϕ_i  = pdf(n01, t1ist)
#                     ϕ_j1 = pdf(n01, t1j1st);  ϕ_j  = pdf(n01, t1jst)

#                     derivNorm1i1 = (T1/pijicfa) * (ϕ_i1*(cdf_j1_i1 - cdf_j_i1))
#                     derivNorm1i  = (T1/pijicfa) * (ϕ_i *(cdf_j_i  - cdf_j1_i))
#                     derivNorm1j1 = (T1/pijicfa) * (ϕ_j1*(cdf_i1_j1 - cdf_i_j1))
#                     derivNorm1j  = (T1/pijicfa) * (ϕ_j *(cdf_i_j  - cdf_i1_j))

#                     aggiorna gradienti soglie CFA
#                     gradtaucfa[xni+1,idx,1] +=  derivNorm1i1 / sd11
#                     gradtaucfa[xni,  idx,1] +=  derivNorm1i  / sd11
#                     gradtaucfa[xnj+1,idx,2] +=  derivNorm1j1 / sd12
#                     gradtaucfa[xnj,  idx,2] +=  derivNorm1j  / sd12

#                     derivate per soglie (EFA)
#                     ϕ2_i1 = pdf(n01, t2i1st);  ϕ2_i = pdf(n01, t2ist)
#                     ϕ2_j1 = pdf(n01, t2j1st);  ϕ2_j = pdf(n01, t2jst)

#                     cdf2_j1_i1 = cdf(n01, (t2j1st - ρ2*t2i1st)/sq2)
#                     cdf2_j_i1  = cdf(n01, (t2jst  - ρ2*t2i1st)/sq2)
#                     cdf2_i1_j1 = cdf(n01, (t2i1st - ρ2*t2j1st)/sq2)
#                     cdf2_i_j1  = cdf(n01, (t2ist  - ρ2*t2j1st)/sq2)
#                     cdf2_j_i   = cdf(n01, (t2jst  - ρ2*t2ist)/sq2)
#                     cdf2_j1_i  = cdf(n01, (t2j1st - ρ2*t2ist)/sq2)

#                     derivNormi1 = (T2/pijiefa) * (ϕ2_i1*(cdf2_j1_i1 - cdf2_j_i1))
#                     derivNormi  = (T2/pijiefa) * (ϕ2_i *(cdf2_j_i  - cdf2_j1_i))
#                     derivNormj1 = (T2/pijiefa) * (ϕ2_j1*(cdf2_i1_j1 - cdf2_i_j1))
#                     derivNormj  = (T2/pijiefa) * (ϕ2_j *(cdf2_i_j  - cdf2_i1_j))

#                     aggiornamenti EFA (thresholds, κ, ψ, L2)
#                     gradtauefa[xni+1,idx,1] +=  derivNormi1 / sd21
#                     gradtauefa[xni,  idx,1] +=  derivNormi  / sd21
#                     gradtauefa[xnj+1,idx,2] +=  derivNormj1 / sd22
#                     gradtauefa[xnj,  idx,2] +=  derivNormj  / sd22

#                     termini "bordo" per rho-derivatives (senza MvNormal)
#                     xnixnjc   = (isfinite(t1ist) && isfinite(t1jst))   ? φ2_std(t1ist,  t1jst,  ρ1) : 0.0
#                     xni1xnj1c = (isfinite(t1i1st) && isfinite(t1j1st)) ? φ2_std(t1i1st, t1j1st, ρ1) : 0.0
#                     xni1xnjc  = (isfinite(t1i1st) && isfinite(t1jst))  ? φ2_std(t1i1st, t1jst,  ρ1) : 0.0
#                     xnixnj1c  = (isfinite(t1ist)  && isfinite(t1j1st)) ? φ2_std(t1ist,  t1j1st, ρ1) : 0.0

#                     dpidrij1 = (T1/pijicfa) * (xni1xnj1c - xni1xnjc - xnixnj1c + xnixnjc)

#                     gradlscfaidx[:,idx,1] .+= dpidrij1 * (Φ * L1j) / (sd11³ * sd12)
#                     gradlscfaidx[:,idx,2] .+= dpidrij1 * (Φ * L1i) / (sd12*sd12*sd12 * sd11)

#                     gradlambci[:,xni+1,idx,1] .+=  derivNorm1i1 * (-(L1i*t1i1st)/s11_1)
#                     gradlambci[:,xnj+1,idx,2] .+=  derivNorm1j1 * (-(L1j*t1j1st)/s22_1)
#                     gradlambci[:,xni,  idx,1] .+=  derivNorm1i  * (-(L1i*t1ist )/s11_1)
#                     gradlambci[:,xnj,  idx,2] .+=  derivNorm1j  * (-(L1j*t1jst )/s22_1)

#                     gradphici[:,:,idx]        .+=  dpidrij1 * (L1i * L1j') / (sd11*sd12)

#                     parte EFA: analoghi "bordo" con rho2
#                     xnixnje   = (isfinite(t2ist) && isfinite(t2jst))   ? φ2_std(t2ist,  t2jst,  ρ2) : 0.0
#                     xni1xnj1e = (isfinite(t2i1st) && isfinite(t2j1st)) ? φ2_std(t2i1st, t2j1st, ρ2) : 0.0
#                     xni1xnje  = (isfinite(t2i1st) && isfinite(t2jst))  ? φ2_std(t2i1st, t2jst,  ρ2) : 0.0
#                     xnixnj1e  = (isfinite(t2ist)  && isfinite(t2j1st)) ? φ2_std(t2ist,  t2j1st, ρ2) : 0.0

#                     dpidrij = (T2/pijiefa) * (xni1xnj1e - xni1xnje - xnixnj1e + xnixnje)

#                     κ, ψ, L2
#                     gradkappacefa[:,idx] .+= derivNormi1 * (.-L2i ./ sd21) .+ derivNormi * (.-L2i ./ sd21) .+
#                                              derivNormj1 * (.-L2j ./ sd22) .+ derivNormj * (.-L2j ./ sd22)

#                     gradpsiefaidx[idx,1] += dpidrij * (-(dot(L2i,L2j)) / (2 * sd21³ * sd22))
#                     gradpsiefaidx[idx,2] += dpidrij * (-(dot(L2j,L2i)) / (2 * sd22*sd22*sd22 * sd21))

#                     gradpsiefa[xni+1,idx,1] += (derivNormi1 * (μi - t2i1)) / (2 * sd21³)
#                     gradpsiefa[xnj+1,idx,2] += (derivNormj1 * (μj - t2j1)) / (2 * sd22*sd22*sd22)
#                     gradpsiefa[xni,  idx,1] += (derivNormi  * (μi - t2i )) / (2 * sd21³)
#                     gradpsiefa[xnj,  idx,2] += (derivNormj  * (μj - t2j )) / (2 * sd22*sd22*sd22)

#                     @inbounds for k in 1:K
#                         L2ik = L2i[k];  L2jk = L2j[k]
#                         termini L2 (semplificati; stessi segni dei tuoi)
#                         gradlsefaidx[k,idx,1] += dpidrij * ( -L2ik*(dot(L2i,L2j))/(sd21³*sd22) + L2jk/(sd21*sd22) )
#                         gradlsefaidx[k,idx,2] += dpidrij * ( -L2jk*(dot(L2j,L2i))/(sd22*sd22*sd22*sd21) + L2ik/(sd22*sd21) )

#                         gradlamb2ci[k,xni+1,idx,1] += derivNormi1 * ( -((L2ik*t2i1st)/s11_1 + 1.0/sd21) )
#                         gradlamb2ci[k,xnj+1,idx,2] += derivNormj1 * ( -((L2jk*t2j1st)/s22_1 + 1.0/sd22) )
#                         gradlamb2ci[k,xni,  idx,1] += derivNormi  * ( -((L2ik*t2ist )/s11_1 + 1.0/sd21) )
#                         gradlamb2ci[k,xnj,  idx,2] += derivNormj  * ( -((L2jk*t2jst )/s22_1 + 1.0/sd22) )
#                     end
#                 end
#             end
#         end

#         riduzioni come nel tuo codice (immutate nella logica, ma con meno allocazioni)
#         vecpair = vec(pairsindex)

#         liktaushared = zeros(p*2); idp = 0
#         @inbounds for j in 2:3, i in 1:p
#             idp += 1
#             gi = (gradtaucfa[j,:,1]'); gj = (gradtaucfa[j,:,2]')
#             ei = (gradtauefa[j,:,1]'); ej = (gradtauefa[j,:,2]')
#             liktaushared[idp] = -sum(vec(vcat(gi,gj))[vecpair .== i]) - sum(vec(vcat(ei,ej))[vecpair .== i])
#         end

#         liktau  = zeros(p*(c-1-2)); idp = 0
#         @inbounds for j in 4:c, i in 1:p
#             idp += 1
#             gi = (gradtaucfa[j,:,1]'); gj = (gradtaucfa[j,:,2]')
#             liktau[idp] = -sum(vec(vcat(gi,gj))[vecpair .== i])
#         end

#         liktau2 = zeros(p*(c-1-2)); idp = 0
#         @inbounds for j in 4:c, i in 1:p
#             idp += 1
#             ei = (gradtauefa[j,:,1]'); ej = (gradtauefa[j,:,2]')
#             liktau2[idp] = -sum(vec(vcat(ei,ej))[vecpair .== i])
#         end

#         llamb1 = sum(gradlambci[:,2:c,:,:], dims=2)[:,1,:,:]
#         sumgradls = sum(gradlamb2ci[:,2:c,:,:], dims=2)[:,1,:,:]
#         sumgradpsi = sum(gradpsiefa[2:c,:,:], dims=1)[1,:,:]

#         likls1 = zeros(p,q); llam1 = zeros(p,q); likpsi = zeros(p); likpsiidx = zeros(p)
#         likls  = zeros(p,K); liklsidx = zeros(p,K)

#         @inbounds for i in 1:p
#             @inbounds for qs in 1:q
#                 llam1[i,qs]   = -sum(vec(llamb1[qs,:,:]')[vecpair .== i])
#                 likls1[i,qs]  = -sum(vec(gradlscfaidx[qs,:,:]')[vecpair .== i])
#             end
#             likpsi[i]    = -sum(vec(sumgradpsi')[vecpair .== i])
#             likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair .== i])
#             @inbounds for k in 1:K
#                 likls[i,k]    = -sum(vec(sumgradls[k,:,:]')[vecpair .== i])
#                 liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair .== i])
#             end
#         end

#         likpsitrue = likpsi .+ likpsiidx
#         liklstrue  = likls  .+ liklsidx
#         liklamtrue = llam1  .+ likls1

#         likp = sum(gradphici, dims=3)[:,:,1]
#         likphida = -likp[triu!(trues(size(likp)), 1)]

#         likappa = -sum(gradkappacefa, dims=2)

#         G .= vec(vcat(liktaushared,
#                       liktau,
#                       sum(liklamtrue .* L_str, dims=2),
#                       likphida,
#                       liktau2,
#                       liklstrue[indices],
#                       likpsitrue,
#                       likappa))
#     end

#     =================== OBIETTIVO (F) ===================
#     if F !== nothing
#         NB: qui lasciamo la tua struttura, ma eliminiamo allocazioni inutili
#         tcfa1 = [t1[j, X[i,j]+1] for i in 1:n, j in 1:p]
#         tcfa  = [t1[j, X[i,j]]   for i in 1:n, j in 1:p]
#         tefa1 = [t2[j, X[i,j]+1] for i in 1:n, j in 1:p]
#         tefa  = [t2[j, X[i,j]]   for i in 1:n, j in 1:p]

#         plimix = zeros(Float64, n)

#         @inbounds for col in 1:size(pairsindex,2)
#             i = pairsindex[1,col]; j = pairsindex[2,col]

#             CFA scalare
#             s11_1 = dot(@view L1[i,:], Φ * @view L1[i,:]) + 1.0
#             s22_1 = dot(@view L1[j,:], Φ * @view L1[j,:]) + 1.0
#             s12_1 = dot(@view L1[i,:], Φ * @view L1[j,:])
#             ρ1 = s12_1 / sqrt(s11_1*s22_1)

#             EFA scalare
#             s11_2 = dot(@view L2[i,:], @view L2[i,:]) + ψvec[i]
#             s22_2 = dot(@view L2[j,:], @view L2[j,:]) + ψvec[j]
#             ρ2    = dot(@view L2[i,:], @view L2[j,:]) / sqrt(s11_2*s22_2)

#             sd11 = sqrt(s11_1); sd12 = sqrt(s22_1)
#             sd21 = sqrt(s11_2); sd22 = sqrt(s22_2)

#             vettori n
#             T1v = @view T[:,1,col]; T2v = @view T[:,2,col]

#             @inbounds for r in 1:n
#                 u1l = tcfa[r,i]/sd11; u1u = tcfa1[r,i]/sd11
#                 v1l = tcfa[r,j]/sd12; v1u = tcfa1[r,j]/sd12
#                 p1  = bvnP(u1l, u1u, v1l, v1u, ρ1, false)

#                 u2l = (tefa[r,i]-μ[i])/sd21; u2u = (tefa1[r,i]-μ[i])/sd21
#                 v2l = (tefa[r,j]-μ[j])/sd22; v2u = (tefa1[r,j]-μ[j])/sd22
#                 p2  = bvnP(u2l, u2u, v2l, v2u, ρ2, false)

#                 plimix[r] += log(p1)*T1v[r] + log(p2)*T2v[r]
#             end
#         end

#         return -sum(plimix)
#     end
# end




function fun_grad!(G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trs = (c-1-2)
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
   

    t2s = x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs*p)]
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+trs*p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))] )

    Psi = diagm(x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    kappa = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p)+K)]    #interc2 = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p)]
    # prop = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]

    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
    br2 = hcat(br1[:,1:2],reshape(t2s,p,trs))
    t2=hcat(fill(-Inf,p),br2,fill(Inf,p))

    mu = L2*kappa

    # Deltas = Diagonal(Diagonal((L1*Phi*L1')) .+ 1);
    

    Sig1=(L1*Phi*L1' + diagm(ones(p))) #(Deltas^(-1/2))*(L1*Phi*L1' + diagm(ones(p)))*(Deltas^(-1/2))
    Sig2=L2*L2'+Psi


        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradkappacefa = zeros(Float64,K,size(pairsindex,2))
        gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)

        gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)


        gradlscfaidx = zeros(Float64,q,size(pairsindex,2),2)
        gradlambci = zeros(Float64,q,c+1,size(pairsindex,2),2)
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
            
        idx = 0
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1

                sig1s = hcat(vcat(Sig1[i,i],Sig1[i,j]),vcat(Sig1[j,i],Sig1[j,j]))
                cormat1 = cov2cor(sig1s)
                rho1 = cormat1[1,2]
                sig1sij = diag(sig1s)   
                sig1ij = sqrt.(sig1sij)
                sig1ij13 = sig1ij[1]^(3)
                sig1ij23 = sig1ij[2]^(3)
                L1i=L1[i,:]
                L1j=L1[j,:]

                sig2s = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sig2s)
                rho2 = cormat2[1,2]
                muij = mu[vcat(i,j)]
                sig2sij = diag(sig2s)
                sig2ij = sqrt.(sig2sij)
                sig2ij13 = sig2ij[1]^(3)
                sig2ij23 = sig2ij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]
                Psii=Psi[i,i]
                Psij=Psi[j,j]

                @inbounds for N in 1:n
                    xni = X[N,i]
                    xnj = X[N,j]
                    t1i = t1[i,xni]
                    t1i1 = t1[i,xni+1]
                    t1j = t1[j,xnj]
                    t1j1 = t1[j,xnj+1]
                    t1ist = (t1i)/sig1ij[1]
                    t1i1st = (t1i1)/sig1ij[1]
                    t1jst = (t1j)/sig1ij[2]
                    t1j1st = (t1j1)/sig1ij[2]
                    T1 = T[N,1,idx][1]
                    T2 = T[N,2,idx][1]
                    t2i = t2[i,xni]
                    t2i1 = t2[i,xni+1]
                    t2j = t2[j,xnj]
                    t2j1 = t2[j,xnj+1]
                    t2ist = (t2i-muij[1])/sig2ij[1]
                    t2i1st = (t2i1-muij[1])/sig2ij[1]
                    t2jst = (t2j-muij[2])/sig2ij[2]
                    t2j1st = (t2j1-muij[2])/sig2ij[2]
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                
                    pijicfa = bvnP(t1ist, t1i1st, t1jst, t1j1st, rho1, false)
                    derivNorm1i1 = (T1/pijicfa)*(pdf(Normal(0,1),t1i1st)*(cdf(Normal(0,1),(t1j1st-rho1*t1i1st)/sqrtrij1)-cdf(Normal(0,1),(t1jst-rho1*t1i1st)/sqrtrij1)))
                    derivNorm1i = (T1/pijicfa)*(pdf(Normal(0,1),t1ist)*(cdf(Normal(0,1),(t1jst-rho1*t1ist)/sqrtrij1)-cdf(Normal(0,1),(t1j1st-rho1*t1ist)/sqrtrij1)))
                    derivNorm1j1 = (T1/pijicfa)*(pdf(Normal(0,1),t1j1st)*(cdf(Normal(0,1),(t1i1st-rho1*t1j1st)/sqrtrij1)-cdf(Normal(0,1),(t1ist-rho1*t1j1st)/sqrtrij1)))
                    derivNorm1j = (T1/pijicfa)*(pdf(Normal(0,1),t1jst)*(cdf(Normal(0,1),(t1ist-rho1*t1jst)/sqrtrij1)-cdf(Normal(0,1),(t1i1st-rho1*t1jst)/sqrtrij1)))
                
                    gradtaucfa[xni+1,idx,1] +=  derivNorm1i1/(sig1ij[1])
                    gradtaucfa[xni,idx,1] += derivNorm1i/(sig1ij[1])
                    gradtaucfa[xnj+1,idx,2] += derivNorm1j1/(sig1ij[2])
                    gradtaucfa[xnj,idx,2] += derivNorm1j/(sig1ij[2])


                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))

                    if isnan(derivNormi1)
                        derivNormi1 = 0.
                    end
                    if isnan(derivNormi)
                        derivNormi = 0.
                    end
                    if isnan(derivNormj1)
                        derivNormj1 = 0.
                    end
                    if isnan(derivNormj)
                        derivNormj = 0.
                    end

                    
                    if !(t1ist == -Inf && t1jst == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1jst))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else                  
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1st == Inf && t1j1st == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1j1st))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end

                    if !(t1i1st == Inf && t1jst == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1jst))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end

                    if !(t1ist == -Inf && t1j1st == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1j1st))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    dpidrij1 = (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)#*vcat((-((L1j'*Phi.*t1j1st)./(2*sig1sij[2]))),(-((L1i'*Phi.*t1i1st)./(2*sig1sij[1]))))#vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)#
                

                    gradlscfaidx[1:q,idx,1] .+= dpidrij1*(Phi*L1j)/(sig1ij13*sig1ij[2])#(-Phi*L1i*(L1i'*L1j)/(sig1ij13*sig1ij[2])+Phi*L1j/((sig1ij[1]*sig1ij[2])))
                    gradlscfaidx[1:q,idx,2] .+= dpidrij1*(Phi*L1i)/(sig1ij23*sig1ij[1])#(-Phi*L1j*(L1j'*L1i)/(sig1ij23*sig1ij[1])+Phi*L1i/((sig1ij[2]*sig1ij[1])))
                    
                    gradlambci[1:q,xni+1,idx,1] .+=  derivNorm1i1 *(-((L1i*t1i1st)./(sig1sij[1])))#(-((L1i*t1i1).+(kappak*Psii))./((sig1ij13)))
                    gradlambci[1:q,xnj+1,idx,2] .+=  derivNorm1j1 *(-((L1j*t1j1st)./(sig1sij[2])))#(-((L1j*t1j1).+(kappak*Psij))./((sig1ij23)))
                    gradlambci[1:q,xni,idx,1]   .+=  derivNorm1i  *(-((L1i*t1ist)./(sig1sij[1])))#(-((L1i*t1i).+(kappak*Psii))./((sig1ij13)))
                    gradlambci[1:q,xnj,idx,2]   .+=  derivNorm1j  *(-((L1j*t1jst)./(sig1sij[2])))#(-((L1j*t1j).+(kappak*Psij))./((sig1ij23)))               
                    
                    gradphici[1:q,1:q,idx] .+= dpidrij1*(L1i*L1j')/(sig1ij[1]*sig1ij[2])

                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sig2ij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sig2ij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sig2ij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sig2ij[2])


                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)

                    # gradkappacefa[xni+1,1:K,idx,1] += (derivNormi1*(-L2[i,:]/sig2ij[1]))
                    # gradkappacefa[xnj+1,1:K,idx,2] += (derivNormj1*(-L2[j,:]/sig2ij[2]))
                    # gradkappacefa[xni,1:K,idx,1] += (derivNormi*(-L2[i,:]/sig2ij[1]))
                    # gradkappacefa[xnj,1:K,idx,2] += (derivNormj*(-L2[j,:]/sig2ij[2]))
                    gradkappacefa[1:K,idx] .+= derivNormi1*(-L2[i,:]/sig2ij[1]) + derivNormi*(-L2[i,:]/sig2ij[1]) + derivNormj1*(-L2[j,:]/sig2ij[2]) + derivNormj*(-L2[j,:]/sig2ij[2])
                    
                    gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sig2ij13)*(sig2ij[2])))
                    gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sig2ij23)*(sig2ij[1])))

                    gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij[1]-t2i1))/(2*(sig2ij13))
                    gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij[2]-t2j1))/(2*(sig2ij23))
                    gradpsiefa[xni,idx,1] += (derivNormi*(muij[1]-t2i))/(2*(sig2ij13))
                    gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij[2]-t2j))/(2*(sig2ij23))
                    
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        kappak=kappa[k]
                        gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sig2ij13*sig2ij[2])+L2jk/((sig2ij[1]*sig2ij[2])))#((L2j*Psi[i,i])/(sig2ij13*sig2ij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sig2ij23*sig2ij[1])+L2ik/((sig2ij[2]*sig2ij[1])))#((L2[i,:].*Psi[j,j])./(sig2ij23*sig2ij[1]))

                        gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(-((L2ik*t2i1st)./(sig2sij[1]).+(kappak)./(sig2ij[1])))#(-((L2ik*t2i1st).+(kappak*Psii))./((sig2ij13)))
                        gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(-((L2jk*t2j1st)./(sig2sij[2]).+(kappak)./(sig2ij[2])))#(-((L2jk*t2j1st).+(kappak*Psij))./((sig2ij23)))
                        gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(-((L2ik*t2ist)./(sig2sij[1]).+(kappak)./(sig2ij[1])))#(-((L2ik*t2ist).+(kappak*Psii))./((sig2ij13)))
                        gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(-((L2jk*t2jst)./(sig2sij[2]).+(kappak)./(sig2ij[2])))#(-((L2jk*t2jst).+(kappak*Psij))./((sig2ij23)))
                    end 
                
                end
            end
        end

        vecpair = vec(pairsindex)
        liktaushared = zeros(p*(2));idpcfa = 0   
        @inbounds for j in 2:3
            @inbounds for i in 1:p
                idpcfa += 1
                liktaushared[idpcfa] = -sum(vec(vcat(gradtaucfa[j,:,1]',gradtaucfa[j,:,2]'))[vecpair.==i])-sum(vec(vcat(gradtauefa[j,:,1]',gradtauefa[j,:,2]'))[vecpair.==i])
            end
        end

        liktau = zeros(p*(c-1-2));idpcfa = 0   
        @inbounds for j in 4:c
            @inbounds for i in 1:p
                idpcfa += 1
                liktau[idpcfa] = -sum(vec(vcat(gradtaucfa[j,:,1]',gradtaucfa[j,:,2]'))[vecpair.==i])
            end
        end

        liktau2 = zeros(p*(c-1-2));idpefa = 0   
        @inbounds for j in (4):c
            @inbounds for i in 1:p
                idpefa += 1
                liktau2[idpefa] = -sum(vec(vcat(gradtauefa[j,:,1]',gradtauefa[j,:,2]'))[vecpair.==i])
            end
        end



        likls1 = zeros(p,q); llamb1 = sum(gradlambci[:,2:c,:,:],dims=2)[:,1,:,:]; llam1 = zeros(p,q); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
        likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
        llam1idx  = zeros(p,q); likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        @inbounds for i in 1:p
            @inbounds for qs in 1:q
                llam1[i,qs] = -sum(vec(llamb1[qs,:,:]')[vecpair.==i])
                llam1idx[i,qs] = -sum(vec(gradlscfaidx[qs,:,:]')[vecpair.==i])
            end
            likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
            likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end

        likpsitrue = likpsi.+likpsiidx
        liklstrue = likls.+liklsidx
        liklamtrue = llam1.+llam1idx

        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]

        likappa = -sum(gradkappacefa,dims=2)

        G .= vec(vcat(liktaushared,liktau,sum(liklamtrue.*L_str,dims=2),likphida,liktau2,liklstrue[indices],likpsitrue,likappa))#likpsitrue[2:end]
end






########################################################################### Binary - Latent means EFA - Millsap's parametrization ##################################################################################################

function startvalslmbin(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
    liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K)+ncov,replix);#p-1
    for i in 1:replix
        println(i)

        tau1 = rand(Uniform(-1.5,1.5),p)
        lambda1 = rand(Uniform(0.5,1.3),p);
        phi=rand(LKJ(q,1))
        phi=phi[tril!(trues(size(phi)), -1)]
        lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(0.15,.90),Int((p*K-(K*(K-1)/2))))))
        lambda2=lamb2[lamb2.!=0.]
        kappa=rand(Uniform(-1.2,1.2),K)



        prop = rand(Uniform(0.65,0.95))
        bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.1,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


        param_initials = vcat(tau1,lambda1,phi,lambda2,kappa,bets)

        # first iteration
        plimix = E_stepMixlmbin(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        params[:,i] = M_parsMixlmbin(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        dt = time() - t
        prrot=params[:,i]
        liksmix[i] = E_stepMixlmbin(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false);
    end
    param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    lik=liksmix[findall(liksmix.==maximum(liksmix))]

    tau1=param_initials0[(1):((c-1)*p)]
    tau1=vec(sort(reshape(tau1,p,c-1),dims=2))
    lambda1=param_initials0[((c-1)*p+1):((c-1)*p+p)]
    phi=param_initials0[((c-1)*p+p+1):((c-1)*p+p+Int(q*(q-1)/2))]
    # lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    # lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    # phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    lambda2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2))))]
    kappa=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+K)]#interc=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    # lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    # lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+K+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+K+ncov)]#    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


    param_initialii = vcat(tau1,lambda1,phi,lambda2,kappa,bets)

    println(param_initialii,liksmix,"first iteration: done")
    return param_initialii
end


function EMparsMixlmbin(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
    local plimix::Array, T,oldlik
    n = size(X,1)
    er=0;
    # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
    it=0;
    dev=Inf; dev1=Inf; dev2=Inf;

    
    oldlik=-1000000.;
    oldlik2=1000000.;
    oldlikA=-1000000.;
    lik=0.
    

    noerror = true; devnostop = true;  sum_devs = 0;
    while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
        it = it+1;
        param=param01;

        T00 = 1
        plimix = E_stepMixlmbin(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        param01 = M_parsMixlmbin(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
        dt = time() - t
        if param01 == false
            return return(vcat(param01,vec(plimix),vec(T),0.,0.))
        end
        lik = E_stepMixlmbin(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false);

        #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
        con = (lik-oldlik)/(oldlik-oldlik2);
        likA = oldlik + (lik-oldlik)/(1-con)
        old_dev = dev1
        dev1 = abs(likA-oldlikA)
        dev2 = abs(lik-oldlik)
        println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
        println(param01)
        dev = abs(lik-oldlik)/abs(oldlik)

        noerror = oldlik<lik

        oldlikA=likA
        oldlik2=oldlik
        oldlik=lik;
        pml_em_log("EMparsMixlmbin", it, lik, dev, dev1, abs(lik - oldlik2), tol, noerror)

        # sum_devs += dev > old_dev
        # devstop = sum_devs <= 3
    end
    
    # classificazioni dopo ultima massimizzazione
    plimix = E_stepMixlmbin(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

    if noerror == true && it<=maxit 
        return vcat(param01,vec(plimix),vec(T),lik,1.)
    else
        return vcat(param01,vec(plimix),vec(T),0.,0.)
    end
end





function E_stepMixlmbin(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep)#::AbstractVector{T},q::Int64,K::Int64,p::Int64,n::Int64,c::Int64,X::Matrix{Int64},ics::Matrix,ncov::Int64,pairsindex::Matrix{Int64},L_str::Matrix,indices,estep::Bool) where T
    # local trs::Int64,trsi::Int64,t1s::Vector{Real}, L1::Matrix{Real}, Phi::Matrix{Real},  t2s::Vector{Real}, L2::Matrix{Real}, Psi::Matrix{Real}, kappa::Vector{Real}, betas::Vector{Real}, propo::Vector{Real}, pg::Matrix{Real}, mu::Vector{Real}, Sigms1::Matrix{Real}, Sigms2::Matrix{Real}, br1::Matrix{Real}, br2::Matrix{Real}, t1::Matrix{Real}, t2::Matrix{Real}, tcfa1::Matrix{Real}, tcfa::Matrix{Real}, tefa1::Matrix{Real}, tefa::Matrix{Real}, plimix::Array{Real}, plimixll::Vector{Real}, Scfa::Matrix{Real}, Rcfa::Matrix{Real}, Bcfa::Matrix{Real}, tcfa3::Vector{Real}, Sefa::Matrix{Real}, Refa::Matrix{Real}, Befa::Matrix{Real}, tefa3::Vector{Real}, cllik::Real
    if sum(param .== 0.) == size(param,1)
        return -Inf
    end
    trsi = (c-1)
    t1s = param[1:trsi*p]
    L1 = param[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    
    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )

    
    kappa = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K)]   # interc2 = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
    betas = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K+ncov)]  #  betas = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)

    mu = L2*kappa
    
    # Sigms=zeros(p,p,2);
    # Deltas = Diagonal(Diagonal((L1*Phi*L1')) .+ 1);
    # Sigms[:,:,1]=(Deltas^(-1/2))*(L1*Phi*L1' + diagm(ones(p)))*(Deltas^(-1/2))
    Sigms1=L1*Phi*L1' + diagm(ones(p))
    Sigms2=L2*L2'+ diagm(ones(p))
    
    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))


    tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
    


    # tlow = cat(tcfa,tefa,dims=3)
    # tup = cat(tcfa1,tefa1,dims=3)
    
    # if estep == true
    #     low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
    #     plimix=ones(n,2,size(pairsindex,2));
    
    #     for idx=1:size(pairsindex,2) #parallelizzabile
    #         for g=1:2
    #             if g==1
    #                 s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
    #                 r2c=cov2cor(s2c);
    #                 up=(tcfa1[:,pairsindex[:,idx]]./sqrt.(diag(s2c)'))
    #                 low=(tcfa[:,pairsindex[:,idx]]./sqrt.(diag(s2c)'))
    #             elseif g==2
    #                 s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
    #                 r2c=cov2cor(s2c);
    #                 up=(tefa1[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
    #                 low=(tefa[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
    #             end
    #             B=hcat(low,up);
    #             t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
    #             plimix[:,g,idx]= t3.*pg[:,g];
    #         end
    #     end
    # # % probabilities in each row must sum to 1 (application of Bayes law)
    #     return plimix
    # elseif estep == false
        # low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        # plimix=ones(n,2,size(pairsindex,2));

        # for idx=1:size(pairsindex,2) #parallelizzabile
        #     for g=1:2
        #         if g==1
        #             s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
        #             r2c=cov2cor(s2c);
        #             up=(tcfa1[:,pairsindex[:,idx]]./sqrt.(diag(s2c)'))
        #             low=(tcfa[:,pairsindex[:,idx]]./sqrt.(diag(s2c)'))
        #         elseif g==2
        #             s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
        #             r2c=cov2cor(s2c);
        #             up=(tefa1[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
        #             low=(tefa[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
        #         end
        #         B=hcat(low,up);
        #         t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
        #         plimix[:,g,idx] = t3
        #     end
        # end
        # lik1=(log.(plimix));
        # lik1=sum(lik1,dims=3);
        # U=exp.(lik1)[:,:,1].*pg; 
        # cllik=sum(log.(sum(U,dims=2)))
    #     return cllik
    # end
    if estep == true
        plimix=zeros(n,2,size(pairsindex,2));
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sigms1[pairsindex[:,idx],pairsindex[:,idx]]
            Rcfa=cov2cor(Scfa);
            Bcfa=hcat(tfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            
            
            #EFA
            Sefa=Sigms2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            
    
            plimix[:,1,idx] = tcfa3 .* pg[:,1]
            plimix[:,2,idx] = tefa3 .* pg[:,2]
        end
    # % probabilities in each row must sum to 1 (application of Bayes law)
        return plimix
    elseif estep == false
        plimixll=zeros(Float64,n);
        # plimix=zeros(n,2,size(pairsindex,2));
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sigms1[pairsindex[:,idx],pairsindex[:,idx]]
            Rcfa=cov2cor(Scfa);
            sdScfa=sqrt.(diag(Scfa)')
            tcfaidx = tfa[:,pairsindex[:,idx]]
            tcfa1idx = tfa1[:,pairsindex[:,idx]]
            Bcfa=hcat(tcfaidx./sdScfa,tcfa1idx./sdScfa);
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            
            
            #EFA
            Sefa=Sigms2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            sdSefa=sqrt.(diag(Sefa)')
            muij=mu[pairsindex[:,idx]]'
            tefaidx = tfa[:,pairsindex[:,idx]]
            tefa1idx = tfa1[:,pairsindex[:,idx]]
            Befa=hcat((tefaidx.-muij)./sdSefa,(tefa1idx.-muij)./sdSefa);
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            
            # plimix[:,1,idx] = tcfa3 .* pg[:,1]
            # plimix[:,2,idx] = tefa3 .* pg[:,2]
            plimixll += log.(tcfa3.*pg[:,1].+tefa3.*pg[:,2])
        end
        # cllik=sum(log.(sum(plimix,dims=2))) 
        cllik=sum(plimixll) 
        return cllik
    end
end
    



function obslikbin(param::AbstractVector{T},q::Int64,K::Int64,p::Int64,n::Int64,c::Int64,X::Matrix{Int64},ics::Matrix,ncov::Int64,pairsindex::Matrix{Int64},L_str::Matrix,indices::Vector,estep::Bool) where T
    local trs::Int64,trsi::Int64,t1s::Vector{Real}, L1::Matrix{Real}, Phi::Matrix{Real},  t2s::Vector{Real}, L2::Matrix{Real}, Psi::Matrix{Real}, kappa::Vector{Real}, betas::Vector{Real}, propo::Vector{Real}, pg::Matrix{Real}, mu::Vector{Real}, Sigms1::Matrix{Real}, Sigms2::Matrix{Real}, br1::Matrix{Real}, br2::Matrix{Real}, t1::Matrix{Real}, t2::Matrix{Real}, tcfa1::Matrix{Real}, tcfa::Matrix{Real}, tefa1::Matrix{Real}, tefa::Matrix{Real}, plimixll::Vector{Real}, Scfa::Matrix{Real}, Rcfa::Matrix{Real}, Bcfa::Matrix{Real}, tcfa3::Vector{Real}, Sefa::Matrix{Real}, Refa::Matrix{Real}, Befa::Matrix{Real}, tefa3::Vector{Real}, cllik::Real
    # if sum(param .== 0.) == size(param,1)
    #     return -Inf
    # end
    trsi = (c-1)*p

    t1s = param[1:trsi]
    L1 = param[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
    
    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )

    
    kappa = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K)]   # interc2 = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
    betas = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K+ncov)]  #  betas = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)

    mu = L2*kappa
    
    # Sigms=zeros(p,p,2);
    # Deltas = Diagonal(Diagonal((L1*Phi*L1')) .+ 1);
    # Sigms[:,:,1]=(Deltas^(-1/2))*(L1*Phi*L1' + diagm(ones(p)))*(Deltas^(-1/2))
    Sigms1=L1*Phi*L1' + diagm(ones(p))
    Sigms2=L2*L2'+ diagm(ones(p))
    
    br1 = reshape(t1s,p,c-1)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))


    tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]

        plimixll=zeros(Float64,n);
        # plimix=zeros(n,2,size(pairsindex,2));
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sigms1[pairsindex[:,idx],pairsindex[:,idx]]
            Rcfa=cov2cor(Scfa);
            sdScfa=sqrt.(diag(Scfa)')
            tcfaidx = tfa[:,pairsindex[:,idx]]
            tcfa1idx=tfa1[:,pairsindex[:,idx]]
            Bcfa=hcat(tcfaidx./sdScfa,tcfa1idx./sdScfa);
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            
            
            #EFA
            Sefa=Sigms2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            sdSefa=sqrt.(diag(Sefa)')
            muij=mu[pairsindex[:,idx]]'
            tefaidx = tfa[:,pairsindex[:,idx]]
            tefa1idx = tfa1[:,pairsindex[:,idx]]
            Befa=hcat((tefaidx.-muij)./sdSefa,(tefa1idx.-muij)./sdSefa);
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            
            # plimix[:,1,idx] = tcfa3 .* pg[:,1]
            # plimix[:,2,idx] = tefa3 .* pg[:,2]
            plimixll += log.(tcfa3.*pg[:,1].+tefa3.*pg[:,2])
        end
        # cllik=sum(log.(sum(plimix,dims=2))) 
        cllik=sum(plimixll) 
        return cllik
end





function M_parsMixlmbin(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)

    trsi = (c-1)*p #attenzione p è già qui

    # PML - CFA & EFA
    theta_init = theta[1:(size(theta,1)-ncov)]#0.000001
    # # -1 loadings
    # lb=vec(vcat(fill(-5.1,trsi),fill(0.00001,p),fill(-1.,Int((q*(q-1))/2)),fill(1.,trs),fill(-1.,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-10,p)));
    # # 0. loadings
    lb=vec(vcat(fill(-Inf,trsi),fill(-Inf,p),fill(-1.,Int((q*(q-1))/2)),fill(-Inf,Int(p*K-K*(K-1)/2)),fill(-Inf,K)));#1.0e-5
    ub=vec(vcat(fill(Inf,trsi),fill(Inf,p),fill(1.,Int((q*(q-1))/2)),fill(Inf,Int(p*K-K*(K-1)/2)),fill(Inf,K)));

    par = zeros(Float64,size(lb))
    if it == 1
        par = optimize(Optim.only_fg!((F,G,x)->fglmbin!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min ,Optim.Options(iterations=1000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0))
    elseif it > 1
        par = optimize(Optim.only_fg!((F,G,x)->fglmbin!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min
    end 

    # func = TwiceDifferentiable((x)->pmllastlmbin(x,X,T,p,q,K,c,n,pairsindex,L_str,indices),(G,x)->fun_grad!(G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices),theta_init; autodiff=:forward);
    # dfc = TwiceDifferentiableConstraints(lb, ub)
    # par = optimize(func,dfc,theta_init,IPNewton(),Optim.Options(iterations=10000,f_reltol=1.0000e-06))


    # par = optimize(x->pmllast(x,X,T,p,q,K,c,n,pairsindex,L_str,indices),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

    # par = zeros(Float64,size(lb))
    # if it == 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 300.0)) #20 min
    # elseif it > 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 1200.0)) #20 min
    # end 
    
    # mixture parameter
    betas = theta[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
    Ez1 = T[:,1,:]
    Ez0 = T[:,2,:]
    
    pairs = size(pairsindex,2)
    betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)

    
    # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[1,1];
    # betasnew = vcat(log(pg/(1 - pg)))
   
    # println(par.minimizer,pg[1])
    if !(Optim.converged(par)) && isnan(betas[1]) 
        println(betas[1])
        if start == true
            return zeros(Float64,size(theta,1))
        elseif start == false
            return Optim.converged(par)
        end
    end
    return vcat(par.minimizer,betasnew)#,pg[1])
end



function pmllastlmbin(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])


    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )

    kappa = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K)]    #interc2 = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
    # prop = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+p+1)]

    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))

    mu = L2*kappa
    
    # Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    
    # Deltas = Diagonal(Diagonal((L1*Phi*L1')) .+ 1);
    # Sig1=(Deltas^(-1/2))*(L1*Phi*L1' + diagm(ones(p)))*(Deltas^(-1/2))
    
    Sig1=L1*Phi*L1' + diagm(ones(p))
    Sig2=L2*L2' + diagm(ones(p))



    plimix=zeros(Float64,n);
    tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
    

    for idx=1:size(pairsindex,2) #parallelizzabile
        # CFA
        Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
        Rcfa=cov2cor(Scfa);
        Bcfa=hcat(tfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
        tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
        Tt1 = T[:,1,idx]
        
        #EFA
        Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
        Refa=cov2cor(Sefa);
        Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
        tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
        Tt2 = T[:,2,idx]

        plimix += log.(tcfa3).*Tt1.+log.(tefa3).*Tt2
    end
    f=-sum(plimix)
return f
end

function fglmbin!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trs = (c-1-2)
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
   

    t2s = x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs*p)]
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )

    kappa = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K)]    #interc2 = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
    # prop = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+p+1)]

    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
   

    mu = L2*kappa

    # Deltas = Diagonal(Diagonal((L1*Phi*L1')) .+ 1);
    

    Sig1=(L1*Phi*L1' + diagm(ones(p))) #(Deltas^(-1/2))*(L1*Phi*L1' + diagm(ones(p)))*(Deltas^(-1/2))
    Sig2=L2*L2' + diagm(ones(p))

    if G !== nothing
      
        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradkappacefa = zeros(Float64,K,size(pairsindex,2))
        

        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)


        gradlscfaidx = zeros(Float64,q,size(pairsindex,2),2)
        gradlambci = zeros(Float64,q,c+1,size(pairsindex,2),2)
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
            
        idx = 0
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1

                sig1s = hcat(vcat(Sig1[i,i],Sig1[i,j]),vcat(Sig1[j,i],Sig1[j,j]))
                cormat1 = cov2cor(sig1s)
                rho1 = cormat1[1,2]
                sig1sij = diag(sig1s)   
                sig1ij = sqrt.(sig1sij)
                sig1ij13 = sig1ij[1]^(3)
                sig1ij23 = sig1ij[2]^(3)
                L1i=L1[i,:]
                L1j=L1[j,:]

                sig2s = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sig2s)
                rho2 = cormat2[1,2]
                muij = mu[vcat(i,j)]
                sig2sij = diag(sig2s)
                sig2ij = sqrt.(sig2sij)
                sig2ij13 = sig2ij[1]^(3)
                sig2ij23 = sig2ij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]

                @inbounds for N in 1:n
                    xni = X[N,i]
                    xnj = X[N,j]
                    t1i = t1[i,xni]
                    t1i1 = t1[i,xni+1]
                    t1j = t1[j,xnj]
                    t1j1 = t1[j,xnj+1]
                    t1ist = (t1i)/sig1ij[1]
                    t1i1st = (t1i1)/sig1ij[1]
                    t1jst = (t1j)/sig1ij[2]
                    t1j1st = (t1j1)/sig1ij[2]
                    T1 = T[N,1,idx][1]
                    T2 = T[N,2,idx][1]
                    t2i = t1[i,xni]
                    t2i1 = t1[i,xni+1]
                    t2j = t1[j,xnj]
                    t2j1 = t1[j,xnj+1]
                    t2ist = (t2i-muij[1])/sig2ij[1]
                    t2i1st = (t2i1-muij[1])/sig2ij[1]
                    t2jst = (t2j-muij[2])/sig2ij[2]
                    t2j1st = (t2j1-muij[2])/sig2ij[2]
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                
                    pijicfa = bvnP(t1ist, t1i1st, t1jst, t1j1st, rho1, false)
                    derivNorm1i1 = (T1/pijicfa)*(pdf(Normal(0,1),t1i1st)*(cdf(Normal(0,1),(t1j1st-rho1*t1i1st)/sqrtrij1)-cdf(Normal(0,1),(t1jst-rho1*t1i1st)/sqrtrij1)))
                    derivNorm1i = (T1/pijicfa)*(pdf(Normal(0,1),t1ist)*(cdf(Normal(0,1),(t1jst-rho1*t1ist)/sqrtrij1)-cdf(Normal(0,1),(t1j1st-rho1*t1ist)/sqrtrij1)))
                    derivNorm1j1 = (T1/pijicfa)*(pdf(Normal(0,1),t1j1st)*(cdf(Normal(0,1),(t1i1st-rho1*t1j1st)/sqrtrij1)-cdf(Normal(0,1),(t1ist-rho1*t1j1st)/sqrtrij1)))
                    derivNorm1j = (T1/pijicfa)*(pdf(Normal(0,1),t1jst)*(cdf(Normal(0,1),(t1ist-rho1*t1jst)/sqrtrij1)-cdf(Normal(0,1),(t1i1st-rho1*t1jst)/sqrtrij1)))
                
                    gradtaucfa[xni+1,idx,1] +=  derivNorm1i1/(sig1ij[1])
                    gradtaucfa[xni,idx,1] += derivNorm1i/(sig1ij[1])
                    gradtaucfa[xnj+1,idx,2] += derivNorm1j1/(sig1ij[2])
                    gradtaucfa[xnj,idx,2] += derivNorm1j/(sig1ij[2])


                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))

                    if isnan(derivNormi1)
                        derivNormi1 = 0.
                    end
                    if isnan(derivNormi)
                        derivNormi = 0.
                    end
                    if isnan(derivNormj1)
                        derivNormj1 = 0.
                    end
                    if isnan(derivNormj)
                        derivNormj = 0.
                    end

                    
                    if !(t1ist == -Inf && t1jst == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1jst))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else                  
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1st == Inf && t1j1st == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1j1st))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end

                    if !(t1i1st == Inf && t1jst == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1jst))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end

                    if !(t1ist == -Inf && t1j1st == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1j1st))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    dpidrij1 = (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)#*vcat((-((L1j'*Phi.*t1j1st)./(2*sig1sij[2]))),(-((L1i'*Phi.*t1i1st)./(2*sig1sij[1]))))#vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)#
                

                    gradlscfaidx[1:q,idx,1] .+= dpidrij1*(Phi*L1j)/(sig1ij13*sig1ij[2])#(-Phi*L1i*(L1i'*L1j)/(sig1ij13*sig1ij[2])+Phi*L1j/((sig1ij[1]*sig1ij[2])))
                    gradlscfaidx[1:q,idx,2] .+= dpidrij1*(Phi*L1i)/(sig1ij23*sig1ij[1])#(-Phi*L1j*(L1j'*L1i)/(sig1ij23*sig1ij[1])+Phi*L1i/((sig1ij[2]*sig1ij[1])))
                    
                    gradlambci[1:q,xni+1,idx,1] .+=  derivNorm1i1 *(-((L1i*t1i1st)./(sig1sij[1])))#(-((L1i*t1i1).+(kappak*Psii))./((sig1ij13)))
                    gradlambci[1:q,xnj+1,idx,2] .+=  derivNorm1j1 *(-((L1j*t1j1st)./(sig1sij[2])))#(-((L1j*t1j1).+(kappak*Psij))./((sig1ij23)))
                    gradlambci[1:q,xni,idx,1]   .+=  derivNorm1i  *(-((L1i*t1ist)./(sig1sij[1])))#(-((L1i*t1i).+(kappak*Psii))./((sig1ij13)))
                    gradlambci[1:q,xnj,idx,2]   .+=  derivNorm1j  *(-((L1j*t1jst)./(sig1sij[2])))#(-((L1j*t1j).+(kappak*Psij))./((sig1ij23)))               
                    
                    gradphici[1:q,1:q,idx] .+= dpidrij1*(L1i*L1j')/(sig1ij[1]*sig1ij[2])

                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sig2ij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sig2ij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sig2ij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sig2ij[2])


                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)

                    # gradkappacefa[xni+1,1:K,idx,1] += (derivNormi1*(-L2[i,:]/sig2ij[1]))
                    # gradkappacefa[xnj+1,1:K,idx,2] += (derivNormj1*(-L2[j,:]/sig2ij[2]))
                    # gradkappacefa[xni,1:K,idx,1] += (derivNormi*(-L2[i,:]/sig2ij[1]))
                    # gradkappacefa[xnj,1:K,idx,2] += (derivNormj*(-L2[j,:]/sig2ij[2]))
                    gradkappacefa[1:K,idx] .+= derivNormi1*(-L2[i,:]/sig2ij[1]) + derivNormi*(-L2[i,:]/sig2ij[1]) + derivNormj1*(-L2[j,:]/sig2ij[2]) + derivNormj*(-L2[j,:]/sig2ij[2])
                    
                   
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        kappak=kappa[k]
                        gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sig2ij13*sig2ij[2])+L2jk/((sig2ij[1]*sig2ij[2])))#((L2j*Psi[i,i])/(sig2ij13*sig2ij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sig2ij23*sig2ij[1])+L2ik/((sig2ij[2]*sig2ij[1])))#((L2[i,:].*Psi[j,j])./(sig2ij23*sig2ij[1]))

                        gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(-((L2ik*t2i1st)./(sig2sij[1]).+(kappak)./(sig2ij[1])))#(-((L2ik*t2i1st).+(kappak*Psii))./((sig2ij13)))
                        gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(-((L2jk*t2j1st)./(sig2sij[2]).+(kappak)./(sig2ij[2])))#(-((L2jk*t2j1st).+(kappak*Psij))./((sig2ij23)))
                        gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(-((L2ik*t2ist)./(sig2sij[1]).+(kappak)./(sig2ij[1])))#(-((L2ik*t2ist).+(kappak*Psii))./((sig2ij13)))
                        gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(-((L2jk*t2jst)./(sig2sij[2]).+(kappak)./(sig2ij[2])))#(-((L2jk*t2jst).+(kappak*Psij))./((sig2ij23)))
                    end 
                
                end
            end
        end

        vecpair = vec(pairsindex)
        liktaushared = zeros(p*(c-1));idpcfa = 0   
        @inbounds for j in 2:c
            @inbounds for i in 1:p
                idpcfa += 1
                liktaushared[idpcfa] = -sum(vec(vcat(gradtaucfa[j,:,1]',gradtaucfa[j,:,2]'))[vecpair.==i])-sum(vec(vcat(gradtauefa[j,:,1]',gradtauefa[j,:,2]'))[vecpair.==i])
            end
        end

        # liktau = zeros(p*(c-1-2));idpcfa = 0   
        # @inbounds for j in 4:c
        #     @inbounds for i in 1:p
        #         idpcfa += 1
        #         liktau[idpcfa] = -sum(vec(vcat(gradtaucfa[j,:,1]',gradtaucfa[j,:,2]'))[vecpair.==i])  
        #     end
        # end

        # liktau2 = zeros(p*(c-1-2));idpefa = 0   
        # @inbounds for j in (4):c
        #     @inbounds for i in 1:p
        #         idpefa += 1
        #         liktau2[idpefa] = -sum(vec(vcat(gradtauefa[j,:,1]',gradtauefa[j,:,2]'))[vecpair.==i])
        #     end
        # end



        likls1 = zeros(p,q); llamb1 = sum(gradlambci[:,2:c,:,:],dims=2)[:,1,:,:]; llam1 = zeros(p,q); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
        llam1idx  = zeros(p,q); likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        @inbounds for i in 1:p
            @inbounds for qs in 1:q
                llam1[i,qs] = -sum(vec(llamb1[qs,:,:]')[vecpair.==i])
                llam1idx[i,qs] = -sum(vec(gradlscfaidx[qs,:,:]')[vecpair.==i])
            end
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end

        liklstrue = likls.+liklsidx
        liklamtrue = llam1.+llam1idx

        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]

        likappa = -sum(gradkappacefa,dims=2)

        G .= vec(vcat(liktaushared,sum(liklamtrue.*L_str,dims=2),likphida,liklstrue[indices],likappa))#likpsitrue[2:end]

    end
    if F !== nothing 
  
    
        plimix=zeros(Float64,n);
        tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
        tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
        
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
            Rcfa=cov2cor(Scfa);
            Bcfa=hcat(tfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            Tt1 = T[:,1,idx]
            
            #EFA
            Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            Tt2 = T[:,2,idx]

            plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
        end
        f=-sum(plimix)
        # @info "Current objective: $f"
        return f
    end
end

















########################################################################### Latent means EFA + invariant thresholds + Psi estimated ##################################################################################################

function startvalsinv(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
    liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)+ncov,replix);#p-1
    for i in 1:replix
        # println(i)

        tau1 = vec(hcat(rand(Uniform(-2.0,-1.5),p).+reshape(repeat(collect(range(0,3.5,c-1)),inner=p),p,(c-1))))
        lambda1 = rand(Uniform(0.3,1.40),p);
        phi=rand(LKJ(q,1))
        phi=phi[tril!(trues(size(phi)), -1)]
        lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(-1.20,1.20),Int((p*K-(K*(K-1)/2))))))
        lambda2=lamb2[lamb2.!=0.]
        psi=rand(Uniform(0.2,2.5),p)
        kappa=rand(Uniform(-1.2,1.2),K)



        prop = rand(Uniform(0.75,0.95))
        bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.1,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


        param_initials = vcat(tau1,lambda1,phi,lambda2,psi,kappa,bets)

        # first iteration
        plimix = E_stepMixinv(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        params[:,i] = M_parsMixinv(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        dt = time() - t
        prrot=params[:,i]
        liksmix[i] = E_stepMixinv(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
    end
    param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    lik=liksmix[findall(liksmix.==maximum(liksmix))]

    tau1=param_initials0[(1):((c-1)*p)]
    tau1=vec(sort(reshape(tau1,p,c-1),dims=2))
    lambda1=param_initials0[((c-1)*p+1):((c-1)*p+p)]
    phi=param_initials0[((c-1)*p+p+1):((c-1)*p+p+Int(q*(q-1)/2))]
    # lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    # lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    # phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    lambda2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2))))]
    psi=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p))]# psi=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1))]
    kappa=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+K)]#interc=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    # lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    # lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+K+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+K+ncov)]#    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


    param_initialii = vcat(tau1,lambda1,phi,lambda2,psi,kappa,bets)

    # println(param_initialii,liksmix,"first iteration: done")
    return param_initialii
end


function EMparsMixinv(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
    local plimix::Array, T,oldlik
    n = size(X,1)
    er=0;
    # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
    it=0;
    dev=Inf; dev1=Inf; dev2=Inf;

    
    oldlik=-1000000.;
    oldlik2=1000000.;
    oldlikA=-1000000.;
    lik=0.
    

    noerror = true; devnostop = true;  sum_devs = 0;
    while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
        it = it+1;
        param=param01;

        T00 = 1
        plimix = E_stepMixinv(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        param01 = M_parsMixinv(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
        dt = time() - t
        if param01 == false
            return return(vcat(param01,vec(plimix),vec(T),0.,0.))
        end
        lik = E_stepMixinv(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);

        #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
        con = (lik-oldlik)/(oldlik-oldlik2);
        likA = oldlik + (lik-oldlik)/(1-con)
        old_dev = dev1
        dev1 = abs(likA-oldlikA)
        dev2 = abs(lik-oldlik)
        # println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
        # println(param01)
        dev = abs(lik-oldlik)/abs(oldlik)

        noerror = oldlik<lik

        oldlikA=likA
        oldlik2=oldlik
        oldlik=lik;
        pml_em_log("EMparsMixinv", it, lik, dev, dev1, abs(lik - oldlik2), tol, noerror)

        # sum_devs += dev > old_dev
        # devstop = sum_devs <= 3
    end
    
    # classificazioni dopo ultima massimizzazione
    plimix = E_stepMixinv(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

    if noerror == true && it<=maxit 
        return vcat(param01,vec(plimix),vec(T),lik,1.)
    else
        return vcat(param01,vec(plimix),vec(T),0.,0.)
    end
end





function E_stepMixinv(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
    if sum(param .== 0.) == size(param,1)
        return -Inf
    end
    
    trsi = (c-1)
    t1s = param[1:trsi*p]
    L1 = param[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    
    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )
    Psi = diagm(param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    kappa = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)]   # interc2 = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
    betas = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+ncov)]  #  betas = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)

    mu = L2*kappa
    
    
    Sigms1=L1*Phi*L1' + diagm(ones(p))
    Sigms2=L2*L2'+Psi 
    
    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
   

    tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
   
    
    
    if estep == true
        plimix=zeros(n,2,size(pairsindex,2));
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sigms1[pairsindex[:,idx],pairsindex[:,idx]]
            Rcfa=cov2cor(Scfa);
            Bcfa=hcat(tfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            
            
            #EFA
            Sefa=Sigms2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            

            plimix[:,1,idx] = tcfa3 .* pg[:,1]
            plimix[:,2,idx] = tefa3 .* pg[:,2]
        end
    # % probabilities in each row must sum to 1 (application of Bayes law)
        return plimix
    elseif estep == false
        plimixll=zeros(Float64,n);
        # plimix=zeros(n,2,size(pairsindex,2));
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sigms1[pairsindex[:,idx],pairsindex[:,idx]]
            Rcfa=cov2cor(Scfa);
            sdScfa=sqrt.(diag(Scfa)')
            tcfaidx = tfa[:,pairsindex[:,idx]]
            tcfa1idx=tfa1[:,pairsindex[:,idx]]
            Bcfa=hcat(tcfaidx./sdScfa,tcfa1idx./sdScfa);
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            
            
            #EFA
            Sefa=Sigms2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            sdSefa=sqrt.(diag(Sefa)')
            muij=mu[pairsindex[:,idx]]'
            tefaidx = tfa[:,pairsindex[:,idx]]
            tefa1idx = tfa1[:,pairsindex[:,idx]]
            Befa=hcat((tefaidx.-muij)./sdSefa,(tefa1idx.-muij)./sdSefa);
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            
            # plimix[:,1,idx] = tcfa3 .* pg[:,1]
            # plimix[:,2,idx] = tefa3 .* pg[:,2]
            plimixll += log.(tcfa3.*pg[:,1].+tefa3.*pg[:,2])
        end
        # cllik=sum(log.(sum(plimix,dims=2))) 
        cllik=sum(plimixll) 
        return cllik
    end
end
    



function obslikinv(param::AbstractVector{T},q::Int64,K::Int64,p::Int64,n::Int64,c::Int64,X::Matrix{Int64},ics::Matrix,ncov::Int64,pairsindex::Matrix{Int64},L_str::Matrix,indices::Vector,estep::Bool) where T
    local trs::Int64,trsi::Int64,t1s::Vector{Real}, L1::Matrix{Real}, Phi::Matrix{Real},   L2::Matrix{Real}, Psi::Matrix{Real}, kappa::Vector{Real}, betas::Vector{Real}, propo::Vector{Real}, pg::Matrix{Real}, mu::Vector{Real}, Sigms1::Matrix{Real}, Sigms2::Matrix{Real}, br1::Matrix{Real}, t1::Matrix{Real}, tfa1::Matrix{Real}, tfa::Matrix{Real}, plimixll::Vector{Real}, Scfa::Matrix{Real}, Rcfa::Matrix{Real}, Bcfa::Matrix{Real}, tcfa3::Vector{Real}, Sefa::Matrix{Real}, Refa::Matrix{Real}, Befa::Matrix{Real}, tefa3::Vector{Real}, cllik::Real
    # if sum(param .== 0.) == size(param,1)
    #     return -Inf
    # end
    trsi = (c-1)
    t1s = param[1:trsi*p]
    L1 = param[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    
    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )
    Psi = diagm(param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    kappa = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)]   # interc2 = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
    betas = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+ncov)]  #  betas = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)

    mu = L2*kappa
    
    
    Sigms1=L1*Phi*L1' + diagm(ones(p))
    Sigms2=L2*L2'+Psi 
    
    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
   

    tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
   
    
        plimixll=zeros(Float64,n);
        # plimix=zeros(n,2,size(pairsindex,2));
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sigms1[pairsindex[:,idx],pairsindex[:,idx]]
            Rcfa=cov2cor(Scfa);
            sdScfa=sqrt.(diag(Scfa)')
            tcfaidx = tfa[:,pairsindex[:,idx]]
            tcfa1idx=tfa1[:,pairsindex[:,idx]]
            Bcfa=hcat(tcfaidx./sdScfa,tcfa1idx./sdScfa);
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            
            
            #EFA
            Sefa=Sigms2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            sdSefa=sqrt.(diag(Sefa)')
            muij=mu[pairsindex[:,idx]]'
            tefaidx = tfa[:,pairsindex[:,idx]]
            tefa1idx = tfa1[:,pairsindex[:,idx]]
            Befa=hcat((tefaidx.-muij)./sdSefa,(tefa1idx.-muij)./sdSefa);
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            
            # plimix[:,1,idx] = tcfa3 .* pg[:,1]
            # plimix[:,2,idx] = tefa3 .* pg[:,2]
            plimixll += log.(tcfa3.*pg[:,1].+tefa3.*pg[:,2])
        end
        # cllik=sum(log.(sum(plimix,dims=2))) 
        cllik=sum(plimixll) 
        return cllik
end



function M_parsMixinv(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)

    trsi = (c-1)*p #attenzione p è già qui
   
    # PML - CFA & EFA
    theta_init = theta[1:(size(theta,1)-ncov)]#0.000001
    # # -1 loadings
    # lb=vec(vcat(fill(-5.1,trsi),fill(0.00001,p),fill(-1.,Int((q*(q-1))/2)),fill(1.,trs),fill(-1.,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-10,p)));
    # # 0. loadings
    lb=vec(vcat(fill(-Inf,trsi),fill(-Inf,p),fill(-1.,Int((q*(q-1))/2)),fill(-Inf,Int(p*K-K*(K-1)/2)),fill(1.0e-6,p),fill(-Inf,K)));#1.0e-5
    ub=vec(vcat(fill(Inf,trsi),fill(Inf,p),fill(1.,Int((q*(q-1))/2)),fill(Inf,Int(p*K-K*(K-1)/2)),fill(Inf,p),fill(Inf,K)));

    par = zeros(Float64,size(lb))
    if it == 1
        par = optimize(Optim.only_fg!((F,G,x)->fginv!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min ,Optim.Options(iterations=1000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0))
    elseif it > 1
        par = optimize(Optim.only_fg!((F,G,x)->fginv!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min
    end 


    # par = optimize(x->pmllast(x,X,T,p,q,K,c,n,pairsindex,L_str,indices),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

    # par = zeros(Float64,size(lb))
    # if it == 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 300.0)) #20 min
    # elseif it > 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 1200.0)) #20 min
    # end 
    
    # mixture parameter
    betas = theta[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
    Ez1 = T[:,1,:]
    Ez0 = T[:,2,:]
    
    pairs = size(pairsindex,2)
    betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)

    
    # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
   
    # println(par.minimizer,pg[1])
    if !(Optim.converged(par)) && isnan(betas[1]) 
        println(betas[1])
        if start == true
            return zeros(Float64,size(theta,1))
        elseif start == false
            return Optim.converged(par)
        end
    end
    return vcat(par.minimizer,betasnew)#,pg[1])
end



function pmllastinv(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])

    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )
    Psi = diagm(x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    kappa = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)]    #interc2 = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p)]
    # prop = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]

    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
    
    mu = L2*kappa
    
    # Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig1 = L1*Phi*L1' + diagm(ones(p))
    Sig2=L2*L2'+Psi



    plimix=zeros(Float64,n);
    tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]

    
    for idx=1:size(pairsindex,2) #parallelizzabile
        # CFA
        Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
        Rcfa=cov2cor(Scfa);
        Bcfa=hcat(tfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
        tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
        Tt1 = T[:,1,idx]
        
        #EFA
        Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
        Refa=cov2cor(Sefa);
        Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
        tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
        Tt2 = T[:,2,idx]

        plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
    end
    f=-sum(plimix)
return f
end





function fginv!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
   
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )

    Psi = diagm(x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    kappa = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)]    #interc2 = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
    # prop = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]

    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
  
    mu = L2*kappa

    Sig1= L1*Phi*L1' + diagm(ones(p))
    Sig2= L2*L2'+Psi


    if G !== nothing

        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradkappacefa = zeros(Float64,K,size(pairsindex,2))
        gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)

        gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)


        gradlscfaidx = zeros(Float64,q,size(pairsindex,2),2)
        gradlambci = zeros(Float64,q,c+1,size(pairsindex,2),2)
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
            
        idx = 0
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1

                sig1s = hcat(vcat(Sig1[i,i],Sig1[i,j]),vcat(Sig1[j,i],Sig1[j,j]))
                cormat1 = cov2cor(sig1s)
                rho1 = cormat1[1,2]
                sig1sij = diag(sig1s)   
                sig1ij = sqrt.(sig1sij)
                sig1ij13 = sig1ij[1]^(3)
                sig1ij23 = sig1ij[2]^(3)
                L1i=L1[i,:]
                L1j=L1[j,:]

                sig2s = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sig2s)
                rho2 = cormat2[1,2]
                muij = mu[vcat(i,j)]
                sig2sij = diag(sig2s)
                sig2ij = sqrt.(sig2sij)
                sig2ij13 = sig2ij[1]^(3)
                sig2ij23 = sig2ij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]
                Psii=Psi[i,i]
                Psij=Psi[j,j]

                @inbounds for N in 1:n
                    xni = X[N,i]
                    xnj = X[N,j]
                    t1i = t1[i,xni]
                    t1i1 = t1[i,xni+1]
                    t1j = t1[j,xnj]
                    t1j1 = t1[j,xnj+1]
                    t1ist = (t1i)/sig1ij[1]
                    t1i1st = (t1i1)/sig1ij[1]
                    t1jst = (t1j)/sig1ij[2]
                    t1j1st = (t1j1)/sig1ij[2]
                    T1 = T[N,1,idx][1]
                    T2 = T[N,2,idx][1]
                    t2i = t1[i,xni]
                    t2i1 = t1[i,xni+1]
                    t2j = t1[j,xnj]
                    t2j1 = t1[j,xnj+1]
                    t2ist = (t2i-muij[1])/sig2ij[1]
                    t2i1st = (t2i1-muij[1])/sig2ij[1]
                    t2jst = (t2j-muij[2])/sig2ij[2]
                    t2j1st = (t2j1-muij[2])/sig2ij[2]
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                
                    pijicfa = bvnP(t1ist, t1i1st, t1jst, t1j1st, rho1, false)
                    derivNorm1i1 = (T1/pijicfa)*(pdf(Normal(0,1),t1i1st)*(cdf(Normal(0,1),(t1j1st-rho1*t1i1st)/sqrtrij1)-cdf(Normal(0,1),(t1jst-rho1*t1i1st)/sqrtrij1)))
                    derivNorm1i = (T1/pijicfa)*(pdf(Normal(0,1),t1ist)*(cdf(Normal(0,1),(t1jst-rho1*t1ist)/sqrtrij1)-cdf(Normal(0,1),(t1j1st-rho1*t1ist)/sqrtrij1)))
                    derivNorm1j1 = (T1/pijicfa)*(pdf(Normal(0,1),t1j1st)*(cdf(Normal(0,1),(t1i1st-rho1*t1j1st)/sqrtrij1)-cdf(Normal(0,1),(t1ist-rho1*t1j1st)/sqrtrij1)))
                    derivNorm1j = (T1/pijicfa)*(pdf(Normal(0,1),t1jst)*(cdf(Normal(0,1),(t1ist-rho1*t1jst)/sqrtrij1)-cdf(Normal(0,1),(t1i1st-rho1*t1jst)/sqrtrij1)))
                
                    gradtaucfa[xni+1,idx,1] +=  derivNorm1i1/(sig1ij[1])
                    gradtaucfa[xni,idx,1] += derivNorm1i/(sig1ij[1])
                    gradtaucfa[xnj+1,idx,2] += derivNorm1j1/(sig1ij[2])
                    gradtaucfa[xnj,idx,2] += derivNorm1j/(sig1ij[2])


                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))

                    if isnan(derivNormi1)
                        derivNormi1 = 0.
                    end
                    if isnan(derivNormi)
                        derivNormi = 0.
                    end
                    if isnan(derivNormj1)
                        derivNormj1 = 0.
                    end
                    if isnan(derivNormj)
                        derivNormj = 0.
                    end

                    
                    if !(t1ist == -Inf && t1jst == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1jst))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else                  
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1st == Inf && t1j1st == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1j1st))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end

                    if !(t1i1st == Inf && t1jst == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1jst))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end

                    if !(t1ist == -Inf && t1j1st == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1j1st))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    dpidrij1 = (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)#*vcat((-((L1j'*Phi.*t1j1st)./(2*sig1sij[2]))),(-((L1i'*Phi.*t1i1st)./(2*sig1sij[1]))))#vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)#
                

                    gradlscfaidx[1:q,idx,1] .+= dpidrij1*(Phi*L1j)/(sig1ij13*sig1ij[2])#(-Phi*L1i*(L1i'*L1j)/(sig1ij13*sig1ij[2])+Phi*L1j/((sig1ij[1]*sig1ij[2])))
                    gradlscfaidx[1:q,idx,2] .+= dpidrij1*(Phi*L1i)/(sig1ij23*sig1ij[1])#(-Phi*L1j*(L1j'*L1i)/(sig1ij23*sig1ij[1])+Phi*L1i/((sig1ij[2]*sig1ij[1])))
                    
                    gradlambci[1:q,xni+1,idx,1] .+=  derivNorm1i1 *(-((L1i*t1i1st)./(sig1sij[1])))#(-((L1i*t1i1).+(kappak*Psii))./((sig1ij13)))
                    gradlambci[1:q,xnj+1,idx,2] .+=  derivNorm1j1 *(-((L1j*t1j1st)./(sig1sij[2])))#(-((L1j*t1j1).+(kappak*Psij))./((sig1ij23)))
                    gradlambci[1:q,xni,idx,1]   .+=  derivNorm1i  *(-((L1i*t1ist)./(sig1sij[1])))#(-((L1i*t1i).+(kappak*Psii))./((sig1ij13)))
                    gradlambci[1:q,xnj,idx,2]   .+=  derivNorm1j  *(-((L1j*t1jst)./(sig1sij[2])))#(-((L1j*t1j).+(kappak*Psij))./((sig1ij23)))               
                    
                    gradphici[1:q,1:q,idx] .+= dpidrij1*(L1i*L1j')/(sig1ij[1]*sig1ij[2])

                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sig2ij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sig2ij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sig2ij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sig2ij[2])


                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)

                    # gradkappacefa[xni+1,1:K,idx,1] += (derivNormi1*(-L2[i,:]/sig2ij[1]))
                    # gradkappacefa[xnj+1,1:K,idx,2] += (derivNormj1*(-L2[j,:]/sig2ij[2]))
                    # gradkappacefa[xni,1:K,idx,1] += (derivNormi*(-L2[i,:]/sig2ij[1]))
                    # gradkappacefa[xnj,1:K,idx,2] += (derivNormj*(-L2[j,:]/sig2ij[2]))
                    gradkappacefa[1:K,idx] .+= derivNormi1*(-L2[i,:]/sig2ij[1]) + derivNormi*(-L2[i,:]/sig2ij[1]) + derivNormj1*(-L2[j,:]/sig2ij[2]) + derivNormj*(-L2[j,:]/sig2ij[2])
                    
                    gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sig2ij13)*(sig2ij[2])))
                    gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sig2ij23)*(sig2ij[1])))

                    gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij[1]-t2i1))/(2*(sig2ij13))
                    gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij[2]-t2j1))/(2*(sig2ij23))
                    gradpsiefa[xni,idx,1] += (derivNormi*(muij[1]-t2i))/(2*(sig2ij13))
                    gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij[2]-t2j))/(2*(sig2ij23))
                    
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        kappak=kappa[k]
                        gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sig2ij13*sig2ij[2])+L2jk/((sig2ij[1]*sig2ij[2])))#((L2j*Psi[i,i])/(sig2ij13*sig2ij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sig2ij23*sig2ij[1])+L2ik/((sig2ij[2]*sig2ij[1])))#((L2[i,:].*Psi[j,j])./(sig2ij23*sig2ij[1]))

                        gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(-((L2ik*t2i1st)./(sig2sij[1]).+(kappak)./(sig2ij[1])))#(-((L2ik*t2i1st).+(kappak*Psii))./((sig2ij13)))
                        gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(-((L2jk*t2j1st)./(sig2sij[2]).+(kappak)./(sig2ij[2])))#(-((L2jk*t2j1st).+(kappak*Psij))./((sig2ij23)))
                        gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(-((L2ik*t2ist)./(sig2sij[1]).+(kappak)./(sig2ij[1])))#(-((L2ik*t2ist).+(kappak*Psii))./((sig2ij13)))
                        gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(-((L2jk*t2jst)./(sig2sij[2]).+(kappak)./(sig2ij[2])))#(-((L2jk*t2jst).+(kappak*Psij))./((sig2ij23)))
                    end 
                
                end
            end
        end

        vecpair = vec(pairsindex)
        liktaushared = zeros(p*(c-1));idpcfa = 0   
        @inbounds for j in 2:c
            @inbounds for i in 1:p
                idpcfa += 1
                liktaushared[idpcfa] = -sum(vec(vcat(gradtaucfa[j,:,1]',gradtaucfa[j,:,2]'))[vecpair.==i])-sum(vec(vcat(gradtauefa[j,:,1]',gradtauefa[j,:,2]'))[vecpair.==i])
            end
        end

        



        likls1 = zeros(p,q); llamb1 = sum(gradlambci[:,2:c,:,:],dims=2)[:,1,:,:]; llam1 = zeros(p,q); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
        likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
        llam1idx  = zeros(p,q); likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        @inbounds for i in 1:p
            @inbounds for qs in 1:q
                llam1[i,qs] = -sum(vec(llamb1[qs,:,:]')[vecpair.==i])
                llam1idx[i,qs] = -sum(vec(gradlscfaidx[qs,:,:]')[vecpair.==i])
            end
            likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
            likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end

        likpsitrue = likpsi.+likpsiidx
        liklstrue = likls.+liklsidx
        liklamtrue = llam1.+llam1idx

        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]

        likappa = -sum(gradkappacefa,dims=2)

        G .= vec(vcat(liktaushared,sum(liklamtrue.*L_str,dims=2),likphida,liklstrue[indices],likpsitrue,likappa))#likpsitrue[2:end]

    end
    if F !== nothing 
  
    
        plimix=zeros(Float64,n);    
        tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
        tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
    
       
       
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
            Rcfa=cov2cor(Scfa);
            Bcfa=hcat(tfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            Tt1 = T[:,1,idx]
            
            #EFA
            Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            Tt2 = T[:,2,idx]
            plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2;
            
        end
        f=-sum(plimix)
        # @info "Current objective: $f"
        return f
    end
end



















# Invariatìnt con delta e theta parameterization
# ########################################################################### Latent means EFA + invariant thresholds + Psi estimated ##################################################################################################

# function startvalsinv(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
#     liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)+ncov,replix);#p-1
#     for i in 1:replix
#         println(i)

#         tau1 = vec(hcat(rand(Uniform(-2.3,-1.),p).+reshape(repeat(collect(range(0,3.5,c-1)),inner=p),p,(c-1))))
#         lambda1 = rand(Uniform(0.5,.90),p);
#         phi=rand(LKJ(q,1))
#         phi=phi[tril!(trues(size(phi)), -1)]
#         lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(0.15,.90),Int((p*K-(K*(K-1)/2))))))
#         lambda2=lamb2[lamb2.!=0.]
#         psi=rand(Uniform(0.2,1.5),p)
#         kappa=rand(Uniform(-0.2,2),K)



#         prop = rand(Uniform(0.65,0.95))
#         bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.1,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


#         param_initials = vcat(tau1,lambda1,phi,lambda2,psi,kappa,bets)

#         # first iteration
#         plimix = E_stepMixinv(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
#         T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
#         t = time()
#         params[:,i] = M_parsMixinv(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
#         dt = time() - t
#         prrot=params[:,i]
#         liksmix[i] = E_stepMixinv(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
#     end
#     param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
#     lik=liksmix[findall(liksmix.==maximum(liksmix))]

#     tau1=param_initials0[(1):((c-1)*p)]
#     tau1=vec(sort(reshape(tau1,p,c-1),dims=2))
#     lambda1=param_initials0[((c-1)*p+1):((c-1)*p+p)]
#     phi=param_initials0[((c-1)*p+p+1):((c-1)*p+p+Int(q*(q-1)/2))]
#     lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
#     lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
#     phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
#     lambda2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2))))]
#     psi=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p))]# psi=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1))]
#     kappa=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+K)]#interc=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
#     lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
#     lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
#     psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
#     bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+K+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+K+ncov)]#    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


#     param_initialii = vcat(tau1,lambda1,phi,lambda2,psi,kappa,bets)

#     println(param_initialii,liksmix,"first iteration: done")
#     return param_initialii
# end


# function EMparsMixinv(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
#     local plimix::Array, T,oldlik
#     n = size(X,1)
#     er=0;
#     # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
#     it=0;
#     dev=Inf; dev1=Inf; dev2=Inf;

    
#     oldlik=-1000000.;
#     oldlik2=1000000.;
#     oldlikA=-1000000.;
#     lik=0.
    

#     noerror = true; devnostop = true;  sum_devs = 0;
#     while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
#         it = it+1;
#         param=param01;

#         T00 = 1
#         plimix = E_stepMixinv(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
#         T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
#         t = time()
#         param01 = M_parsMixinv(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
#         dt = time() - t
#         if param01 == false
#             return return(vcat(param01,vec(plimix),vec(T),0.,0.))
#         end
#         lik = E_stepMixinv(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);

#         #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
#         con = (lik-oldlik)/(oldlik-oldlik2);
#         likA = oldlik + (lik-oldlik)/(1-con)
#         old_dev = dev1
#         dev1 = abs(likA-oldlikA)
#         dev2 = abs(lik-oldlik)
#         println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
#         println(param01)
#         dev = abs(lik-oldlik)/abs(oldlik)

#         noerror = oldlik<lik

#         oldlikA=likA
#         oldlik2=oldlik
#         oldlik=lik;

#         # sum_devs += dev > old_dev
#         # devstop = sum_devs <= 3
#     end
    
#     # classificazioni dopo ultima massimizzazione
#     plimix = E_stepMixinv(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
#     T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

#     if noerror == true && it<=maxit 
#         return vcat(param01,vec(plimix),vec(T),lik,1.)
#     else
#         return vcat(param01,vec(plimix),vec(T),0.,0.)
#     end
# end





# function E_stepMixinv(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
#     if sum(param .== 0.) == size(param,1)
#         return -Inf
#     end
    
#     trsi = (c-1)
#     t1s = param[1:trsi*p]
#     L1 = param[(trsi*p+1):(trsi*p+p)].*L_str 
#     Phi = reconstruct_correlation_matrix(param[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    
#     L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )
#     Psi = diagm(param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p-1))]))#K
#     kappa = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)]   # interc2 = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
#     betas = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+ncov)]  #  betas = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

#     propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
#     pg = hcat(propo,1. .-propo)

#     mu = L2*kappa
    
#     Sigms=ones(p,p,2);
#     Sigms[:,:,1]=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))#L1*Phi*L1' + diagm(ones(p))
#     Sigms[:,:,2]=L2*L2'+Psi 
    
#     br1 = reshape(t1s,p,trsi)
#     t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
   

#     tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
   
    
#     if estep == true
#         low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#         plimix=ones(n,2,size(pairsindex,2));
    
#         for idx=1:size(pairsindex,2) #parallelizzabile
#             for g=1:2
#                 if g==1
#                     s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                     r2c=cov2cor(s2c);
#                     up=(tfa1[:,pairsindex[:,idx]])#./sqrt.(diag(s2c)'))
#                     low=(tfa[:,pairsindex[:,idx]])#./sqrt.(diag(s2c)'))
#                 elseif g==2
#                     s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                     r2c=cov2cor(s2c);
#                     up=(tfa1[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                     low=(tfa[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 end
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix[:,g,idx]= t3.*pg[:,g];
#             end
#         end
#     # % probabilities in each row must sum to 1 (application of Bayes law)
#         return plimix
#     elseif estep == false
#         low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#         plimix=ones(n,2,size(pairsindex,2));

#         for idx=1:size(pairsindex,2) #parallelizzabile
#             for g=1:2
#                 if g==1
#                     s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                     r2c=cov2cor(s2c);
#                     up=(tfa1[:,pairsindex[:,idx]])#./sqrt.(diag(s2c)'))
#                     low=(tfa[:,pairsindex[:,idx]])#./sqrt.(diag(s2c)'))
#                 elseif g==2
#                     s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                     r2c=cov2cor(s2c);
#                     up=(tfa1[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                     low=(tfa[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 end
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix[:,g,idx] = t3
#             end
#         end
#         lik1=(log.(plimix));
#         lik1=sum(lik1,dims=3);
#         U=exp.(lik1)[:,:,1].*pg; 
#         cllik=sum(log.(sum(U,dims=2)))
#         return cllik
#     end
# end
    

    


# function M_parsMixinv(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)

#     trsi = (c-1)*p #attenzione p è già qui
   
#     # PML - CFA & EFA
#     theta_init = theta[1:(size(theta,1)-ncov)]#0.000001
#     # # -1 loadings
#     # lb=vec(vcat(fill(-5.1,trsi),fill(0.00001,p),fill(-1.,Int((q*(q-1))/2)),fill(1.,trs),fill(-1.,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-10,p)));
#     # # 0. loadings
#     lb=vec(vcat(fill(-5.1,trsi),fill(-1,p),fill(-1.,Int((q*(q-1))/2)),fill(-Inf,Int(p*K-K*(K-1)/2)),fill(1.0e-36,p),fill(-5,K)));#1.0e-5
#     ub=vec(vcat(fill(5.1,trsi),fill(1.,p),fill(1.,Int((q*(q-1))/2)),fill(Inf,Int(p*K-K*(K-1)/2)),fill(Inf,p),fill(10.,K)));

#     par = zeros(Float64,size(lb))
#     if it == 1
#         par = optimize(Optim.only_fg!((F,G,x)->fginv!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min ,Optim.Options(iterations=1000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0))
#     elseif it > 1
#         par = optimize(Optim.only_fg!((F,G,x)->fginv!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min
#     end 


#     # par = optimize(x->pmllast(x,X,T,p,q,K,c,n,pairsindex,L_str,indices),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

#     # par = zeros(Float64,size(lb))
#     # if it == 1
#     #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 300.0)) #20 min
#     # elseif it > 1
#     #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 1200.0)) #20 min
#     # end 
    
#     # mixture parameter
#     betas = theta[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
#     Ez1 = T[:,1,:]
#     Ez0 = T[:,2,:]
    
#     pairs = size(pairsindex,2)
#     betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)

    
#     # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
   
#     # println(par.minimizer,pg[1])
#     if !(Optim.converged(par)) && isnan(betas[1]) 
#         println(betas[1])
#         if start == true
#             return zeros(Float64,size(theta,1))
#         elseif start == false
#             return Optim.converged(par)
#         end
#     end
#     return vcat(par.minimizer,betasnew)#,pg[1])
# end



# function pmllastinv(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
#     trsi = (c-1)
#     t1s = x[1:trsi*p]
#     L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
#     Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])

#     L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )
#     Psi = diagm(x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p-1))]))#K
#     kappa = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)]    #interc2 = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p)]
#     # prop = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]

#     br1 = reshape(t1s,p,trsi)
#     t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
    
#     mu = L2*kappa
    
#     Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
#     # Sig1 = L1*Phi*L1' + diagm(ones(p))
#     Sig2=L2*L2'+Psi



#     plimix=zeros(Float64,n);
#     tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]

    
#     for idx=1:size(pairsindex,2) #parallelizzabile
#         # CFA
#         Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
#         Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
#         # Bcfa=hcat(tfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
#         tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
#         Tt1 = T[:,1,idx]
        
#         #EFA
#         Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
#         Refa=cov2cor(Sefa);
#         Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
#         tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
#         Tt2 = T[:,2,idx]

#         plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
#     end
#     f=-sum(plimix)
# return f
# end

# function fginv!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
#     # do common computations here
#     # ...
#     trsi = (c-1)
#     t1s = x[1:trsi*p]
#     L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
#     Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
   
#     L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )

#     Psi = diagm(x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p-1))]))#K
#     kappa = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)]    #interc2 = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
#     # prop = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]

#     br1 = reshape(t1s,p,trsi)
#     t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
  
#     mu = L2*kappa

#     Sig1= L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))#L1*Phi*L1' + diagm(ones(p))
#     Sig2= L2*L2'+Psi


#     if G !== nothing

#         # gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
#         # gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
#         # gradkappacefa = zeros(Float64,K,size(pairsindex,2))
#         # gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)

#         # gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
#         # gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)


#         # gradlscfaidx = zeros(Float64,q,size(pairsindex,2),2)
#         # gradlambci = zeros(Float64,q,c+1,size(pairsindex,2),2)
#         # gradphici = zeros(Float64,q,q,size(pairsindex,2))
#         # gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
            
#         # idx = 0
#         # @inbounds for i in 1:(p-1)
#         #     @inbounds for j in (i+1):p
#         #         idx += 1

#         #         sig1s = hcat(vcat(Sig1[i,i],Sig1[i,j]),vcat(Sig1[j,i],Sig1[j,j]))
#         #         cormat1 = cov2cor(sig1s)
#         #         rho1 = cormat1[1,2]
#         #         sig1sij = diag(sig1s)   
#         #         sig1ij = sqrt.(sig1sij)
#         #         sig1ij13 = sig1ij[1]^(3)
#         #         sig1ij23 = sig1ij[2]^(3)
#         #         L1i=L1[i,:]
#         #         L1j=L1[j,:]

#         #         sig2s = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
#         #         cormat2 = cov2cor(sig2s)
#         #         rho2 = cormat2[1,2]
#         #         muij = mu[vcat(i,j)]
#         #         sig2sij = diag(sig2s)
#         #         sig2ij = sqrt.(sig2sij)
#         #         sig2ij13 = sig2ij[1]^(3)
#         #         sig2ij23 = sig2ij[2]^(3)
#         #         L2i=L2[i,:]
#         #         L2j=L2[j,:]
#         #         Psii=Psi[i,i]
#         #         Psij=Psi[j,j]

#         #         @inbounds for N in 1:n
#         #             xni = X[N,i]
#         #             xnj = X[N,j]
#         #             t1i = t1[i,xni]
#         #             t1i1 = t1[i,xni+1]
#         #             t1j = t1[j,xnj]
#         #             t1j1 = t1[j,xnj+1]
#         #             t1ist = (t1i)/sig1ij[1]
#         #             t1i1st = (t1i1)/sig1ij[1]
#         #             t1jst = (t1j)/sig1ij[2]
#         #             t1j1st = (t1j1)/sig1ij[2]
#         #             T1 = T[N,1,idx][1]
#         #             T2 = T[N,2,idx][1]
#         #             t2i = t1[i,xni]
#         #             t2i1 = t1[i,xni+1]
#         #             t2j = t1[j,xnj]
#         #             t2j1 = t1[j,xnj+1]
#         #             t2ist = (t2i-muij[1])/sig2ij[1]
#         #             t2i1st = (t2i1-muij[1])/sig2ij[1]
#         #             t2jst = (t2j-muij[2])/sig2ij[2]
#         #             t2j1st = (t2j1-muij[2])/sig2ij[2]
#         #             sqrtrij1 = sqrt(1-rho1^2)
#         #             sqrtrij2 = sqrt(1-rho2^2)
                                
#         #             pijicfa = bvnP(t1ist, t1i1st, t1jst, t1j1st, rho1, false)
#         #             derivNorm1i1 = (T1/pijicfa)*(pdf(Normal(0,1),t1i1st)*(cdf(Normal(0,1),(t1j1st-rho1*t1i1st)/sqrtrij1)-cdf(Normal(0,1),(t1jst-rho1*t1i1st)/sqrtrij1)))
#         #             derivNorm1i = (T1/pijicfa)*(pdf(Normal(0,1),t1ist)*(cdf(Normal(0,1),(t1jst-rho1*t1ist)/sqrtrij1)-cdf(Normal(0,1),(t1j1st-rho1*t1ist)/sqrtrij1)))
#         #             derivNorm1j1 = (T1/pijicfa)*(pdf(Normal(0,1),t1j1st)*(cdf(Normal(0,1),(t1i1st-rho1*t1j1st)/sqrtrij1)-cdf(Normal(0,1),(t1ist-rho1*t1j1st)/sqrtrij1)))
#         #             derivNorm1j = (T1/pijicfa)*(pdf(Normal(0,1),t1jst)*(cdf(Normal(0,1),(t1ist-rho1*t1jst)/sqrtrij1)-cdf(Normal(0,1),(t1i1st-rho1*t1jst)/sqrtrij1)))
                
#         #             pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
#         #             derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
#         #             derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
#         #             derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
#         #             derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))

#         #             if isnan(derivNormi1)
#         #                 derivNormi1 = 0.
#         #             end
#         #             if isnan(derivNormi)
#         #                 derivNormi = 0.
#         #             end
#         #             if isnan(derivNormj1)
#         #                 derivNormj1 = 0.
#         #             end
#         #             if isnan(derivNormj)
#         #                 derivNormj = 0.
#         #             end

                    
#         #             if !(t1ist == -Inf && t1jst == -Inf) 
#         #                 xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1jst))
#         #                 xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
#         #             else                  
#         #                 xnixnjc = 0.0
#         #                 xnixnje = 0.0
#         #             end
                        
#         #             if !(t1i1st == Inf && t1j1st == Inf) 
#         #                 xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1j1st))
#         #                 xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
#         #             else
#         #                 xni1xnj1c = 0.0
#         #                 xni1xnj1e = 0.0
#         #             end

#         #             if !(t1i1st == Inf && t1jst == -Inf) 
#         #                 xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1jst))
#         #                 xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
#         #             else
#         #                 xni1xnjc = 0.0
#         #                 xni1xnje = 0.0
#         #             end

#         #             if !(t1ist == -Inf && t1j1st == Inf) 
#         #                 xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1j1st))
#         #                 xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
#         #             else
#         #                 xnixnj1c = 0.0
#         #                 xnixnj1e = 0.0
#         #             end
#         #             dpidrij1 = (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)#*vcat((-((L1j'*Phi.*t1j1st)./(2*sig1sij[2]))),(-((L1i'*Phi.*t1i1st)./(2*sig1sij[1]))))#vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)#
                

#         #             gradlscfaidx[1:q,idx,1] .+= dpidrij1*(Phi*L1j)/(sig1ij13*sig1ij[2])#(-Phi*L1i*(L1i'*L1j)/(sig1ij13*sig1ij[2])+Phi*L1j/((sig1ij[1]*sig1ij[2])))
#         #             gradlscfaidx[1:q,idx,2] .+= dpidrij1*(Phi*L1i)/(sig1ij23*sig1ij[1])#(-Phi*L1j*(L1j'*L1i)/(sig1ij23*sig1ij[1])+Phi*L1i/((sig1ij[2]*sig1ij[1])))
                    
#         #             gradlambci[1:q,xni+1,idx,1] .+=  derivNorm1i1 *(-((L1i*t1i1st)./(sig1sij[1])))#(-((L1i*t1i1).+(kappak*Psii))./((sig1ij13)))
#         #             gradlambci[1:q,xnj+1,idx,2] .+=  derivNorm1j1 *(-((L1j*t1j1st)./(sig1sij[2])))#(-((L1j*t1j1).+(kappak*Psij))./((sig1ij23)))
#         #             gradlambci[1:q,xni,idx,1]   .+=  derivNorm1i  *(-((L1i*t1ist)./(sig1sij[1])))#(-((L1i*t1i).+(kappak*Psii))./((sig1ij13)))
#         #             gradlambci[1:q,xnj,idx,2]   .+=  derivNorm1j  *(-((L1j*t1jst)./(sig1sij[2])))#(-((L1j*t1j).+(kappak*Psij))./((sig1ij23)))               
                    
#         #             gradphici[1:q,1:q,idx] .+= dpidrij1*(L1i*L1j')#/(sig1ij[1]*sig1ij[2])

#         #             gradtaucfa[xni+1,idx,1] +=  derivNorm1i1/(sig1ij[1])
#         #             gradtaucfa[xni,idx,1] += derivNorm1i/(sig1ij[1])
#         #             gradtaucfa[xnj+1,idx,2] += derivNorm1j1/(sig1ij[2])
#         #             gradtaucfa[xnj,idx,2] += derivNorm1j/(sig1ij[2])
                    
#         #             gradtauefa[xni+1,idx,1] += (derivNormi1)/(sig2ij[1])
#         #             gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sig2ij[2])
#         #             gradtauefa[xni,idx,1] += (derivNormi)/(sig2ij[1])
#         #             gradtauefa[xnj,idx,2] += (derivNormj)/(sig2ij[2])


#         #             dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)

#         #             # gradkappacefa[xni+1,1:K,idx,1] += (derivNormi1*(-L2[i,:]/sig2ij[1]))
#         #             # gradkappacefa[xnj+1,1:K,idx,2] += (derivNormj1*(-L2[j,:]/sig2ij[2]))
#         #             # gradkappacefa[xni,1:K,idx,1] += (derivNormi*(-L2[i,:]/sig2ij[1]))
#         #             # gradkappacefa[xnj,1:K,idx,2] += (derivNormj*(-L2[j,:]/sig2ij[2]))
#         #             gradkappacefa[1:K,idx] .+= derivNormi1*(-L2[i,:]/sig2ij[1]) + derivNormi*(-L2[i,:]/sig2ij[1]) + derivNormj1*(-L2[j,:]/sig2ij[2]) + derivNormj*(-L2[j,:]/sig2ij[2])
                    
#         #             gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sig2ij13)*(sig2ij[2])))
#         #             gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sig2ij23)*(sig2ij[1])))

#         #             gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij[1]-t2i1))/(2*(sig2ij13))
#         #             gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij[2]-t2j1))/(2*(sig2ij23))
#         #             gradpsiefa[xni,idx,1] += (derivNormi*(muij[1]-t2i))/(2*(sig2ij13))
#         #             gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij[2]-t2j))/(2*(sig2ij23))
                    
#         #             @inbounds for k in 1:K
#         #                 L2ik=L2[i,k]
#         #                 L2jk=L2[j,k]
#         #                 kappak=kappa[k]
#         #                 gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sig2ij13*sig2ij[2])+L2jk/((sig2ij[1]*sig2ij[2])))#((L2j*Psi[i,i])/(sig2ij13*sig2ij[2]))
#         #                 gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sig2ij23*sig2ij[1])+L2ik/((sig2ij[2]*sig2ij[1])))#((L2[i,:].*Psi[j,j])./(sig2ij23*sig2ij[1]))

#         #                 gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(-((L2ik*t2i1st)./(sig2sij[1]).+(kappak)./(sig2ij[1])))#(-((L2ik*t2i1st).+(kappak*Psii))./((sig2ij13)))
#         #                 gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(-((L2jk*t2j1st)./(sig2sij[2]).+(kappak)./(sig2ij[2])))#(-((L2jk*t2j1st).+(kappak*Psij))./((sig2ij23)))
#         #                 gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(-((L2ik*t2ist)./(sig2sij[1]).+(kappak)./(sig2ij[1])))#(-((L2ik*t2ist).+(kappak*Psii))./((sig2ij13)))
#         #                 gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(-((L2jk*t2jst)./(sig2sij[2]).+(kappak)./(sig2ij[2])))#(-((L2jk*t2jst).+(kappak*Psij))./((sig2ij23)))
#         #             end 
                
#         #         end
#         #     end
#         # end


#         # vecpair = vec(pairsindex)
#         # liktaushared = zeros(p*(c-1));idpcfa = 0   
#         # @inbounds for j in 2:c
#         #     @inbounds for i in 1:p
#         #         idpcfa += 1
#         #         liktaushared[idpcfa] = -sum(vec(vcat(gradtaucfa[j,:,1]',gradtaucfa[j,:,2]'))[vecpair.==i])-sum(vec(vcat(gradtauefa[j,:,1]',gradtauefa[j,:,2]'))[vecpair.==i])
#         #     end
#         # end


#         # likls1 = zeros(p,q); llamb1 = sum(gradlambci[:,2:c,:,:],dims=2)[:,1,:,:]; llam1 = zeros(p,q); ids = find_ones_columns(L_str);
#         # likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
#         # likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
#         # llam1idx  = zeros(p,q); likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
#         # @inbounds for i in 1:p
#         #     @inbounds for qs in 1:q
#         #         llam1[i,qs] = -sum(vec(llamb1[qs,:,:]')[vecpair.==i])
#         #         llam1idx[i,qs] = -sum(vec(gradlscfaidx[qs,:,:]')[vecpair.==i])
#         #     end
#         #     likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
#         #     likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
#         #     @inbounds for k in 1:K
#         #         likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
#         #         liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
#         #     end
#         # end

    

#         # likpsitrue = likpsi.+likpsiidx
#         # liklstrue = likls.+liklsidx
#         # liklamtrue = llam1.+llam1idx

#         # likp = sum(gradphici,dims=(3))[:,:,1]
#         # likphida = -likp[triu!(trues(size(likp)), 1)]

#         # likappa = -sum(gradkappacefa,dims=2)

#         # G .= vec(vcat(liktaushared,sum(liklamtrue.*L_str,dims=2),likphida,liklstrue[indices],likpsitrue,likappa))#likpsitrue[2:end]


        # gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        # gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        # gradkappacefa = zeros(Float64,K,size(pairsindex,2))
        # gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)
    
        # gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
        # gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)
    
        # gradlambci = zeros(Float64,2,q,size(pairsindex,2))
        # gradphici = zeros(Float64,q,q,size(pairsindex,2))
        # gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
        
        # idx = 0
        # @inbounds for i in 1:(p-1)
        #     @inbounds for j in (i+1):p
        #         idx += 1
        #         rho1 = Sig1[i,j]
        #         cormat1 = ones(Float64,2,2); cormat1[1,2]=rho1; cormat1[2,1]=rho1
        #         sigs = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
        #         cormat2 = cov2cor(sigs)
        #         rho2 = cormat2[1,2]
        #         muij = mu[vcat(i,j)]
        #         sigsij = diag(sigs)
        #         sigij = sqrt.(sigsij)
        #         sigij13 = sigij[1]^(3)
        #         sigij23 = sigij[2]^(3)
        #         L2i=L2[i,:]
        #         L2j=L2[j,:]
        #         Psii=Psi[i,i]
        #         Psij=Psi[j,j]
        #         @inbounds for N in 1:n
        #             xni = X[N,i]
        #             xnj = X[N,j]
        #             t1i = t1[i,xni]
        #             t1i1 = t1[i,xni+1]
        #             t1j = t1[j,xnj]
        #             t1j1 = t1[j,xnj+1]
        #             T1 = T[N,1,idx][1]
        #             T2 = T[N,2,idx][1]
        #             t2i = t1[i,xni]
        #             t2i1 = t1[i,xni+1]
        #             t2j = t1[j,xnj]
        #             t2j1 = t1[j,xnj+1]
        #             t2ist = (t2i-muij[1])/sigij[1]
        #             t2i1st = (t2i1-muij[1])/sigij[1]
        #             t2jst = (t2j-muij[2])/sigij[2]
        #             t2j1st = (t2j1-muij[2])/sigij[2]
        #             sqrtrij1 = sqrt(1-rho1^2)
        #             sqrtrij2 = sqrt(1-rho2^2)
                                   
        #             pijicfa = bvnP(t1i, t1i1, t1j, t1j1, rho1, false)
        #             gradtaucfa[xni+1,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i1)*(cdf(Normal(0,1),(t1j1-rho1*t1i1)/sqrtrij1)-cdf(Normal(0,1),(t1j-rho1*t1i1)/sqrtrij1))) 
        #             gradtaucfa[xni,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i)*(cdf(Normal(0,1),(t1j-rho1*t1i)/sqrtrij1)-cdf(Normal(0,1),(t1j1-rho1*t1i)/sqrtrij1))) 
        #             gradtaucfa[xnj+1,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j1)*(cdf(Normal(0,1),(t1i1-rho1*t1j1)/sqrtrij1)-cdf(Normal(0,1),(t1i-rho1*t1j1)/sqrtrij1))) 
        #             gradtaucfa[xnj,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j)*(cdf(Normal(0,1),(t1i-rho1*t1j)/sqrtrij1)-cdf(Normal(0,1),(t1i1-rho1*t1j)/sqrtrij1))) 
                    
    
        #             pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
        #             derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
        #             derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
        #             derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
        #             derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))
    
        #             if isnan(derivNormi1)
        #                 derivNormi1 = 0.
        #             end
        #             if isnan(derivNormi)
        #                 derivNormi = 0.
        #             end
        #             if isnan(derivNormj1)
        #                 derivNormj1 = 0.
        #             end
        #             if isnan(derivNormj)
        #                 derivNormj = 0.
        #             end
    
                    
        #             if !(t1i == -Inf && t1j == -Inf) 
        #                 xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j))
        #                 xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
        #             else                  
        #                 xnixnjc = 0.0
        #                 xnixnje = 0.0
        #             end
                        
        #             if !(t1i1 == Inf && t1j1 == Inf) 
        #                 xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j1))
        #                 xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
        #             else
        #                 xni1xnj1c = 0.0
        #                 xni1xnj1e = 0.0
        #             end
        
        #             if !(t1i1 == Inf && t1j == -Inf) 
        #                 xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j))
        #                 xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
        #             else
        #                 xni1xnjc = 0.0
        #                 xni1xnje = 0.0
        #             end
        
        #             if !(t1i == -Inf && t1j1 == Inf) 
        #                 xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j1))
        #                 xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
        #             else
        #                 xnixnj1c = 0.0
        #                 xnixnj1e = 0.0
        #             end
        #             gradlambci[1:2,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)
        #             gradphici[1:q,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*(L1[i,:]*L1[j,:]')
                    
        #             #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
        #             gradtauefa[xni+1,idx,1] += (derivNormi1)/(sigij[1])
        #             gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sigij[2])
        #             gradtauefa[xni,idx,1] += (derivNormi)/(sigij[1])
        #             gradtauefa[xnj,idx,2] += (derivNormj)/(sigij[2])
    
    
        #             dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)

        #             gradkappacefa[1:K,idx] .+= derivNormi1*(-L2[i,:]/sigij[1]) + derivNormi*(-L2[i,:]/sigij[1]) + derivNormj1*(-L2[j,:]/sigij[2]) + derivNormj*(-L2[j,:]/sigij[2])
                    
        #             gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sigij13)*(sigij[2])))
        #             gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sigij23)*(sigij[1])))
    
        #             gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij[1]-t2i1))/(2*(sigij13))
        #             gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij[2]-t2j1))/(2*(sigij23))
        #             gradpsiefa[xni,idx,1] += (derivNormi*(muij[1]-t2i))/(2*(sigij13))
        #             gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij[2]-t2j))/(2*(sigij23))
                    
        #             @inbounds for k in 1:K
        #                 L2ik=L2[i,k]
        #                 L2jk=L2[j,k]
        #                 kappak=kappa[k]
        #                 gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sigij13*sigij[2])+L2jk/((sigij[1]*sigij[2])))#((L2j*Psi[i,i])/(sigij13*sigij[2]))
        #                 gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sigij23*sigij[1])+L2ik/((sigij[2]*sigij[1])))#((L2[i,:].*Psi[j,j])./(sigij23*sigij[1]))
    
        #                 gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(-((L2ik*t2i1st)./(sigsij[1]).+(kappak)./(sigij[1])))#(-((L2ik*t2i1).+(kappak*Psii))./((sigij13)))
        #                 gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(-((L2jk*t2j1st)./(sigsij[2]).+(kappak)./(sigij[2])))#(-((L2jk*t2j1).+(kappak*Psij))./((sigij23)))
        #                 gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(-((L2ik*t2ist)./(sigsij[1]).+(kappak)./(sigij[1])))#(-((L2ik*t2i).+(kappak*Psii))./((sigij13)))
        #                 gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(-((L2jk*t2jst)./(sigsij[2]).+(kappak)./(sigij[2])))#(-((L2jk*t2j).+(kappak*Psij))./((sigij23)))
        #             end 
                
        #         end
        #     end
        # end

        # vecpair = vec(pairsindex)
        # liktaushared = zeros(p*(c-1));idpcfa = 0   
        # @inbounds for j in 2:c
        #     @inbounds for i in 1:p
        #         idpcfa += 1
        #         liktaushared[idpcfa] = -sum(vec(vcat(gradtaucfa[j,:,1]',gradtaucfa[j,:,2]'))[vecpair.==i])-sum(vec(vcat(gradtauefa[j,:,1]',gradtauefa[j,:,2]'))[vecpair.==i])
        #     end
        # end
        
        
    
        
        # llamb1 = permutedims(gradlambci,(2,3,1)); llam1 = zeros(p); ids = find_ones_columns(L_str);
        # likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
        # likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
        # likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        # @inbounds for i in 1:p
        #     llam1[i] = -sum(vec(llamb1[ids[i],:,:]')[vecpair.==i])
        #     likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
        #     likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
        #     @inbounds for k in 1:K
        #         likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
        #         liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
        #     end
        # end
    
        # likpsitrue = likpsi.+likpsiidx
        # liklstrue = likls.+liklsidx
    
        # likp = sum(gradphici,dims=(3))[:,:,1]
        # likphida = -likp[triu!(trues(size(likp)), 1)]
    
        # likappa = -sum(gradkappacefa,dims=2)

#         G .= vec(vcat(liktaushared,llam1,likphida,liklstrue[indices],likpsitrue,likappa))#likpsitrue[2:end]
      
#     end
#     if F !== nothing 
  
    
#         plimix=zeros(Float64,n);    
#         tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
#         tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
    
       
       
#         for idx=1:size(pairsindex,2) #parallelizzabile
#             # CFA
#             Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
#             Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
#             # Bcfa=hcat(tfa[:,pairsindex[:,idx]]./sqrt.(diag(Sefa)'),tfa1[:,pairsindex[:,idx]]./sqrt.(diag(Sefa)'));
#             tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
#             Tt1 = T[:,1,idx]
            
#             #EFA
#             Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
#             Refa=cov2cor(Sefa);
#             Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
#             tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
#             Tt2 = T[:,2,idx]
#             plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2;
            
#         end
#         f=-sum(plimix)
#         # @info "Current objective: $f"
#         return f
#     end
# end







########################################################################### Latent means EFA + ONLY ONE SET OF invariant thresholds + Psi estimated ##################################################################################################


function startvalsthinv(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
    liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)+ncov,replix);#p-1
    for i in 1:replix
        println(i)

        tau1 = rand(Uniform(-1.9,-1.)).+repeat(collect(range(0,3.5,c-1)))
        lambda1 = rand(Uniform(0.5,.90),p);
        phi=rand(LKJ(q,1))
        phi=phi[tril!(trues(size(phi)), -1)]
        lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(-1.20,1.20),Int((p*K-(K*(K-1)/2))))))
        lambda2=lamb2[lamb2.!=0.]
        psi=rand(Uniform(0.2,2.5),p)
        kappa=rand(Uniform(-1.2,1.2),K)



        prop = rand(Uniform(0.65,0.97))
        bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.1,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


        param_initials = vcat(tau1,lambda1,phi,lambda2,psi,kappa,bets)

        # first iteration
        plimix = E_stepMixthinv(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        params[:,i] = M_parsMixthinv(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        dt = time() - t
        prrot=params[:,i]
        liksmix[i] = E_stepMixthinv(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
    end
    param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    lik=liksmix[findall(liksmix.==maximum(liksmix))]

    tau1=param_initials0[(1):((c-1))]
    tau1=vec(sort(reshape(tau1,1,c-1),dims=2))
    lambda1=param_initials0[((c-1)+1):((c-1)+p)]
    phi=param_initials0[((c-1)+p+1):((c-1)+p+Int(q*(q-1)/2))]
    # lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    # lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    # phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    lambda2=param_initials0[((c-1)+p+Int(q*(q-1)/2)+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2))))]
    psi=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p))]# psi=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1))]
    kappa=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+K)]#interc=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    # lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    # lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
    bets = param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+K+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+K+ncov)]#    bets = param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


    param_initialii = vcat(tau1,lambda1,phi,lambda2,psi,kappa,bets)

    println(param_initialii,liksmix,"first iteration: done")
    return param_initialii
end


function EMparsMixthinv(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
    local plimix::Array, T,oldlik
    n = size(X,1)
    er=0;
    # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
    it=0;
    dev=Inf; dev1=Inf; dev2=Inf;

    
    oldlik=-1000000.;
    oldlik2=1000000.;
    oldlikA=-1000000.;
    lik=0.
    

    noerror = true; devnostop = true;  sum_devs = 0;
    while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
        it = it+1;
        param=param01;

        T00 = 1
        plimix = E_stepMixthinv(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        param01 = M_parsMixthinv(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
        dt = time() - t
        if param01 == false
            return return(vcat(param01,vec(plimix),vec(T),0.,0.))
        end
        lik = E_stepMixthinv(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);

        #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
        con = (lik-oldlik)/(oldlik-oldlik2);
        likA = oldlik + (lik-oldlik)/(1-con)
        old_dev = dev1
        dev1 = abs(likA-oldlikA)
        dev2 = abs(lik-oldlik)
        println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
        println(param01)
        dev = abs(lik-oldlik)/abs(oldlik)

        noerror = oldlik<lik

        oldlikA=likA
        oldlik2=oldlik
        oldlik=lik;
        pml_em_log("EMparsMixthinv", it, lik, dev, dev1, abs(lik - oldlik2), tol, noerror)

        # sum_devs += dev > old_dev
        # devstop = sum_devs <= 3
    end
    
    # classificazioni dopo ultima massimizzazione
    plimix = E_stepMixthinv(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

    if noerror == true && it<=maxit 
        return vcat(param01,vec(plimix),vec(T),lik,1.)
    else
        return vcat(param01,vec(plimix),vec(T),0.,0.)
    end
end





function E_stepMixthinv(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
    if sum(param .== 0.) == size(param,1)
        return -Inf
    end
    
    trsi = (c-1)
    t1s = param[1:trsi]
    L1 = param[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
    
    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )
    Psi = diagm(param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    kappa = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)]   # interc2 = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
    betas = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+ncov)]  #  betas = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)

    mu = L2*kappa
    
    Sigms1=L1*Phi*L1' + diagm(ones(p))
    Sigms2=L2*L2'+Psi  
    
    t1=vcat(-Inf,t1s,Inf)
    tfa1=t1[X.+1]
    tfa=t1[X]
   
    
    if estep == true
        plimix=zeros(n,2,size(pairsindex,2));
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sigms1[pairsindex[:,idx],pairsindex[:,idx]]
            Rcfa=cov2cor(Scfa);
            Bcfa=hcat(tfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            
            
            #EFA
            Sefa=Sigms2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            

            plimix[:,1,idx] = tcfa3 .* pg[:,1]
            plimix[:,2,idx] = tefa3 .* pg[:,2]
        end
    # % probabilities in each row must sum to 1 (application of Bayes law)
        return plimix
    elseif estep == false
        plimixll=zeros(Float64,n);
        # plimix=zeros(n,2,size(pairsindex,2));
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sigms1[pairsindex[:,idx],pairsindex[:,idx]]
            Rcfa=cov2cor(Scfa);
            sdScfa=sqrt.(diag(Scfa)')
            tcfaidx = tfa[:,pairsindex[:,idx]]
            tcfa1idx=tfa1[:,pairsindex[:,idx]]
            Bcfa=hcat(tcfaidx./sdScfa,tcfa1idx./sdScfa);
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            
            
            #EFA
            Sefa=Sigms2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            sdSefa=sqrt.(diag(Sefa)')
            muij=mu[pairsindex[:,idx]]'
            tefaidx = tfa[:,pairsindex[:,idx]]
            tefa1idx = tfa1[:,pairsindex[:,idx]]
            Befa=hcat((tefaidx.-muij)./sdSefa,(tefa1idx.-muij)./sdSefa);
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            
            # plimix[:,1,idx] = tcfa3 .* pg[:,1]
            # plimix[:,2,idx] = tefa3 .* pg[:,2]
            plimixll += log.(tcfa3.*pg[:,1].+tefa3.*pg[:,2])
        end
        # cllik=sum(log.(sum(plimix,dims=2))) 
        cllik=sum(plimixll) 
        return cllik
    end
end
    

    


function M_parsMixthinv(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)

    trsi = (c-1) #attenzione p è già qui
   
    # PML - CFA & EFA
    theta_init = theta[1:(size(theta,1)-ncov)]#0.000001
    # # -1 loadings
    # lb=vec(vcat(fill(-5.1,trsi),fill(0.00001,p),fill(-1.,Int((q*(q-1))/2)),fill(1.,trs),fill(-1.,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-10,p)));
    # # 0. loadings
    lb=vec(vcat(fill(-Inf,trsi),fill(-Inf,p),fill(-1.,Int((q*(q-1))/2)),fill(-Inf,Int(p*K-K*(K-1)/2)),fill(1.0e-6,p),fill(-Inf,K)));#1.0e-5
    ub=vec(vcat(fill(Inf,trsi),fill(Inf,p),fill(1.,Int((q*(q-1))/2)),fill(Inf,Int(p*K-K*(K-1)/2)),fill(Inf,p),fill(Inf,K)));

    par = zeros(Float64,size(lb))
    if it == 1
        par = optimize(Optim.only_fg!((F,G,x)->fgthinv!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min ,Optim.Options(iterations=1000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0))
    elseif it > 1
        par = optimize(Optim.only_fg!((F,G,x)->fgthinv!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min
    end 


    # par = optimize(x->pmllast(x,X,T,p,q,K,c,n,pairsindex,L_str,indices),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

    # par = zeros(Float64,size(lb))
    # if it == 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 300.0)) #20 min
    # elseif it > 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 1200.0)) #20 min
    # end 
    
    # mixture parameter
    betas = theta[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
    Ez1 = T[:,1,:]
    Ez0 = T[:,2,:]
    
    pairs = size(pairsindex,2)
    betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)

    
    # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
   
    # println(par.minimizer,pg[1])
    if !(Optim.converged(par)) && isnan(betas[1]) 
        println(betas[1])
        if start == true
            return zeros(Float64,size(theta,1))
        elseif start == false
            return Optim.converged(par)
        end
    end
    return vcat(par.minimizer,betasnew)#,pg[1])
end



function pmllastthinv(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    trsi = (c-1)
    t1s = x[1:trsi]
    L1 = x[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])

    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )
    Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    kappa = x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)]    #interc2 = x[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p)]
    # prop = x[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p+1)]

    mu = L2*kappa
    
    Sig1=L1*Phi*L1' + diagm(ones(p))
    Sig2=L2*L2'+Psi

    plimix=zeros(Float64,n);
    t1=vcat(-Inf,t1s,Inf)
    tfa1=t1[X.+1]
    tfa=t1[X]
   
    
    for idx=1:size(pairsindex,2) #parallelizzabile
        # CFA
        Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
        Rcfa=cov2cor(Scfa);
        Bcfa=hcat(tfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
        tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
        Tt1 = T[:,1,idx]
        
        #EFA
        Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
        Refa=cov2cor(Sefa);
        Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
        tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
        Tt2 = T[:,2,idx]

        plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
    end
    f=-sum(plimix)
return f
end

function fgthinv!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trsi = (c-1)
    t1s = x[1:trsi]
    L1 = x[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
   
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )

    Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p-1))]))#K
    kappa = x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)]    #interc2 = x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
    # prop = x[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p+1)]

    t1=vcat(-Inf,t1s,Inf)
   
    mu = L2*kappa

    Sig1= L1*Phi*L1' + diagm(ones(p))
    Sig2= L2*L2'+Psi

    if G !== nothing

        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradkappacefa = zeros(Float64,K,size(pairsindex,2))
        gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)
        
        gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)
        
        
        gradlscfaidx = zeros(Float64,q,size(pairsindex,2),2)
        gradlambci = zeros(Float64,q,c+1,size(pairsindex,2),2)
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
            
        idx = 0
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1
        
                sig1s = hcat(vcat(Sig1[i,i],Sig1[i,j]),vcat(Sig1[j,i],Sig1[j,j]))
                cormat1 = cov2cor(sig1s)
                rho1 = cormat1[1,2]
                sig1sij = diag(sig1s)   
                sig1ij = sqrt.(sig1sij)
                sig1ij13 = sig1ij[1]^(3)
                sig1ij23 = sig1ij[2]^(3)
                L1i=L1[i,:]
                L1j=L1[j,:]
        
                sig2s = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sig2s)
                rho2 = cormat2[1,2]
                muij = mu[vcat(i,j)]
                sig2sij = diag(sig2s)
                sig2ij = sqrt.(sig2sij)
                sig2ij13 = sig2ij[1]^(3)
                sig2ij23 = sig2ij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]
                Psii=Psi[i,i]
                Psij=Psi[j,j]
        
                @inbounds for N in 1:n
                    xni = X[N,i]
                    xnj = X[N,j]
                    t1i = t1[xni]
                    t1i1 = t1[xni+1]
                    t1j = t1[xnj]
                    t1j1 = t1[xnj+1]
                    t1ist = (t1i)/sig1ij[1]
                    t1i1st = (t1i1)/sig1ij[1]
                    t1jst = (t1j)/sig1ij[2]
                    t1j1st = (t1j1)/sig1ij[2]
                    T1 = T[N,1,idx][1]
                    T2 = T[N,2,idx][1]
                    t2i = t1[xni]
                    t2i1 = t1[xni+1]
                    t2j = t1[xnj]
                    t2j1 = t1[xnj+1]
                    t2ist = (t2i-muij[1])/sig2ij[1]
                    t2i1st = (t2i1-muij[1])/sig2ij[1]
                    t2jst = (t2j-muij[2])/sig2ij[2]
                    t2j1st = (t2j1-muij[2])/sig2ij[2]
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                
                    pijicfa = bvnP(t1ist, t1i1st, t1jst, t1j1st, rho1, false)
                    derivNorm1i1 = (T1/pijicfa)*(pdf(Normal(0,1),t1i1st)*(cdf(Normal(0,1),(t1j1st-rho1*t1i1st)/sqrtrij1)-cdf(Normal(0,1),(t1jst-rho1*t1i1st)/sqrtrij1)))
                    derivNorm1i = (T1/pijicfa)*(pdf(Normal(0,1),t1ist)*(cdf(Normal(0,1),(t1jst-rho1*t1ist)/sqrtrij1)-cdf(Normal(0,1),(t1j1st-rho1*t1ist)/sqrtrij1)))
                    derivNorm1j1 = (T1/pijicfa)*(pdf(Normal(0,1),t1j1st)*(cdf(Normal(0,1),(t1i1st-rho1*t1j1st)/sqrtrij1)-cdf(Normal(0,1),(t1ist-rho1*t1j1st)/sqrtrij1)))
                    derivNorm1j = (T1/pijicfa)*(pdf(Normal(0,1),t1jst)*(cdf(Normal(0,1),(t1ist-rho1*t1jst)/sqrtrij1)-cdf(Normal(0,1),(t1i1st-rho1*t1jst)/sqrtrij1)))
                
                    gradtaucfa[xni+1,idx,1] +=  derivNorm1i1/(sig1ij[1])
                    gradtaucfa[xni,idx,1] += derivNorm1i/(sig1ij[1])
                    gradtaucfa[xnj+1,idx,2] += derivNorm1j1/(sig1ij[2])
                    gradtaucfa[xnj,idx,2] += derivNorm1j/(sig1ij[2])
        
        
                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))
        
                    if isnan(derivNormi1)
                        derivNormi1 = 0.
                    end
                    if isnan(derivNormi)
                        derivNormi = 0.
                    end
                    if isnan(derivNormj1)
                        derivNormj1 = 0.
                    end
                    if isnan(derivNormj)
                        derivNormj = 0.
                    end
        
                    
                    if !(t1ist == -Inf && t1jst == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1jst))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else                  
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1st == Inf && t1j1st == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1j1st))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end
        
                    if !(t1i1st == Inf && t1jst == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1jst))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end
        
                    if !(t1ist == -Inf && t1j1st == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1j1st))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    dpidrij1 = (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)#*vcat((-((L1j'*Phi.*t1j1st)./(2*sig1sij[2]))),(-((L1i'*Phi.*t1i1st)./(2*sig1sij[1]))))#vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)#
                
        
                    gradlscfaidx[1:q,idx,1] .+= dpidrij1*(Phi*L1j)/(sig1ij13*sig1ij[2])#(-Phi*L1i*(L1i'*L1j)/(sig1ij13*sig1ij[2])+Phi*L1j/((sig1ij[1]*sig1ij[2])))
                    gradlscfaidx[1:q,idx,2] .+= dpidrij1*(Phi*L1i)/(sig1ij23*sig1ij[1])#(-Phi*L1j*(L1j'*L1i)/(sig1ij23*sig1ij[1])+Phi*L1i/((sig1ij[2]*sig1ij[1])))
                    
                    gradlambci[1:q,xni+1,idx,1] .+=  derivNorm1i1 *(-((L1i*t1i1st)./(sig1sij[1])))#(-((L1i*t1i1).+(kappak*Psii))./((sig1ij13)))
                    gradlambci[1:q,xnj+1,idx,2] .+=  derivNorm1j1 *(-((L1j*t1j1st)./(sig1sij[2])))#(-((L1j*t1j1).+(kappak*Psij))./((sig1ij23)))
                    gradlambci[1:q,xni,idx,1]   .+=  derivNorm1i  *(-((L1i*t1ist)./(sig1sij[1])))#(-((L1i*t1i).+(kappak*Psii))./((sig1ij13)))
                    gradlambci[1:q,xnj,idx,2]   .+=  derivNorm1j  *(-((L1j*t1jst)./(sig1sij[2])))#(-((L1j*t1j).+(kappak*Psij))./((sig1ij23)))               
                    
                    gradphici[1:q,1:q,idx] .+= dpidrij1*(L1i*L1j')/(sig1ij[1]*sig1ij[2])
        
                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sig2ij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sig2ij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sig2ij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sig2ij[2])
        
        
                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)
        
                    # gradkappacefa[xni+1,1:K,idx,1] += (derivNormi1*(-L2[i,:]/sig2ij[1]))
                    # gradkappacefa[xnj+1,1:K,idx,2] += (derivNormj1*(-L2[j,:]/sig2ij[2]))
                    # gradkappacefa[xni,1:K,idx,1] += (derivNormi*(-L2[i,:]/sig2ij[1]))
                    # gradkappacefa[xnj,1:K,idx,2] += (derivNormj*(-L2[j,:]/sig2ij[2]))
                    gradkappacefa[1:K,idx] .+= derivNormi1*(-L2[i,:]/sig2ij[1]) + derivNormi*(-L2[i,:]/sig2ij[1]) + derivNormj1*(-L2[j,:]/sig2ij[2]) + derivNormj*(-L2[j,:]/sig2ij[2])
                    
                    gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sig2ij13)*(sig2ij[2])))
                    gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sig2ij23)*(sig2ij[1])))
        
                    gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij[1]-t2i1))/(2*(sig2ij13))
                    gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij[2]-t2j1))/(2*(sig2ij23))
                    gradpsiefa[xni,idx,1] += (derivNormi*(muij[1]-t2i))/(2*(sig2ij13))
                    gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij[2]-t2j))/(2*(sig2ij23))
                    
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        kappak=kappa[k]
                        gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sig2ij13*sig2ij[2])+L2jk/((sig2ij[1]*sig2ij[2])))#((L2j*Psi[i,i])/(sig2ij13*sig2ij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sig2ij23*sig2ij[1])+L2ik/((sig2ij[2]*sig2ij[1])))#((L2[i,:].*Psi[j,j])./(sig2ij23*sig2ij[1]))
        
                        gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(-((L2ik*t2i1st)./(sig2sij[1]).+(kappak)./(sig2ij[1])))#(-((L2ik*t2i1st).+(kappak*Psii))./((sig2ij13)))
                        gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(-((L2jk*t2j1st)./(sig2sij[2]).+(kappak)./(sig2ij[2])))#(-((L2jk*t2j1st).+(kappak*Psij))./((sig2ij23)))
                        gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(-((L2ik*t2ist)./(sig2sij[1]).+(kappak)./(sig2ij[1])))#(-((L2ik*t2ist).+(kappak*Psii))./((sig2ij13)))
                        gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(-((L2jk*t2jst)./(sig2sij[2]).+(kappak)./(sig2ij[2])))#(-((L2jk*t2jst).+(kappak*Psij))./((sig2ij23)))
                    end 
                
                end
            end
        end
        
        vecpair = vec(pairsindex)
        liktaushared = -sum(gradtaucfa,dims=(2,3))[2:(c)] -sum(gradtauefa,dims=(2,3))[2:(c)]
        
        
        
        
        
        likls1 = zeros(p,q); llamb1 = sum(gradlambci[:,2:c,:,:],dims=2)[:,1,:,:]; llam1 = zeros(p,q); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
        likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
        llam1idx  = zeros(p,q); likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        @inbounds for i in 1:p
            @inbounds for qs in 1:q
                llam1[i,qs] = -sum(vec(llamb1[qs,:,:]')[vecpair.==i])
                llam1idx[i,qs] = -sum(vec(gradlscfaidx[qs,:,:]')[vecpair.==i])
            end
            likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
            likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end
        
        likpsitrue = likpsi.+likpsiidx
        liklstrue = likls.+liklsidx
        liklamtrue = llam1.+llam1idx
        
        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]
        
        likappa = -sum(gradkappacefa,dims=2)
        
        G .= vec(vcat(liktaushared,sum(liklamtrue.*L_str,dims=2),likphida,liklstrue[indices],likpsitrue,likappa))#likpsitrue[2:end]    
      
    end
    if F !== nothing 
  
    
        plimix=zeros(Float64,n);    
        tfa1=t1[X.+1]
        tfa=t1[X]
       
       
       
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
            Rcfa=cov2cor(Scfa);
            Bcfa=hcat(tfa[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'),tfa1[:,pairsindex[:,idx]]./sqrt.(diag(Scfa)'));
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Rcfa[1,2],false);
            Tt1 = T[:,1,idx]
            
            #EFA
            Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            Tt2 = T[:,2,idx]
            plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2;
            
        end
        f=-sum(plimix)
        # @info "Current objective: $f"
        return f
    end
end




# mixed theta and delta parameterization
# ########################################################################### Latent means EFA + ONLY ONE SET OF invariant thresholds + Psi estimated ##################################################################################################


# function startvalsthinv(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
#     liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)+ncov,replix);#p-1
#     for i in 1:replix
#         println(i)

#         tau1 = rand(Uniform(-2.3,-1.)).+repeat(collect(range(0,3.5,c-1)))
#         lambda1 = rand(Uniform(0.5,.90),p);
#         phi=rand(LKJ(q,1))
#         phi=phi[tril!(trues(size(phi)), -1)]
#         lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(0.15,.90),Int((p*K-(K*(K-1)/2))))))
#         lambda2=lamb2[lamb2.!=0.]
#         psi=rand(Uniform(0.2,1.5),p)
#         kappa=rand(Uniform(-0.2,2),K)



#         prop = rand(Uniform(0.65,0.95))
#         bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.1,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


#         param_initials = vcat(tau1,lambda1,phi,lambda2,psi,kappa,bets)

#         # first iteration
#         plimix = E_stepMixthinv(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
#         T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
#         t = time()
#         params[:,i] = M_parsMixthinv(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
#         dt = time() - t
#         prrot=params[:,i]
#         liksmix[i] = E_stepMixthinv(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
#     end
#     param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
#     lik=liksmix[findall(liksmix.==maximum(liksmix))]

#     tau1=param_initials0[(1):((c-1))]
#     tau1=vec(sort(reshape(tau1,1,c-1),dims=2))
#     lambda1=param_initials0[((c-1)+1):((c-1)+p)]
#     phi=param_initials0[((c-1)+p+1):((c-1)+p+Int(q*(q-1)/2))]
#     lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
#     lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
#     phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
#     lambda2=param_initials0[((c-1)+p+Int(q*(q-1)/2)+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2))))]
#     psi=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p))]# psi=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1))]
#     kappa=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+K)]#interc=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
#     lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
#     lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
#     psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
#     bets = param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+K+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+K+ncov)]#    bets = param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


#     param_initialii = vcat(tau1,lambda1,phi,lambda2,psi,kappa,bets)

#     println(param_initialii,liksmix,"first iteration: done")
#     return param_initialii
# end


# function EMparsMixthinv(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
#     local plimix::Array, T,oldlik
#     n = size(X,1)
#     er=0;
#     # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
#     it=0;
#     dev=Inf; dev1=Inf; dev2=Inf;

    
#     oldlik=-1000000.;
#     oldlik2=1000000.;
#     oldlikA=-1000000.;
#     lik=0.
    

#     noerror = true; devnostop = true;  sum_devs = 0;
#     while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
#         it = it+1;
#         param=param01;

#         T00 = 1
#         plimix = E_stepMixthinv(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
#         T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
#         t = time()
#         param01 = M_parsMixthinv(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
#         dt = time() - t
#         if param01 == false
#             return return(vcat(param01,vec(plimix),vec(T),0.,0.))
#         end
#         lik = E_stepMixthinv(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);

#         #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
#         con = (lik-oldlik)/(oldlik-oldlik2);
#         likA = oldlik + (lik-oldlik)/(1-con)
#         old_dev = dev1
#         dev1 = abs(likA-oldlikA)
#         dev2 = abs(lik-oldlik)
#         println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
#         println(param01)
#         dev = abs(lik-oldlik)/abs(oldlik)

#         noerror = oldlik<lik

#         oldlikA=likA
#         oldlik2=oldlik
#         oldlik=lik;

#         # sum_devs += dev > old_dev
#         # devstop = sum_devs <= 3
#     end
    
#     # classificazioni dopo ultima massimizzazione
#     plimix = E_stepMixthinv(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
#     T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

#     if noerror == true && it<=maxit 
#         return vcat(param01,vec(plimix),vec(T),lik,1.)
#     else
#         return vcat(param01,vec(plimix),vec(T),0.,0.)
#     end
# end





# function E_stepMixthinv(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
#     if sum(param .== 0.) == size(param,1)
#         return -Inf
#     end
    
#     trsi = (c-1)
#     t1s = param[1:trsi]
#     L1 = param[(trsi+1):(trsi+p)].*L_str 
#     Phi = reconstruct_correlation_matrix(param[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
    
#     L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )
#     Psi = diagm(param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p-1))]))#K
#     kappa = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)]   # interc2 = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
#     betas = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+ncov)]  #  betas = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

#     propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
#     pg = hcat(propo,1. .-propo)

#     mu = L2*kappa
    
#     Sigms=ones(p,p,2);
#     Sigms[:,:,1]=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
#     Sigms[:,:,2]=L2*L2'+Psi 
    
#     t1=vcat(-Inf,t1s,Inf)
#     tfa1=t1[X.+1]
#     tfa=t1[X]
   
    
#     if estep == true
#         low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#         plimix=ones(n,2,size(pairsindex,2));
    
#         for idx=1:size(pairsindex,2) #parallelizzabile
#             for g=1:2
#                 if g==1
#                     s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                     r2c=cov2cor(s2c);
#                     up=(tfa1[:,pairsindex[:,idx]])
#                     low=(tfa[:,pairsindex[:,idx]])
#                 elseif g==2
#                     s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                     r2c=cov2cor(s2c);
#                     up=(tfa1[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                     low=(tfa[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 end
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix[:,g,idx]= t3.*pg[:,g];
#             end
#         end
#     # % probabilities in each row must sum to 1 (application of Bayes law)
#         return plimix
#     elseif estep == false
#         low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#         plimix=ones(n,2,size(pairsindex,2));

#         for idx=1:size(pairsindex,2) #parallelizzabile
#             for g=1:2
#                 if g==1
#                     s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                     r2c=cov2cor(s2c);
#                     up=(tfa1[:,pairsindex[:,idx]])
#                     low=(tfa[:,pairsindex[:,idx]])
#                 elseif g==2
#                     s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                     r2c=cov2cor(s2c);
#                     up=(tfa1[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                     low=(tfa[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 end
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix[:,g,idx] = t3
#             end
#         end
#         lik1=(log.(plimix));
#         lik1=sum(lik1,dims=3);
#         U=exp.(lik1)[:,:,1].*pg; 
#         cllik=sum(log.(sum(U,dims=2)))
#         return cllik
#     end
# end
    

    


# function M_parsMixthinv(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)

#     trsi = (c-1) #attenzione p è già qui
   
#     # PML - CFA & EFA
#     theta_init = theta[1:(size(theta,1)-ncov)]#0.000001
#     # # -1 loadings
#     # lb=vec(vcat(fill(-5.1,trsi),fill(0.00001,p),fill(-1.,Int((q*(q-1))/2)),fill(1.,trs),fill(-1.,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-10,p)));
#     # # 0. loadings
#     lb=vec(vcat(fill(-5.1,trsi),fill(-1,p),fill(-1.,Int((q*(q-1))/2)),fill(-Inf,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-5,K)));#1.0e-5
#     ub=vec(vcat(fill(5.1,trsi),fill(1.,p),fill(1.,Int((q*(q-1))/2)),fill(Inf,Int(p*K-K*(K-1)/2)),fill(Inf,p),fill(10.,K)));

#     par = zeros(Float64,size(lb))
#     if it == 1
#         par = optimize(Optim.only_fg!((F,G,x)->fgthinv!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min ,Optim.Options(iterations=1000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0))
#     elseif it > 1
#         par = optimize(Optim.only_fg!((F,G,x)->fgthinv!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min
#     end 


#     # par = optimize(x->pmllast(x,X,T,p,q,K,c,n,pairsindex,L_str,indices),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

#     # par = zeros(Float64,size(lb))
#     # if it == 1
#     #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 300.0)) #20 min
#     # elseif it > 1
#     #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 1200.0)) #20 min
#     # end 
    
#     # mixture parameter
#     betas = theta[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
#     Ez1 = T[:,1,:]
#     Ez0 = T[:,2,:]
    
#     pairs = size(pairsindex,2)
#     betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)

    
#     # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
   
#     # println(par.minimizer,pg[1])
#     if !(Optim.converged(par)) && isnan(betas[1]) 
#         println(betas[1])
#         if start == true
#             return zeros(Float64,size(theta,1))
#         elseif start == false
#             return Optim.converged(par)
#         end
#     end
#     return vcat(par.minimizer,betasnew)#,pg[1])
# end



# function pmllastthinv(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
#     trsi = (c-1)
#     t1s = x[1:trsi]
#     L1 = x[(trsi+1):(trsi+p)].*L_str 
#     Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])

#     L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )
#     Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p-1))]))#K
#     kappa = x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)]    #interc2 = x[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p)]
#     # prop = x[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p+1)]

#     mu = L2*kappa
    
#     Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
#     Sig2=L2*L2'+Psi

#     plimix=zeros(Float64,n);
#     t1=vcat(-Inf,t1s,Inf)
#     tfa1=t1[X.+1]
#     tfa=t1[X]
   
    
#     for idx=1:size(pairsindex,2) #parallelizzabile
#         # CFA
#         Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
#         Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
#         tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
#         Tt1 = T[:,1,idx]
        
#         #EFA
#         Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
#         Refa=cov2cor(Sefa);
#         Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
#         tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
#         Tt2 = T[:,2,idx]

#         plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
#     end
#     f=-sum(plimix)
# return f
# end

# function fgthinv!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
#     # do common computations here
#     # ...
#     trsi = (c-1)
#     t1s = x[1:trsi]
#     L1 = x[(trsi+1):(trsi+p)].*L_str 
#     Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
   
#     L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )

#     Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p))])#Psi = diagm(vcat(1.,x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+(p-1))]))#K
#     kappa = x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+K)]    #interc2 = x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
#     # prop = x[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p+1)]

#     t1=vcat(-Inf,t1s,Inf)
   
#     mu = L2*kappa

#     Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
#     Sig2=L2*L2'+Psi


#     if G !== nothing
#         gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
#         gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
#         gradkappacefa = zeros(Float64,K,size(pairsindex,2))
#         gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)
    
#         gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
#         gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)
    
#         gradlambci = zeros(Float64,2,q,size(pairsindex,2))
#         gradphici = zeros(Float64,q,q,size(pairsindex,2))
#         gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
        
#         idx = 0
#         @inbounds for i in 1:(p-1)
#             @inbounds for j in (i+1):p
#                 idx += 1
#                 rho1 = Sig1[i,j]
#                 cormat1 = ones(Float64,2,2); cormat1[1,2]=rho1; cormat1[2,1]=rho1
#                 sigs = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
#                 cormat2 = cov2cor(sigs)
#                 rho2 = cormat2[1,2]
#                 muij = mu[vcat(i,j)]
#                 sigsij = diag(sigs)
#                 sigij = sqrt.(sigsij)
#                 sigij13 = sigij[1]^(3)
#                 sigij23 = sigij[2]^(3)
#                 L2i=L2[i,:]
#                 L2j=L2[j,:]
#                 Psii=Psi[i,i]
#                 Psij=Psi[j,j]
#                 @inbounds for N in 1:n
#                     xni = X[N,i]
#                     xnj = X[N,j]
#                     t1i = t1[xni]
#                     t1i1 = t1[xni+1]
#                     t1j = t1[xnj]
#                     t1j1 = t1[xnj+1]
#                     T1 = T[N,1,idx][1]
#                     T2 = T[N,2,idx][1]
#                     t2i = t1[xni]
#                     t2i1 = t1[xni+1]
#                     t2j = t1[xnj]
#                     t2j1 = t1[xnj+1]
#                     t2ist = (t2i-muij[1])/sigij[1]
#                     t2i1st = (t2i1-muij[1])/sigij[1]
#                     t2jst = (t2j-muij[2])/sigij[2]
#                     t2j1st = (t2j1-muij[2])/sigij[2]
#                     sqrtrij1 = sqrt(1-rho1^2)
#                     sqrtrij2 = sqrt(1-rho2^2)
                                   
#                     pijicfa = bvnP(t1i, t1i1, t1j, t1j1, rho1, false)
#                     gradtaucfa[xni+1,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i1)*(cdf(Normal(0,1),(t1j1-rho1*t1i1)/sqrtrij1)-cdf(Normal(0,1),(t1j-rho1*t1i1)/sqrtrij1))) 
#                     gradtaucfa[xni,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i)*(cdf(Normal(0,1),(t1j-rho1*t1i)/sqrtrij1)-cdf(Normal(0,1),(t1j1-rho1*t1i)/sqrtrij1))) 
#                     gradtaucfa[xnj+1,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j1)*(cdf(Normal(0,1),(t1i1-rho1*t1j1)/sqrtrij1)-cdf(Normal(0,1),(t1i-rho1*t1j1)/sqrtrij1))) 
#                     gradtaucfa[xnj,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j)*(cdf(Normal(0,1),(t1i-rho1*t1j)/sqrtrij1)-cdf(Normal(0,1),(t1i1-rho1*t1j)/sqrtrij1))) 
                    
    
#                     pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
#                     derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
#                     derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
#                     derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
#                     derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))
    
#                     if isnan(derivNormi1)
#                         derivNormi1 = 0.
#                     end
#                     if isnan(derivNormi)
#                         derivNormi = 0.
#                     end
#                     if isnan(derivNormj1)
#                         derivNormj1 = 0.
#                     end
#                     if isnan(derivNormj)
#                         derivNormj = 0.
#                     end
    
                    
#                     if !(t1i == -Inf && t1j == -Inf) 
#                         xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j))
#                         xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
#                     else                  
#                         xnixnjc = 0.0
#                         xnixnje = 0.0
#                     end
                        
#                     if !(t1i1 == Inf && t1j1 == Inf) 
#                         xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j1))
#                         xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
#                     else
#                         xni1xnj1c = 0.0
#                         xni1xnj1e = 0.0
#                     end
        
#                     if !(t1i1 == Inf && t1j == -Inf) 
#                         xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j))
#                         xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
#                     else
#                         xni1xnjc = 0.0
#                         xni1xnje = 0.0
#                     end
        
#                     if !(t1i == -Inf && t1j1 == Inf) 
#                         xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j1))
#                         xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
#                     else
#                         xnixnj1c = 0.0
#                         xnixnj1e = 0.0
#                     end
#                     gradlambci[1:2,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)
#                     gradphici[1:q,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*(L1[i,:]*L1[j,:]')
                    
#                     #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
#                     gradtauefa[xni+1,idx,1] += (derivNormi1)/(sigij[1])
#                     gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sigij[2])
#                     gradtauefa[xni,idx,1] += (derivNormi)/(sigij[1])
#                     gradtauefa[xnj,idx,2] += (derivNormj)/(sigij[2])
    
    
#                     dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)

#                     gradkappacefa[1:K,idx] .+= derivNormi1*(-L2[i,:]/sigij[1]) + derivNormi*(-L2[i,:]/sigij[1]) + derivNormj1*(-L2[j,:]/sigij[2]) + derivNormj*(-L2[j,:]/sigij[2])
                    
#                     gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sigij13)*(sigij[2])))
#                     gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sigij23)*(sigij[1])))
    
#                     gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij[1]-t2i1))/(2*(sigij13))
#                     gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij[2]-t2j1))/(2*(sigij23))
#                     gradpsiefa[xni,idx,1] += (derivNormi*(muij[1]-t2i))/(2*(sigij13))
#                     gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij[2]-t2j))/(2*(sigij23))
                    
#                     @inbounds for k in 1:K
#                         L2ik=L2[i,k]
#                         L2jk=L2[j,k]
#                         kappak=kappa[k]
#                         gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sigij13*sigij[2])+L2jk/((sigij[1]*sigij[2])))#((L2j*Psi[i,i])/(sigij13*sigij[2]))
#                         gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sigij23*sigij[1])+L2ik/((sigij[2]*sigij[1])))#((L2[i,:].*Psi[j,j])./(sigij23*sigij[1]))
    
#                         gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(-((L2ik*t2i1st)./(sigsij[1]).+(kappak)./(sigij[1])))#(-((L2ik*t2i1).+(kappak*Psii))./((sigij13)))
#                         gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(-((L2jk*t2j1st)./(sigsij[2]).+(kappak)./(sigij[2])))#(-((L2jk*t2j1).+(kappak*Psij))./((sigij23)))
#                         gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(-((L2ik*t2ist)./(sigsij[1]).+(kappak)./(sigij[1])))#(-((L2ik*t2i).+(kappak*Psii))./((sigij13)))
#                         gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(-((L2jk*t2jst)./(sigsij[2]).+(kappak)./(sigij[2])))#(-((L2jk*t2j).+(kappak*Psij))./((sigij23)))
#                     end 
                
#                 end
#             end
#         end

#         vecpair = vec(pairsindex)
#         liktaushared = -sum(gradtaucfa,dims=(2,3))[2:(c)] -sum(gradtauefa,dims=(2,3))[2:(c)]
        
        
        
    
        
#         llamb1 = permutedims(gradlambci,(2,3,1)); llam1 = zeros(p); ids = find_ones_columns(L_str);
#         likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
#         likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
#         likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
#         @inbounds for i in 1:p
#             llam1[i] = -sum(vec(llamb1[ids[i],:,:]')[vecpair.==i])
#             likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
#             likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
#             @inbounds for k in 1:K
#                 likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
#                 liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
#             end
#         end
    
#         likpsitrue = likpsi.+likpsiidx
#         liklstrue = likls.+liklsidx
    
#         likp = sum(gradphici,dims=(3))[:,:,1]
#         likphida = -likp[triu!(trues(size(likp)), 1)]
    
#         likappa = -sum(gradkappacefa,dims=2)

#         G .= vec(vcat(liktaushared,llam1,likphida,liklstrue[indices],likpsitrue,likappa))#likpsitrue[2:end]
      
#     end
#     if F !== nothing 
  
    
#         plimix=zeros(Float64,n);    
#         tfa1=t1[X.+1]
#         tfa=t1[X]
       
       
       
#         for idx=1:size(pairsindex,2) #parallelizzabile
#             # CFA
#             Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
#             Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
#             tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
#             Tt1 = T[:,1,idx]
            
#             #EFA
#             Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
#             Refa=cov2cor(Sefa);
#             Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
#             tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
#             Tt2 = T[:,2,idx]
#             plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2;
            
#         end
#         f=-sum(plimix)
#         # @info "Current objective: $f"
#         return f
#     end
# end






########################################################################### Latent means EFA no psi and invariant thersholds ##################################################################################################

function startvalslminv(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
    liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K)+ncov,replix);#p-1
    for i in 1:replix
        println(i)

        tau1 = vec(hcat(rand(Uniform(-2.3,-1.),p).+reshape(repeat(collect(range(0,3.5,c-1)),inner=p),p,(c-1))))
        lambda1 = rand(Uniform(0.5,.90),p);
        phi=rand(LKJ(q,1))
        phi=phi[tril!(trues(size(phi)), -1)]
        lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(0.15,.90),Int((p*K-(K*(K-1)/2))))))
        lambda2=lamb2[lamb2.!=0.]
        kappa=rand(Uniform(0.5,3),K)



        prop = rand(Uniform(0.65,0.95))
        bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.1,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


        param_initials = vcat(tau1,lambda1,phi,lambda2,kappa,bets)

        # first iteration
        plimix = E_stepMixlminv(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        params[:,i] = M_parsMixlminv(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        dt = time() - t
        prrot=params[:,i]
        liksmix[i] = E_stepMixlminv(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
    end
    param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    lik=liksmix[findall(liksmix.==maximum(liksmix))]

    tau1=param_initials0[(1):((c-1)*p)]
    tau1=vec(sort(reshape(tau1,p,c-1),dims=2))
    lambda1=param_initials0[((c-1)*p+1):((c-1)*p+p)]
    phi=param_initials0[((c-1)*p+p+1):((c-1)*p+p+Int(q*(q-1)/2))]
    lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    lambda2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2))))]
    kappa=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+K)]#interc=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+K+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+K+ncov)]#    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


    param_initialii = vcat(tau1,lambda1,phi,lambda2,kappa,bets)

    println(param_initialii,liksmix,"first iteration: done")
    return param_initialii
end


function EMparsMixlminv(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
    local plimix::Array, T,oldlik
    n = size(X,1)
    er=0;
    # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
    it=0;
    dev=Inf; dev1=Inf; dev2=Inf;

    
    oldlik=-1000000.;
    oldlik2=1000000.;
    oldlikA=-1000000.;
    lik=0.
    

    noerror = true; devnostop = true;  sum_devs = 0;
    while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
        it = it+1;
        param=param01;

        T00 = 1
        plimix = E_stepMixlminv(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        param01 = M_parsMixlminv(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
        dt = time() - t
        if param01 == false
            return return(vcat(param01,vec(plimix),vec(T),0.,0.))
        end
        lik = E_stepMixlminv(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);

        #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
        con = (lik-oldlik)/(oldlik-oldlik2);
        likA = oldlik + (lik-oldlik)/(1-con)
        old_dev = dev1
        dev1 = abs(likA-oldlikA)
        dev2 = abs(lik-oldlik)
        println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
        println(param01)
        dev = abs(lik-oldlik)/abs(oldlik)

        noerror = oldlik<lik

        oldlikA=likA
        oldlik2=oldlik
        oldlik=lik;
        pml_em_log("EMparsMixlminv", it, lik, dev, dev1, abs(lik - oldlik2), tol, noerror)

        # sum_devs += dev > old_dev
        # devstop = sum_devs <= 3
    end
    
    # classificazioni dopo ultima massimizzazione
    plimix = E_stepMixlminv(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

    if noerror == true && it<=maxit 
        return vcat(param01,vec(plimix),vec(T),lik,1.)
    else
        return vcat(param01,vec(plimix),vec(T),0.,0.)
    end
end





function E_stepMixlminv(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
    if sum(param .== 0.) == size(param,1)
        return -Inf
    end
    trsi = (c-1)
    t1s = param[1:trsi*p]
    L1 = param[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    
    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )

    kappa = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K)]   # interc2 = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
    betas = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K+ncov)]  #  betas = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)

    mu = L2*kappa
    
    Sigms=ones(p,p,2);
    Sigms[:,:,1]=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sigms[:,:,2]=L2*L2'+ Diagonal(1 .-Diagonal((L2*L2')))
    
    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))

    tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]

    if estep == true
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));
    
        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]])
                    low=(tfa[:,pairsindex[:,idx]])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')
                    low=(tfa[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx] = t3.*pg[:,g];
            end
        end
    # % probabilities in each row must sum to 1 (application of Bayes law)
        return plimix
    elseif estep == false
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));

        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]])
                    low=(tfa[:,pairsindex[:,idx]])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')
                    low=(tfa[:,pairsindex[:,idx]].-(mu)[pairsindex[:,idx]]')
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx] = t3
            end
        end
        lik1=(log.(plimix));
        lik1=sum(lik1,dims=3);
        U=exp.(lik1)[:,:,1].*pg; 
        cllik=sum(log.(sum(U,dims=2)))
        return cllik
    end
end
    

    


function M_parsMixlminv(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)

    trsi = (c-1)*p #attenzione p è già qui


    # PML - CFA & EFA
    theta_init = theta[1:(size(theta,1)-ncov)]#0.000001
    # # -1 loadings
    # lb=vec(vcat(fill(-5.1,trsi),fill(0.00001,p),fill(-1.,Int((q*(q-1))/2)),fill(1.,trs),fill(-1.,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-10,p)));
    # # 0. loadings
    lb=vec(vcat(fill(-5.1,trsi),fill(-1,p),fill(-1.,Int((q*(q-1))/2)),fill(-1,Int(p*K-K*(K-1)/2)),fill(-5,K)));#1.0e-5
    ub=vec(vcat(fill(5.1,trsi),fill(1.,p),fill(1.,Int((q*(q-1))/2)),fill(1,Int(p*K-K*(K-1)/2)),fill(10.,K)));

    par = zeros(Float64,size(lb))
    if it == 1
        par = optimize(Optim.only_fg!((F,G,x)->fglminv!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min ,Optim.Options(iterations=1000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0))
    elseif it > 1
        par = optimize(Optim.only_fg!((F,G,x)->fglminv!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min
    end 


    # par = optimize(x->pmllast(x,X,T,p,q,K,c,n,pairsindex,L_str,indices),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

    # par = zeros(Float64,size(lb))
    # if it == 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 300.0)) #20 min
    # elseif it > 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 1200.0)) #20 min
    # end 
    
    # mixture parameter
    betas = theta[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
    Ez1 = T[:,1,:]
    Ez0 = T[:,2,:]
    
    pairs = size(pairsindex,2)
    betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)

    
    # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
   
    # println(par.minimizer,pg[1])
    if !(Optim.converged(par)) && isnan(betas[1]) 
        println(betas[1])
        if start == true
            return zeros(Float64,size(theta,1))
        elseif start == false
            return Optim.converged(par)
        end
    end
    return vcat(par.minimizer,betasnew)#,pg[1])
end



function pmllastlminv(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )
    kappa = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K)]    #interc2 = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p)]
    # prop = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]

    

    mu = L2*kappa
    
    Sig1=L1*Phi*L1' + Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2' + Diagonal(1 .-Diagonal((L2*L2')))

   
    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))

    tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
    
    plimix=zeros(Float64,n);
    for idx=1:size(pairsindex,2) #parallelizzabile
        # CFA
        Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
        Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
        tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
        Tt1 = T[:,1,idx]
        
        #EFA
        Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
        Refa=cov2cor(Sefa);
        Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]'));
        tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
        Tt2 = T[:,2,idx]

        plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
    end
    f=-sum(plimix)
return f
end

function fglminv!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))])
    kappa = x[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+K)]    #interc2 = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p)]
    # prop = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]

    mu = L2*kappa

    Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+ Diagonal(1 .-Diagonal((L2*L2')))


    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))

   

    if G !== nothing
        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradkappacefa = zeros(Float64,K,size(pairsindex,2))
        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)
    
        gradlambci = zeros(Float64,2,q,size(pairsindex,2))
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        # gradlamb2ci = zeros(Float64,2,K,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)

        idx = 0
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1
                rho1 = Sig1[i,j]
                cormat1 = ones(Float64,2,2); cormat1[1,2]=rho1; cormat1[2,1]=rho1
                sigs = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sigs)
                rho2 = cormat2[1,2]
                muij = mu[vcat(i,j)]
                sigsij = diag(sigs)
                sigij = sqrt.(sigsij)
                sigij13 = sigij[1]^(3)
                sigij23 = sigij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]
                @inbounds for N in 1:n
                    xni = X[N,i]
                    xnj = X[N,j]
                    t1i = t1[i,xni]
                    t1i1 = t1[i,xni+1]
                    t1j = t1[j,xnj]
                    t1j1 = t1[j,xnj+1]
                    T1 = T[N,1,idx][1]
                    T2 = T[N,2,idx][1]
                    t2i = t1[i,xni]
                    t2i1 = t1[i,xni+1]
                    t2j = t1[j,xnj]
                    t2j1 = t1[j,xnj+1]
                    t2ist = (t2i-muij[1])/sigij[1]
                    t2i1st = (t2i1-muij[1])/sigij[1]
                    t2jst = (t2j-muij[2])/sigij[2]
                    t2j1st = (t2j1-muij[2])/sigij[2]
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                   
                    pijicfa = bvnP(t1i, t1i1, t1j, t1j1, rho1, false)
                    gradtaucfa[xni+1,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i1)*(cdf(Normal(0,1),(t1j1-rho1*t1i1)/sqrtrij1)-cdf(Normal(0,1),(t1j-rho1*t1i1)/sqrtrij1))) 
                    gradtaucfa[xni,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i)*(cdf(Normal(0,1),(t1j-rho1*t1i)/sqrtrij1)-cdf(Normal(0,1),(t1j1-rho1*t1i)/sqrtrij1))) 
                    gradtaucfa[xnj+1,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j1)*(cdf(Normal(0,1),(t1i1-rho1*t1j1)/sqrtrij1)-cdf(Normal(0,1),(t1i-rho1*t1j1)/sqrtrij1))) 
                    gradtaucfa[xnj,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j)*(cdf(Normal(0,1),(t1i-rho1*t1j)/sqrtrij1)-cdf(Normal(0,1),(t1i1-rho1*t1j)/sqrtrij1))) 
                    
    
                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))
    
                    if isnan(derivNormi1)
                        derivNormi1 = 0.
                    end
                    if isnan(derivNormi)
                        derivNormi = 0.
                    end
                    if isnan(derivNormj1)
                        derivNormj1 = 0.
                    end
                    if isnan(derivNormj)
                        derivNormj = 0.
                    end
    
                    
                    if !(t1i == -Inf && t1j == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else                  
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1 == Inf && t1j1 == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j1))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end
        
                    if !(t1i1 == Inf && t1j == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end
        
                    if !(t1i == -Inf && t1j1 == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j1))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    gradlambci[1:2,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)
                    gradphici[1:q,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*(L1[i,:]*L1[j,:]')
                    
                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sigij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sigij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sigij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sigij[2])
    
    
                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)
    
                    # gradkappacefa[xni+1,1:K,idx,1] += (derivNormi1*(-L2[i,:]/sigij[1]))
                    # gradkappacefa[xnj+1,1:K,idx,2] += (derivNormj1*(-L2[j,:]/sigij[2]))
                    # gradkappacefa[xni,1:K,idx,1] += (derivNormi*(-L2[i,:]/sigij[1]))
                    # gradkappacefa[xnj,1:K,idx,2] += (derivNormj*(-L2[j,:]/sigij[2]))
                    gradkappacefa[1:K,idx] .+= derivNormi1*(-L2[i,:]/sigij[1]) + derivNormi*(-L2[i,:]/sigij[1]) + derivNormj1*(-L2[j,:]/sigij[2]) + derivNormj*(-L2[j,:]/sigij[2])
                    
                    # gradlamb2ci[1:2,1:K,idx] .+=  dpidrij*vcat(L2[j,:]',L2[i,:]')   .+  kappa .* vcat(derivNormi1 + derivNormi, derivNormj1 + derivNormj)      
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        kappak=kappa[k]
                        gradlsefaidx[k,idx,1] += dpidrij*(L2jk)#((L2j*Psi[i,i])/(sigij13*sigij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(L2ik)#((L2[i,:].*Psi[j,j])./(sigij23*sigij[1]))
    
                        gradlamb2ci[k,xni+1,idx,1] -=  derivNormi1 * kappak
                        gradlamb2ci[k,xnj+1,idx,2] -=  derivNormj1 * kappak
                        gradlamb2ci[k,xni,idx,1]   -=  derivNormi  * kappak
                        gradlamb2ci[k,xnj,idx,2]   -=  derivNormj  * kappak
                    end 
                          
                end
            end
        end
    
        vecpair = vec(pairsindex)
        liktaushared = zeros(p*(c-1));idpcfa = 0   
        @inbounds for j in 2:c
            @inbounds for i in 1:p
                idpcfa += 1
                liktaushared[idpcfa] = -sum(vec(vcat(gradtaucfa[j,:,1]',gradtaucfa[j,:,2]'))[vecpair.==i])-sum(vec(vcat(gradtauefa[j,:,1]',gradtauefa[j,:,2]'))[vecpair.==i])
            end
        end
        
                
        llamb1 = permutedims(gradlambci,(2,3,1)); llam1 = zeros(p); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:]; liklsidx  = zeros(p,K)
        @inbounds for i in 1:p
            llam1[i] = -sum(vec(llamb1[ids[i],:,:]')[vecpair.==i])
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end
        # llamb2 = permutedims(gradlamb2ci,(2,3,1)); llam2 = zeros(p);
        # liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        # @inbounds for i in 1:p
        #     llam1[i] = -sum(vec(llamb1[ids[i],:,:]')[vecpair.==i])
        #     llam2[i] = -sum(vec(llamb2)[vecpair.==i])
        # end
    
       
    
        liklstrue = likls.+liklsidx
        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]
    
        likappa = -sum(gradkappacefa,dims=2)


        G .= vec(vcat(liktaushared,llam1,likphida,liklstrue[indices],likappa))#likpsitrue[2:end]
      
    end
    if F !== nothing 
  
    
        plimix=zeros(Float64,n);
        tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
        tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
       
       
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
            Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
            Tt1 = T[:,1,idx]
            
            #EFA
            Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tfa[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-mu[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            Tt2 = T[:,2,idx]
            plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2;
            
        end
        f=-sum(plimix)
        # @info "Current objective: $f"
        return f
    end
end



########################################################################### Only Lambda2: No latent means EFA, no psi and invariant thersholds ##################################################################################################

function startvalsl2inv(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
    liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+ncov,replix);#p-1
    for i in 1:replix
        println(i)

        tau1 = vec(hcat(rand(Uniform(-2.3,-1.),p).+reshape(repeat(collect(range(0,3.5,c-1)),inner=p),p,(c-1))))
        lambda1 = rand(Uniform(0.5,.90),p);
        phi=rand(LKJ(q,1))
        phi=phi[tril!(trues(size(phi)), -1)]
        lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(0.15,.90),Int((p*K-(K*(K-1)/2))))))
        lambda2=lamb2[lamb2.!=0.]



        prop = rand(Uniform(0.65,0.95))
        bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.1,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


        param_initials = vcat(tau1,lambda1,phi,lambda2,bets)

        # first iteration
        plimix = E_stepMixl2inv(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        params[:,i] = M_parsMixl2inv(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        dt = time() - t
        prrot=params[:,i]
        liksmix[i] = E_stepMixl2inv(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
    end
    param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    lik=liksmix[findall(liksmix.==maximum(liksmix))]

    tau1=param_initials0[(1):((c-1)*p)]
    tau1=vec(sort(reshape(tau1,p,c-1),dims=2))
    lambda1=param_initials0[((c-1)*p+1):((c-1)*p+p)]
    phi=param_initials0[((c-1)*p+p+1):((c-1)*p+p+Int(q*(q-1)/2))]
    lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    lambda2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2))))]
    lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+ncov)]#    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)*p+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


    param_initialii = vcat(tau1,lambda1,phi,lambda2,bets)

    println(param_initialii,liksmix,"first iteration: done")
    return param_initialii
end


function EMparsMixl2inv(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
    local plimix::Array, T,oldlik
    n = size(X,1)
    er=0;
    # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
    it=0;
    dev=Inf; dev1=Inf; dev2=Inf;

    
    oldlik=-1000000.;
    oldlik2=1000000.;
    oldlikA=-1000000.;
    lik=0.
    

    noerror = true; devnostop = true;  sum_devs = 0;
    while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
        it = it+1;
        param=param01;

        T00 = 1
        plimix = E_stepMixl2inv(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        param01 = M_parsMixl2inv(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
        dt = time() - t
        if param01 == false
            return return(vcat(param01,vec(plimix),vec(T),0.,0.))
        end
        lik = E_stepMixl2inv(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);

        #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
        con = (lik-oldlik)/(oldlik-oldlik2);
        likA = oldlik + (lik-oldlik)/(1-con)
        old_dev = dev1
        dev1 = abs(likA-oldlikA)
        dev2 = abs(lik-oldlik)
        println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
        println(param01)
        dev = abs(lik-oldlik)/abs(oldlik)

        noerror = oldlik<lik

        oldlikA=likA
        oldlik2=oldlik
        oldlik=lik;
        pml_em_log("EMparsMixl2inv", it, lik, dev, dev1, abs(lik - oldlik2), tol, noerror)

        # sum_devs += dev > old_dev
        # devstop = sum_devs <= 3
    end
    
    # classificazioni dopo ultima massimizzazione
    plimix = E_stepMixl2inv(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

    if noerror == true && it<=maxit 
        return vcat(param01,vec(plimix),vec(T),lik,1.)
    else
        return vcat(param01,vec(plimix),vec(T),0.,0.)
    end
end





function E_stepMixl2inv(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
    if sum(param .== 0.) == size(param,1)
        return -Inf
    end
    trsi = (c-1)
    t1s = param[1:trsi*p]
    L1 = param[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    
    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )
   # interc2 = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p-1)+p)]
    betas = param[Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+ncov)]  #  betas = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)

    
    Sigms=ones(p,p,2);
    Sigms[:,:,1]=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sigms[:,:,2]=L2*L2'+ Diagonal(1 .-Diagonal((L2*L2')))
    
    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))

    tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]

    if estep == true
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));
    
        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]])
                    low=(tfa[:,pairsindex[:,idx]])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]])
                    low=(tfa[:,pairsindex[:,idx]])
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx] = t3.*pg[:,g];
            end
        end
    # % probabilities in each row must sum to 1 (application of Bayes law)
        return plimix
    elseif estep == false
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));

        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]])
                    low=(tfa[:,pairsindex[:,idx]])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]])
                    low=(tfa[:,pairsindex[:,idx]])
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx] = t3
            end
        end
        lik1=(log.(plimix));
        lik1=sum(lik1,dims=3);
        U=exp.(lik1)[:,:,1].*pg; 
        cllik=sum(log.(sum(U,dims=2)))
        return cllik
    end
end
    

    


function M_parsMixl2inv(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)

    trsi = (c-1)*p #attenzione p è già qui


    # PML - CFA & EFA
    theta_init = theta[1:(size(theta,1)-ncov)]#0.000001
    # # -1 loadings
    # lb=vec(vcat(fill(-5.1,trsi),fill(0.00001,p),fill(-1.,Int((q*(q-1))/2)),fill(1.,trs),fill(-1.,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-10,p)));
    # # 0. loadings
    lb=vec(vcat(fill(-5.1,trsi),fill(-1,p),fill(-1.,Int((q*(q-1))/2)),fill(-1,Int(p*K-K*(K-1)/2))));#1.0e-5
    ub=vec(vcat(fill(5.1,trsi),fill(1.,p),fill(1.,Int((q*(q-1))/2)),fill(1,Int(p*K-K*(K-1)/2))));

    par = zeros(Float64,size(lb))
    if it == 1
        par = optimize(Optim.only_fg!((F,G,x)->fgl2inv!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min ,Optim.Options(iterations=1000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0))
    elseif it > 1
        par = optimize(Optim.only_fg!((F,G,x)->fgl2inv!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06, time_limit = 300.0)) #20 min
    end 


    # par = optimize(x->pmllast(x,X,T,p,q,K,c,n,pairsindex,L_str,indices),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

    # par = zeros(Float64,size(lb))
    # if it == 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 300.0)) #20 min
    # elseif it > 1
    #     par = optimize(Optim.only_fg!((F,G,x)->fg!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),theta_init,LBFGS(),Optim.Options(time_limit = 1200.0)) #20 min
    # end 
    
    # mixture parameter
    betas = theta[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
    Ez1 = T[:,1,:]
    Ez0 = T[:,2,:]
    
    pairs = size(pairsindex,2)
    betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)

    
    # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
   
    # println(par.minimizer,pg[1])
    if !(Optim.converged(par)) && isnan(betas[1]) 
        println(betas[1])
        if start == true
            return zeros(Float64,size(theta,1))
        elseif start == false
            return Optim.converged(par)
        end
    end
    return vcat(par.minimizer,betasnew)#,pg[1])
end



function pmllastl2inv(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )
    
    
    
    Sig1=L1*Phi*L1' + Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2' + Diagonal(1 .-Diagonal((L2*L2')))

   
    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))

    tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
    
    plimix=zeros(Float64,n);
    for idx=1:size(pairsindex,2) #parallelizzabile
        # CFA
        Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
        Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
        tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
        Tt1 = T[:,1,idx]
        
        #EFA
        Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
        Refa=cov2cor(Sefa);
        Befa=hcat((tfa[:,pairsindex[:,idx]]),(tfa1[:,pairsindex[:,idx]]));
        tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
        Tt2 = T[:,2,idx]

        plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
    end
    f=-sum(plimix)
return f
end

function fgl2inv!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))])
    # prop = x[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]


    Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+ Diagonal(1 .-Diagonal((L2*L2')))


    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))

   

    if G !== nothing
        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)
    
        gradlambci = zeros(Float64,2,q,size(pairsindex,2))
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        # gradlamb2ci = zeros(Float64,2,K,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)

        idx = 0
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1
                rho1 = Sig1[i,j]
                cormat1 = ones(Float64,2,2); cormat1[1,2]=rho1; cormat1[2,1]=rho1
                sigs = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sigs)
                rho2 = cormat2[1,2]
                sigsij = diag(sigs)
                sigij = sqrt.(sigsij)
                sigij13 = sigij[1]^(3)
                sigij23 = sigij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]
                @inbounds for N in 1:n
                    xni = X[N,i]
                    xnj = X[N,j]
                    t1i = t1[i,xni]
                    t1i1 = t1[i,xni+1]
                    t1j = t1[j,xnj]
                    t1j1 = t1[j,xnj+1]
                    T1 = T[N,1,idx][1]
                    T2 = T[N,2,idx][1]
                    t2i = t1[i,xni]
                    t2i1 = t1[i,xni+1]
                    t2j = t1[j,xnj]
                    t2j1 = t1[j,xnj+1]
                    t2ist = t2i
                    t2i1st = t2i1
                    t2jst = t2j
                    t2j1st = t2j1
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                   
                    pijicfa = bvnP(t1i, t1i1, t1j, t1j1, rho1, false)
                    gradtaucfa[xni+1,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i1)*(cdf(Normal(0,1),(t1j1-rho1*t1i1)/sqrtrij1)-cdf(Normal(0,1),(t1j-rho1*t1i1)/sqrtrij1))) 
                    gradtaucfa[xni,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i)*(cdf(Normal(0,1),(t1j-rho1*t1i)/sqrtrij1)-cdf(Normal(0,1),(t1j1-rho1*t1i)/sqrtrij1))) 
                    gradtaucfa[xnj+1,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j1)*(cdf(Normal(0,1),(t1i1-rho1*t1j1)/sqrtrij1)-cdf(Normal(0,1),(t1i-rho1*t1j1)/sqrtrij1))) 
                    gradtaucfa[xnj,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j)*(cdf(Normal(0,1),(t1i-rho1*t1j)/sqrtrij1)-cdf(Normal(0,1),(t1i1-rho1*t1j)/sqrtrij1))) 
                    
    
                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))
    
                    if isnan(derivNormi1)
                        derivNormi1 = 0.
                    end
                    if isnan(derivNormi)
                        derivNormi = 0.
                    end
                    if isnan(derivNormj1)
                        derivNormj1 = 0.
                    end
                    if isnan(derivNormj)
                        derivNormj = 0.
                    end
    
                    
                    if !(t1i == -Inf && t1j == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else                  
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1 == Inf && t1j1 == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j1))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end
        
                    if !(t1i1 == Inf && t1j == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end
        
                    if !(t1i == -Inf && t1j1 == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j1))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    gradlambci[1:2,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)
                    gradphici[1:q,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*(L1[i,:]*L1[j,:]')
                    
                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sigij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sigij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sigij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sigij[2])
    
    
                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)
         # gradlamb2ci[1:2,1:K,idx] .+=  dpidrij*vcat(L2[j,:]',L2[i,:]')   .+  kappa .* vcat(derivNormi1 + derivNormi, derivNormj1 + derivNormj)      
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        
                        gradlsefaidx[k,idx,1] += dpidrij*(L2jk)#((L2j*Psi[i,i])/(sigij13*sigij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(L2ik)#((L2[i,:].*Psi[j,j])./(sigij23*sigij[1]))
    
                        # gradlamb2ci[k,xni+1,idx,1] -=  derivNormi1 
                        # gradlamb2ci[k,xnj+1,idx,2] -=  derivNormj1 
                        # gradlamb2ci[k,xni,idx,1]   -=  derivNormi  
                        # gradlamb2ci[k,xnj,idx,2]   -=  derivNormj  
                    end 
                          
                end
            end
        end
    
        vecpair = vec(pairsindex)
        liktaushared = zeros(p*(c-1));idpcfa = 0   
        @inbounds for j in 2:c
            @inbounds for i in 1:p
                idpcfa += 1
                liktaushared[idpcfa] = -sum(vec(vcat(gradtaucfa[j,:,1]',gradtaucfa[j,:,2]'))[vecpair.==i])-sum(vec(vcat(gradtauefa[j,:,1]',gradtauefa[j,:,2]'))[vecpair.==i])
            end
        end
        
                
        llamb1 = permutedims(gradlambci,(2,3,1)); llam1 = zeros(p); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:]; liklsidx  = zeros(p,K)
        @inbounds for i in 1:p
            llam1[i] = -sum(vec(llamb1[ids[i],:,:]')[vecpair.==i])
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end
        # llamb2 = permutedims(gradlamb2ci,(2,3,1)); llam2 = zeros(p);
        # liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        # @inbounds for i in 1:p
        #     llam1[i] = -sum(vec(llamb1[ids[i],:,:]')[vecpair.==i])
        #     llam2[i] = -sum(vec(llamb2)[vecpair.==i])
        # end
    
       
    
        liklstrue = likls.+liklsidx
        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]
    


        G .= vec(vcat(liktaushared,llam1,likphida,liklstrue[indices]))#likpsitrue[2:end]
      
    end
    if F !== nothing 
  
    
        plimix=zeros(Float64,n);
        tfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
        tfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
       
       
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
            Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
            Tt1 = T[:,1,idx]
            
            #EFA
            Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tfa[:,pairsindex[:,idx]])./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]])./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            Tt2 = T[:,2,idx]
            plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2;
            
        end
        f=-sum(plimix)
        # @info "Current objective: $f"
        return f
    end
end




########################################################################### Functions for no intercept CFA estimated ##############################################################


# function PMLcfaefa(X,p,q,K,c,tol,maxit, maxNumRestart::Int64)
#     local res::Vector{Float64}, err_0_count
   
    
#     iter_err = 1; code_error = 0.0; 
#     while code_error != 1.0
#         println("Restart: n. ", iter_err)

#         # Starting values
#         paramcfa = vcat(rand(Uniform(2.3,5),p),rand(Uniform(0.4,1),p),rand(Uniform(0.6,1.5),p),rand(Uniform(0.2,0.8),Int(q*(q-1)/2)))
#         paramcfaefa = vcat(paramcfa,vec(rand(Uniform(1.2,4),p,1)),rand(Uniform(0.1,1.),p*K),rand(Uniform(0.1,2),p),rand(Uniform(0,6),p))
#         param_init = vcat(paramcfaefa,rand(Uniform(0.82,0.96)))
        
#         # EM
#         res =  EMparsMix(param_init,tol,X,p,q,K,c,maxit)
        
#         # Controls
#         code_error = res[size(res,1)]
        
        
#         if code_error == 0.0  && iter_err < maxNumRestart
#             err_0_count += 1.
#             iter_err += 1.
#         end
#     end 
#     return vcat(res,err_0_count) #numero parametri + codice errore + numero restart
# end





# function EMparsMix(param01,tol,X,p,q,K,c,maxit)
#     local plimix::Array, T,lik,oldlik
#     n = size(X,1)
#     er=0;
#     pairsindex=mapreduce(permutedims,hcat,collect(combinations(1:p,2))');
#     # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
#     it=0;
#     dev=Inf;

#     A1 = diagm(ones(q));
#     L_str = kron(A1,ones(Int64(p/q)));
#     A2 = diagm(ones(K));
#     L2_str = kron(A2,ones(Int64(p)));
#     # T_init = T_init.*ones(Float64,n,2,size(pairsindex,2))
#     # param = M_parsMix(T_init,param_init,X,n,p,q,K,pairsindex,L_str,L2_str); #solo primo argomento di output
#     # T=T_init;
    
#     oldlik=-1000000.;
#     oldlik2=1000000.;
#     oldlikA=-1000000.;
#     lik=0.
#     noerror = true
#     while (dev>tol) && it <= maxit && noerror
        
#         it = it+1;
#         param=param01;

#         T00 = 1
#         plimix = E_stepMix(param,q,K,p,n,c,X,pairsindex,L_str,L2_str,true);
#         T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
#         param01 = M_parsMix(T,param,X,n,p,q,K,c,pairsindex,L_str,L2_str); #solo primo argomento di output
#         lik = E_stepMix(param01,q,K,p,n,c,X,pairsindex,L_str,L2_str,false,T);

#         #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
#         con = (lik-oldlik)/(oldlik-oldlik2);
#         likA = oldlik + (lik-oldlik)/(1-con)
#         dev = abs(likA-oldlikA)
#         # dev = abs(lik-oldlik)
#         println(lik," ",dev," ",oldlik," ",oldlik2," ",oldlik<lik )

#         noerror = oldlik<lik

#         oldlikA=likA
#         oldlik2=oldlik
#         oldlik=lik;

#     end
#     if  oldlik<lik 
#         return(vcat(param01,vec(T),lik,0.))
#     end

#     return vcat(param01,vec(T),lik,1.)
# end



# function E_stepMix(param,q,K,p,n,c,X,pairsindex,L_str,L2_str,estep,T=false)
#     trs = (c-1-2)
#     t1 = param[1:trs*p]
#     L1 = param[(trs*p+1):(trs*p+p)].*L_str 
#     Thet = diagm(param[(trs*p+p+1):(trs*p+p+p)])
#     phi = param[Int(trs*p+p+p+1):Int(trs*p+p+p+(q*(q-1))/2)]
#     Phi = vcat(hcat(vcat(0.,0.),reshape(vcat(phi[1:2],0.,phi[3]),2,2)'),zeros(Float64,3)')'.+vcat(hcat(vcat(1.,0.),reshape(vcat(phi[1:2],1,phi[3]),2,2)'),vcat(zeros(Float64,2),1.)')
   

#     t2 = param[Int(trs*p+p+p+(q*(q-1))/2+1):Int(trs*p+p+p+(q*(q-1))/2+trs*p)]
#     L2 = param[Int(trs*p+p+p+(q*(q-1))/2+trs*p+1):Int(trs*p+p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))] #.*idL #K
#     # L2 = hcat(Lam2[1:p],vcat(0.,Lam2[(p+1):(p*2)]))
#     # L2 = kron(ones(K),diagm(ones(p)))'*(Lam2 .* L2_str)
#     # L2 = hcat(Lam2[1:p],vcat(0.,Lam2[(p+1):p*K]))
#     if K == 3
#         L2 = hcat(L2[1:p],vcat(0.,L2[Int(p+1):(p*2-1)]),vcat(0.,0.,L2[Int(p*2):(p*3-3)]))
#     end

#     Psi = diagm(param[Int(trs*p+p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trs*p+p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+p)])#K
#     # prop = param[Int((trs*p+p+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+p)+1]#K
#     interc2 = param[Int(trs*p+p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+1):Int(trs*p+p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p)]
#     # nu = param[Int(((trs*p+p+p+(q*(q-1))/2+trs*p+p+p*K-(K*(K-1)/2))+p)+1):(Int((trs*p+p+p+(q*(q-1))/2+trs*p+p+p*K-(K*(K-1)/2))+p+K))]
#     prop = param[Int(trs*p+p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]

#     pg = vcat(prop,1. -prop)
#     #type==8
#     Sigms=ones(p,p,2);
#     Sigms[:,:,1]=L1*Phi*L1'+ Thet #+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
#     Sigms[:,:,2]=L2*L2'+Psi 
#     br1 = hcat(zeros(Float64,p),ones(Float64,p),reshape(t1,p,c-1-2))
#     vari1=hcat(fill(-Inf,p),br1,fill(Inf,p))
#     br2 = hcat(zeros(Float64,p),ones(Float64,p),reshape(t2,p,c-1-2))
#     vari2=hcat(fill(-Inf,p),br2,fill(Inf,p))
#     t11=[vari1[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     t12=[vari1[j,X[i,j]] for i in 1:n, j in 1:p]
#     t21=[vari2[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     t22=[vari2[j,X[i,j]] for i in 1:n, j in 1:p]
#     tuni = cat(t11,t21,dims=3)
#     tdui = cat(t12,t22,dims=3)
    
#     if estep == true
#         low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#         plimix=ones(n,2,size(pairsindex,2));
    
#         for idx=1:size(pairsindex,2) #parallelizzabile
#             for g=1:2
#                 s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                 if g==1
#                     # s2c=s2c+0.0001*Matrix(I,2,2)
#                     r2c=cov2cor(s2c);
#                     up=(tuni[:,pairsindex[:,idx],g])./sqrt.(diag(s2c)');
#                     low=(tdui[:,pairsindex[:,idx],g])./sqrt.(diag(s2c)');
#                 else 
#                     # s2c=s2c+0.0001*Matrix(I,2,2)
#                     r2c=cov2cor(s2c);
#                     up=(tuni[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                     low=(tdui[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 end
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix[:,g,idx]= t3*pg[g];
#             end
#         end
#     # % probabilities in each row must sum to 1 (application of Bayes law)
#         return plimix
#     elseif estep == false
#         low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#         plimix=ones(n,2,size(pairsindex,2));

#         for idx=1:size(pairsindex,2) #parallelizzabile
#             for g=1:2
#                 s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                 if g==1
#                     # s2c=s2c+0.0001*Matrix(I,2,2)
#                     r2c=cov2cor(s2c);
#                     up=(tuni[:,pairsindex[:,idx],g])./sqrt.(diag(s2c)');
#                     low=(tdui[:,pairsindex[:,idx],g])./sqrt.(diag(s2c)');
#                 else 
#                     # s2c=s2c+0.0001*Matrix(I,2,2)
#                     r2c=cov2cor(s2c);
#                     up=(tuni[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                     low=(tdui[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 end
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix[:,g,idx]= t3
#             end
#         end
#         lik1=(log.(plimix));
#         lik1=sum(lik1,dims=3);
#         U=exp.(lik1)[:,:,1]*diagm(pg); 
#         cllik=sum(log.(sum(U,dims=2)))
#         return cllik
#     end
# end
    


# function M_parsMix(T,theta,X,n,p,q,K,c,pairsindex,L_str,L2_str)
#     trsi = (c-1-2)*p
#     # pn = p*trs
    
#     # V=reshape(1:(p*q),p,q);
#     # Vu=(triu(V[1:q,1:q],1))
#     # idL=ones(p*q); idL[Vu[Vu.>0]].=0.0

    
#     lb=vec(vcat(fill(-1,trsi),fill(0.00001,p),fill(0.0001,p),fill(-1.,Int((q*(q-1))/2)),fill(-1.,trsi),fill(0.000001,Int(p*K-K*(K-1)/2)),fill(0.0001,p),fill(0.1,p)));
#     ub=vec(vcat(fill(7.,trsi),fill(1.,p),fill(3.,p),fill(1.,Int((q*(q-1))/2)),fill(7.,trsi),fill(1,Int(p*K-K*(K-1)/2)),fill(2.,p),fill(8.,p)));
#     # func = TwiceDifferentiable(param -> pmlmix(param,X,T,p,q,K,c,n,pairsindex,L_str,L2_str),theta[1:(size(theta,1)-1)]; autodiff=:forward);
#     # dfc = TwiceDifferentiableConstraints(lb, ub)
#     # par = optimize(func,dfc,theta[1:(size(theta,1)-1)],IPNewton(),Optim.Options(iterations=10000,f_reltol=1.0000e-06))
    
#     par = optimize(param ->pmlmix(param,X,T,p,q,K,c,n,pairsindex,L_str,L2_str),lb,ub,theta[1:(size(theta,1)-1)],Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

#     pg = 0.
#     if size(T,3) > 1
#         pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
#     else
#         pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
#     end
#     return vcat(par.minimizer,pg[1])
# end





# function pmlmix(param,X,T,p,q,K,c,n,pairsindex,L_str,L2_str)

#     trs = (c-1-2)
#     t1 = param[1:trs*p]
#     L1 = param[(trs*p+1):(trs*p+p)].*L_str 
#     # L1 = hcat(vcat(1.,L1s[1:2],zeros(6)),vcat(zeros(3),1.,L1s[3:4],zeros(3)),vcat(zeros(6),1.,L1s[5:6]))
#     Thet = diagm(param[(trs*p+p+1):(trs*p+p+p)])
#     phi = param[Int(trs*p+p+p+1):Int(trs*p+p+p+(q*(q-1))/2)]
#     Phi = vcat(hcat(vcat(0.,0.),reshape(vcat(phi[1:2],0.,phi[3]),2,2)'),zeros(Float64,3)')'.+vcat(hcat(vcat(1.,0.),reshape(vcat(phi[1:2],1,phi[3]),2,2)'),vcat(zeros(Float64,2),1.)')
#     # Phi = hcat(vcat(1.,phi[1:4]),vcat(phi[1],1.,phi[5:7]),vcat(phi[vcat(2,5)],1.,phi[8:9]),vcat(phi[vcat(3,6,8)],1.,phi[10]),vcat(phi[vcat(4,7,9,10)],1.))
    

#     t2 = param[Int(trs*p+p+p+(q*(q-1))/2+1):Int(trs*p+p+p+(q*(q-1))/2+trs*p)]
#     L2 = param[Int(trs*p+p+p+(q*(q-1))/2+trs*p+1):Int(trs*p+p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))] #.*idL #K-1
#     # L2 = hcat(Lam2[1:p],vcat(0.,Lam2[(p+1):(p*2-1)]))
#     # L2 = kron(ones(K),diagm(ones(p)))'*(Lam2 .* L2_str)
#     if K == 3
#         L2 = hcat(L2[1:p],vcat(0.,L2[Int(p+1):(p*2-1)]),vcat(0.,0.,L2[Int(p*2):(p*3-3)]))
#     end
#     Psi = diagm(param[Int(trs*p+p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+1):Int((trs*p+p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2)))+p)])#K
#     interc2 = param[Int(trs*p+p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+1):Int(trs*p+p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+p+p)]
   



#     Sig1=L1*Phi*L1'+ Thet#Diagonal(1 .-Diagonal((L1*Phi*L1')))
#     Sig2=L2*L2'+Psi

#     plimix1=zeros(n);plimix2=zeros(n);
   
#     br1 = hcat(zeros(Float64,p),ones(Float64,p),reshape(t1,p,c-1-2))
#     vari1=hcat(fill(-Inf,p),br1,fill(Inf,p))
#     br2 = hcat(zeros(Float64,p),ones(Float64,p),reshape(t2,p,c-1-2))
#     vari2=hcat(fill(-Inf,p),br2,fill(Inf,p))
    
#     t11=[vari1[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     t12=[vari1[j,X[i,j]] for i in 1:n, j in 1:p]
#     t21=[vari2[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     t22=[vari2[j,X[i,j]] for i in 1:n, j in 1:p]
   
   
#     tuni = cat(t11,t21,dims=3)
#     tdui = cat(t12,t22,dims=3)
#     low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#     plimix=ones(n,2,size(pairsindex,2));
#     for idx=1:size(pairsindex,2) #parallelizzabile
#         for g=1:2
#             if g==1
#                 s2c=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
#                 # s2c=s2c+0.0001*Matrix(I,2,2)
#                 r2c=cov2cor(s2c);
            
#                 up=(tuni[:,pairsindex[:,idx],g])./sqrt.(diag(s2c)');
#                 low=(tdui[:,pairsindex[:,idx],g])./sqrt.(diag(s2c)');
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix1 += log.(t3).*T[:,1,idx];
#             else 
#                 s2c=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
#                 # s2c=s2c+0.0001*Matrix(I,2,2)
#                 r2c=cov2cor(s2c);
#                 up=(tuni[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 low=(tdui[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix2 += log.(t3).*T[:,2,idx];
#             end
            
#         end
#     end
#     f=-sum(plimix1+plimix2)
#     # @info "Current objective: $f"
#     return f
# end
  




########################################################################### Functions for intercept CFA estimated ##############################################################


# function EMparsMix(param01,tol,X,p,q,K,c,maxit)
#     local plimix::Array, T,lik,oldlik
#     n = size(X,1)
#     er=0;
#     pairsindex=mapreduce(permutedims,hcat,collect(combinations(1:p,2))');
#     # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
#     it=0;
#     dev=Inf;

#     A1 = diagm(ones(q));
#     L_str = kron(A1,ones(Int64(p/q)));
#     A2 = diagm(ones(K));
#     L2_str = kron(A2,ones(Int64(p)));
#     # T_init = T_init.*ones(Float64,n,2,size(pairsindex,2))
#     # param = M_parsMix(T_init,param_init,X,n,p,q,K,pairsindex,L_str,L2_str); #solo primo argomento di output
#     # T=T_init;
    
#     oldlik=-1000000.;
#     oldlik2=1000000.;
#     oldlikA=-1000000.;
#     lik=0.
#     noerror = true
#     while (dev>tol) && it <= maxit && noerror
        
#         it = it+1;
#         param=param01;

#         T00 = 1
#         plimix = E_stepMix(param,q,K,p,n,c,X,pairsindex,L_str,L2_str,true);
#         T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
#         param01 = M_parsMix(T,param,X,n,p,q,K,c,pairsindex,L_str,L2_str); #solo primo argomento di output
#         lik = E_stepMix(param01,q,K,p,n,c,X,pairsindex,L_str,L2_str,false,T);

#         #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
#         con = (lik-oldlik)/(oldlik-oldlik2);
#         likA = oldlik + (lik-oldlik)/(1-con)
#         dev = abs(likA-oldlikA)
#         # dev = abs(lik-oldlik)
#         println(lik," ",dev," ",oldlik," ",oldlik2," ",oldlik<lik )

#         noerror = oldlik<lik

#         oldlikA=likA
#         oldlik2=oldlik
#         oldlik=lik;

#     end
#     if  oldlik<lik 
#         return(vcat(param01,vec(T),lik,0.))
#     end

#     return vcat(param01,vec(T),lik,1.)
# end



# function E_stepMix(param,q,K,p,n,c,X,pairsindex,L_str,L2_str,estep,T=false)
#     trs = (c-1-2)
#     t1 = param[1:trs*p]
#     L1 = param[(trs*p+1):(trs*p+p)].*L_str 
#     Thet = diagm(param[(trs*p+p+1):(trs*p+p+p)])
#     phi = param[Int(trs*p+p+p+1):Int(trs*p+p+p+(q*(q-1))/2)]
#     Phi = vcat(hcat(vcat(0.,0.),reshape(vcat(phi[1:2],0.,phi[3]),2,2)'),zeros(Float64,3)')'.+vcat(hcat(vcat(1.,0.),reshape(vcat(phi[1:2],1,phi[3]),2,2)'),vcat(zeros(Float64,2),1.)')
#     interc1 = param[Int(trs*p+p+p+(q*(q-1))/2+1):Int(trs*p+p+p+(q*(q-1))/2+p)]

#     t2 = param[Int(trs*p+p+p+(q*(q-1))/2+p+1):Int(trs*p+p+p+(q*(q-1))/2+p+trs*p)]
#     L2 = param[Int(trs*p+p+p+(q*(q-1))/2+p+trs*p+1):Int(trs*p+p+p+(q*(q-1))/2+p+trs*p+(p*K-(K*(K-1)/2)))] #.*idL #K
#     # L2 = hcat(Lam2[1:p],vcat(0.,Lam2[(p+1):(p*2)]))
#     # L2 = kron(ones(K),diagm(ones(p)))'*(Lam2 .* L2_str)
#     # L2 = hcat(Lam2[1:p],vcat(0.,Lam2[(p+1):p*K]))
#     if K == 3
#         L2 = hcat(L2[1:p],vcat(0.,L2[Int(p+1):(p*2-1)]),vcat(0.,0.,L2[Int(p*2):(p*3-3)]))
#     end

#     Psi = diagm(param[Int(trs*p+p+p+(q*(q-1))/2+p+trs*p+(p*K-(K*(K-1)/2))+1):Int((trs*p+p+p+(q*(q-1))/2+p+trs*p+(p*K-(K*(K-1)/2)))+p)])#K
#     # prop = param[Int((trs*p+p+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2)))+p)+1]#K
#     interc2 = param[Int(trs*p+p+p+(q*(q-1))/2+p+trs*p+(p*K-(K*(K-1)/2))+p+1):Int(trs*p+p+p+(q*(q-1))/2+p+trs*p+(p*K-(K*(K-1)/2))+p+p)]
#     # nu = param[Int(((trs*p+p+p+(q*(q-1))/2+p+trs*p+p+p*K-(K*(K-1)/2))+p)+1):(Int((trs*p+p+p+(q*(q-1))/2+p+trs*p+p+p*K-(K*(K-1)/2))+p+K))]
#     prop = param[Int(trs*p+p+p+(q*(q-1))/2+p+trs*p+(p*K-(K*(K-1)/2))+p+p+1)]

#     pg = vcat(prop,1. -prop)
#     #type==8
#     Sigms=ones(p,p,2);
#     Sigms[:,:,1]=L1*Phi*L1'+ Thet #+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
#     Sigms[:,:,2]=L2*L2'+Psi 
#     br1 = hcat(zeros(Float64,p),ones(Float64,p),reshape(t1,p,c-1-2))
#     vari1=hcat(fill(-Inf,p),br1,fill(Inf,p))
#     br2 = hcat(zeros(Float64,p),ones(Float64,p),reshape(t2,p,c-1-2))
#     vari2=hcat(fill(-Inf,p),br2,fill(Inf,p))
#     t11=[vari1[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     t12=[vari1[j,X[i,j]] for i in 1:n, j in 1:p]
#     t21=[vari2[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     t22=[vari2[j,X[i,j]] for i in 1:n, j in 1:p]
#     tuni = cat(t11,t21,dims=3)
#     tdui = cat(t12,t22,dims=3)
    
#     if estep == true
#         low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#         plimix=ones(n,2,size(pairsindex,2));
    
#         for idx=1:size(pairsindex,2) #parallelizzabile
#             for g=1:2
#                 s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                 if g==1
#                     # s2c=s2c+0.0001*Matrix(I,2,2)
#                     r2c=cov2cor(s2c);
#                     up=(tuni[:,pairsindex[:,idx],g].-(interc1)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
#                     low=(tdui[:,pairsindex[:,idx],g].-(interc1)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
#                 else 
#                     # s2c=s2c+0.0001*Matrix(I,2,2)
#                     r2c=cov2cor(s2c);
#                     up=(tuni[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                     low=(tdui[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 end
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix[:,g,idx]= t3*pg[g];
#             end
#         end
#     # % probabilities in each row must sum to 1 (application of Bayes law)
#         return plimix
#     elseif estep == false
#         low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#         plimix=ones(n,2,size(pairsindex,2));

#         for idx=1:size(pairsindex,2) #parallelizzabile
#             for g=1:2
#                 s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
#                 if g==1
#                     # s2c=s2c+0.0001*Matrix(I,2,2)
#                     r2c=cov2cor(s2c);
#                     up=(tuni[:,pairsindex[:,idx],g].-(interc1)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
#                     low=(tdui[:,pairsindex[:,idx],g].-(interc1)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
#                 else 
#                     # s2c=s2c+0.0001*Matrix(I,2,2)
#                     r2c=cov2cor(s2c);
#                     up=(tuni[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                     low=(tdui[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 end
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix[:,g,idx]= t3
#             end
#         end
#         lik1=(log.(plimix));
#         lik1=sum(lik1,dims=3);
#         U=exp.(lik1)[:,:,1]*diagm(pg); 
#         cllik=sum(log.(sum(U,dims=2)))
#         return cllik
#     end
# end
    


# function M_parsMix(T,theta,X,n,p,q,K,c,pairsindex,L_str,L2_str)
#     trsi = (c-1-2)*p
#     # pn = p*trs
    
#     # V=reshape(1:(p*q),p,q);
#     # Vu=(triu(V[1:q,1:q],1))
#     # idL=ones(p*q); idL[Vu[Vu.>0]].=0.0

    
#     lb=vec(vcat(fill(-1,trsi),fill(0.00001,p),fill(0.0001,p),fill(-1.,Int((q*(q-1))/2)),fill(-1.,p),fill(-1.,trsi),fill(0.000001,Int(p*K-K*(K-1)/2)),fill(0.0001,p),fill(0.1,p)));
#     ub=vec(vcat(fill(7.,trsi),fill(1.,p),fill(3.,p),fill(1.,Int((q*(q-1))/2)),fill(8,p),fill(7.,trsi),fill(1,Int(p*K-K*(K-1)/2)),fill(2.,p),fill(8.,p)));
#     # func = TwiceDifferentiable(param -> pmlmix(param,X,T,p,q,K,c,n,pairsindex,L_str,L2_str),theta[1:(size(theta,1)-1)]; autodiff=:forward);
#     # dfc = TwiceDifferentiableConstraints(lb, ub)
#     # par = optimize(func,dfc,theta[1:(size(theta,1)-1)],IPNewton(),Optim.Options(iterations=10000,f_reltol=1.0000e-06))
    
#     par = optimize(param ->pmlmix(param,X,T,p,q,K,c,n,pairsindex,L_str,L2_str),lb,ub,theta[1:(size(theta,1)-1)],Fminbox(LBFGS()),Optim.Options(iterations=10000,f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06);autodiff=:forward)

#     pg = 0.
#     if size(T,3) > 1
#         pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
#     else
#         pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
#     end
#     return vcat(par.minimizer,pg[1])
# end





# function pmlmix(param,X,T,p,q,K,c,n,pairsindex,L_str,L2_str)

#     trs = (c-1-2)
#     t1 = param[1:trs*p]
#     L1 = param[(trs*p+1):(trs*p+p)].*L_str 
#     # L1 = hcat(vcat(1.,L1s[1:2],zeros(6)),vcat(zeros(3),1.,L1s[3:4],zeros(3)),vcat(zeros(6),1.,L1s[5:6]))
#     Thet = diagm(param[(trs*p+p+1):(trs*p+p+p)])
#     phi = param[Int(trs*p+p+p+1):Int(trs*p+p+p+(q*(q-1))/2)]
#     Phi = vcat(hcat(vcat(0.,0.),reshape(vcat(phi[1:2],0.,phi[3]),2,2)'),zeros(Float64,3)')'.+vcat(hcat(vcat(1.,0.),reshape(vcat(phi[1:2],1,phi[3]),2,2)'),vcat(zeros(Float64,2),1.)')
#     # Phi = hcat(vcat(1.,phi[1:4]),vcat(phi[1],1.,phi[5:7]),vcat(phi[vcat(2,5)],1.,phi[8:9]),vcat(phi[vcat(3,6,8)],1.,phi[10]),vcat(phi[vcat(4,7,9,10)],1.))
#     interc1 = param[Int(trs*p+p+p+(q*(q-1))/2+1):Int(trs*p+p+p+(q*(q-1))/2+p)]

#     t2 = param[Int(trs*p+p+p+(q*(q-1))/2+p+1):Int(trs*p+p+p+(q*(q-1))/2+p+trs*p)]
#     L2 = param[Int(trs*p+p+p+(q*(q-1))/2+p+trs*p+1):Int(trs*p+p+p+(q*(q-1))/2+p+trs*p+(p*K-(K*(K-1)/2)))] #.*idL #K-1
#     # L2 = hcat(Lam2[1:p],vcat(0.,Lam2[(p+1):(p*2-1)]))
#     # L2 = kron(ones(K),diagm(ones(p)))'*(Lam2 .* L2_str)
#     if K == 3
#         L2 = hcat(L2[1:p],vcat(0.,L2[Int(p+1):(p*2-1)]),vcat(0.,0.,L2[Int(p*2):(p*3-3)]))
#     end
#     Psi = diagm(param[Int(trs*p+p+p+(q*(q-1))/2+p+trs*p+(p*K-(K*(K-1)/2))+1):Int((trs*p+p+p+(q*(q-1))/2+p+trs*p+(p*K-(K*(K-1)/2)))+p)])#K
#     interc2 = param[Int(trs*p+p+p+(q*(q-1))/2+p+trs*p+(p*K-(K*(K-1)/2))+p+1):Int(trs*p+p+p+(q*(q-1))/2+p+trs*p+(p*K-(K*(K-1)/2))+p+p)]
   



#     Sig1=L1*Phi*L1'+ Thet#Diagonal(1 .-Diagonal((L1*Phi*L1')))
#     Sig2=L2*L2'+Psi

#     plimix1=zeros(n);plimix2=zeros(n);
   
#     br1 = hcat(zeros(Float64,p),ones(Float64,p),reshape(t1,p,c-1-2))
#     vari1=hcat(fill(-Inf,p),br1,fill(Inf,p))
#     br2 = hcat(zeros(Float64,p),ones(Float64,p),reshape(t2,p,c-1-2))
#     vari2=hcat(fill(-Inf,p),br2,fill(Inf,p))
    
#     t11=[vari1[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     t12=[vari1[j,X[i,j]] for i in 1:n, j in 1:p]
#     t21=[vari2[j,X[i,j]+1] for i in 1:n, j in 1:p]
#     t22=[vari2[j,X[i,j]] for i in 1:n, j in 1:p]
   
   
#     tuni = cat(t11,t21,dims=3)
#     tdui = cat(t12,t22,dims=3)
#     low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
#     plimix=ones(n,2,size(pairsindex,2));
#     for idx=1:size(pairsindex,2) #parallelizzabile
#         for g=1:2
#             if g==1
#                 s2c=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
#                 # s2c=s2c+0.0001*Matrix(I,2,2)
#                 r2c=cov2cor(s2c);
            
#                 up=(tuni[:,pairsindex[:,idx],g].-(interc1)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
#                 low=(tdui[:,pairsindex[:,idx],g].-(interc1)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix1 += log.(t3).*T[:,1,idx];
#             else 
#                 s2c=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
#                 # s2c=s2c+0.0001*Matrix(I,2,2)
#                 r2c=cov2cor(s2c);
#                 up=(tuni[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 low=(tdui[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
#                 B=hcat(low,up);
#                 t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
#                 plimix2 += log.(t3).*T[:,2,idx];
#             end
            
#         end
#     end
#     f=-sum(plimix1+plimix2)
#     # @info "Current objective: $f"
#     return f
# end
  











########################################################################## Functions for CDF bivariate Normal ###################################################################
function bvnP( xl, xu, yl, yu, r, vect=true)
    # %BVN
    # %  A function for computing bivariate normal probabilities.
    # %  bvn calculates the probability that 
    # %    xl < x < xu and yl < y < yu, 
    # %  with correlation coefficient r.
    # %   p = bvn( xl, xu, yl, yu, r )
    # %
    
    # %   Author
    # %       Alan Genz, Department of Mathematics
    # %       Washington State University, Pullman, Wa 99164-3113
    # %       Email : alangenz@wsu.edu
    # %
    if vect == true
        pip=ones(length(xl));
        for i=1:length(xl)
        pip[i] = bvnlr(xl[i],yl[i],r) - bvnlr(xu[i],yl[i],r) - bvnlr(xl[i],yu[i],r) + bvnlr(xu[i],yu[i],r); 
        pip[i] = max( 2.2204e-16, min( pip[i], 1 ) );
        end

        return pip

    elseif vect == false
        p=ones(length(xl));
        p = bvnlr(xl,yl,r) - bvnlr(xu,yl,r) - bvnlr(xl,yu,r) + bvnlr(xu,yu,r); 
        p = max( 2.2204e-16, min( p, 1 ) )
        
        return p
    end 
    
end




function bvnlr(dh, dk, r)
    # %BVNLR
    # %  A function for computing bivariate normal probabilities.
    # %  bvnlr calculates the probability that x < dh and y < dk.
    # %    parameters
    # %      dh 1st upper integration limit
    # %      dk 2nd upper integration limit
    # %      r   correlation coefficient
    # %  Example:
    # %    p = bvnlr( 3, 1, .35 )
    # %
    
    # %   Author
    # %       Alan Genz
    # %       Department of Mathematics
    # %       Washington State University
    # %       Pullman, Wa 99164-3113
    # %       Email : alangenz@wsu.edu
    # %   This function is based on the method described by
    # %        Drezner, Z and G.O. Wesolowsky, (1989),
    # %        On the computation of the bivariate normal inegral,
    # %        Journal of Statist. Comput. Simul. 35, pp. 101-107,
    # %    with major modifications for double precision, for |r| close to 1,
    # %    and for matlab by Alan Genz - last modifications 7/98.
    # %
    dh=-dh;
    dk=-dk;
    if dh ==  Inf || dk ==  Inf
         p = 0;
    elseif dh == -Inf
        if dk == -Inf
            p = 1;
        else
            p = phid(-dk); 
        end
    elseif dk == -Inf
        p = phid(-dh);
    else
        if abs(r) < 0.3
            # %       Gauss Legendre points and weights, n =  6
            w = [0.1713244923791705 0.3607615730481384 0.4679139345726904]';
            x = [0.9324695142031522 0.6612093864662647 0.2386191860831970]';
        elseif abs(r) < 0.75
            # %       Gauss Legendre points and weights, n = 12
            w = [.04717533638651177 0.1069393259953183 0.1600783285433464 0.2031674267230659 0.2334925365383547 0.2491470458134029]';
            x = [0.9815606342467191 0.9041172563704750 0.7699026741943050 0.5873179542866171 0.3678314989981802 0.1252334085114692]';
        else
            # %       Gauss Legendre points and weights, n = 20
            w = [.01761400713915212 .04060142980038694 .06267204833410906 .08327674157670475 0.1019301198172404 0.1181945319615184 0.1316886384491766 0.1420961093183821 0.1491729864726037 0.1527533871307259]';
            x = [0.9931285991850949 0.9639719272779138 0.9122344282513259 0.8391169718222188 0.7463319064601508 0.6360536807265150 0.5108670019508271 0.3737060887154196 0.2277858511416451 0.07652652113349733]';
        end
        h = dh; k = dk; hk = h*k; bvn = 0;
        if abs(r) < 0.925 
            hs = ( h*h + k*k )/2; 
            asr = asin(r);
            sn = sin.( asr.*( 1 .- x)./2 );
            bvn = w'*exp.( ( sn*hk .- hs )./( 1 .- sn.*sn ) );
            sn = sin.( asr.*( 1 .+ x )./2 );
            bvn = bvn + w'*exp.( ( sn*hk .- hs )./( 1 .- sn.*sn ) );
            bvn = bvn.*asr/( 4*pi );
            bvn = bvn .+ phid(-h)*phid(-k);
        else
            twopi = 2*pi;
            if r < 0
                k = -k; 
                hk = -hk;
            end
            if abs(r) < 1
                as = ( 1 - r )*( 1 + r ); a = sqrt(as); 
                bs = ( h - k )^2;
                c = ( 4 - hk )/8 ; d = ( 12 - hk )/16; 
                asr = -( bs/as + hk )/2;
                if asr > -100
                    bvn = a*exp(asr)*( 1 - c*(bs-as)*(1-d*bs/5)/3 + c*d*as*as/5 );
                end
                if hk > -100
                    b = sqrt(bs); sp = sqrt(twopi)*phid(-b/a);
                    bvn = bvn[1,1] - exp(-hk/2)*sp*b*( 1 - c*bs*( 1 - d*bs/5 )/3 );
                end
                a = a/2;
                xs = ( a .- a*x ).^2;
                asr = -( bs./xs .+ hk )./2;
                ind = asr .> -100;
                xs=xs[ind];
                rs = sqrt.( 1 .- xs );
                sp = ( 1 .+ c*xs.*( 1 .+ d*xs ) );
                ep = exp.( -hk.*( 1 .- rs )./( 2 .*( 1 .+ rs ) ) )./rs;
                bvn = bvn + a*w[ind]'*(exp.(asr[ind]).*( ep - sp ));
                xs = ( a .+ a*x ).^2;
                asr = -( bs./xs .+ hk )./2;
                ind = asr .> -100;
                xs=xs[ind];
                rs = sqrt.( 1 .- xs );
                sp = ( 1 .+ c*xs.*( 1 .+ d*xs ) );
                ep = exp.( -hk*( 1 .- rs )./( 2 .*( 1 .+ rs ) ) )./rs;
                bvn = bvn + a*w[ind]'*(exp.(asr[ind]).*( ep - sp ));
                bvn = -bvn/twopi;
        end
        if r > 0
            bvn =  bvn .+ phid( -max( h, k ) ); 
        end
        if r < 0
            bvn = -bvn .+ max( 0, phid(-h)-phid(-k) ); 
        end
    end
    p = max( 0, min( 1, bvn[1,1] ) );
end
end
    

function phid(z)
    p = erfc( -z/sqrt(2) )/2;# % Normal cdf
end






function pmlcfa(param::AbstractVector,p::Int64,q::Int64,c::Int64,X::Matrix{Int64},L_str)
    t = reshape(param[1:((c-1)*p)],p,c-1)
    t = hcat(fill(-Inf,p),t,fill(Inf,p))
    Lambda = param[((c-1)*p+1):((c-1)*p+p)].*L_str #reshape(param[((c-1)*p+1):((c-1)*p+p*q)],p,q)
    Phi = reconstruct_correlation_matrix(param[Int((c-1)*p+p+1):Int((c-1)*p+p+(q*(q-1))/2)])
   
    nlp = 0.
    for i in 1:(p-1)
        for j in (i+1):p
            rho = Lambda[i,:]'*Phi*Lambda[j,:]  
            nij = Matrix{Float64}(contingency_table(X[:,i], X[:,j], c))
            for ci in 1:(c)
                for cj in 1:(c)
                    # println(t[ci]," ",t[ci+1]," ",t[cj]," ",t[cj+1]," ",nij[ci,cj])
                    nlp += nij[ci,cj]*log(bvnP(t[i,ci], t[i,ci+1], t[j,cj], t[j,cj+1], rho,false))
                
                end
            end
        end
    end 
    f = -sum(nlp)
    # @info "Current objective: $f"
    return f
end



# function pmlcfa(param::AbstractVector,p::Int64,q::Int64,c::Int64,X::Matrix{Int64},L_str::Matrix{Float64})
#     t = reshape(param[1:((c-1)*p)],p,c-1)
#     t = hcat(fill(-Inf,p),t,fill(Inf,p))
#     Lambda = param[((c-1)*p+1):((c-1)*p+p)].*L_str #reshape(param[((c-1)*p+1):((c-1)*p+p*q)],p,q)
#     Phi = reconstruct_correlation_matrix(param[Int((c-1)*p+p+1):Int((c-1)*p+p+(q*(q-1))/2)])
   
#     N = size(X,1)
#     nlp = 0.
#     Sig1=Lambda*Phi*Lambda'+ Diagonal(1 .-Diagonal((Lambda*Phi*Lambda')))
#     tcfa1=[t[j,X[i,j]+1] for i in 1:N, j in 1:p]
#     tcfa=[t[j,X[i,j]] for i in 1:N, j in 1:p]
   
#     for n in 1:N
#         R"""
#         probs = mvtnorm::pmvnorm($tcfa[$n,], $tcfa1[$n,], rep(0,$p), $Sig1)
#         """
#         probss = rcopy(R"probs")
#         nlp += log(probss)
#     end 
#     f = -sum(nlp)
#     # @info "Current objective: $f"
#     return f
# end
function count_and_reduce_rows(arr::Matrix{T}) where T
    # Convert each row to a tuple for easier comparison
    row_tuples = [tuple(row...) for row in eachrow(arr)]
    
    # Create a dictionary to count occurrences of each row
    row_counts = Dict()
    for row in row_tuples
        row_counts[row] = get(row_counts, row, 0) + 1
    end
    
    # Extract unique rows and their counts
    unique_rows = hcat([collect(v) for v in keys(row_counts)]...)
    counts = collect(values(row_counts))
    
    # Convert counts to a column vector
    counts_col = reshape(counts, :, 1)
    
    # Concatenate the unique rows with counts column
    return hcat(unique_rows', counts_col)
end


# function pmlcfa(param::AbstractVector, p::Int64, q::Int64, c::Int64, X::Matrix{Int64}, L_str::Matrix{Float64})
#     t = reshape(param[1:((c-1)*p)], p, c-1)
#     t = hcat(fill(-Inf, p), t, fill(Inf, p))
#     Lambda = param[((c-1)*p+1):((c-1)*p+p)] .* L_str
#     Phi = reconstruct_correlation_matrix(param[Int((c-1)*p+p+1):Int((c-1)*p+p+(q*(q-1))/2)])

#     Xcom = count_and_reduce_rows(X)
#     R = size(Xcom, 1)
#     countss = Xcom[:,p+1]
#     Xred = Xcom[:,1:p]
#     nlp = 0.0
#     Sig1 = Lambda * Phi * Lambda' + Diagonal(1.0 .- Diagonal(Lambda * Phi * Lambda'))
#     tcfa1 = [t[j, Xred[i,j]+1] for i in 1:R, j in 1:p]
#     tcfa = [t[j, Xred[i,j]] for i in 1:R, j in 1:p]

#     for r in 1:R
#         tcfa_r = tcfa[r, :]
#         tcfa1_r = tcfa1[r, :]
#         Sig1_r = Sig1
#         # tcfa_r = ForwardDiff.value.(tcfa[r, :])
#         # tcfa1_r = ForwardDiff.value.(tcfa1[r, :])
#         # Sig1_r = ForwardDiff.value.(Sig1)

#         probss = mvnormcdf(Sig1, tcfa[r,:],tcfa1[r,:],m=10000*p)
#         # @rput tcfa_r tcfa1_r Sig1_r p

#         # R"""
#         # probs <- mvtnorm::pmvnorm(tcfa_r, tcfa1_r, rep(0,p), Sig1_r, algorithm = "GenzBretz")#, algorithm="Miwa")
#         # """
#         # probss = rcopy(R"probs")
#         nlp += countss[r]*log(probss[1])
#     end 

#     f = -nlp
#     @info "Current objective: $f"
#     return f
# end

function pmlefaint(param::AbstractVector,p::Int64,K::Int64,c::Int64,X::Matrix{Int64},indices::Vector{Int64})
    t = reshape(param[1:(c-1-2)*p],p,c-1-2)
    ts = hcat(fill(-Inf,p),zeros(Float64,p),ones(Float64,p),t,fill(Inf,p))
    L2 = reconstruct_L2_matrix(p, K, indices,param[Int((c-1-2)*p+1):Int((c-1-2)*p+(p*K-(K*(K-1)/2)))])
    Psi = diagm(param[Int((c-1-2)*p+(p*K-(K*(K-1)/2))+1):Int((c-1-2)*p+(p*K-(K*(K-1)/2))+p)])
    interc = param[Int((c-1-2)*p+(p*K-(K*(K-1)/2))+p+1):Int((c-1-2)*p+(p*K-(K*(K-1)/2))+p+p)]
    
    Sig = L2*L2' + Psi
    nlp = 0.
    for i in 1:(p-1)
        for j in (i+1):p
            s2c=Sig[i,j]
            sigs = hcat(vcat(Sig[i,i],Sig[i,j]),vcat(Sig[j,i],Sig[j,j]))
            rhos=cov2cor(sigs);
            rho=rhos[1,2];
            nu = interc[vcat(i,j)]
            nij = Matrix{Float64}(contingency_table(X[:,i], X[:,j], c))
            sigmaij = sqrt.(diag(sigs))
            for ci in 1:(c)
                for cj in 1:(c)
                    nlp += nij[ci,cj]*log(bvnP((ts[i,ci]-nu[1])/sigmaij[1], (ts[i,ci+1]-nu[1])/sigmaij[1], (ts[j,cj]-nu[2])/sigmaij[2], (ts[j,cj+1]-nu[2])/sigmaij[2], rho,false))
                                                               
                end
            end
        end
    end 
    f = -sum(nlp)
    # @info "Current objective: $f"
    return f
end



function pmlefalm(param::AbstractVector,p::Int64,K::Int64,c::Int64,X::Matrix{Int64},indices::Vector{Int64})
    t = reshape(param[1:(c-1-2)*p],p,c-1-2)
    ts = hcat(fill(-Inf,p),zeros(Float64,p),ones(Float64,p),t,fill(Inf,p))
    L2 = reconstruct_L2_matrix(p, K, indices,param[Int((c-1-2)*p+1):Int((c-1-2)*p+(p*K-(K*(K-1)/2)))])
    Psi = diagm(param[Int((c-1-2)*p+(p*K-(K*(K-1)/2))+1):Int((c-1-2)*p+(p*K-(K*(K-1)/2))+p)])
    kap = param[Int((c-1-2)*p+(p*K-(K*(K-1)/2))+p+1):Int((c-1-2)*p+(p*K-(K*(K-1)/2))+p+K)]
    
    mu = L2*kap
    Sig = L2*L2' + Psi
    nlp = 0.
    for i in 1:(p-1)
        for j in (i+1):p
            s2c=Sig[i,j]
            sigs = hcat(vcat(Sig[i,i],Sig[i,j]),vcat(Sig[j,i],Sig[j,j]))
            rhos=cov2cor(sigs);
            rho=rhos[1,2];
            nu = mu[vcat(i,j)]
            nij = Matrix{Float64}(contingency_table(X[:,i], X[:,j], c))
            sigmaij = sqrt.(diag(sigs))
            for ci in 1:(c)
                for cj in 1:(c)
                    nlp += nij[ci,cj]*log(bvnP((ts[i,ci]-nu[1])/sigmaij[1], (ts[i,ci+1]-nu[1])/sigmaij[1], (ts[j,cj]-nu[2])/sigmaij[2], (ts[j,cj+1]-nu[2])/sigmaij[2], rho,false))
                                                               
                end
            end
        end
    end 
    f = -sum(nlp)
    # @info "Current objective: $f"
    return f
end


function contingency_table(v1::Vector{Int}, v2::Vector{Int}, c::Int)
    # Initialize the contingency table with zeros
    contingency = zeros(Int, c, c)
    
    # Fill the contingency table
    for i in 1:length(v1)
        contingency[v1[i], v2[i]] += 1
    end
    
    return contingency
end

# function care_truegeneration(design,B)
#     local TruePar::Matrix{Float64}, Lambda0::Matrix{Float64}, Theta_delta0::Diagonal{Float64, Vector{Float64}}, Phi0::Matrix{Float64}
    
#     seed = vcat(4352746652536,8401675438322,4933401147650,3759925580001,4726536670900,5703165301676,7833702708054,7792120499866, 648698602044,6859789824488, 814784558462,1263536034076)
#     numpar = 200
#     TruePar = Matrix{Float64}(undef,numpar,size(design,1))
#     for d in 1:size(design,1)
#         n = 1000
#         p = 9
#         q = 3
#         K = design[d,1]
#         propo = design[d,2]
#         careless = design[d,3]
#         c = 4
        

        
#         # True CFA 
#         t10 = vcat(-Inf,0, 1, 3, Inf)
#         A = reshape(Diagonal(ones(q)),q,q);
#         L_str = kron(A,ones(Int64(p/q)));
#         Lambda0 = [0.69, 0.43, 0.89, 0.36, 0.78, 0.96, 0.56, 0.42, 0.67] .* L_str;
#         Phi0 = [1 0.6 0.4;0.6 1 0.2;.4 0.2 1]; #Diagonal(ones(q,q))
#         Theta_delta0 =  Diagonal(1 .-Diagonal((Lambda0*Phi0*Lambda0')));

#         npartot= p*(c-1) + p + size(vec(Phi0),1) #per avere vettori output con dimensioni uniformi nonostante la differenza tra disegni sperimentali
#         TruePar[:,d] = vcat(vec(Lambda0), vec(Phi0), vec(zeros(numpar-npartot))) 

#         Random.seed!(seed[d])
#         # Generation
#         for b in 1:B
#             Eta0 = rand(MvNormal(zeros(q),Phi0),n);
#             Delta0 = rand(MvNormal(zeros(p),Theta_delta0),n);
#             Ystarcfa = Eta0'*Lambda0' + Delta0';
#             Ycfa = Matrix{Int64}(undef,n,p);
#             for j in 1:size(Ystarcfa,2) Ycfa[:,j] = cut(Ystarcfa[:,j],t10,labels=1:c) end
            
#             carel = Matrix{Int64}(undef,n,p)
#             if careless == 0. #random
#                 carel = reshape(wsample([1,2,3,4],[0.25,0.25,0.25,0.25],n*p,replace=true),n,p)
#             elseif careless == 1. #complete straightlining 
#                 carel = reshape(repeat(wsample([1,2,3,4],[0.25,0.25,0.25,0.25],n,replace=true),outer=p),n,p)
#             end 

#             z = rand(Binomial(1, propo),n)
#             X=zeros(Int64,n,p)
#             for i=1:n
#                 if z[i]==1
#                     X[i,:] = vec(Ycfa[i,:]);
#                 elseif z[i]==0
#                     X[i,:] = vec(carel[i,:]);
#                 end
#             end
#         end
#         jldsave(string("/home/nc/MEGA/LatentVariables/data/data_",d,"_",b,".jld2"),data=X)  
#     end
#     jldsave("/home/nc/MEGA/LatentVariables/data/TruePar.jld2",TruePars=TruePar) 
#     return TruePar
# end

function find_ones_columns(matrix)
    # Initialize the result vector for the column indices
    ones_positions = Vector{Int}(undef, size(matrix, 1))
    
    # For each row, find the column where the 1 is located
    @inbounds for i in 1:size(matrix, 1)
        ones_positions[i] = argmax(matrix[i, :])
    end
    
    return ones_positions
end

function classifications(z,Ez)
    z = Int.(z); Ez = Int.(Ez);
    TrueNegatives = sum(Ez.== 1 .&& z.==  1)
    FalseNegatives = sum(Ez.!= 1 .&& z .== 1)
    FalsePositives = sum(Ez.== 1 .&& z .!= 1)
    TruePositives = sum(Ez.!= 1 .&& z .!= 1)


    SP = TrueNegatives/(TrueNegatives+FalsePositives)
    SE = TruePositives/(TruePositives+FalseNegatives)
    #MCC è riscalato a [0,1]
    MCC = (1+(TruePositives*TrueNegatives-FalsePositives*FalseNegatives)/sqrt((TruePositives+FalsePositives)*(TruePositives+FalseNegatives)*(TrueNegatives+FalsePositives)*(TrueNegatives+FalseNegatives)))/2
    BACC = (TruePositives/(TruePositives+FalseNegatives)+(TrueNegatives/(TrueNegatives+FalsePositives))) / 2
    # Cohen = Metrics.cohen_kappa(Int.(Ez), z)
    ARI = Clustering.randindex(z,Int.(Ez))[1]

    return vcat(hcat("BACC","MCC","SE","SP","ARI"),hcat(BACC,MCC,SE,SP,ARI))
end





########################################## Fixed thresholds ##################################################
function startvalsth(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
    liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)+p+(q*(q-1))/2+(c-1-2)+(p*K-(K*(K-1)/2))+(p)+p)+ncov,replix);#p-1
    for i in 1:replix
        println(i)

        tau1 = rand(Uniform(-2.3,-1.)).+repeat(collect(range(0,3.5,c-1)))
        lambda1 = rand(Uniform(0.5,.90),p);
        phi=rand(LKJ(q,1))
        phi=phi[tril!(trues(size(phi)), -1)]
        tau2=vec(reshape(rand(Uniform(1.2,4),1*(c-1-2)),1,(c-1-2)).+reshape(repeat(collect(0:(c-2-2)),inner=1),1,(c-1-2)))
        tau2=sort(tau2)
        lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(0.15,.90),Int((p*K-(K*(K-1)/2))))))
        lambda2=lamb2[lamb2.!=0.]
        psi=rand(Uniform(0.5,1.5),p)
        interc=rand(Uniform(1.2,4),p)



        prop = rand(Uniform(0.65,0.95))
        bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.03,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


        param_initials = vcat(tau1,lambda1,phi,tau2,lambda2,psi,interc,bets)

        # first iteration
        plimix = E_stepMixth(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        params[:,i] = M_parsMixth(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        dt = time() - t
        prrot=params[:,i]
        liksmix[i] = E_stepMixth(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
    end
    param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    lik=liksmix[findall(liksmix.==maximum(liksmix))]

    tau1=param_initials0[(1):((c-1))]
    tau1=vec(sort(reshape(tau1,1,c-1),dims=2))
    lambda1=param_initials0[((c-1)+1):((c-1)+p)]
    phi=param_initials0[((c-1)+p+1):((c-1)+p+Int(q*(q-1)/2))]
    lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    tau2=param_initials0[((c-1)+p+Int(q*(q-1)/2)+1):((c-1)+p+Int(q*(q-1)/2)+(c-1-2))]
    tau2=vec(sort(reshape(tau2,1,c-1-2),dims=2))
    tau2[tau2.<=1.01].=rand(Uniform(1.05,3),sum(tau2.<=1.01))
    lambda2=param_initials0[((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+1):((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2))))]
    psi=param_initials0[((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p))]# psi=param_initials0[((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p-1))]
    interc=param_initials0[((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p)+1):((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p)+p)]#interc=param_initials0[((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
    bets = param_initials0[((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p)+p+1):((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p)+p+ncov)]#    bets = param_initials0[((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1)+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


    param_initialii = vcat(tau1,lambda1,phi,tau2,lambda2,psi,interc,bets)

    println(param_initialii,liksmix,"first iteration: done")
    return param_initialii
end


function pmlefath(param::AbstractVector,p::Int64,K::Int64,c::Int64,X::Matrix{Int64},indices::Vector{Int64})
    t = param[1:(c-1-2)]
    ts = hcat(fill(-Inf,p),zeros(Float64,p),ones(Float64,p),t,fill(Inf,p))
    L2 = reconstruct_L2_matrix(p, K, indices,param[Int((c-1-2)+1):Int((c-1-2)+(p*K-(K*(K-1)/2)))])
    Psi = diagm(param[Int((c-1-2)+(p*K-(K*(K-1)/2))+1):Int((c-1-2)+(p*K-(K*(K-1)/2))+p)])
    interc = param[Int((c-1-2)+(p*K-(K*(K-1)/2))+p+1):Int((c-1-2)+(p*K-(K*(K-1)/2))+p+p)]
    
    Sig = L2*L2' + Psi
    nlp = 0.
    for i in 1:(p-1)
        for j in (i+1):p
            s2c=Sig[i,j]
            sigs = hcat(vcat(Sig[i,i],Sig[i,j]),vcat(Sig[j,i],Sig[j,j]))
            rhos=cov2cor(sigs);
            rho=rhos[1,2];
            nu = interc[vcat(i,j)]
            nij = Matrix{Float64}(contingency_table(X[:,i], X[:,j], c))
            sigmaij = sqrt.(diag(sigs))
            for ci in 1:(c)
                for cj in 1:(c)
                    nlp += nij[ci,cj]*log(bvnP((ts[ci]-nu[1])/sigmaij[1], (ts[ci+1]-nu[1])/sigmaij[1], (ts[cj]-nu[2])/sigmaij[2], (ts[cj+1]-nu[2])/sigmaij[2], rho,false))
                                                               
                end
            end
        end
    end 
    f = -sum(nlp)
    # @info "Current objective: $f"
    return f
end


function pmlefathfixed(param::AbstractVector,p::Int64,K::Int64,c::Int64,X::Matrix{Int64},indices::Vector{Int64},thresholds)
    thresholds = vcat(-Inf,thresholds,Inf)
    L2 = reconstruct_L2_matrix(p, K, indices,param[Int(1):Int((p*K-(K*(K-1)/2)))])
    Psi = diagm(param[Int((p*K-(K*(K-1)/2))+1):Int((p*K-(K*(K-1)/2))+p)])
    interc = param[Int((p*K-(K*(K-1)/2))+p+1):Int((p*K-(K*(K-1)/2))+p+p)]
    
    Sig = L2*L2' + Psi
    nlp = 0.
    for i in 1:(p-1)
        for j in (i+1):p
            s2c=Sig[i,j]
            sigs = hcat(vcat(Sig[i,i],Sig[i,j]),vcat(Sig[j,i],Sig[j,j]))
            rhos=cov2cor(sigs);
            rho=rhos[1,2];
            nu = interc[vcat(i,j)]
            nij = Matrix{Float64}(contingency_table(X[:,i], X[:,j], c))
            sigmaij = sqrt.(diag(sigs))
            for ci in 1:(c)
                for cj in 1:(c)
                    nlp += nij[ci,cj]*log(bvnP((thresholds[ci]-nu[1])/sigmaij[1], (thresholds[ci+1]-nu[1])/sigmaij[1], (thresholds[cj]-nu[2])/sigmaij[2], (thresholds[cj+1]-nu[2])/sigmaij[2], rho,false))
                                                               
                end
            end
        end
    end 
    f = -sum(nlp)
    # @info "Current objective: $f"
    return f
end



function pmlcfath(param::AbstractVector,p::Int64,q::Int64,c::Int64,X::Matrix{Int64},L_str)
    t = param[1:((c-1))]
    t = vcat(-Inf,t,Inf)
    Lambda = param[((c-1)+1):((c-1)+p)].*L_str #reshape(param[((c-1)+1):((c-1)+p*q)],p,q)
    Phi = reconstruct_correlation_matrix(param[Int((c-1)+p+1):Int((c-1)+p+(q*(q-1))/2)])
   
    nlp = 0.
    for i in 1:(p-1)
        for j in (i+1):p
            rho = Lambda[i,:]'*Phi*Lambda[j,:]  
            nij = Matrix{Float64}(contingency_table(X[:,i], X[:,j], c))
            for ci in 1:(c)
                for cj in 1:(c)
                    # println(t[ci]," ",t[ci+1]," ",t[cj]," ",t[cj+1]," ",nij[ci,cj])
                    nlp += nij[ci,cj]*log(bvnP(t[ci], t[ci+1], t[cj], t[cj+1], rho,false))
                
                end
            end
        end
    end 
    f = -sum(nlp)
    # @info "Current objective: $f"
    return f
end



function EMparsMixth(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
    local plimix::Array, T,oldlik
    n = size(X,1)
    er=0;
    # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
    it=0;
    dev=Inf; dev1=Inf; dev2=Inf;

    
    oldlik=-1000000.;
    oldlik2=1000000.;
    oldlikA=-1000000.;
    lik=0.
    

    noerror = true; devnostop = true;  sum_devs = 0;
    while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
        it = it+1;
        param=param01;

        T00 = 1
        plimix = E_stepMixth(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        param01 = M_parsMixth(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
        dt = time() - t
        if param01 == false
            return return(vcat(param01,vec(plimix),vec(T),0.,0.))
        end
        lik = E_stepMixth(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);

        #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
        con = (lik-oldlik)/(oldlik-oldlik2);
        likA = oldlik + (lik-oldlik)/(1-con)
        old_dev = dev1
        dev1 = abs(likA-oldlikA)
        dev2 = abs(lik-oldlik)
        println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
        println(param01)
        dev = abs(lik-oldlik)/abs(oldlik)

        noerror = oldlik<lik

        oldlikA=likA
        oldlik2=oldlik
        oldlik=lik;
        pml_em_log("EMparsMixth", it, lik, dev, dev1, abs(lik - oldlik2), tol, noerror)

        # sum_devs += dev > old_dev
        # devstop = sum_devs <= 3
    end
    
    # classificazioni dopo ultima massimizzazione
    plimix = E_stepMixth(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

    if noerror == true && it<=maxit 
        return vcat(param01,vec(plimix),vec(T),lik,1.)
    else
        return vcat(param01,vec(plimix),vec(T),0.,0.)
    end
end

function E_stepMixth(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
    if sum(param .== 0.) == size(param,1)
        return -Inf
    end
    trs = (c-1-2)
    trsi = (c-1)
    t1s = param[1:trsi]
    L1 = param[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])

    t2s = param[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+trs)]
    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi+p+(q*(q-1))/2+trs+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(param[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = param[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p)]
    # prop = param[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p+1)]
    betas = param[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+p+ncov)]  #  betas = param[Int(trsi+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    # br1 = vcat(0.,1.,t1s)
    t1=vcat(-Inf,t1s,Inf)
    br2 = vcat(0.,1.,t2s)
    t2=vcat(-Inf,br2,Inf)


    Sigms=ones(p,p,2);
    Sigms[:,:,1]=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sigms[:,:,2]=L2*L2'+Psi 

    plimix1=zeros(Float64,n);plimix2=zeros(Float64,n);
    tcfa1=t1[X.+1]
    tcfa=t1[X]
    tefa1=t2[X.+1]
    tefa=t2[X]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)
    
    tuni = cat(tcfa,tefa,dims=3)
    tdui = cat(tcfa1,tefa1,dims=3)
    
    if estep == true
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));
    
        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    # s2c=s2c+0.0001*Matrix(I,2,2)
                    r2c=cov2cor(s2c);
                    up=(tuni[:,pairsindex[:,idx],g])
                    low=(tdui[:,pairsindex[:,idx],g])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    # s2c=s2c+0.0001*Matrix(I,2,2)
                    r2c=cov2cor(s2c);
                    up=(tuni[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                    low=(tdui[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3.*pg[:,g];
            end
        end
    # % probabilities in each row must sum to 1 (application of Bayes law)
        return plimix
    elseif estep == false
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));

        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    # s2c=s2c+0.0001*Matrix(I,2,2)
                    r2c=cov2cor(s2c);
                    up=(tuni[:,pairsindex[:,idx],g])
                    low=(tdui[:,pairsindex[:,idx],g])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    # s2c=s2c+0.0001*Matrix(I,2,2)
                    r2c=cov2cor(s2c);
                    up=(tuni[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                    low=(tdui[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3
            end
        end
        lik1=(log.(plimix));
        lik1=sum(lik1,dims=3);
        U=exp.(lik1)[:,:,1].*pg; 
        cllik=sum(log.(sum(U,dims=2)))
        return cllik
    end
end
    


function pmlth(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trs = (c-1-2)
    trsi = (c-1)
    t1s = x[1:trsi]
    L1 = x[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
   

    t2s = x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+trs)]
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+trs+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = x[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p)]
    # prop = x[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p+1)]

    br1 = t1s
    t1=vcat(-Inf,br1,Inf)
    br2 = vcat(0.,1.,t2s)
    t2=vcat(-Inf,br2,Inf)


    Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi


    plimix1=zeros(Float64,n);plimix2=zeros(Float64,n);
    tcfa1=t1[X.+1]
    tcfa=t1[X]
    tefa1=t2[X.+1]
    tefa=t2[X]

    low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
    for idx=1:size(pairsindex,2) #parallelizzabile
        for g=1:2
            if g==1
                s2c=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
                # s2c=s2c+0.0001*Matrix(I,2,2)
                r2c=cov2cor(s2c);
                B=hcat(tcfa[:,pairsindex[:,idx]],tcfa1[:,pairsindex[:,idx]]);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                Tt1 = T[:,1,idx]
                plimix1 += log.(t3).*Tt1;
            else 
                s2c=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
                # s2c=s2c+0.0001*Matrix(I,2,2)
                r2c=cov2cor(s2c);
                B=hcat((tefa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(s2c)'),(tefa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(s2c)'));
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                Tt2 = T[:,2,idx]
                plimix2 += log.(t3).*Tt2;
            end
            
        end
    end
    f=-sum(plimix1+plimix2)
end

function fgth!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trs = (c-1-2)
    trsi = (c-1)
    t1s = x[1:trsi]
    L1 = x[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
   

    t2s = x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+trs)]
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+trs+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = x[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p)]
    # prop = x[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p+1)]

    br1 = t1s
    t1=vcat(-Inf,br1,Inf)
    br2 = vcat(0.,1.,t2s)
    t2=vcat(-Inf,br2,Inf)


    Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi


    if G !== nothing
        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradintercefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)
    
        gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)
    
        gradlambci = zeros(Float64,2,q,size(pairsindex,2))
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
        
        idx = 0
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1
                rho1 = Sig1[i,j]
                cormat1 = ones(Float64,2,2); cormat1[1,2]=rho1; cormat1[2,1]=rho1
                sigs = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sigs)
                rho2 = cormat2[1,2]
                muij = interc2[vcat(i,j)]
                sigsij = diag(sigs)
                sigij = sqrt.(sigsij)
                sigij13 = sigij[1]^(3)
                sigij23 = sigij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]
                @inbounds for N in 1:n
                    xni = X[N,i]
                    xnj = X[N,j]
                    t1i = t1[xni]
                    t1i1 = t1[xni+1]
                    t1j = t1[xnj]
                    t1j1 = t1[xnj+1]
                    T1 = T[N,1,idx][1]
                    T2 = T[N,2,idx][1]
                    t2i = t2[xni]
                    t2i1 = t2[xni+1]
                    t2j = t2[xnj]
                    t2j1 = t2[xnj+1]
                    t2ist = (t2i-muij[1])/sigij[1]
                    t2i1st = (t2i1-muij[1])/sigij[1]
                    t2jst = (t2j-muij[2])/sigij[2]
                    t2j1st = (t2j1-muij[2])/sigij[2]
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                   
                    pijicfa = bvnP(t1i, t1i1, t1j, t1j1, rho1, false)
                    gradtaucfa[xni+1,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i1)*(cdf(Normal(0,1),(t1j1-rho1*t1i1)/sqrtrij1)-cdf(Normal(0,1),(t1j-rho1*t1i1)/sqrtrij1))) 
                    gradtaucfa[xni,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i)*(cdf(Normal(0,1),(t1j-rho1*t1i)/sqrtrij1)-cdf(Normal(0,1),(t1j1-rho1*t1i)/sqrtrij1))) 
                    gradtaucfa[xnj+1,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j1)*(cdf(Normal(0,1),(t1i1-rho1*t1j1)/sqrtrij1)-cdf(Normal(0,1),(t1i-rho1*t1j1)/sqrtrij1))) 
                    gradtaucfa[xnj,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j)*(cdf(Normal(0,1),(t1i-rho1*t1j)/sqrtrij1)-cdf(Normal(0,1),(t1i1-rho1*t1j)/sqrtrij1))) 
                    
    
                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))
    
    
                    
                    if !(t1i == -Inf && t1j == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1 == Inf && t1j1 == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j1))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end
        
                    if !(t1i1 == Inf && t1j == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end
        
                    if !(t1i == -Inf && t1j1 == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j1))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    gradlambci[1:2,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)
                    gradphici[1:q,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*(L1[i,:]*L1[j,:]')
                    
                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sigij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sigij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sigij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sigij[2])
    
                    gradintercefa[xni+1,idx,1] += (derivNormi1/(-sigij[1]))
                    gradintercefa[xnj+1,idx,2] += (derivNormj1/(-sigij[2]))
                    gradintercefa[xni,idx,1] += (derivNormi/(-sigij[1]))
                    gradintercefa[xnj,idx,2] += (derivNormj/(-sigij[2]))
    
                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)
    
                    
                    gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sigij13)*(sigij[2])))
                    gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sigij23)*(sigij[1])))
    
                    gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij[1]-t2i1))/(2*(sigij13))
                    gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij[2]-t2j1))/(2*(sigij23))
                    gradpsiefa[xni,idx,1] += (derivNormi*(muij[1]-t2i))/(2*(sigij13))
                    gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij[2]-t2j))/(2*(sigij23))
                    
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sigij13*sigij[2])+L2jk/((sigij[1]*sigij[2])))#((L2j*Psi[i,i])/(sigij13*sigij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sigij23*sigij[1])+L2ik/((sigij[2]*sigij[1])))#((L2[i,:].*Psi[j,j])./(sigij23*sigij[1]))
    
                        gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(((muij[1].-t2i1).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(((muij[2].-t2j1).*L2jk)./((sigij23)))
                        gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(((muij[1].-t2i).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(((muij[2].-t2j).*L2jk)./((sigij23)))
                    end 
                
                end
            end
        end
    
        vecpair = vec(pairsindex)    
        liktau = -sum(gradtaucfa,dims=(2,3))[2:(c)]
        liktau2 = -sum(gradtauefa,dims=(2,3))[4:(c)]
        
    
        llamb1 = permutedims(gradlambci,(2,3,1)); llam1 = zeros(p); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
        likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
        likinterc2 = zeros(p); sumgradinter = sum(gradintercefa[2:c,:,:],dims=1)[1,:,:]
        likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        @inbounds for i in 1:p
            llam1[i] = -sum(vec(llamb1[ids[i],:,:]')[vecpair.==i])
            likinterc2[i] = -sum(vec(sumgradinter')[vecpair.==i])
            likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
            likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end
    
        likpsitrue = likpsi.+likpsiidx
        liklstrue = likls.+liklsidx
    
        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]
    
        G .= vec(vcat(liktau,llam1,likphida,liktau2,liklstrue[indices],likpsitrue,likinterc2))
      
    end
    if F !== nothing 
  
    
        plimix1=zeros(Float64,n);plimix2=zeros(Float64,n);
        tcfa1=t1[X.+1]
        tcfa=t1[X]
        tefa1=t2[X.+1]
        tefa=t2[X]
       
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
                    r2c=cov2cor(s2c);
                    B=hcat(tcfa[:,pairsindex[:,idx]],tcfa1[:,pairsindex[:,idx]]);
                    t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                    Tt1 = T[:,1,idx]
                    plimix1 += log.(t3).*Tt1;
                else 
                    s2c=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
                    r2c=cov2cor(s2c);
                    B=hcat((tefa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(s2c)'),(tefa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(s2c)'));
                    t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                    Tt2 = T[:,2,idx]
                    plimix2 += log.(t3).*Tt2;
                end
                
            end
        end
        f=-sum(plimix1+plimix2)
        # @info "Current objective: $f"
        return f
    end
end



function M_parsMixth(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)
    trsi = (c-1) #attenzione p è già qui
    trs = (c-1-2)

    # PML - CFA & EFA
    theta_init = theta[1:(size(theta,1)-ncov)]

    lb=vec(vcat(fill(-5.1,trsi),fill(-1,p),fill(-1.,Int((q*(q-1))/2)),fill(1.,trs),fill(-Inf,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-5,p)));#1.0e-5
    ub=vec(vcat(fill(5.1,trsi),fill(1.,p),fill(1.,Int((q*(q-1))/2)),fill(7.,trs),fill(Inf,Int(p*K-K*(K-1)/2)),fill(Inf,p),fill(10.,p)));    
    par = zeros(Float64,size(lb))
    if it == 1
        par = optimize(Optim.only_fg!((F,G,x)->fgth!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06,time_limit = 300.0)) #20 min
    elseif it > 1
        par = optimize(Optim.only_fg!((F,G,x)->fgth!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06,time_limit = 1200.0)) #20 min
    end 

    
    # mixture parameter
    betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+p+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
    Ez1 = T[:,1,:]
    Ez0 = T[:,2,:]
    
    pairs = size(pairsindex,2)
    betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)
    
    
    # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
   
    # println(par.minimizer,pg[1])
    if !(Optim.converged(par)) && isnan(betas[1]) 
        println(betas[1])
        if start == true
            return zeros(Float64,size(theta,1))
        elseif start == false
            return Optim.converged(par)
        end
    end
    return vcat(par.minimizer,betasnew)#,pg[1])
end














########################################## Fixed and invariant thresholds ##################################################
function startvalsonlyth(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
    liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p)+ncov,replix);#p-1
    for i in 1:replix
        println(i)

        tau1 = rand(Uniform(-2.3,-1.)).+repeat(collect(range(0,3.5,c-1)))
        lambda1 = rand(Uniform(0.5,.90),p);
        phi=rand(LKJ(q,1))
        phi=phi[tril!(trues(size(phi)), -1)]
        lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(0.15,.90),Int((p*K-(K*(K-1)/2))))))
        lambda2=lamb2[lamb2.!=0.]
        psi=rand(Uniform(0.5,1.25),p)
        interc=rand(Uniform(0.2,2),p)



        prop = rand(Uniform(0.65,0.95))
        bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.03,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


        param_initials = vcat(tau1,lambda1,phi,lambda2,psi,interc,bets)

        # first iteration
        plimix = E_stepMixonlyth(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        params[:,i] = M_parsMixonlyth(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        dt = time() - t
        prrot=params[:,i]
        liksmix[i] = E_stepMixonlyth(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
    end
    param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    lik=liksmix[findall(liksmix.==maximum(liksmix))]

    tau1=param_initials0[(1):((c-1))]
    tau1=vec(sort(reshape(tau1,1,c-1),dims=2))
    lambda1=param_initials0[((c-1)+1):((c-1)+p)]
    phi=param_initials0[((c-1)+p+1):((c-1)+p+Int(q*(q-1)/2))]
    lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    lambda2=param_initials0[((c-1)+p+Int(q*(q-1)/2)+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2))))]
    psi=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p))]# psi=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1))]
    interc=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+p)]#interc=param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
    bets = param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+p+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+p+ncov)]#    bets = param_initials0[((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


    param_initialii = vcat(tau1,lambda1,phi,lambda2,psi,interc,bets)

    println(param_initialii,liksmix,"first iteration: done")
    return param_initialii
end



function EMparsMixonlyth(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
    local plimix::Array, T,oldlik
    n = size(X,1)
    er=0;
    # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
    it=0;
    dev=Inf; dev1=Inf; dev2=Inf;

    
    oldlik=-1000000.;
    oldlik2=1000000.;
    oldlikA=-1000000.;
    lik=0.
    

    noerror = true; devnostop = true;  sum_devs = 0;
    while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
        it = it+1;
        param=param01;

        T00 = 1
        plimix = E_stepMixonlyth(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        param01 = M_parsMixonlyth(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
        dt = time() - t
        if param01 == false
            return return(vcat(param01,vec(plimix),vec(T),0.,0.))
        end
        lik = E_stepMixonlyth(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);

        #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
        con = (lik-oldlik)/(oldlik-oldlik2);
        likA = oldlik + (lik-oldlik)/(1-con)
        old_dev = dev1
        dev1 = abs(likA-oldlikA)
        dev2 = abs(lik-oldlik)
        println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
        println(param01)
        dev = abs(lik-oldlik)/abs(oldlik)

        noerror = oldlik<lik

        oldlikA=likA
        oldlik2=oldlik
        oldlik=lik;
        pml_em_log("EMparsMixonlyth", it, lik, dev, dev1, abs(lik - oldlik2), tol, noerror)

        # sum_devs += dev > old_dev
        # devstop = sum_devs <= 3
    end
    
    # classificazioni dopo ultima massimizzazione
    plimix = E_stepMixonlyth(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

    if noerror == true && it<=maxit 
        return vcat(param01,vec(plimix),vec(T),lik,1.)
    else
        return vcat(param01,vec(plimix),vec(T),0.,0.)
    end
end

function E_stepMixonlyth(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
    if sum(param .== 0.) == size(param,1)
        return -Inf
    end 
    trsi = (c-1)
    t1s = param[1:trsi]
    L1 = param[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])

    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+p)]
    # prop = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+p+1)]
    betas = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p+ncov)]  #  betas = param[Int(trsi+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    t1=vcat(-Inf,t1s,Inf)

    Sigms=ones(p,p,2);
    Sigms[:,:,1]=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sigms[:,:,2]=L2*L2'+Psi 

    plimix1=zeros(Float64,n);plimix2=zeros(Float64,n);
    tfa1=t1[X.+1]
    tfa=t1[X]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)
    
    if estep == true
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));
    
        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]])
                    low=(tfa[:,pairsindex[:,idx]])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
                    low=(tfa[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3.*pg[:,g];
            end
        end
    # % probabilities in each row must sum to 1 (application of Bayes law)
        return plimix
    elseif estep == false
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));

        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]])
                    low=(tfa[:,pairsindex[:,idx]])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                    low=(tfa[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3
            end
        end
        lik1=(log.(plimix));
        lik1=sum(lik1,dims=3);
        U=exp.(lik1)[:,:,1].*pg; 
        cllik=sum(log.(sum(U,dims=2)))
        return cllik
    end
end
    

function pmllastonlyth(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    trsi = (c-1)
    t1s = x[1:trsi]
    L1 = x[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
   
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+p)]
   
    Sig1=L1*Phi*L1' + Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi

   
    t1=vcat(-Inf,t1s,Inf)
    tfa1=t1[X.+1]
    tfa=t1[X]
    
    plimix=zeros(Float64,n);
    for idx=1:size(pairsindex,2) #parallelizzabile
        # CFA
        Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
        Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
        tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
        Tt1 = T[:,1,idx]
        
        #EFA
        Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
        Refa=cov2cor(Sefa);
        Befa=hcat((tfa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
        tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
        Tt2 = T[:,2,idx]

        plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
    end
    f=-sum(plimix)
return f
end




function fgonlyth!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trsi = (c-1)
    t1s = x[1:trsi]
    L1 = x[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
   

    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+p)]
    # prop = x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+p+1)]

    t1=vcat(-Inf,t1s,Inf)
 

    Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi


    if G !== nothing
        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradintercefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)
    
        gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)
    
        gradlambci = zeros(Float64,2,q,size(pairsindex,2))
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
        
        idx = 0
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1
                rho1 = Sig1[i,j]
                cormat1 = ones(Float64,2,2); cormat1[1,2]=rho1; cormat1[2,1]=rho1
                sigs = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sigs)
                rho2 = cormat2[1,2]
                muij = interc2[vcat(i,j)]
                sigsij = diag(sigs)
                sigij = sqrt.(sigsij)
                sigij13 = sigij[1]^(3)
                sigij23 = sigij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]
                @inbounds for N in 1:n
                    xni = X[N,i]
                    xnj = X[N,j]
                    t1i = t1[xni]
                    t1i1 = t1[xni+1]
                    t1j = t1[xnj]
                    t1j1 = t1[xnj+1]
                    T1 = T[N,1,idx][1]
                    T2 = T[N,2,idx][1]
                    t2i = t1[xni]
                    t2i1 = t1[xni+1]
                    t2j = t1[xnj]
                    t2j1 = t1[xnj+1]
                    t2ist = (t2i-muij[1])/sigij[1]
                    t2i1st = (t2i1-muij[1])/sigij[1]
                    t2jst = (t2j-muij[2])/sigij[2]
                    t2j1st = (t2j1-muij[2])/sigij[2]
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                   
                    pijicfa = bvnP(t1i, t1i1, t1j, t1j1, rho1, false)
                    gradtaucfa[xni+1,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i1)*(cdf(Normal(0,1),(t1j1-rho1*t1i1)/sqrtrij1)-cdf(Normal(0,1),(t1j-rho1*t1i1)/sqrtrij1))) 
                    gradtaucfa[xni,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i)*(cdf(Normal(0,1),(t1j-rho1*t1i)/sqrtrij1)-cdf(Normal(0,1),(t1j1-rho1*t1i)/sqrtrij1))) 
                    gradtaucfa[xnj+1,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j1)*(cdf(Normal(0,1),(t1i1-rho1*t1j1)/sqrtrij1)-cdf(Normal(0,1),(t1i-rho1*t1j1)/sqrtrij1))) 
                    gradtaucfa[xnj,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j)*(cdf(Normal(0,1),(t1i-rho1*t1j)/sqrtrij1)-cdf(Normal(0,1),(t1i1-rho1*t1j)/sqrtrij1))) 
                    
    
                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))
    
    
                    
                    if !(t1i == -Inf && t1j == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1 == Inf && t1j1 == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j1))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end
        
                    if !(t1i1 == Inf && t1j == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end
        
                    if !(t1i == -Inf && t1j1 == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j1))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    gradlambci[1:2,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)
                    gradphici[1:q,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*(L1[i,:]*L1[j,:]')
                    
                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sigij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sigij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sigij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sigij[2])
    
                    gradintercefa[xni+1,idx,1] += (derivNormi1/(-sigij[1]))
                    gradintercefa[xnj+1,idx,2] += (derivNormj1/(-sigij[2]))
                    gradintercefa[xni,idx,1] += (derivNormi/(-sigij[1]))
                    gradintercefa[xnj,idx,2] += (derivNormj/(-sigij[2]))
    
                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)
    
                    
                    gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sigij13)*(sigij[2])))
                    gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sigij23)*(sigij[1])))
    
                    gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij[1]-t2i1))/(2*(sigij13))
                    gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij[2]-t2j1))/(2*(sigij23))
                    gradpsiefa[xni,idx,1] += (derivNormi*(muij[1]-t2i))/(2*(sigij13))
                    gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij[2]-t2j))/(2*(sigij23))
                    
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sigij13*sigij[2])+L2jk/((sigij[1]*sigij[2])))#((L2j*Psi[i,i])/(sigij13*sigij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sigij23*sigij[1])+L2ik/((sigij[2]*sigij[1])))#((L2[i,:].*Psi[j,j])./(sigij23*sigij[1]))
    
                        gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(((muij[1].-t2i1).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(((muij[2].-t2j1).*L2jk)./((sigij23)))
                        gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(((muij[1].-t2i).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(((muij[2].-t2j).*L2jk)./((sigij23)))
                    end 
                
                end
            end
        end
    
        vecpair = vec(pairsindex)    
        likshared = -sum(gradtaucfa,dims=(2,3))[2:(c)] -sum(gradtauefa,dims=(2,3))[2:(c)]
        
    
        llamb1 = permutedims(gradlambci,(2,3,1)); llam1 = zeros(p); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
        likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
        likinterc2 = zeros(p); sumgradinter = sum(gradintercefa[2:c,:,:],dims=1)[1,:,:]
        likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        @inbounds for i in 1:p
            llam1[i] = -sum(vec(llamb1[ids[i],:,:]')[vecpair.==i])
            likinterc2[i] = -sum(vec(sumgradinter')[vecpair.==i])
            likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
            likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end
    
        likpsitrue = likpsi.+likpsiidx
        liklstrue = likls.+liklsidx
    
        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]
    
        G .= vec(vcat(likshared,llam1,likphida,liklstrue[indices],likpsitrue,likinterc2))
      
    end
    if F !== nothing 
        plimix=zeros(Float64,n);
        tfa1=t1[X.+1]
        tfa=t1[X]
        
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
            Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
            Tt1 = T[:,1,idx]
            
            # EFA
            Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tfa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            Tt2 = T[:,2,idx]
    
            plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
        end
        f=-sum(plimix)
        return f
    end
end



function M_parsMixonlyth(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)
    trsi = (c-1) #attenzione p è già qui
 
    # PML - CFA & EFA
    theta_init = theta[1:(size(theta,1)-ncov)]

    lb=vec(vcat(fill(-5.1,trsi),fill(-1,p),fill(-1.,Int((q*(q-1))/2)),fill(-Inf,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-5,p)));#1.0e-5
    ub=vec(vcat(fill(5.1,trsi),fill(1.,p),fill(1.,Int((q*(q-1))/2)),fill(Inf,Int(p*K-K*(K-1)/2)),fill(Inf,p),fill(10.,p)));    
    par = zeros(Float64,size(lb))
    if it == 1
        par = optimize(Optim.only_fg!((F,G,x)->fgonlyth!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06,time_limit = 300.0)) #20 min
    elseif it > 1
        par = optimize(Optim.only_fg!((F,G,x)->fgonlyth!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06,time_limit = 1200.0)) #20 min
    end 

    
    # mixture parameter
    betas = theta[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
    Ez1 = T[:,1,:]
    Ez0 = T[:,2,:]
    
    pairs = size(pairsindex,2)
    betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)
    
    
    # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
   
    # println(par.minimizer,pg[1])
    if !(Optim.converged(par)) && isnan(betas[1]) 
        println(betas[1])
        if start == true
            return zeros(Float64,size(theta,1))
        elseif start == false
            return Optim.converged(par)
        end
    end
    return vcat(par.minimizer,betasnew)#,pg[1])
end















########################################## Fixed and invariant thresholds con 0 e 1 ##################################################
function startvalsonlyth01(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
    local tau1
    liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1-2)+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p)+ncov,replix);#p-1
    for i in 1:replix
        println(i)

        # if c==4 tau1 = rand(Uniform(1,1.2)).+repeat(collect(range(1,1,c-1-2))) else rand(Uniform(1,1.2)).+repeat(collect(range(0.5,1.5,c-1-2))) end
        tau1 = rand(Uniform(1,2.2)).+repeat(collect(range(0.5,1.5,c-1-2)))
        lambda1 = rand(Uniform(0.5,.90),p);
        phi=rand(LKJ(q,1))
        phi=phi[tril!(trues(size(phi)), -1)]
        lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(0.15,.90),Int((p*K-(K*(K-1)/2))))))
        lambda2=lamb2[lamb2.!=0.]
        psi=rand(Uniform(0.5,1.25),p)
        interc=rand(Uniform(0.2,2),p)



        prop = rand(Uniform(0.65,0.95))
        bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.03,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


        param_initials = vcat(tau1,lambda1,phi,lambda2,psi,interc,bets)

        # first iteration
        plimix = E_stepMixonlyth01(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        params[:,i] = M_parsMixonlyth01(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        dt = time() - t
        prrot=params[:,i]
        liksmix[i] = E_stepMixonlyth01(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
    end
    param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    lik=liksmix[findall(liksmix.==maximum(liksmix))]

    tau1=param_initials0[(1):((c-1-2))]
    tau1=vec(sort(reshape(tau1,1,c-1-2),dims=2))
    lambda1=param_initials0[((c-1-2)+1):((c-1-2)+p)]
    phi=param_initials0[((c-1-2)+p+1):((c-1-2)+p+Int(q*(q-1)/2))]
    lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    lambda2=param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+1):((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2))))]
    psi=param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p))]# psi=param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+1):((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1))]
    interc=param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+1):((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+p)]#interc=param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
    bets = param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+p+1):((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p)+p+ncov)]#    bets = param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


    param_initialii = vcat(tau1,lambda1,phi,lambda2,psi,interc,bets)

    println(param_initialii,liksmix,"first iteration: done")
    return param_initialii
end



function EMparsMixonlyth01(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
    local plimix::Array, T,oldlik
    n = size(X,1)
    er=0;
    # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
    it=0;
    dev=Inf; dev1=Inf; dev2=Inf;

    
    oldlik=-1000000.;
    oldlik2=1000000.;
    oldlikA=-1000000.;
    lik=0.
    

    noerror = true; devnostop = true;  sum_devs = 0;
    while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
        it = it+1;
        param=param01;

        T00 = 1
        plimix = E_stepMixonlyth01(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        param01 = M_parsMixonlyth01(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
        dt = time() - t
        if param01 == false
            return return(vcat(param01,vec(plimix),vec(T),0.,0.))
        end
        lik = E_stepMixonlyth01(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);

        #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
        con = (lik-oldlik)/(oldlik-oldlik2);
        likA = oldlik + (lik-oldlik)/(1-con)
        old_dev = dev1
        dev1 = abs(likA-oldlikA)
        dev2 = abs(lik-oldlik)
        println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
        println(param01)
        dev = abs(lik-oldlik)/abs(oldlik)

        noerror = oldlik<lik

        oldlikA=likA
        oldlik2=oldlik
        oldlik=lik;
        pml_em_log("EMparsMixonlyth01", it, lik, dev, dev1, abs(lik - oldlik2), tol, noerror)

        # sum_devs += dev > old_dev
        # devstop = sum_devs <= 3
    end
    
    # classificazioni dopo ultima massimizzazione
    plimix = E_stepMixonlyth01(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

    if noerror == true && it<=maxit 
        return vcat(param01,vec(plimix),vec(T),lik,1.)
    else
        return vcat(param01,vec(plimix),vec(T),0.,0.)
    end
end

function E_stepMixonlyth01(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
    if sum(param .== 0.) == size(param,1)
        return -Inf
    end 
    trsi = (c-1-2)
    t1s = param[1:trsi]
    L1 = param[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])

    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+p)]
    # prop = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+p+1)]
    betas = param[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p+ncov)]  #  betas = param[Int(trsi+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    t1=vcat(-Inf,0,1,t1s,Inf)

    Sigms=ones(p,p,2);
    Sigms[:,:,1]=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sigms[:,:,2]=L2*L2'+Psi 

    plimix1=zeros(Float64,n);plimix2=zeros(Float64,n);
    tfa1=t1[X.+1]
    tfa=t1[X]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)
    
    if estep == true
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));
    
        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]])
                    low=(tfa[:,pairsindex[:,idx]])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
                    low=(tfa[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3.*pg[:,g];
            end
        end
    # % probabilities in each row must sum to 1 (application of Bayes law)
        return plimix
    elseif estep == false
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));

        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]])
                    low=(tfa[:,pairsindex[:,idx]])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                    low=(tfa[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3
            end
        end
        lik1=(log.(plimix));
        lik1=sum(lik1,dims=3);
        U=exp.(lik1)[:,:,1].*pg; 
        cllik=sum(log.(sum(U,dims=2)))
        return cllik
    end
end
    

function pmllastonlyth01(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    trsi = (c-1-2)
    t1s = x[1:trsi]
    L1 = x[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
   
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+p)]
   
    Sig1=L1*Phi*L1' + Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi

   
    t1=vcat(-Inf,0,1,t1s,Inf)
    tfa1=t1[X.+1]
    tfa=t1[X]
    
    plimix=zeros(Float64,n);
    for idx=1:size(pairsindex,2) #parallelizzabile
        # CFA
        Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
        Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
        tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
        Tt1 = T[:,1,idx]
        
        #EFA
        Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
        Refa=cov2cor(Sefa);
        Befa=hcat((tfa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
        tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
        Tt2 = T[:,2,idx]

        plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
    end
    f=-sum(plimix)
return f
end




function fgonlyth01!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trsi = (c-1-2)
    t1s = x[1:trsi]
    L1 = x[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
   

    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+p)]
    # prop = x[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+p+p+1)]

    t1=vcat(-Inf,0,1,t1s,Inf)
 

    Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi


    if G !== nothing
        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradintercefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)
    
        gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)
    
        gradlambci = zeros(Float64,2,q,size(pairsindex,2))
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
        
        idx = 0
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1
                rho1 = Sig1[i,j]
                cormat1 = ones(Float64,2,2); cormat1[1,2]=rho1; cormat1[2,1]=rho1
                sigs = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sigs)
                rho2 = cormat2[1,2]
                muij = interc2[vcat(i,j)]
                sigsij = diag(sigs)
                sigij = sqrt.(sigsij)
                sigij13 = sigij[1]^(3)
                sigij23 = sigij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]
                @inbounds for N in 1:n
                    xni = X[N,i]
                    xnj = X[N,j]
                    t1i = t1[xni]
                    t1i1 = t1[xni+1]
                    t1j = t1[xnj]
                    t1j1 = t1[xnj+1]
                    T1 = T[N,1,idx][1]
                    T2 = T[N,2,idx][1]
                    t2i = t1[xni]
                    t2i1 = t1[xni+1]
                    t2j = t1[xnj]
                    t2j1 = t1[xnj+1]
                    t2ist = (t2i-muij[1])/sigij[1]
                    t2i1st = (t2i1-muij[1])/sigij[1]
                    t2jst = (t2j-muij[2])/sigij[2]
                    t2j1st = (t2j1-muij[2])/sigij[2]
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                   
                    pijicfa = bvnP(t1i, t1i1, t1j, t1j1, rho1, false)
                    gradtaucfa[xni+1,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i1)*(cdf(Normal(0,1),(t1j1-rho1*t1i1)/sqrtrij1)-cdf(Normal(0,1),(t1j-rho1*t1i1)/sqrtrij1))) 
                    gradtaucfa[xni,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i)*(cdf(Normal(0,1),(t1j-rho1*t1i)/sqrtrij1)-cdf(Normal(0,1),(t1j1-rho1*t1i)/sqrtrij1))) 
                    gradtaucfa[xnj+1,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j1)*(cdf(Normal(0,1),(t1i1-rho1*t1j1)/sqrtrij1)-cdf(Normal(0,1),(t1i-rho1*t1j1)/sqrtrij1))) 
                    gradtaucfa[xnj,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j)*(cdf(Normal(0,1),(t1i-rho1*t1j)/sqrtrij1)-cdf(Normal(0,1),(t1i1-rho1*t1j)/sqrtrij1))) 
                    
    
                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))
    
    
                    
                    if !(t1i == -Inf && t1j == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1 == Inf && t1j1 == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j1))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end
        
                    if !(t1i1 == Inf && t1j == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end
        
                    if !(t1i == -Inf && t1j1 == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j1))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    gradlambci[1:2,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)
                    gradphici[1:q,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*(L1[i,:]*L1[j,:]')
                    
                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sigij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sigij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sigij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sigij[2])
    
                    gradintercefa[xni+1,idx,1] += (derivNormi1/(-sigij[1]))
                    gradintercefa[xnj+1,idx,2] += (derivNormj1/(-sigij[2]))
                    gradintercefa[xni,idx,1] += (derivNormi/(-sigij[1]))
                    gradintercefa[xnj,idx,2] += (derivNormj/(-sigij[2]))
    
                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)
    
                    
                    gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sigij13)*(sigij[2])))
                    gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sigij23)*(sigij[1])))
    
                    gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij[1]-t2i1))/(2*(sigij13))
                    gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij[2]-t2j1))/(2*(sigij23))
                    gradpsiefa[xni,idx,1] += (derivNormi*(muij[1]-t2i))/(2*(sigij13))
                    gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij[2]-t2j))/(2*(sigij23))
                    
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sigij13*sigij[2])+L2jk/((sigij[1]*sigij[2])))#((L2j*Psi[i,i])/(sigij13*sigij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sigij23*sigij[1])+L2ik/((sigij[2]*sigij[1])))#((L2[i,:].*Psi[j,j])./(sigij23*sigij[1]))
    
                        gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(((muij[1].-t2i1).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(((muij[2].-t2j1).*L2jk)./((sigij23)))
                        gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(((muij[1].-t2i).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(((muij[2].-t2j).*L2jk)./((sigij23)))
                    end 
                
                end
            end
        end
    
        vecpair = vec(pairsindex)    
        likshared = -sum(gradtaucfa,dims=(2,3))[4:(c)] -sum(gradtauefa,dims=(2,3))[4:(c)]
        
    
        llamb1 = permutedims(gradlambci,(2,3,1)); llam1 = zeros(p); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
        likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
        likinterc2 = zeros(p); sumgradinter = sum(gradintercefa[2:c,:,:],dims=1)[1,:,:]
        likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        @inbounds for i in 1:p
            llam1[i] = -sum(vec(llamb1[ids[i],:,:]')[vecpair.==i])
            likinterc2[i] = -sum(vec(sumgradinter')[vecpair.==i])
            likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
            likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end
    
        likpsitrue = likpsi.+likpsiidx
        liklstrue = likls.+liklsidx
    
        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]
    
        G .= vec(vcat(likshared,llam1,likphida,liklstrue[indices],likpsitrue,likinterc2))
      
    end
    if F !== nothing 
        plimix=zeros(Float64,n);
        tfa1=t1[X.+1]
        tfa=t1[X]
        
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
            Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
            Tt1 = T[:,1,idx]
            
            # EFA
            Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tfa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            Tt2 = T[:,2,idx]
    
            plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
        end
        f=-sum(plimix)
        return f
    end
end



function M_parsMixonlyth01(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)
    trsi = (c-1-2) #attenzione p è già qui
 
    # PML - CFA & EFA
    theta_init = theta[1:(size(theta,1)-ncov)]

    lb=vec(vcat(fill(-5.1,trsi),fill(-1,p),fill(-1.,Int((q*(q-1))/2)),fill(-Inf,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-5,p)));#1.0e-5
    ub=vec(vcat(fill(5.1,trsi),fill(1.,p),fill(1.,Int((q*(q-1))/2)),fill(Inf,Int(p*K-K*(K-1)/2)),fill(Inf,p),fill(10.,p)));    
    par = zeros(Float64,size(lb))
    if it == 1
        par = optimize(Optim.only_fg!((F,G,x)->fgonlyth01!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06,time_limit = 300.0)) #20 min
    elseif it > 1
        par = optimize(Optim.only_fg!((F,G,x)->fgonlyth01!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06,time_limit = 1200.0)) #20 min
    end 

    
    # mixture parameter
    betas = theta[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
    Ez1 = T[:,1,:]
    Ez0 = T[:,2,:]
    
    pairs = size(pairsindex,2)
    betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)
    
    
    # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
   
    # println(par.minimizer,pg[1])
    if !(Optim.converged(par)) && isnan(betas[1]) 
        println(betas[1])
        if start == true
            return zeros(Float64,size(theta,1))
        elseif start == false
            return Optim.converged(par)
        end
    end
    return vcat(par.minimizer,betasnew)#,pg[1])
end













########################################## intercetta CFA + Fixed and invariant thresholds con 0 e 1 ##################################################
function startvalsonlyth01intercCFA(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
    local tau1
    liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1-2)+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+(p)+p)+ncov,replix);#p-1
    for i in 1:replix
        println(i)

        # if c==4 tau1 = rand(Uniform(1,1.2)).+repeat(collect(range(1,1,c-1-2))) else tau1 = rand(Uniform(1,1.2)).+repeat(collect(range(0.5,1.5,c-1-2))) end
        tau1 = rand(Uniform(1.5,2.5))#.+repeat(collect(range(0.5,1.5,c-1-2)))
        lambda1 = rand(Uniform(0.5,.90),p);
        phi=rand(LKJ(q,1))
        phi=phi[tril!(trues(size(phi)), -1)]
        interc1=rand(Uniform(0.2,2),p)
        lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(0.15,.90),Int((p*K-(K*(K-1)/2))))))
        lambda2=lamb2[lamb2.!=0.]
        psi=rand(Uniform(0.5,1.25),p)
        interc2=rand(Uniform(-2,2),p)



        prop = rand(Uniform(0.65,0.95))
        bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.03,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


        param_initials = vcat(tau1,lambda1,phi,interc1,lambda2,psi,interc2,bets)

        # first iteration
        plimix = E_stepMixonlyth01intercCFA(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        params[:,i] = M_parsMixonlyth01intercCFA(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        dt = time() - t
        prrot=params[:,i]
        liksmix[i] = E_stepMixonlyth01intercCFA(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
    end
    param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    lik=liksmix[findall(liksmix.==maximum(liksmix))]

    tau1=param_initials0[(1):((c-1-2))]
    tau1=vec(sort(reshape(tau1,1,c-1-2),dims=2))
    lambda1=param_initials0[((c-1-2)+1):((c-1-2)+p)]
    phi=param_initials0[((c-1-2)+p+1):((c-1-2)+p+Int(q*(q-1)/2))]
    lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    interc1=param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+1):((c-1-2)+p+Int(q*(q-1)/2)+p)]#interc=param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1-2)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    lambda2=param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+p+1):((c-1-2)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2))))]
    psi=param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+1):((c-1-2)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p))]# psi=param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+1):((c-1-2)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p-1))]
    interc2=param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p)+1):((c-1-2)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p)+p)]#interc=param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1-2)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
    bets = param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p)+p+1):((c-1-2)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p)+p+ncov)]#    bets = param_initials0[((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1-2)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


    param_initialii = vcat(tau1,lambda1,phi,interc1,lambda2,psi,interc2,bets)

    println(param_initialii,liksmix,"first iteration: done")
    return param_initialii
end



function EMparsMixonlyth01intercCFA(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
    local plimix::Array, T,oldlik
    n = size(X,1)
    er=0;
    # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
    it=0;
    dev=Inf; dev1=Inf; dev2=Inf;

    
    oldlik=-1000000.;
    oldlik2=1000000.;
    oldlikA=-1000000.;
    lik=0.
    

    noerror = true; devnostop = true;  sum_devs = 0;
    while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
        it = it+1;
        param=param01;

        T00 = 1
        plimix = E_stepMixonlyth01intercCFA(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        param01 = M_parsMixonlyth01intercCFA(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
        dt = time() - t
        if param01 == false
            return return(vcat(param01,vec(plimix),vec(T),0.,0.))
        end
        lik = E_stepMixonlyth01intercCFA(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);

        #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
        con = (lik-oldlik)/(oldlik-oldlik2);
        likA = oldlik + (lik-oldlik)/(1-con)
        old_dev = dev1
        dev1 = abs(likA-oldlikA)
        dev2 = abs(lik-oldlik)
        println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
        println(param01)
        dev = abs(lik-oldlik)/abs(oldlik)

        noerror = oldlik<lik

        oldlikA=likA
        oldlik2=oldlik
        oldlik=lik;
        pml_em_log("EMparsMixonlyth01intercCFA", it, lik, dev, dev1, abs(lik - oldlik2), tol, noerror)

        # sum_devs += dev > old_dev
        # devstop = sum_devs <= 3
    end
    
    # classificazioni dopo ultima massimizzazione
    plimix = E_stepMixonlyth01intercCFA(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

    if noerror == true && it<=maxit 
        return vcat(param01,vec(plimix),vec(T),lik,1.)
    else
        return vcat(param01,vec(plimix),vec(T),0.,0.)
    end
end

function E_stepMixonlyth01intercCFA(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
    if sum(param .== 0.) == size(param,1)
        return -Inf
    end 
    trsi = (c-1-2)
    t1s = param[1:trsi]
    L1 = param[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
    interc1 = param[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+p)]

    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi+p+(q*(q-1))/2+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(param[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = param[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+p)]
    # prop = param[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+p+1)]
    betas = param[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+(p)+p+ncov)]  #  betas = param[Int(trsi+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    t1=vcat(-Inf,0,1,t1s,Inf)

    Sigms=ones(p,p,2);
    Sigms[:,:,1]=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sigms[:,:,2]=L2*L2'+Psi 

    plimix1=zeros(Float64,n);plimix2=zeros(Float64,n);
    tfa1=t1[X.+1]
    tfa=t1[X]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)
    
    if estep == true
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));
    
        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(interc1)[pairsindex[:,idx]]')
                    low=(tfa[:,pairsindex[:,idx]].-(interc1)[pairsindex[:,idx]]')
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
                    low=(tfa[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3.*pg[:,g];
            end
        end
    # % probabilities in each row must sum to 1 (application of Bayes law)
        return plimix
    elseif estep == false
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));

        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(interc1)[pairsindex[:,idx]]')
                    low=(tfa[:,pairsindex[:,idx]].-(interc1)[pairsindex[:,idx]]')
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                    low=(tfa[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3
            end
        end
        lik1=(log.(plimix));
        lik1=sum(lik1,dims=3);
        U=exp.(lik1)[:,:,1].*pg; 
        cllik=sum(log.(sum(U,dims=2)))
        return cllik
    end
end
    

function pmllastonlyth01intercCFA(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    trsi = (c-1-2)
    t1s = x[1:trsi]
    L1 = x[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
    interc1 = x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+p)]

    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = x[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+p)]
   
    Sig1=L1*Phi*L1' + Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi

   
    t1=vcat(-Inf,0,1,t1s,Inf)
    tfa1=t1[X.+1]
    tfa=t1[X]
    
    plimix=zeros(Float64,n);
    for idx=1:size(pairsindex,2) #parallelizzabile
        # CFA
        Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
        Bcfa=hcat(tfa[:,pairsindex[:,idx]].-interc1[pairsindex[:,idx]]',tfa1[:,pairsindex[:,idx]].-interc1[pairsindex[:,idx]]');
        tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
        Tt1 = T[:,1,idx]
        
        #EFA
        Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
        Refa=cov2cor(Sefa);
        Befa=hcat((tfa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
        tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
        Tt2 = T[:,2,idx]

        plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
    end
    f=-sum(plimix)
return f
end




function fgonlyth01intercCFA!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trsi = (c-1-2)
    t1s = x[1:trsi]
    L1 = x[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
    interc1 = x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+p)]

    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2)))] )
    Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = x[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+p)]
   
    t1=vcat(-Inf,0,1,t1s,Inf)
 

    Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi


    if G !== nothing
        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradinterccfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradintercefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)
    
        gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)
    
        gradlambci = zeros(Float64,2,q,size(pairsindex,2))
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
        
        idx = 0
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1
                rho1 = Sig1[i,j]
                cormat1 = ones(Float64,2,2); cormat1[1,2]=rho1; cormat1[2,1]=rho1
                sigs = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sigs)
                rho2 = cormat2[1,2]
                muij1 = interc1[vcat(i,j)]
                muij2 = interc2[vcat(i,j)]
                sigsij = diag(sigs)
                sigij = sqrt.(sigsij)
                sigij13 = sigij[1]^(3)
                sigij23 = sigij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]
                @inbounds for N in 1:n
                    xni = X[N,i]
                    xnj = X[N,j]
                    t1i = t1[xni]
                    t1i1 = t1[xni+1]
                    t1j = t1[xnj]
                    t1j1 = t1[xnj+1]
                    t1ist = (t1i-muij1[1])
                    t1i1st = (t1i1-muij1[1])
                    t1jst = (t1j-muij1[2])
                    t1j1st = (t1j1-muij1[2])
                    T1 = T[N,1,idx][1]
                    T2 = T[N,2,idx][1]
                    t2i = t1[xni]
                    t2i1 = t1[xni+1]
                    t2j = t1[xnj]
                    t2j1 = t1[xnj+1]
                    t2ist = (t2i-muij2[1])/sigij[1]
                    t2i1st = (t2i1-muij2[1])/sigij[1]
                    t2jst = (t2j-muij2[2])/sigij[2]
                    t2j1st = (t2j1-muij2[2])/sigij[2]
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                   
                    pijicfa = bvnP(t1ist, t1i1st, t1jst, t1j1st, rho1, false)
                    derivcfai1 = (T1/pijicfa)*(pdf(Normal(0,1),t1i1st)*(cdf(Normal(0,1),(t1j1st-rho1*t1i1st)/sqrtrij1)-cdf(Normal(0,1),(t1jst-rho1*t1i1st)/sqrtrij1))) 
                    derivcfaj1 = (T1/pijicfa)*(pdf(Normal(0,1),t1j1st)*(cdf(Normal(0,1),(t1i1st-rho1*t1j1st)/sqrtrij1)-cdf(Normal(0,1),(t1ist-rho1*t1j1st)/sqrtrij1))) 
                    derivcfai =  (T1/pijicfa)*(pdf(Normal(0,1),t1ist)*(cdf(Normal(0,1),(t1jst-rho1*t1ist)/sqrtrij1)-cdf(Normal(0,1),(t1j1st-rho1*t1ist)/sqrtrij1))) 
                    derivcfaj = (T1/pijicfa)*(pdf(Normal(0,1),t1jst)*(cdf(Normal(0,1),(t1ist-rho1*t1jst)/sqrtrij1)-cdf(Normal(0,1),(t1i1st-rho1*t1jst)/sqrtrij1))) 
                    
                    gradtaucfa[xni+1,idx,1] += derivcfai1
                    gradtaucfa[xnj+1,idx,2] += derivcfaj1
                    gradtaucfa[xni,idx,1] += derivcfai
                    gradtaucfa[xnj,idx,2] += derivcfaj

                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))
    
    
                    
                    if !(t1ist == -Inf && t1jst == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1jst))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1st == Inf && t1j1st == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1j1st))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end
        
                    if !(t1i1st == Inf && t1j == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1jst))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end
        
                    if !(t1i == -Inf && t1j1 == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1j1st))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    gradlambci[1:2,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)
                    gradphici[1:q,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*(L1[i,:]*L1[j,:]')
                    
                    gradinterccfa[xni+1,idx,1]  -= derivcfai1
                    gradinterccfa[xnj+1,idx,2]  -= derivcfaj1
                    gradinterccfa[xni,idx,1] -=  derivcfai
                    gradinterccfa[xnj,idx,2] -=  derivcfaj    

                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sigij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sigij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sigij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sigij[2])
    
                    gradintercefa[xni+1,idx,1] += (derivNormi1/(-sigij[1]))
                    gradintercefa[xnj+1,idx,2] += (derivNormj1/(-sigij[2]))
                    gradintercefa[xni,idx,1] += (derivNormi/(-sigij[1]))
                    gradintercefa[xnj,idx,2] += (derivNormj/(-sigij[2]))
    
                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)
    
                    
                    gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sigij13)*(sigij[2])))
                    gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sigij23)*(sigij[1])))
    
                    gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij2[1]-t2i1))/(2*(sigij13))
                    gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij2[2]-t2j1))/(2*(sigij23))
                    gradpsiefa[xni,idx,1] += (derivNormi*(muij2[1]-t2i))/(2*(sigij13))
                    gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij2[2]-t2j))/(2*(sigij23))
                    
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sigij13*sigij[2])+L2jk/((sigij[1]*sigij[2])))#((L2j*Psi[i,i])/(sigij13*sigij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sigij23*sigij[1])+L2ik/((sigij[2]*sigij[1])))#((L2[i,:].*Psi[j,j])./(sigij23*sigij[1]))
    
                        gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(((muij2[1].-t2i1).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(((muij2[2].-t2j1).*L2jk)./((sigij23)))
                        gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(((muij2[1].-t2i).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(((muij2[2].-t2j).*L2jk)./((sigij23)))
                    end 
                
                end
            end
        end
    
        vecpair = vec(pairsindex)    
        likshared = -sum(gradtaucfa,dims=(2,3))[4:(c)] -sum(gradtauefa,dims=(2,3))[4:(c)]
        
    
        llamb1 = permutedims(gradlambci,(2,3,1)); llam1 = zeros(p); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
        likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
        likinterc1 = zeros(p); sumgradinter1 = sum(gradinterccfa[2:c,:,:],dims=1)[1,:,:]
        likinterc2 = zeros(p); sumgradinter2 = sum(gradintercefa[2:c,:,:],dims=1)[1,:,:]
        likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        @inbounds for i in 1:p
            llam1[i] = -sum(vec(llamb1[ids[i],:,:]')[vecpair.==i])
            likinterc1[i] = -sum(vec(sumgradinter1')[vecpair.==i])
            likinterc2[i] = -sum(vec(sumgradinter2')[vecpair.==i])
            likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
            likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end
    
        likpsitrue = likpsi.+likpsiidx
        liklstrue = likls.+liklsidx
    
        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]
    
        G .= vec(vcat(likshared,llam1,likphida,likinterc1,liklstrue[indices],likpsitrue,likinterc2))
      
    end
    if F !== nothing 
        plimix=zeros(Float64,n);
        tfa1=t1[X.+1]
        tfa=t1[X]
        
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
            Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
            Tt1 = T[:,1,idx]
            
            # EFA
            Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tfa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            Tt2 = T[:,2,idx]
    
            plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
        end
        f=-sum(plimix)
        return f
    end
end



function M_parsMixonlyth01intercCFA(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)
    trsi = (c-1-2) #attenzione p è già qui
 
    # PML - CFA & EFA
    theta_init = theta[1:(size(theta,1)-ncov)]

    lb=vec(vcat(fill(-5.1,trsi),fill(-1,p),fill(-1.,Int((q*(q-1))/2)),fill(-5,p),fill(-Inf,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-5,p)));#1.0e-5
    ub=vec(vcat(fill(5.1,trsi),fill(1.,p),fill(1.,Int((q*(q-1))/2)),fill(5.,p),fill(Inf,Int(p*K-K*(K-1)/2)),fill(Inf,p),fill(10.,p)));    
    par = zeros(Float64,size(lb))
    if it == 1
        par = optimize(Optim.only_fg!((F,G,x)->fgonlyth01intercCFA!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06,time_limit = 300.0)) #20 min
    elseif it > 1
        par = optimize(Optim.only_fg!((F,G,x)->fgonlyth01intercCFA!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06,time_limit = 1200.0)) #20 min
    end 

    
    # mixture parameter
    betas = theta[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
    Ez1 = T[:,1,:]
    Ez0 = T[:,2,:]
    
    pairs = size(pairsindex,2)
    betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)
    
    
    # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
   
    # println(par.minimizer,pg[1])
    if !(Optim.converged(par)) && isnan(betas[1]) 
        println(betas[1])
        if start == true
            return zeros(Float64,size(theta,1))
        elseif start == false
            return Optim.converged(par)
        end
    end
    return vcat(par.minimizer,betasnew)#,pg[1])
end







########################################## intercetta CFA + Fixed and invariant thresholds con 0 ##################################################
function startvalsonlyth0intercCFA(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
    local tau1
    liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1-1)+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+(p)+p)+ncov,replix);#p-1
    for i in 1:replix
        println(i)

        # if c==4 tau1 = rand(Uniform(0.2,.7)).+repeat(collect(range(1,3,c-1-1))) else tau1 = rand(Uniform(1,1.2)).+repeat(collect(range(0.5,1.5,c-1-1))) end
        tau1 = rand(Uniform(1,1.2)).+repeat(collect(range(0.5,1.5,c-1-1)))
        lambda1 = rand(Uniform(0.5,.90),p);
        phi=rand(LKJ(q,1))
        phi=phi[tril!(trues(size(phi)), -1)]
        interc1=rand(Uniform(0.2,2),p)
        lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(0.15,.90),Int((p*K-(K*(K-1)/2))))))
        lambda2=lamb2[lamb2.!=0.]
        psi=rand(Uniform(0.5,1.25),p)
        interc2=rand(Uniform(0.2,2),p)



        prop = rand(Uniform(0.65,0.95))
        bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.03,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


        param_initials = vcat(tau1,lambda1,phi,interc1,lambda2,psi,interc2,bets)

        # first iteration
        plimix = E_stepMixonlyth0intercCFA(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        params[:,i] = M_parsMixonlyth0intercCFA(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        dt = time() - t
        prrot=params[:,i]
        liksmix[i] = E_stepMixonlyth0intercCFA(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
    end
    param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    lik=liksmix[findall(liksmix.==maximum(liksmix))]

    tau1=param_initials0[(1):((c-1-1))]
    tau1=vec(sort(reshape(tau1,1,c-1-1),dims=2))
    lambda1=param_initials0[((c-1-1)+1):((c-1-1)+p)]
    phi=param_initials0[((c-1-1)+p+1):((c-1-1)+p+Int(q*(q-1)/2))]
    lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    interc1=param_initials0[((c-1-1)+p+Int(q*(q-1)/2)+1):((c-1-1)+p+Int(q*(q-1)/2)+p)]#interc=param_initials0[((c-1-1)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1-1)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    lambda2=param_initials0[((c-1-1)+p+Int(q*(q-1)/2)+p+1):((c-1-1)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2))))]
    psi=param_initials0[((c-1-1)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+1):((c-1-1)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p))]# psi=param_initials0[((c-1-1)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+1):((c-1-1)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p-1))]
    interc2=param_initials0[((c-1-1)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p)+1):((c-1-1)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p)+p)]#interc=param_initials0[((c-1-1)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1-1)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
    bets = param_initials0[((c-1-1)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p)+p+1):((c-1-1)+p+Int(q*(q-1)/2)+p+Int((p*K-(K*(K-1)/2)))+(p)+p+ncov)]#    bets = param_initials0[((c-1-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1-1)+p+Int(q*(q-1)/2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


    param_initialii = vcat(tau1,lambda1,phi,interc1,lambda2,psi,interc2,bets)

    println(param_initialii,liksmix,"first iteration: done")
    return param_initialii
end



function EMparsMixonlyth0intercCFA(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
    local plimix::Array, T,oldlik
    n = size(X,1)
    er=0;
    # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
    it=0;
    dev=Inf; dev1=Inf; dev2=Inf;

    
    oldlik=-1000000.;
    oldlik2=1000000.;
    oldlikA=-1000000.;
    lik=0.
    

    noerror = true; devnostop = true;  sum_devs = 0;
    while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
        it = it+1;
        param=param01;

        T00 = 1
        plimix = E_stepMixonlyth0intercCFA(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        param01 = M_parsMixonlyth0intercCFA(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
        dt = time() - t
        if param01 == false
            return return(vcat(param01,vec(plimix),vec(T),0.,0.))
        end
        lik = E_stepMixonlyth0intercCFA(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);

        #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
        con = (lik-oldlik)/(oldlik-oldlik2);
        likA = oldlik + (lik-oldlik)/(1-con)
        old_dev = dev1
        dev1 = abs(likA-oldlikA)
        dev2 = abs(lik-oldlik)
        println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
        println(param01)
        dev = abs(lik-oldlik)/abs(oldlik)

        noerror = oldlik<lik

        oldlikA=likA
        oldlik2=oldlik
        oldlik=lik;
        pml_em_log("EMparsMixonlyth0intercCFA", it, lik, dev, dev1, abs(lik - oldlik2), tol, noerror)

        # sum_devs += dev > old_dev
        # devstop = sum_devs <= 3
    end
    
    # classificazioni dopo ultima massimizzazione
    plimix = E_stepMixonlyth0intercCFA(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

    if noerror == true && it<=maxit 
        return vcat(param01,vec(plimix),vec(T),lik,1.)
    else
        return vcat(param01,vec(plimix),vec(T),0.,0.)
    end
end

function E_stepMixonlyth0intercCFA(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
    if sum(param .== 0.) == size(param,1)
        return -Inf
    end 
    trsi = (c-1-1)
    t1s = param[1:trsi]
    L1 = param[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
    interc1 = param[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+p)]

    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi+p+(q*(q-1))/2+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(param[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = param[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+p)]
    # prop = param[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+p+1)]
    betas = param[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+(p)+p+ncov)]  #  betas = param[Int(trsi+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    t1=vcat(-Inf,0,t1s,Inf)

    Sigms=ones(p,p,2);
    Sigms[:,:,1]=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sigms[:,:,2]=L2*L2'+Psi 

    plimix1=zeros(Float64,n);plimix2=zeros(Float64,n);
    tfa1=t1[X.+1]
    tfa=t1[X]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)
    
    if estep == true
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));
    
        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(interc1)[pairsindex[:,idx]]')
                    low=(tfa[:,pairsindex[:,idx]].-(interc1)[pairsindex[:,idx]]')
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
                    low=(tfa[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3.*pg[:,g];
            end
        end
    # % probabilities in each row must sum to 1 (application of Bayes law)
        return plimix
    elseif estep == false
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));

        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(interc1)[pairsindex[:,idx]]')
                    low=(tfa[:,pairsindex[:,idx]].-(interc1)[pairsindex[:,idx]]')
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    r2c=cov2cor(s2c);
                    up=(tfa1[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                    low=(tfa[:,pairsindex[:,idx]].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3
            end
        end
        lik1=(log.(plimix));
        lik1=sum(lik1,dims=3);
        U=exp.(lik1)[:,:,1].*pg; 
        cllik=sum(log.(sum(U,dims=2)))
        return cllik
    end
end
    

function pmllastonlyth0intercCFA(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    trsi = (c-1-1)
    t1s = x[1:trsi]
    L1 = x[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
    interc1 = x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+p)]

    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = x[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+p)]
   
    Sig1=L1*Phi*L1' + Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi

   
    t1=vcat(-Inf,0,t1s,Inf)
    tfa1=t1[X.+1]
    tfa=t1[X]
    
    plimix=zeros(Float64,n);
    for idx=1:size(pairsindex,2) #parallelizzabile
        # CFA
        Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
        Bcfa=hcat(tfa[:,pairsindex[:,idx]].-interc1[pairsindex[:,idx]]',tfa1[:,pairsindex[:,idx]].-interc1[pairsindex[:,idx]]');
        tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
        Tt1 = T[:,1,idx]
        
        #EFA
        Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
        Refa=cov2cor(Sefa);
        Befa=hcat((tfa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
        tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
        Tt2 = T[:,2,idx]

        plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
    end
    f=-sum(plimix)
return f
end




function fgonlyth0intercCFA!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trsi = (c-1-1)
    t1s = x[1:trsi]
    L1 = x[(trsi+1):(trsi+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi+p+1):Int(trsi+p+(q*(q-1))/2)])
    interc1 = x[Int(trsi+p+(q*(q-1))/2+1):Int(trsi+p+(q*(q-1))/2+p)]

    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi+p+(q*(q-1))/2+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2)))] )
    Psi = diagm(x[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+1):Int((trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = x[Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+1):Int(trsi+p+(q*(q-1))/2+p+(p*K-(K*(K-1)/2))+p+p)]
   
    t1=vcat(-Inf,0,t1s,Inf)
 

    Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi


    if G !== nothing
        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradinterccfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradintercefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)
    
        gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)
    
        gradlambci = zeros(Float64,2,q,size(pairsindex,2))
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
        
        idx = 0
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1
                rho1 = Sig1[i,j]
                cormat1 = ones(Float64,2,2); cormat1[1,2]=rho1; cormat1[2,1]=rho1
                sigs = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sigs)
                rho2 = cormat2[1,2]
                muij1 = interc1[vcat(i,j)]
                muij2 = interc2[vcat(i,j)]
                sigsij = diag(sigs)
                sigij = sqrt.(sigsij)
                sigij13 = sigij[1]^(3)
                sigij23 = sigij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]
                @inbounds for N in 1:n
                    xni = X[N,i]
                    xnj = X[N,j]
                    t1i = t1[xni]
                    t1i1 = t1[xni+1]
                    t1j = t1[xnj]
                    t1j1 = t1[xnj+1]
                    t1ist = (t1i-muij1[1])
                    t1i1st = (t1i1-muij1[1])
                    t1jst = (t1j-muij1[2])
                    t1j1st = (t1j1-muij1[2])
                    T1 = T[N,1,idx][1]
                    T2 = T[N,2,idx][1]
                    t2i = t1[xni]
                    t2i1 = t1[xni+1]
                    t2j = t1[xnj]
                    t2j1 = t1[xnj+1]
                    t2ist = (t2i-muij2[1])/sigij[1]
                    t2i1st = (t2i1-muij2[1])/sigij[1]
                    t2jst = (t2j-muij2[2])/sigij[2]
                    t2j1st = (t2j1-muij2[2])/sigij[2]
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                   
                    pijicfa = bvnP(t1ist, t1i1st, t1jst, t1j1st, rho1, false)
                    derivcfai1 = (T1/pijicfa)*(pdf(Normal(0,1),t1i1st)*(cdf(Normal(0,1),(t1j1st-rho1*t1i1st)/sqrtrij1)-cdf(Normal(0,1),(t1jst-rho1*t1i1st)/sqrtrij1))) 
                    derivcfaj1 = (T1/pijicfa)*(pdf(Normal(0,1),t1j1st)*(cdf(Normal(0,1),(t1i1st-rho1*t1j1st)/sqrtrij1)-cdf(Normal(0,1),(t1ist-rho1*t1j1st)/sqrtrij1))) 
                    derivcfai =  (T1/pijicfa)*(pdf(Normal(0,1),t1ist)*(cdf(Normal(0,1),(t1jst-rho1*t1ist)/sqrtrij1)-cdf(Normal(0,1),(t1j1st-rho1*t1ist)/sqrtrij1))) 
                    derivcfaj = (T1/pijicfa)*(pdf(Normal(0,1),t1jst)*(cdf(Normal(0,1),(t1ist-rho1*t1jst)/sqrtrij1)-cdf(Normal(0,1),(t1i1st-rho1*t1jst)/sqrtrij1))) 
                    
                    gradtaucfa[xni+1,idx,1] += derivcfai1
                    gradtaucfa[xnj+1,idx,2] += derivcfaj1
                    gradtaucfa[xni,idx,1] += derivcfai
                    gradtaucfa[xnj,idx,2] += derivcfaj

                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))
    
    
                    
                    if !(t1ist == -Inf && t1jst == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1jst))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1st == Inf && t1j1st == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1j1st))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end
        
                    if !(t1i1st == Inf && t1j == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1st,t1jst))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end
        
                    if !(t1i == -Inf && t1j1 == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1ist,t1j1st))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    gradlambci[1:2,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)
                    gradphici[1:q,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*(L1[i,:]*L1[j,:]')
                    
                    gradinterccfa[xni+1,idx,1]  -= derivcfai1
                    gradinterccfa[xnj+1,idx,2]  -= derivcfaj1
                    gradinterccfa[xni,idx,1] -=  derivcfai
                    gradinterccfa[xnj,idx,2] -=  derivcfaj    

                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sigij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sigij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sigij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sigij[2])
    
                    gradintercefa[xni+1,idx,1] += (derivNormi1/(-sigij[1]))
                    gradintercefa[xnj+1,idx,2] += (derivNormj1/(-sigij[2]))
                    gradintercefa[xni,idx,1] += (derivNormi/(-sigij[1]))
                    gradintercefa[xnj,idx,2] += (derivNormj/(-sigij[2]))
    
                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)
    
                    
                    gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sigij13)*(sigij[2])))
                    gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sigij23)*(sigij[1])))
    
                    gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij2[1]-t2i1))/(2*(sigij13))
                    gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij2[2]-t2j1))/(2*(sigij23))
                    gradpsiefa[xni,idx,1] += (derivNormi*(muij2[1]-t2i))/(2*(sigij13))
                    gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij2[2]-t2j))/(2*(sigij23))
                    
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sigij13*sigij[2])+L2jk/((sigij[1]*sigij[2])))#((L2j*Psi[i,i])/(sigij13*sigij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sigij23*sigij[1])+L2ik/((sigij[2]*sigij[1])))#((L2[i,:].*Psi[j,j])./(sigij23*sigij[1]))
    
                        gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(((muij2[1].-t2i1).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(((muij2[2].-t2j1).*L2jk)./((sigij23)))
                        gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(((muij2[1].-t2i).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(((muij2[2].-t2j).*L2jk)./((sigij23)))
                    end 
                
                end
            end
        end
    
        vecpair = vec(pairsindex)    
        likshared = -sum(gradtaucfa,dims=(2,3))[3:(c)] -sum(gradtauefa,dims=(2,3))[3:(c)]
        
    
        llamb1 = permutedims(gradlambci,(2,3,1)); llam1 = zeros(p); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
        likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
        likinterc1 = zeros(p); sumgradinter1 = sum(gradinterccfa[2:c,:,:],dims=1)[1,:,:]
        likinterc2 = zeros(p); sumgradinter2 = sum(gradintercefa[2:c,:,:],dims=1)[1,:,:]
        likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        @inbounds for i in 1:p
            llam1[i] = -sum(vec(llamb1[ids[i],:,:]')[vecpair.==i])
            likinterc1[i] = -sum(vec(sumgradinter1')[vecpair.==i])
            likinterc2[i] = -sum(vec(sumgradinter2')[vecpair.==i])
            likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
            likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end
    
        likpsitrue = likpsi.+likpsiidx
        liklstrue = likls.+liklsidx
    
        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]
    
        G .= vec(vcat(likshared,llam1,likphida,likinterc1,liklstrue[indices],likpsitrue,likinterc2))
      
    end
    if F !== nothing 
        plimix=zeros(Float64,n);
        tfa1=t1[X.+1]
        tfa=t1[X]
        
        for idx=1:size(pairsindex,2) #parallelizzabile
            # CFA
            Scfa=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
            Bcfa=hcat(tfa[:,pairsindex[:,idx]],tfa1[:,pairsindex[:,idx]]);
            tcfa3=bvnP.(Bcfa[:,1],Bcfa[:,3],Bcfa[:,2],Bcfa[:,4],Scfa[1,2],false);
            Tt1 = T[:,1,idx]
            
            # EFA
            Sefa=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
            Refa=cov2cor(Sefa);
            Befa=hcat((tfa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'),(tfa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(Sefa)'));
            tefa3=bvnP.(Befa[:,1],Befa[:,3],Befa[:,2],Befa[:,4],Refa[1,2],false);
            Tt2 = T[:,2,idx]
    
            plimix += log.(tcfa3).*Tt1 + log.(tefa3).*Tt2; 
        end
        f=-sum(plimix)
        return f
    end
end



function M_parsMixonlyth0intercCFA(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)
    trsi = (c-1-1) #attenzione p è già qui
 
    # PML - CFA & EFA
    theta_init = theta[1:(size(theta,1)-ncov)]

    lb=vec(vcat(fill(-5.1,trsi),fill(-1,p),fill(-1.,Int((q*(q-1))/2)),fill(-5,p),fill(-Inf,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-5,p)));#1.0e-5
    ub=vec(vcat(fill(5.1,trsi),fill(1.,p),fill(1.,Int((q*(q-1))/2)),fill(5.,p),fill(Inf,Int(p*K-K*(K-1)/2)),fill(Inf,p),fill(10.,p)));    
    par = zeros(Float64,size(lb))
    if it == 1
        par = optimize(Optim.only_fg!((F,G,x)->fgonlyth0intercCFA!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06,time_limit = 300.0)) #20 min
    elseif it > 1
        par = optimize(Optim.only_fg!((F,G,x)->fgonlyth0intercCFA!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06,time_limit = 1200.0)) #20 min
    end 

    
    # mixture parameter
    betas = theta[Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi+p+(q*(q-1))/2+(p*K-(K*(K-1)/2))+(p)+p+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
    Ez1 = T[:,1,:]
    Ez0 = T[:,2,:]
    
    pairs = size(pairsindex,2)
    betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)
    
    
    # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
   
    # println(par.minimizer,pg[1])
    if !(Optim.converged(par)) && isnan(betas[1]) 
        println(betas[1])
        if start == true
            return zeros(Float64,size(theta,1))
        elseif start == false
            return Optim.converged(par)
        end
    end
    return vcat(par.minimizer,betasnew)#,pg[1])
end













########################################################### Fixed Thresholds only for EFA ##################################################################



function EMparsMixthefa(param01,tol,X,p,q,K,c,ics,ncov,maxit,L_str,indices,pairsindex)
    local plimix::Array, T,oldlik
    n = size(X,1)
    er=0;
    # if P>3 pairsindex=pairsindex[end:-1:1,:]; end 
    it=0;
    dev=Inf; dev1=Inf; dev2=Inf;

    
    oldlik=-1000000.;
    oldlik2=1000000.;
    oldlikA=-1000000.;
    lik=0.
    

    noerror = true; devnostop = true;  sum_devs = 0;
    while (dev>tol) && (dev1>tol)  && (dev2>tol) && it <= maxit && noerror #&& devnostop
        
        it = it+1;
        param=param01;

        T00 = 1
        plimix = E_stepMixthefa(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        param01 = M_parsMixthefa(T,param,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,t,plimix,false); 
        dt = time() - t
        if param01 == false
            return return(vcat(param01,vec(plimix),vec(T),0.,0.))
        end
        lik = E_stepMixthefa(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);

        #Aitken acceleration based stopping criterion p.142 EM algorithm and extensions
        con = (lik-oldlik)/(oldlik-oldlik2);
        likA = oldlik + (lik-oldlik)/(1-con)
        old_dev = dev1
        dev1 = abs(likA-oldlikA)
        dev2 = abs(lik-oldlik)
        println(lik," ",dev," ",dev2," ",dev1," ",oldlik<lik," ",dt)
        println(param01)
        dev = abs(lik-oldlik)/abs(oldlik)

        noerror = oldlik<lik

        oldlikA=likA
        oldlik2=oldlik
        oldlik=lik;
        pml_em_log("EMparsMixthefa", it, lik, dev, dev1, abs(lik - oldlik2), tol, noerror)

        # sum_devs += dev > old_dev
        # devstop = sum_devs <= 3
    end
    
    # classificazioni dopo ultima massimizzazione
    plimix = E_stepMixthefa(param01,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
    T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;

    if noerror == true && it<=maxit 
        return vcat(param01,vec(plimix),vec(T),lik,1.)
    else
        return vcat(param01,vec(plimix),vec(T),0.,0.)
    end
end

function E_stepMixthefa(param,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,estep,T=false)
    if sum(param .== 0.) == size(param,1)
        return -Inf
    end
    trs = (c-1-2)
    trsi = (c-1)
    t1s = param[1:trsi*p]
    L1 = param[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(param[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])

    t2s = param[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs)]
    L2 = reconstruct_L2_matrix(p, K, indices, param[Int(trsi*p+p+(q*(q-1))/2+trs+1):Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))] )


    Psi = diagm(param[Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = param[Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+1):Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p)]
    # prop = param[Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p+1)]
    betas = param[Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+p+ncov)]  #  betas = param[Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi*p+p+(q*(q-1))/2+trs*p+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]

    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
    br2 = vcat(0.,1.,t2s)
    t2=vcat(-Inf,br2,Inf)


    Sigms=ones(p,p,2);
    Sigms[:,:,1]=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sigms[:,:,2]=L2*L2'+Psi 

    plimix1=zeros(Float64,n);plimix2=zeros(Float64,n);
    tcfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tcfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
    tefa1=t2[X.+1]
    tefa=t2[X]

    propo = exp.(ics*betas)./(exp.(ics*betas).+1)
    
    pg = hcat(propo,1. .-propo)
    
    tuni = cat(tcfa1,tefa1,dims=3)
    tdui = cat(tcfa,tefa,dims=3)
    
    if estep == true
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));
    
        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    # s2c=s2c+0.0001*Matrix(I,2,2)
                    r2c=cov2cor(s2c);
                    up=(tuni[:,pairsindex[:,idx],g])
                    low=(tdui[:,pairsindex[:,idx],g])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    # s2c=s2c+0.0001*Matrix(I,2,2)
                    r2c=cov2cor(s2c);
                    up=(tuni[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                    low=(tdui[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3.*pg[:,g];
            end
        end
    # % probabilities in each row must sum to 1 (application of Bayes law)
        return plimix
    elseif estep == false
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        plimix=ones(n,2,size(pairsindex,2));

        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    # s2c=s2c+0.0001*Matrix(I,2,2)
                    r2c=cov2cor(s2c);
                    up=(tuni[:,pairsindex[:,idx],g])
                    low=(tdui[:,pairsindex[:,idx],g])
                elseif g==2
                    s2c=Sigms[pairsindex[:,idx],pairsindex[:,idx],g]
                    # s2c=s2c+0.0001*Matrix(I,2,2)
                    r2c=cov2cor(s2c);
                    up=(tuni[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                    low=(tdui[:,pairsindex[:,idx],g].-(interc2)[pairsindex[:,idx]]')./sqrt.(diag(s2c)');#./sqrt.(diag(s2c)');
                end
                B=hcat(low,up);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                plimix[:,g,idx]= t3
            end
        end
        lik1=(log.(plimix));
        lik1=sum(lik1,dims=3);
        U=exp.(lik1)[:,:,1].*pg; 
        cllik=sum(log.(sum(U,dims=2)))
        return cllik
    end
end
    

function pmlthefa(x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...

    trs = (c-1-2)
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])

    t2s = x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs)]
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+trs+1):Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))] )
    Psi = diagm(x[Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = x[Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+1):Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p)]
    # prop = x[Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p+1)]

    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
    br2 = vcat(0.,1.,t2s)
    t2=vcat(-Inf,br2,Inf)


    Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi


    plimix1=zeros(Float64,n);plimix2=zeros(Float64,n);
    tcfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
    tcfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
    tefa1=t2[X.+1]
    tefa=t2[X]

    low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
    for idx=1:size(pairsindex,2) #parallelizzabile
        for g=1:2
            if g==1
                s2c=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
                # s2c=s2c+0.0001*Matrix(I,2,2)
                r2c=cov2cor(s2c);
                B=hcat(tcfa[:,pairsindex[:,idx]],tcfa1[:,pairsindex[:,idx]]);
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                Tt1 = T[:,1,idx]
                plimix1 += log.(t3).*Tt1;
            else 
                s2c=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
                # s2c=s2c+0.0001*Matrix(I,2,2)
                r2c=cov2cor(s2c);
                B=hcat((tefa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(s2c)'),(tefa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(s2c)'));
                t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                Tt2 = T[:,2,idx]
                plimix2 += log.(t3).*Tt2;
            end
            
        end
    end
    f=-sum(plimix1+plimix2)
end

function fgthefa!(F, G, x,X,T,p,q,K,c,n,pairsindex,L_str,indices)
    # do common computations here
    # ...
    trs = (c-1-2)
    trsi = (c-1)
    t1s = x[1:trsi*p]
    L1 = x[(trsi*p+1):(trsi*p+p)].*L_str 
    Phi = reconstruct_correlation_matrix(x[Int(trsi*p+p+1):Int(trsi*p+p+(q*(q-1))/2)])
    
    t2s = x[Int(trsi*p+p+(q*(q-1))/2+1):Int(trsi*p+p+(q*(q-1))/2+trs)]
    L2 = reconstruct_L2_matrix(p, K, indices, x[Int(trsi*p+p+(q*(q-1))/2+trs+1):Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))] )
    Psi = diagm(x[Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+1):Int((trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2)))+p)])#K
    interc2 = x[Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+1):Int(trsi*p+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p)]
   
    # prop = x[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+p+p+1)]

    br1 = reshape(t1s,p,trsi)
    t1=hcat(fill(-Inf,p),br1,fill(Inf,p))
    br2 = vcat(0.,1.,t2s)
    t2=vcat(-Inf,br2,Inf)


    Sig1=L1*Phi*L1'+ Diagonal(1 .-Diagonal((L1*Phi*L1')))
    Sig2=L2*L2'+Psi


    if G !== nothing
        gradtaucfa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradtauefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradintercefa = zeros(Float64,c+1,size(pairsindex,2),2)
        gradpsiefa = zeros(Float64,c+1,size(pairsindex,2),2)
    
        gradpsiefaidx = zeros(Float64,size(pairsindex,2),2)
        gradlsefaidx = zeros(Float64,K,size(pairsindex,2),2)
    
        gradlambci = zeros(Float64,2,q,size(pairsindex,2))
        gradphici = zeros(Float64,q,q,size(pairsindex,2))
        gradlamb2ci = zeros(Float64,K,c+1,size(pairsindex,2),2)
        
        idx = 0
        @inbounds for i in 1:(p-1)
            @inbounds for j in (i+1):p
                idx += 1
                rho1 = Sig1[i,j]
                cormat1 = ones(Float64,2,2); cormat1[1,2]=rho1; cormat1[2,1]=rho1
                sigs = hcat(vcat(Sig2[i,i],Sig2[i,j]),vcat(Sig2[j,i],Sig2[j,j]))
                cormat2 = cov2cor(sigs)
                rho2 = cormat2[1,2]
                muij = interc2[vcat(i,j)]
                sigsij = diag(sigs)
                sigij = sqrt.(sigsij)
                sigij13 = sigij[1]^(3)
                sigij23 = sigij[2]^(3)
                L2i=L2[i,:]
                L2j=L2[j,:]
                @inbounds for N in 1:n
                    xni = X[N,i]
                    xnj = X[N,j]
                    t1i = t1[i,xni]
                    t1i1 = t1[i,xni+1]
                    t1j = t1[j,xnj]
                    t1j1 = t1[j,xnj+1]
                    T1 = T[N,1,idx][1]
                    T2 = T[N,2,idx][1]
                    t2i = t2[xni]
                    t2i1 = t2[xni+1]
                    t2j = t2[xnj]
                    t2j1 = t2[xnj+1]
                    t2ist = (t2i-muij[1])/sigij[1]
                    t2i1st = (t2i1-muij[1])/sigij[1]
                    t2jst = (t2j-muij[2])/sigij[2]
                    t2j1st = (t2j1-muij[2])/sigij[2]
                    sqrtrij1 = sqrt(1-rho1^2)
                    sqrtrij2 = sqrt(1-rho2^2)
                                   
                    pijicfa = bvnP(t1i, t1i1, t1j, t1j1, rho1, false)
                    gradtaucfa[xni+1,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i1)*(cdf(Normal(0,1),(t1j1-rho1*t1i1)/sqrtrij1)-cdf(Normal(0,1),(t1j-rho1*t1i1)/sqrtrij1))) 
                    gradtaucfa[xni,idx,1] += (T1/pijicfa)*(pdf(Normal(0,1),t1i)*(cdf(Normal(0,1),(t1j-rho1*t1i)/sqrtrij1)-cdf(Normal(0,1),(t1j1-rho1*t1i)/sqrtrij1))) 
                    gradtaucfa[xnj+1,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j1)*(cdf(Normal(0,1),(t1i1-rho1*t1j1)/sqrtrij1)-cdf(Normal(0,1),(t1i-rho1*t1j1)/sqrtrij1))) 
                    gradtaucfa[xnj,idx,2] += (T1/pijicfa)*(pdf(Normal(0,1),t1j)*(cdf(Normal(0,1),(t1i-rho1*t1j)/sqrtrij1)-cdf(Normal(0,1),(t1i1-rho1*t1j)/sqrtrij1))) 
                    
    
                    pijiefa = bvnP((t2ist), (t2i1st), (t2jst), (t2j1st), rho2, false)
                    derivNormi1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2i1st))*(cdf(Normal(0,1),((t2j1st)-rho2*(t2i1st))/sqrtrij2)-cdf(Normal(0,1),((t2jst)-rho2*(t2i1st))/sqrtrij2)))
                    derivNormj1 = (T2/pijiefa)*(pdf(Normal(0,1),(t2j1st))*(cdf(Normal(0,1),((t2i1st)-rho2*(t2j1st))/sqrtrij2)-cdf(Normal(0,1),((t2ist)-rho2*(t2j1st))/sqrtrij2)))
                    derivNormi = (T2/pijiefa)*(pdf(Normal(0,1),(t2ist))*(cdf(Normal(0,1),((t2jst)-rho2*(t2ist))/sqrtrij2)-cdf(Normal(0,1),((t2j1st)-rho2*(t2ist))/sqrtrij2)))
                    derivNormj = (T2/pijiefa)*(pdf(Normal(0,1),(t2jst))*(cdf(Normal(0,1),((t2ist)-rho2*(t2jst))/sqrtrij2)-cdf(Normal(0,1),((t2i1st)-rho2*(t2jst))/sqrtrij2)))
    
    
                    
                    if !(t1i == -Inf && t1j == -Inf) 
                        xnixnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j))
                        xnixnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2jst)))
                    else
                        xnixnjc = 0.0
                        xnixnje = 0.0
                    end
                        
                    if !(t1i1 == Inf && t1j1 == Inf) 
                        xni1xnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j1))
                        xni1xnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2j1st)))
                    else
                        xni1xnj1c = 0.0
                        xni1xnj1e = 0.0
                    end
        
                    if !(t1i1 == Inf && t1j == -Inf) 
                        xni1xnjc = pdf(MvNormal(zeros(2),cormat1),vcat(t1i1,t1j))
                        xni1xnje = pdf(MvNormal(zeros(2),cormat2),vcat((t2i1st),(t2jst)))
                    else
                        xni1xnjc = 0.0
                        xni1xnje = 0.0
                    end
        
                    if !(t1i == -Inf && t1j1 == Inf) 
                        xnixnj1c = pdf(MvNormal(zeros(2),cormat1),vcat(t1i,t1j1))
                        xnixnj1e = pdf(MvNormal(zeros(2),cormat2),vcat((t2ist),(t2j1st)))
                    else
                        xnixnj1c = 0.0
                        xnixnj1e = 0.0
                    end
                    gradlambci[1:2,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*vcat(L1[j,:]'*Phi,L1[i,:]'*Phi)
                    gradphici[1:q,1:q,idx] .+= (T1/pijicfa)*(xni1xnj1c-xni1xnjc-xnixnj1c+xnixnjc)*(L1[i,:]*L1[j,:]')
                    
                    #CHAIN RULE PER TRASFORMAZIONE DI THRESHOLDS
                    gradtauefa[xni+1,idx,1] += (derivNormi1)/(sigij[1])
                    gradtauefa[xnj+1,idx,2] += (derivNormj1)/(sigij[2])
                    gradtauefa[xni,idx,1] += (derivNormi)/(sigij[1])
                    gradtauefa[xnj,idx,2] += (derivNormj)/(sigij[2])
    
                    gradintercefa[xni+1,idx,1] += (derivNormi1/(-sigij[1]))
                    gradintercefa[xnj+1,idx,2] += (derivNormj1/(-sigij[2]))
                    gradintercefa[xni,idx,1] += (derivNormi/(-sigij[1]))
                    gradintercefa[xnj,idx,2] += (derivNormj/(-sigij[2]))
    
                    dpidrij=(T2/pijiefa)*(xni1xnj1e-xni1xnje-xnixnj1e+xnixnje)
    
                    
                    gradpsiefaidx[idx,1] += dpidrij*(-(L2i'*L2j)/(2*(sigij13)*(sigij[2])))
                    gradpsiefaidx[idx,2] += dpidrij*(-(L2j'*L2i)/(2*(sigij23)*(sigij[1])))
    
                    gradpsiefa[xni+1,idx,1] += (derivNormi1*(muij[1]-t2i1))/(2*(sigij13))
                    gradpsiefa[xnj+1,idx,2] +=  (derivNormj1*(muij[2]-t2j1))/(2*(sigij23))
                    gradpsiefa[xni,idx,1] += (derivNormi*(muij[1]-t2i))/(2*(sigij13))
                    gradpsiefa[xnj,idx,2] +=  (derivNormj*(muij[2]-t2j))/(2*(sigij23))
                    
                    @inbounds for k in 1:K
                        L2ik=L2[i,k]
                        L2jk=L2[j,k]
                        gradlsefaidx[k,idx,1] += dpidrij*(-L2ik*(L2i'*L2j)/(sigij13*sigij[2])+L2jk/((sigij[1]*sigij[2])))#((L2j*Psi[i,i])/(sigij13*sigij[2]))
                        gradlsefaidx[k,idx,2] += dpidrij*(-L2jk*(L2j'*L2i)/(sigij23*sigij[1])+L2ik/((sigij[2]*sigij[1])))#((L2[i,:].*Psi[j,j])./(sigij23*sigij[1]))
    
                        gradlamb2ci[k,xni+1,idx,1] +=  derivNormi1 *(((muij[1].-t2i1).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj+1,idx,2] +=  derivNormj1 *(((muij[2].-t2j1).*L2jk)./((sigij23)))
                        gradlamb2ci[k,xni,idx,1]   +=  derivNormi  *(((muij[1].-t2i).*L2ik)./((sigij13)))
                        gradlamb2ci[k,xnj,idx,2]   +=  derivNormj  *(((muij[2].-t2j).*L2jk)./((sigij23)))
                    end 
                
                end
            end
        end
    
    
        vecpair = vec(pairsindex)
        liktau = zeros(p*(c-1));idpcfa = 0   
        @inbounds for j in 2:c
            @inbounds for i in 1:p
                idpcfa += 1
                liktau[idpcfa] = -sum(vec(vcat(gradtaucfa[j,:,1]',gradtaucfa[j,:,2]'))[vecpair.==i])
            end
        end
    
        
        liktau2 = -sum(gradtauefa,dims=(2,3))[4:(c)]
        
    
        llamb1 = permutedims(gradlambci,(2,3,1)); llam1 = zeros(p); ids = find_ones_columns(L_str);
        likls = zeros(p,K); sumgradls = sum(gradlamb2ci[:,2:c,:,:],dims=2)[:,1,:,:] #ATTENZIONE PER K>1 BISOGNA RIVEDERE QUELL'1 IN [1,2:c,...]#sum(gradlsefa[2:c,:,:],dims=1)[1,:,:]
        likpsi = zeros(p); sumgradpsi = sum(gradpsiefa[2:c,:,:],dims=1)[1,:,:]
        likinterc2 = zeros(p); sumgradinter = sum(gradintercefa[2:c,:,:],dims=1)[1,:,:]
        likpsiidx  = zeros(p); liklsidx  = zeros(p,K)#Int((p*K-(K*(K-1)/2))));
        @inbounds for i in 1:p
            llam1[i] = -sum(vec(llamb1[ids[i],:,:]')[vecpair.==i])
            likinterc2[i] = -sum(vec(sumgradinter')[vecpair.==i])
            likpsi[i] = -sum(vec(sumgradpsi')[vecpair.==i])
            likpsiidx[i] = -sum(vec(gradpsiefaidx')[vecpair.==i])
            @inbounds for k in 1:K
                likls[i,k] = -sum(vec(sumgradls[k,:,:]')[vecpair.==i]) 
                liklsidx[i,k] = -sum(vec(gradlsefaidx[k,:,:]')[vecpair.==i]) 
            end
        end
    
        likpsitrue = likpsi.+likpsiidx
        liklstrue = likls.+liklsidx
    
        likp = sum(gradphici,dims=(3))[:,:,1]
        likphida = -likp[triu!(trues(size(likp)), 1)]
    
        G .= vec(vcat(liktau,llam1,likphida,liktau2,liklstrue[indices],likpsitrue,likinterc2))
      
    end
    if F !== nothing 
  
    
        plimix1=zeros(Float64,n);plimix2=zeros(Float64,n);
        tcfa1=[t1[j,X[i,j]+1] for i in 1:n, j in 1:p]
        tcfa=[t1[j,X[i,j]] for i in 1:n, j in 1:p]
        tefa1=t2[X.+1]
        tefa=t2[X]
       
        low = Matrix{Float64}(undef,n,2); up = Matrix{Float64}(undef,n,2)
        for idx=1:size(pairsindex,2) #parallelizzabile
            for g=1:2
                if g==1
                    s2c=Sig1[pairsindex[:,idx],pairsindex[:,idx]]
                    # s2c=s2c+0.0001*Matrix(I,2,2)
                    r2c=cov2cor(s2c);
                    B=hcat(tcfa[:,pairsindex[:,idx]],tcfa1[:,pairsindex[:,idx]]);
                    t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                    Tt1 = T[:,1,idx]
                    plimix1 += log.(t3).*Tt1;
                else 
                    s2c=Sig2[pairsindex[:,idx],pairsindex[:,idx]]
                    # s2c=s2c+0.0001*Matrix(I,2,2)
                    r2c=cov2cor(s2c);
                    B=hcat((tefa[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(s2c)'),(tefa1[:,pairsindex[:,idx]].-interc2[pairsindex[:,idx]]')./sqrt.(diag(s2c)'));
                    t3=bvnP.(B[:,1],B[:,3],B[:,2],B[:,4],r2c[1,2],false);
                    Tt2 = T[:,2,idx]
                    plimix2 += log.(t3).*Tt2;
                end
                
            end
        end
        f=-sum(plimix1+plimix2)
        # @info "Current objective: $f"
        return f
    end
end



function M_parsMixthefa(T,theta,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,it,time,plimix,start)
    trsi = (c-1)*p #attenzione p è già qui
    trs = (c-1-2)

    # PML - CFA & EFA
    theta_init = theta[1:(size(theta,1)-ncov)]

    lb=vec(vcat(fill(-5.1,trsi),fill(-1,p),fill(-1.,Int((q*(q-1))/2)),fill(1.,trs),fill(-Inf,Int(p*K-K*(K-1)/2)),fill(1.0e-3,p),fill(-5,p)));#1.0e-5
    ub=vec(vcat(fill(5.1,trsi),fill(1.,p),fill(1.,Int((q*(q-1))/2)),fill(7.,trs),fill(Inf,Int(p*K-K*(K-1)/2)),fill(Inf,p),fill(10.,p)));    
    par = zeros(Float64,size(lb))
    if it == 1
        par = optimize(Optim.only_fg!((F,G,x)->fgthefa!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06,time_limit = 300.0)) #20 min
    elseif it > 1
        par = optimize(Optim.only_fg!((F,G,x)->fgthefa!(F,G,x,X,T,p,q,K,c,n,pairsindex,L_str,indices)),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(f_reltol=1.0000e-06,outer_f_reltol=1.0000e-06,time_limit = 1200.0)) #20 min
    end 
    # if it == 1
    #     par = optimize((x)->pmlthefa(x,X,T,p,q,K,c,n,pairsindex,L_str,indices),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(time_limit = 300.0);autodiff=:forward)#20 min
    # elseif it > 1
    #     par = optimize((x)->pmlthefa(x,X,T,p,q,K,c,n,pairsindex,L_str,indices),lb,ub,theta_init,Fminbox(LBFGS()),Optim.Options(time_limit = 1200.0);autodiff=:forward)#20 min
    # end 
    
    
    # mixture parameter
    betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p)+p+ncov)]   # betas = theta[Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+1):Int(trsi+p+(q*(q-1))/2+trs+(p*K-(K*(K-1)/2))+(p-1)+p+ncov)]
    Ez1 = T[:,1,:]
    Ez0 = T[:,2,:]
    
    pairs = size(pairsindex,2)
    betasnew = NewtonRaphson(betas_gradient, betas, ics, Ez1, Ez0, ncov,pairs)
    
    
    # pg=mean(sum(T,dims=1)[1,:,:]./n,dims=2)[:,1];
   
    # println(par.minimizer,pg[1])
    if !(Optim.converged(par)) && isnan(betas[1]) 
        println(betas[1])
        if start == true
            return zeros(Float64,size(theta,1))
        elseif start == false
            return Optim.converged(par)
        end
    end
    return vcat(par.minimizer,betasnew)#,pg[1])
end


function startvalsthefa(replix,X,n,p,c,q,K,ncov,ics,pairsindex,L_str,indices)
    liksmix = Vector{Float64}(undef,replix); params = Matrix{Float64}(undef,Int((c-1)*p+p+(q*(q-1))/2+(c-1-2)+(p*K-(K*(K-1)/2))+(p)+p)+ncov,replix);#p-1
    for i in 1:replix
        println(i)

        tau1= vec(hcat(rand(Uniform(-2.3,-1.),p).+reshape(repeat(collect(range(0,3.5,c-1)),inner=p),p,(c-1))))
        lambda1 = rand(Uniform(0.5,.90),p);
        phi=rand(LKJ(q,1))
        phi=phi[tril!(trues(size(phi)), -1)]
        tau2=vec(reshape(rand(Uniform(1.2,4),1*(c-1-2)),1,(c-1-2)).+reshape(repeat(collect(0:(c-2-2)),inner=1),1,(c-1-2)))
        tau2=sort(tau2)
        lamb2= vec(reconstruct_L2_matrix(p, K, indices, rand(Uniform(0.15,.90),Int((p*K-(K*(K-1)/2))))))
        lambda2=lamb2[lamb2.!=0.]
        psi=rand(Uniform(0.5,2.5),p)
        interc=rand(Uniform(1.1,4),p)



        prop = rand(Uniform(0.65,0.95))# prop = rand(Uniform(0.65,0.95))
        bets = vcat(log(prop/(1 - prop)),rand(Uniform(0.03,2),ncov-1))#,rand(Uniform(-1.5,1.5),ncov-1))


        param_initials = vcat(tau1,lambda1,phi,tau2,lambda2,psi,interc,bets)

        # first iteration
        plimix = E_stepMixthefa(param_initials,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,true);
        T = repeat(1 ./sum(plimix,dims=2),1,2).*plimix;
        t = time()
        params[:,i] = M_parsMixthefa(T,param_initials,X,n,p,q,K,c,ics,ncov,pairsindex,L_str,indices,1,t,plimix,false); #solo primo argomento di output
        dt = time() - t
        prrot=params[:,i]
        liksmix[i] = E_stepMixthefa(prrot,q,K,p,n,c,X,ics,ncov,pairsindex,L_str,indices,false,T);
    end
    param_initials0 = params[:,findall(liksmix.==maximum(liksmix))]
    lik=liksmix[findall(liksmix.==maximum(liksmix))]

    tau1=param_initials0[(1):((c-1)*p)]
    tau1=vec(sort(reshape(tau1,p,c-1),dims=2))
    lambda1=param_initials0[((c-1)*p+1):((c-1)*p+p)]
    phi=param_initials0[((c-1)*p+p+1):((c-1)*p+p+Int(q*(q-1)/2))]
    lambda1[lambda1.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda1.<=0.1))
    lambda1[lambda1.>=0.99].=rand(Uniform(0.3,0.9),sum(lambda1.>=0.99))
    phi[phi.<=0.1].=rand(Uniform(0.02,0.9),sum(phi.<=0.1))
    tau2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2))]
    tau2=vec(sort(reshape(tau2,1,c-1-2),dims=2))
    tau2[tau2.<=1.01].=rand(Uniform(1.05,3),sum(tau2.<=1.01))
    lambda2=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2))))]
    psi=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p))]# psi=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p-1))]
    interc=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p)+p)]#interc=param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p-1)+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p)]
    lambda2[lambda2.<=0.1].=rand(Uniform(0.1,0.7),sum(lambda2.<=0.1))
    lambda2[lambda2.>=0.96].=rand(Uniform(0.3,0.9),sum(lambda2.>=0.96))
    psi[psi.<=0.1].=rand(Uniform(0.3,1.7),sum(psi.<=0.1))
    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p)+p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p)+p+ncov)]#    bets = param_initials0[((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+1):((c-1)*p+p+Int(q*(q-1)/2)+(c-1-2)+Int((p*K-(K*(K-1)/2)))+(p-1)+p+ncov)]


    param_initialii = vcat(tau1,lambda1,phi,tau2,lambda2,psi,interc,bets)

    println(param_initialii,liksmix,"first iteration: done")
    return param_initialii
end



end