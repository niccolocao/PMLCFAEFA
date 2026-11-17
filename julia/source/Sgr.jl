module Sgr


using Distributions, SortingAlgorithms, SpecialFunctions

export partition_replacement

function partition_replacement(Dx,PM,Pparm,Q=nothing,fake_model=nothing,p=nothing)
  
  if Q === nothing Q = round(Int,maximum(Dx)) end
  
  # Pparm["p"] = reshape(Pparm["p"],:,2) end
  # Pparm["gam"] = reshape(Pparm["gam"],:,2) end
  # Pparm["del"] = reshape(Pparm["del"],:,2) end
  

  Fx = Matrix(Dx)
  Dxx = Matrix(Dx)
  for h in round.(Int,deleteat!(sort(unique(PM)),1))
    righe = 0; colonne = 0
    for i in 1:size(PM,1) 
        for j in 1:size(PM,2)
          if PM[i,j]==h
            righe = hcat(righe,i)
            colonne = hcat(colonne,j)
          end
        end
    end
    righe = unique(righe)[2:end]
    colonne = unique(colonne)[2:end]
    K = reshape(Dx[righe,colonne],length(righe),length(colonne))
    
    R = replacement_matrix(Q,Pparm["p"][h,:],Pparm["gam"][h,:],Pparm["del"][h,:])
    Df = rdatarepl(K,R,false)["Fx"]
    
    for i in 1:length(righe) 
      for j in 1:length(colonne) 
        Fx[righe[i],colonne[j]] = Df[i,j]
      end
    end
    
  end
  
  Delta = Dxx-Fx
  Delta[findall(Delta.!=0)] .= 1
  Fperc = sum(Delta)/(prod(size(Delta)))*100
  println(join([round(Fperc,digits=2),"% of data replaced."]))
  return Fx 
end


function rdatarepl(Dx,RM,printfp=true) 
  Fx = Matrix(Dx)
  for j in 1:size(Dx,2)
    for i in 1:size(RM,1) 
      quali = findall(Dx[:,j].==i) 
      Fx[quali,j] .= Distributions.wsample(collect(1:size(RM,2)),Vector{Float64}(RM[i,:]),length(quali),replace=true)
    end
  end
  
  Delta = Dx-Fx
  Delta[findall(Delta.!=0)] .= 1
  Fperc = sum(Delta)/(prod(size(Delta)))*100
  
  if (printfp) println(join([round(Fperc,2),"% of data replaced."])) end
  return Dict("Fx"=>Fx,"Fperc"=>Fperc)
end


function replacement_matrix(Q,p=[0,0],gam=[1,1],del=[1,1])  #togli o gem e del o fake_model ,fake_model=["uninformative","average","slight","extreme"]  
  if (length(p)<2) 
    p = repeat(p,2)
  end
  if (length(gam)<2) 
    gam = repeat(gam,2)
  end
  if (length(del)<2) 
    del = repeat(del,2)
  end
  
  # if (fake_model!="uninformative") 
  #   MF = model_fake_par(fake_model)
  #   gam = MF["gam"]; del = MF["del"]
  # end
  
  fake = "both"
  if ((p[1]==0)&&(p[2]>0)) fake = "bad" end
  if ((p[2]==0)&&(p[1]>0)) fake = "good" end
    
  R = Matrix(undef,Q,Q)
  if (fake=="good")  
    p = p[1]
    gam = gam[1]
    del = del[1]
    for i in 1:Q
      for j in 1:Q
        # println([j i])
        R[i,j] = pfakegood(j,i,p,Q,gam,del)
      end
    end
  end
  if (fake=="bad")   
    p = p[2]
    gam = gam[2]
    del = del[2] 
    for i in 1:Q
      for j in 1:Q
        # println([j i])
        R[i,j] = pfakebad(j,i,p,Q,gam,del)
      end
    end
  end
  if (fake=="both") 
    for i in 1:Q
      for j in 1:Q
        # println([j i])
        R[i,j] = pfake(j,i,p,Q,gam,del)
      end
    end
  end
  
  
  #return(list(gam,del))
  return R
end


# function replacement_matrix(Q,p=[0,0],gam=[1,1],del=[1,1])  #togli o gem e del o fake_model ,fake_model=["uninformative","average","slight","extreme"]  
#   if (length(p)<2) 
#     p = repeat(p,2)
#   end
#   if (length(gam)<2) 
#     gam = repeat(gam,2)
#   end
#   if (length(del)<2) 
#     del = repeat(del,2)
#   end
  
#   # if (fake_model!="uninformative") 
#   #   MF = model_fake_par(fake_model)
#   #   gam = MF["gam"]; del = MF["del"]
#   # end
  
#   fake = "both"
#   if ((p[1]==0)&&(p[2]>0)) fake = "bad" end
#   if ((p[2]==0)&&(p[1]>0)) fake = "good" end
    
#   if (fake=="good")  
#     FUN = "pfakegood"
#     p = p[1]
#     gam = gam[1]
#     del = del[1]
#   end
#   if (fake=="bad")  
#     FUN = "pfakebad"  
#     p = p[2]
#     gam = gam[2]
#     del = del[2]
#   end
#   if (fake=="both") 
#     FUN = "pfake"
#   end
  
#   R = Matrix(undef,Q,Q)
#   for i in 1:Q
#     for j in 1:Q
#       # println([j i])
#       R[i,j] = getfield(Main,Symbol(FUN))(j,i,p,Q,gam,del)
#     end
#   end
  
#   #return(list(gam,del))
#   return R
# end


 
function pfake(k,h=k,p=[0,0],Q=5,gam=[1,1],del=[1,1]) 
  if (fake_model!="uninformative") 
    MF = model_fake_par(fake_model)
    gam = MF["gam"]; del = MF["del"]
  end
  
  if ((k>Q)||(h>Q)) 
    println("value/s out of scale boundaries")
    outP = 0
  end

  if ((1<=k)&&(k<h)&&(h<=Q)) 
    outP = p[2]*dgBetaD(k,1,h-1,gam[2],del[2]) 
  end

  if ((1<=k)&&(k==h)&&(h<=Q)) 
    outP = 1-(p[1]+p[2])
    if (k==1) 
      outP = 1-p[1] 
    end
    if (k==Q) 
      outP = 1-p[2] 
    end
  end

  if ((1<=h)&&(h<k)&&(k<=Q)) 
    outP = p[1]*dgBetaD(k,h+1,Q,gam[1],del[1]) 
  end

  return outP 
end 




function model_fake_par(fake_model) 
  gam = Vector{Real}(undef,2); del = Vector{Real}(undef,2); 
  if (fake_model=="uninformative") 
    gam = [1,1]; del = [1,1]
  end
  if (fake_model=="average") 
    gam = [3,3]; del = [3,3]
  end
  if (fake_model=="slight") 
    gam = [1.5,4]; del = [4,1.5]
  end
  if (fake_model=="extreme") 
    gam = [4,1.5]; del = [1.5,4]
  end
  
  return Dict("gam"=>gam,"del"=>del)
end

function pfakegood(k,h=k,p=0,Q=5,gam=1,del=1) 
  if ((k>Q) || (h>Q)) 
     outP = 0
  end
  gB = dgBetaD(k,h+1,Q,gam,del)
  if ((h==k)&&(h==Q)) 
    outP = 1
  else 
    if ((1<=h)&&(h<k)&&(k<=Q)) outP = p*gB end
    if ((1<=k)&&(k==h)&&(h<Q)) outP = 1-p end
    if ((1<=k)&&(k<h)&&(h<=Q)) outP = 0 end
  end
  return outP
  #return(list(outP=outP,gBetaD=gB,Gd=gB1))
end


function isinsorted(a, x)
  pos = searchsortedfirst(a, x)
  (pos <= length(a)) && (a[pos] == x)
end

function isin(x, y)
  sorted_y = sort(y, alg=RadixSort)
  isinsorted.(Ref(sorted_y), x)
end

function dgBetaD(x,a=minimum(x),b=maximum(x),gam=1,del=1,ct=1) 
  X = collect(a:b)
  t1 = (beta(gam,del))^(-1)
  t2 = ((b-a+2*ct)^(gam+del-1))^(-1)
  t3 = (X.-a.+ct).^(gam-1)
  t4 = (b.-X.+ct).^(del-1)
  l = t1.*t2.*t3.*t4
  L = sum(l)
  
  if (x-a+ct) < 0 && modf(gam)[1] > 0 #for negative numbers use let ... (no parentesi)
    t3 = 0
  else
    t3 = (x-a+ct)^(gam-1)
  end
  t4 = (b-x+ct)^(del-1)
  l = t1*t2*t3*t4
  L = l/L
  
  if isin(x,X)
    y = L
  else
    y = 0
  end
  return y
end


function pfakebad(k,h=k,p=0,Q=5,gam=1,del=1) 
  if ((k>Q)||(h>Q)) 
    outP = 0
  end
  gB = dgBetaD(k,1,h-1,gam,del)
  if ((h==k)&&(h==1)) 
    outP = 1
  else 
    if ((1<=k)&&(k<h)&&(h<=Q)) outP = p*gB end
    if ((1<k)&&(k==h)&&(h<=Q)) outP = 1-p end
    if ((1<=h)&&(h<k)&&(k<=Q)) outP = 0 end
  end
  return outP
  #return(list(outP=outP,gBetaD=gB,Gd=gB1))
end


function dgBeta(x,a=minimum(x),b=maximum(x),gam=1,del=1) 
  e1 = 1/(beta(gam,del)*(b-a)^(gam+del-1))
  d = e1*(x-a)^(gam-1)*(b-x)^(del-1)
  return d
end



end