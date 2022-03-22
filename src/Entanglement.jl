function reshapelist(subA,ED) #subA is the set of sites in subsystem A e.g. subA=[1,3,5]
    list=Dict();
    d=ED.d
    N=ED.N
    NA=length(subA)
    fullsys=Array{Int64}(1:N)
    subB=setdiff(fullsys,subA)
    NB=N-NA;
    dimA=d^NA;
    dimB=d^NB;
    for i=1:dimA, j=1:dimB
        basisA=codetobasis(i-1,NA,d);
        basisB=codetobasis(j-1,NB,d);
        basis=zeros(Int64,N);
        basis[subA]=basisA
        basis[subB]=basisB
        code=basistocode(basis,d)
        list[code]=(i,j)
    end
    return list,dimA,dimB
end

function reshapeψ(ψ,subA,ED) #subA is the set of sites in subsystem A e.g. subA=[1,3,5]
    list,dimA,dimB=reshapelist(subA,ED);
    Mψ=zeros(dimA,dimB)
    for code in keys(list)
        Mψ[list[code]...]=ψ[ED.index[code]];
    end
    return Mψ
end

function Schmidtdecomposition(ψ,subA,ED)
    Mψ=reshapeψ(ψ,subA,ED) 
    U,S,V=svd(Mψ)
    return U,S,V 
end

function entanglemententropy(ψ,subA,ED)
    _,spec,_ = Schmidtdecomposition(ψ,subA,ED)
    return -sum(spec .^2 .* (log.(spec.^2)))
end




#Not yet modified ------------------------------------------
#functions for entanglement calculation
function Svon_sym(LA,ψ,ED,symmetry) #calculate von Neumann entropy with subsys LA of a pure state ψ
    base=ED.base
    L=ED.L
    LB=L-LA
    M=zeros(eltype(ψ),base^LA,base^LB)
    for i in eachindex(ψ)
        G,gofr,χg=symmetry(ED.state[i])
        G==1 ? Sn=1 : Sn=ED.Snorm[i]
        for n=1:G
           mi=div(gofr[n],base^LB)+1
           mj=rem(gofr[n],base^LB)+1 
           M[mi,mj]+=ψ[i]*χg[n]/sqrt(G)/Sn
        end
    end
    
    U,S,V=svd(M)
    EE=0
    for i=1:length(S)
        abs(S[i])>1.0E-30 ? EE+=-S[i]^2*log(S[i]^2) : nothing
    end
    return EE
end



#functions for entanglement calculation
function Svon(LA,ψ,ED) #calculate von Neumann entropy with subsys LA of a pure state ψ
    base=ED.base
    L=ED.L
    LB=L-LA
    M=zeros(eltype(ψ),base^LA,base^LB)
    for i=1:base^LA, j=1:base^LB
        intA=num2basis(i-1,LA,base);
        intB=num2basis(j-1,LB,base);
        append!(intA,intB)
        sta=basis2num(intA,base)
        
        if haskey(ED.index,sta)
            M[i,j]=ψ[ED.index[sta]]
        end
    end
    U,S,V=svd(M)
    EE=0
    for i=1:length(S)
        abs(S[i])>1.0E-30 ? EE+=-S[i]^2*log(S[i]^2) : nothing
    end
    return EE
end

function Sparsity(LA,ψ,ED,thres=1E-12) #calculate the sparsity a pure state ψ
    base=ED.base
    L=ED.L
    LB=L-LA
    M=zeros(eltype(ψ),base^LA,base^LB)
    for i=1:base^LA, j=1:base^LB
        intA=num2basis(i-1,LA,base);
        intB=num2basis(j-1,LB,base);
        append!(intA,intB)
        sta=basis2num(intA,base)
        
        if haskey(ED.index,sta)
            M[i,j]=ψ[ED.index[sta]]
        end
    end
    
    S=svdvals(M)
    rk=length(S)
    for i=1:length(S)
       if S[i]<thres
            rk=i-1
            break
        end
    end
    return rk

end

function EntSpectrum(LA,ψ,ED) #calculate the sparsity a pure state ψ
    base=ED.base
    L=ED.L
    LB=L-LA
    M=zeros(eltype(ψ),base^LA,base^LB)
    for i=1:base^LA, j=1:base^LB
        intA=num2basis(i-1,LA,base);
        intB=num2basis(j-1,LB,base);
        append!(intA,intB)
        sta=basis2num(intA,base)
        
        if haskey(ED.index,sta)
            M[i,j]=ψ[ED.index[sta]]
        end
    end
    
    return svdvals(M)

end