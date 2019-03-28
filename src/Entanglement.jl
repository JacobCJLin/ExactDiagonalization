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

function reshapelist(LA,ED)
    list=Dict();
    base=ED.base
    L=ED.L
    LB=L-LA 
    for i=1:base^LA, j=1:base^LB
        intA=num2basis(i-1,LA,base);
        intB=num2basis(j-1,LB,base);
        append!(intA,intB)
        sta=basis2num(intA,base)
        list[sta]=(i,j)
    end
    return list
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