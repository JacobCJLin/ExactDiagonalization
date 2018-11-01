#functions for entanglement calculation
function Svon(LA,ψ,ED,base=2) #calculate von Neumann entropy with subsys LA of a pure state ψ
    L=ED.L
    LB=L-LA
    M=zeros(base^LA,base^LB)
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

function Sparsity(LA,ψ,ED,base=2,thres=1E-12) #calculate the sparsity a pure state ψ
    L=ED.L
    LB=L-LA
    M=zeros(base^LA,base^LB)
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

function EntSpectrum(LA,ψ,ED,base=2) #calculate the sparsity a pure state ψ
    L=ED.L
    LB=L-LA
    M=zeros(base^LA,base^LB)
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