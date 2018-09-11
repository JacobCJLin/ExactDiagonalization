#functions for entanglement calculation
function Svon(LA,ψ,ED) #calculate von Neumann entropy with subsys LA of a pure state ψ
    L=ED.L
    LB=L-LA
    M=zeros(2^LA,2^LB)
    for i=1:2^LA, j=1:2^LB
        intA=num2basis(i-1,LA);
        intB=num2basis(j-1,LB);
        append!(intA,intB)
        sta=basis2num(intA)
        
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

function Sparsity(LA,ψ,thres=1E-12,ED) #calculate the sparsity a pure state ψ
    L=ED.L
    LB=L-LA
    M=zeros(2^LA,2^LB)
    for i=1:2^LA, j=1:2^LB
        intA=num2basis(i-1,LA);
        intB=num2basis(j-1,LB);
        append!(intA,intB)
        sta=basis2num(intA)
        
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
    return rk/2^LA

end

function EntSpectrum(LA,ψ,ED) #calculate the sparsity a pure state ψ
    L=ED.L
    LB=L-LA
    M=zeros(2^LA,2^LB)
    for i=1:2^LA, j=1:2^LB
        intA=num2basis(i-1,LA);
        intB=num2basis(j-1,LB);
        append!(intA,intB)
        sta=basis2num(intA)
        
        if haskey(ED.index,sta)
            M[i,j]=ψ[ED.index[sta]]
        end
    end
    
    return svdvals(M)

end