#commonly used symmtry operations in 1d
#---------------------
#translation operation
function translation_code(a::Int64,L::Int64,base=2)
    basis=codetobasis(a,L,base)
    Tbasis=zeros(Int64,L);
    Tbasis[1]=basis[L];
    Tbasis[2:L]=basis[1:L-1];
    return basistocode(Tbasis,base)
end

function translation(vec,ED)
    outvec=zeros(eltype(vec),ED.dim)
    for i=1:ED.dim
       a=ED.state[i]
       Ta=translation_code(a,ED.L,ED.base)
       outvec[ED.index[Ta]]+=vec[i]  
    end
    return outvec
end

#space inversion operation
function inversion_code(a::Int64,L::Int64,base=2)
   oldarray=codetobasis(a,L,base)
   newarray=zeros(Int64,L)
    for i=1:L
    newarray[i]=oldarray[L-i+1]
    end
    return basistocode(newarray,base)
end

function inversion(vec,ED)
    outvec=zeros(eltype(vec),ED.dim)
    for i=1:ED.dim
       a=ED.state[i]
       Ia=inversion_code(a,ED.L,ED.base)
       outvec[ED.index[Ia]]+=vec[i]  
    end
    return outvec 
end
#--------------------------------------------











