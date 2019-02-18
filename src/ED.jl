#type for ED calculation
struct EDtable
    state::Dict
    index::Dict
    Snorm;
    dim::Int64
    base::Int64
    L
end
#some tools
function printstate(state,ED::EDtable;tol=1E-12)
dim=ED.dim
L=ED.L
base=ED.base
for i=1:dim
    if abs(state[i])>tol
        println(state[i],"\t",num2basis(ED.state[i],L,base))
    end
end    
end

#Loading EDdata if exist; if not, diagonalize it
function EDdata(H,EDfilename)
    if isfile(EDfilename)
    EDfile=open(EDfilename,"r")
    readdata=deserialize(EDfile)
    E=readdata["E"]
    S=readdata["S"]
    readdata=0;
    close(EDfile)
else
    println("no file: Doing ED calculation.")
sol=eigen(H);
E=sol.values;
S=sol.vectors;
EDdata=Dict()
EDdata["E"]=E;
EDdata["S"]=S;
EDfile=open(EDfilename,"w")
serialize(EDfile,EDdata)
close(EDfile)
end
   return E,S 
end


#functions for ED
#functions
function code2basis(a::Int64,L::Int64,base=2) #convert a state code a to basis, e.g. 1=0000001
    basis=zeros(Int64,L);
    temp=a;
    for i=1:L
        basis[i]=rem(temp,base);
        temp=div(temp,base);
    end
    return basis;    
end

num2basis(a::Int64,L::Int64,base=2)=code2basis(a,L,base)

#basis to number converter
function basis2code(basis::Array{Int64},base=2)
    temp=0;
    for i=1:length(basis)
        temp+=basis[i]*base^(i-1);
    end    
    return temp;  
end

basis2num(basis::Array{Int64},base=2)=basis2code(basis,base)

#translation operation
function translation_code(a::Int64,L::Int64,base=2)
    basis=code2basis(a,L,base)
    Tbasis=zeros(Int64,L);
    Tbasis[1]=basis[L];
    Tbasis[2:L]=basis[1:L-1];
    return basis2code(Tbasis,base)
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
   oldarray=code2basis(a,L,base)
   newarray=zeros(Int64,L)
    for i=1:L
    newarray[i]=oldarray[L-i+1]
    end
    return basis2code(newarray,base)
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

#function for calculate ED evolution
function evolution(E,S,Δt)
    eiEt=exp.(-1.0im*Δt*E)
    dim1,dim2=size(S)
    U=zeros(complex(Float64),dim1,dim1)
    for i=1:dim1, j=1:dim1
       for k=1:dim2
        U[i,j]+= S[i,k]*eiEt[k]*conj(S[j,k])
        end
    end
    return U
end



##################################################################################
#### Old
##################################################################################
#generate states without any symmetry consideration
function GenerateED(L::Int64,base=2;statecheck= i -> true)
    fulldim=base^L
    counter=0;
    state=Dict()
    index=Dict()
    for i=0:fulldim-1
        if statecheck(i)
            counter+=1;
            state[counter]=i;
            index[i]=counter;  
        end
    end
    dim=counter
    ED=EDtable(state,index,nothing,dim,base,L)
    return ED
end

#generate states with momentum resolution
function GenerateEDK(L::Int64,K::Int64,base=2;tol=1E-14,statecheck= i -> true)
    fulldim=base^L
    counter=0;
    state=Dict()
    index=Dict()
    Snorm=Dict()
    for i=0:fulldim-1
      if statecheck(i)  
        leaststate=i;
        tempstate=i;
        for r=1:L-1   #loop over the translational operation; if the state is the smallest then it's representative
            tempstate=translation_code(tempstate,L,base);
            tempstate<leaststate ? leaststate=tempstate : nothing
        end
         if leaststate==i
            #determine it's statenorm
            tempnorm=1;
            tempstate=i;
            for r=1:L-1
                tempstate=translation_code(tempstate,L,base);
                tempstate==i ? tempnorm+=exp(2im*π*K*r/L) : nothing
            end
            
            if abs(tempnorm)>tol  #count the state if statenorm!=0
            counter+=1
            state[counter]=i;
            index[i]=counter;
            Snorm[counter]=sqrt(tempnorm)
            end
          end
        end  #end statecheck  
    end
    dim=counter
    ED=EDtable(state,index,Snorm,dim,base,L)
    return ED
end

#generate states with momentum resolution + inversion symmetry resolution
function GenerateEDKI(L::Int64,K::Int64,Inv::Bool,base=2;tol=1E-12,statecheck= i -> true)
    fulldim=base^L
    Inv ? Isgn=1 : Isgn=-1;  #sign of the inversion symmetry
    counter=0;
    state=Dict()
    index=Dict()
    Snorm=Dict()
    for i=0:fulldim-1
        if statecheck(i)
        leaststate=i;
        tempstate=i;
        Itempstate=inversion_code(tempstate,L,base)
        Itempstate<leaststate ? leaststate=Itempstate : nothing
        for r=1:L-1   #loop over the translational operation; if the state is the smallest then it's representative
            tempstate=translation_code(tempstate,L,base);
            tempstate<leaststate ? leaststate=tempstate : nothing
            Itempstate=inversion_code(tempstate,L,base)   #checked the inversion symmetry
            Itempstate<leaststate ? leaststate=Itempstate : nothing
        end
         if leaststate==i
            #determine it's statenorm
            tempnorm=1;
            tempstate=i;
            Itempstate=inversion_code(tempstate,L,base)
            Itempstate==i ? tempnorm+=Isgn : nothing
            for r=1:L-1
                tempstate=translation_code(tempstate,L,base);
                tempstate==i ? tempnorm+=exp(2im*π*K*r/L) : nothing
                Itempstate=inversion_code(tempstate,L,base)
                Itempstate==i ? tempnorm+=Isgn*exp(2im*π*K*r/L) : nothing
            end
            
            if abs(tempnorm)>tol  #count the state if statenorm!=0
            counter+=1
            state[counter]=i;
            index[i]=counter;
            Snorm[counter]=sqrt(tempnorm)
            end
          end
        end
    end
    dim=counter
    ED=EDtable(state,index,Snorm,dim,base,L)
    return ED
end

function findrep(a::Int64,K::Int64,ED::EDtable,base=2)
    L=ED.L
    tempstate=a;
    for r=0:L-1
        if haskey(ED.index,tempstate)
            return tempstate,exp(2im*π*r*K/L) 
        else
            tempstate=translation_code(tempstate,L,base)
        end
    end
    return -1,0
end

function findrep(a::Int64,K::Int64,Inv::Bool,ED::EDtable,base=2)
    L=ED.L
    Inv ? Isgn=1 : Isgn=-1;  #sign of the inversion symmetry
    tempstate=a;
   
    for r=0:L-1
            Itempstate=inversion_code(tempstate,L,base)
        if haskey(ED.index,tempstate)
            return tempstate,exp(2im*π*r*K/L) 
        elseif haskey(ED.index,Itempstate)
            return Itempstate,Isgn*exp(2im*π*r*K/L)    
        else
            tempstate=translation_code(tempstate,L,base)
        end
    end
    return -1,0
end











