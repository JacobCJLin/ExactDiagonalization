#type for ED calculation
struct EDtable
    state::Dict
    index::Dict
    normsq;
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
        println(state[i],"\t",codetobasis(ED.state[i],L,base))
    end
end    
end

#Loading EDtable data if exist; if not, construct it
function EDtabledata(L,symf,statecheck,filename)
    if isfile(filename)
    EDfile=open(filename,"r")
    ED=deserialize(EDfile)
    close(EDfile)  
    else
    println("no file: Construct EDtable.")
    ED=generateEDsym(L,symf,statecheck=statecheck)
    EDfile=open(filename,"w")
    serialize(EDfile,ED)
    close(EDfile)
    end    
    return ED
end

#Loading EDtabledata if exist
function loadEDtabledata(EDfilename)
    if isfile(filename)
    EDfile=open(filename,"r")
    ED=deserialize(EDfile)
    close(EDfile)
    return ED 
    else
    println("no file: ",EDfilename)
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

#Loading EDdata if exist
function loadEDdata(EDfilename)
    if isfile(EDfilename)
    EDfile=open(EDfilename,"r")
    readdata=deserialize(EDfile)
    E=readdata["E"]
    S=readdata["S"]
    readdata=0;
    close(EDfile)
    return E,S 
    else
    println("no file: ",EDfilename)
    end
end

#functions for ED
function codetobasis(a::Int64,L::Int64,base=2) #convert a state code a to basis, e.g. 1=0000001
    basis=zeros(Int64,L);
    temp=a;
    for i=1:L
        basis[i]=rem(temp,base);
        temp=div(temp,base);
    end
    return basis;    
end

num2basis(a::Int64,L::Int64,base=2)=codetobasis(a,L,base)

#basis to number converter
function basistocode(basis::Array{Int64},base=2)
    temp=0;
    for i=1:length(basis)
        temp+=basis[i]*base^(i-1);
    end    
    return temp;  
end

basis2num(basis::Array{Int64},base=2)=basistocode(basis,base)

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

#function for calculate ED evolution
function evolution(E,S,Δt)
    eiEt=exp.(-1.0im*Δt*E)
    U=S*Diagonal(eiEt)*S'
    #dim1,dim2=size(S)
    #U=zeros(complex(Float64),dim1,dim1)
    #for i=1:dim1, j=1:dim1
    #   for k=1:dim2
    #    U[i,j]+= S[i,k]*eiEt[k]*conj(S[j,k])
    #    end
    #end
    return U
end

function generateEDsym(L::Int64,symf,base=2;tol=1E-14,statecheck= i -> true)
   #L: total number of sites
   #symf: function for the symmetry operation, requiring input ("code", L, base), output ("codelist","χlist")  
    fulldim=base^L
    counter=0;
    state=Dict();index=Dict();normsq=Dict(); 
    for i=0:fulldim-1
        if statecheck(i) 
            codelist,χlist=symf(i,L,base)
            minimum(codelist)!= i ? continue : nothing
            indlist=findall(x -> x==i , codelist)
            snormsq=sum(χlist[indlist])
            if abs(snormsq) > tol
            counter+=1;
            state[counter]=i;
            index[i]=counter;
            normsq[counter]=snormsq    
            end
        end 
    end
    dim=counter
    ED=EDtable(state,index,normsq,dim,base,L)
    return ED
end

function findrepsym(a,L,symf,base=2)
    codelist,χlist=symf(a,L,base)
    ind=findmin(codelist)
    return ind[1], χlist[ind[2]]
end

#generate states without any symmetry consideration
function generateED(L::Int64,base=2;statecheck= i -> true)
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













