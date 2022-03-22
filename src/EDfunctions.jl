#type for ED calculation
struct EDtable
    code::Dict   #given an index, output the state code
    index::Dict  #given a state code, output its index
    normsq;      #The norm square of the state
    dim::Int64   #dimension of the total Hilbert space 
    N::Int64     #number of sites
    d::Int64     #dimension of the local Hilbert space

end

#convert a state code a to computational basis, e.g. 1=0000001
function codetobasis(a::Int64,N::Int64,d=2) 
    basis=zeros(Int64,N);
    temp=a;
    for i=1:N
        basis[i]=rem(temp,d);
        temp=div(temp,d);
    end
    return basis;    
end

#basis to number converter
function basistocode(basis::Array{Int64},d=2)
    temp=0;
    for i=1:length(basis)
        temp+=basis[i]*d^(i-1);
    end    
    return temp;  
end

#generate states without any symmetry consideration
function generateED(N::Int64,d=2;statecheck= i -> true)
    fulldim=d^N
    counter=0;
    code=Dict()
    index=Dict()
    for i=0:fulldim-1
        if statecheck(i)
            counter+=1;
            code[counter]=i;
            index[i]=counter;  
        end
    end
    dim=counter
    ED=EDtable(code,index,nothing,dim,N,d)
    return ED
end


function generateHmat(Hcode,ED::EDtable)
    #Hcode: takes a computational basis |a>  and generate H|a>
    dim=ED.dim
    Hmat=zeros(dim,dim)
    for i=1:dim
        outvec=Hcode(ED.code[i],ED)
        Hmat[:,i]=outvec
    end
    return Hmat
end


#some tools
function printψ(ψ,ED::EDtable;tol=1E-16) #printour the state with nonzero amplitudes
    dim=ED.dim
    N=ED.N
    d=ED.d
    for i=1:dim
        if abs(ψ[i])>tol
            println(ψ[i],"\t",codetobasis(ED.code[i],N,d))
        end
    end    
end

#function for calculate evolution operator 
function generateUt(E,S,Δt)
    eiEt=exp.(-1.0im*Δt*E)
    Ut=S*Diagonal(eiEt)*S'
    return Ut
end

####not yet modified below--------------------------------------------------------------

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
function loadEDtabledata(filename)
    if isfile(filename)
    EDfile=open(filename,"r")
    ED=deserialize(EDfile)
    close(EDfile)
    return ED 
    else
    println("no file: ",filename)
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


function generateEDsym(L::Int64,symf,base=2;tol=1E-14,statecheck= i -> true,symfile = nothing)
   #L: total number of sites
   #symf: function for the symmetry operation, requiring input ("code", L, base), output ("codelist","χlist")  
    fulldim=base^L
    counter=0;
    state=Dict();index=Dict();normsq=Dict();symlist=Dict() 
    for i=0:fulldim-1
        if statecheck(i) 
            codelist,χlist=symf(i,L,base)
            if !isnothing(symfile)
            symlist[i,"code"]=codelist    
            symlist[i,"χ"]=χlist
            end    
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
    if !isnothing(symfile)
    return ED,symlist
    else
    return ED
    end
end

function findrepsym(a,L,symf,base=2)
    codelist,χlist=symf(a,L,base)
    ind=findmin(codelist)
    return ind[1], χlist[ind[2]]
end



function savesymfile(symlist,symfilename)
    file=open(symfilename,"w")
    serialize(file,symlist)
    close(file)
end

function loadsymfile(symfilename)
    if isfile(symfilename)
    file=open(symfilename,"r")
    symlist=deserialize(file)
    close(file)
    return symlist
    else
        error("no symmetry file!")
    end
end












