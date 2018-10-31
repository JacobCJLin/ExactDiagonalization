module ExactDiagonalization

using LinearAlgebra

include("ED.jl")
include("Entanglement.jl")

export EDtable,printstate,num2basis,basis2num,tran,Op_translation, inversion,evolution,EDdata
export GenerateED,GenerateEDK,GenerateEDKI,findrep
export Svon,Sparsity,EntSpectrum



end