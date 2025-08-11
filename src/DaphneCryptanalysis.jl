module DaphneCryptanalysis
using DaphneCipher,Base.Threads,OffsetArrays,CairoMakie,Printf
import DaphneCipher:stepp
import DaphneCipher:invStep
import DaphneCipher:mul257
import DaphneCipher:mulOdd
import DaphneCipher:left
import DaphneCipher:right
import OffsetArrays:Origin
export stepRow,interstep,nonlinearity,sameness,concoctShiftRegister,decryptOne
export avalanche
export plotNonlinearity,plotSameness,chosenCiphertext16M

function hadamard(buf::OffsetVector{<:Real})
  tmp0=copy(buf)
  tmp1=copy(buf)
  sz=length(buf)
  h=sz÷2
  if ispow2(sz) && Origin(buf)==Origin((0,))
    while h>0
      for i in 0:sz-1
	j=i⊻h
	if i>j
	  @inbounds tmp1[i]=tmp0[j]-tmp0[i]
	else
	  @inbounds tmp1[i]=tmp0[j]+tmp0[i]
	end
      end
      tmp0,tmp1=tmp1,tmp0
      if h%3==1
	tmp0./=2
      end
      h÷=2
    end
  end
  if length(buf)%3==2
    tmp0.*=√2
  end
  tmp0
end

function avalanche(buf::OffsetVector{<:Integer},nbits::Integer=0)
  if nbits<=0
    nbits=count_zeros(zero(eltype(buf))) # this will error if BigInt
  end
  sz=length(buf)
  szbits=count_ones(sz-1)
  @assert ispow2(sz)
  @assert sz>2 # 2 is too small to satisfy the strict avalanche criterion
  ret=OffsetMatrix(fill(-sz÷4,nbits,szbits),-1,-1)
  for i in 0:szbits-1
    h=1<<i
    for n in eachindex(buf)
      if n⊻h>n
        x=buf[n⊻h]⊻buf[n]
        for j in 0:nbits-1
          ret[j,i]+=(x>>j)&1
        end
      end
    end
  end
  ret
end

# Tests of the S-box, preceded and followed by multiplications,
# and of the multiplications by themselves

function stepRow(left::UInt8,right::UInt8)
  row=OffsetVector(UInt8[],-1)
  for i in 0x00:0xff
    push!(row,stepp(i,left,right))
  end
  row
end

function interstep(left::UInt8,right::UInt8)
  row=OffsetVector(UInt8[],-1)
  for i in 0x00:0xff
    push!(row,mul257[mulOdd[i,right],left])
  end
  row
end

function nonlinearity(bytes::OffsetVector{UInt8})
  buf=OffsetVector(Float64[],-1)
  for b in bytes
    for i in 0:7
      push!(buf,((b>>i)&1)*-2+1)
    end
  end
  had=hadamard(buf)
  maxNonlin=√(length(had))
  (maxNonlin-maximum(abs.(had)))/maxNonlin
end

function sameness(a::OffsetVector{T},b::OffsetVector{T}) where T
  count(map(==,a,b))
end

"""
    plotNonlinearity()

Plot a histogram of the nonlinearity of step rows. A step row is the result of
fixing l and r in stepp(x,l,r) and running x through all 256 bytes.

Most of the nonlinearities are between 7/8 (which is the nonlinearity of the S-box)
and 15/16. The maximum conceivable nonlinearity is 0.9558, which is unattainable
because the number of bits (2048) is not a power of 4.
"""
function plotNonlinearity()
  data=OffsetVector(fill(0.,65536),-1)
  @threads for i in 0:65535
    data[i]=nonlinearity(stepRow(UInt8(i÷256),UInt8(i%256)))
  end # Data are multiples of 1/512
  minbar=round(minimum(data);digits=9,base=2)
  maxbar=round(maximum(data);digits=9,base=2)
  nl=Figure(size=(1189,841))
  nlax=Axis(nl[1,1],
    title="Daphne Stepping Function Nonlinearity",
    xlabel="Nonlinearity",
    ylabel="Number of stepping function rows")
  hist!(nlax,OffsetArrays.no_offset_view(data),bins=minbar-1/1024:1/512:maxbar+1/1024,
        color=:values,bar_labels=:values,
        label_formatter=x->x<1000 ? string(round(Int,x)) : string(round(Int,x/1000))*'k')
  save("daphne-nonlinearity.svg",nl)
end

"""
    plotSameness()

Plot histograms of the number of same bytes in different rows of both stepp(x,l,r)
= mulOdd[sbox[mul257[x,l]],r] and interstep = mul257[mulOdd[i,r],l], which is what
happens between successive uses of the S-box. The plot of stepp looks like a Poisson
distribution; the maximum is 7. The plot of interstep tails off to 34, which is only
2/15 of 255, which shows that the two multiplications do not resemble each other.
While mul257 is a cyclic group, mulOdd is not, and their identities differ by 1.

If thinRatio is 1, this uses 8 GB just for the arrays, and some more gigs for plotting.
"""
function plotSameness()
  thinRatio=257 # 1 for all data (slow), 257 to thin the data
  rows=OffsetVector(fill(OffsetVector(UInt8[],-1),65536),-1)
  tweens=copy(rows)
  @threads for i in 0:65535
    rows[i]=stepRow(UInt8(i÷256),UInt8(i%256))
    tweens[i]=interstep(UInt8(i÷256),UInt8(i%256))
  end
  stp=UInt16[]
  intr=UInt16[]
  k=0
  for i in 1:65535
    for j in 0:i-1
      if k%thinRatio==0
        push!(stp,i)
        push!(intr,j)
      end
      k+=1
    end
  end
  @threads for n in eachindex(stp)
    i=stp[n]
    j=intr[n]
    stp[n]=sameness(rows[i],rows[j])
    intr[n]=sameness(tweens[i],tweens[j])
    if intr[n]>=35 # max is 34, attained about 68619 times
      @printf "%04x %04x %d\n" i j intr[n]
    end
  end
  isam=Figure(size=(1189,841))
  isamax=Axis(isam[1,1],
    title="Daphne Interstep Sameness",
    xlabel="Number of same bytes",
    ylabel="Number of pairs of rows")
  hist!(isamax,intr,bins=-0.5:maximum(intr)+.5)
  save("daphne-interstep-sameness.svg",isam)
  sam=Figure(size=(1189,841))
  samax=Axis(sam[1,1],
    title="Daphne Stepping Function Sameness")
  hist!(samax,stp,bins=-0.5:maximum(stp)+.5)
  save("daphne-step-sameness.svg",sam)
  maximum(intr)
end

# Chosen-ciphertext cryptanalysis, where the chosen ciphertext consists of
# 0x00 and 0x01, and enough data are accumulated to get all accumulator values

function concoctShiftRegister(bits::Integer)
  ret=UInt8[]
  for i in 1:16
    push!(ret,bits&1)
    bits÷=2
  end
  ret
end

function decryptOne(key::Vector{UInt8},accBits::Integer)
  daph=Daphne()
  daph.key=key
  daph.sreg=concoctShiftRegister(accBits)
  daph.acc=accBits>>16
  l=left(daph)
  r=right(daph)
  plain0=invStep(0,l,r)
  plain1=invStep(1,l,r)
  plain1*256+plain0
end

function chosenCiphertext16M(key::Vector{UInt8})
  ret=OffsetVector(fill(0x0000,16777216),-1)
  @threads for i in 0:16777215
    ret[i]=decryptOne(key,i)
  end
  ret
end

struct PlainAcc
  acc		::UInt8
  pt		::UInt8
  ct		::UInt8
  shiftreg	::UInt
end

keepRunning::Bool=false

plainAccs=Channel{PlainAcc}(256)

function chosen16MWorker(daph::Daphne)
  shiftreg=0
  i=0
  acc=0x0
  while keepRunning
    ct=UInt8(rand(Bool))
    pt=decrypt!(daph,ct)
    if (i>16)
      put!(plainAccs,PlainAcc(acc,pt,ct,shiftreg))
    end
    i+=1
    acc+=pt
    shiftreg=(shiftreg+ct*65536)>>1
  end
end

function chosenCiphertext16M(daph::Daphne)
  ret=OffsetVector(fill(0x0000,16777216),-1) # decryptOne never returns 0x0000
  full=0
  shiftreg=0
  i=0
  acc=0x0
  while full<10000000
    ct=UInt8(rand(Bool))
    pt=decrypt!(daph,ct)
    if (i>16)
      inx=acc*65536+shiftreg
      if ret[inx]==0
        ret[inx]=0x101*pt
      elseif ret[inx]%0x101==0
        if ct>0
          newval=ret[inx]&0xff+pt*256
        else
          newval=ret[inx]&0xff00+pt
        end
        @printf "%04x %04x\n" newval ret[inx]
        if newval!=ret[inx]
          full+=1
        end
        ret[inx]=newval
      else
        @assert ret[inx]&0xff==pt || ret[inx]>>8==pt
      end
    end
    i+=1
    acc+=pt
    shiftreg=(shiftreg+ct*65536)>>1
  end
  ret
end

end # module DaphneCryptanalysis
