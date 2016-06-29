using Seismic
using Base.Test

# test of Conjugate Gradients
nd = 1000
nm = 300
L = randn(nd,nm)
m = 10000*randn(nm);
d = L*m;
d = d + 0.2*rand(nd);

m1,cost = ConjugateGradients(d,[MatrixMultiplyOp],[Dict(:matrix=>L)],Niter=nm,mu=0.2)
ext = Seismic.Extent(size(d,1),size(d,2),1,1,1,
0,0,0,1,1,
1,1,1,1,1,
"","","","","",
"","","","","",
"")
h = Header[]
for ix = 1:size(d,2)
		push!(h,Seismic.InitSeisHeader())
		h[ix] = Seismic.InitSeisHeader();
		h[ix].tracenum = ix;
		h[ix].n1 = size(d,1);
		h[ix].d1 = 1;
end

println(size(d))

SeisWrite("tmp_d.seis",d,h,ext)

ConjugateGradients("tmp_m.seis","tmp_d.seis",[MatrixMultiplyOp],[Dict(:matrix=>L)],"tmp_cost.txt",Niter=nm,mu=0.2)
m2,ext = SeisRead("tmp_m.seis")

# test that quality factor between disk and memory based CG 
# is greater than 50 Decibels
quality_factor = 10*log10(norm(m1[:],2)/norm(m2[:]-m1[:],2))
println("Quality factor = ",quality_factor)
@test quality_factor > 10.
