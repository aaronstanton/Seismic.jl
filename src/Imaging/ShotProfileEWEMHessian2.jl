function ShotProfileEWEMHessian2(w::Array{ASCIIString,1},m::Array{ASCIIString,1};pspi=true,nref=5,vp="vp.seis",vs="vs.seis",angx="angx.seis",angy="angy.seis",wav="wav.seis",sz=0.,gz=0.,nangx=1,oangx=0,dangx=1,nangy=1,oangy=0,dangy=1,fmin=0,fmax=80,padt=1,padx=1,verbose=false,sx=[0],sy=[0],kz_eps=0.0005,eps=0.01,Niter=5)

	# approximate Hessian for shot profile elastic wave equation migration operator
	# using random realizations
	# reference: http://library.seg.org/doi/pdf/10.1190/1.2792973
	
	rand_string = string(round(Int,rand()*100000))
	mref = [join(["tmp_Hess_m1_",rand_string]);join(["tmp_Hess_m2_",rand_string]);join(["tmp_Hess_m3_",rand_string])]
	mref_fwd = [join(["tmp_Hess_m1_fwd_",rand_string]);join(["tmp_Hess_m2_fwd_",rand_string]);join(["tmp_Hess_m3_fwd_",rand_string])]
	mref_fwd_adj = [join(["tmp_Hess_m1_fwd_adj_",rand_string]);join(["tmp_Hess_m2_fwd_adj_",rand_string]);join(["tmp_Hess_m3_fwd_adj_",rand_string])]
	wraw = [join(["tmp_Hess_w1_raw_",rand_string]);join(["tmp_Hess_w2_raw_",rand_string]);join(["tmp_Hess_w3_raw_",rand_string])]
	@compat param = Dict(:pspi=>pspi,:nref=>nref,:vp=>vp,:vs=>vs,:angx=>angx,:angy=>angy,:wav=>wav,:sz=>sz,:gz=>gz,:nangx=>nangx,:oangx=>oangx,:dangx=>dangx,:nangy=>nangy,:oangy=>oangy,:dangy=>dangy,:fmin=>fmin,:fmax=>fmax,:padt=>padt,:padx=>padx,:verbose=>verbose,:sx=>sx,:sy=>sy,:kz_eps=>kz_eps)	

	m1,h,ext = SeisRead(m[1])
	m2,h,ext = SeisRead(m[2])
	m3,h,ext = SeisRead(m[3])
	w1 = zeros(m1)
	w2 = zeros(m2)
	w3 = zeros(m3)
	
	for iter = 1 : Niter
		m1 = randn(size(m1))
		SeisWrite(mref[1],m1,h,ext)
		m2 = randn(size(m2))
		SeisWrite(mref[2],0.*m2,h,ext)
		m3 = randn(size(m3))
		SeisWrite(mref[3],m3,h,ext)
		ShotProfileEWEM(mref,mref_fwd,false;param...)
		ShotProfileEWEM(mref_fwd_adj,mref_fwd,true;param...)
		m1a,h,ext = SeisRead(mref[1])
		m1b,h,ext = SeisRead(mref_fwd_adj[1])
		m2a,h,ext = SeisRead(mref[2])
		m2b,h,ext = SeisRead(mref_fwd_adj[2])
		m3a,h,ext = SeisRead(mref[3])
		m3b,h,ext = SeisRead(mref_fwd_adj[3])
		w1 = w1 + m1a.*m1b 
		w2 = w2 + m2a.*m2b 
		w3 = w3 + m3a.*m3b 
	end
	w1 = w1/Niter
	w2 = w2/Niter
	w3 = w3/Niter
	SeisWrite(wraw[1],w1,h,ext)
	SeisWrite(wraw[2],w2,h,ext)
	SeisWrite(wraw[3],w3,h,ext)

	# smooth the weights
	SeisSmooth2(wraw[1],w[1];key=["iaz";"iang"],L1=11,L2=11,Nrepeat=2)
	SeisSmooth2(wraw[2],w[2];key=["iaz";"iang"],L1=11,L2=11,Nrepeat=2)
	SeisSmooth2(wraw[3],w[3];key=["iaz";"iang"],L1=11,L2=11,Nrepeat=2)

	SeisRemove(mref)
	SeisRemove(mref_fwd)
	SeisRemove(mref_fwd_adj)
	SeisRemove(wraw)

end


