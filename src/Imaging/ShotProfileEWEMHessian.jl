function ShotProfileEWEMHessian(w::Array{ASCIIString,1},d::Array{ASCIIString,1};pspi=true,nref=5,vp="vp.seis",vs="vs.seis",angx="angx.seis",angy="angy.seis",wav="wav.seis",sz=0.,gz=0.,nangx=1,oangx=0,dangx=1,nangy=1,oangy=0,dangy=1,fmin=0,fmax=80,padt=1,padx=1,verbose=false,sx=[0],sy=[0],kz_eps=0.0005)

	# approximate Hessian for shot profile elastic wave equation migration operator
	rand_string = string(round(Int,rand()*100000))
	mref = [join(["tmp_Hess_m1_",rand_string]);join(["tmp_Hess_m2_",rand_string]);join(["tmp_Hess_m3_",rand_string])]
	mref_env = [join(["tmp_Hess_m1_env_",rand_string]);join(["tmp_Hess_m2_env_",rand_string]);join(["tmp_Hess_m3_env_",rand_string])]
	mref_env_s = [join(["tmp_Hess_m1_env_s_",rand_string]);join(["tmp_Hess_m2_env_s_",rand_string]);join(["tmp_Hess_m3_env_s_",rand_string])]
	mref_fwd = [join(["tmp_Hess_m1_fwd_",rand_string]);join(["tmp_Hess_m2_fwd_",rand_string]);join(["tmp_Hess_m3_fwd_",rand_string])]
	mref_fwd_wd = [join(["tmp_Hess_m1_fwd_wd_",rand_string]);join(["tmp_Hess_m2_fwd_wd_",rand_string]);join(["tmp_Hess_m3_fwd_wd_",rand_string])]
	mref_fwd_adj = [join(["tmp_Hess_m1_fwd_adj_",rand_string]);join(["tmp_Hess_m2_fwd_adj_",rand_string]);join(["tmp_Hess_m3_fwd_adj_",rand_string])]
	mref_fwd_adj_env = [join(["tmp_Hess_m1_fwd_adj_env_",rand_string]);join(["tmp_Hess_m2_fwd_adj_env_",rand_string]);join(["tmp_Hess_m3_fwd_adj_env_",rand_string])]
	mref_fwd_adj_env_s = [join(["tmp_Hess_m1_fwd_adj_env_s_",rand_string]);join(["tmp_Hess_m2_fwd_adj_env_s_",rand_string]);join(["tmp_Hess_m3_fwd_adj_env_s_",rand_string])]
	@compat param = Dict(:pspi=>pspi,:nref=>nref,:vp=>vp,:vs=>vs,:angx=>angx,:angy=>angy,:wav=>wav,:sz=>sz,:gz=>gz,:nangx=>nangx,:oangx=>oangx,:dangx=>dangx,:nangy=>nangy,:oangy=>oangy,:dangy=>dangy,:fmin=>fmin,:fmax=>fmax,:padt=>padt,:padx=>padx,:verbose=>verbose,:sx=>sx,:sy=>sy,:kz_eps=>kz_eps)	
	# adjoint
	ShotProfileEWEM(mref,d,true;param...)
	# forward applied to adjoint
	ShotProfileEWEM(mref,mref_fwd,false;param...)
	# zero any traces that were zero on input
	rand_string = string(round(Int,rand()*100000))
	wd = join(["tmp_Hess_wd_",rand_string])
	CalculateSampling(d[1],wd,cutoff=1e-15)
	WeightingOp(mref_fwd,mref_fwd_wd,false;w=wd)
	# adjoint applied to forward of adjoint
	ShotProfileEWEM(mref_fwd_adj,mref_fwd_wd,true;param...)
	# calculate envelopes and smooth them
	SeisEnvelope(mref[1],mref_env[1])
	SeisEnvelope(mref[2],mref_env[2])
	SeisEnvelope(mref[3],mref_env[3])
	SeisSmooth2(mref_env[1],mref_env_s[1];key=["iaz";"iang"],L1=11,L2=11)
	SeisSmooth2(mref_env[2],mref_env_s[2];key=["iaz";"iang"],L1=11,L2=11)
	SeisSmooth2(mref_env[3],mref_env_s[3];key=["iaz";"iang"],L1=11,L2=11)
	SeisEnvelope(mref_fwd_adj[1],mref_fwd_adj_env[1])
	SeisEnvelope(mref_fwd_adj[2],mref_fwd_adj_env[2])
	SeisEnvelope(mref_fwd_adj[3],mref_fwd_adj_env[3])
	SeisSmooth2(mref_fwd_adj_env[1],mref_fwd_adj_env_s[1];key=["iaz";"iang"],L1=11,L2=11)
	SeisSmooth2(mref_fwd_adj_env[2],mref_fwd_adj_env_s[2];key=["iaz";"iang"],L1=11,L2=11)
	SeisSmooth2(mref_fwd_adj_env[3],mref_fwd_adj_env_s[3];key=["iaz";"iang"],L1=11,L2=11)
	# normalize
	m1,h,e = SeisRead(mref_env_s[1])
	m2,h,e = SeisRead(mref_fwd_adj_env_s[1])
	if maximum(abs(m1)) > 0. 
		m1 = m1/maximum(abs(m1))
	end
	if maximum(abs(m2)) > 0. 
		m2 = m2/maximum(abs(m2))
	end
	# calculate weight
	w1 = m1./(m2 + 0.01)
	SeisWrite(w[1],w1,h,e);

	# normalize
	m1,h,e = SeisRead(mref_env_s[2])
	m2,h,e = SeisRead(mref_fwd_adj_env_s[2])
	if maximum(abs(m1)) > 0. 
		m1 = m1/maximum(abs(m1))
	end
	if maximum(abs(m2)) > 0. 
		m2 = m2/maximum(abs(m2))
	end
	# calculate weight
	w1 = m1./(m2 + 0.01)
	SeisWrite(w[2],w1,h,e);

	# normalize
	m1,h,e = SeisRead(mref_env_s[3])
	m2,h,e = SeisRead(mref_fwd_adj_env_s[3])
	if maximum(abs(m1)) > 0. 
		m1 = m1/maximum(abs(m1))
	end
	if maximum(abs(m2)) > 0. 
		m2 = m2/maximum(abs(m2))
	end
	# calculate weight
	w1 = m1./(m2 + 0.01)
	SeisWrite(w[3],w1,h,e);

	#SeisRemove(mref)
	#SeisRemove(mref_fwd)
	#SeisRemove(mref_fwd_wd)
	#SeisRemove(wd)
	#SeisRemove(mref_fwd_adj)
	#SeisRemove(mref_env)
	#SeisRemove(mref_env_s)
	#SeisRemove(mref_fwd_adj_env)
	#SeisRemove(mref_fwd_adj_env_s)

end


