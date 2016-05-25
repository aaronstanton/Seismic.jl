function ShotProfileLSEWEM(m::Array{ASCIIString,1},d::Array{ASCIIString,1};Niter=10,mu=[0.;0.;0.],cost="cost.txt",precon=false,Nsmooth=5,tmute=0.6,vmute=1500.,wm="NULL",pspi=true,nref=5,vp="vp.seis",vs="vs.seis",angx="angx.seis",angy="angy.seis",wav="wav.seis",sz=0.,gz=0.,nangx=1,oangx=0,dangx=1,nangy=1,oangy=0,dangy=1,fmin=0,fmax=80,padt=1,padx=1,verbose=false,sx=[0],sy=[0],kz_eps=0.0005)

	# least squares shot profile wave equation migration of isotropic 3C data. 

	rand_string = string(round(Int,rand()*100000))
	wd = join(["tmp_LSM_wd_",rand_string])
	CalculateSampling(d[1],wd,cutoff=1e-15)
	@compat param_wd = Dict(:w=>wd)
	@compat param_mute = Dict(:tmute=>tmute,:vmute=>vmute)
	@compat param_mig = Dict(:pspi=>pspi,:nref=>nref,:vp=>vp,:vs=>vs,:angx=>angx,:angy=>angy,:wav=>wav,:sz=>sz,:gz=>gz,:nangx=>nangx,:oangx=>oangx,:dangx=>dangx,:nangy=>nangy,:oangy=>oangy,:dangy=>dangy,:fmin=>fmin,:fmax=>fmax,:padt=>padt,:padx=>padx,:verbose=>verbose,:sx=>sx,:sy=>sy,:kz_eps=>kz_eps)
	@compat param_wm = Dict(:w=>wm)
	if (precon)
		@compat param_smooth = Dict(:Nsmooth=>Nsmooth,:Nrepeat=>1)
		if (wm != "NULL")
			operators = [WeightingOp ;  SeisMute ; ShotProfileEWEM ; WeightingOp ; SmoothGathers]
			parameters = [param_wd ; param_mute ; param_mig; param_wm ; param_smooth]
		else
			operators = [WeightingOp ;  SeisMute ; ShotProfileEWEM ; SmoothGathers]
			parameters = [param_wd ; param_mute ; param_mig ; param_smooth]		
		end	
	else
		if (wm != "NULL")
			operators = [WeightingOp ;  SeisMute ; ShotProfileEWEM ; WeightingOp]
			parameters = [param_wd ; param_mute ; param_mig; param_wm]
		else
			operators = [WeightingOp ;  SeisMute ; ShotProfileEWEM]
			parameters = [param_wd ; param_mute ; param_mig]
		end	
	end
	ConjugateGradients(m,d,operators,parameters,cost,Niter=Niter,mu=mu)

	if (wm != "NULL")
		m_w = [join(["tmp_LSM_m1_w_",rand_string]);join(["tmp_LSM_m2_w_",rand_string]);join(["tmp_LSM_m3_w_",rand_string])]
		WeightingOp(m,m_w,false;w=wm)
		SeisCopy(m_w,m)
		SeisRemove(m_w)	
	end

	SeisRemove(wd)

end
