function ShotProfileLSWEM_suboff(m::ASCIIString,d::ASCIIString;Niter=10,mu=0.,cost="cost.txt",precon=false,Nsmooth=5,pspi=true,nref=5,vel="vp.seis",wav="wav.seis",sz=0.,gz=0.,npx=1,opx=0,dpx=1,nhx=1,ohx=0,dhx=1,fmin=0,fmax=80,padt=1,padx=1,verbose=false,sx=[0],sy=[0],kz_eps=0.0005)

	# least squares shot profile wave equation migration of isotropic 1C data. 

	rand_string = string(round(Int,rand()*100000))
	wd = join(["tmp_LSM_wd_",rand_string])
	CalculateSampling(d,wd,cutoff=1e-15)
	@compat param_wd = Dict(:w=>wd)
	@compat param_mig = Dict(:pspi=>pspi,:nref=>nref,:vel=>vel,:nhx=>nhx,:ohx=>ohx,:dhx=>dhx,:wav=>wav,:sz=>sz,:gz=>gz,:fmin=>fmin,:fmax=>fmax,:padt=>padt,:padx=>padx,:verbose=>verbose,:sx=>sx,:sy=>sy,:kz_eps=>kz_eps)
	if (precon)
		operators = [WeightingOp ; ShotProfileWEM_suboff ; SmoothGathers]
		@compat param_smooth = Dict(:Nsmooth=>Nsmooth,:Nrepeat=>1)
		parameters = [param_wd ; param_mig ; param_smooth]
	else
		operators = [WeightingOp ; ShotProfileWEM_suboff]
		parameters = [param_wd ; param_mig]
	end
	ConjugateGradients(m,d,operators,parameters,cost,Niter=Niter,mu=mu)
	SeisRemove(wd)

end
