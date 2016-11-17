function ShotProfileLSEWEM_suboff(m::Array{ASCIIString,1},d::Array{ASCIIString,1};Niter=10,mu=[0.;0.],cost="cost.txt",precon=false,Nsmooth=5,wm="NULL",wd="NULL",pspi=true,nref=5,vp="vp.seis",vs="vs.seis",wav="wav.seis",sz=0.,gz=0.,npx=1,opx=0,dpx=1,nhx=1,ohx=0,dhx=1,fmin=0,fmax=80,padt=1,padx=1,verbose=false,sx=[0],kz_eps=0.5)

	# least squares shot profile wave equation migration of isotropic 2C data. 

	rand_string = string(round(Int,rand()*100000))
	wd_flag = true
	if (wd == "NULL")
	  wd = join(["tmp_LSM_wd_",rand_string])
	  CalculateSampling(d[1],wd,cutoff=1e-15)
	  wd_flag = false
        end  
	@compat param_wd = Dict(:w=>wd)
	@compat param_mig = Dict(:pspi=>pspi,:nref=>nref,:vp=>vp,:vs=>vs,:wav=>wav,:sz=>sz,:gz=>gz,:nhx=>nhx,:ohx=>ohx,:dhx=>dhx,:fmin=>fmin,:fmax=>fmax,:padt=>padt,:padx=>padx,:verbose=>verbose,:sx=>sx,:kz_eps=>kz_eps)
	@compat param_wm = Dict(:w=>wm)
	@compat param_ang = Dict(:npx=>npx,:opx=>opx,:dpx=>dpx,:nhx=>nhx,:ohx=>ohx,:dhx=>dhx)
	if (precon)
		@compat param_smooth = Dict(:Nsmooth=>Nsmooth,:Nrepeat=>1)
		if (wm != "NULL")
			operators = [WeightingOp ;  ShotProfileEWEM_suboff ; OffsetToAngle ; WeightingOp ; SmoothGathers]
			parameters = [param_wd ; param_mig ; param_ang ; param_wm ; param_smooth]
		else
			operators = [WeightingOp ;  ShotProfileEWEM_suboff ; OffsetToAngle ; SmoothGathers]
			parameters = [param_wd ; param_mig ; param_ang ; param_smooth]		
		end	
	else
		if (wm != "NULL")
			operators = [WeightingOp ;  ShotProfileEWEM_suboff ; OffsetToAngle ; WeightingOp]
			parameters = [param_wd ; param_mig ; param_ang ; param_wm]
		else
			operators = [WeightingOp ;  ShotProfileEWEM_suboff ; OffsetToAngle]
			parameters = [param_wd ; param_mig ; param_ang]
		end	
	end
	ConjugateGradients(m,d,operators,parameters,cost,Niter=Niter,mu=mu)

	if (wm != "NULL")
		m_w = [join(["tmp_LSM_m1_w_",rand_string]);join(["tmp_LSM_m2_w_",rand_string])]
		WeightingOp(m,m_w,false;w=wm)
		SeisCopy(m_w,m)
		SeisRemove(m_w)	
	end
        if (wd_flag == false)
	  SeisRemove(wd)
        end

end
