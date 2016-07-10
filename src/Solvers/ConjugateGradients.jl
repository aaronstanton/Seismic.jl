function ConjugateGradients(d,operators,parameters;Niter=10,mu=0,tol=1.0e-15)
	# Conjugate Gradients following Algorithm 2 from Scales, 1987. 
	# The user provides an array of linear operators. Ensure linear operator(s) pass the dot product.

	# approximate largest eigenvalue of the sequence of operators
	z = LinearOperator(d,operators,parameters,adj=true)
	z = ones(size(z))
	Lz = LinearOperator(z,operators,parameters,adj=false)
	LTLz = LinearOperator(Lz,operators,parameters,adj=true) 
	eig_val = InnerProduct(LTLz,LTLz)/InnerProduct(z,z)
	println("eig_val=",eig_val)

	cost = Float64[]
	r = copy(d)
	s = LinearOperator(r,operators,parameters,adj=true)
	p = copy(s)
	m = zeros(s)
	norm_s0 = norm(s)
	gamma = norm_s0^2
	norm_m = 0.0
	m_max = norm_m
	for iter = 1 : Niter
		q = LinearOperator(p,operators,parameters,adj=false)
		delta = norm(q)^2 + mu*norm(p)^2
		if delta <= 0 
			error("delta <= 0")
		end
		alpha = gamma / delta
		m = m + alpha*p
		r = r - alpha*q
		s = LinearOperator(r,operators,parameters,adj=true)
		s = s - mu*m
		norm_s = norm(s)
		gamma1 = gamma
		gamma = norm_s^2
		beta = gamma / gamma1
		p = s + beta*p	
		norm_m = norm(m)
		m_max = max(m_max,norm_m)
		push!(cost,norm_s / norm_s0)
		if (norm_s <= norm_s0 * tol) || (norm_m * tol >= 1)
			println("tolerance reached, ending at iteration ",iter)
			break;
		end
	end

	return m, cost
end

function ConjugateGradients2(d,operators,parameters;Niter=10,mu=0)
	# Conjugate Gradients following Algorithm 2 from Scales, 1987. 
	# The user provides an array of linear operators. Ensure linear operator(s) pass the dot product.

	mu = sqrt(mu)
	cost = Float64[]
	r = copy(d)
	push!(cost,InnerProduct(r,r))
	g = LinearOperator(r,operators,parameters,adj=true)
	m = zeros(g)
	rr = zeros(g)
	s = copy(g)
	gamma_old = InnerProduct(g,g)
	for iter = 1 : Niter
		t = LinearOperator(s,operators,parameters,adj=false)
		tt = mu.*s
		delta = InnerProduct(t,t) + InnerProduct(tt,tt)
		alpha = gamma_old/delta
		m = m + alpha*s
		r = r - alpha*t
		rr = rr - alpha*tt
		push!(cost,InnerProduct(r,r) + InnerProduct(rr,rr))
		g = LinearOperator(r,operators,parameters,adj=true)
		g = g + mu*rr
		gamma = InnerProduct(g,g)
		beta = gamma/gamma_old
		gamma_old = copy(gamma)
		s = beta*s + g
	end

	return m, cost
end



function ConjugateGradients(m::ASCIIString,d::ASCIIString,operators,parameters,cost_file::ASCIIString;Niter=10,mu=0,tol=1.0e-15)
	# Conjugate Gradients following Algorithm 2 from Scales, 1987. 
	# The user provides an array of linear operators. Ensure linear operator(s) pass the dot product.


#=
	
	cost = Float64[]
	rand_string = string(round(Int,rand()*100000))
	g = join(["tmp_CG_g_",rand_string])
	s = join(["tmp_CG_s_",rand_string])
	rr = join(["tmp_CG_rr_",rand_string])
	tt = join(["tmp_CG_tt_",rand_string])
	r = join(["tmp_CG_r_",rand_string])
	t = join(["tmp_CG_t_",rand_string])
	SeisCopy(d,r)
	push!(cost,InnerProduct(r,r))
	fp = open(cost_file,"w")
	@compat write(fp,join(["began execution at: ",Libc.strftime(time()),"\n"]))
	write(fp,join([string(cost[1]),"\n"]))
	close(fp)
	LinearOperator(g,r,operators,parameters,adj=true)
	SeisCopy(g,s)
	SeisCopy(g,rr)
	SeisCopy(s,tt)
	CGStep(tt,s,a=0.,b=0.)
	CGStep(rr,g,a=0.,b=0.)
	SeisCopy(tt,m)
	gamma_old = InnerProduct(g,g)
	for iter = 1 : Niter	
		LinearOperator(s,t,operators,parameters,adj=false)
		CGStep(tt,s,a=0.,b=mu)
		delta = InnerProduct(t,t) + InnerProduct(tt,tt)
		alpha = gamma_old/delta
		CGStep(m,s,a=1.,b=alpha)  
		CGStep(r,t,a=1.,b=-alpha)
		CGStep(rr,tt,a=1.,b=-alpha)
		push!(cost,InnerProduct(r,r) + InnerProduct(rr,rr))
		fp = open(cost_file,"a")
		write(fp,join([string(cost[iter+1]),"\n"]))
		close(fp)
		LinearOperator(g,r,operators,parameters,adj=true)
		CGStep(g,rr,a=1.,b=mu)
		gamma = InnerProduct(g,g)
		beta = gamma/gamma_old
		gamma_old = copy(gamma)
		CGStep(s,g,a=beta,b=1.)
	end
	SeisRemove(g);
	SeisRemove(s);
	SeisRemove(rr);
	SeisRemove(tt)
	SeisRemove(r);
	SeisRemove(t);

=#


	cost = Float64[]
	rand_string = string(round(Int,rand()*100000))


	# approximate largest eigenvalue of the sequence of operators
	z = join(["tmp_CG_z_",rand_string])
	Lz = join(["tmp_CG_Lz_",rand_string])
	LTLz = join(["tmp_CG_LTLz_",rand_string])
	LinearOperator(z,d,operators,parameters,adj=true)
	tmp,h,ext = SeisRead(z)
	SeisWrite(z,ones(size(tmp)),h,ext)
	LinearOperator(z,Lz,operators,parameters,adj=false)
	LinearOperator(LTLz,Lz,operators,parameters,adj=true) 
	eig_val = InnerProduct(LTLz,LTLz)/InnerProduct(z,z)
	println("eig_val=",eig_val)



	s = join(["tmp_CG_s_",rand_string])
	p = join(["tmp_CG_p_",rand_string])
	r = join(["tmp_CG_r_",rand_string])
	q = join(["tmp_CG_q_",rand_string])
	SeisCopy(d,r)
	fp = open(cost_file,"w")
	@compat write(fp,join(["began execution at: ",Libc.strftime(time()),"\n"]))
	close(fp)
	LinearOperator(s,r,operators,parameters,adj=true)
	SeisCopy(s,p)
	SeisCopy(s,m)
	CGStep(m,s,a=0.,b=0.)
	norm_s0 = sqrt(InnerProduct(s,s))
#println("norm_s0=",norm_s0)	
	gamma = norm_s0^2
#println("gamma=",gamma)	
	norm_m = 0.0
	m_max = norm_m
	for iter = 1 : Niter	
		LinearOperator(p,q,operators,parameters,adj=false)
		delta = InnerProduct(q,q) + mu*InnerProduct(p,p)
#println("delta=",delta)	
		if delta <= 0 
			error("delta <= 0")
		end
		alpha = gamma / delta
#println("alpha=",alpha)	
		CGStep(m,p,a=1.,b=alpha)  
		CGStep(r,q,a=1.,b=-alpha)
		LinearOperator(s,r,operators,parameters,adj=true)
		CGStep(s,m,a=1.,b=-mu)
#println("mu=",mu)	
		norm_s = sqrt(InnerProduct(s,s))
#println("norm_s=",norm_s)	
		gamma1 = gamma
#println("gamma1=",gamma1)	
		gamma = norm_s^2
#println("gamma=",gamma)	
		beta = gamma / gamma1
#println("beta=",beta)	
		CGStep(p,s,a=beta,b=1.)
		norm_m = sqrt(InnerProduct(m,m))
#println("norm_m=",norm_m)	
		m_max = max(m_max,norm_m)
#println("m_max=",m_max)	
		push!(cost,norm_s / norm_s0)
#println("cost=",norm_s / norm_s0)	
		fp = open(cost_file,"a")
		write(fp,join([string(cost[iter]),"\n"]))
		close(fp)
		
#println("(norm_s <= norm_s0 * tol) || (norm_m * tol >= 1) ")	
#println(norm_s,"   vs.    ",norm_s0 * tol,"     ",norm_m * tol,"   vs.    ",1)	
		
		
		
		if (norm_s <= norm_s0 * tol) || (norm_m * tol >= 1)
			println("tolerance reached, ending at iteration ",iter)
			break;
		end	
	end
	
#println("shrink=",norm_m/m_max)	
	
	
	
	SeisRemove(p);
	SeisRemove(s);
	SeisRemove(r);
	SeisRemove(q)
	
	
	
	
end

function ConjugateGradients(m::Array{ASCIIString,1},d::Array{ASCIIString,1},operators,parameters,cost_file::ASCIIString;Niter=10,mu=[0 0],tol=1.0e-15)
	# Conjugate Gradients following Algorithm 2 from Scales, 1987. 
	# The user provides an array of linear operators. Ensure linear operator(s) pass the dot product.

	cost = Float64[]
	rand_string = string(round(Int,rand()*100000))
	s = [join(["tmp_CG_s1_",rand_string]);join(["tmp_CG_s2_",rand_string])]
	p = [join(["tmp_CG_p1_",rand_string]);join(["tmp_CG_p2_",rand_string])]
	r = [join(["tmp_CG_r1_",rand_string]);join(["tmp_CG_r2_",rand_string])]
	q = [join(["tmp_CG_q1_",rand_string]);join(["tmp_CG_q2_",rand_string])]
	SeisCopy(d,r)
	fp = open(cost_file,"w")
	@compat write(fp,join(["began execution at: ",Libc.strftime(time()),"\n"]))
	close(fp)
	LinearOperator(s,r,operators,parameters,adj=true)
	SeisCopy(s,p)
	SeisCopy(s,m)
	CGStep(m,s,a=[0.;0.],b=[0.;0.])
	norm_s0 = sqrt(InnerProduct(s,s))
	gamma = norm_s0^2
	norm_m = 0.0
	m_max = norm_m
	for iter = 1 : Niter	
		LinearOperator(p,q,operators,parameters,adj=false)
		delta = InnerProduct(q,q) + mu[1]*InnerProduct(p[1],p[1]) + mu[2]*InnerProduct(p[2],p[2])
		if delta <= 0 
			error("delta <= 0")
		end
		alpha = gamma / delta
		CGStep(m,p,a=[1.0;1.0],b=[alpha;alpha])  
		CGStep(r,q,a=[1.0;1.0],b=[-alpha;-alpha]) 
		LinearOperator(s,r,operators,parameters,adj=true)
		CGStep(s,m,a=[1.0;1.0],b=-mu) 
		norm_s = sqrt(InnerProduct(s,s))
		gamma1 = gamma
		gamma = norm_s^2
		beta = gamma / gamma1
		CGStep(p,s,a=[beta;beta],b=[1.0;1.0])
		norm_m = sqrt(InnerProduct(m,m))
		m_max = max(m_max,norm_m)
		push!(cost,norm_s / norm_s0)
		fp = open(cost_file,"a")
		write(fp,join([string(cost[iter]),"\n"]))
		close(fp)
		if (norm_s <= norm_s0 * tol) || (norm_m * tol >= 1)
			println("tolerance reached, ending at iteration ",iter)
			break;
		end	
	end
	SeisRemove(p);
	SeisRemove(s);
	SeisRemove(r);
	SeisRemove(q)
	
end
