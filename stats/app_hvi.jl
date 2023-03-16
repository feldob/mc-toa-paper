
function dominates(a, b)
	d = false
	for i = 1:length(a)
		if (a[i] > b[i])
			return false
		elseif (a[i] < b[i])
			d = true
		end
	end
    return d
end

function reduce_to_non_dominating(existing_pf, new_candidates)
	front = copy(existing_pf)
	for i in 1:size(new_candidates, 2)
		nondominated = true
		a = new_candidates[:,i]
		for j in 1:size(existing_pf, 2)
			b = existing_pf[:,j]
			if dominates(a,b)
				nondominated = false
				break
			end
		end
		for j in 1:size(new_candidates, 2)
			if i == j
				continue
			end
			b = new_candidates[:,j]
			if dominates(a,b)
				nondominated = false
				break
			end
		end
		if nondominated
			front = hcat(front, a)
		end
	end
	return front
end

function reduce_to_non_dominating(partprep)
	front = Array{Float64}(undef, size(partprep, 1), 0)
	for i in 1:size(partprep, 2)
		nondominated = true
		a = partprep[:,i]
		for j in 1:size(partprep, 2)
			if i == j
				continue
			end
			b = partprep[:,j]
			if dominates(a,b)
				nondominated = false
				break
			end
		end
		if nondominated
			front = hcat(front, a)
		end
	end
	return front
end

function approximate_hypervolume_ms(F, lb = zeros(size(F, 2)), ub=maximum(F, dims=2), samples=10000)
    n_entries = size(F, 2)
    n_dims = size(F, 1)

	F_samples = repeat(lb,1,samples) + rand(n_dims,samples) .* repeat((ub - lb),1,samples)
	is_dominated_count = 0
	for i = 1:samples
		for j = 1:n_entries
			if (dominates(F[:,j], F_samples[:,i]))
				is_dominated_count = is_dominated_count + 1
				break
			end
		end
	end
	return prod(ub - lb) * (is_dominated_count / samples)
end

function approximate_hypervolume_ms_new(F, lb = zeros(size(F, 2)), ub=maximum(F, dims=2), samples=10_000)
    n_entries = size(F, 2)
    n_dims = size(F, 1)

	is_dominated_count = 0
	for i = 1:samples
        C = lb + rand(n_dims) .* (ub - lb)
		for j = 1:n_entries
			if (dominates(F[:,j], C))
				is_dominated_count = is_dominated_count + 1
				break
			end
		end
	end
	return prod(ub - lb) * (is_dominated_count / samples)
end