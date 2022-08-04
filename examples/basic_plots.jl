### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d990be32-1343-11ed-0acb-0fdad105c87c
begin
	import Pkg
	Pkg.activate()
	using JSServe
	using WGLMakie
	using MRSimulator
	using PlutoUI
    Page()
end

# ╔═╡ ecac5fd3-313c-407d-a299-6ec1825dd9c8
@bind TE PlutoUI.Slider(5:50, default=20, show_value=true)

# ╔═╡ 9512c409-850c-4b7d-ac5a-ecbacf706f9d
begin
	sequence = Sequence([RFPulse(flip_angle=90), RFPulse(flip_angle=180, time=TE/2.)], 1.5 * TE)
	MRSimulator.sequenceplot(sequence)
end

# ╔═╡ ee4b67ce-95a1-4a89-968a-210e1d19f3ed
begin
	snap = Snapshot([
		Spin(transverse=0.3), 
		Spin(position=[0, 1, 2.], transverse=1.),
		Spin(position=[2., 1, 2.], transverse=0.)], 0.3)
	plot(snap)
end

# ╔═╡ Cell order:
# ╠═d990be32-1343-11ed-0acb-0fdad105c87c
# ╠═ecac5fd3-313c-407d-a299-6ec1825dd9c8
# ╠═9512c409-850c-4b7d-ac5a-ecbacf706f9d
# ╠═ee4b67ce-95a1-4a89-968a-210e1d19f3ed
