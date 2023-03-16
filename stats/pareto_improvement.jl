include("hv_indicator.jl")

_sysid=1
#_expid = "original_system_$(_sysid)_objective_small"
_expid = "original_system_$(_sysid)_objective"
_objectives_labels = ["size", "exercises", "failrate"]
_objectives_tomax = [(x) -> -x, identity, identity]
_runtime = 300
_caseid = "System $(_sysid)_1_1.0"

_p = plot_hypervolumes(_expid, _objectives_labels, _objectives_tomax, _runtime, _caseid)
title!(_p, "System $(_sysid) (objective)")
savefig(_p, "sys1objall.pdf")
#_p
