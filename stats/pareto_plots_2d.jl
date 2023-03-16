include("../io.jl")

sys1_obj = "original_system_1_objective"
_p = plot_2ds(sys1_obj, "sys1_obj_2d.pdf", "System 1, objective", 107, 2, 32, 17, 92)

sys1_constr = "original_system_1_constraint"
_p = plot_2ds(sys1_constr, "sys1_constr_2d.pdf", "System 1, constraints", 136, 181, 166, 151, 197)

sys2_obj = "original_system_2_objective"
_p = plot_2ds(sys2_obj, "sys2_obj_2d.pdf", "System 2, objectives", 111, 2, 32, 17, 92)

sys2_constr = "original_system_2_constraint"
_p = plot_2ds(sys2_constr, "sys2_constr_2d.pdf", "System 2, constraints", 136, 181, 166, 151, 197)
