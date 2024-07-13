#-------------------------
#   save
#-------------------------
includet("../src/sv.jl")
using .sv


foldername = "20240706_1"


path = pwd()*"/demo_scripts/data/"*foldername
mkdir(path)

notes = "\nNotes:"
save_variables(joinpath(path, "var.txt"), kx, δ, A, cloud_r, cloud_T,n,notes)  
save_plots(path, pPhase,p2d,p1d_1,p1d_2)

print("\nSaved!")