using TightlyBound
using Plots

#setup chosen Plots backend, i recommend pyplot
#pyplot()


#make the crystal object
#here we choose Sc P rocksalt

types=["Sc", "P"];

#positions, crystal units
pos = [0 0 0 ; 0.5000000000  0.5000000000  0.5000000000]
#lattice vectors, in Bohr units currently
A=[ [4.7 4.7 0]; [4.7 0 4.7 ]; [ 0 4.7 4.7]];

#makes the crystal
c=makecrys(A, pos, types)

println("Starting crystal is")
println(c)


energy, tbc, flag = scf_energy(c)

println("Starting energy is $energy")
println("--------------------------")

cfinal, tbc, energy, forces, stress = relax_structure(c)

println("Final Energy $energy ; DFT energy is -0.9836659797393565 Ryd")

println("Final crystal is ")
println(cfinal)
print()
println("forces")
println(forces)
print()
println("stress")
println(stress)
print()


#plot band structure
kpath=[0.0 0.0 0.0;0.5 0 0; 0.5 0.25 0.75; 0.0 0.0 0.0]
knames=["Γ", "X", "W", "Γ"]

plot_bandstr(tbc, kpath=kpath, names=knames, proj_types=["P"])

savefig("ScP_rocksalt.pdf")



