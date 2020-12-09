# Running tight-binding calculations

How to run tight-binding calculations using the pre-fit tight-binding
coefficients. Note, only elemental and binary systems are currently
supported.

!!! note

    Running a julia function for the first time will compile the function. Future runs will be *much* faster.

## **Create a crystal object**

Consists of lattice vectors, atomic positions, and atom types. Current units are Bohr only, this will change.

```@example 1
using TightlyBound
A = [4.0 4.0 0.0;4.0 0.0 4.0;0.0 4.0 4.0];
pos = [0.0 0.0 0.0];
types =        ["Al"];
fcc_al = makecrys(A, pos, types)
```

Alternatively, you can read the positions from a simple POSCAR or Quantum Espresso input file.

```@example 1
rbcl = makecrys("../src/POSCAR_rbcl")
```

## **Do an self-consistent calculation.**

Gets the energy and charge density.

```@example 1
alp = makecrys("../src/POSCAR_alp")
energy, tbc_alp = scf_energy(alp); 
println("The energy is $energy Ryd")
```
This returns the non-magnetic atomization energy, and a tight-binding object with the SCF electron density calculated.

## **Plot the band structure.**

Using the tight-binding object from above. Note: SCF must be done first to get meaningful results.

```@example 1
using Plots #hide
gr() #hide
plot_bandstr(tbc_alp); 
savefig("alp.png"); #hide
```

![AlP plot](alp.png)

The default just picks some random kpoints, but you can add your own kpath. We also project onto the *s* orbital of Al.

```@example 1
kpath=[0.0 0.0 0.0; 0.5 0.5 0.5; 0.0 0.5 0.5];
knames=["Γ", "X", "V"];
plot_bandstr(tbc_alp, kpath=kpath, names=knames, npts=100, proj_orbs=[:s], proj_types=["Al"]);
savefig("alp2.png"); #hide
```

![AlP plot](alp2.png)

## **Calculate force / stress**

```@example 1
energy, force, stress, tbc = scf_energy_force_stress(tbc_alp);

println("energy $energy")
println()
println("Forces")
show(stdout, "text/plain", force)
println()
println("Stress")
show(stdout, "text/plain", stress)
nothing #hide
```
Can also be called directly on a new crystal structure instead of a tb_crys object.

## **Relax structure**

```@example 1
crys_new, tbc_updated, energy, force, stress = relax_structure(alp);

println("Energy new $energy")
println()
println("Force")
show(stdout, "text/plain", force)
println()
println("Stress")
show(stdout, "text/plain", stress)
nothing #hide
```
Energy is lower, stress is near zero, forces are zero by symmetry in Zinc Blende structure.


