var documenterSearchIndex = {"docs":
[{"location":"B/#BBB","page":"BB","title":"BBB","text":"","category":"section"},{"location":"B/#BBB1","page":"BB","title":"BBB1","text":"","category":"section"},{"location":"B/","page":"BB","title":"BB","text":"stuff\nstuff2","category":"page"},{"location":"B/#BBB2","page":"BB","title":"BBB2","text":"","category":"section"},{"location":"B/","page":"BB","title":"BB","text":"stuffx\nstuffy","category":"page"},{"location":"theindex/#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"theindex/","page":"Index","title":"Index","text":"","category":"page"},{"location":"#TightlyBound.jl-Documentation","page":"Home","title":"TightlyBound.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A three-body tight binding program written in Julia","category":"page"},{"location":"","page":"Home","title":"Home","text":"Primary Author: Kevin F. Garrity, NIST","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nThis package currently barely works. Needs much more testing.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Run tight-binding calculations with near DFT level accuracy (PBEsol).\nGet results in seconds based on pre-fit parameters from across periodic table.\nCalculate bandstructures and total energies.\nGet forces, stresses, and relax structures.\nParameters based on two- and three-body interactions.\nIncludes self-consistent treatment of long-range Coulomb interaction.\nPlotting based on interface of Plots.jl","category":"page"},{"location":"#Outline","page":"Home","title":"Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Stuff\nOther stuff","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = TightlyBound","category":"page"},{"location":"","page":"Home","title":"Home","text":"TightlyBound\nmakecrys\nrelax_structure","category":"page"},{"location":"#TightlyBound.TightlyBound","page":"Home","title":"TightlyBound.TightlyBound","text":"holds three body tight binding important stuff\n\n\n\n\n\n","category":"module"},{"location":"#TightlyBound.CrystalMod.makecrys","page":"Home","title":"TightlyBound.CrystalMod.makecrys","text":"makecrys(A,coords,types)\n\nReturn a crystal object from 3×3 lattice A (Bohr), nat × 3 coords in crystal units, and nat element strings\n\nNote: also export-ed directly from TightlyBound for convenience\n\n#Examples\n\njulia> makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0.0 0.0 0.0], [\"H\"])\nA1=     10.00000  0.00000  0.00000\nA2=     0.00000  10.00000  0.00000\nA3=     0.00000  0.00000  10.00000\n\nH    0.00000  0.00000  0.00000\n\n\n\n\n\nmakecrys(filename::String)\n\nRead filename, return crystal object. File can be POSCAR, or simple quantum espresso inputfile\n\n\n\n\n\nmakecrys(lines::Array{String,1})\n\nRead string array, return crystal object. string array can be POSCAR, or simple quantum espresso inputfile\n\n\n\n\n\n","category":"function"},{"location":"#TightlyBound.relax_structure","page":"Home","title":"TightlyBound.relax_structure","text":"relax_structure(c::crystal; mode=\"vc-relax\")\n\nFind the lowest energy atomic configuration of crystal c.\n\n...\n\nArguments\n\nc::crystal: the structure to relax, only required argument\nmode=\"vc-relax\": Default (variable-cell relax) will relax structure and cell, anything else will relax structure only.\ndatabase=missing: coefficent database, default is to use the pre-fit pbesol database\nsmearing=0.01: smearing temperature (ryd), default = 0.01\ngrid=missing: k-point grid, e.g. [10,10,10], default chosen automatically\nnsteps=100: maximum iterations\nupdate_grid=true: update automatic k-point grid during relaxation\n\n...\n\n\n\n\n\n","category":"function"},{"location":"A/#AAA","page":"AA","title":"AAA","text":"","category":"section"},{"location":"A/#AAA1","page":"AA","title":"AAA1","text":"","category":"section"},{"location":"A/","page":"AA","title":"AA","text":"stuff\nstuff2","category":"page"},{"location":"A/#AAA2","page":"AA","title":"AAA2","text":"","category":"section"},{"location":"A/","page":"AA","title":"AA","text":"stuffx\nstuffy","category":"page"}]
}
