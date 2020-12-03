using TightlyBound
using Plots

#setup chosen backend
#pyplot()

# load tb_crys from file

filname = "../test/data_forces/znse.in_vnscf_vol_2/projham.xml.gz"
tbc = read_tb_crys(filname)

println(tbc)

#load dft calculation
filname_dft = "../test/data_forces/znse.in_vnscf_vol_2/qe.save/"
dft = TightlyBound.QE.loadXML(filname_dft)

println(dft)

#compare to dft calculation
plot_compare_dft(tbc, dft)

savefig("mgs_compare_dft.pdf")



