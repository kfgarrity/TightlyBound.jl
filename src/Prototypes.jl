using LinearAlgebra
using ..Atomdata:atom_prefered_oxidation
using ..Atomdata:min_dimer_dist_dict
using ..QE:loadXML
struct proto_data

    CalcD::Dict
    core_mono::Array{String}
    core_binary::Array{String}

    A0::Array{String}
    A1B1::Array{String}
    A1B2::Array{String}
    A1B3::Array{String}
    A1B4::Array{String}
    A1B5::Array{String}
    A1B6::Array{String}
    A2B3::Array{String}
    A2B5::Array{String}
    metals::Array{String}
    short_bonds::Array{String}
    all_ternary::Array{String}

end

function setup_proto_data()

    CalcD = Dict()


    CalcD["atom"] = ["../reference_structures/atom.in", "scf", "all", "scf", "nscf"]
    CalcD["sc"] = ["../reference_structures/sc.in.up", "vc-relax", "all", "vol-big", "nscf"]
    CalcD["sc_inv"] = ["../reference_structures/sc.in.up", "vc-relax", "all", "break_inv", "nscf"]
    CalcD["bcc"] = ["../reference_structures/bcc.in.up", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["bcc_inv"] = ["../reference_structures/bcc_atom2.in", "vc-relax", "all", "break_inv", "nscf"]
    CalcD["fcc"] = ["../reference_structures/fcc.in.up", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["line"] = ["../reference_structures/line.in.up", "vc-relax", "z", "vol", "nscf"]
    CalcD["line_rumple"] = ["../reference_structures/line.in.rumple.up", "vc-relax", "z", "vol", "nscf"]


    CalcD["sc_verydense"] = ["../reference_structures/sc.in.up", "vc-relax", "all", "vol-verydense", "nscf"]

    CalcD["fcc_verydense"] = ["../reference_structures/fcc.in.up", "vc-relax", "all", "vol-verydense", "nscf"]
    CalcD["bcc_verydense"] = ["../reference_structures/bcc.in.up", "vc-relax", "all", "vol-verydense", "nscf"]
    CalcD["diamond_verydense"] = ["../reference_structures/diamond.in.up", "vc-relax", "all", "vol-verydense", "nscf"]


    CalcD["fcc_dense"] = ["../reference_structures/fcc.in.up", "vc-relax", "all", "vol-dense", "nscf"]
    CalcD["bcc_dense"] = ["../reference_structures/bcc.in.up", "vc-relax", "all", "vol-dense", "nscf"]

    CalcD["hcp"] = ["../reference_structures/hcp.in.up", "vc-relax", "all", "vol", "nscf"]
    CalcD["diamond"] = ["../reference_structures/diamond.in.up", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["graphene"] = ["../reference_structures/fake_graphene.in", "vc-relax", "2Dxy", "2D-mid", "nscf"]
    CalcD["hex"] = ["../reference_structures/hex.in.up", "vc-relax", "2Dxy", "2D-mid", "nscf"]
    CalcD["hex_short"] = ["../reference_structures/hex.in.up", "vc-relax", "2Dxy", "2D-short", "nscf"]
    CalcD["square"] = ["../reference_structures/square.in.up", "vc-relax", "2Dxy", "2D", "nscf"]
    CalcD["dimer"] =       ["../reference_structures/dimer.in", "relax", "2Dxy", "coords", "nscf"]
    CalcD["dimer_short"] = ["../reference_structures/dimer.in", "relax", "2Dxy", "coords-short", "nscf"]
    CalcD["dimer_small"] = ["../reference_structures/dimer_small.in", "relax", "2Dxy", "coords", "nscf"]
    CalcD["dimer_super"] = ["../reference_structures/dimer.in", "relax", "2Dxy", "coords_super", "nscf"]
    CalcD["hcp_shape"] = ["../reference_structures/hcp.in.up", "vc-relax", "all", "shape", "nscf"]

    CalcD["hex_2lay"] = ["../reference_structures/hex_2layers.in.up", "vc-relax", "2Dxy", "2D", "nscf"]
    CalcD["bcc_2lay"] = ["../reference_structures/bcc_2layers.in.up", "vc-relax", "2Dxy", "2D", "nscf"]

    CalcD["fcc_huge"] = ["../reference_structures/fcc.in.up", "vc-relax", "all", "vol-huge", "nscf"]

    CalcD["bcc_tet"] = ["../reference_structures/bcc_tet.in", "vc-relax", "all", "vol", "nscf"]


    CalcD["trimer"] =       ["../reference_structures/trimer.in", "none", "2Dxy", "coords_trimer", "nscf"]
    CalcD["trimer_dense"] =       ["../reference_structures/trimer.in", "none", "2Dxy", "coords_trimer_dense", "nscf"]
    CalcD["trimer2"] =       ["../reference_structures/trimer.in2", "none", "2Dxy", "coords_trimer2", "nscf"]

    CalcD["trimer_ab2"] =       ["../reference_structures/binary/trimer.in.ab2", "none", "2Dxy", "coords_trimer_ab", "nscf"]
    CalcD["trimer2_ab2"] =       ["../reference_structures/binary/trimer.in2.ab2", "none", "2Dxy", "coords_trimer_ab", "nscf"]

    CalcD["trimer_ba2"] =       ["../reference_structures/binary/trimer.in.ba2", "none", "2Dxy", "coords_trimer_ab", "nscf"]
    CalcD["trimer2_ba2"] =       ["../reference_structures/binary/trimer.in2.ba2", "none", "2Dxy", "coords_trimer_ab", "nscf"]


    CalcD["trimer_ab2_dense"] =       ["../reference_structures/binary/trimer.in.ab2", "none", "2Dxy", "coords_trimer_ab_dense", "nscf"]
    CalcD["trimer2_ab2_dense"] =       ["../reference_structures/binary/trimer.in2.ab2", "none", "2Dxy", "coords_trimer_ab_dense", "nscf"]

    CalcD["trimer_ba2_dense"] =       ["../reference_structures/binary/trimer.in.ba2", "none", "2Dxy", "coords_trimer_ab_dense", "nscf"]
    CalcD["trimer2_ba2_dense"] =       ["../reference_structures/binary/trimer.in2.ba2", "none", "2Dxy", "coords_trimer_ab_dense", "nscf"]


    CalcD["trimer_ba2_big"] =       ["../reference_structures/binary/trimer.in.ba2.big", "none", "2Dxy", "coords_trimer_ab_big", "nscf"]


    CalcD["as_221"] = ["../reference_structures/POSCAR_As_221", "vc-relax", "all", "vol", "nscf"]
    CalcD["as_orth"] = ["../reference_structures/POSCAR_As_ortho", "vc-relax", "all", "vol", "nscf"]
    CalcD["ga_tet"] = ["../reference_structures/POSCAR_ga_tet", "vc-relax", "all", "vol", "nscf"]
    CalcD["ge_wurtz"] = ["../reference_structures/POSCAR_ge_wurtz", "vc-relax", "all", "vol", "nscf"]
    CalcD["pb_r3m"] = ["../reference_structures/POSCAR_pb_r3m_2atom", "vc-relax", "all", "vol", "nscf"]
    CalcD["beta_sn"] = ["../reference_structures/beta_sn.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["n"] = ["../reference_structures/n.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["f"] = ["../reference_structures/f.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["ga"] = ["../reference_structures/POSCAR_ga", "vc-relax", "all", "vol", "nscf"]
    CalcD["bi"] = ["../reference_structures/bi.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["te"] = ["../reference_structures/te.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["in"] = ["../reference_structures/in.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["i2"] = ["../reference_structures/i2.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["li_p6mmm"] = ["../reference_structures/POSCAR_li_p6mmm", "vc-relax", "all", "vol", "nscf"]

    CalcD["sc_shape"] = ["../reference_structures/sc.in.up", "vc-relax", "all", "shape", "nscf"]
    CalcD["diamond_shear"] = ["../reference_structures/diamond.in.up", "vc-relax", "all", "shear", "nscf"]

    CalcD["hh_mono"] = ["../reference_structures/POSCAR_hh_mono", "vc-relax", "all", "vol-mid", "nscf"]


    CalcD["cscl"] = ["../reference_structures/binary/cscl.in", "vc-relax", "all", "vol-big", "nscf"]
    CalcD["cscl_layers"] = ["../reference_structures/binary/cscl_layers.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["cscl_inv"] = ["../reference_structures/binary/cscl.in", "vc-relax", "all", "break_inv", "nscf"]
    CalcD["rocksalt"] = ["../reference_structures/binary/rocksalt.in", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["rocksalt_inv"] = ["../reference_structures/binary/rocksalt.in", "vc-relax", "all", "break_inv", "nscf"]

    CalcD["rocksalt_dense"] = ["../reference_structures/binary/rocksalt.in", "vc-relax", "all", "vol-dense", "nscf"]
    CalcD["hcp_v2_dense"] = ["../reference_structures/binary/hcp.in2", "vc-relax", "all", "vol-dense", "nscf"]
    CalcD["znse_dense"] = ["../reference_structures/binary/znse.in", "vc-relax", "all", "vol-dense", "nscf"]

    CalcD["nias"] = ["../reference_structures/binary/POSCAR_nias_hex", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["227"] = ["../reference_structures/binary/POSCAR_227", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["laga"] = ["../reference_structures/binary/POSCAR_laga_63", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["distort"] = ["../reference_structures/binary/POSCAR_distort", "vc-relax", "all", "vol", "nscf"]
    CalcD["distort_ab2"] = ["../reference_structures/binary/POSCAR_distort_ab2", "vc-relax", "all", "vol", "nscf"]
    CalcD["distort_ba2"] = ["../reference_structures/binary/POSCAR_distort_ba2", "vc-relax", "all", "vol", "nscf"]

##    CalcD["227_BA"] = ["../reference_structures/binary/POSCAR_227_BA", "vc-relax", "all", "vol-mid", "nscf"]


    CalcD["znse_shear"] = ["../reference_structures/binary/znse.in", "vc-relax", "all", "shear", "nscf"]
    CalcD["hbn"] = ["../reference_structures/binary/hbn.in", "vc-relax", "2Dxy", "2D-mid", "nscf"]
    CalcD["znse"] = ["../reference_structures/binary/znse.in", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["dimer2"] = ["../reference_structures/binary/dimer.in", "relax", "all", "coords", "nscf"]
    CalcD["dimer2_rev"] = ["../reference_structures/binary/dimer_rev.in", "relax", "all", "coords", "nscf"]
    CalcD["square2"] = ["../reference_structures/binary/square.in", "vc-relax", "2Dxy", "2D", "nscf"]
    CalcD["caf2"] = ["../reference_structures/binary/POSCAR_caf2", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["co2"] = ["../reference_structures/binary/co2.in", "relax", "all", "coords-small", "nscf"]
    CalcD["square_ab2"] = ["../reference_structures/binary/square_ab2.in", "vc-relax", "2Dxy", "2D", "nscf"]
    CalcD["mgf2"] = ["../reference_structures/binary/POSCAR_mgf2", "vc-relax", "all", "vol", "nscf"]
    CalcD["hcp_v2"] = ["../reference_structures/binary/hcp.in2", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["p2ca3"] = ["../reference_structures/binary/POSCAR_p2ca3", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["triangle"] = ["../reference_structures/binary/triangle.in", "relax", "all", "coords-small", "nscf"]
    CalcD["triangle2"] = ["../reference_structures/binary/triangle2.in", "relax", "all", "coords-small", "nscf"]


    CalcD["beta_sn2"] = ["../reference_structures/binary/beta_sn.in.up", "vc-relax", "all", "vol", "nscf"]
    CalcD["hbn_real"] = ["../reference_structures/binary/hbn_real.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["znseAAAB"] = ["../reference_structures/binary/znse.in.super.AAAB", "vc-relax", "all", "vol", "nscf"]
    CalcD["rocksaltAAAB"] = ["../reference_structures/binary/rocksalt.in.super.AAAB", "vc-relax", "all", "vol", "nscf"]
    CalcD["znseABBB"] = ["../reference_structures/binary/znse.in.super.ABBB", "vc-relax", "all", "vol", "nscf"]
    CalcD["rocksaltABBB"] = ["../reference_structures/binary/rocksalt.in.super.ABBB", "vc-relax", "all", "vol", "nscf"]
    CalcD["al2o3"] = ["../reference_structures/binary/POSCAR_al2o3", "vc-relax", "all", "vol", "nscf"]
    CalcD["sis2"] = ["../reference_structures/binary/POSCAR_sis2", "vc-relax", "all", "vol", "nscf"]
    CalcD["tio2_rutile"] = ["../reference_structures/binary/POSCAR_tio2_rutile", "vc-relax", "all", "vol", "nscf"]
    CalcD["bi2se3"] = ["../reference_structures/binary/POSCAR_bi2se3", "vc-relax", "all", "vol", "nscf"]
    CalcD["sis_2d"] = ["../reference_structures/binary/sis_2d.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["sio2_224"] = ["../reference_structures/binary/sio2.224.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["ges"] = ["../reference_structures/binary/POSCAR_ges", "vc-relax", "all", "vol", "nscf"]
    CalcD["alf3"] = ["../reference_structures/binary/POSCAR_alf3", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["mgf2_v2"] = ["../reference_structures/binary/POSCAR_mgf2_v2", "vc-relax", "all", "vol", "nscf"]
    CalcD["sns"] = ["../reference_structures/binary/POSCAR_sns", "vc-relax", "all", "vol", "nscf"]
    CalcD["sns2"] = ["../reference_structures/binary/POSCAR_sns2", "vc-relax", "all", "vol", "nscf"]
    CalcD["y2o3"] = ["../reference_structures/binary/POSCAR_y2o3", "vc-relax", "all", "vol", "nscf"]
    CalcD["wurtz"] = ["../reference_structures/binary/POSCAR_wurtz", "vc-relax", "all", "vol", "nscf"]
    CalcD["mgcl2"] = ["../reference_structures/binary/POSCAR_mgcl2", "vc-relax", "all", "vol", "nscf"]
    CalcD["mgcl2_tet"] = ["../reference_structures/binary/POSCAR_mgcl2_tet", "vc-relax", "all", "vol", "nscf"]
    CalcD["asna3_2d"] = ["../reference_structures/binary/POSCAR_asna3", "vc-relax", "2Dxy", "2D", "nscf"]
    CalcD["gain3"] = ["../reference_structures/binary/POSCAR_gain3", "vc-relax", "all", "vol", "nscf"]

    CalcD["li3n_hex"] = ["../reference_structures/binary/POSCAR_li3n", "vc-relax", "all", "vol", "nscf"]
    CalcD["nan3"] = ["../reference_structures/binary/POSCAR_nan3", "vc-relax", "all", "vol", "nscf"]
    CalcD["rbo2"] = ["../reference_structures/binary/POSCAR_rbo2", "vc-relax", "all", "vol", "nscf"]



    CalcD["nb2o5"] = ["../reference_structures/binary/POSCAR_nb2o5", "vc-relax", "all", "vol", "nscf"]
    CalcD["sif4"] = ["../reference_structures/binary/POSCAR_sif4", "vc-relax", "all", "vol", "nscf"]
    CalcD["ticl2"] = ["../reference_structures/binary/POSCAR_ticl2", "vc-relax", "all", "vol", "nscf"]

    CalcD["snf4"] = ["../reference_structures/binary/POSCAR_snf4", "vc-relax", "all", "vol", "nscf"]
    CalcD["ga2s3"] = ["../reference_structures/binary/POSCAR_ga2s3", "vc-relax", "all", "vol", "nscf"]

    CalcD["bcc_13"] = ["../reference_structures/binary/POSCAR_bcc_13", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["bcc_31"] = ["../reference_structures/binary/POSCAR_bcc_31", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["mg2si"] = ["../reference_structures/binary/POSCAR_mg2si", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["simg2"] = ["../reference_structures/binary/POSCAR_simg2", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["fcc_12"] = ["../reference_structures/binary/fcc_12.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["fcc_21"] = ["../reference_structures/binary/fcc_21.in", "vc-relax", "all", "vol", "nscf"]

    CalcD["fcc_conv_ABBB"] = ["../reference_structures/binary/fcc_conv_ABBB.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["fcc_conv_BAAA"] = ["../reference_structures/binary/fcc_conv_BAAA.in", "vc-relax", "all", "vol", "nscf"]

    CalcD["mgb2_12"] = ["../reference_structures/binary/POSCAR_mgb2", "vc-relax", "all", "vol", "nscf"]
    CalcD["mgb2_21"] = ["../reference_structures/binary/POSCAR_mgb2_21", "vc-relax", "all", "vol", "nscf"]

    CalcD["mgb2_AB"] = ["../reference_structures/binary/POSCAR_mgb2_AB", "vc-relax", "all", "vol", "nscf"]
    CalcD["mgb2_BA"] = ["../reference_structures/binary/POSCAR_mgb2_BA", "vc-relax", "all", "vol", "nscf"]

    CalcD["ab2_71"] = ["../reference_structures/binary/POSCAR_ab2_71", "vc-relax", "all", "vol", "nscf"]
    CalcD["ba2_71"] = ["../reference_structures/binary/POSCAR_ba2_71", "vc-relax", "all", "vol", "nscf"]

    CalcD["irn2_38"] = ["../reference_structures/binary/POSCAR_irn2_38", "vc-relax", "all", "vol", "nscf"]
    CalcD["cao2_12"] = ["../reference_structures/binary/POSCAR_cao2_12", "vc-relax", "all", "vol", "nscf"]

    CalcD["p2o5"] = ["../reference_structures/binary/POSCAR_p2o5", "vc-relax", "all", "vol", "nscf"]


    CalcD["znseAABB"] = ["../reference_structures/binary/znse.AABB.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["squareAABB"] = ["../reference_structures/binary/square.in.AABB", "vc-relax", "2Dxy", "2D", "nscf"]


    CalcD["dimer_pair"] = ["../reference_structures/binary/POSCAR_dimer_pair", "vc-relax", "all", "vol", "nscf"]


    CalcD["mgb2_AB-mid"] = ["../reference_structures/binary/POSCAR_mgb2_AB", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["mgb2_BA-mid"] = ["../reference_structures/binary/POSCAR_mgb2_BA", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["mgb2_12-mid"] = ["../reference_structures/binary/POSCAR_mgb2", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["mgb2_21-mid"] = ["../reference_structures/binary/POSCAR_mgb2_21", "vc-relax", "all", "vol-mid", "nscf"]


#    CalcD["i4mmm_tet_AB"] = ["../reference_structures/binary/POSCAR_i4mmm_4atom_AB", "vc-relax", "all", "vol", "nscf"]
#    CalcD["i4mmm_tet_BA"] = ["../reference_structures/binary/POSCAR_i4mmm_4atom_BA", "vc-relax", "all", "vol", "nscf"]



    CalcD["gei2"] = ["../reference_structures/binary/POSCAR_GeI2", "vc-relax", "all", "vol", "nscf"]
    CalcD["gei4_mol"] = ["../reference_structures/binary/gei4_molecule.in", "relax", "all", "coords-small", "nscf"]

    CalcD["ab3_mol"] = ["../reference_structures/binary/ab3_molecule.in", "relax", "all", "coords-small", "nscf"]


    CalcD["tet"] = ["../reference_structures/binary/tet.in", "vc-relax", "all", "vol", "nscf"]


    CalcD["mof6"] = ["../reference_structures/binary/POSCAR_mof6", "vc-relax",  "all", "vol", "nscf"]


    CalcD["ascl5"] = ["../reference_structures/binary/POSCAR_ascl5", "vc-relax",  "2Dxy", "2D", "nscf"]
    CalcD["rocksalt_shape"] = ["../reference_structures/binary/rocksalt.in", "vc-relax", "all", "shape", "nscf"]
    CalcD["rocksalt_2lay"] = ["../reference_structures/binary/rocksalt.in.2lay", "vc-relax", "2Dxy", "2D", "nscf"]

    CalcD["gas"] = ["../reference_structures/binary/POSCAR_gas", "vc-relax", "all", "vol", "nscf"]

    CalcD["hex12"] = ["../reference_structures/binary/hex_trim_12.in", "vc-relax",  "2Dxy", "2D", "nscf"]
    CalcD["hex21"] = ["../reference_structures/binary/hex_trim_21.in", "vc-relax",  "2Dxy", "2D", "nscf"]

    CalcD["hex12a"] = ["../reference_structures/binary/hex_trim_12.in", "vc-relax",  "2Dxy", "2D-mid", "nscf"]
    CalcD["hex21a"] = ["../reference_structures/binary/hex_trim_21.in", "vc-relax",  "2Dxy", "2D-mid", "nscf"]

    #ternary
    CalcD["abc_line"] = ["../reference_structures/ternary/abc_line.in", "relax", "all", "coords-small", "nscf"]
    CalcD["bac_line"] = ["../reference_structures/ternary/bac_line.in", "relax", "all", "coords-small", "nscf"]
    CalcD["cab_line"] = ["../reference_structures/ternary/cab_line.in", "relax", "all", "coords-small", "nscf"]

    CalcD["fcc_tern"] = ["../reference_structures/ternary/fcc_tern.in", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["hex_trim"] = ["../reference_structures/ternary/hex_trim_3.in", "vc-relax",  "2Dxy", "2D-mid", "nscf"]

    CalcD["hh1"] = ["../reference_structures/ternary/POSCAR_hh1", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["hh2"] = ["../reference_structures/ternary/POSCAR_hh2", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["hh3"] = ["../reference_structures/ternary/POSCAR_hh3", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["stuffhex_1"] = ["../reference_structures/ternary/POSCAR_mgb2_1", "vc-relax", "all", "vol", "nscf"]
    CalcD["stuffhex_2"] = ["../reference_structures/ternary/POSCAR_mgb2_2", "vc-relax", "all", "vol", "nscf"]
    CalcD["stuffhex_3"] = ["../reference_structures/ternary/POSCAR_mgb2_3", "vc-relax", "all", "vol", "nscf"]

    CalcD["stuffhex_z_1"] = ["../reference_structures/ternary/POSCAR_z_mgb2_1", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["stuffhex_z_2"] = ["../reference_structures/ternary/POSCAR_z_mgb2_2", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["stuffhex_z_3"] = ["../reference_structures/ternary/POSCAR_z_mgb2_3", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["rocksalt_2lay_abo2"] = ["../reference_structures/ternary/rocksalt.in.2lay_abo2", "vc-relax", "2Dxy", "2D", "nscf"]
    CalcD["caf2_abc"] = ["../reference_structures/ternary/POSCAR_caf2_ABC", "vc-relax", "all", "vol", "nscf"]

    CalcD["perov"] = ["../reference_structures/ternary/POSCAR_abo3", "vc-relax", "all", "vol", "nscf"]
    CalcD["perov2"] = ["../reference_structures/ternary/POSCAR_abo3_2", "vc-relax", "all", "vol", "nscf"]
    CalcD["perov3"] = ["../reference_structures/ternary/POSCAR_abo3_3", "vc-relax", "all", "vol", "nscf"]
    CalcD["perov4"] = ["../reference_structures/ternary/POSCAR_abo3_4", "vc-relax", "all", "vol", "nscf"]
    CalcD["perov5"] = ["../reference_structures/ternary/POSCAR_abo3_5", "vc-relax", "all", "vol", "nscf"]
    CalcD["perov6"] = ["../reference_structures/ternary/POSCAR_abo3_6", "vc-relax", "all", "vol", "nscf"]


    CalcD["simple_hex"] = ["../reference_structures/simple_hex.in", "vc-relax", "all", "flyaway", "nscf"]







    core_mono = [     "sc", "atom",     "sc_inv",     "bcc",     "bcc_inv",     "fcc",     "hcp", "hcp_shape",      "diamond",     "graphene",     "hex",     "square",     "dimer", "hex_2lay", "bcc_2lay", "fcc_dense", "bcc_dense", "znse_dense"]




    core_binary = [    "cscl",     "znse_shear",     "hbn",     "rocksalt",  "rocksalt_inv",   "znse",     "dimer2",     "square2", "hcp_v2", "rocksalt_shape", "rocksalt_2lay", "cscl_layers", "fcc_12", "fcc_21", "bcc_13", "bcc_31", "rocksalt_dense", "hcp_v2_dense", "znse_dense",  "distort", "227", "znseAABB", "trimer_ab2", "trimer2_ab2", "trimer_ba2", "trimer2_ba2"]

    A0 = [   "as_orth",    "ga_tet",    "ge_wurtz",    "pb_r3m",    "beta_sn",   "ga",    "bi",    "te",    "in",    "i2",    "li_p6mmm",    "hcp_shape",     "diamond_shear",   "n" ] #"bcc_tet.in",  same as POSCAR_ga_tet   #"as_221",  is simple cubic
    
    A1B1 = ["hbn_real", "sis_2d", "ges", "sns", "wurtz"] #nias #
    A1B2 = ["mgcl2", "mgcl2_tet", "caf2", "sis2", "tio2_rutile", "co2",  "mgf2",  "ticl2", "gei2"]   # "sns2" duplictes ticl2 #"mgf2_v2" is rutile
    A1B3 = ["alf3", "asna3_2d", "ab3_mol", "gain3", "li3n_hex"]
    A1B4 = ["sif4", "snf4", "gei4_mol"]
    A1B5 = ["ascl5"]
    A1B6 = ["mof6"]
    A2B3 = ["y2o3", "p2ca3", "al2o3", "bi2se3", "ga2s3", "gas"]
    A2B5 = ["nb2o5", "p2o5"]

    short_bonds = ["dimer_pair", "cao2_12", "irn2_38"]
    

#    metals = [ "mg2si", "simg2",  "mgb2_12", "mgb2_21", "fcc_conv_ABBB","fcc_conv_BAAA", "ab2_71", "ba2_71"]
#    metals = [ "mg2si", "simg2",  "mgb2_12", "mgb2_21",  "fcc_conv_ABBB","fcc_conv_BAAA",  "ab2_71", "ba2_71"]
    metals = [ "mg2si", "simg2",  "mgb2_12", "mgb2_21",    "ab2_71", "ba2_71"]

    all_ternary = ["abc_line", "bac_line", "cab_line", "fcc_tern", "hex_trim", "hh1", "hh2", "hh3", "stuffhex_1", "stuffhex_2", "stuffhex_3","stuffhex_z_1", "stuffhex_z_2", "stuffhex_z_3", "rocksalt_2lay_abo2", "caf2_abc", "perov", "perov2", "perov3",  "perov4",  "perov5",  "perov6"  ]

    pd = proto_data(CalcD, core_mono, core_binary, A0, A1B1, A1B2, A1B3, A1B4, A1B5, A1B6, A2B3, A2B5, metals, short_bonds, all_ternary)

    return pd

end

function  do_run(pd, T1, T2, T3, tmpname, dir, procs, torun; nscf_only = false, only_kspace = false, check_only = false)

    TORUN = []
    for t in torun
        if t == :core_mono
            TORUN = [TORUN;pd.core_mono]
        elseif t == :core_binary
            TORUN = [TORUN;pd.core_binary]
        elseif t == :A0
            TORUN = [TORUN;pd.A0]
        elseif t == :A1B1
            TORUN = [TORUN;pd.A1B1]
        elseif t == :A1B2
            TORUN = [TORUN;pd.A1B2]
        elseif t == :A1B3
            TORUN = [TORUN;pd.A1B3]
        elseif t == :A1B4
            TORUN = [TORUN;pd.A1B4]
        elseif t == :A1B5
            TORUN = [TORUN;pd.A1B5]
        elseif t == :A1B6
            TORUN = [TORUN;pd.A1B6]
        elseif t == :A2B3
            TORUN = [TORUN;pd.A2B3]
        elseif t == :A2B5
            TORUN = [TORUN;pd.A2B5]
        elseif t == :metals
            TORUN = [TORUN;pd.metals]
        elseif t == :short_bonds
            TORUN = [TORUN;pd.short_bonds]
        elseif t == :all_ternary
            TORUN = [TORUN;pd.all_ternary]
        else
            TORUN = [TORUN; t]
        end
    end

#    println("TORUN")
#    for t in TORUN
#        println(t)
#    end
#    println("--")
#    sleep(1)

    already_done = []
    not_done = []

    for st in TORUN

        file, scf, free, newst, calc_mode = pd.CalcD[st]

        if nscf_only 
            if calc_mode != "nscf"
                println("nscf_only $nscf_only skip $st")
                continue
            end
        end

        arr = split(file, '/')
        name = arr[end]
#        println("st $st $scf")

#############
        #check if already done

        if newst == "vol"
            ncalc = length([ 0.95 1.0 1.05])
        elseif newst == "vol-mid"
            ncalc = length( [ 0.9 0.95 1.0 1.05 1.1 ])
        elseif newst == "vol-dense"
            ncalc = length( [0.8 0.83 0.87 ])
        elseif newst == "vol-verydense"
            ncalc = length( [0.77 0.82 0.75 0.70 0.65 0.60 0.55 0.50])
        elseif newst == "vol-big"
            ncalc = length( [0.80 0.85 0.9 0.95 1.0 1.05 1.1 1.2 1.3 1.5 ])
        elseif newst == "vol-huge"
            ncalc = length( [0.9 0.95 1.0 1.05 1.1 1.2 1.3 1.5 2.0 2.5 3.0 3.5 4.0 5.0])
        elseif newst == "2D"
            ncalc = length( [0.90 0.95 1.0 1.05 1.10])
        elseif newst == "2D-mid"
            ncalc = length( [0.86 0.88 0.91 0.96 1.0 1.05 1.10])
        elseif newst == "2D-short"
            ncalc = length( [0.80 0.83])
        elseif newst == "shape"
            ncalc = length( [-0.06 -0.03 0.03 0.06])
        elseif newst == "coords"
            ncalc = length( [-0.20 -0.17 -0.14 -0.10 -0.07 -0.03 0.0 0.03 0.07 0.10 0.15 0.2 0.25 0.35 0.5])
        elseif newst == "coords-short"
            ncalc = length( [-0.3 -0.27 -0.25])
        elseif newst == "coords-small"
            ncalc = length( [ -0.15  -0.10 -0.05  0.0 0.05 0.10  0.15  ])
        elseif newst == "coords_super"
            ncalc = length( 0.08:.01:0.25)
        elseif newst == "break_inv"
            ncalc = length( [0.01 0.02 0.05 0.07 ])
        elseif newst == "shear"
            ncalc = length([0.01 0.02 ])
        elseif newst == "flyaway"
            ncalc = length( [-0.05, 0.0, 0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 3.0, 4.0])
        elseif newst == "scf"
            ncalc = 1
        elseif (newst == "coords_trimer"  ||  newst == "coords_trimer2" || newst == "coords_trimer_ab" || newst == "coords_trimer_ab_big" )
#            ncalc = length([1.05, 1.1, 1.15, 1.2, 1.25, 1.3])
            ncalc = length([1.05, 1.1,  1.2,  1.3])
        elseif (newst == "coords_trimer_dense" || newst == "coords_trimer_ab_dense" )
            ncalc = length([1.0, 0.95])
        else
            println("error newst: ", newst)
            ncalc = 1
        end

        for n in 1:ncalc
            d="$dir/$name"*"_vnscf_"*"$newst"*"_"*"$n"        
            if isfile(d*"/projham_K.xml") ||  isfile(d*"/projham_K.xml.gz")
                push!(already_done, d)
            else
                push!(not_done, d)
            end
        end
        if check_only
            continue
        end

        d="$dir/$name"*"_vnscf_"*"$newst"*"_"*"$ncalc"        
        
        if isfile(d*"/projham_K.xml") ||  isfile(d*"/projham_K.xml.gz")
            println("everything is already done! $d")
            println("we can move on")
            continue
        else
            println("not done yet")
        end
#######################


        randi = Int64(round(rand()*1000000))
        try
            println("read crys")
            c = CrystalMod.makecrys(file)
            
            for i in 1:c.nat
                if c.types[i] == "A"
                    c.types[i] = T1
                elseif c.types[i] == "B"
                    c.types[i] = T2
                elseif c.types[i] == "C"
                    c.types[i] = T3
                else
                    c.types[i] = T2
                end
            end

            println(c)
            println("did read crys, run dft")

            if scf != "none"
            

                #preadjust vol
                avg_rad = 0.0
                for t in c.types
                    avg_rad += Atomdata.atom_radius[t] / 100.0 / 0.529177
                end
                avg_rad = avg_rad / c.nat
                
                vol_peratom = abs(det(c.A)) / c.nat
                
                ratio = (avg_rad^3 * (4 * pi / 3)) / vol_peratom
                println("ratio $ratio")
                if ratio > 1.1
                    c.A = c.A * (ratio^(1.0/3.0))
                    println("new starting c.A")
                    println(c.A)
                end
                

                println("START DFT.runSCF")
                
                dft_ref = DFT.runSCF(c, inputstr=name, nprocs=procs, prefix="$name.qe.relax", directory="$dir", tmpdir="/$tmpname/$name.$randi", wannier=false, code="QE", skip=false, calculation=scf, dofree=free, cleanup=true)

                println("did dft, get new struct")

                if (scf == "relax" || scf == "vc-relax") && maximum(abs.(dft_ref.stress)) > 2e-5
                    println()
                    println("structure convergence not reached, relax again")
                    c2 = dft_ref.crys
                    println(c2)
                    println()

                    dft_ref = DFT.runSCF(c2, inputstr=name, nprocs=procs, prefix="$name.qe.relax", directory="$dir", tmpdir="/$tmpname/$name.$randi", wannier=false, code="QE", skip=false, calculation=scf, dofree=free, cleanup=true)
                end

                println("END DFT.runSCF")
                println(dft_ref)

                cnew = deepcopy(dft_ref.crys)
                #            cnew = deepcopy(c)

            else
                cnew = deepcopy(c)
            end

            torun = []
            if newst == "vol"
                for x in [ 0.95 1.0 1.05]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "vol-mid"
                for x in [ 0.9 0.95 1.0 1.05 1.1 ]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "vol-dense"
                for x in [0.8 0.83 0.87 ]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "vol-verydense"
                for x in [0.77 0.82 0.75 0.70 0.65 0.60 0.55 0.50]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "vol-big"
                for x in [0.80 0.85 0.9 0.95 1.0 1.05 1.1 1.2 1.3 1.5 ]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "vol-huge"
                for x in [0.9 0.95 1.0 1.05 1.1 1.2 1.3 1.5 2.0 2.5 3.0 3.5 4.0 5.0]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "2D"
                for x in [0.90 0.95 1.0 1.05 1.10]
                    c = deepcopy(cnew)
                    c.A[1:2,:] = c.A[1:2,:] * x
                    push!(torun, deepcopy(c))
                end 
            elseif newst == "2D-mid"
                for x in [0.86 0.88 0.91 0.96 1.0 1.05 1.10]
                    c = deepcopy(cnew)
                    c.A[1:2,:] = c.A[1:2,:] * x
                    push!(torun, deepcopy(c))
                end 
            elseif newst == "2D-short"
                for x in [0.80 0.83]
                    c = deepcopy(cnew)
                    c.A[1:2,:] = c.A[1:2,:] * x
                    push!(torun, deepcopy(c))
                end 
            elseif newst == "shape"
                for x in [-0.06 -0.03 0.03 0.06]
                    c = deepcopy(cnew)
                    c.A[1:2,:] = c.A[1:2,:] * (1+x)
                    c.A[3,:] = c.A[3,:] * (1- 2.0 * x)
                    push!(torun, deepcopy(c))
                end 
            elseif newst == "coords_trimer_ab"

                a = min_dimer_dist_dict[ cnew.types[1]]
                println("coords_trimer_ab $a ", cnew.types[1], " " , cnew.types[3] )
                ab = 1.0
                try
                    name_t="dimer.in"
                    ncalc_t = 1
                    newst_t = "coords"
                    ab_dir="$dir/$name_t"*"_vnscf_"*"$newst_t"*"_"*"$ncalc_t"        
                    println("try $ab_dir")
                    dft = QE.loadXML(ab_dir*"/qe.save")

                    ab = -dft.crys.coords[1,3] * dft.crys.A[3,3] * 2.0
                    println("ab $ab loaded")

                catch
                    ab = (min_dimer_dist_dict[T1] + min_dimer_dist_dict[T2]) / 2.0 
                    println("ab $ab estimated")

                end
                
#                for x in [1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
                for x in [1.05, 1.1,  1.2,  1.3]

                    c = deepcopy(cnew)
                    if c.coords[1,3] == 0.8
                        c.A[1,1] = ab / (0.4^2 + 0.2^2)^0.5 * x
                    else
                        c.A[1,1] = ab / 0.4 * x
                    end

                    c.A[2,2] = ab * 2.0 * x
                    c.A[3,3] = a / 0.4 * x
                    push!(torun, deepcopy(c))
                end


            elseif newst == "coords_trimer_ab_dense"

                a = min_dimer_dist_dict[ cnew.types[1]]
                println("coords_trimer_ab $a ", cnew.types[1], " " , cnew.types[3] )
                ab = 1.0
                try
                    name_t="dimer.in"
                    ncalc_t = 1
                    newst_t = "coords"
                    ab_dir="$dir/$name_t"*"_vnscf_"*"$newst_t"*"_"*"$ncalc_t"        
                    println("try $ab_dir")
                    dft = QE.loadXML(ab_dir*"/qe.save")

                    ab = -dft.crys.coords[1,3] * dft.crys.A[3,3] * 2.0
                    println("ab $ab loaded")

                catch
                    ab = (min_dimer_dist_dict[T1] + min_dimer_dist_dict[T2]) / 2.0 
                    println("ab $ab estimated")

                end
                
#                for x in [1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
                for x in [1.0, 0.95]

                    c = deepcopy(cnew)
                    if c.coords[1,3] == 0.8
                        c.A[1,1] = ab / (0.4^2 + 0.2^2)^0.5 * x
                    else
                        c.A[1,1] = ab / 0.4 * x
                    end

                    c.A[2,2] = ab * 2.0 * x
                    c.A[3,3] = a / 0.4 * x
                    push!(torun, deepcopy(c))
                end


            elseif newst == "coords_trimer_ab_big"

                a = min_dimer_dist_dict[ cnew.types[1]]
                println("coords_trimer_ab big $a ", cnew.types[1], " " , cnew.types[3] )
                ab = 1.0
                try
                    name_t="dimer.in"
                    ncalc_t = 1
                    newst_t = "coords"
                    ab_dir="$dir/$name_t"*"_vnscf_"*"$newst_t"*"_"*"$ncalc_t"        
                    println("try $ab_dir")
                    dft = QE.loadXML(ab_dir*"/qe.save")

                    ab = -dft.crys.coords[1,3] * dft.crys.A[3,3] * 2.0
                    println("ab $ab loaded")

                catch
                    ab = (min_dimer_dist_dict[T1] + min_dimer_dist_dict[T2]) / 2.0 
                    println("ab $ab estimated")

                end
                
                for x in [1.05, 1.15, 1.25, 1.35, 1.45, 1.55]
#                for x in [1.05, 1.1,  1.2,  1.3]

                    c = deepcopy(cnew)
                    c.A[1,1] = ab / (0.4) * x

                    c.A[2,2] = ab * 3.0 * x
                    c.A[3,3] = a / 0.2 * x
                    push!(torun, deepcopy(c))
                end


            elseif newst == "coords_trimer"
                a = min_dimer_dist_dict[T1]
                for x in [1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
                    c = deepcopy(cnew)
                    c.A[1,1] = a / (0.4^2 + 0.2^2)^0.5 * x
                    c.A[2,2] = a * 1.5 * x
                    c.A[3,3] = a / 0.4 * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "coords_trimer_dense"
                a = min_dimer_dist_dict[T1]
                for x in [1.0, 0.95]
                    c = deepcopy(cnew)
                    c.A[1,1] = a / (0.4^2 + 0.2^2)^0.5 * x
                    c.A[2,2] = a * 1.5 * x
                    c.A[3,3] = a / 0.4 * x
                    push!(torun, deepcopy(c))
                end

            elseif newst == "coords_trimer2"
                a = min_dimer_dist_dict[T1]
                for x in [1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
                    c = deepcopy(cnew)
                    c.A[1,1] = a / 0.4 * x
                    c.A[2,2] = a * 1.5 * x
                    c.A[3,3] = a / 0.4 * x
                    push!(torun, deepcopy(c))
                end

            elseif newst == "coords"
                for x in [-0.20 -0.17 -0.14 -0.10 -0.07 -0.03 0.0 0.03 0.07 0.10 0.15 0.2 0.25 0.35 0.5]
                    c = deepcopy(cnew)
                    c.coords = c.coords * (1+x)
                    push!(torun, deepcopy(c))
                end
            elseif newst == "coords-short"
                for x in [-0.3 -0.27 -0.25]
                    c = deepcopy(cnew)
                    c.coords = c.coords * (1+x)
                    push!(torun, deepcopy(c))
                end
            elseif newst == "coords-small"
                for x in [ -0.15  -0.10 -0.05  0.0 0.05 0.10  0.15  ]
                    c = deepcopy(cnew)
                    c.coords = c.coords * (1+x)
                    push!(torun, deepcopy(c))
                end
            elseif newst == "coords_super"
                for x in 0.08:.01:0.25
                    c = deepcopy(cnew)
                    c.coords[1,3] = -1.0*x
                    c.coords[2,3] = x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "break_inv"
                for x in [0.01 0.02 0.05 0.07 ]
                    c = deepcopy(cnew)
		    if c.nat == 1
		        c = c * [2 1 1]
                    end
                    c.A = c.A * 1.02
                    c.coords[1,1] = c.coords[1,1] + x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "shear"
                for x in [0.01 0.02 ]
                    c = deepcopy(cnew)
                    c = c 
                    c.A = c.A * (I +  [0 x 0; x 0 0; 0 0 0])

                    push!(torun, deepcopy(c))
                end 

            elseif newst == "flyaway"
                for x in [-0.05, 0.0, 0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 3.0, 4.0]
                    c = deepcopy(cnew)
                    c = c 
                    c.A = c.A * (I +  [0 0 0; 0 0 0; 0 0 x])

                    push!(torun, deepcopy(c))
                end 
                
            elseif newst == "scf"
                push!(torun, deepcopy(cnew))
            else
                println("error newst: ", newst)
            end
            

            for (i,c) in enumerate(torun)
                println("start torun $name  $i")
                try

                    d="$dir/$name"*"_vnscf_"*"$newst"*"_"*"$i"        
                    println(d)
                    println(c)
                    if isfile(d*"/projham_K.xml") ||  isfile(d*"/projham_K.xml.gz")
                        println("continue")
                        continue
                    end

                    dft = DFT.runSCF(c, nprocs=procs, prefix="qe", directory="$d", tmpdir="$d", wannier=false, code="QE", skip=true, cleanup=true)
                    if calc_mode == "nscf"

                        try
                            tbc, tbck = AtomicProj.projwfx_workf(dft, nprocs=procs, directory=d, skip_og=true, skip_proj=true, freeze=true, localized_factor = 0.15, cleanup=true, only_kspace=only_kspace)
                        catch err3
                            println("err3")
                            println(err3)
                            println("skip nscf $i ")##
                        end
                    end

                catch err2
                    println("err2")
                    println(err2)
                    println("skip $i dft $i ")##
                   
                end
            end
        catch err
            println(err)
            println("skip everything  $st ")
        end
    end
    
    if check_only
        return already_done, not_done
    end

end

pd = setup_proto_data()

function oxidation_guess(atom1, atom2)

    if atom1 == atom2
        keep = [[atom1, atom1, :core_mono]]
        keep = push!(keep , [atom1, atom1, :A0])
        return keep
    end


    possible_configs = zeros(Int64, 0, 5)

    for (c1, o1) = enumerate(atom_prefered_oxidation[atom1])
        for (c2,o2) = enumerate(atom_prefered_oxidation[atom2])

            if c1 == 1
                score1=1
            elseif o1 == 0
                score1=4
            else
                score1=2
            end

            if c2 == 1
                score2=1
            elseif o2 == 0
                score2=4
            else
                score2=2
            end


            for n1 = 1:6
                for n2 = 1:6
                    
                    if gcd(n1,n2) > 1
                        continue
                    end

                    if o1 * n1 + o2 * n2 == 0
                        possible_configs = [possible_configs; n1 n2 score1+score2   score1  score2]
                    elseif  abs(o1 * n1 + o2 * n2) == 1
                        possible_configs = [possible_configs; n1 n2 score1+score2+10  score1  score2] 
                    end
                end
            end
        end
    end

    possible_configs = possible_configs[sortperm(possible_configs[:, 3]), :]

#    println("possible configs $atom1 $atom2 score score1 score2")
#    for p in 1:min(size(possible_configs,1), 6)
#        println(possible_configs[p,:])
#    end

    n_config = 0
    keep = []
    
    keep = [[atom1, atom2, :core_binary]]

    for p in 1:size(possible_configs,1)
        use = false
        if possible_configs[p,3] == 2 || possible_configs[p,3] == 3
#            println("a ",  possible_configs[p,:])
            use = true
            n_config += 1
        elseif possible_configs[p,3] == 4 && n_config <= 2
            use = true
#            println("b ",  possible_configs[p,:])
            n_config += 1
#        else
#            println("R ",  possible_configs[p,:])

        end
        if use
            if possible_configs[p,1] == 1 && possible_configs[p,2] == 1
                push!(keep, [atom1, atom2, :A1B1])
            elseif possible_configs[p,1] == 1 && possible_configs[p,2] == 2
                push!(keep, [atom1, atom2, :A1B2])
            elseif possible_configs[p,1] == 1 && possible_configs[p,2] ==3
                push!(keep, [atom1, atom2, :A1B3])
            elseif possible_configs[p,1] == 1 && possible_configs[p,2] == 4
                push!(keep, [atom1, atom2, :A1B4])
            elseif possible_configs[p,1] == 1 && possible_configs[p,2] == 5
                push!(keep, [atom1, atom2, :A1B5])
            elseif possible_configs[p,1] == 1 && possible_configs[p,2] == 6
                push!(keep, [atom1, atom2, :A1B6])
            elseif possible_configs[p,1] == 2 && possible_configs[p,2] == 3
                push!(keep, [atom1, atom2, :A2B3])
            elseif possible_configs[p,1] == 2 && possible_configs[p,2] == 5
                push!(keep, [atom1, atom2, :A2B5])


            elseif possible_configs[p,2] == 1 && possible_configs[p,1] == 2
                push!(keep, [atom2, atom1, :A1B2])
            elseif possible_configs[p,2] == 1 && possible_configs[p,1] == 3
                push!(keep, [atom2, atom1, :A1B3])
            elseif possible_configs[p,2] == 1 && possible_configs[p,1] == 4
                push!(keep, [atom2, atom1, :A1B4])
            elseif possible_configs[p,2] == 1 && possible_configs[p,1] == 5
                push!(keep, [atom2, atom1, :A1B5])
            elseif possible_configs[p,2] == 1 && possible_configs[p,1] == 6
                push!(keep, [atom2, atom1, :A1B6])
            elseif possible_configs[p,2] == 2 && possible_configs[p,1] == 3
                push!(keep, [atom2, atom1, :A2B3])
            elseif possible_configs[p,2] == 2 && possible_configs[p,1] == 5
                push!(keep, [atom2, atom1, :A2B5])
            end


        end
    end
    if length(keep) == 1
        if maximum(atom_prefered_oxidation[atom1]) > 0 || maximum(atom_prefered_oxidation[atom2]) > 0
            push!(keep, [atom1, atom2, :metals])
        end
    end

#    for k in keep
#        println(k)
#    end
    return keep 
end

#    A1B1 = ["hbn_real", "sis_2d", "ges", "sns", "wurtz"]
#    A1B2 = ["mgcl2", "mgcl2_tet", "caf2", "sis2", "tio2_rutile", "co2", "square_ab2", "mgf2",  "ticl2", "gei2"]   # "sns2" duplictes ticl2 #"mgf2_v2" is rutile
#    A1B3 = ["alf3", "asna3_2d", "ab3_mol"]
#    A1B4 = ["sif4", "snf4", "gei4_mol"]
#    A1B5 = ["ascl5"]
#    A2B3 = ["y2o3", "p2ca3", "al2o3", "bi2se3", "ga2s3", "gas"]
#    A2B5 = ["nb2o5"]