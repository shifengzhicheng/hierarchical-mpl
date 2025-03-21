source load_design.tcl

rtl_macro_placer -halo_width 5 -halo_height 5 -target_util 0.7

write_def ${case}_fp_HierRTLMP.def