source load_design.tcl

rtl_macro_placer -halo_width 10 -halo_height 10 -target_util 0.7

write_def ${case}_fp_HierRTLMP.def