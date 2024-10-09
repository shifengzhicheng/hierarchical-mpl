# Case in which we treat the root as a macro cluster and jump
# from initializing the physical tree to macro placement
source "helpers.tcl"

read_lef "./BP/lef/NangateOpenCellLibrary.tech.lef"
read_lef "./BP/lef/NangateOpenCellLibrary.macro.rect.lef" 
read_lef "./BP/lef/fakeram45_32x32.lef"
read_lef "./BP/lef/fakeram45_64x62.lef"
read_lef "./BP/lef/fakeram45_64x124.lef"
read_lef "./BP/lef/fakeram45_128x116.lef"
read_lef "./BP/lef/fakeram45_256x48.lef"
read_lef "./BP/lef/fakeram45_512x64.lef"


read_liberty "./BP/lib/NangateOpenCellLibrary_typical.lib"
read_liberty "./BP/lib/fakeram45_32x32.lib"
read_liberty "./BP/lib/fakeram45_64x62.lib"
read_liberty "./BP/lib/fakeram45_64x124.lib"
read_liberty "./BP/lib/fakeram45_128x116.lib"
read_liberty "./BP/lib/fakeram45_256x48.lib"
read_liberty "./BP/lib/fakeram45_512x64.lib"

put "reading verilog"
read_verilog "./BP/1.v"
put "done reading verilog"
put "reading def"
read_def "./BP/1.def"
put "done reading def"
if 0 {
    
    set_thread_count 12
    rtl_macro_placer -report_directory results/1 -halo_width 10

    set def_file [make_result_file 1.def]
    write_def $def_file

    diff_files 1.defok $def_file
}
