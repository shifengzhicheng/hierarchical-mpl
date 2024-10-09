# Case in which we treat the root as a macro cluster and jump
# from initializing the physical tree to macro placement
source "helpers.tcl"

read_lef ./SW/lef/asap7_tech_1x_201209.lef
set lef_files [glob "./SW/lef/*.lef"]
foreach lef_file $lef_files {
    put "reading $lef_file"
	read_lef $lef_file
}

set lib_files [glob "./SW/lib/*.lib"]
foreach lib_file $lib_files {
    put "reading $lib_file"
	read_liberty $lib_file
}

put "reading verilog"
read_verilog "./SW/sw.v"
put "done reading verilog"
put "reading def"
read_def "./SW/sw.def"
put "done reading def"
