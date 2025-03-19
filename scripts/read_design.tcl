# Setting lef files
# check if pwd have case1.odb
set tech_dir "/OpenROAD-flow-scripts/flow/platforms/nangate45"
set tech_lef "${tech_dir}/lef/NangateOpenCellLibrary.tech.lef"
set std_lef "${tech_dir}/lef/NangateOpenCellLibrary.macro.mod.lef"

set lefs "
    ${tech_dir}/lef/NangateOpenCellLibrary.macro.lef \
    ${tech_dir}/lef/NangateOpenCellLibrary.macro.rect.lef \
    ./fakeram45_256x16.lef
"
# Setting lib files
set libs "
    ${tech_dir}/lib/NangateOpenCellLibrary_typical.lib \
    ./fakeram45_256x16.lib \
"

read_lef $tech_lef
read_lef $std_lef

foreach lef_file ${lefs} {
  read_lef $lef_file
}

foreach lib_file ${libs} {
  read_liberty $lib_file
}

puts "read ariane_CMP.def"
read_verilog ../ariane.v
link_design ariane
read_sdc ../ariane.sdc
read_def -floorplan_initialize ./ariane_fp_HierRTLMP.def
# Save the design
write_db ariane_CMP.odb

# Liberty units are fF,kOhm
set_layer_rc -layer metal1 -resistance 5.4286e-03 -capacitance 7.41819E-02
set_layer_rc -layer metal2 -resistance 3.5714e-03 -capacitance 6.74606E-02
set_layer_rc -layer metal3 -resistance 3.5714e-03 -capacitance 8.88758E-02
set_layer_rc -layer metal4 -resistance 1.5000e-03 -capacitance 1.07121E-01
set_layer_rc -layer metal5 -resistance 1.5000e-03 -capacitance 1.08964E-01
set_layer_rc -layer metal6 -resistance 1.5000e-03 -capacitance 1.02044E-01
set_layer_rc -layer metal7 -resistance 1.8750e-04 -capacitance 1.10436E-01
set_layer_rc -layer metal8 -resistance 1.8750e-04 -capacitance 9.69714E-02
# No calibration data available for metal9 and metal10
#set_layer_rc -layer metal9 -resistance 3.7500e-05 -capacitance 3.6864e-02
#set_layer_rc -layer metal10 -resistance 3.7500e-05 -capacitance 2.8042e-02

set_wire_rc -signal -layer metal3
set_wire_rc -clock  -layer metal5

