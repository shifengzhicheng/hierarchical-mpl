source ./.synopsys_dc.setup
source ./pre_dc.tcl
compile_ultra -retime
source ./post_dc.tcl

report_area
report_power
report_timing