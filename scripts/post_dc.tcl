set VERILOG_OUTPUT_FILE "${DESIGN_NAME}_netlist.v"
set DDC_OUTPUT_FILE "${DESIGN_NAME}_post.ddc"

write -format verilog -hierarchy -output ${OUTPUT_DIR}/${DESIGN_NAME}_netlist.v
write -format ddc     -hierarchy -output ${OUTPUT_DIR}/${DESIGN_NAME}_post.ddc
write_sdf ${OUTPUT_DIR}/${DESIGN_NAME}.sdf
write_sdc ${OUTPUT_DIR}/${DESIGN_NAME}.sdc
set_svf -off
rt > ${OUTPUT_DIR}/report_timing.rpt
