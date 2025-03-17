set DESIGN_NAME        "bsg_chip"
set CONSTRAIN_FILE     "./bsg_chip.sdc"
set DES_LIB  "./work/"
set OUTPUT_DIR "./DC/"


if {![file exists $DES_LIB]} {
    echo "Generate dir design lib"
    sh mkdir $DES_LIB
}
if {![file exists ${OUTPUT_DIR}]} {
    echo "Generate output directory"
    sh mkdir ${OUTPUT_DIR}
}
define_design_lib MY_LIB -path $DES_LIB

# set_svf ${OUTPUT_DIR}${DESIGN_NAME}.svf
# saif_map -start
# analyze -format verilog ${RTL_SOURCE_FILES} -library MY_LIB
# elaborate ${DESIGN_NAME} -library MY_LIB
# write -hierarchy -format ddc -output ./${OUTPUT_DIR}${DESIGN_NAME}_pre.ddc
# source -echo -verbose ${CONSTRAIN_FILE}

set_svf ${OUTPUT_DIR}${DESIGN_NAME}.svf
saif_map -start
foreach rtl_file $rtl_list {
    analyze -format sverilog -define $HYPER_DEFINE $rtl_file -library MY_LIB
}
elaborate ${DESIGN_NAME} -library MY_LIB
write -hierarchy -format ddc -output ${OUTPUT_DIR}${DESIGN_NAME}_pre.ddc
source -echo -verbose ${CONSTRAIN_FILE}
