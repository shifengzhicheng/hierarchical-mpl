# include all rtl files here
set rtl_dir "./rtl"
set rtl_list {
    ./rtl/fakeram45_32x32_dp.v
    ./rtl/bsg_chip_block.sv2v.v
}

set HYPER_DEFINE {} 

proc getAllSubdirs {path} {
    set result {}
    foreach item [glob -nocomplain -directory $path *] {
        if {[file isdirectory $item]} {
            set norm [file normalize $item]
            lappend result $norm
            set subs [getAllSubdirs $item]
            foreach subdir $subs {
                lappend result $subdir
            }
        }
    }
    return $result
}

if {![info exists search_path]} {
    set search_path {}
}

set additionalPaths [getAllSubdirs $rtl_dir]

foreach p $additionalPaths {
    if {[lsearch -exact $search_path $p] == -1} {
        lappend search_path $p
    }
}