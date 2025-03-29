# include all rtl files here

set_app_var target_library  [glob ../libs/*.db]
set_app_var link_library $target_library

set HYPER_DEFINE {}

set sdc_file ariane.sdc

set rtl_list {
./rtl/ariane.v
}

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

set additionalPaths [getAllSubdirs .]

foreach p $additionalPaths {
    if {[lsearch -exact $search_path $p] == -1} {
        lappend search_path $p
    }
}
