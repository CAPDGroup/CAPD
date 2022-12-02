#m4_define([CAPD_VERSION_NUMBER], m4_esyscmd(tr -d '\n' < capdMake_PATH[/capd_version_number] ))
#m4_define([CAPD_LIBRARY_VERSION_NUMBER], m4_esyscmd(tr -d '\n' < capdMake_PATH[/capd_library_version_number]))
m4_include(capdMake_PATH[/capd_library_version_number.m4])
m4_include(capdMake_PATH[/capd_version_number.m4])
