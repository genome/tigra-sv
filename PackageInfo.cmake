cmake_minimum_required(VERSION 2.8)

# .deb packaging
set(ARCH "i686")
if(${CMAKE_SIZEOF_VOID_P} MATCHES 8)
    set(ARCH "x86_64")
endif ()

# The format of the description field is a short summary line followed by a
# longer paragraph indented by a single space on each line
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
"TIGRA_SV is a program that conducts targeted local assembly of structural variants.
 TIGRA_SV is a program that conducts targeted local assembly of structural
 variants (SV) using the iterative graph routing assembly (TIGRA) algorithm (L.
 Chen et al, unpublished). It takes as input a list of putative SV calls and a
 set of bam files that contain reads mapped to a reference genome such as NCBI
 build36.  For each SV call, it assembles the set of reads that were mapped or
 partially mapped to the region of interest (ROI) in the corresponding bam
 files. Instead of outputing a single consensus sequence, TIGRA_SV attempts to
 construct all the alternative alleles in the ROI as long as they received
 sufficient sequence coverage (usually >= 2x).  It also utilizes the variant
 type information in the input files to select reads for assembly.  TIGRA_SV is
 effective at improving the SV prediction accuracy and resolution in short reads
 analysis and can produce accurate breakpoint sequences that are useful to
 understand the origin, mechanism and pathology underlying the SVs.")

set(CPACK_PACKAGE_NAME "tigra-sv${EXE_VERSION_SUFFIX}")
set(CPACK_PACKAGE_VENDOR "wugc")
set(CPACK_PACKAGE_VERSION ${FULL_VERSION})
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Ken Chen <kchenken@gmail.com>")
set(CPACK_SYSTEM_NAME "Linux-${ARCH}")
set(CPACK_TOPLEVEL_TAG "Linux-${ARCH}")
set(CPACK_DEBIAN_PACKAGE_SECTION science)
set(CPACK_DEBIAN_PACKAGE_PRIORITY optional)
set(CPACK_DEBIAN_PACKAGE_REPLACES "")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.4), libgcc1 (>= 1:4.1.1-21), libstdc++6 (>= 4.2.1-4), zlib1g (>= 1:1.1.4)")
if (CMAKE_BUILD_TYPE MATCHES package)
    set(CPACK_GENERATOR "DEB")
else(CMAKE_BUILD_TYPE MATCHES package)
    set(CPACK_GENERATOR "TGZ")
endif(CMAKE_BUILD_TYPE MATCHES package)

configure_file(debian/postinst.in debian/postinst @ONLY)
configure_file(debian/prerm.in debian/prerm @ONLY)
set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "debian/postinst;debian/prerm")

include(CPack)
