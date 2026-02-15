fn main() {
    // C/C++ PGO support:
    //   QZ_PGO_GENERATE=/path/to/dir  → instrument for profile collection
    //   QZ_PGO_USE=/path/to/dir       → optimize using collected profiles
    let pgo_generate = std::env::var("QZ_PGO_GENERATE").ok();
    let pgo_use = std::env::var("QZ_PGO_USE").ok();

    // Compile libbsc C++ library together with libsais in one build
    // This ensures proper linking of dependencies
    let mut bsc_build = cc::Build::new();
    bsc_build
        .cpp(true)
        .files(&[
            // libsais (C source - dependency of libbsc)
            "third_party/libbsc/libbsc/bwt/libsais/libsais.c",
            // libbsc C++ sources
            "third_party/libbsc/libbsc/libbsc/libbsc.cpp",
            "third_party/libbsc/libbsc/adler32/adler32.cpp",
            "third_party/libbsc/libbsc/bwt/bwt.cpp",
            "third_party/libbsc/libbsc/coder/coder.cpp",
            "third_party/libbsc/libbsc/coder/qlfc/qlfc.cpp",
            "third_party/libbsc/libbsc/coder/qlfc/qlfc_model.cpp",
            "third_party/libbsc/libbsc/filters/detectors.cpp",
            "third_party/libbsc/libbsc/filters/preprocessing.cpp",
            "third_party/libbsc/libbsc/lzp/lzp.cpp",
            "third_party/libbsc/libbsc/platform/platform.cpp",
            "third_party/libbsc/libbsc/st/st.cpp",
        ])
        .flag_if_supported("-O3")
        .flag_if_supported("-march=native")     // Use AVX2/SSE for faster BWT
        .flag_if_supported("-std=c++11")
        .flag_if_supported("-fopenmp") // For multithreading support
        .define("LIBBSC_OPENMP_SUPPORT", None) // Enable BSC's OpenMP code paths
        .define("LIBSAIS_OPENMP", None)        // Enable libsais parallel BWT
        .include("third_party/libbsc/libbsc")
        .include("third_party/libbsc/libbsc/bwt/libsais") // For libsais headers
        .warnings(false);

    if let Some(ref dir) = pgo_generate {
        bsc_build.flag(&format!("-fprofile-generate={}", dir));
    }
    if let Some(ref dir) = pgo_use {
        bsc_build.flag(&format!("-fprofile-use={}", dir));
        bsc_build.flag("-fprofile-correction"); // Handle slightly mismatched counts
    }
    bsc_build.compile("libbsc");

    println!("cargo:rerun-if-changed=third_party/libbsc");

    // Compile htscodecs (fqzcomp_qual + utils) as static C lib
    let mut hts_build = cc::Build::new();
    hts_build
        .files(&[
            "third_party/htscodecs/htscodecs/fqzcomp_qual.c",
            "third_party/htscodecs/htscodecs/utils.c",
        ])
        .flag_if_supported("-O3")
        .include("third_party/htscodecs/htscodecs")
        .warnings(false);

    if let Some(ref dir) = pgo_generate {
        hts_build.flag(&format!("-fprofile-generate={}", dir));
    }
    if let Some(ref dir) = pgo_use {
        hts_build.flag(&format!("-fprofile-use={}", dir));
        hts_build.flag("-fprofile-correction");
    }
    hts_build.compile("htscodecs");

    println!("cargo:rerun-if-changed=third_party/htscodecs/htscodecs");

    // Link OpenMP for multithreading
    // Pass both -L and -lgomp as link-args so they appear together
    // after -Bdynamic in the link command (rust-lld needs this)
    for entry in std::fs::read_dir("/usr/lib/gcc/x86_64-linux-gnu").into_iter().flatten() {
        if let Ok(e) = entry {
            if e.path().join("libgomp.so").exists() {
                println!("cargo:rustc-link-arg=-L{}", e.path().display());
            }
        }
    }
    println!("cargo:rustc-link-arg=-lgomp");
    println!("cargo:rustc-link-arg=-lstdc++");
    // Link math lib for htscodecs (log, etc.)
    println!("cargo:rustc-link-lib=m");

    // Link gcov runtime for C/C++ PGO instrumentation
    if pgo_generate.is_some() {
        // Force whole-archive inclusion since gcov symbols are referenced
        // from instrumented object files that are already linked
        println!("cargo:rustc-link-arg=-Wl,--whole-archive");
        println!("cargo:rustc-link-arg=/usr/lib/gcc/x86_64-linux-gnu/13/libgcov.a");
        println!("cargo:rustc-link-arg=-Wl,--no-whole-archive");
    }
}
