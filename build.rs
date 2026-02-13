fn main() {
    // Compile libbsc C++ library together with libsais in one build
    // This ensures proper linking of dependencies
    cc::Build::new()
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
        .flag_if_supported("-std=c++11")
        .flag_if_supported("-fopenmp") // For multithreading support
        .include("third_party/libbsc/libbsc")
        .include("third_party/libbsc/libbsc/bwt/libsais") // For libsais headers
        .warnings(false) // Suppress warnings from third-party code
        .compile("libbsc");

    println!("cargo:rerun-if-changed=third_party/libbsc");

    // Link OpenMP for multithreading (if available)
    println!("cargo:rustc-link-lib=gomp"); // GNU OpenMP
}
