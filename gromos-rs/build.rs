use std::env;
use std::process::Command;

#[cfg(feature = "use-cuda")]
use std::path::PathBuf;

fn main() {
    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();

    // Generate C header file using cbindgen (with config file)
    // Note: Binary I/O modules are Rust-only and excluded from C bindings
    // Skip cbindgen on stable Rust to avoid -Zunpretty errors
    if env::var("RUSTUP_TOOLCHAIN")
        .map(|t| t.contains("nightly"))
        .unwrap_or(false)
    {
        let config = cbindgen::Config::from_file("cbindgen.toml").unwrap_or_default();
        match cbindgen::Builder::new()
            .with_crate(&crate_dir)
            .with_config(config)
            .generate()
        {
            Ok(bindings) => {
                bindings.write_to_file("gromos_rs.h");
                println!("cargo:warning=C bindings generated successfully");
            },
            Err(e) => {
                println!("cargo:warning=Unable to generate C bindings: {:?}", e);
                println!("cargo:warning=Binary I/O modules are Rust-only and excluded from C API");
            },
        }
    } else {
        println!("cargo:warning=Skipping C bindings generation (requires nightly Rust)");
        println!("cargo:warning=Binary I/O modules are Rust-only and excluded from C API");
    }

    // Auto-detect FFTW3 library (but don't auto-enable - user must opt-in via features)
    if is_fftw3_available() {
        println!("cargo:warning=FFTW3 detected - enable with --features use-fftw for high-performance FFT");
    } else {
        println!("cargo:warning=FFTW3 not found - using rustfft fallback");
    }

    // Compile CUDA kernels if feature enabled
    #[cfg(feature = "use-cuda")]
    compile_cuda_kernels();

    println!("cargo:rerun-if-changed=src/");
    println!("cargo:rerun-if-changed=src/gpu_kernels.cu");
}

/// Check if FFTW3 library is available on the system
fn is_fftw3_available() -> bool {
    // Try pkg-config first (Linux, macOS with pkg-config)
    if Command::new("pkg-config")
        .args(["--exists", "fftw3"])
        .output()
        .map(|output| output.status.success())
        .unwrap_or(false)
    {
        return true;
    }

    // Try to find library manually (for systems without pkg-config)
    // Check common library paths
    let lib_paths = [
        "/usr/lib",
        "/usr/local/lib",
        "/opt/homebrew/lib",         // macOS ARM
        "/usr/lib/x86_64-linux-gnu", // Debian/Ubuntu
    ];

    for path in &lib_paths {
        let lib_file_so = format!("{}/libfftw3.so", path);
        let lib_file_dylib = format!("{}/libfftw3.dylib", path);
        let lib_file_a = format!("{}/libfftw3.a", path);

        if std::path::Path::new(&lib_file_so).exists()
            || std::path::Path::new(&lib_file_dylib).exists()
            || std::path::Path::new(&lib_file_a).exists()
        {
            return true;
        }
    }

    false
}

/// Compile CUDA kernels to PTX
#[cfg(feature = "use-cuda")]
fn compile_cuda_kernels() {
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let cuda_file = "src/gpu_kernels.cu";
    let ptx_file = out_dir.join("gpu_kernels.ptx");

    // Check if CUDA file exists
    if !std::path::Path::new(cuda_file).exists() {
        println!("cargo:warning=CUDA kernel file not found: {}", cuda_file);
        return;
    }

    // Find nvcc compiler
    let nvcc = match find_nvcc() {
        Some(path) => path,
        None => {
            println!("cargo:warning=nvcc not found. CUDA kernels will not be compiled.");
            println!("cargo:warning=Install CUDA Toolkit to enable GPU acceleration.");
            println!("cargo:warning=Falling back to placeholder kernels.");
            return;
        },
    };

    println!("cargo:warning=Found nvcc at: {}", nvcc);
    println!("cargo:warning=Compiling CUDA kernels to PTX...");

    // Compile CUDA to PTX
    let status = Command::new(&nvcc)
        .args([
            "-ptx",               // Compile to PTX (portable)
            "-O3",                // Optimization level 3
            "--use_fast_math",    // Fast math (slight precision trade-off)
            "-arch=sm_50",        // Minimum: Maxwell (CUDA 6.0+, GTX 700+)
            "--ptxas-options=-v", // Verbose PTX assembly
            "-o",
            ptx_file.to_str().unwrap(),
            cuda_file,
        ])
        .status()
        .expect("Failed to execute nvcc");

    if !status.success() {
        println!("cargo:warning=CUDA kernel compilation failed!");
        println!("cargo:warning=GPU acceleration will use placeholder kernels.");
        return;
    }

    println!(
        "cargo:warning=CUDA kernels compiled successfully: {}",
        ptx_file.display()
    );

    // Embed PTX path in the binary
    println!("cargo:rustc-env=CUDA_PTX_PATH={}", ptx_file.display());
}

/// Find nvcc CUDA compiler
#[cfg(feature = "use-cuda")]
fn find_nvcc() -> Option<String> {
    // Try CUDA_PATH environment variable
    if let Ok(cuda_path) = env::var("CUDA_PATH") {
        let nvcc_path = PathBuf::from(cuda_path).join("bin").join("nvcc");
        if nvcc_path.exists() {
            return Some(nvcc_path.to_string_lossy().to_string());
        }
    }

    // Try common installation paths
    let paths = vec![
        "/usr/local/cuda/bin/nvcc",
        "/usr/bin/nvcc",
        "/opt/cuda/bin/nvcc",
        "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v11.8\\bin\\nvcc.exe",
        "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.0\\bin\\nvcc.exe",
    ];

    for path in paths {
        if std::path::Path::new(path).exists() {
            return Some(path.to_string());
        }
    }

    // Try PATH environment variable
    if let Ok(output) = Command::new("which").arg("nvcc").output() {
        if output.status.success() {
            let path = String::from_utf8_lossy(&output.stdout).trim().to_string();
            if !path.is_empty() {
                return Some(path);
            }
        }
    }

    // Try Windows 'where' command
    if let Ok(output) = Command::new("where").arg("nvcc").output() {
        if output.status.success() {
            let path = String::from_utf8_lossy(&output.stdout)
                .lines()
                .next()
                .unwrap_or("")
                .trim()
                .to_string();
            if !path.is_empty() {
                return Some(path);
            }
        }
    }

    None
}
