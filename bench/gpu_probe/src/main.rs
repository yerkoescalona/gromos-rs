//! Feasibility probe for BENCHMARKING.md "GPU": pair-arithmetic throughput of the
//! integrated GPU (wgpu/Vulkan, f32 and f64) vs the CPU (f64, 1 and 8 cores).
//! Standalone: `cd bench/gpu_probe && cargo run --release`.

use wgpu::util::DeviceExt;
const N: u32 = 1 << 20;   // independent "pairs" per dispatch
const ITERS: u32 = 64;    // dependent LJ+CRF evaluations per pair, to amortise dispatch overhead
const REPS: u32 = 20;

fn kernel(f64: bool) -> String {
    let (ty, cv) = if f64 { ("f64", "f64") } else { ("f32", "f32") };
    format!(r#"
        @group(0) @binding(0) var<storage, read> x: array<f32>;
        @group(0) @binding(1) var<storage, read_write> y: array<f32>;
        @compute @workgroup_size(256)
        fn main(@builtin(global_invocation_id) id: vec3<u32>) {{
            let i = id.x;
            if (i >= arrayLength(&x)) {{ return; }}
            var r2: {ty} = {cv}(x[i]) + {cv}(0.5);
            let c6: {ty} = {cv}(0.0026); let c12: {ty} = {cv}(0.0000026); let q: {ty} = {cv}(0.17);
            let k1: {ty} = {cv}(0.3644); let k2: {ty} = {cv}(0.1822); let cut: {ty} = {cv}(1.2);
            var acc: {ty} = {cv}(0.0);
            for (var it = 0u; it < {iters}u; it = it + 1u) {{
                let inv_r2 = {cv}(1.0) / r2;
                let inv_r6 = inv_r2 * inv_r2 * inv_r2;
                let e_lj = (c12 * inv_r6 - c6) * inv_r6;
                let f_lj = ({cv}(12.0) * c12 * inv_r6 - {cv}(6.0) * c6) * inv_r6 * inv_r2;
                let inv_r = sqrt(inv_r2);
                let e_crf = q * (inv_r - k1 * r2 - cut);
                let f_crf = q * (inv_r * inv_r2 + k2);
                acc = acc + f_lj + f_crf + e_lj + e_crf;
                r2 = r2 + {cv}(1e-3) * (acc - acc) + {cv}(1e-4);   // keep the chain dependent
            }}
            y[i] = f32(acc);
        }}"#, ty = ty, cv = cv, iters = ITERS)
}

fn cpu_eval(x: &[f32]) -> f64 {
    let mut total = 0.0f64;
    for &xi in x {
        let mut r2 = xi as f64 + 0.5;
        let (c6, c12, q, k1, k2, cut) = (0.0026f64, 0.0000026f64, 0.17f64, 0.3644f64, 0.1822f64, 1.2f64);
        let mut acc = 0.0f64;
        for _ in 0..ITERS {
            let inv_r2 = 1.0 / r2;
            let inv_r6 = inv_r2 * inv_r2 * inv_r2;
            let e_lj = (c12 * inv_r6 - c6) * inv_r6;
            let f_lj = (12.0 * c12 * inv_r6 - 6.0 * c6) * inv_r6 * inv_r2;
            let inv_r = inv_r2.sqrt();
            let e_crf = q * (inv_r - k1 * r2 - cut);
            let f_crf = q * (inv_r * inv_r2 + k2);
            acc += f_lj + f_crf + e_lj + e_crf;
            r2 += 1e-3 * (acc - acc) + 1e-4;
        }
        total += acc;
    }
    total
}

fn main() {
    let x: Vec<f32> = (0..N).map(|i| (i % 1000) as f32 * 1e-3).collect();
    let evals = (N as f64) * (ITERS as f64);

    // CPU, 1 core and 8 cores (std threads), f64.
    let t = std::time::Instant::now(); let s1 = cpu_eval(&x); let cpu1 = t.elapsed().as_secs_f64();
    let t = std::time::Instant::now();
    let s8: f64 = std::thread::scope(|sc| {
        let hs: Vec<_> = x.chunks(x.len() / 8).map(|c| sc.spawn(move || cpu_eval(c))).collect();
        hs.into_iter().map(|h| h.join().unwrap()).sum()
    });
    let cpu8 = t.elapsed().as_secs_f64();
    println!("CPU f64  1 core : {:7.1} M pair-evals/s   (checksum {:.3e})", evals / cpu1 / 1e6, s1);
    println!("CPU f64  8 cores: {:7.1} M pair-evals/s   (checksum {:.3e})", evals / cpu8 / 1e6, s8);

    let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor { backends: wgpu::Backends::VULKAN, ..Default::default() });
    let adapter = instance.enumerate_adapters(wgpu::Backends::VULKAN).into_iter()
        .find(|a| a.get_info().device_type == wgpu::DeviceType::IntegratedGpu).expect("iGPU");
    println!("GPU: {}", adapter.get_info().name);
    let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
        label: None, required_features: wgpu::Features::SHADER_F64, required_limits: wgpu::Limits::default(), memory_hints: Default::default(),
    }, None)).expect("device with f64");

    for (label, is_f64) in [("GPU f32", false), ("GPU f64", true)] {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: None, source: wgpu::ShaderSource::Wgsl(kernel(is_f64).into()) });
        let xb = device.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: None, contents: bytemuck::cast_slice(&x), usage: wgpu::BufferUsages::STORAGE });
        let yb = device.create_buffer(&wgpu::BufferDescriptor { label: None, size: (N as u64) * 4, usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC, mapped_at_creation: false });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor { label: None, layout: None, module: &shader, entry_point: Some("main"), compilation_options: Default::default(), cache: None });
        let bg = device.create_bind_group(&wgpu::BindGroupDescriptor { label: None, layout: &pipeline.get_bind_group_layout(0), entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: xb.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: yb.as_entire_binding() } ] });
        // warm-up
        let mut enc = device.create_command_encoder(&Default::default());
        { let mut p = enc.begin_compute_pass(&Default::default()); p.set_pipeline(&pipeline); p.set_bind_group(0, &bg, &[]); p.dispatch_workgroups(N / 256, 1, 1); }
        queue.submit([enc.finish()]); device.poll(wgpu::Maintain::Wait);
        let t = std::time::Instant::now();
        for _ in 0..REPS {
            let mut enc = device.create_command_encoder(&Default::default());
            { let mut p = enc.begin_compute_pass(&Default::default()); p.set_pipeline(&pipeline); p.set_bind_group(0, &bg, &[]); p.dispatch_workgroups(N / 256, 1, 1); }
            queue.submit([enc.finish()]);
        }
        device.poll(wgpu::Maintain::Wait);
        let per = t.elapsed().as_secs_f64() / REPS as f64;
        println!("{label}       : {:7.1} M pair-evals/s   ({:.2} ms per dispatch of {} pairs x {} evals)", evals / per / 1e6, per * 1e3, N, ITERS);
    }
}
