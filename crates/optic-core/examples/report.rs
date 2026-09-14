//! Print first-order data and a spot summary for the sample systems.

use optic_core::{launch, lines, samples, trace, Field, Paraxial, System};

fn rms_spot(
    sys: &System<f64>,
    par: &Paraxial<f64>,
    field: Field,
    wl: f64,
    n: usize,
) -> (f64, usize) {
    let mut pts = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let px = -1.0 + 2.0 * (i as f64 + 0.5) / n as f64;
            let py = -1.0 + 2.0 * (j as f64 + 0.5) / n as f64;
            if px * px + py * py > 1.0 {
                continue;
            }
            let r = trace(sys, wl, launch(sys, par, field, px, py));
            if let Some(p) = r.image_point() {
                pts.push((p.x, p.y));
            }
        }
    }
    if pts.is_empty() {
        return (f64::NAN, 0);
    }
    let cx = pts.iter().map(|p| p.0).sum::<f64>() / pts.len() as f64;
    let cy = pts.iter().map(|p| p.1).sum::<f64>() / pts.len() as f64;
    let var = pts
        .iter()
        .map(|p| (p.0 - cx).powi(2) + (p.1 - cy).powi(2))
        .sum::<f64>()
        / pts.len() as f64;
    (var.sqrt(), pts.len())
}

fn report(name: &str, sys: &System<f64>) {
    let par = Paraxial::compute(sys, lines::D);
    println!("\n=== {name} : {} ===", sys.title);
    println!("  EFL          {:.6} mm", par.efl);
    println!("  BFD          {:.6} mm", par.bfd);
    println!("  EPD          {:.6} mm", par.epd);
    println!("  EP z         {:.6} mm (from surface 1)", par.ep_z);
    println!("  XP z         {:.6} mm (from image)", par.xp_z);
    println!("  working f/#  {:.6}", par.fno);
    let fields: Vec<Field> = sys.fields.clone();
    for f in fields {
        let (rms, n) = rms_spot(sys, &par, f, lines::D, 24);
        let h = par.image_height(sys, lines::D, f);
        println!(
            "  field {:>5.1}: rms spot {:9.4} um  ({n:3} rays)  paraxial image y {:9.4} mm",
            f.radius(),
            rms * 1000.0,
            h
        );
    }
}

fn main() {
    report("singlet", &samples::singlet());
    report("cooke", &samples::cooke_triplet());

    println!("\n=== catalog ===");
    for (name, m, nd, vd) in optic_core::catalog::PUBLISHED {
        println!(
            "  {name:9} nd {:.6} (published {nd:.5}, d {:+.2e})   Vd {:.4} (published {vd:.2}, d {:+.3})",
            m.nd(),
            m.nd() - nd,
            m.abbe(),
            m.abbe() - vd
        );
    }
}
