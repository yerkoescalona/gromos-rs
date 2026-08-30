//! gromos++ `utils::Time`: `@time <t0> <dt>` overrides the times in the trajectory (t0 + n·dt);
//! without it the frame's own time is used. Printed fixed with nine decimals in a 15-wide field.

use crate::args::Arguments;

#[derive(Debug, Clone)]
pub struct Time {
    t0: f64,
    dt: f64,
    read: bool,
    steps: usize,
    current: f64,
}

impl Time {
    pub fn from_args(args: &Arguments) -> Result<Self, String> {
        let v = args.values("time");
        let read = v.is_empty();
        let t0 = v
            .first()
            .map(|s| {
                s.parse::<f64>()
                    .map_err(|_| "@time: starting time is not numeric".to_string())
            })
            .transpose()?
            .unwrap_or(0.0);
        let dt = v
            .get(1)
            .map(|s| {
                s.parse::<f64>()
                    .map_err(|_| "@time: dt is not numeric".to_string())
            })
            .transpose()?
            .unwrap_or(1.0);
        Ok(Time {
            t0,
            dt,
            read,
            steps: 0,
            current: 0.0,
        })
    }

    /// Advance with a frame: its own time, or t0 + n·dt when `@time` was given.
    pub fn next(&mut self, frame_time: f64) -> f64 {
        self.current = if self.read {
            frame_time
        } else {
            self.t0 + self.steps as f64 * self.dt
        };
        self.steps += 1;
        self.current
    }

    pub fn current(&self) -> f64 {
        self.current
    }

    /// gromos++ `Time::print`: `setw(15) << setprecision(9)`, fixed.
    pub fn fmt(t: f64) -> String {
        format!("{:15.9}", t)
    }
}
