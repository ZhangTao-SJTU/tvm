// SPDX-License-Identifier: GPL-3.0-or-later
//! MT19937-64 and Marsaglia polar normals. Each accepted pair returns the
//! y sample and caches the x sample for the next draw.
#[derive(Clone, Debug)]
pub struct Rng {
    words: [u64; 312],
    index: usize,
    saved: Option<f64>,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        let mut words = [0; 312];
        words[0] = seed;
        for i in 1..312 {
            words[i] = 6364136223846793005u64
                .wrapping_mul(words[i - 1] ^ (words[i - 1] >> 62))
                .wrapping_add(i as u64);
        }
        Self {
            words,
            index: 312,
            saved: None,
        }
    }

    pub fn next_u64(&mut self) -> u64 {
        if self.index == 312 {
            for i in 0..312 {
                let x =
                    (self.words[i] & 0xffffffff80000000) | (self.words[(i + 1) % 312] & 0x7fffffff);
                self.words[i] = self.words[(i + 156) % 312]
                    ^ (x >> 1)
                    ^ if x & 1 != 0 { 0xb5026f5aa96619e9 } else { 0 };
            }
            self.index = 0;
        }
        let mut x = self.words[self.index];
        self.index += 1;
        x ^= (x >> 29) & 0x5555555555555555;
        x ^= (x << 17) & 0x71d67fffeda60000;
        x ^= (x << 37) & 0xfff7eee000000000;
        x ^= x >> 43;
        x
    }

    fn uniform(&mut self) -> f64 {
        let value = self.next_u64() as f64 / 18446744073709551616.0;
        value.min(f64::from_bits(1.0f64.to_bits() - 1))
    }

    pub fn normal(&mut self) -> f64 {
        if let Some(x) = self.saved.take() {
            return x;
        }
        loop {
            let x = 2.0 * self.uniform() - 1.0;
            let y = 2.0 * self.uniform() - 1.0;
            let r2 = x * x + y * y;
            if r2 > 0.0 && r2 <= 1.0 {
                let factor = (-2.0 * r2.ln() / r2).sqrt();
                self.saved = Some(x * factor);
                return y * factor;
            }
        }
    }

    /// Final random state for audit; this is not a simulation checkpoint.
    pub fn state(&self) -> String {
        let mut result = String::from("mt19937_64\n");
        for word in self.words {
            result.push_str(&format!("{word} "));
        }
        result.push_str(&format!("{}\nnormal_saved {:?}\n", self.index, self.saved));
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn standard_mt19937_64_reference() {
        let mut rng = Rng::new(5489);
        // Standard MT19937-64 reference sequence.
        let expected = [
            14514284786278117030,
            4620546740167642908,
            13109570281517897720,
            17462938647148434322,
            355488278567739596,
        ];
        for x in expected {
            assert_eq!(rng.next_u64(), x);
        }
    }
    #[test]
    fn cloned_normal_state_continues_exactly() {
        let mut a = Rng::new(20261028);
        a.normal();
        let mut b = a.clone();
        for _ in 0..1000 {
            assert_eq!(a.normal(), b.normal());
        }
        assert_eq!(a.state(), b.state());
    }
}
