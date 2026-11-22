//! Geometry primitives used by the CPU reference mesher.

use std::ops::{Add, Div, Mul, Sub};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);

    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(self, other: Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    pub fn length(self) -> f64 {
        self.x.hypot(self.y).hypot(self.z)
    }

    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    pub fn max_abs_component(self) -> f64 {
        self.x.abs().max(self.y.abs()).max(self.z.abs())
    }
}

impl Add for Vec3 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl Sub for Vec3 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, factor: f64) -> Self {
        Self::new(self.x * factor, self.y * factor, self.z * factor)
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;

    fn div(self, divisor: f64) -> Self {
        Self::new(self.x / divisor, self.y / divisor, self.z / divisor)
    }
}

/// A point on the plane and a unit normal pointing away from the retained half-space.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane {
    point: Vec3,
    normal: Vec3,
}

impl Plane {
    /// Returns `None` for a nonfinite point or a zero or nonfinite normal.
    pub fn new(point: Vec3, normal: Vec3) -> Option<Self> {
        let scale = normal.max_abs_component();
        if !point.is_finite() || !normal.is_finite() || scale == 0.0 {
            return None;
        }
        let scaled = normal / scale;
        Some(Self {
            point,
            normal: scaled / scaled.length(),
        })
    }

    pub fn normal(self) -> Vec3 {
        self.normal
    }

    /// Negative distances are on the retained side.
    pub fn signed_distance(self, point: Vec3) -> f64 {
        (point - self.point).dot(self.normal)
    }
}

/// A box with finite bounds and strictly positive, representable extent on every axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisAlignedBox {
    min: Vec3,
    max: Vec3,
}

impl AxisAlignedBox {
    pub fn new(min: Vec3, max: Vec3) -> Option<Self> {
        let extent = max - min;
        if !extent.is_finite() || extent.x <= 0.0 || extent.y <= 0.0 || extent.z <= 0.0 {
            return None;
        }
        Some(Self { min, max })
    }

    pub fn min(self) -> Vec3 {
        self.min
    }

    pub fn max(self) -> Vec3 {
        self.max
    }

    pub fn extent(self) -> Vec3 {
        self.max - self.min
    }

    pub fn contains_interior(self, point: Vec3) -> bool {
        self.min.x < point.x
            && point.x < self.max.x
            && self.min.y < point.y
            && point.y < self.max.y
            && self.min.z < point.z
            && point.z < self.max.z
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
pub enum BoxSide {
    XMin,
    XMax,
    YMin,
    YMax,
    ZMin,
    ZMax,
}

/// Index into the surface domain's named patch table.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PatchId(pub usize);

/// Index into the surface domain's boundary triangle table.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct FacetId(pub usize);

impl BoxSide {
    pub const ALL: [Self; 6] = [
        Self::XMin,
        Self::XMax,
        Self::YMin,
        Self::YMax,
        Self::ZMin,
        Self::ZMax,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Self::XMin => "xmin",
            Self::XMax => "xmax",
            Self::YMin => "ymin",
            Self::YMax => "ymax",
            Self::ZMin => "zmin",
            Self::ZMax => "zmax",
        }
    }

    pub fn axis(self) -> usize {
        self as usize / 2
    }

    pub fn is_max(self) -> bool {
        matches!(self, Self::XMax | Self::YMax | Self::ZMax)
    }
}
