//! A small HDF5 writer: groups, hard and soft links, attributes and
//! contiguous datasets of fixed-size elements, streamed in one pass.
//!
//! Written from the HDF5 File Format Specification version 3.0 (The HDF
//! Group): superblock version 2 (section II.A), version 2 object headers
//! (IV.A.1.b) with the messages Dataspace v2 (IV.A.2.b), Datatype v1
//! (IV.A.2.d), Fill Value v3 (IV.A.2.f), Link v1 (IV.A.2.g), Data Layout v3
//! (IV.A.2.i), Attribute v3 (IV.A.2.m), Link Info v0 (IV.A.2.c) and Group
//! Info v0 (IV.A.2.j). Links and attributes are stored compactly in the
//! object headers; the checksums are Bob Jenkins' lookup3 `hashlittle`
//! with initial value 0, as the specification prescribes (its appendix on
//! the Jenkins hash).
//!
//! The file is laid out as the superblock, every object header, then the
//! data of every dataset in the order the datasets were made. The headers
//! hold the data addresses, so all shapes are declared first ([`H5File`]),
//! the metadata is written, and then the data is streamed dataset by
//! dataset ([`DataWriter`]) without seeking: the file can go straight into
//! a hashing stream.

use std::io::{self, Write};

/// The "undefined address" of the specification (all bits set).
const UNDEFINED: u64 = u64::MAX;
const SUPERBLOCK_SIZE: u64 = 48;

/// The element type of an attribute or dataset (little-endian throughout).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Element {
    I64,
    U32,
    U8,
    F64,
    /// An ASCII string of fixed length, NUL-padded (HDF5 "null
    /// terminate": a string that fills the length has no terminator).
    Str(usize),
}

impl Element {
    pub fn size(self) -> usize {
        match self {
            Self::I64 | Self::F64 => 8,
            Self::U32 => 4,
            Self::U8 => 1,
            Self::Str(n) => n,
        }
    }

    /// The Datatype message (IV.A.2.d, version 1).
    fn encode(self, out: &mut Vec<u8>) {
        let size = self.size() as u32;
        match self {
            Self::I64 | Self::U32 | Self::U8 => {
                // class 0 fixed-point; bit 3 of the class bits: signed
                let signed = matches!(self, Self::I64);
                out.extend_from_slice(&[0x10, if signed { 0x08 } else { 0 }, 0, 0]);
                out.extend_from_slice(&size.to_le_bytes());
                out.extend_from_slice(&0u16.to_le_bytes()); // bit offset
                out.extend_from_slice(&((8 * size) as u16).to_le_bytes()); // precision
            }
            Self::F64 => {
                // class 1 floating-point: mantissa normalisation "implied
                // msb" (bits 4-5 = 2), sign bit 63; IEEE 754 binary64
                out.extend_from_slice(&[0x11, 0x20, 63, 0]);
                out.extend_from_slice(&size.to_le_bytes());
                out.extend_from_slice(&0u16.to_le_bytes()); // bit offset
                out.extend_from_slice(&64u16.to_le_bytes()); // precision
                out.extend_from_slice(&[52, 11, 0, 52]); // exponent and mantissa
                out.extend_from_slice(&1023u32.to_le_bytes()); // exponent bias
            }
            Self::Str(_) => {
                // class 3 string: padding 0 (null terminate), ASCII
                out.extend_from_slice(&[0x13, 0, 0, 0]);
                out.extend_from_slice(&size.to_le_bytes());
            }
        }
    }
}

/// An object of an [`H5File`]: a group or a dataset.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Obj(usize);

enum Target {
    Hard(usize),
    Soft(String),
}

struct Link {
    name: String,
    target: Target,
}

struct Attribute {
    name: String,
    element: Element,
    /// Empty for a scalar.
    shape: Vec<u64>,
    data: Vec<u8>,
}

enum Kind {
    Group(Vec<Link>),
    Dataset { element: Element, shape: Vec<u64> },
}

struct Object {
    kind: Kind,
    attributes: Vec<Attribute>,
}

/// The structure of a file: groups, links, attributes and dataset shapes.
pub struct H5File {
    objects: Vec<Object>,
}

impl Default for H5File {
    fn default() -> Self {
        Self::new()
    }
}

impl H5File {
    /// A file with an empty root group.
    pub fn new() -> Self {
        Self {
            objects: vec![Object {
                kind: Kind::Group(Vec::new()),
                attributes: Vec::new(),
            }],
        }
    }

    pub fn root(&self) -> Obj {
        Obj(0)
    }

    fn add(&mut self, parent: Obj, name: &str, kind: Kind) -> Obj {
        let id = self.objects.len();
        self.objects.push(Object {
            kind,
            attributes: Vec::new(),
        });
        self.link(parent, name, Target::Hard(id));
        Obj(id)
    }

    fn link(&mut self, parent: Obj, name: &str, target: Target) {
        if let Kind::Group(links) = &mut self.objects[parent.0].kind {
            links.push(Link {
                name: name.to_string(),
                target,
            });
        }
    }

    /// A new group `name` in `parent` (a group).
    pub fn group(&mut self, parent: Obj, name: &str) -> Obj {
        self.add(parent, name, Kind::Group(Vec::new()))
    }

    /// A new contiguous dataset; its data follows in [`DataWriter`].
    pub fn dataset(&mut self, parent: Obj, name: &str, element: Element, shape: &[u64]) -> Obj {
        self.add(
            parent,
            name,
            Kind::Dataset {
                element,
                shape: shape.to_vec(),
            },
        )
    }

    /// A soft link `name` in `parent` to the absolute path `path`.
    pub fn soft_link(&mut self, parent: Obj, name: &str, path: &str) {
        self.link(parent, name, Target::Soft(path.to_string()));
    }

    fn attribute(&mut self, obj: Obj, name: &str, element: Element, shape: &[u64], data: Vec<u8>) {
        self.objects[obj.0].attributes.push(Attribute {
            name: name.to_string(),
            element,
            shape: shape.to_vec(),
            data,
        });
    }

    /// A scalar 64-bit integer attribute.
    pub fn attr_i64(&mut self, obj: Obj, name: &str, value: i64) {
        self.attribute(obj, name, Element::I64, &[], value.to_le_bytes().to_vec());
    }

    /// A one-dimensional 64-bit integer attribute.
    pub fn attr_i64s(&mut self, obj: Obj, name: &str, values: &[i64]) {
        let data = values.iter().flat_map(|v| v.to_le_bytes()).collect();
        self.attribute(obj, name, Element::I64, &[values.len() as u64], data);
    }

    /// A scalar string attribute, stored with the length of the string.
    pub fn attr_str(&mut self, obj: Obj, name: &str, value: &str) {
        // a zero-sized type is not allowed; an empty string keeps one NUL
        let data = if value.is_empty() {
            vec![0]
        } else {
            value.as_bytes().to_vec()
        };
        self.attribute(obj, name, Element::Str(data.len()), &[], data);
    }

    /// Write the superblock and every object header to `out`, and return
    /// the writer of the datasets' data, which must follow at once.
    pub fn write<'a>(self, out: &'a mut dyn Write) -> io::Result<DataWriter<'a>> {
        // Header sizes do not depend on addresses (each address is eight
        // bytes), so one encoding with placeholder addresses sizes them.
        let zeros = vec![0u64; self.objects.len()];
        let mut sizes = Vec::with_capacity(self.objects.len());
        for object in &self.objects {
            sizes.push(object_header(object, &zeros, 0)?.len() as u64);
        }
        let mut header_at = Vec::with_capacity(sizes.len());
        let mut at = SUPERBLOCK_SIZE;
        for size in &sizes {
            header_at.push(at);
            at += size;
        }
        let mut data_at = vec![UNDEFINED; self.objects.len()];
        let mut datasets = Vec::new();
        for (i, object) in self.objects.iter().enumerate() {
            if let Kind::Dataset { element, shape } = &object.kind {
                let bytes = shape
                    .iter()
                    .try_fold(element.size() as u64, |n, &d| n.checked_mul(d))
                    .ok_or_else(|| io::Error::other("an HDF5 dataset is too large"))?;
                if bytes > 0 {
                    data_at[i] = at;
                    at += bytes;
                }
                datasets.push((i, bytes));
            }
        }
        out.write_all(&superblock(at, header_at[0]))?;
        for (object, &at) in self.objects.iter().zip(&data_at) {
            out.write_all(&object_header(object, &header_at, at)?)?;
        }
        Ok(DataWriter {
            out,
            buf: Vec::with_capacity(DataWriter::SPILL + 64),
            datasets,
            next: 0,
            current: None,
            written: 0,
        })
    }
}

/// The data of the datasets, written in the order they were made: call
/// [`DataWriter::begin`] for each dataset (also an empty one), then its
/// elements in row-major order, and [`DataWriter::finish`] at the end. A
/// dataset given more or fewer bytes than its shape holds is an error.
pub struct DataWriter<'a> {
    out: &'a mut dyn Write,
    buf: Vec<u8>,
    /// (object, bytes) of every dataset, in order.
    datasets: Vec<(usize, u64)>,
    next: usize,
    /// The dataset being written and the byte count it needs.
    current: Option<(usize, u64)>,
    written: u64,
}

impl DataWriter<'_> {
    const SPILL: usize = 1 << 16;

    fn close(&mut self) -> io::Result<()> {
        if let Some((obj, bytes)) = self.current.take()
            && self.written != bytes
        {
            return Err(io::Error::other(format!(
                "HDF5 object {obj}: {} bytes written for a dataset of {bytes}",
                self.written
            )));
        }
        Ok(())
    }

    /// Start the data of dataset `ds`, the next one in order.
    pub fn begin(&mut self, ds: Obj) -> io::Result<()> {
        self.close()?;
        match self.datasets.get(self.next) {
            Some(&(obj, bytes)) if obj == ds.0 => {
                self.current = Some((obj, bytes));
                self.next += 1;
                self.written = 0;
                Ok(())
            }
            _ => Err(io::Error::other(format!(
                "HDF5 object {} is not the next dataset",
                ds.0
            ))),
        }
    }

    fn put(&mut self, bytes: &[u8]) -> io::Result<()> {
        self.buf.extend_from_slice(bytes);
        self.written += bytes.len() as u64;
        if self.buf.len() >= Self::SPILL {
            self.out.write_all(&self.buf)?;
            self.buf.clear();
        }
        Ok(())
    }

    pub fn i64(&mut self, v: i64) -> io::Result<()> {
        self.put(&v.to_le_bytes())
    }

    /// A non-negative count or id as a 64-bit signed integer.
    pub fn id(&mut self, v: usize) -> io::Result<()> {
        let v = i64::try_from(v).map_err(|_| io::Error::other("an id exceeds 2^63"))?;
        self.i64(v)
    }

    pub fn u32(&mut self, v: u32) -> io::Result<()> {
        self.put(&v.to_le_bytes())
    }

    pub fn u8(&mut self, v: u8) -> io::Result<()> {
        self.put(&[v])
    }

    pub fn f64(&mut self, v: f64) -> io::Result<()> {
        self.put(&v.to_le_bytes())
    }

    /// One element of a string dataset of `width` bytes, NUL-padded.
    pub fn str(&mut self, s: &str, width: usize) -> io::Result<()> {
        if s.len() > width {
            return Err(io::Error::other("a string is longer than its dataset type"));
        }
        self.put(s.as_bytes())?;
        for _ in s.len()..width {
            self.put(&[0])?;
        }
        Ok(())
    }

    /// Check that every dataset got its data and hand on the rest.
    pub fn finish(mut self) -> io::Result<()> {
        self.close()?;
        if self.next != self.datasets.len() {
            return Err(io::Error::other(format!(
                "{} HDF5 datasets got no data",
                self.datasets.len() - self.next
            )));
        }
        self.out.write_all(&self.buf)
    }
}

fn superblock(end_of_file: u64, root: u64) -> Vec<u8> {
    let mut s = Vec::with_capacity(SUPERBLOCK_SIZE as usize);
    s.extend_from_slice(b"\x89HDF\r\n\x1a\n");
    // version 2, offsets and lengths of eight bytes, no consistency flags
    s.extend_from_slice(&[2, 8, 8, 0]);
    s.extend_from_slice(&0u64.to_le_bytes()); // base address
    s.extend_from_slice(&UNDEFINED.to_le_bytes()); // no superblock extension
    s.extend_from_slice(&end_of_file.to_le_bytes());
    s.extend_from_slice(&root.to_le_bytes());
    let sum = lookup3(&s);
    s.extend_from_slice(&sum.to_le_bytes());
    s
}

/// One header message: type, data size, flags (none), data.
fn message(out: &mut Vec<u8>, kind: u8, data: &[u8]) -> io::Result<()> {
    let size = u16::try_from(data.len())
        .map_err(|_| io::Error::other("an HDF5 header message exceeds 64 KiB"))?;
    out.push(kind);
    out.extend_from_slice(&size.to_le_bytes());
    out.push(0);
    out.extend_from_slice(data);
    Ok(())
}

/// Dataspace message version 2: scalar (no dimensions) or simple.
fn dataspace(shape: &[u64], out: &mut Vec<u8>) {
    let kind = if shape.is_empty() { 0 } else { 1 };
    out.extend_from_slice(&[2, shape.len() as u8, 0, kind]);
    for d in shape {
        out.extend_from_slice(&d.to_le_bytes());
    }
}

/// A version 2 object header with all its messages in chunk 0.
/// `data_at` is the address of a dataset's data (undefined when empty).
fn object_header(object: &Object, header_at: &[u64], data_at: u64) -> io::Result<Vec<u8>> {
    let mut m = Vec::new();
    let mut d = Vec::new();
    match &object.kind {
        Kind::Group(links) => {
            // Link Info v0: creation order tracked and indexed (VTK lists
            // the blocks of a collection in creation order), the next
            // creation index, and compact links (no heap, no indices)
            d.extend_from_slice(&[0, 0x03]);
            d.extend_from_slice(&(links.len() as u64).to_le_bytes());
            d.extend_from_slice(&UNDEFINED.to_le_bytes());
            d.extend_from_slice(&UNDEFINED.to_le_bytes());
            d.extend_from_slice(&UNDEFINED.to_le_bytes());
            message(&mut m, 0x02, &d)?;
            // Group Info v0 with the default phase changes
            message(&mut m, 0x0A, &[0, 0])?;
            for (order, link) in links.iter().enumerate() {
                d.clear();
                let name = link.name.as_bytes();
                let long = name.len() > usize::from(u8::MAX);
                // flags: name length in 1 or 2 bytes; bit 2 the creation
                // order; bit 3 the link type
                let mut flags = u8::from(long) | 0x04;
                if matches!(link.target, Target::Soft(_)) {
                    flags |= 0x08;
                }
                d.extend_from_slice(&[1, flags]);
                if let Target::Soft(_) = link.target {
                    d.push(1);
                }
                d.extend_from_slice(&(order as u64).to_le_bytes());
                if long {
                    let n = u16::try_from(name.len())
                        .map_err(|_| io::Error::other("an HDF5 link name is too long"))?;
                    d.extend_from_slice(&n.to_le_bytes());
                } else {
                    d.push(name.len() as u8);
                }
                d.extend_from_slice(name);
                match &link.target {
                    Target::Hard(id) => d.extend_from_slice(&header_at[*id].to_le_bytes()),
                    Target::Soft(path) => {
                        let n = u16::try_from(path.len())
                            .map_err(|_| io::Error::other("an HDF5 soft link is too long"))?;
                        d.extend_from_slice(&n.to_le_bytes());
                        d.extend_from_slice(path.as_bytes());
                    }
                }
                message(&mut m, 0x06, &d)?;
            }
        }
        Kind::Dataset { element, shape } => {
            d.clear();
            dataspace(shape, &mut d);
            message(&mut m, 0x01, &d)?;
            d.clear();
            element.encode(&mut d);
            message(&mut m, 0x03, &d)?;
            // Fill Value v3: allocation late, written if set, none set
            message(&mut m, 0x05, &[3, 0x0A])?;
            // Data Layout v3, contiguous
            d.clear();
            d.extend_from_slice(&[3, 1]);
            let size = shape
                .iter()
                .fold(element.size() as u64, |n, &x| n.saturating_mul(x));
            d.extend_from_slice(&data_at.to_le_bytes());
            d.extend_from_slice(&size.to_le_bytes());
            message(&mut m, 0x08, &d)?;
        }
    }
    for a in &object.attributes {
        d.clear();
        let mut datatype = Vec::new();
        a.element.encode(&mut datatype);
        let mut space = Vec::new();
        dataspace(&a.shape, &mut space);
        let name_size = u16::try_from(a.name.len() + 1)
            .map_err(|_| io::Error::other("an HDF5 attribute name is too long"))?;
        d.extend_from_slice(&[3, 0]);
        d.extend_from_slice(&name_size.to_le_bytes());
        d.extend_from_slice(&(datatype.len() as u16).to_le_bytes());
        d.extend_from_slice(&(space.len() as u16).to_le_bytes());
        d.push(0); // ASCII name
        d.extend_from_slice(a.name.as_bytes());
        d.push(0);
        d.extend_from_slice(&datatype);
        d.extend_from_slice(&space);
        d.extend_from_slice(&a.data);
        message(&mut m, 0x0C, &d)?;
    }
    let chunk = u32::try_from(m.len())
        .map_err(|_| io::Error::other("an HDF5 object header exceeds 4 GiB"))?;
    let mut h = Vec::with_capacity(m.len() + 14);
    // version 2; flags: the size of chunk 0 in four bytes
    h.extend_from_slice(b"OHDR");
    h.extend_from_slice(&[2, 0x02]);
    h.extend_from_slice(&chunk.to_le_bytes());
    h.extend_from_slice(&m);
    let sum = lookup3(&h);
    h.extend_from_slice(&sum.to_le_bytes());
    Ok(h)
}

/// Bob Jenkins' lookup3 `hashlittle` (public domain, 2006) with initial
/// value 0, byte by byte: the checksum of HDF5 metadata.
pub fn lookup3(key: &[u8]) -> u32 {
    fn mix(a: &mut u32, b: &mut u32, c: &mut u32) {
        *a = a.wrapping_sub(*c);
        *a ^= c.rotate_left(4);
        *c = c.wrapping_add(*b);
        *b = b.wrapping_sub(*a);
        *b ^= a.rotate_left(6);
        *a = a.wrapping_add(*c);
        *c = c.wrapping_sub(*b);
        *c ^= b.rotate_left(8);
        *b = b.wrapping_add(*a);
        *a = a.wrapping_sub(*c);
        *a ^= c.rotate_left(16);
        *c = c.wrapping_add(*b);
        *b = b.wrapping_sub(*a);
        *b ^= a.rotate_left(19);
        *a = a.wrapping_add(*c);
        *c = c.wrapping_sub(*b);
        *c ^= b.rotate_left(4);
        *b = b.wrapping_add(*a);
    }
    fn last(a: &mut u32, b: &mut u32, c: &mut u32) {
        *c ^= *b;
        *c = c.wrapping_sub(b.rotate_left(14));
        *a ^= *c;
        *a = a.wrapping_sub(c.rotate_left(11));
        *b ^= *a;
        *b = b.wrapping_sub(a.rotate_left(25));
        *c ^= *b;
        *c = c.wrapping_sub(b.rotate_left(16));
        *a ^= *c;
        *a = a.wrapping_sub(c.rotate_left(4));
        *b ^= *a;
        *b = b.wrapping_sub(a.rotate_left(14));
        *c ^= *b;
        *c = c.wrapping_sub(b.rotate_left(24));
    }
    fn word(k: &[u8]) -> u32 {
        k.iter()
            .enumerate()
            .fold(0u32, |w, (i, &b)| w.wrapping_add(u32::from(b) << (8 * i)))
    }
    let init = 0xdead_beef_u32.wrapping_add(key.len() as u32);
    let (mut a, mut b, mut c) = (init, init, init);
    let mut k = key;
    while k.len() > 12 {
        a = a.wrapping_add(word(&k[0..4]));
        b = b.wrapping_add(word(&k[4..8]));
        c = c.wrapping_add(word(&k[8..12]));
        mix(&mut a, &mut b, &mut c);
        k = &k[12..];
    }
    if k.is_empty() {
        return c;
    }
    // the last 1 to 12 bytes, zero-extended
    a = a.wrapping_add(word(&k[..k.len().min(4)]));
    if k.len() > 4 {
        b = b.wrapping_add(word(&k[4..k.len().min(8)]));
    }
    if k.len() > 8 {
        c = c.wrapping_add(word(&k[8..]));
    }
    last(&mut a, &mut b, &mut c);
    c
}

/// A reader of exactly the files [`H5File`] writes, for tests: every
/// object by path, with checksums verified.
#[cfg(test)]
pub(crate) mod read {
    use std::collections::BTreeMap;

    use super::lookup3;

    #[derive(Debug, Default)]
    pub(crate) struct Item {
        pub attrs: BTreeMap<String, Vec<u8>>,
        /// Datasets: the shape and the raw little-endian data.
        pub shape: Option<Vec<u64>>,
        pub data: Vec<u8>,
        /// Links of a group, in creation order.
        pub links: Vec<String>,
        pub soft: Option<String>,
    }

    impl Item {
        pub fn i64s(&self) -> Vec<i64> {
            let (words, _) = self.data.as_chunks::<8>();
            words.iter().map(|c| i64::from_le_bytes(*c)).collect()
        }

        pub fn f64s(&self) -> Vec<f64> {
            let (words, _) = self.data.as_chunks::<8>();
            words.iter().map(|c| f64::from_le_bytes(*c)).collect()
        }

        pub fn u32s(&self) -> Vec<u32> {
            let (words, _) = self.data.as_chunks::<4>();
            words.iter().map(|c| u32::from_le_bytes(*c)).collect()
        }
    }

    fn u64_at(b: &[u8], at: usize) -> u64 {
        u64::from_le_bytes(b[at..at + 8].try_into().unwrap())
    }

    fn u16_at(b: &[u8], at: usize) -> usize {
        usize::from(u16::from_le_bytes(b[at..at + 2].try_into().unwrap()))
    }

    pub(crate) fn read(file: &[u8]) -> BTreeMap<String, Item> {
        assert_eq!(&file[..8], b"\x89HDF\r\n\x1a\n");
        assert_eq!(file[8], 2, "superblock version");
        let sum = u32::from_le_bytes(file[44..48].try_into().unwrap());
        assert_eq!(sum, lookup3(&file[..44]), "superblock checksum");
        assert_eq!(u64_at(file, 28) as usize, file.len(), "end of file address");
        let mut items = BTreeMap::new();
        object(file, u64_at(file, 36) as usize, String::new(), &mut items);
        items
    }

    fn object(file: &[u8], at: usize, path: String, items: &mut BTreeMap<String, Item>) {
        let h = &file[at..];
        assert_eq!(&h[..4], b"OHDR");
        assert_eq!(h[4], 2);
        assert_eq!(h[5], 0x02, "chunk size in four bytes, nothing optional");
        let chunk = u32::from_le_bytes(h[6..10].try_into().unwrap()) as usize;
        let end = 10 + chunk;
        let sum = u32::from_le_bytes(h[end..end + 4].try_into().unwrap());
        assert_eq!(sum, lookup3(&h[..end]), "object header checksum at {at}");
        let mut item = Item::default();
        let mut children = Vec::new();
        let mut layout = None;
        let mut k = 10;
        while k < end {
            let (kind, size) = (h[k], u16_at(h, k + 1));
            let d = &h[k + 4..k + 4 + size];
            k += 4 + size;
            match kind {
                0x01 => {
                    assert_eq!(d[0], 2);
                    let rank = usize::from(d[1]);
                    item.shape = Some((0..rank).map(|i| u64_at(d, 4 + 8 * i)).collect());
                }
                0x06 => {
                    let flags = d[1];
                    let mut i = 2;
                    let soft = flags & 0x08 != 0 && d[i] == 1;
                    if flags & 0x08 != 0 {
                        i += 1;
                    }
                    if flags & 0x04 != 0 {
                        i += 8;
                    }
                    let n = if flags & 0x03 == 0 {
                        i += 1;
                        usize::from(d[i - 1])
                    } else {
                        i += 2;
                        u16_at(d, i - 2)
                    };
                    let name = String::from_utf8(d[i..i + n].to_vec()).unwrap();
                    i += n;
                    item.links.push(name.clone());
                    let child = format!("{path}/{name}");
                    if soft {
                        let len = u16_at(d, i);
                        let target = String::from_utf8(d[i + 2..i + 2 + len].to_vec()).unwrap();
                        items.insert(
                            child,
                            Item {
                                soft: Some(target),
                                ..Item::default()
                            },
                        );
                    } else {
                        children.push((u64_at(d, i) as usize, child));
                    }
                }
                0x08 => {
                    assert_eq!(&d[..2], &[3, 1], "contiguous layout");
                    layout = Some((u64_at(d, 2), u64_at(d, 10) as usize));
                }
                0x0C => {
                    let name_size = u16_at(d, 2);
                    let (type_size, space_size) = (u16_at(d, 4), u16_at(d, 6));
                    let name = String::from_utf8(d[9..9 + name_size - 1].to_vec()).unwrap();
                    let data = d[9 + name_size + type_size + space_size..].to_vec();
                    item.attrs.insert(name, data);
                }
                _ => {}
            }
        }
        if let Some((address, size)) = layout {
            item.data = if size == 0 {
                assert_eq!(address, u64::MAX, "an empty dataset has no address");
                Vec::new()
            } else {
                file[address as usize..address as usize + size].to_vec()
            };
        }
        items.insert(path, item);
        for (address, child) in children {
            object(file, address, child, items);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Invariant: the checksum is lookup3 `hashlittle` with initial value
    /// 0, which HDF5 checks on every version 2 superblock and object header.
    /// Breaks if: a rotation, the `> 12` loop bound or the tail handling
    /// changes. No other test: a wrong checksum only shows when a real
    /// HDF5 reader opens the file. Values from the comments of lookup3.c.
    #[test]
    fn lookup3_matches_the_published_values() {
        assert_eq!(lookup3(b""), 0xdead_beef);
        assert_eq!(lookup3(b"Four score and seven years ago"), 0x1777_0551);
    }

    /// Invariant: the data of every dataset is checked against its shape,
    /// and the datasets come in the order they were made.
    /// Breaks if: `begin` or `finish` stop counting. No other test: a short
    /// dataset gives a file whose last bytes are missing, which only a
    /// reader notices.
    #[test]
    fn datasets_must_get_exactly_their_bytes_in_order() {
        let build = || {
            let mut f = H5File::new();
            let root = f.root();
            let a = f.dataset(root, "a", Element::I64, &[2]);
            let b = f.dataset(root, "b", Element::U8, &[0]);
            (f, a, b)
        };
        let (f, a, b) = build();
        let mut out = Vec::new();
        let mut w = f.write(&mut out).unwrap();
        w.begin(a).unwrap();
        w.i64(1).unwrap();
        w.i64(2).unwrap();
        w.begin(b).unwrap();
        w.finish().unwrap();
        // the file ends with the two values
        assert_eq!(&out[out.len() - 16..out.len() - 8], &1i64.to_le_bytes());

        let (f, a, b) = build();
        let mut out = Vec::new();
        let mut w = f.write(&mut out).unwrap();
        w.begin(a).unwrap();
        w.i64(1).unwrap();
        assert!(w.begin(b).is_err(), "a short dataset is an error");

        let (f, _, b) = build();
        let mut out = Vec::new();
        let mut w = f.write(&mut out).unwrap();
        assert!(w.begin(b).is_err(), "datasets come in order");
    }
}
