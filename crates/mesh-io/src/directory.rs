use std::ffi::OsStr;
use std::fmt;
use std::io::{self, Read, Write};
use std::path::{Component, Path};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::thread::{self, JoinHandle};

use cap_std::ambient_authority;
use cap_std::fs::{Dir, File, OpenOptions};

/// Access explicitly granted to one working directory.
pub struct MeshDirectory(Dir);

/// Why a streamed file was not written. The file is removed either way.
#[derive(Debug)]
pub enum StreamError {
    /// The file system refused: creating, writing or syncing the file, or
    /// a name that is not a single file name.
    Io(io::Error),
    /// The writer of the content failed.
    Content(io::Error),
}

impl fmt::Display for StreamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(error) | Self::Content(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for StreamError {}

impl MeshDirectory {
    pub fn open(path: &Path) -> io::Result<Self> {
        Dir::open_ambient_dir(path, ambient_authority()).map(Self)
    }

    pub fn read_regular(&self, name: &OsStr, limit: usize) -> io::Result<Vec<u8>> {
        let path = single_entry(name)?;
        if !self.0.symlink_metadata(path)?.is_file() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "expected a regular file",
            ));
        }
        let file = self.0.open(path)?;
        let mut bytes = Vec::new();
        file.take((limit as u64).saturating_add(1))
            .read_to_end(&mut bytes)?;
        if bytes.len() > limit {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "file exceeds the size limit",
            ));
        }
        Ok(bytes)
    }

    pub fn write_new(&self, name: &OsStr, bytes: &[u8]) -> io::Result<()> {
        let path = single_entry(name)?;
        let mut file = self
            .0
            .open_with(path, OpenOptions::new().write(true).create_new(true))?;
        if let Err(error) = file.write_all(bytes).and_then(|_| file.sync_all()) {
            drop(file);
            return Err(with_cleanup(error, self.0.remove_file(path)));
        }
        Ok(())
    }

    /// Create the new file `name` and stream `content` into it. A second
    /// thread hashes and writes each block of [`BLOCK`] bytes while the
    /// next one is produced, so a file of any size needs a few blocks of
    /// memory. Returns the BLAKE3 hash of the file as lower-case hex
    /// ([`FILE_HASH`]). On any
    /// error the file is removed.
    pub fn stream_new(
        &self,
        name: &OsStr,
        content: impl FnOnce(&mut dyn Write) -> io::Result<()>,
    ) -> Result<String, StreamError> {
        let path = single_entry(name).map_err(StreamError::Io)?;
        let file = self
            .0
            .open_with(path, OpenOptions::new().write(true).create_new(true))
            .map_err(StreamError::Io)?;
        let mut out = HashedFile::new(file);
        let result = match content(&mut out) {
            Ok(()) => out.close(true).map_err(StreamError::Io),
            // a failed disk shows in the content's error too; report the
            // disk's own error then
            Err(error) => Err(match out.close(false) {
                Err(disk) => StreamError::Io(disk),
                Ok(_) => StreamError::Content(error),
            }),
        };
        result.map_err(|error| match error {
            StreamError::Io(e) => StreamError::Io(with_cleanup(e, self.0.remove_file(path))),
            StreamError::Content(e) => {
                StreamError::Content(with_cleanup(e, self.0.remove_file(path)))
            }
        })
    }

    /// Start the new directory `name` for a bundle of files written one by
    /// one with [`Bundle::stream`].
    pub fn bundle(&self, name: &OsStr) -> io::Result<Bundle<'_>> {
        let path = single_entry(name)?;
        self.0.create_dir(path)?;
        match self.0.open_dir(path) {
            Ok(dir) => Ok(Bundle {
                parent: self,
                name: name.to_os_string(),
                dir: Self(dir),
                files: Vec::new(),
                done: false,
            }),
            Err(error) => Err(with_cleanup(error, self.0.remove_dir(path))),
        }
    }
}

/// A new directory of files written one at a time. Unless [`Bundle::finish`]
/// is called, the directory and its files are removed again: when a file
/// fails, at once (and its error says if the cleanup failed too), and when
/// the bundle is dropped, for errors that happen between files.
pub struct Bundle<'a> {
    parent: &'a MeshDirectory,
    name: std::ffi::OsString,
    dir: MeshDirectory,
    files: Vec<(String, String)>,
    done: bool,
}

impl Bundle<'_> {
    /// Stream the new file `name` into the bundle (see
    /// [`MeshDirectory::stream_new`]) and record its BLAKE3 hash.
    pub fn stream(
        &mut self,
        name: &str,
        content: impl FnOnce(&mut dyn Write) -> io::Result<()>,
    ) -> Result<(), StreamError> {
        match self.dir.stream_new(OsStr::new(name), content) {
            Ok(digest) => {
                self.files.push((name.to_string(), digest));
                Ok(())
            }
            Err(error) => {
                let cleanup = self.remove();
                self.done = true;
                Err(match error {
                    StreamError::Io(e) => StreamError::Io(with_cleanup(e, cleanup)),
                    StreamError::Content(e) => StreamError::Content(with_cleanup(e, cleanup)),
                })
            }
        }
    }

    /// The files written so far with their BLAKE3 hash, in order.
    pub fn files(&self) -> &[(String, String)] {
        &self.files
    }

    /// Keep the bundle; returns its files with their BLAKE3 hash.
    pub fn finish(mut self) -> Vec<(String, String)> {
        self.done = true;
        std::mem::take(&mut self.files)
    }

    fn remove(&mut self) -> io::Result<()> {
        let mut result = Ok(());
        for (name, _) in self.files.drain(..) {
            if let Err(e) = self.dir.0.remove_file(&name) {
                result = Err(e);
            }
        }
        let removed = single_entry(&self.name).and_then(|path| self.parent.0.remove_dir(path));
        result.and(removed)
    }
}

impl Drop for Bundle<'_> {
    fn drop(&mut self) {
        if !self.done {
            let _ = self.remove();
        }
    }
}

/// The bytes gathered before a block goes to the thread that hashes and
/// writes it.
pub const BLOCK: usize = 1 << 20;

/// A new file fed in blocks to a thread that hashes and writes them, so
/// producing, hashing and the disk overlap.
struct HashedFile {
    block: Vec<u8>,
    blocks: Option<SyncSender<Vec<u8>>>,
    spare: Receiver<Vec<u8>>,
    failed: Arc<AtomicBool>,
    worker: Option<JoinHandle<(String, io::Result<()>)>>,
}

impl HashedFile {
    fn new(mut file: File) -> Self {
        let (blocks, incoming) = mpsc::sync_channel::<Vec<u8>>(2);
        let (give_back, spare) = mpsc::channel();
        let failed = Arc::new(AtomicBool::new(false));
        let flag = Arc::clone(&failed);
        let worker = thread::spawn(move || {
            let mut hash = blake3::Hasher::new();
            let mut result = Ok(());
            for mut block in incoming {
                if result.is_ok() {
                    hash.update(&block);
                    result = file.write_all(&block);
                    if result.is_err() {
                        flag.store(true, Ordering::Relaxed);
                    }
                }
                block.clear();
                let _ = give_back.send(block);
            }
            if result.is_ok() {
                result = file.sync_all();
            }
            (hash.finalize().to_hex().to_string(), result)
        });
        Self {
            block: Vec::with_capacity(BLOCK),
            blocks: Some(blocks),
            spare,
            failed,
            worker: Some(worker),
        }
    }

    fn send(&mut self) -> io::Result<()> {
        if self.failed.load(Ordering::Relaxed) {
            return Err(io::Error::other("the file could not be written"));
        }
        let next = self
            .spare
            .try_recv()
            .unwrap_or_else(|_| Vec::with_capacity(BLOCK));
        let full = std::mem::replace(&mut self.block, next);
        match &self.blocks {
            Some(blocks) => blocks
                .send(full)
                .map_err(|_| io::Error::other("the thread writing the file stopped")),
            None => Err(io::Error::other("the file is closed")),
        }
    }

    /// Hand over the last block if `keep`, wait for the thread and return
    /// the digest, or the error of the file system.
    fn close(mut self, keep: bool) -> io::Result<String> {
        let rest = std::mem::take(&mut self.block);
        if let Some(blocks) = self.blocks.take()
            && keep
            && !rest.is_empty()
        {
            // a refused send means the thread is gone; its result says why
            let _ = blocks.send(rest);
        }
        match self.worker.take().map(JoinHandle::join) {
            Some(Ok((digest, Ok(())))) => Ok(digest),
            Some(Ok((_, Err(error)))) => Err(error),
            _ => Err(io::Error::other("the thread writing the file stopped")),
        }
    }
}

impl Write for HashedFile {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        self.write_all(data)?;
        Ok(data.len())
    }

    fn write_all(&mut self, mut data: &[u8]) -> io::Result<()> {
        while !data.is_empty() {
            let n = (BLOCK - self.block.len()).min(data.len());
            self.block.extend_from_slice(&data[..n]);
            data = &data[n..];
            if self.block.len() == BLOCK {
                self.send()?;
            }
        }
        Ok(())
    }

    /// Blocks go out when full and at the end; nothing to do here.
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl Drop for HashedFile {
    fn drop(&mut self) {
        self.blocks = None;
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}

/// The hash of every file digest KoflCraft writes: in `mesh-evidence.json`
/// (`file_hash`, format 1.1) and as the key of each entry of a tool
/// result's `files`. BLAKE3 hashes a large file several times faster than
/// SHA-256 and on many threads (the maintainer's choice, 2026-11-24:
/// SHA-256 of a 5 GB mesh took 17 s on a CPU without SHA instructions).
pub const FILE_HASH: &str = "blake3";

/// The BLAKE3 hash ([`FILE_HASH`]) of `bytes`, lower-case hex.
pub fn file_digest(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}

fn single_entry(name: &OsStr) -> io::Result<&Path> {
    let path = Path::new(name);
    let mut components = path.components();
    if let Some(Component::Normal(component)) = components.next()
        && component == name
        && components.next().is_none()
    {
        return Ok(path);
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        "expected a single filename within the authorized directory",
    ))
}

fn with_cleanup(error: io::Error, cleanup: io::Result<()>) -> io::Error {
    match cleanup {
        Ok(()) => error,
        Err(cleanup) => io::Error::new(error.kind(), format!("{error}. Cleanup failed: {cleanup}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Scratch(std::path::PathBuf);

    impl Scratch {
        fn new(name: &str) -> Self {
            let path =
                std::env::temp_dir().join(format!("kofl-stream-{}-{name}", std::process::id()));
            let _ = std::fs::remove_dir_all(&path);
            std::fs::create_dir_all(&path).unwrap();
            Self(path)
        }

        fn dir(&self) -> MeshDirectory {
            MeshDirectory::open(&self.0).unwrap()
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    /// Bytes that differ from block to block, over several blocks and a
    /// partial last one.
    fn content() -> Vec<u8> {
        (0..3 * BLOCK + 12_345)
            .map(|i| (i * 7 + i / BLOCK) as u8)
            .collect()
    }

    /// Invariant: a streamed file holds exactly the bytes written, and its
    /// digest is the BLAKE3 hash of those bytes, with writes of any size across
    /// block borders. Breaks if a block is lost, reordered or hashed
    /// differently from what reaches the disk, or the partial last block is
    /// dropped; the evidence digests of every bundle depend on it, and no
    /// other test writes more than one block.
    #[test]
    fn a_streamed_file_and_its_digest_match_the_bytes_written() {
        let scratch = Scratch::new("digest");
        let bytes = content();
        let digest = scratch
            .dir()
            .stream_new(OsStr::new("f.bin"), |out| {
                let mut rest = &bytes[..];
                let mut k = 1;
                while !rest.is_empty() {
                    let n = (k * 4099).min(rest.len());
                    out.write_all(&rest[..n])?;
                    rest = &rest[n..];
                    k = k * 3 % 1000 + 1;
                }
                Ok(())
            })
            .unwrap();
        assert_eq!(std::fs::read(scratch.0.join("f.bin")).unwrap(), bytes);
        assert_eq!(digest, file_digest(&bytes));
    }

    /// Invariant: a writer that fails leaves no file behind and is reported
    /// as a content error; a name that exists is a file-system error and
    /// the existing file is kept. Breaks if a half-written mesh file stays
    /// (a later reader would take it for a result, AG-08).
    #[test]
    fn a_failed_stream_removes_its_file_and_keeps_existing_ones() {
        let scratch = Scratch::new("fail");
        let dir = scratch.dir();
        let error = dir
            .stream_new(OsStr::new("f.bin"), |out| {
                out.write_all(&content())?;
                Err(io::Error::other("writer failed"))
            })
            .unwrap_err();
        assert!(matches!(error, StreamError::Content(_)), "{error}");
        assert!(!scratch.0.join("f.bin").exists());

        std::fs::write(scratch.0.join("old.bin"), b"old").unwrap();
        let error = dir
            .stream_new(OsStr::new("old.bin"), |out| out.write_all(b"new"))
            .unwrap_err();
        assert!(matches!(error, StreamError::Io(_)), "{error}");
        assert_eq!(std::fs::read(scratch.0.join("old.bin")).unwrap(), b"old");
    }

    /// Invariant: a bundle is either complete or gone: a failing file, or a
    /// bundle dropped without `finish`, removes the files written so far and
    /// the directory. Breaks if cleanup forgets the earlier files; the
    /// generators rely on it to never leave a partial bundle.
    #[test]
    fn a_bundle_is_complete_or_removed() {
        let scratch = Scratch::new("bundle");
        let dir = scratch.dir();
        let mut bundle = dir.bundle(OsStr::new("b")).unwrap();
        bundle.stream("a.txt", |out| out.write_all(b"a")).unwrap();
        let error = bundle
            .stream("b.txt", |_| Err(io::Error::other("writer failed")))
            .unwrap_err();
        assert!(matches!(error, StreamError::Content(_)));
        drop(bundle);
        assert!(!scratch.0.join("b").exists());

        let mut bundle = dir.bundle(OsStr::new("b")).unwrap();
        bundle.stream("a.txt", |out| out.write_all(b"a")).unwrap();
        drop(bundle);
        assert!(!scratch.0.join("b").exists());

        let mut bundle = dir.bundle(OsStr::new("b")).unwrap();
        bundle.stream("a.txt", |out| out.write_all(b"a")).unwrap();
        let files = bundle.finish();
        assert_eq!(files[0].1, file_digest(b"a"));
        assert_eq!(std::fs::read(scratch.0.join("b/a.txt")).unwrap(), b"a");
    }
}
