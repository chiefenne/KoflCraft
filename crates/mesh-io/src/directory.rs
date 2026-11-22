use std::ffi::OsStr;
use std::io::{self, Read, Write};
use std::path::{Component, Path};

use cap_std::ambient_authority;
use cap_std::fs::{Dir, OpenOptions};

/// Access explicitly granted to one working directory.
pub struct MeshDirectory(Dir);

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

    pub fn write_bundle(&self, name: &OsStr, files: &[(&str, &[u8])]) -> io::Result<()> {
        let path = single_entry(name)?;
        self.0.create_dir(path)?;
        let target = match self.0.open_dir(path) {
            Ok(dir) => Self(dir),
            Err(error) => return Err(with_cleanup(error, self.0.remove_dir(path))),
        };
        for (index, (name, bytes)) in files.iter().enumerate() {
            if let Err(mut error) = target.write_new(OsStr::new(name), bytes) {
                for (name, _) in &files[..index] {
                    error = with_cleanup(error, target.0.remove_file(name));
                }
                drop(target);
                return Err(with_cleanup(error, self.0.remove_dir(path)));
            }
        }
        Ok(())
    }
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
