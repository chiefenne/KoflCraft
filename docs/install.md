# Install and update

KoflCraft is built from source; release downloads are planned. This page
lists every tool to install on macOS, Linux and Windows, with the commands,
so that a fresh machine gets to a working `kofl` in one pass.

Contents: [What you need](#what-you-need) · [macOS](#macos) ·
[Linux](#linux) · [Windows](#windows) · [Build](#build) ·
[Optional: Gmsh](#optional-gmsh) ·
[Optional: Foehncraft](#optional-foehncraft) · [Update](#update)

## What you need

| Tool | Needed for | Required |
|---|---|---|
| Rust 1.88 or newer (rustup) | building KoflCraft | yes |
| C compiler and linker | linking the Rust programs | yes |
| Git | cloning the repository | yes |
| Gmsh program (`gmsh` on `PATH`) | the extra check `cargo m1 --require-gmsh` | no |
| Foehncraft's `foehn` | the round trip `cargo m1 --require-foehn` | no |

KoflCraft needs **no system libraries, no CMake and no Gmsh** to build,
mesh and export: it writes Gmsh `.msh`, VTU and the Foehncraft formats
with its own code. The Gmsh program is only used by the acceptance run to
open the exported `.msh` file once more (`gmsh cli.msh -check`).

CI builds and tests on Ubuntu (`ubuntu-latest`, with `apt install gmsh`);
macOS and Windows builds are not checked by CI yet.

## macOS

```
xcode-select --install                 # C compiler, linker and Git
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Open a new terminal, then `rustc --version` (1.88 or newer; else
`rustup update stable`). Continue with [Build](#build).

## Linux

```
# Ubuntu and Debian
sudo apt update
sudo apt install build-essential git curl
# Fedora, RHEL, Rocky and Alma
sudo dnf install gcc make git curl

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
rustc --version                        # 1.88 or newer; else: rustup update stable
```

Continue with [Build](#build).

## Windows

In PowerShell, with `winget` (part of Windows 10 and 11):

```
winget install --id Microsoft.VisualStudio.2022.BuildTools -e --override "--quiet --wait --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
winget install --id Rustlang.Rustup -e
winget install --id Git.Git -e
```

The first command installs the MSVC compiler and linker that Rust needs
(`rustup-init.exe` offers them too). Close PowerShell, open a new one, and
check `rustc --version`. With [Chocolatey](https://chocolatey.org):
`choco install rustup.install git visualstudio2022buildtools visualstudio2022-workload-vctools`.

## Build

On all three systems (PowerShell on Windows):

```
git clone https://github.com/chiefenne/KoflCraft.git
cd KoflCraft
cargo build --release --locked -p mesh-cli -p mesh-mcp
target/release/kofl generate test-fixtures/box-4.json box.mesh.json
target/release/kofl check box.mesh.json
target/release/kofl export box.mesh.json gmsh box.msh
```

This builds `kofl` (the command line; `kofl mcp` is the MCP server),
`mesh-cli` (the same commands under the older name) and `mesh-mcp` (the MCP
server on its own) into `target/release/` (`.exe` on Windows). To call
`kofl` from any folder, install it into `~/.cargo/bin`, which rustup puts
on your `PATH`:

```
cargo install --locked --path crates/mesh-cli
kofl commands
```

To build a release instead of `main`, check out its tag first (`git tag`
lists them, for example `git checkout v0.1.0`).

Tests: `cargo test --release --workspace`.

## Optional: Gmsh

Only for `cargo m1 --require-gmsh`, or to look at an exported `.msh` in the
Gmsh window.

| System | Command |
|---|---|
| macOS | `brew install gmsh` |
| Ubuntu, Debian | `sudo apt install gmsh` (what CI uses) |
| Fedora | `sudo dnf install gmsh` |
| Windows | download the Windows zip from <https://gmsh.info>, unpack it (for example to `C:\gmsh`) and add that folder to `PATH`: `setx PATH "$env:PATH;C:\gmsh"`, then open a new PowerShell |
| any, with conda | `conda install -c conda-forge gmsh` (puts `gmsh` on `PATH` in the active environment) |

On macOS and Linux, `pip install gmsh` in an active virtual environment
also gives a working `gmsh` command (a small Python launcher; on Linux it
needs `sudo apt install libglu1-mesa libxcursor1 libxft2 libxinerama1
libgomp1`). On Windows that launcher is not found as `gmsh`; use the zip.

Check and run:

```
gmsh -version
cargo m1 --require-gmsh
```

## Optional: Foehncraft

For `cargo m1 --require-foehn` (the export read back by Foehncraft), build
Foehncraft as in its [install guide](https://github.com/chiefenne/Foehncraft/blob/main/docs/install.md)
(it needs CMake in addition) and either put `foehn` on `PATH` or name it:

```
# macOS and Linux
FOEHN=$HOME/Foehncraft/target/release/foehn cargo m1 --require-gmsh --require-foehn
# Windows (PowerShell)
$env:FOEHN="$HOME\Foehncraft\target\release\foehn.exe"; cargo m1 --require-gmsh --require-foehn
```

## Update

```
git checkout main && git pull --ff-only
cargo build --release --locked -p mesh-cli -p mesh-mcp
```

If you installed `kofl` with `cargo install`, run that again.

## Problems on a first build

| Message | Cause and fix |
|---|---|
| `linker 'cc' not found` (Linux), `link.exe not found` (Windows) | No C toolchain: `build-essential` / `gcc`, or the Visual Studio C++ Build Tools. |
| `requires rustc 1.88` | `rustup update stable`. |
| `gmsh must be on PATH` | Only with `--require-gmsh`: install Gmsh as above, open a new terminal. |
