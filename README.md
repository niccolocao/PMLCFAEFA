# PMLCFAEFA: Julia algorithms behind an R interface

This draft converts the three provided Julia modules into a package that can be
called **from R without installing Julia on the end user's machine**. The
source code is retained in `julia/source/`; the preparation script adapts the
module import paths and adds a restart limit to the generated compilation copy.

**Status: prototype, not yet compiled or numerically validated.** The original
algorithms contain a large number of model branches. The initial build targets
64-bit Linux; Windows and macOS would need separate compiled distributions
and testing. Do not publish or rely on this package for statistical results
until numerical equivalence has been checked against the Julia originals.

## Repository structure

- `julia/source/*.jl`: original files, kept unchanged
- `julia/PMLRBridge/src/PMLRBridge.jl`: C-callable Julia wrapper
- `scripts/prepare.jl`: copies the source into the Julia package, adapts module imports, and installs dependencies
- `scripts/build.jl`: compiles the Julia project to a relocatable library bundle
- `Rpkg`: R package with a registered `.Call` native interface
- `.github/workflows/build-linux.yml`: manual or tag-triggered Linux bundle build

## Building (on a Linux developer or CI machine)

Install **Julia 1.11**, **R**, a C compiler, and a usable compiler toolchain
first. Then execute from the repository root:

```bash
./scripts/build-linux.sh
R CMD INSTALL PMLCFAEFA_0.1.0.tar.gz
Rscript scripts/smoke.R
```

The first command installs Julia package dependencies, compiles the library
using PackageCompiler, and builds an R source tarball that **includes** the
Julia runtime. It may require substantial time, disk, and memory.

## Installation for an R user without Julia

Transfer `PMLCFAEFA_0.1.0.tar.gz`, generated above or downloaded from the
GitHub Actions build artifacts, to a compatible **64-bit Linux** machine.
Install it using:

```r
install.packages("PMLCFAEFA_0.1.0.tar.gz", repos = NULL, type = "source")
library(PMLCFAEFA)
PMLR_status()
```

The user does **not** need a separate Julia installation. Because this first
package release is installed from an R source archive and has a C extension,
the target machine still needs an R-compatible C toolchain. A future Linux
binary build would eliminate that requirement as well. The Julia runtime is
bundled, not eliminated.

## Example R call

```r
library(PMLCFAEFA)
set.seed(123)
X <- matrix(sample(1:4, 60 * 6, replace = TRUE), 60, 6)
ics <- cbind(intercept = 1, x1 = rnorm(60))

fit <- PMLcfaefa(
  X = X, q = 2, K = 2, c = 4, ics = ics,
  tol = 1e-6, maxit = 100, maxNumRestart = 3,
  replix = 3, fixthrehsolds = 0, latentmeans = 1
)
print(fit)
head(fit$raw)
```

This is an **API illustration**, not a validated example. The original
optimizer may fail to converge on random toy responses. The R function
returns the original flattened numeric result in `fit$raw`, plus a
`status_code` and a `restart_count` taken from the last two values.
Parameter extraction should be implemented only after each model branch has
been checked for consistent output indexing.

## Build through GitHub Actions

The repository lives at https://github.com/niccolocao/PMLCFAEFA and retains
its existing GPL license. Use **Actions → Build Julia library and R package
(Linux) → Run workflow**. After a successful run, download the
`PMLCFAEFA-linux-r-package` artifact and install the bundled R archive on a
compatible 64-bit Linux system. The bundle includes the Julia runtime and
does not require a separate Julia installation.

## Design limitations

- The R interface supports the main `PMLcfaefa` function, not every export from
  the three Julia modules.
- The optional `parametr` custom initial vector is not yet exposed.
- Calls are sequential within an R process. R and Julia share one process.
- The Julia runtime is initialized once and kept active for the R session.
- The bridge returns the raw flat vector; validated statistical post-processing
  is a separate step.
- Installation via GitHub source with `remotes::install_github()` **alone**
  does not work because the Julia bundle is generated during CI and excluded
  from version control. Install a built archive instead.
- Separate platform-specific builds are required for Windows, macOS and
  different processor architectures. Linux binary compatibility also depends
  on system-level libraries and build environment.
- Some original simulation helpers contain hardcoded absolute paths. The
  current bridge does not expose these helpers.

## Upstream code modifications during preparation

`PMLCFAEFA.jl` is preserved verbatim under `julia/source/`. The generated
module changes `using MixCFAEFA, Sgr, ...` to relative imports and replaces
`while code_error == 0.0` with an iteration-limited condition. When the
maximum restart count is reached without a nonzero status, a Julia exception
is returned to R instead of a potentially infinite loop.

## References

- PackageCompiler libraries: https://julialang.github.io/PackageCompiler.jl/stable/libs.html
- PackageCompiler API: https://julialang.github.io/PackageCompiler.jl/stable/refs.html
- Writing R Extensions: https://cran.r-project.org/doc/manuals/r-release/R-exts.html