This is a quick Chapel port of https://omlc.org/software/mc/small_mc.c.

### Different versions
- `small_mc.c`: the base version
- `small_mc.chpl`: more or less direct translation into Chapel
- `small_mc_v2.chpl`: code improvements to prepare for parallel execution
- `small_mc_v3.chpl`: first take for parallel implementation. This uses
  `atomic`s heavily for synchornization. It certainly kills performance. We need
  per-task arrays for the histogram-like pattern. That should be relatively easy
  implementation.
- `small_mc_v3_racy.chpl`: This removes `atomic`s causing racy execution. This
  is _not correct_ but represents the upper bound of parallel performance in
  this implementation. I think per-task arrays for histogram can come pretty
  close to this.
- `small_mc_v4.chpl`: This is the naive GPU implementation. The code uses
  array-of-structs representation for photons. This almost certainly makes the
  application memory bound. A next version should use struct of arrays. Note
  that this is a typical optimization for GPUs regardless of the programming
  language. This will reduce readability, but hopefully will improve performance
  significantly.
- `small_mc_v5: A less memory intensive version of the GPU implementation

### Quick performance experiment

Below are results from my workstation. While it is relatively idle, it is not an
ideal testbed.

- CPU: i5-11400 (6 cores)
- GPU: RTX A2000
- Chapel 2.2 (compiled with only `--fast`)
- gcc 13.2 (compiled with `-O3 -lm`)

All numbers are seconds (lower is better).

|numPhotons  | C     | v1     | v2     | v3 (parallel+naive sync)| v3 (parallel+race)| v4 (naive GPU) | v5 (less memory GPU)
|------------|-------|--------|--------|-------------------------|-------------------|----------------|---------------------
|100_000     | 0.625 | 0.276  | 0.320  | 1.550                   | 0.100             | 1.091          |
|1_000_000   | 6.131 | 2.351  | 2.715  | 15.047                  | 0.627             | 1.375          |
|10_000_000  | 60.98 | 23.434 | 24.014 | 46.931                  | 5.761             | 4.102          | 4.087
|100_000_000 |       |        |        |                         | 65.33             | OOM            | 31.087
