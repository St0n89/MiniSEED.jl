# MiniSEED.jl

Pure-Julia miniSEED 2.x reader and decoder. One file, zero packages — only
`Printf` and `Dates` from the Julia standard library. All binary parsing is
hand-rolled, all Steim decoding is a faithful transcription of the reference
algorithms in [libmseed](https://github.com/EarthScope/libmseed)'s
`unpackdata.c`.

Written in the same stdlib-only spirit as
[Newton_Fractals](https://github.com/St0n89/Newton_Fractals),
[INSAR](https://github.com/St0n89/INSAR) and
[bmp2png.jl](https://github.com/St0n89/bmp2png.jl).

Decodes data-only SEED (miniSEED v2) records with:

* **Steim-1** and **Steim-2** compression (the workhorses of seismic data
  distribution), big- and little-endian
* INT16, INT32, FLOAT32, FLOAT64, ASCII/TEXT
* fixed section of data header (FSDH) parsing, Blockette 1000/1001/100
* automatic record length and byte order detection (also for files without
  Blockette 1000)
* time correction handling, leap-second-safe microsecond timestamps

## Usage

```
julia miniseed.jl [options] FILE [FILE ...]
```

* (no flags) — one summary line per record
* `--verbose` — full record details (header, blockettes, flag bits)
* `--sample N` — also print the first N decoded samples per record
* `--channel GLOB` — only decode channels matching, e.g. `'*.?HZ'`
* `--csv FILE` — append a per-record summary CSV
* `--no-correct` — do not apply the time correction field
* `--conv MODE` — encoding convention: `auto` (default), `seed24`, `libmseed3`
* `--quiet` — suppress per-record warnings
* `--selftest` — run the built-in self test and exit

```
$ julia miniseed.jl reference-testdata-steim2.mseed2
000000 XX.TEST..BHZ  2012-05-12T00:00:00.000000  40 Hz  247 samples  STEIM2  [Xn ok]
000000 XX.TEST..BHZ  2012-05-12T00:00:06.175000  40 Hz  104 samples  STEIM2  [Xn ok]
...
```

## The encoding byte problem

SEED 2.4 and libmseed 3 disagree on Blockette 1000 encoding values, and both
conventions exist in real v2 files:

| byte | SEED 2.4 (ObsPy, IRIS, dataselect) | libmseed 3 internal |
|------|-----------------------------------|---------------------|
| 2    | INT32                             | —                   |
| 3    | FLOAT32                           | INT32               |
| 4    | FLOAT64                           | FLOAT32             |
| 5    | STEIM1                            | FLOAT64             |
| 10   | STEIM2                            | STEIM1              |
| 11   | GEOSCOPE (legacy)                 | STEIM2              |

libmseed 3 writes its internal values into v2 Blockette 1000 without
remapping (verified against its sources and its own reference files). In
`auto` mode this decoder tries SEED 2.4 first and falls back to the libmseed 3
meaning, validated by the Steim Xn integrity check, payload-size limits and
float plausibility (non-finite values, subnormal/normal mixes, absurd
magnitudes). `--conv seed24` or `--conv libmseed3` forces one convention.

## What's inside

* **FSDH parser** — 48-byte fixed header, station/network/channel/location,
  B1000/B1001/B1000 sample-rate override, activity/IO/data-quality flag bits,
  time correction (applied unless the "time applied" flag is set).
* **Record scanner** — sticky record length detection (B1000 exponent, header
  scan on a 64-byte grid, power-of-two remainder), byte order detection by
  year/day sanity, non-SEED gap skipping, truncated final records.
* **Steim decoders** — frame word layout (16 words, nibble word, X0/Xn),
  frame-0 first-difference skip, wrap-around 32-bit accumulation, Xn
  integrity warning, Steim1 (4×8/2×16/1×32) and Steim2 (4×8, 1×30, 2×15,
  3×10, 5×6, 6×5, 7×4) with file-byte-order sub-values.
* **Self test** — 146 checks: round-trips through a transcribed libmseed
  Steim encoder, synthetic records for every encoding and both byte orders,
  auto-convention fallbacks, and a deterministic 600-trial fuzz loop
  (xorshift64, 6 difference-magnitude classes) that re-encodes random sample
  streams and compares exact decode plus Xn.

## Validation

Verified against the official libmseed reference test data (packed by
libmseed's own `lm_pack` example, expected samples from `test/testdata.h`),
plus an independently written Python transcription of the libmseed decoders:

* 146/146 built-in self test checks.
* All 10 reference files decode with exact sample equality: steim1 BE/LE
  (500), steim2 BE/LE (499), int16 (220), int32 (500), float32 (500),
  float64 (500), text (234 chars), no-Blockette-1000 steim1 (2 records,
  7312 samples, Xn ok).
* Independent cross-check: `test/py_steim_check.py` (pure Python stdlib
  transcription of `msr_decode_steim1/2`) reproduces the Julia output
  sample-for-sample on all four Steim reference files.

## Test harness

```
julia miniseed.jl --selftest        # built-in self test (146 checks)
bash testdata/fetch.sh              # fetch official libmseed reference files
julia test/validate.jl              # exact-sample corpus validation
python3 test/py_steim_check.py       # independent Python cross-check
```

## References

* IRIS/EarthScope SEED manual, chapter 8 (Steim-1/Steim-2 compression)
* libmseed (Apache-2.0): https://github.com/EarthScope/libmseed

## License

MIT — see [LICENSE](LICENSE).
