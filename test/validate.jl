#!/usr/bin/env julia
# =============================================================================
# test/validate.jl — full corpus validation for miniseed.jl
#
# Decodes every reference file in testdata/ and compares the samples exactly
# against the expected files in test/. The reference files are the official
# libmseed test data (fetch them first:  bash testdata/fetch.sh).
#
# Expected samples are derived from the dsinedata table in libmseed's
# test/testdata.h, packed by libmseed's own example program (lm_pack.c):
#   int32 / float32 / float64 / steim1: 500 samples
#   steim2: 499 (the last difference needs 31 bits)
#   int16: 220 (the first 220 values fit into 16 bits)
#   text: 234 ASCII characters
#
# Run from the repository root:  julia test/validate.jl
# Exit status: 0 if everything matches, 1 on any failure.
# =============================================================================

const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "miniseed.jl"))

function decoded(file)
    buf = read(file)
    out = Tuple{MS2Record,DecodeResult}[]
    for (rec, payload, _) in scan_ms2_records(buf)
        push!(out, (rec, decode_record(rec, payload)))
    end
    return out
end

function expect_int(path, ::Type{T}) where T <: Integer
    return T[parse(T, strip(l)) for l in eachline(path)]
end

function main()::Int
    ok = true
    td = joinpath(ROOT, "testdata")

    function check(label, got, expected; xn = true)
        pass = (got == expected) && xn
        ok = ok && pass
        println((pass ? "PASS" : "FAIL") * "  " * label *
                "  ($(length(got)) samples)")
    end

    # Steim-1 and Steim-2, big- and little-endian
    for f in ("reference-testdata-steim1.mseed2", "reference-testdata-steim1-LE.mseed2")
        s = Int32[]; xn = true
        for (rec, res) in decoded(joinpath(td, f))
            append!(s, res.samples)
            xn = xn && (res.xn_ok === nothing || res.xn_ok)
        end
        check(f, s, expect_int(joinpath(ROOT, "test/expected_steim1.txt"), Int32); xn)
    end
    for f in ("reference-testdata-steim2.mseed2", "reference-testdata-steim2-LE.mseed2")
        s = Int32[]; xn = true
        for (rec, res) in decoded(joinpath(td, f))
            append!(s, res.samples)
            xn = xn && (res.xn_ok === nothing || res.xn_ok)
        end
        check(f, s, expect_int(joinpath(ROOT, "test/expected_steim2.txt"), Int32); xn)
    end

    # INT16, INT32
    for (f, ef) in (("reference-testdata-int16.mseed2", "expected_int16.txt"),
                    ("reference-testdata-int32.mseed2", "expected_int32.txt"))
        s = Int32[]
        for (rec, res) in decoded(joinpath(td, f))
            append!(s, res.samples)
        end
        T = occursin("int16", f) ? Int32 : Int32
        check(f, s, expect_int(joinpath(ROOT, "test", ef), Int32))
    end

    # FLOAT64 (repr doubles), FLOAT32 (hex words, little-endian file bytes)
    s = Float64[]
    for (rec, res) in decoded(joinpath(td, "reference-testdata-float64.mseed2"))
        append!(s, res.samples)
    end
    check("reference-testdata-float64.mseed2", s,
          [parse(Float64, strip(l)) for l in eachline(joinpath(ROOT, "test/expected_float64.txt"))])

    s = Float32[]
    for (rec, res) in decoded(joinpath(td, "reference-testdata-float32.mseed2"))
        append!(s, res.samples)
    end
    expected32 = Float32[]
    for l in eachline(joinpath(ROOT, "test/expected_float32.hex"))
        push!(expected32, reinterpret(Float32, parse(UInt32, strip(l); base = 16)))
    end
    check("reference-testdata-float32.mseed2", s, expected32)

    # TEXT
    recs = decoded(joinpath(td, "reference-testdata-text.mseed2"))
    t = String(recs[1][2].samples)
    expected_text = read(joinpath(ROOT, "test/expected_text.txt"), String)
    check("reference-testdata-text.mseed2", [t], [expected_text])

    # no-Blockette-1000 file: Steim-1 assumed, Xn integrity check
    recs = decoded(joinpath(td, "testdata-no-blockette1000-steim1.mseed2"))
    n = 0; xn = true
    for (rec, res) in recs
        n += length(res.samples)
        xn = xn && (res.xn_ok === nothing || res.xn_ok)
    end
    pass = xn && n > 0
    ok = ok && pass
    println((pass ? "PASS" : "FAIL") * "  testdata-no-blockette1000-steim1.mseed2  " *
            "($(length(recs)) records, $n samples, Xn ok)")

    println(ok ? "ALL PASS" : "FAILURES")
    return ok ? 0 : 1
end

exit(main())
