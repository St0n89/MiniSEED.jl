#!/usr/bin/env julia
# =============================================================================
# miniseed.jl — Pure-Julia miniSEED 2.x reader and decoder
#
# Reads and decodes miniSEED v2 data records (data-only SEED):
#   Steim-1, Steim-2, INT16, INT32, FLOAT32, FLOAT64, ASCII/TEXT
#   with big- and little-endian headers and payloads.
#
# No external packages: only Printf and Dates from the Julia standard library.
# All binary parsing is hand-rolled, all Steim decoding is a faithful
# transcription of the reference algorithms in libmseed (unpackdata.c).
#
# Usage:
#   julia miniseed.jl [options] FILE [FILE ...]
#
#   (no flags)        One summary line per record
#   --info            Same as default
#   --verbose         Full record details (header, blockettes, flags)
#   --sample N        Also print the first N decoded samples per record
#   --channel GLOB    Only decode channels matching GLOB (* and ?),
#                     e.g. 'XX.TEST..BHZ' or '*.?HZ'
#   --csv FILE        Append a per-record summary CSV to FILE
#   --no-correct      Do not apply the time correction field
#   --conv MODE       Encoding convention for Blockette 1000 values:
#                     auto (default), seed24, libmseed3 (see note below)
#   --quiet           Suppress per-record warnings on stdout
#   --selftest        Run the built-in self test and exit
#   -h, --help        This text
#
# Encoding conventions (important):
#   SEED 2.4 and libmseed 3 disagree on Blockette 1000 encoding values.
#   Real-world v2 files (ObsPy, IRIS/DMC, dataselect, qlibmseed) follow
#   SEED 2.4: 2=INT32, 3=FLOAT32, 4=FLOAT64, 5=STEIM1, 10=STEIM2.
#   libmseed 3 writes its internal FDSN values into v2 Blockette 1000
#   (verified against its sources and reference files):
#   3=INT32, 4=FLOAT32, 5=FLOAT64, 10=STEIM1, 11=STEIM2.
#   In auto mode this decoder tries SEED 2.4 first and falls back to the
#   libmseed 3 meaning, validated by the Steim Xn integrity check and
#   payload-size plausibility:
#     0 -> TEXT        1 -> INT16       2 -> INT32
#     3 -> FLOAT32, fallback INT32 if the float interpretation is not
#          plausible (non-finite values or a subnormal/normal mix)
#     4 -> FLOAT64, fallback FLOAT32 if nsamp*8 does not fit the payload
#          or the float64 interpretation is implausible (non-finite,
#          subnormal, or absurd magnitudes)
#     5 -> STEIM1, fallback FLOAT64 on Steim error or Xn mismatch
#     10 -> STEIM2, fallback STEIM1 on Steim error or Xn mismatch
#     11 -> STEIM2 (SEED 11 is legacy GEOSCOPE, unsupported here)
#     (no B1000) -> STEIM1 assumed, as libmseed does
#   Use --conv seed24 or --conv libmseed3 to force one convention.
#
# Exit status: 0 on success, 1 on error or failed self test.
#
# Reference algorithms transcribed from libmseed
# (https://github.com/EarthScope/libmseed, Apache-2.0):
#   unpackdata.c msr_decode_steim1/2, msr_encode_steim1/2 (self test),
#   packdata.c BITWIDTH and packing order, ms_nomsamprate.
#
# Copyright (c) 2026 Dennis Kunde. MIT licensed, see LICENSE.
# =============================================================================

using Printf
using Dates

# =============================================================================
# Section 1: binary helpers (0-based offsets, 1-based buffer indexing)
# =============================================================================

u8(buf, o::Int) = buf[o + 1]

function u16(buf, o::Int, be::Bool)::UInt16
    b1 = UInt16(buf[o + 1]); b2 = UInt16(buf[o + 2])
    return be ? (b1 << 8) | b2 : b1 | (b2 << 8)
end

function u32(buf, o::Int, be::Bool)::UInt32
    if be
        return (UInt32(buf[o + 1]) << 24) | (UInt32(buf[o + 2]) << 16) |
               (UInt32(buf[o + 3]) << 8) | UInt32(buf[o + 4])
    else
        return UInt32(buf[o + 1]) | (UInt32(buf[o + 2]) << 8) |
               (UInt32(buf[o + 3]) << 16) | (UInt32(buf[o + 4]) << 24)
    end
end

function u64(buf, o::Int, be::Bool)::UInt64
    v = UInt64(0)
    if be
        for k in 0:7
            v |= UInt64(buf[o + 1 + k]) << (8 * (7 - k))
        end
    else
        for k in 0:7
            v |= UInt64(buf[o + 1 + k]) << (8 * k)
        end
    end
    return v
end

i16(buf, o::Int, be::Bool)::Int32 = Int32(reinterpret(Int16, u16(buf, o, be)))
i32(buf, o::Int, be::Bool)::Int32 = reinterpret(Int32, u32(buf, o, be))
f32(buf, o::Int, be::Bool)::Float32 = reinterpret(Float32, u32(buf, o, be))
f64(buf, o::Int, be::Bool)::Float64 = reinterpret(Float64, u64(buf, o, be))

"""Write a 32-bit word at 0-based offset o in the given byte order.
Julia's UInt8() conversion is checked, so every byte must be masked."""
function put32!(buf, o::Int, v::UInt32, be::Bool)
    if be
        buf[o + 1] = UInt8((v >> 24) & 0xFF)
        buf[o + 2] = UInt8((v >> 16) & 0xFF)
        buf[o + 3] = UInt8((v >> 8) & 0xFF)
        buf[o + 4] = UInt8(v & 0xFF)
    else
        buf[o + 1] = UInt8(v & 0xFF)
        buf[o + 2] = UInt8((v >> 8) & 0xFF)
        buf[o + 3] = UInt8((v >> 16) & 0xFF)
        buf[o + 4] = UInt8(v >> 24)
    end
    return buf
end

function put16!(buf, o::Int, v::UInt16, be::Bool)
    if be
        buf[o + 1] = UInt8((v >> 8) & 0xFF)
        buf[o + 2] = UInt8(v & 0xFF)
    else
        buf[o + 1] = UInt8(v & 0xFF)
        buf[o + 2] = UInt8(v >> 8)
    end
    return buf
end

"""Fixed-length string field: stop at NUL, then strip padding spaces."""
function fstr(buf, o::Int, len::Int)::String
    out = IOBuffer()
    for k in 0:(len - 1)
        b = buf[o + 1 + k]
        b == 0x00 && break
        print(out, Char(b))
    end
    return strip(String(take!(out)))
end

"""Wrap-around 32-bit addition/subtraction (Steim accumulation)."""
wadd(a::Int32, b::Int32)::Int32 =
    reinterpret(Int32, reinterpret(UInt32, a) + reinterpret(UInt32, b))
wsub(a::Int32, b::Int32)::Int32 =
    reinterpret(Int32, reinterpret(UInt32, a) - reinterpret(UInt32, b))

"""Sign-extend the low `bits` bits of v (UInt32) to Int32."""
function signext(v::UInt32, bits::Int)::Int32
    m = UInt32(1) << (bits - 1)          # sign bit
    x = Int64(v & (m - UInt32(1)))       # magnitude bits (Int64, no reinterpret)
    (v & m) != UInt32(0) && (x -= Int64(m))
    return Int32(x)
end

"""Low `bits` of a signed value as unsigned."""
bmask(v::Int32, bits::Int)::UInt32 = reinterpret(UInt32, v) & (UInt32(1) << bits - UInt32(1))

"""Sample rate from SEED rate factor/multiplier, as libmseed ms_nomsamprate."""
function nomsamprate(factor::Int, mult::Int)::Float64
    rate = 0.0
    if factor > 0
        rate = Float64(factor)
    elseif factor < 0
        rate = -1.0 / Float64(factor)
    end
    if mult > 0
        rate *= Float64(mult)
    elseif mult < 0
        rate = -rate / Float64(mult)
    end
    return rate
end

# =============================================================================
# Section 2: encoding names and conventions
# =============================================================================

const ENC_NAMES_SEED24 = Dict(
    0  => "ASCII",  1 => "INT16",   2 => "INT32", 3 => "FLOAT32",
    4  => "FLOAT64", 5 => "STEIM1", 10 => "STEIM2",
    11 => "GEOSCOPE16-3 (legacy)", 12 => "GEOSCOPE24 (legacy)",
    13 => "GEOSCOPE16-4 (legacy)", 16 => "CDSN16-4 (legacy)",
    30 => "SRO16-3 (legacy)", 31 => "SRO16-4 (legacy)",
    32 => "DWWSSN16-4 (legacy)")

const ENC_NAMES_LIBMSEED3 = Dict(
    0  => "TEXT", 1 => "INT16", 3 => "INT32", 4 => "FLOAT32", 5 => "FLOAT64",
    10 => "STEIM1", 11 => "STEIM2", 12 => "GEOSCOPE24 (legacy)",
    13 => "GEOSCOPE16-3 (legacy)", 14 => "GEOSCOPE16-4 (legacy)",
    16 => "CDSN16-4 (legacy)", 30 => "SRO16-3 (legacy)", 32 => "DWWSSN16-4 (legacy)")

encname(enc::Int, conv::Symbol)::String =
    (enc == -1 ? "none" :
     get(conv == :libmseed3 ? ENC_NAMES_LIBMSEED3 : ENC_NAMES_SEED24, enc,
         "unknown ($enc)"))

# =============================================================================
# Section 3: Steim decoders (transcriptions of libmseed unpackdata.c)
# =============================================================================

struct SteimError <: Exception
    msg::String
end

"""Decode Steim-2 compressed data. Returns (samples, xn_mismatch).
`data` must start at the record payload (data offset) and be >= 64*frames bytes.
Faithful transcription of msr_decode_steim2, including the first-frame
first-difference skip and the Xn integrity check (only meaningful when the
expected sample count is reached)."""
function decode_steim2(data::AbstractVector{UInt8}, nsamples::Int, be::Bool)
    maxframes = length(data) ÷ 64
    nsamples <= 0 || maxframes == 0 && return (Int32[], false)
    out = Vector{Int32}(undef, nsamples)
    diff = Vector{Int32}(undef, 105)   # max 15 words x 7 differences
    Xn = Int32(0)
    nout = 0
    for frameidx in 0:(maxframes - 1)
        (nout >= nsamples) && break
        fo = frameidx * 64
        w0 = u32(data, fo, be)         # nibble word, file byte order
        startw = frameidx == 0 ? 3 : 1
        nd = 0                         # diffs collected in this frame
        if frameidx == 0
            out[1] = i32(data, fo + 4, be)   # X0, output directly
            nout = 1
            Xn = i32(data, fo + 8, be)       # reverse integration constant
        end
        for widx in startw:15
            nibble = (w0 >> (30 - 2 * widx)) & UInt32(3)
            nibble == 0x0 && continue
            wo = fo + 4 * widx
            if nibble == 0x1           # 4 x 8-bit, file byte order
                for k in 0:3
                    b = (be ? (u32(data, wo, be) >> (24 - 8 * k)) :
                              (u32(data, wo, be) >> (8 * k))) & UInt32(0xFF)
                    nd += 1
                    diff[nd] = signext(b, 8)
                end
            else
                w = u32(data, wo, be)  # logical word for bit extraction
                dnib = (w >> 30) & UInt32(3)
                if nibble == 0x2       # 1x30, 2x15 or 3x10
                    if dnib == 0x1
                        nd += 1; diff[nd] = signext(w & UInt32(0x3FFFFFFF), 30)
                    elseif dnib == 0x2
                        nd += 1; diff[nd] = signext((w >> 15) & UInt32(0x7FFF), 15)
                        nd += 1; diff[nd] = signext(w & UInt32(0x7FFF), 15)
                    elseif dnib == 0x3
                        nd += 1; diff[nd] = signext((w >> 20) & UInt32(0x3FF), 10)
                        nd += 1; diff[nd] = signext((w >> 10) & UInt32(0x3FF), 10)
                        nd += 1; diff[nd] = signext(w & UInt32(0x3FF), 10)
                    else
                        throw(SteimError("impossible Steim2 dnib=00 for nibble=10"))
                    end
                else                    # nibble == 0x3: 5x6, 6x5 or 7x4
                    if dnib == 0x0
                        for k in 0:4
                            nd += 1
                            diff[nd] = signext((w >> (24 - 6 * k)) & UInt32(0x3F), 6)
                        end
                    elseif dnib == 0x1
                        for k in 0:5
                            nd += 1
                            diff[nd] = signext((w >> (25 - 5 * k)) & UInt32(0x1F), 5)
                        end
                    elseif dnib == 0x2
                        for k in 0:6
                            nd += 1
                            diff[nd] = signext((w >> (24 - 4 * k)) & UInt32(0xF), 4)
                        end
                    else
                        throw(SteimError("impossible Steim2 dnib=11 for nibble=11"))
                    end
                end
            end
        end
        # Apply differences; frame 0 skips its first difference
        di = frameidx == 0 ? 2 : 1
        while di <= nd && nout < nsamples
            out[nout + 1] = wadd(out[nout], diff[di])
            nout += 1
            di += 1
        end
    end
    xn_bad = (nout == nsamples && nsamples > 0 && out[nout] != Xn)
    return (nout == nsamples ? out : out[1:nout], xn_bad)
end

"""Decode Steim-1 compressed data. Returns (samples, xn_mismatch).
Transcription of msr_decode_steim1: 4x8-bit file byte order, 2x16-bit
file 16-bit order, 1x32-bit."""
function decode_steim1(data::AbstractVector{UInt8}, nsamples::Int, be::Bool)
    maxframes = length(data) ÷ 64
    nsamples <= 0 || maxframes == 0 && return (Int32[], false)
    out = Vector{Int32}(undef, nsamples)
    diff = Vector{Int32}(undef, 60)   # max 15 words x 4 differences
    Xn = Int32(0)
    nout = 0
    for frameidx in 0:(maxframes - 1)
        (nout >= nsamples) && break
        fo = frameidx * 64
        w0 = u32(data, fo, be)
        startw = frameidx == 0 ? 3 : 1
        nd = 0
        if frameidx == 0
            out[1] = i32(data, fo + 4, be)
            nout = 1
            Xn = i32(data, fo + 8, be)
        end
        for widx in startw:15
            nibble = (w0 >> (30 - 2 * widx)) & UInt32(3)
            nibble == 0x0 && continue
            wo = fo + 4 * widx
            w = u32(data, wo, be)
            if nibble == 0x1           # 4 x 8-bit, file byte order
                for k in 0:3
                    b = (be ? (w >> (24 - 8 * k)) : (w >> (8 * k))) & UInt32(0xFF)
                    nd += 1
                    diff[nd] = signext(b, 8)
                end
            elseif nibble == 0x2       # 2 x 16-bit, file 16-bit order
                nd += 1
                diff[nd] = Int32(reinterpret(Int16, be ? UInt16((w >> 16) & 0xFFFF) :
                                                        UInt16(w & 0xFFFF)))
                nd += 1
                diff[nd] = Int32(reinterpret(Int16, be ? UInt16(w & 0xFFFF) :
                                                        UInt16((w >> 16) & 0xFFFF)))
            else                        # 1 x 32-bit
                nd += 1
                diff[nd] = reinterpret(Int32, w)
            end
        end
        di = frameidx == 0 ? 2 : 1
        while di <= nd && nout < nsamples
            out[nout + 1] = wadd(out[nout], diff[di])
            nout += 1
            di += 1
        end
    end
    xn_bad = (nout == nsamples && nsamples > 0 && out[nout] != Xn)
    return (nout == nsamples ? out : out[1:nout], xn_bad)
end

# =============================================================================
# Section 4: simple payload decoders
# =============================================================================

function decode_int16(data, nsamples::Int, be::Bool)::Vector{Int32}
    n = min(nsamples, length(data) ÷ 2)
    return [i16(data, 2 * k, be) for k in 0:(n - 1)]
end

function decode_int32(data, nsamples::Int, be::Bool)::Vector{Int32}
    n = min(nsamples, length(data) ÷ 4)
    return [i32(data, 4 * k, be) for k in 0:(n - 1)]
end

function decode_float32(data, nsamples::Int, be::Bool)::Vector{Float32}
    n = min(nsamples, length(data) ÷ 4)
    return [f32(data, 4 * k, be) for k in 0:(n - 1)]
end

function decode_float64(data, nsamples::Int, be::Bool)::Vector{Float64}
    n = min(nsamples, length(data) ÷ 8)
    return [f64(data, 8 * k, be) for k in 0:(n - 1)]
end

function decode_text(data, nsamples::Int)::String
    n = min(nsamples, length(data))
    s = String(copy(@view data[1:n]))
    return rstrip(s, ['\0', ' ', '\n', '\r', '\t'])
end

"""Plausibility of interpreting raw words as float32: no non-finite values
and no mix of subnormal and normal magnitudes (integers misread as floats
produce subnormals)."""
function plausible_f32(vals::Vector{Float32})::Bool
    sub = false; big = false
    for v in vals
        !isfinite(v) && return false
        a = abs(v)
        a > 0 && a < Float32(1.17549435e-38) && (sub = true)
        a > Float32(1e-30) && (big = true)
    end
    return !(sub && big)
end

"""Plausibility of interpreting raw words as float64: no non-finite values
and no absurd magnitudes (float32 pairs misread as float64 produce values
around 1e60)."""
function plausible_f64(vals::Vector{Float64})::Bool
    for v in vals
        isfinite(v) || return false
        a = abs(v)
        (a > 0 && a < 1e-300) && return false   # subnormal
        a > 1e30 && return false
    end
    return true
end

# =============================================================================
# Section 5: fixed section of data header (FSDH) and blockettes
# =============================================================================

mutable struct MS2Record
    offset::Int               # byte offset in file
    reclen::Int               # detected record length (may exceed available bytes)
    truncated::Bool
    be::Bool                  # header byte order
    seq::String
    quality::UInt8
    station::String
    location::String
    channel::String
    network::String
    year::Int; day::Int; hour::Int; minute::Int; second::Int
    fsec::Int                 # 0.0001 s units
    nsamples::Int
    rate_factor::Int; rate_mult::Int
    actflags::UInt8; ioflags::UInt8; dqflags::UInt8
    numblockettes::Int
    timecorr::Int32           # 0.0001 s units, raw
    dataoffset::Int
    blocketteoffset::Int
    encoding::Int             # -1 if absent
    b1000_byteorder::Union{Int,Nothing}
    b1000_reclen::Union{Int,Nothing}
    b1001_tq::Union{Int,Nothing}
    b1001_usec::Union{Int,Nothing}
    b1001_framecount::Union{Int,Nothing}
    b100_rate::Union{Float32,Nothing}
    blockettes::Vector{Tuple{Int,Int}}  # (type, offset from record start)
    warnings::Vector{String}
end

MS2Record() = MS2Record(0, 0, false, true, "", 0x00, "", "", "", "",
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0x00, 0x00, 0x00,
                        0, Int32(0), 0, 0, -1, nothing, nothing,
                        nothing, nothing, nothing, nothing, Tuple{Int,Int}[],
                        String[])

"""Minimal check that bytes at `o` look like a v2 FSDH:
sequence number digits/spaces and a plausible year/day in either byte order."""
function is_ms2_header(buf, o::Int, len::Int)::Bool
    o + 48 > len && return false
    for k in 0:5
        b = buf[o + 1 + k]
        (b < 0x20 || b > 0x7e) && return false
    end
    _, ok = detect_byteorder(buf, o, len)
    return ok
end

"""Detect header byte order from year/day sanity. Returns (be, valid)."""
function detect_byteorder(buf, o::Int, len::Int)
    o + 24 > len && return (true, false)
    be = true
    year = Int(u16(buf, o + 20, be))
    day = Int(u16(buf, o + 22, be))
    valid(x, d) = 1900 <= x <= 2100 && 1 <= d <= 366
    if valid(year, day)
        return (true, true)
    end
    year = Int(u16(buf, o + 20, false))
    day = Int(u16(buf, o + 22, false))
    if valid(year, day)
        return (false, true)
    end
    return (true, false)
end

"""Parse a v2 record header at `o`. `avail` is how many bytes are actually
available (reclen or the truncated remainder of the file)."""
function parse_ms2_record(buf, o::Int, reclen::Int, avail::Int, be::Bool)::MS2Record
    r = MS2Record()
    r.offset = o
    r.reclen = reclen
    r.truncated = avail < reclen
    r.be = be
    r.seq = fstr(buf, o, 6)
    r.quality = u8(buf, o + 6)
    r.station = String(fstr(buf, o + 8, 5))
    r.location = String(fstr(buf, o + 13, 2))
    r.channel = String(fstr(buf, o + 15, 3))
    r.network = String(fstr(buf, o + 18, 2))
    r.year = Int(u16(buf, o + 20, be))
    r.day = Int(u16(buf, o + 22, be))
    r.hour = Int(u8(buf, o + 24))
    r.minute = Int(u8(buf, o + 25))
    r.second = Int(u8(buf, o + 26))
    r.fsec = Int(u16(buf, o + 28, be))
    r.nsamples = Int(u16(buf, o + 30, be))
    r.rate_factor = Int(reinterpret(Int16, u16(buf, o + 32, be)))
    r.rate_mult = Int(reinterpret(Int16, u16(buf, o + 34, be)))
    r.actflags = u8(buf, o + 36)
    r.ioflags = u8(buf, o + 37)
    r.dqflags = u8(buf, o + 38)
    r.numblockettes = Int(u8(buf, o + 39))
    r.timecorr = i32(buf, o + 40, be)
    r.dataoffset = Int(u16(buf, o + 44, be))
    r.blocketteoffset = Int(u16(buf, o + 46, be))

    if !(1900 <= r.year <= 2100)
        push!(r.warnings, "implausible year $(r.year)")
    end
    if !(1 <= r.day <= 366)
        push!(r.warnings, "implausible day-of-year $(r.day)")
    end
    if r.hour > 23 || r.minute > 59 || r.second > 60
        push!(r.warnings, "implausible time $(r.hour):$(r.minute):$(r.second)")
    end

    # Walk the blockette chain
    boff = r.blocketteoffset
    seen = 0
    if boff != 0
        while boff != 0 && seen < 64
            if boff < 48 || boff + 4 > avail
                push!(r.warnings, "blockette chain broken at offset $boff")
                break
            end
            btype = Int(u16(buf, o + boff, be))
            next = Int(u16(buf, o + boff + 2, be))
            push!(r.blockettes, (btype, boff))
            if btype == 1000
                if boff + 8 <= avail
                    r.encoding = Int(u8(buf, o + boff + 4))
                    bo_ = Int(u8(buf, o + boff + 5))
                    exp_ = Int(u8(buf, o + boff + 6))
                    if bo_ == 0 || bo_ == 1
                        r.b1000_byteorder = bo_
                    else
                        push!(r.warnings, "blockette 1000 byte order flag invalid: $bo_")
                    end
                    if 0 <= exp_ <= 30
                        r.b1000_reclen = 2^exp_
                        if !r.truncated && 2^exp_ != reclen
                            push!(r.warnings,
                                   "blockette 1000 record length $(2^exp_) != detected $reclen")
                        end
                    else
                        push!(r.warnings, "blockette 1000 record length exponent invalid: $exp_")
                    end
                else
                    push!(r.warnings, "blockette 1000 truncated")
                end
            elseif btype == 1001
                if boff + 8 <= avail
                    r.b1001_tq = Int(u8(buf, o + boff + 4))
                    r.b1001_usec = Int(reinterpret(Int8, u8(buf, o + boff + 5)))
                    r.b1001_framecount = Int(u8(buf, o + boff + 7))
                end
            elseif btype == 100
                if boff + 8 <= avail
                    r.b100_rate = f32(buf, o + boff + 4, be)
                end
            end
            (next != 0 && (next <= boff)) &&
                push!(r.warnings, "blockette chain not increasing at offset $boff")
            boff = next
            seen += 1
        end
    end
    if r.dataoffset != 0 && (r.dataoffset < 48 || r.dataoffset > reclen - 1) &&
       !r.truncated
        push!(r.warnings, "data offset $((r.dataoffset)) outside record")
    end
    return r
end

"""Detect record length for the record at `o`, as libmseed does:
blockette 1000 exponent first, then scan a 64-byte grid for the next
header, then a power-of-two remainder as last resort."""
function detect_reclen(buf, o::Int, len::Int)::Int
    # Try blockette 1000
    if is_ms2_header(buf, o, len)
        be, ok = detect_byteorder(buf, o, len)
        if ok
            boff = Int(u16(buf, o + 46, be))
            seen = 0
            while boff != 0 && seen < 64 && boff >= 48 && boff + 8 <= len - o
                btype = Int(u16(buf, o + boff, be))
                if btype == 1000
                    exp_ = Int(u8(buf, o + boff + 6))
                    0 <= exp_ <= 30 && return 2^exp_
                    break
                end
                boff = Int(u16(buf, o + boff + 2, be))
                seen += 1
            end
        end
    end
    # Scan 64-byte grid for the next header
    maxscan = min(len - o, 2^20)
    for r in 64:64:maxscan
        if is_ms2_header(buf, o + r, len)
            return r
        end
    end
    # Power-of-two remainder
    rem = len - o
    if rem >= 128 && rem <= 2^20 && (rem & (rem - 1)) == 0
        return rem
    end
    return 0   # undetectable
end

const MINRECLEN = 128
const MAXRECLEN = 2^20

"""Scan a buffer for v2 records. Returns Vector{Tuple{MS2Record,SubArray}}
where the SubArray is the record payload (data offset .. record end).
Skips non-SEED bytes on the 64-byte grid, keeps a sticky record length."""
function scan_ms2_records(buf::Vector{UInt8})
    results = Tuple{MS2Record,SubArray{UInt8,1,Vector{UInt8}},Vector{UInt8}}[]
    len = length(buf)
    o = 0
    reclen = 0
    while o + 48 <= len
        if !is_ms2_header(buf, o, len)
            o += 64
            continue
        end
        be, ok = detect_byteorder(buf, o, len)
        ok || (o += 64; continue)
        # Sticky record length: reuse if the next record fits the grid
        if reclen == 0 || o + reclen > len || !is_ms2_header(buf, o, len) == false
            r2 = detect_reclen(buf, o, len)
            if r2 >= MINRECLEN || (r2 > 0 && o + r2 == len)
                reclen = r2
            end
        end
        rl = reclen
        if rl == 0
            o += 64
            continue
        end
        avail = min(rl, len - o)
        rec = parse_ms2_record(buf, o, rl, avail, be)
        if avail < rl
            rec.truncated = true
            push!(rec.warnings, "record truncated: $avail of $rl bytes")
        end
        doff = rec.dataoffset
        pstart = clamp(doff, 48, avail)
        payload = @view buf[(o + pstart + 1):(o + avail)]
        push!(results, (rec, payload, buf))
        o += avail
    end
    return results
end

# =============================================================================
# Section 6: channel identification, sample rate, time
# =============================================================================

channelid(r::MS2Record) = "$(r.network).$(r.station).$(r.location).$(r.channel)"

samprate(r::MS2Record)::Float64 =
    r.b100_rate !== nothing ? Float64(r.b100_rate) :
    nomsamprate(r.rate_factor, r.rate_mult)

"""Start time in microseconds since the Unix epoch.
Applies the time correction (0.0001 s units) when `corrected` is true and
the "time correction applied" activity flag (0x02) is clear, as SEED defines.
Second 60 (leap second) rolls monotonically into the next second."""
function start_us(r::MS2Record; corrected::Bool=true)::Int64
    sec = min(r.second, 59)
    dt = DateTime(r.year, 1, 1) + Day(r.day - 1) + Hour(r.hour) +
         Minute(r.minute) + Second(sec)
    us = (r.second == 60 ? 1_000_000 : 0) + r.fsec * 100
    if r.b1001_usec !== nothing
        us += r.b1001_usec
    end
    if corrected && (r.actflags & 0x02) == 0x00
        us += Int64(r.timecorr) * 100
    end
    total = Int64(Dates.value(dt - DateTime(1970, 1, 1))) * 1000 + us
    return total
end

"""ISO 8601 with microseconds: yyyy-mm-ddTHH:MM:SS.ffffff"""
function iso_us(total_us::Int64)::String
    sec, remus = fldmod(total_us, 1_000_000)
    dt = DateTime(1970, 1, 1) + Dates.Second(sec)
    return @sprintf("%s.%06d", Dates.format(dt, "yyyy-mm-ddTHH:MM:SS"), remus)
end

iso_us(r::MS2Record; corrected::Bool=true) = iso_us(start_us(r; corrected))

"""Payload byte order: blockette 1000 byte order, else header byte order."""
payload_be(r::MS2Record)::Bool = r.b1000_byteorder === nothing ? r.be : (r.b1000_byteorder == 1)

# =============================================================================
# Section 7: record decoding with dual encoding convention
# =============================================================================

struct DecodeResult
    samples::Union{Vector{Int32},Vector{Float32},Vector{Float64},String,Nothing}
    nsamples::Int             # samples decoded
    encoding_used::String
    xn_ok::Union{Bool,Nothing} # Steim integrity check (nothing if n/a)
    warnings::Vector{String}
end

const CONVS = ("auto", "seed24", "libmseed3")

"""Try Steim `kind` (1 or 2). Success requires the full expected sample
count to be decoded and the Xn integrity check to pass."""
function try_steim(payload, nsamp::Int, be::Bool, kind::Int)
    try
        s, xnbad = (kind == 2 ? decode_steim2(payload, nsamp, be) :
                               decode_steim1(payload, nsamp, be))
        if length(s) < nsamp
            return (s, false, "only $(length(s)) of $nsamp samples decodable")
        end
        return (s, true, xnbad ? "Xn mismatch" : "")
    catch e
        e isa SteimError || rethrow()
        return (Int32[], false, e.msg)
    end
end

"""Decode one record payload following the encoding convention.
conv in (:auto, :seed24, :libmseed3)."""
function decode_record(r::MS2Record, payload::AbstractVector{UInt8};
                       conv::Symbol=:auto)
    warn = copy(r.warnings)
    enc = r.encoding
    be = payload_be(r)
    nsamp = r.nsamples
    local samples, used
    xn_ok = nothing

    function short(n::Int, why::String)
        if length(payload) ÷ (n == 0 ? 1 : n) < nsamp
            push!(warn, why)
        end
    end

    if enc == -1
        push!(warn, "no blockette 1000, assuming Steim-1")
        s, ok, err = try_steim(payload, nsamp, be, 1)
        samples = s; used = "STEIM1"; xn_ok = ok && err == ""
        err != "" && push!(warn, "Steim-1 decode failed: $err")
    elseif enc == 0
        samples = decode_text(payload, nsamp); used = "TEXT"
    elseif enc == 1
        samples = decode_int16(payload, nsamp, be); used = "INT16"
    elseif enc == 2
        samples = decode_int32(payload, nsamp, be); used = "INT32"
    elseif enc == 3
        f = decode_float32(payload, nsamp, be)
        i = decode_int32(payload, nsamp, be)
        if conv == :libmseed3
            samples = i; used = "INT32"
        elseif conv == :seed24 || plausible_f32(f)
            samples = f; used = "FLOAT32"
        else
            push!(warn, "encoding 3: float interpretation implausible, " *
                        "using INT32 (libmseed 3 convention)")
            samples = i; used = "INT32"
        end
    elseif enc == 4
        fits64 = nsamp * 8 <= length(payload)
        if conv == :libmseed3
            samples = decode_float32(payload, nsamp, be); used = "FLOAT32"
        elseif conv == :seed24
            samples = decode_float64(payload, nsamp, be); used = "FLOAT64"
        elseif fits64
            d = decode_float64(payload, nsamp, be)
            if plausible_f64(d)
                samples = d; used = "FLOAT64"
            else
                push!(warn, "encoding 4: float64 interpretation implausible, " *
                            "using FLOAT32 (libmseed 3 convention)")
                samples = decode_float32(payload, nsamp, be); used = "FLOAT32"
            end
        else
            push!(warn, "encoding 4: data does not fit as FLOAT64, " *
                        "using FLOAT32 (libmseed 3 convention)")
            samples = decode_float32(payload, nsamp, be); used = "FLOAT32"
        end
    elseif enc == 5
        fits64 = nsamp * 8 <= length(payload)
        steimok = false
        if conv != :libmseed3
            s, ok, err = try_steim(payload, nsamp, be, 1)
            if ok && err == ""
                samples = s; used = "STEIM1"; xn_ok = true; steimok = true
            end
        end
        if !steimok
            if conv == :seed24
                push!(warn, "Steim-1 decode failed: $err")
                samples = Int32[]; used = "STEIM1"; xn_ok = false
            elseif fits64
                push!(warn, "encoding 5: Steim-1 failed, " *
                            "using FLOAT64 (libmseed 3 convention)")
                samples = decode_float64(payload, nsamp, be); used = "FLOAT64"
            else
                push!(warn, "Steim-1 decode failed: $err")
                samples = Int32[]; used = "STEIM1"; xn_ok = false
            end
        end
    elseif enc == 10
        if conv == :libmseed3
            s, ok, err = try_steim(payload, nsamp, be, 1)
            samples = s; used = "STEIM1"; xn_ok = ok && err == ""
            err != "" && push!(warn, "Steim-1 decode failed: $err")
        elseif conv == :seed24
            s, ok, err = try_steim(payload, nsamp, be, 2)
            samples = s; used = "STEIM2"; xn_ok = ok && err == ""
            err != "" && push!(warn, "Steim-2 decode failed: $err")
        else
            s, ok, err = try_steim(payload, nsamp, be, 2)
            if ok && err == ""
                samples = s; used = "STEIM2"; xn_ok = true
            else
                s2, xnbad = decode_steim1(payload, nsamp, be)
                samples = s2; used = "STEIM1"; xn_ok = !xnbad
                push!(warn, "encoding 10: Steim-2 failed ($err), " *
                            "using Steim-1 (libmseed 3 convention)")
            end
        end
    elseif enc == 11
        if conv == :seed24
            push!(warn, "encoding 11 (SEED: GEOSCOPE 16-3) not supported")
            samples = nothing; used = "GEOSCOPE16-3"
        else
            s, ok, err = try_steim(payload, nsamp, be, 2)
            samples = s; used = "STEIM2"; xn_ok = ok && err == ""
            err != "" && push!(warn, "Steim-2 decode failed: $err")
        end
    else
        push!(warn, "unsupported encoding $enc")
        samples = nothing
        used = encname(enc, :seed24)
    end
    if samples isa AbstractVector && length(samples) < nsamp
        push!(warn, "decoded $(length(samples)) of $nsamp samples")
    end
    return DecodeResult(samples, samples isa AbstractVector ? length(samples) : 0,
                        used, xn_ok, warn)
end

# =============================================================================
# Section 8: command line interface
# =============================================================================

function glob_to_re(g::String)
    out = ""
    for c in g
        if c == '*'
            out *= ".*"
        elseif c == '?'
            out *= "."
        elseif c == '.'
            out *= "\\."
        else
            out *= c
        end
    end
    return Regex("^" * out * "\$")
end

function print_flag_bits(byte::UInt8, names::Vector{String})
    on = String[]
    for k in 1:length(names)
        (byte & (UInt8(1) << (k - 1))) != 0 && push!(on, names[k])
    end
    return on
end

const ACT_FLAG_NAMES = ["calibration", "time correction applied", "event begin",
                        "event end", "event in progress", "negative leap second",
                        "positive leap second"]
const IO_FLAG_NAMES = ["station volume parity error", "long record", "short record",
                       "start of time series", "end of time series", "clock locked"]
const DQ_FLAG_NAMES = ["amplifier saturation", "digitizer clipping", "spikes",
                       "glitches", "missing/padded data", "telemetry sync error",
                       "filter charging", "time tag questionable"]

function print_record(rec::MS2Record, res::DecodeResult, verbose::Bool,
                      nsample::Int, corrected::Bool)
    println("$(rpad(rec.seq, 6)) $(rpad(channelid(rec), 20)) " *
            "$(rpad(iso_us(rec; corrected), 27)) " *
            "$(rpad(@sprintf("%.6g", samprate(rec)), 10)) Hz  " *
            "$(rpad(rec.nsamples, 5)) samples  $(res.encoding_used)" *
            (res.xn_ok === nothing ? "" : (res.xn_ok ? "  [Xn ok]" : "  [Xn MISMATCH]")))
    if !verbose
        return
    end
    println("    header byte order:    $(rec.be ? "big" : "little") endian")
    println("    payload byte order:   $(payload_be(rec) ? "big" : "little") endian")
    println("    record length:        $(rec.reclen)$(rec.truncated ? " (truncated)" : "")")
    println("    FSDH quality/flags:   $(Char(rec.quality))  act=$(rec.actflags) io=$(rec.ioflags) dq=$(rec.dqflags)")
    act = print_flag_bits(rec.actflags, ACT_FLAG_NAMES)
    io_ = print_flag_bits(rec.ioflags, IO_FLAG_NAMES)
    dq = print_flag_bits(rec.dqflags, DQ_FLAG_NAMES)
    println("    activity flags:       $(isempty(act) ? "-" : join(act, ", "))")
    println("    I/O flags:            $(isempty(io_) ? "-" : join(io_, ", "))")
    println("    data quality flags:   $(isempty(dq) ? "-" : join(dq, ", "))")
    println("    rate factor/mult:     $(rec.rate_factor)/$(rec.rate_mult)" *
            (rec.b100_rate !== nothing ? "  (blockette 100: $(rec.b100_rate) Hz)" : ""))
    println("    encoding byte:        $(rec.encoding)$(rec.encoding >= 0 ? "" : " (none)")" *
            "  [SEED 2.4: $(encname(rec.encoding, :seed24)), libmseed 3: $(encname(rec.encoding, :libmseed3))]")
    if rec.b1001_tq !== nothing
        println("    timing quality:        $(rec.b1001_tq)%  usec=$(rec.b1001_usec)  framecount=$(rec.b1001_framecount)")
    end
    if !isempty(rec.blockettes)
        bs = join(["B$(t)@$(o)" for (t, o) in rec.blockettes], ", ")
        println("    blockettes:            $bs")
    end
    println("    time correction:       $(Int(rec.timecorr) / 10000.0) s" *
            ((rec.actflags & 0x02) != 0 ? " (already applied)" : ""))
    if nsample > 0 && res.samples isa AbstractVector && !isempty(res.samples)
        n = min(nsample, length(res.samples))
        println("    first $n samples:      $(join(res.samples[1:n], ", "))")
    end
    if !isempty(res.warnings)
        for w in res.warnings
            println("    warning: $w")
        end
    end
end

function usage()
    println("""
    miniseed.jl — pure-Julia miniSEED 2.x reader/decoder

    usage: julia miniseed.jl [options] FILE [FILE ...]

      (no flags)      one summary line per record
      --info          same as default
      --verbose      full record details
      --sample N     print first N decoded samples per record
      --channel GLOB filter channels, e.g. 'XX.TEST..BHZ' or '*.?HZ'
      --csv FILE     append per-record summary CSV to FILE
      --no-correct   do not apply the time correction field
      --conv MODE    encoding convention: auto (default), seed24, libmseed3
      --quiet        suppress per-record warnings
      --selftest     run the built-in self test and exit
      -h, --help     this text""")
end

function main(args)::Int
    verbose = false
    quiet = false
    selftest = false
    nsample = 0
    corrected = true
    conv = :auto
    channelpat = nothing
    csvfile = nothing
    files = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--selftest"
            selftest = true
        elseif a == "--verbose" || a == "-v"
            verbose = true
        elseif a == "--info"
            nothing
        elseif a == "--quiet" || a == "-q"
            quiet = true
        elseif a == "--no-correct"
            corrected = false
        elseif a == "--conv"
            i += 1
            (i > length(args) || !(args[i] in CONVS)) &&
                (println(stderr, "error: --conv needs one of: $(join(CONVS, ", "))"); return 1)
            conv = Symbol(args[i])
        elseif a == "--sample"
            i += 1
            (i > length(args) || tryparse(Int, args[i]) === nothing) &&
                (println(stderr, "error: --sample needs a number"); return 1)
            nsample = parse(Int, args[i])
        elseif a == "--channel"
            i += 1
            i > length(args) && (println(stderr, "error: --channel needs a pattern"); return 1)
            channelpat = glob_to_re(args[i])
        elseif a == "--csv"
            i += 1
            i > length(args) && (println(stderr, "error: --csv needs a file name"); return 1)
            csvfile = args[i]
        elseif a == "-h" || a == "--help"
            usage(); return 0
        elseif startswith(a, "-") && a != "-"
            println(stderr, "error: unknown option $a"); usage(); return 1
        else
            push!(files, a)
        end
        i += 1
    end
    if selftest
        return run_selftest() ? 0 : 1
    end
    isempty(files) && (usage(); return 1)

    csv_io = nothing
    if csvfile !== nothing
        csv_io = open(csvfile, append=true)
        println(csv_io, "file,record,seq,channel,start_iso,samprate_hz,nsamples," *
                       "encoding_byte,encoding_used,decoded,xn_ok,warnings")
    end
    status = 0
    recidx = 0
    for f in files
        isfile(f) || (println(stderr, "error: cannot read $f"); status = 1; continue)
        buf = read(f)
        records = scan_ms2_records(buf)
        if isempty(records)
            println(stderr, "warning: $f: no miniSEED records found")
        end
        for (rec, payload, _) in records
            recidx += 1
            if channelpat !== nothing && !occursin(channelpat, channelid(rec))
                continue
            end
            res = decode_record(rec, payload; conv=conv)
            print_record(rec, res, verbose, nsample, corrected)
            if !quiet
                for w in res.warnings
                    println(stderr, "warning: $(basename(f)): $(channelid(rec)): $w")
                end
            end
            if csv_io !== nothing
                xn = res.xn_ok === nothing ? "" : (res.xn_ok ? "ok" : "mismatch")
                row = (basename(f), recidx, rec.seq, channelid(rec),
                       iso_us(rec; corrected), samprate(rec), rec.nsamples,
                       rec.encoding, res.encoding_used, res.nsamples, xn,
                       length(res.warnings))
                println(csv_io, join([string(x) for x in row], ","))
            end
        end
    end
    csv_io !== nothing && close(csv_io)
    return status
end

# =============================================================================
# Section 9: self test
#   - unit tests for all binary helpers and both Steim decoders
#   - Steim encoders transcribed from libmseed packdata.c, used for
#     round-trip and fuzz tests (decoder against independent encoder)
#   - synthetic records for header parsing, time handling and the
#     dual encoding-convention logic
# =============================================================================

const TC = Ref(0)      # total checks
const TF = Ref(0)      # failed checks

function check(cond::Bool, name::String)
    TC[] += 1
    if cond
        @printf("ok %2d - %s\n", TC[], name)
    else
        TF[] += 1
        @printf("not ok %2d - %s\n", TC[], name)
    end
end

function checkeq(a, b, name::String)
    check(a == b, "$name (got $(repr(a)), expected $(repr(b)))")
end

"""libmseed BITWIDTH: minimal bits from (4,5,6,8,10,15,16,30,32) holding v."""
function bitwidth(v::Int32)::Int
    for w in (4, 5, 6, 8, 10, 15, 16, 30, 32)
        if -(Int64(1) << (w - 1)) <= Int64(v) <= (Int64(1) << (w - 1)) - 1
            return w
        end
    end
    return 64
end

"""Encode samples as Steim-2 frames (transcription of msr_encode_steim2).
diff0 is the (skipped) difference before the first sample. Returns a byte
vector of whole frames, or nothing if the data is not representable
(differences needing more than 30 bits) or does not fit maxframes."""
function encode_steim2(samples::Vector{Int32}, diff0::Int32, maxframes::Int, be::Bool)
    nsamp = length(samples)
    nsamp == 0 && return UInt8[]
    maxframes == 0 && return nothing
    frames = Vector{UInt8}(undef, 64 * maxframes)
    fill!(frames, 0x00)
    # Difference queue: bitwidth must be representable to pack
    bw0 = bitwidth(diff0)
    bw0 > 30 && return nothing
    queue = Int32[diff0]
    bws = Int[bw0]
    qi = 1                       # next sample diff to push (samples[qi+1]-samples[qi])
    outputsamples = 0
    for frameidx in 0:(maxframes - 1)
        outputsamples >= nsamp && break
        fo = 64 * frameidx
        nibbles = zeros(UInt32, 16)
        startw = frameidx == 0 ? 3 : 1
        if frameidx == 0
            put32!(frames, fo + 4, reinterpret(UInt32, samples[1]), be)
        end
        for widx in startw:15
            outputsamples >= nsamp && break
            # Refill queue to 7
            while length(queue) < 7 && qi <= nsamp - 1
                d = wsub(samples[qi + 1], samples[qi])
                bw = bitwidth(d)
                bw > 30 && return nothing
                push!(queue, d); push!(bws, bw)
                qi += 1
            end
            isempty(queue) && break
            local n, wv
            # Greedy packing order, as libmseed
            if length(queue) == 7 && all(bws[1:7] .<= 4)
                n = 7
                wv = UInt32(0x2) << 30                   # dnib=10: 7x4
                for k in 0:6
                    wv |= bmask(queue[k + 1], 4) << (24 - 4 * k)
                end
                nibbles[widx + 1] = 0x3
            elseif length(queue) >= 6 && all(bws[1:6] .<= 5)
                n = 6
                wv = UInt32(0x1) << 30                   # dnib=01: 6x5
                for k in 0:5
                    wv |= bmask(queue[k + 1], 5) << (25 - 5 * k)
                end
                nibbles[widx + 1] = 0x3
            elseif length(queue) >= 5 && all(bws[1:5] .<= 6)
                n = 5
                wv = UInt32(0x0)                          # dnib=00: 5x6
                for k in 0:4
                    wv |= bmask(queue[k + 1], 6) << (24 - 6 * k)
                end
                nibbles[widx + 1] = 0x3
            elseif length(queue) >= 4 && all(bws[1:4] .<= 8)
                n = 4
                wv = UInt32(0)                            # nibble 01: 4x8, file bytes
                for k in 0:3
                    sh = be ? (24 - 8 * k) : (8 * k)
                    wv |= bmask(queue[k + 1], 8) << sh
                end
                nibbles[widx + 1] = 0x1
            elseif length(queue) >= 3 && all(bws[1:3] .<= 10)
                n = 3
                wv = UInt32(0x3) << 30                    # dnib=11: 3x10
                for k in 0:2
                    wv |= bmask(queue[k + 1], 10) << (20 - 10 * k)
                end
                nibbles[widx + 1] = 0x2
            elseif length(queue) >= 2 && all(bws[1:2] .<= 15)
                n = 2
                wv = UInt32(0x2) << 30                    # dnib=10: 2x15
                wv |= bmask(queue[1], 15) << 15
                wv |= bmask(queue[2], 15)
                nibbles[widx + 1] = 0x2
            elseif length(queue) >= 1 && bws[1] <= 30
                n = 1
                wv = UInt32(0x1) << 30                    # dnib=01: 1x30
                wv |= bmask(queue[1], 30)
                nibbles[widx + 1] = 0x2
            else
                return nothing
            end
            put32!(frames, fo + 4 * widx, wv, be)
            deleteat!(queue, 1:n)
            deleteat!(bws, 1:n)
            outputsamples += n
        end
        # Nibble word
        w0 = UInt32(0)
        for widx in 0:15
            w0 |= nibbles[widx + 1] << (30 - 2 * widx)
        end
        put32!(frames, fo, w0, be)
    end
    outputsamples >= nsamp || return nothing
    # Xn = last sample, frame 0 word 2
    put32!(frames, 4 + 4, reinterpret(UInt32, samples[nsamp]), be)
    return frames
end

"""Encode samples as Steim-1 frames (transcription of msr_encode_steim1)."""
function encode_steim1(samples::Vector{Int32}, diff0::Int32, maxframes::Int, be::Bool)
    nsamp = length(samples)
    nsamp == 0 && return UInt8[]
    maxframes == 0 && return nothing
    frames = Vector{UInt8}(undef, 64 * maxframes)
    fill!(frames, 0x00)
    bw0 = bitwidth(diff0)
    bw0 > 32 && return nothing
    queue = Int32[diff0]
    bws = Int[bw0]
    qi = 1
    outputsamples = 0
    for frameidx in 0:(maxframes - 1)
        outputsamples >= nsamp && break
        fo = 64 * frameidx
        nibbles = zeros(UInt32, 16)
        startw = frameidx == 0 ? 3 : 1
        if frameidx == 0
            put32!(frames, fo + 4, reinterpret(UInt32, samples[1]), be)
        end
        for widx in startw:15
            outputsamples >= nsamp && break
            while length(queue) < 4 && qi <= nsamp - 1
                d = wsub(samples[qi + 1], samples[qi])
                push!(queue, d); push!(bws, bitwidth(d))
                qi += 1
            end
            isempty(queue) && break
            local n, wv
            if length(queue) >= 4 && all(bws[1:4] .<= 8)
                n = 4
                wv = UInt32(0)
                for k in 0:3
                    sh = be ? (24 - 8 * k) : (8 * k)
                    wv |= bmask(queue[k + 1], 8) << sh
                end
                nibbles[widx + 1] = 0x1
            elseif length(queue) >= 2 && all(bws[1:2] .<= 16)
                n = 2
                # 2x16 in file 16-bit order; mask sign extension, keep bits
                u1 = reinterpret(UInt32, queue[1]) & UInt32(0xFFFF)
                u2 = reinterpret(UInt32, queue[2]) & UInt32(0xFFFF)
                wv = be ? (u1 << 16) | u2 : u1 | (u2 << 16)
                nibbles[widx + 1] = 0x2
            else
                n = 1
                wv = reinterpret(UInt32, queue[1])
                nibbles[widx + 1] = 0x3
            end
            put32!(frames, fo + 4 * widx, wv, be)
            deleteat!(queue, 1:n)
            deleteat!(bws, 1:n)
            outputsamples += n
        end
        w0 = UInt32(0)
        for widx in 0:15
            w0 |= nibbles[widx + 1] << (30 - 2 * widx)
        end
        put32!(frames, fo, w0, be)
    end
    outputsamples >= nsamp || return nothing
    put32!(frames, 4 + 4, reinterpret(UInt32, samples[nsamp]), be)
    return frames
end

"""Test helper: build one 64-byte frame from nibble/word lists.
nibbles: 16 entries (0..3), words: Dict(word index => logical UInt32)."""
function mkframe(nibbles::Vector{Int}, words::Dict{Int,UInt32}, be::Bool; startw::Int=3)
    fr = Vector{UInt8}(undef, 64)
    fill!(fr, 0x00)
    w0 = UInt32(0)
    for widx in 0:15
        w0 |= (UInt32(nibbles[widx + 1]) << (30 - 2 * widx))
    end
    put32!(fr, 0, w0, be)
    for (widx, w) in words
        put32!(fr, 4 * widx, w, be)
    end
    return fr
end

"""Pack 4 x 8-bit diffs into a logical word (file byte order)."""
function pack4x8(d0::Int32, d1::Int32, d2::Int32, d3::Int32, be::Bool)::UInt32
    w = UInt32(0)
    for (k, d) in enumerate((d0, d1, d2, d3))
        sh = be ? (24 - 8 * (k - 1)) : (8 * (k - 1))
        w |= (reinterpret(UInt32, d) & UInt32(0xFF)) << sh
    end
    return w
end

"""Build a synthetic 512-byte v2 record. Steim payloads are encoded with the
built-in encoder. FSDH fields are chosen to exercise parser paths."""
function build_record(; samples::Vector{Int32}=Int32[],
                       encoding::Int=10, be::Bool=true, nsamp::Union{Int,Nothing}=nothing,
                       station="TEST", location="  ", channel="BHZ", network="XX",
                       year=2012, day=133, hour=0, minute=0, second=0, fsec=0,
                       rate_factor=40, rate_mult=1, actflags=0x00, ioflags=0x00,
                       dqflags=0x00, timecorr=Int32(0), usec=0, tq=100,
                       framecount=7, b100_rate=nothing, with_b1000=true,
                       with_b1001=true, seq="000001", payload_override=nothing)
    ns = nsamp === nothing ? length(samples) : nsamp
    buf = Vector{UInt8}(undef, 512)
    fill!(buf, 0x00)
    for (k, c) in enumerate(seq)
        buf[k] = UInt8(c)
    end
    buf[7] = UInt8('R')
    buf[8] = UInt8(' ')
    for (k, c) in enumerate(station)
        buf[8 + k] = UInt8(c)
    end
    for (k, c) in enumerate(location)
        buf[13 + k] = UInt8(c)
    end
    for (k, c) in enumerate(channel)
        buf[15 + k] = UInt8(c)
    end
    for (k, c) in enumerate(network)
        buf[18 + k] = UInt8(c)
    end
    put16!(buf, 20, UInt16(year), be)
    put16!(buf, 22, UInt16(day), be)
    buf[25] = UInt8(hour); buf[26] = UInt8(minute); buf[27] = UInt8(second)
    put16!(buf, 28, UInt16(fsec), be)
    put16!(buf, 30, UInt16(ns), be)
    put16!(buf, 32, reinterpret(UInt16, Int16(rate_factor)), be)
    put16!(buf, 34, reinterpret(UInt16, Int16(rate_mult)), be)
    buf[37] = actflags; buf[38] = ioflags; buf[39] = dqflags
    nblk = 0
    with_b1000 && (nblk += 1)
    with_b1001 && (nblk += 1)
    b100_rate !== nothing && (nblk += 1)
    buf[40] = UInt8(nblk)
    doff = nblk > 0 ? 48 + 8 * nblk : 64
    put32!(buf, 40, reinterpret(UInt32, timecorr), be)   # bytes 40-43
    put16!(buf, 44, UInt16(doff), be)                    # data offset
    put16!(buf, 46, UInt16(48), be)                       # blockette offset
    boff = 48
    if with_b1000
        put16!(buf, boff, UInt16(1000), be)
        put16!(buf, boff + 2, UInt16(with_b1001 ? boff + 8 : (b100_rate !== nothing ? boff + 8 : 0)), be)
        buf[boff + 5] = UInt8(encoding)
        buf[boff + 6] = UInt8(be ? 1 : 0)
        buf[boff + 7] = UInt8(9)                             # 512 = 2^9
        boff += 8
    end
    if with_b1001
        put16!(buf, boff, UInt16(1001), be)
        put16!(buf, boff + 2, UInt16(b100_rate !== nothing ? boff + 8 : 0), be)
        buf[boff + 5] = UInt8(tq)
        buf[boff + 6] = reinterpret(UInt8, Int8(usec))
        buf[boff + 8] = UInt8(framecount)
        boff += 8
    end
    if b100_rate !== nothing
        put16!(buf, boff, UInt16(100), be)
        put16!(buf, boff + 2, UInt16(0), be)
        put32!(buf, boff + 4, reinterpret(UInt32, b100_rate), be)
        boff += 8
    end
    if payload_override !== nothing
        n = min(length(payload_override), 512 - doff)
        for k in 1:n
            buf[doff + k] = payload_override[k]
        end
    elseif encoding == 10 || encoding == 5 || encoding == -1
        enc = encode_steim1(samples, samples[1], 7, be)
        enc === nothing && error("test data not Steim-1 encodable")
        for k in 1:min(length(enc), 512 - doff)
            buf[doff + k] = enc[k]
        end
    elseif encoding == 11
        enc = encode_steim2(samples, samples[1], 7, be)
        enc === nothing && error("test data not Steim-2 encodable")
        for k in 1:min(length(enc), 512 - doff)
            buf[doff + k] = enc[k]
        end
    elseif encoding == 1
        for k in 1:length(samples)
            put16!(buf, doff + 2 * (k - 1), reinterpret(UInt16, Int16(samples[k])), be)
        end
    elseif encoding == 2 || encoding == 3
        for k in 1:length(samples)
            put32!(buf, doff + 4 * (k - 1), reinterpret(UInt32, samples[k]), be)
        end
    end
    return buf
end

# ---- xorshift64 fuzz RNG (deterministic) ----
const FZ = Ref(UInt64(0x9E3779B97F4A7C15))
function rnd()::UInt64
    x = FZ[]
    x ⊻= x << 13
    x ⊻= x >> 7
    x ⊻= x << 17
    FZ[] = x
    return x
end
r32(lo::Int32, hi::Int32)::Int32 =
    Int32(Int64(lo) + Int64(rnd() % UInt64(UInt64(Int64(hi) - Int64(lo) + 1))))

function run_selftest()::Bool
    TC[] = 0; TF[] = 0
    println("miniseed.jl self test")

    # --- binary helpers -----------------------------------------------------
    buf = Vector{UInt8}([0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08])
    checkeq(u16(buf, 0, true), UInt16(0x0102), "u16 BE")
    checkeq(u16(buf, 0, false), UInt16(0x0201), "u16 LE")
    checkeq(u32(buf, 0, true), UInt32(0x01020304), "u32 BE")
    checkeq(u32(buf, 0, false), UInt32(0x04030201), "u32 LE")
    checkeq(i32(buf, 0, true), Int32(0x01020304), "i32 BE")
    checkeq(i32(buf, 0, false), reinterpret(Int32, UInt32(0x04030201)), "i32 LE")
    checkeq(f32(buf, 0, false), reinterpret(Float32, UInt32(0x04030201)), "f32 LE")
    checkeq(f64(buf, 0, true), reinterpret(Float64, UInt64(0x0102030405060708)), "f64 BE")
    sbuf = Vector{UInt8}([0x41, 0x42, 0x00, 0x43])
    checkeq(fstr(sbuf, 0, 4), "AB", "fstr stops at NUL")
    pbuf = Vector{UInt8}(undef, 8)
    put32!(pbuf, 0, UInt32(0xDEADBEEF), true)
    checkeq(u32(pbuf, 0, true), UInt32(0xDEADBEEF), "put32/u32 BE roundtrip")
    put32!(pbuf, 0, UInt32(0xDEADBEEF), false)
    checkeq(u32(pbuf, 0, false), UInt32(0xDEADBEEF), "put32/u32 LE roundtrip")
    checkeq(pbuf[1], UInt8(0xEF), "put32 LE masks bytes (UInt8 check)")
    checkeq(signext(UInt32(0x0F), 4), Int32(-1), "signext 4-bit negative")
    checkeq(signext(UInt32(0x07), 4), Int32(7), "signext 4-bit positive")
    checkeq(signext(UInt32(0x1FFFFFFF), 30), Int32(536870911), "signext 30-bit max")
    checkeq(signext(UInt32(0x3FFFFFFF), 30), Int32(-1), "signext 30-bit all ones")
    checkeq(signext(UInt32(0x20000000), 30), Int32(-536870912), "signext 30-bit min")
    checkeq(wadd(Int32(0x7FFFFFFF), Int32(1)), typemin(Int32), "wadd wraps")
    checkeq(wsub(typemin(Int32), Int32(1)), typemax(Int32), "wsub wraps")
    checkeq(bmask(Int32(-1), 8), UInt32(0xFF), "bmask")

    # --- nomsamprate ---------------------------------------------------------
    checkeq(nomsamprate(40, 1), 40.0, "nomsamprate 40/1")
    checkeq(nomsamprate(-40, 1), 0.025, "nomsamprate -40/1")
    checkeq(nomsamprate(1, -2), 0.5, "nomsamprate 1/-2")
    checkeq(nomsamprate(0, 0), 0.0, "nomsamprate 0/0")
    checkeq(nomsamprate(1, 40), 40.0, "nomsamprate 1/40")

    # --- Steim-2 decoder unit tests (hand-built frames) -----------------------
    # Frame 0 semantics: X0 is output directly, diff[1] of the frame (the
    # first collected difference) is skipped, Xn must match the last sample.
    for be in (true, false)
        # w3: nibble 10/dnib 10 (2x15: 99 thrown away, +3);
        # w4: nibble 01 (4x8: +1, -1, +2, 0)
        w3 = UInt32(0x2) << 30 | (UInt32(99) << 15) | UInt32(3)
        w4 = pack4x8(Int32(1), Int32(-1), Int32(2), Int32(0), be)
        fr = mkframe([0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     Dict(1 => UInt32(100), 2 => UInt32(105), 3 => w3, 4 => w4), be)
        s, xnbad = decode_steim2(fr, 5, be)
        # apply +3, +1, -1, +2 -> 100, 103, 104, 103, 105
        checkeq(s, Int32[100, 103, 104, 103, 105], "steim2 unit ($(be ? "BE" : "LE"))")
        check(!xnbad, "steim2 unit Xn ok ($(be ? "BE" : "LE"))")
        # Xn mismatch detection
        fr2 = mkframe([0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                      Dict(1 => UInt32(100), 2 => UInt32(999), 3 => w3, 4 => w4), be)
        _, xnbad2 = decode_steim2(fr2, 5, be)
        check(xnbad2, "steim2 unit Xn mismatch ($(be ? "BE" : "LE"))")
        # nibble 11 / dnib 01: 6x5
        w6 = UInt32(0x1) << 30
        d5 = (Int32(-1), Int32(2), Int32(-3), Int32(4), Int32(-5), Int32(6))
        for (k, d) in enumerate(d5)
            w6 |= bmask(d, 5) << (25 - 5 * (k - 1))
        end
        fr3 = mkframe([0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                      Dict(1 => UInt32(0), 2 => UInt32(4), 3 => w6), be)
        # skip -1, apply 2, -3, 4, -5, 6 -> 0, 2, -1, 3, -2, 4
        s3, xn3 = decode_steim2(fr3, 6, be)
        checkeq(s3, Int32[0, 2, -1, 3, -2, 4], "steim2 6x5 ($(be ? "BE" : "LE"))")
        check(!xn3, "steim2 6x5 Xn ($(be ? "BE" : "LE"))")
        # nibble 11 / dnib 10: 7x4
        w7 = UInt32(0x2) << 30
        d4 = (Int32(7), Int32(-8), Int32(1), Int32(2), Int32(-3), Int32(4), Int32(-5))
        for (k, d) in enumerate(d4)
            w7 |= bmask(d, 4) << (24 - 4 * (k - 1))
        end
        fr4 = mkframe([0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                      Dict(1 => UInt32(100), 2 => UInt32(91), 3 => w7), be)
        # skip +7, apply -8, +1, +2, -3, +4, -5 -> 100, 92, 93, 95, 92, 96, 91
        s4, xn4 = decode_steim2(fr4, 7, be)
        checkeq(s4, Int32[100, 92, 93, 95, 92, 96, 91], "steim2 7x4 ($(be ? "BE" : "LE"))")
        check(!xn4, "steim2 7x4 Xn ($(be ? "BE" : "LE"))")
        # nibble 10 / dnib 00: undefined, must throw
        frb = mkframe([0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                      Dict(1 => UInt32(0), 2 => UInt32(0), 3 => UInt32(1)), be)
        threw = try
            decode_steim2(frb, 3, be); false
        catch e
            e isa SteimError
        end
        check(threw, "steim2 impossible dnib throws ($(be ? "BE" : "LE"))")
    end

    # --- Steim-1 decoder unit tests -------------------------------------------
    for be in (true, false)
        u1 = reinterpret(UInt32, Int32(2)) & UInt32(0xFFFF)
        u2 = reinterpret(UInt32, Int32(-2)) & UInt32(0xFFFF)
        w2x16 = be ? (u1 << 16) | u2 : u1 | (u2 << 16)   # +2, -2 in file order
        fr = mkframe([0, 0, 0, 2, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     Dict(1 => UInt32(50), 2 => UInt32(62),
                          3 => w2x16,
                          4 => reinterpret(UInt32, Int32(10)),
                          5 => pack4x8(Int32(1), Int32(1), Int32(1), Int32(1), be)), be)
        s, xnbad = decode_steim1(fr, 7, be)
        # frame0 diffs: (+2, -2, +10, +1, +1, +1, +1); first diff skipped:
        # out: 50, 48, 58, 59, 60, 61, 62
        checkeq(s, Int32[50, 48, 58, 59, 60, 61, 62], "steim1 unit ($(be ? "BE" : "LE"))")
        check(!xnbad, "steim1 unit Xn ok ($(be ? "BE" : "LE"))")
    end

    # --- encoder round-trips: every packing form -------------------------------
    for be in (true, false)
        # 7x4 (diffs within +-7)
        s = Int32[0; accumulate(wadd, Int32[1, -2, 3, -4, 5, -6, 7, 1, 2, -1])]
        enc = encode_steim2(s, s[2], 7, be)
        d, xnbad = decode_steim2(enc, length(s), be)
        checkeq(d, s, "steim2 7x4 roundtrip ($(be ? "BE" : "LE"))")
        check(!xnbad, "steim2 7x4 roundtrip Xn ($(be ? "BE" : "LE"))")
        # 6x5 / 5x6 / 4x8 / 3x10 / 2x15 / 1x30 by diff magnitude
        mags = Dict("6x5" => Int32(15), "5x6" => Int32(31), "4x8" => Int32(127),
                    "3x10" => Int32(511), "2x15" => Int32(16383), "1x30" => Int32(2^29 - 1))
        for (form, m) in mags
            diffs = Int32[m, -m, m, -m, m]
            s2 = Int32[7; accumulate(wadd, diffs)]
            enc2 = encode_steim2(s2, s2[2], 7, be)
            check(enc2 !== nothing, "steim2 $form encodable ($(be ? "BE" : "LE"))")
            enc2 === nothing && continue
            d2, xn2 = decode_steim2(enc2, length(s2), be)
            checkeq(d2, s2, "steim2 $form roundtrip ($(be ? "BE" : "LE"))")
            check(!xn2, "steim2 $form roundtrip Xn ($(be ? "BE" : "LE"))")
        end
        # Steim-1 forms
        s3 = Int32[100; accumulate(wadd, Int32[1, -1, 100, -100, 30000, -30000, 5])]
        enc3 = encode_steim1(s3, s3[2], 7, be)
        d3, xn3 = decode_steim1(enc3, length(s3), be)
        checkeq(d3, s3, "steim1 roundtrip ($(be ? "BE" : "LE"))")
        check(!xn3, "steim1 roundtrip Xn ($(be ? "BE" : "LE"))")
        # unencodable: diff needing 31+ bits
        s4 = Int32[0, Int32(0x40000000)]
        checkeq(encode_steim2(s4, s4[2], 7, be), nothing, "steim2 rejects 31-bit diff ($(be ? "BE" : "LE"))")
    end

    # --- fuzz: encoder vs decoder, 600 trials per (kind, steim, endianness) ---
    kinds = [(Int32(-1), Int32(1)), (Int32(-7), Int32(7)), (Int32(-31), Int32(31)),
             (Int32(-127), Int32(127)), (Int32(-511), Int32(511)),
             (Int32(-16383), Int32(16383))]
    skipped = 0
    failures = 0
    for be in (true, false), steim2 in (true, false), (lo, hi) in kinds
        for _ in 1:600
            n = Int(rnd() % 64) + 8
            base = r32(Int32(-2^20), Int32(2^20))
            diffs = [r32(lo, hi) for _ in 1:n]
            samples = Int32[base; accumulate(wadd, diffs)]
            enc = steim2 ? encode_steim2(samples, samples[2], 7, be) :
                          encode_steim1(samples, samples[2], 7, be)
            if enc === nothing
                skipped += 1
                continue
            end
            d, xnbad = steim2 ? decode_steim2(enc, length(samples), be) :
                                decode_steim1(enc, length(samples), be)
            if d != samples || xnbad
                failures += 1
            end
        end
    end
    checkeq(failures, 0, "fuzz round-trips exact (skipped $skipped unencodable)")
    check(skipped <= 24, "fuzz skip count plausible (skipped=$skipped)")

    # --- record parsing -------------------------------------------------------
    samples = Int32[100; accumulate(wadd, [r32(Int32(-127), Int32(127)) for _ in 1:120])]
    rec_be = build_record(samples=samples, encoding=10, be=true, usec=-25,
                          timecorr=Int32(1234), actflags=0x00)
    recs = scan_ms2_records(rec_be)
    checkeq(length(recs), 1, "scan finds 1 record")
    r = recs[1][1]
    checkeq(r.seq, "000001", "FSDH seq")
    checkeq(channelid(r), "XX.TEST..BHZ", "FSDH channel id")
    checkeq(r.encoding, 10, "B1000 encoding")
    checkeq(r.b1000_byteorder, 1, "B1000 byte order")
    checkeq(r.b1000_reclen, 512, "B1000 record length")
    checkeq(r.b1001_usec, -25, "B1001 usec")
    checkeq(r.b1001_tq, 100, "B1001 timing quality")
    checkeq(r.b1001_framecount, 7, "B1001 frame count")
    checkeq(r.rate_factor, 40, "rate factor")
    checkeq(samprate(r), 40.0, "samprate from factor/mult")
    checkeq(r.timecorr, Int32(1234), "timecorr")
    checkeq(r.dataoffset, 64, "data offset")
    checkeq(r.nsamples, length(samples), "nsamples")
    # BE/LE header detection
    rec_le = build_record(samples=samples, encoding=5, be=false)
    rle = scan_ms2_records(rec_le)[1][1]
    check(!rle.be, "LE header byte order detected")
    checkeq(rle.year, 2012, "LE header year")
    # B100 sample rate override
    rec100 = build_record(samples=samples, encoding=10, be=true, b100_rate=Float32(99.5))
    r100 = scan_ms2_records(rec100)[1][1]
    checkeq(samprate(r100), 99.5, "B100 sample rate override")
    # no B1000: grid-scan record length + assumed Steim-1
    nob = build_record(samples=samples, encoding=10, be=true, with_b1000=false)
    file2 = vcat(nob, nob)
    recs2 = scan_ms2_records(file2)
    checkeq(length(recs2), 2, "2 records scanned")
    rn = recs2[1][1]
    checkeq(rn.reclen, 512, "grid-scan record length")
    checkeq(rn.encoding, -1, "no B1000 -> encoding -1")
    res = decode_record(rn, recs2[1][2])
    checkeq(res.samples, samples, "no B1000 -> Steim-1 assumed, decode ok")
    # truncated record
    trunc = rec_be[1:300]
    rt = scan_ms2_records(trunc)
    check(length(rt) == 1 && rt[1][1].truncated, "truncated record flagged")
    # gap before a record is skipped
    gap = vcat(fill(0xAB, 128), rec_be)
    rg = scan_ms2_records(gap)
    checkeq(length(rg), 1, "non-SEED gap skipped on 64-byte grid")
    checkeq(rg[1][1].seq, "000001", "record after gap found")

    # --- time handling ----------------------------------------------------------
    t0 = scan_ms2_records(rec_be)[1][1]
    # 2012 day 133 = May 12; fsec 0, usec -25, timecorr 0.1234 s not applied (flag clear)
    us = start_us(t0)
    expect = (Int64(Dates.value(DateTime(2012, 5, 12) - DateTime(1970, 1, 1))) * 1000
              - 25 + 123400)
    checkeq(us, expect, "start_us: fsec/usec/timecorr")
    checkeq(iso_us(t0), "2012-05-12T00:00:00.123375", "iso_us format")
    # time correction already applied (actflags 0x02)
    rec_applied = build_record(samples=samples, encoding=10, be=true, usec=-25,
                               timecorr=Int32(1234), actflags=0x02)
    ra = scan_ms2_records(rec_applied)[1][1]
    us2 = start_us(ra)
    checkeq(us2, us - 123400, "start_us: correction applied flag honored")
    us3 = start_us(ra, corrected=false)
    checkeq(us3, us - 123400, "start_us: --no-correct equivalence")
    # leap second rolls into next second
    leap = build_record(samples=samples, encoding=10, be=true, second=60, minute=59)
    rl = scan_ms2_records(leap)[1][1]
    checkeq(iso_us(rl), "2012-05-12T01:00:00.000000", "leap second monotonic")

    # --- dual convention decode ---------------------------------------------------
    # enc=10 with Steim-1 data: auto falls back to Steim-1 with warning
    res = decode_record(r, recs[1][2]; conv=:auto)
    checkeq(res.encoding_used, "STEIM1", "enc=10 auto: falls back to STEIM1")
    check(res.samples == samples && res.xn_ok, "enc=10 auto: Steim-1 samples correct")
    check(any(w -> occursin("libmseed 3 convention", w), res.warnings),
          "enc=10 auto: fallback warning present")
    res = decode_record(r, recs[1][2]; conv=:libmseed3)
    checkeq(res.encoding_used, "STEIM1", "enc=10 libmseed3: direct STEIM1")
    res = decode_record(r, recs[1][2]; conv=:seed24)
    checkeq(res.encoding_used, "STEIM2", "enc=10 seed24: direct STEIM2")
    check(!res.xn_ok, "enc=10 seed24: Steim-1 data as Steim-2 fails Xn")
    # enc=11 with Steim-2 data decodes in auto and libmseed3
    rec11 = build_record(samples=samples, encoding=11, be=true)
    r11, p11 = scan_ms2_records(rec11)[1]
    res11 = decode_record(r11, p11; conv=:auto)
    checkeq(res11.encoding_used, "STEIM2", "enc=11 auto: STEIM2")
    check(res11.samples == samples, "enc=11 auto: samples correct")
    # enc=5 with float64 data: Steim fails -> FLOAT64 fallback
    payload = UInt8[]
    for v in Float64[1.0, 2.0, 3.0, -4.5, 10.0]
        w = reinterpret(UInt64, v)
        append!(payload, [UInt8((w >> 56) & 0xFF), UInt8((w >> 48) & 0xFF),
                          UInt8((w >> 40) & 0xFF), UInt8((w >> 32) & 0xFF),
                          UInt8((w >> 24) & 0xFF), UInt8((w >> 16) & 0xFF),
                          UInt8((w >> 8) & 0xFF), UInt8(w & 0xFF)])
    end
    rec5 = build_record(encoding=5, be=true, nsamp=5, payload_override=payload)
    r5, p5 = scan_ms2_records(rec5)[1]
    res5 = decode_record(r5, p5; conv=:auto)
    checkeq(res5.encoding_used, "FLOAT64", "enc=5 auto: FLOAT64 fallback")
    check(res5.samples == [1.0, 2.0, 3.0, -4.5, 10.0], "enc=5 auto: float64 samples")
    # enc=3 with int data whose float interpretation is implausible
    rec3 = build_record(samples=samples[1:50], encoding=3, be=true)
    r3, p3 = scan_ms2_records(rec3)[1]
    res3 = decode_record(r3, p3; conv=:auto)
    checkeq(res3.encoding_used, "INT32", "enc=3 auto: implausible float -> INT32")
    # enc=4, payload too small for float64 -> FLOAT32
    f32data = Float32[1.5, -2.5, 3.25]
    payload4 = UInt8[]
    for v in f32data
        w = reinterpret(UInt32, v)
        append!(payload4, [UInt8((w >> 24) & 0xFF), UInt8((w >> 16) & 0xFF),
                           UInt8((w >> 8) & 0xFF), UInt8(w & 0xFF)])
    end
    rec4 = build_record(encoding=4, be=true, nsamp=100, payload_override=payload4)
    r4, p4 = scan_ms2_records(rec4)[1]
    res4 = decode_record(r4, p4; conv=:auto)
    checkeq(res4.encoding_used, "FLOAT32", "enc=4 auto: size forces FLOAT32")
    check(res4.samples[1:3] == f32data, "enc=4 auto: float32 samples")
    # enc=4, float32 data that would fit as float64 -> plausibility forces FLOAT32
    f32big = Float32[-1.8728645e8, -2.0122045e8, -1.3588125e8]
    payload4b = UInt8[]
    for v in f32big
        w = reinterpret(UInt32, v)
        append!(payload4b, [UInt8((w >> 24) & 0xFF), UInt8((w >> 16) & 0xFF),
                            UInt8((w >> 8) & 0xFF), UInt8(w & 0xFF)])
    end
    append!(payload4b, fill(0x00, 12))                # 24 bytes: fits 3 float64
    rec4b = build_record(encoding=4, be=true, nsamp=3, payload_override=payload4b)
    r4b, p4b = scan_ms2_records(rec4b)[1]
    res4b = decode_record(r4b, p4b; conv=:auto)
    checkeq(res4b.encoding_used, "FLOAT32", "enc=4 auto: implausible float64 -> FLOAT32")
    check(res4b.samples == f32big, "enc=4 auto: implausible float64 samples")
    # enc=1 int16 and enc=2 int32
    rec1 = build_record(samples=samples[1:20], encoding=1, be=true)
    r1, p1 = scan_ms2_records(rec1)[1]
    res1 = decode_record(r1, p1)
    checkeq(res1.encoding_used, "INT16", "enc=1: INT16")
    check(res1.samples == samples[1:20], "enc=1: samples")
    rec2 = build_record(samples=samples[1:20], encoding=2, be=true)
    r2b, p2b = scan_ms2_records(rec2)[1]
    res2 = decode_record(r2b, p2b)
    checkeq(res2.encoding_used, "INT32", "enc=2: INT32")
    check(res2.samples == samples[1:20], "enc=2: samples")
    # text
    text = collect(codeunits("Hello miniSEED"))
    rect = build_record(encoding=0, be=true, nsamp=length(text),
                        payload_override=text)
    rt2, pt2 = scan_ms2_records(rect)[1]
    rest = decode_record(rt2, pt2)
    checkeq(strip(String(res2 isa DecodeResult ? rest.samples : "")), "Hello miniSEED", "enc=0: TEXT")

    # --- glob filter -----------------------------------------------------------
    checkeq(match(glob_to_re("XX.TEST..BHZ"), "XX.TEST..BHZ") !== nothing, true, "glob exact")
    checkeq(match(glob_to_re("*.?HZ"), "XX.TEST..BHZ") !== nothing, true, "glob wildcards")
    checkeq(match(glob_to_re("XX.*.*.BHZ"), "XX.TEST..BHZ") !== nothing, true, "glob star")
    checkeq(match(glob_to_re("XX.TEST..BHZ"), "XX.TEST..BHN") === nothing, true, "glob nonmatch")

    @printf("\n%d checks, %d failed\n", TC[], TF[])
    return TF[] == 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
