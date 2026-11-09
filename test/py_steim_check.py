#!/usr/bin/env python3
"""Independent Python transcription of libmseed's msr_decode_steim1 and
msr_decode_steim2 (libmseed unpackdata.c). Pure stdlib, no dependencies.

Used to cross-check miniseed.jl against a second, independently written
decoder on the official libmseed reference test data.

Byte-order semantics follow the C code exactly:
- The nibble word and all multi-byte sub-values are read as logical 32-bit
  words in the data byte order (B1000 byte order flag).
- 4x8-bit differences are read in FILE byte order (the C code reads d8[]
  memory bytes without swapping).
- 2x16-bit differences (Steim1) are logical 16-bit values in file order.
- Frame 0: word 1 is X0 (output directly), word 2 is Xn (integrity check),
  data words start at 3; the first decoded difference of frame 0 is ignored.
- Accumulation is modulo 2^32.
- Xn is compared only when the decoded count equals the expected count.
"""
import struct


def u32(buf, off, be):
    return struct.unpack_from('>I' if be else '<I', buf, off)[0]


def signext(v, bits):
    m = (1 << bits) - 1
    x = v & m
    return x - (1 << bits) if x & (1 << (bits - 1)) else x


def _word_diffs_steim2(w, nib, be):
    """Differences encoded in one 32-bit word (logical value) for Steim2."""
    if nib == 0:
        return []
    if nib == 1:  # four 8-bit differences, file byte order
        bs = [(w >> 24) & 0xFF, (w >> 16) & 0xFF, (w >> 8) & 0xFF, w & 0xFF] if be \
            else [w & 0xFF, (w >> 8) & 0xFF, (w >> 16) & 0xFF, (w >> 24) & 0xFF]
        return [signext(b, 8) for b in bs]
    if nib == 2:  # consult dnib (high two bits of logical word)
        dnib = (w >> 30) & 0x3
        if dnib == 1:
            return [signext(w & 0x3FFFFFFF, 30)]
        if dnib == 2:
            return [signext((w >> 15) & 0x7FFF, 15), signext(w & 0x7FFF, 15)]
        if dnib == 3:
            return [signext((w >> 20) & 0x3FF, 10),
                    signext((w >> 10) & 0x3FF, 10),
                    signext(w & 0x3FF, 10)]
        raise ValueError('Steim2: undefined dnib=00 for nibble=10')
    if nib == 3:  # consult dnib
        dnib = (w >> 30) & 0x3
        if dnib == 0:
            return [signext((w >> (24 - 6 * k)) & 0x3F, 6) for k in range(5)]
        if dnib == 1:
            return [signext((w >> (25 - 5 * k)) & 0x1F, 5) for k in range(6)]
        if dnib == 2:
            return [signext((w >> (24 - 4 * k)) & 0xF, 4) for k in range(7)]
        raise ValueError('Steim2: undefined dnib=11 for nibble=11')
    raise ValueError('Steim2: bad nibble')


def _word_diffs_steim1(w, nib, be):
    """Differences encoded in one 32-bit word (logical value) for Steim1."""
    if nib == 0:
        return []
    if nib == 1:  # four 8-bit differences, file byte order
        bs = [(w >> 24) & 0xFF, (w >> 16) & 0xFF, (w >> 8) & 0xFF, w & 0xFF] if be \
            else [w & 0xFF, (w >> 8) & 0xFF, (w >> 16) & 0xFF, (w >> 24) & 0xFF]
        return [signext(b, 8) for b in bs]
    if nib == 2:  # two 16-bit differences, logical 16-bit values in file order
        return [signext((w >> 16) & 0xFFFF, 16), signext(w & 0xFFFF, 16)] if be \
            else [signext(w & 0xFFFF, 16), signext((w >> 16) & 0xFFFF, 16)]
    if nib == 3:
        return [signext(w, 32)]
    raise ValueError('Steim1: bad nibble')


def _decode_steim(data, nsamples, be, worddiffs):
    maxframes = len(data) // 64
    samples = []
    xn = None
    for frameidx in range(maxframes):
        if len(samples) >= nsamples and frameidx > 0:
            break
        base = frameidx * 64
        nibble_word = u32(data, base, be)
        start = 3 if frameidx == 0 else 1
        diffs = []
        for widx in range(start, 16):
            nib = (nibble_word >> (30 - 2 * widx)) & 0x3
            if nib == 0:
                continue
            diffs += worddiffs(u32(data, base + 4 * widx, be), nib, be)
        if frameidx == 0:
            x0 = signext(u32(data, base + 4, be), 32)
            xn = signext(u32(data, base + 8, be), 32)
            samples.append(x0)
            diffs = diffs[1:]  # first difference of frame 0 is ignored
        prev = samples[-1]
        for d in diffs:
            if len(samples) >= nsamples:
                break
            v = (prev + d) & 0xFFFFFFFF
            v = v - (1 << 32) if v & 0x80000000 else v
            samples.append(v)
            prev = v
    xn_ok = (len(samples) == nsamples and samples[-1] == xn)
    return samples, xn_ok


def decode_steim1(data, nsamples, be):
    return _decode_steim(data, nsamples, be, _word_diffs_steim1)


def decode_steim2(data, nsamples, be):
    return _decode_steim(data, nsamples, be, _word_diffs_steim2)


def records(path, fixed_reclen=None):
    """Minimal miniSEED v2 record splitter: detects header byte order,
    walks the blockette chain for B1000, returns (enc, nsamp, payload, be)."""
    data = open(path, 'rb').read()
    recs = []
    off = 0
    while off + 48 <= len(data):
        rec = data[off:off + (fixed_reclen or 512)]
        if len(rec) < 48:
            break
        year_be = struct.unpack_from('>H', rec, 20)[0]
        hbe = 1900 <= year_be <= 2100
        fmt = '>' if hbe else '<'
        nsamp = struct.unpack_from(fmt + 'H', rec, 30)[0]
        doff = struct.unpack_from(fmt + 'H', rec, 44)[0]
        enc, dbe = -1, hbe
        boff = struct.unpack_from(fmt + 'H', rec, 46)[0]
        while boff:
            btype = struct.unpack_from(fmt + 'H', rec, boff)[0]
            nxt = struct.unpack_from(fmt + 'H', rec, boff + 2)[0]
            if btype == 1000:
                enc = rec[boff + 4]
                dbe = (rec[boff + 5] == 1)
            boff = nxt
        payload = rec[doff:]
        recs.append((enc, nsamp, payload, dbe))
        off += (fixed_reclen or 512)
    return recs


def main():
    import os
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ok = True
    for f, kind, expfile in [
            ('reference-testdata-steim1.mseed2', 1, 'expected_steim1.txt'),
            ('reference-testdata-steim1-LE.mseed2', 1, 'expected_steim1.txt'),
            ('reference-testdata-steim2.mseed2', 2, 'expected_steim2.txt'),
            ('reference-testdata-steim2-LE.mseed2', 2, 'expected_steim2.txt')]:
        expected = [int(l) for l in open(os.path.join(base, 'test', expfile)) if l.strip()]
        allsamples = []
        xnallok = True
        for enc, nsamp, payload, dbe in records(os.path.join(base, 'testdata', f)):
            dec = decode_steim1 if kind == 1 else decode_steim2
            s, xn_ok = dec(payload, nsamp, dbe)
            allsamples += s
            xnallok = xnallok and xn_ok
        match = allsamples == expected
        ok = ok and match and xnallok
        print(('PASS' if (match and xnallok) else 'FAIL'), f,
              len(allsamples), 'samples, Xn ok:', xnallok)
    print('ALL PASS' if ok else 'FAILURES')
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
