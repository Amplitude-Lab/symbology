"""Exact unrestricted weight-two word basis used by the example workflows."""
import struct


def write_word_basis(path, n):
    # B[k,i,j] = 1 when k = n*i+j, zero otherwise (zero-based indices).
    # Use the documented uncompressed WXF SparseArray CSR representation.
    def varint(v):
        out=bytearray()
        while v > 127:
            out.append((v & 127) | 128); v >>= 7
        out.append(v)
        return bytes(out)
    def symbol(s): return b's' + varint(len(s)) + s.encode('ascii')
    def function(s,k): return b'f' + varint(k) + symbol(s)
    def array(dims, values):
        return b'\xc1\x02' + varint(len(dims)) + b''.join(varint(d) for d in dims) + b''.join(struct.pack('<i', v) for v in values)
    count=n*n
    path.write_bytes(b'8:' + function('SparseArray',4) + symbol('Automatic') +
        array([3],[count,n,n]) + b'C\0' + function('List',3) + b'C\x01' +
        function('List',2) + array([count+1],range(count+1)) +
        array([count,2],(v for i in range(n) for j in range(n) for v in (i+1,j+1))) +
        array([count],[1]*count))
