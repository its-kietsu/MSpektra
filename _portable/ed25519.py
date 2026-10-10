"""Ed25519 signatures (RFC 8032) in plain Python: the updater checks that an update was signed with the
MSpektra release key (updater.PUBLIC_KEY) before anything is installed; the release script signs with
the private key, which never leaves the release computer. Slow (about 50 ms per check) but one check per
update is all it does."""
import hashlib

_P = 2 ** 255 - 19
_L = 2 ** 252 + 27742317777372353535851937790883648493
_D = -121665 * pow(121666, _P - 2, _P) % _P
_I = pow(2, (_P - 1) // 4, _P)


def _inv(x):
    return pow(x, _P - 2, _P)


def _xrecover(y):
    xx = (y * y - 1) * _inv(_D * y * y + 1) % _P
    x = pow(xx, (_P + 3) // 8, _P)
    if (x * x - xx) % _P:
        x = x * _I % _P
    if (x * x - xx) % _P:
        raise ValueError("not a point of the curve")
    return x


_BY = 4 * _inv(5) % _P
_BX = _xrecover(_BY)
if _BX & 1:
    _BX = _P - _BX
_B = (_BX, _BY, 1, _BX * _BY % _P)  # extended coordinates


def _add(p, q):
    x1, y1, z1, t1 = p
    x2, y2, z2, t2 = q
    a = (y1 - x1) * (y2 - x2) % _P
    b = (y1 + x1) * (y2 + x2) % _P
    c = t1 * 2 * _D * t2 % _P
    d = z1 * 2 * z2 % _P
    e, f, g, h = b - a, d - c, d + c, b + a
    return (e * f % _P, g * h % _P, f * g % _P, e * h % _P)


def _mul(s, p):
    q = (0, 1, 1, 0)
    while s:
        if s & 1:
            q = _add(q, p)
        p = _add(p, p)
        s >>= 1
    return q


def _equal(p, q):
    return (p[0] * q[2] - q[0] * p[2]) % _P == 0 and (p[1] * q[2] - q[1] * p[2]) % _P == 0


def _encode(p):
    zi = _inv(p[2])
    x, y = p[0] * zi % _P, p[1] * zi % _P
    return (y | ((x & 1) << 255)).to_bytes(32, "little")


def _decode(s):
    if len(s) != 32:
        raise ValueError("a point has 32 bytes")
    y = int.from_bytes(s, "little")
    sign = y >> 255
    y &= (1 << 255) - 1
    if y >= _P:
        raise ValueError("not a point of the curve")
    x = _xrecover(y)
    if x == 0 and sign:
        raise ValueError("not a point of the curve")
    if (x & 1) != sign:
        x = _P - x
    return (x, y, 1, x * y % _P)


def _expand(seed):
    h = hashlib.sha512(seed).digest()
    a = int.from_bytes(h[:32], "little")
    a &= (1 << 254) - 8
    a |= 1 << 254
    return a, h[32:]


def public_key(seed):
    """The 32 byte public key of a 32 byte private key (seed)."""
    return _encode(_mul(_expand(seed)[0], _B))


def sign(seed, msg):
    a, prefix = _expand(seed)
    pub = _encode(_mul(a, _B))
    r = int.from_bytes(hashlib.sha512(prefix + msg).digest(), "little") % _L
    rb = _encode(_mul(r, _B))
    h = int.from_bytes(hashlib.sha512(rb + pub + msg).digest(), "little") % _L
    return rb + ((r + h * a) % _L).to_bytes(32, "little")


def verify(pub, msg, sig):
    """True only for a valid signature of msg made with the private key of pub."""
    try:
        if len(sig) != 64:
            return False
        a = _decode(pub)
        r = _decode(sig[:32])
        s = int.from_bytes(sig[32:], "little")
        if s >= _L:
            return False
        h = int.from_bytes(hashlib.sha512(sig[:32] + pub + msg).digest(), "little") % _L
        return _equal(_mul(s, _B), _add(r, _mul(h, a)))
    except (ValueError, TypeError):
        return False
