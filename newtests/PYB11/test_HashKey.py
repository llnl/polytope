import polytope


def test_hash_keys():
    point2 = polytope.CoordPoint2d(17, 31)
    key2 = polytope.HashKey2d.hash(point2)
    assert isinstance(key2, int)
    roundtrip2 = polytope.HashKey2d.unhash(key2)
    assert (roundtrip2.x, roundtrip2.y) == (point2.x, point2.y)

    big = 1 << 40
    point3 = polytope.CoordPoint3d(big, big + 1, big + 2)
    key3 = polytope.HashKey3d.hash(point3)
    assert isinstance(key3, int)
    assert key3.bit_length() > 64
    roundtrip3 = polytope.HashKey3d.unhash(key3)
    assert (roundtrip3.x, roundtrip3.y, roundtrip3.z) == (point3.x, point3.y, point3.z)

    hash_point3 = polytope.HashPoint3d(key3, key3 + 1, key3 + 2)
    assert isinstance(hash_point3.x, int)
    assert hash_point3.x == key3
    assert hash_point3[2] == key3 + 2


if __name__ == "__main__":
    test_hash_keys()
