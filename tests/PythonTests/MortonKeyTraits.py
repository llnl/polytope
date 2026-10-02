import polytope


def test_morton_keys():
    point2 = polytope.CoordinatePoint2d(17, 31)
    key2 = polytope.MortonKeyTraits2d.encode(point2)
    assert isinstance(key2, int)
    roundtrip2 = polytope.MortonKeyTraits2d.decode(key2)
    assert (roundtrip2.x, roundtrip2.y) == (point2.x, point2.y)

    big = 1 << 40
    point3 = polytope.CoordinatePoint3d(big, big + 1, big + 2)
    key3 = polytope.MortonKeyTraits3d.encode(point3)
    assert isinstance(key3, int)
    assert key3.bit_length() > 64
    roundtrip3 = polytope.MortonKeyTraits3d.decode(key3)
    assert (roundtrip3.x, roundtrip3.y, roundtrip3.z) == (point3.x, point3.y, point3.z)

    encode_point3 = polytope.KeyPoint3d(key3, key3 + 1, key3 + 2)
    assert isinstance(encode_point3.x, int)
    assert encode_point3.x == key3
    assert encode_point3[2] == key3 + 2


if __name__ == "__main__":
    test_morton_keys()
