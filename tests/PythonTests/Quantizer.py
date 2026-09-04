import polytope

def test_quantizer():
    q2 = polytope.Quantizer2d.instance()
    q2.init([(0.0, 0.0), (1.0, 1.0)])
    assert q2.m_init

    point = polytope.Point2d(0.5, 0.5)
    qpoint = q2.quantize(point)
    assert isinstance(qpoint.x, int)
    assert q2.inBounds(point)
    assert q2.inQBounds(qpoint)

    ckey = q2.encode(qpoint)
    assert isinstance(ckey, int)
    uckey = q2.decode(ckey)
    assert (uckey.x, uckey.y) == (qpoint.x, qpoint.y)
    rpoint = q2.decodeAndDequantize(ckey)
    assert isinstance(rpoint.x, float)


if __name__ == "__main__":
    test_quantizer()
