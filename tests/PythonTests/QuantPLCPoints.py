import polytope


def test_caller_owned_points():
    quantizer = polytope.Quantizer2d()
    quantizer.init([(0.0, 0.0), (4.0, 4.0)], 0.1)
    points = polytope.QuantPLC2d.quantizePoints(
        quantizer,
        [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (2.0, 2.0)],
    )
    plc = polytope.PLC2d()
    plc.facets = [[0, 1], [1, 2], [2, 3], [3, 0]]
    qplc = polytope.QuantPLC2d(quantizer)
    qplc.init(points, plc)

    assert len(points) == 5
    assert qplc.within((1.0, 1.0), points)
    qplc.reduce(points)
    assert len(points) == 4
    assert qplc.within((1.0, 1.0), points)


if __name__ == "__main__":
    test_caller_owned_points()
