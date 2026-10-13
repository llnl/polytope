import polytope

def test_plc():
    plc = polytope.PLC2d()
    plc.facets = ((0, 1), [1, 2], (2, 3), [3, 0])
    plc.holes = [
        ((4, 5), [5, 6], (6, 7), [7, 4]),
    ]
    assert plc.valid()
    assert len(plc.facets) == 4
    assert len(plc.holes) == 1
    assert len(plc.holes[0]) == 4

if __name__ == "__main__":
    test_plc()
