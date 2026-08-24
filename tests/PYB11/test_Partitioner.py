import polytope


def _coordinates(points):
    return [(point.x, point.y) for point in points]


def _quantized_generators():
    quantizer = polytope.Quantizer2d.instance()
    real_points = [(0.10, 0.10), (0.30, 0.70),
                   (0.55, 0.25), (0.90, 0.90)]
    quantizer.init(real_points)
    return [quantizer.quantize(polytope.Point2d(x, y)) for x, y in real_points]


def test_random_partitioner_is_deterministic():
    points = _quantized_generators()
    partitioner = polytope.RandomPartitioner2d(123456789)

    first = partitioner.computePartition(points)
    second = partitioner.computePartition(points)
    assert _coordinates(first) == _coordinates(second)

    if polytope.Communicator.getNProcs() == 1:
        assert _coordinates(first) == _coordinates(points)


def test_lattice_partitioner_uses_quantizer_bounds():
    points = _quantized_generators()
    rank = polytope.Communicator.getRank()
    nranks = polytope.Communicator.getNProcs()
    quantizer = polytope.Quantizer2d.instance()

    partitioner = polytope.LatticePartitioner2d([nranks, 1])
    local_points = partitioner.computePartition(points)

    lower = quantizer.minBound
    upper = quantizer.maxBound

    def owner(point):
        if point.x == upper.x:
            return nranks - 1
        return min((point.x - lower.x)*nranks // (upper.x - lower.x), nranks - 1)

    expected = [point for point in points if owner(point) == rank]
    assert _coordinates(local_points) == _coordinates(expected)


if __name__ == "__main__":
    test_random_partitioner_is_deterministic()
    test_lattice_partitioner_uses_quantizer_bounds()
