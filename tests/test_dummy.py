from sopa.utils.data import blobs, uniform


def test_make_uniform():
    assert uniform() is not None


def test_make_blobs():
    assert blobs() is not None
