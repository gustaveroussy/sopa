from sopa.utils.data import blobs, uniform


def test_make_uniform():
    assert uniform(length=512) is not None


def test_make_blobs():
    assert blobs(length=512) is not None
