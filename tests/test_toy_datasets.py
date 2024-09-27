from sopa.utils.data import blobs, toy_dataset


def test_make_toy_dataset():
    assert toy_dataset(length=512) is not None


def test_make_blobs():
    assert blobs(length=512) is not None
