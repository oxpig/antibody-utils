import pytest

from antibody_utils import _data


def test_blosum62_is_symmetric_ncbi_matrix():
    blosum62 = _data.load_blosum62()
    assert len(blosum62) == 24 * 24
    assert blosum62["W", "W"] == 11
    assert all(blosum62[a, b] == blosum62[b, a] for a, b in blosum62)


def test_missing_blosum62_explains_how_to_fetch_it(monkeypatch):
    def missing(name):
        raise FileNotFoundError(name)

    monkeypatch.setattr(_data, "_read_bytes", missing)
    _data.load_blosum62.cache_clear()
    try:
        with pytest.raises(FileNotFoundError, match=r"fetch_blosum62\.py"):
            _data.load_blosum62()
    finally:
        _data.load_blosum62.cache_clear()
